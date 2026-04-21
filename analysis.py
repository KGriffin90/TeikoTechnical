"""
analysis.py
   --
Core data-access and analysis functions for the Loblaw Bio Miraclib trial.
"""

import functools
import os
import sqlite3

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold, cross_val_score
from statsmodels.stats.multitest import multipletests

from config import DB_PATH, POPULATIONS

# Cache helpers
def _db_mtime(db_path: str) -> float:
    """Return the file-modification time; 0.0 if the file does not exist yet."""
    try:
        return os.path.getmtime(db_path)
    except FileNotFoundError:
        return 0.0

# Helper Methods
def get_connection(db_path: str = DB_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    return conn

def fetch_df(query: str, db_path: str = DB_PATH, params=()) -> pd.DataFrame:
    with get_connection(db_path) as conn:
        return pd.read_sql_query(query, conn, params=params)

# Compute Frequency Table
@functools.lru_cache(maxsize=8)
def _compute_frequency_table_cached(db_path: str, _mtime: float) -> pd.DataFrame:
    sql = """
        SELECT external_id AS sample,
               b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte,
               (b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte) AS total_count
        FROM sample
    """
    raw = fetch_df(sql, db_path)

    records = []
    for _, row in raw.iterrows():
        for pop in POPULATIONS:
            records.append({
                "sample":      row["sample"],
                "total_count": row["total_count"],
                "population":  pop,
                "count":       int(row[pop]),
                "percentage":  round(row[pop] / row["total_count"] * 100, 4)
                               if row["total_count"] else 0,
            })
    return pd.DataFrame(records)


def compute_frequency_table(db_path: str = DB_PATH) -> pd.DataFrame:
    """
    Returns a tidy DataFrame: sample | total_count | population | count | percentage.
    Result is cached until the DB file is modified.
    """
    return _compute_frequency_table_cached(db_path, _db_mtime(db_path))

# Responder Frequencies
def get_melanoma_miraclib_pbmc(db_path: str = DB_PATH) -> pd.DataFrame:
    """Pull melanoma PBMC samples treated with miraclib (all timepoints)."""
    sql = """
        SELECT s.external_id AS sample,
               st.response,
               s.time_from_treatment_start,
               s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte,
               (s.b_cell + s.cd8_t_cell + s.cd4_t_cell + s.nk_cell + s.monocyte) AS total_count
        FROM sample s
        JOIN subject           sub ON s.subject_id  = sub.subject_id
        JOIN subject_treatment st  ON sub.subject_id = st.subject_id
        JOIN treatment_arm     ta  ON st.arm_id      = ta.arm_id
        WHERE sub.condition   = 'melanoma'
          AND ta.treatment    = 'miraclib'
          AND s.sample_type   = 'PBMC'
    """
    return fetch_df(sql, db_path)

@functools.lru_cache(maxsize=8)
def _compute_responder_frequencies_cached(db_path: str, _mtime: float) -> pd.DataFrame:
    raw = get_melanoma_miraclib_pbmc(db_path)
    if raw.empty:
        return pd.DataFrame()

    records = []
    for _, row in raw.iterrows():
        for pop in POPULATIONS:
            pct = row[pop] / row["total_count"] * 100 if row["total_count"] else 0
            records.append({
                "sample":     row["sample"],
                "response":   row["response"],
                "population": pop,
                "percentage": round(pct, 4),
            })
    return pd.DataFrame(records)

def compute_responder_frequencies(db_path: str = DB_PATH) -> pd.DataFrame:
    """
    Returns tidy DataFrame: sample | response | population | percentage
    for melanoma / miraclib / PBMC samples.
    Cached until the DB file is modified.
    """
    return _compute_responder_frequencies_cached(db_path, _db_mtime(db_path))

# Statistical Tests
def run_statistical_tests(freq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Mann-Whitney U test per population (responders vs non-responders).
    Returns: population | p_value | p_value_adj | significant | effect_size
    """
    results = []

    for pop in POPULATIONS:
        sub = freq_df[freq_df["population"] == pop]
        yes = sub[sub["response"] == "yes"]["percentage"].values
        no  = sub[sub["response"] == "no"]["percentage"].values

        if len(yes) < 2 or len(no) < 2:
            results.append({
                "population":       pop,
                "p_value":          np.nan,
                "p_value_adj":      np.nan,
                "significant":      False,
                "effect_size":      np.nan,
                "n_responders":     len(yes),
                "n_non_responders": len(no),
            })
            continue

        stat, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
        n1, n2  = len(yes), len(no)
        r       = 1 - (2 * stat) / (n1 * n2)   # rank biserial correlation

        results.append({
            "population":       pop,
            "p_value":          float(p),
            "p_value_adj":      np.nan,
            "significant":      False,
            "effect_size":      round(r, 4),
            "n_responders":     n1,
            "n_non_responders": n2,
        })

    # Benjamini Hochberg FDR correction
    valid = [i for i, r in enumerate(results) if not np.isnan(r["p_value"])]
    pvals = [results[i]["p_value"] for i in valid]
    if pvals:
        _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
        for idx, p_adj in zip(valid, pvals_adj):
            results[idx]["p_value_adj"] = round(p_adj, 6)
            results[idx]["significant"] = p_adj < 0.05

    for r in results:
        if not np.isnan(r["p_value"]):
            r["p_value"] = round(r["p_value"], 6)

    return pd.DataFrame(results).sort_values("p_value_adj")

# Longitudinal Figure
def make_longitudinal_figure():
    """Line chart of immune-cell trajectories over t = 0 / 7 / 14."""
    df = get_melanoma_miraclib_pbmc()
    if df.empty:
        return go.Figure()

    records = []
    for _, row in df.iterrows():
        total = row["total_count"]
        for pop in POPULATIONS:
            records.append({
                "sample":     row["sample"],
                "response":   row["response"],
                "time":       row["time_from_treatment_start"],
                "population": pop,
                "percentage": (row[pop] / total * 100) if total else 0,
            })

    long_df = pd.DataFrame(records)

    fig = px.line(
        long_df,
        x="time",
        y="percentage",
        color="population",
        line_group="sample",
        facet_col="response",
        markers=True,
        title="Immune Cell Dynamics Over Time",
    )
    fig.update_layout(height=500, margin=dict(t=50, l=40, r=20, b=40))
    return fig


# Subset Analysis
def avg_b_cells_melanoma_males_responders_t0(db_path: str = DB_PATH) -> float:
    """Average B-cell count: melanoma + male + responders + baseline (t=0)."""
    sql = """
        SELECT AVG(s.b_cell) AS avg_b_cells
        FROM sample s
        JOIN subject           sub ON s.subject_id  = sub.subject_id
        JOIN subject_treatment st  ON sub.subject_id = st.subject_id
        JOIN treatment_arm     ta  ON st.arm_id      = ta.arm_id
        WHERE sub.condition               = 'melanoma'
          AND sub.sex                     = 'M'
          AND st.response                 = 'yes'
          AND ta.treatment                = 'miraclib'
          AND s.sample_type               = 'PBMC'
          AND s.time_from_treatment_start = 0
    """
    df = fetch_df(sql, db_path)
    return round(df["avg_b_cells"].iloc[0], 2)


def get_baseline_melanoma_miraclib(db_path: str = DB_PATH) -> pd.DataFrame:
    """
    Baseline melanoma / miraclib / PBMC subset (t=0).
    Includes raw cell counts so build_response_model can compute frequencies.
    """
    sql = """
        SELECT s.external_id AS sample,
               sub.external_id AS subject,
               p.name AS project,
               sub.condition, sub.sex,
               ta.treatment,
               st.response,
               s.sample_type,
               s.time_from_treatment_start,
               s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte,
               (s.b_cell + s.cd8_t_cell + s.cd4_t_cell + s.nk_cell + s.monocyte)
                   AS total_count
        FROM sample s
        JOIN subject           sub ON s.subject_id  = sub.subject_id
        JOIN project           p   ON sub.project_id = p.project_id
        JOIN subject_treatment st  ON sub.subject_id = st.subject_id
        JOIN treatment_arm     ta  ON st.arm_id      = ta.arm_id
        WHERE sub.condition               = 'melanoma'
          AND s.sample_type               = 'PBMC'
          AND s.time_from_treatment_start = 0
          AND ta.treatment                = 'miraclib'
    """
    return fetch_df(sql, db_path)

def compute_subset_summaries(db_path: str = DB_PATH) -> dict:
    """
    Part 4 sub-questions on the baseline melanoma/miraclib/PBMC subset.
    Returns dict: baseline_df | samples_per_project | response_counts | sex_counts
    """
    df = get_baseline_melanoma_miraclib(db_path)

    if df.empty:
        return {
            "baseline_df":         df,
            "samples_per_project": pd.DataFrame(),
            "response_counts":     pd.DataFrame(),
            "sex_counts":          pd.DataFrame(),
        }

    samples_per_project = (
        df.groupby("project")["sample"]
          .count()
          .reset_index()
          .rename(columns={"sample": "n_samples"})
    )

    unique_subjects = df.drop_duplicates("subject")
    response_counts = (
        unique_subjects.groupby("response")["subject"]
        .count()
        .reset_index()
        .rename(columns={"subject": "n_subjects"})
    )
    sex_counts = (
        unique_subjects.groupby("sex")["subject"]
        .count()
        .reset_index()
        .rename(columns={"subject": "n_subjects"})
    )

    return {
        "baseline_df":         df,
        "samples_per_project": samples_per_project,
        "response_counts":     response_counts,
        "sex_counts":          sex_counts,
    }

# Logistic Regression  (Predictive Modelling tab)
@functools.lru_cache(maxsize=8)
def _build_response_model_cached(db_path: str, _mtime: float):
    df = get_baseline_melanoma_miraclib(db_path)
    if df.empty:
        return None

    data = df.drop_duplicates("sample").copy()
    data["total"] = data[POPULATIONS].sum(axis=1)
    for pop in POPULATIONS:
        data[pop] = data[pop] / data["total"] * 100

    data["y"] = data["response"].map({"yes": 1, "no": 0})
    data = data.dropna(subset=["y"])

    X = data[POPULATIONS].values
    y = data["y"].values

    if len(np.unique(y)) < 2:
        # Cannot train a classifier with only one class present
        return None

    model = LogisticRegression(max_iter=1000, random_state=0)

    # Cross validated AUC
    n_splits = min(5, int(np.bincount(y.astype(int)).min()))
    n_splits = max(n_splits, 2)
    cv       = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    cv_aucs  = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
    model.fit(X, y)

    importances = pd.DataFrame({
        "population":  POPULATIONS,
        "coefficient": model.coef_[0],
    }).sort_values("coefficient", ascending=False).reset_index(drop=True)

    # In sample ROC curve (for illustrative purposes only)
    y_score         = model.predict_proba(X)[:, 1]
    fpr, tpr, _     = roc_curve(y, y_score)
    insample_auc    = float(auc(fpr, tpr))

    return {
        "auc_mean":          float(np.mean(cv_aucs)),
        "auc_std":           float(np.std(cv_aucs)),
        "cv_aucs":           cv_aucs.tolist(),
        "n_splits":          n_splits,
        "insample_auc":      insample_auc,
        "fpr":               fpr.tolist(),
        "tpr":               tpr.tolist(),
        "feature_importance": importances.to_dict(orient="records"),
        "n_samples":         len(data),
        "n_responders":      int(y.sum()),
        "n_non_responders":  int((y == 0).sum()),
    }

def build_response_model(db_path: str = DB_PATH):
    """
    Fit a logistic regression model and return serialisable results.
    Cached until the DB file is modified.
    """
    return _build_response_model_cached(db_path, _db_mtime(db_path))

# CLI smoke-test
if __name__ == "__main__":
    print("=== Part 2: Frequency table (first 10 rows) ===")
    ft = compute_frequency_table()
    print(ft.head(10).to_string(index=False))

    print("\n=== Part 3: Statistical tests ===")
    freq     = compute_responder_frequencies()
    stats_df = run_statistical_tests(freq)
    print(stats_df.to_string(index=False))

    print("\n=== Part 4: Subset summaries ===")
    summaries = compute_subset_summaries()
    print(f"Baseline samples: {len(summaries['baseline_df'])}")
    print(summaries["samples_per_project"].to_string(index=False))

    print("\n=== Melanoma Male Responders (Time 0) ===")
    avg = avg_b_cells_melanoma_males_responders_t0()
    print(f"Average B cells: {avg:.2f}")