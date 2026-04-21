"""
load_data.py
    
Initialises the SQLite database schema and loads cell-count.csv.

Run:
    python3 load_data.py
"""

import os
import sys
import sqlite3
import csv
import argparse

from config import DB_PATH

DEFAULT_CSV = "data/cell-count.csv"

# Schema
SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS project (
    project_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    name         TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS subject (
    subject_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    external_id  TEXT NOT NULL UNIQUE,
    project_id   INTEGER NOT NULL REFERENCES project(project_id),
    condition    TEXT,
    age          INTEGER,
    sex          TEXT CHECK(sex IN ('M','F'))
);

CREATE TABLE IF NOT EXISTS treatment_arm (
    arm_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    treatment    TEXT NOT NULL,
    UNIQUE(treatment)
);

CREATE TABLE IF NOT EXISTS subject_treatment (
    subject_id   INTEGER NOT NULL REFERENCES subject(subject_id),
    arm_id       INTEGER NOT NULL REFERENCES treatment_arm(arm_id),
    response     TEXT CHECK(response IN ('yes','no')),
    PRIMARY KEY (subject_id, arm_id)
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    external_id                TEXT NOT NULL UNIQUE,
    subject_id                 INTEGER NOT NULL REFERENCES subject(subject_id),
    sample_type                TEXT,
    time_from_treatment_start  INTEGER,
    b_cell                     INTEGER,
    cd8_t_cell                 INTEGER,
    cd4_t_cell                 INTEGER,
    nk_cell                    INTEGER,
    monocyte                   INTEGER
);

CREATE VIEW IF NOT EXISTS v_sample_full AS
SELECT
    s.external_id            AS sample,
    sub.external_id          AS subject,
    p.name                   AS project,
    sub.condition,
    sub.age,
    sub.sex,
    ta.treatment,
    st.response,
    s.sample_type,
    s.time_from_treatment_start,
    s.b_cell,
    s.cd8_t_cell,
    s.cd4_t_cell,
    s.nk_cell,
    s.monocyte,
    (s.b_cell + s.cd8_t_cell + s.cd4_t_cell + s.nk_cell + s.monocyte) AS total_count
FROM sample s
JOIN subject           sub ON s.subject_id  = sub.subject_id
JOIN project           p   ON sub.project_id = p.project_id
JOIN subject_treatment st  ON sub.subject_id = st.subject_id
JOIN treatment_arm     ta  ON st.arm_id      = ta.arm_id;
"""

# Loader
def load_csv(csv_path: str, conn: sqlite3.Connection) -> int:
    """Parse csv_path and insert all rows; returns row count inserted."""
    cur = conn.cursor()

    project_cache: dict[str, int] = {}
    subject_cache: dict[str, int] = {}
    arm_cache:     dict[str, int] = {}

    def get_or_create_project(name: str) -> int:
        if name not in project_cache:
            cur.execute("INSERT OR IGNORE INTO project(name) VALUES(?)", (name,))
            cur.execute("SELECT project_id FROM project WHERE name=?", (name,))
            project_cache[name] = cur.fetchone()[0]
        return project_cache[name]

    def get_or_create_arm(treatment: str) -> int:
        if treatment not in arm_cache:
            cur.execute(
                "INSERT OR IGNORE INTO treatment_arm(treatment) VALUES(?)", (treatment,))
            cur.execute(
                "SELECT arm_id FROM treatment_arm WHERE treatment=?", (treatment,))
            arm_cache[treatment] = cur.fetchone()[0]
        return arm_cache[treatment]

    def get_or_create_subject(ext_id, proj_id, condition, age, sex,
                              treatment, response) -> int:
        if ext_id not in subject_cache:
            cur.execute(
                """
                INSERT OR IGNORE INTO subject(external_id, project_id, condition, age, sex)
                VALUES(?,?,?,?,?)
                """,
                (ext_id, proj_id, condition,
                 int(age) if age else None,
                 sex or None),
            )
            cur.execute(
                "SELECT subject_id FROM subject WHERE external_id=?", (ext_id,))
            sid = cur.fetchone()[0]
            subject_cache[ext_id] = sid

            arm_id = get_or_create_arm(treatment)
            cur.execute(
                """
                INSERT OR IGNORE INTO subject_treatment(subject_id, arm_id, response)
                VALUES(?,?,?)
                """,
                (sid, arm_id, response or None),
            )
        return subject_cache[ext_id]

    delimiter    = "\t"
    skipped      = 0
    rows_inserted = 0

    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        sample_line = fh.readline()
        fh.seek(0)
        if "," in sample_line and "\t" not in sample_line:
            delimiter = ","

        reader = csv.DictReader(fh, delimiter=delimiter)

        for row in reader:
            proj_id = get_or_create_project(row["project"])
            subj_id = get_or_create_subject(
                row["subject"],
                proj_id,
                row.get("condition", row.get("indication", "")),
                row.get("age", ""),
                row.get("sex", row.get("gender", "")),
                row["treatment"],
                row["response"],
            )

            cur.execute(
                """
                INSERT OR IGNORE INTO sample(
                    external_id, subject_id, sample_type,
                    time_from_treatment_start,
                    b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte
                )
                VALUES(?,?,?,?,?,?,?,?,?)
                """,
                (
                    row["sample"],
                    subj_id,
                    row.get("sample_type", ""),
                    int(row["time_from_treatment_start"])
                    if row.get("time_from_treatment_start") else None,
                    int(row["b_cell"]),
                    int(row["cd8_t_cell"]),
                    int(row["cd4_t_cell"]),
                    int(row["nk_cell"]),
                    int(row["monocyte"]),
                ),
            )

            if cur.rowcount == 0:
                skipped += 1
            else:
                rows_inserted += 1

    conn.commit()
    print(f"Rows inserted: {rows_inserted}")
    print(f"Duplicate/skipped rows: {skipped}")
    return rows_inserted

# Entry Point
def main():
    parser = argparse.ArgumentParser(
        description="Initialise DB and load cell-count CSV.")
    parser.add_argument("--csv", default=DEFAULT_CSV,
                        help="Path to cell-count CSV/TSV file")
    parser.add_argument("--db",  default=DB_PATH,
                        help="Output SQLite database path")
    args = parser.parse_args()

    if not os.path.exists(args.csv):
        print(f"ERROR: CSV file not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    print(f"Creating database: {args.db}")
    conn = sqlite3.connect(args.db)
    conn.executescript(SCHEMA)

    print(f"Loading data from: {args.csv}")
    n = load_csv(args.csv, conn)
    conn.close()
    print(f"Done. {n} samples loaded into {args.db}")

if __name__ == "__main__":
    main()