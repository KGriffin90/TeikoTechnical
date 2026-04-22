"""dashboard.py 
Interactive dashboard for the Miraclib trial. 
To view, run in terminal with make dashboard or python3 dashboard.py."""

import json
import os

import pandas as pd
import dash
import numpy as np
from dash import Input, Output, dcc, html, dash_table
import plotly.graph_objects as go

from analysis import (
    POPULATIONS,
    avg_b_cells_melanoma_males_responders_t0,
    build_response_model,
    compute_frequency_table,
    compute_responder_frequencies,
    compute_subset_summaries,
    make_longitudinal_figure,
    run_statistical_tests,
)

COLORS = {
    "bg":           "#f4f5f7",
    "panel":        "#ffffff",
    "border":       "#dde1e7",
    "text":         "#1a202c",
    "muted":        "#718096",
    "accent":       "#2c5282",
    "responder":    "#2c5282",
    "nonresponder": "#E07B39",
    "grid":         "#edf2f7",
    "sig":          "#276749",
}

POP_COLORS = ["#2c5282", "#2d6a4f", "#b7791f", "#6b46c1", "#c53030"]

LAYOUT_BASE = dict(
    paper_bgcolor="white",
    plot_bgcolor="white",
    font=dict(family="IBM Plex Sans, sans-serif", color=COLORS["text"], size=12),
    margin=dict(l=50, r=20, t=40, b=40),
)

def axis_style(**kwargs):
    return dict(
        showgrid=True,
        gridcolor=COLORS["grid"],
        gridwidth=1,
        zeroline=False,
        linecolor=COLORS["border"],
        tickfont=dict(size=11, color=COLORS["muted"]),
        **kwargs,
    )

def card(children, style=None):
    base = {
        "background":    COLORS["panel"],
        "border":        f"1px solid {COLORS['border']}",
        "borderRadius":  "6px",
        "padding":       "20px 24px",
        "marginBottom":  "20px",
    }
    if style:
        base.update(style)
    return html.Div(children, style=base)


def section_label(text):
    return html.P(text, style={
        "fontFamily":    "'IBM Plex Mono', monospace",
        "fontSize":      "11px",
        "color":         COLORS["accent"],
        "letterSpacing": "0.12em",
        "textTransform": "uppercase",
        "marginBottom":  "4px",
    })


def section_heading(text):
    return html.H2(text, style={
        "fontSize":     "22px",
        "fontWeight":   "600",
        "color":        COLORS["text"],
        "marginBottom": "4px",
    })


def section_sub(text):
    return html.P(text, style={
        "fontSize":     "13px",
        "color":        COLORS["muted"],
        "marginBottom": "24px",
    })


def stat_card(label, value, desc):
    return html.Div([
        html.P(label, style={
            "fontFamily":    "'IBM Plex Mono', monospace",
            "fontSize":      "10px",
            "color":         COLORS["muted"],
            "textTransform": "uppercase",
            "letterSpacing": "0.08em",
            "marginBottom":  "6px",
        }),
        html.P(str(value), style={
            "fontFamily":   "'IBM Plex Mono', monospace",
            "fontSize":     "28px",
            "fontWeight":   "600",
            "color":        COLORS["accent"],
            "marginBottom": "4px",
        }),
        html.P(desc, style={
            "fontSize": "12px",
            "color":    COLORS["muted"],
        }),
    ], style={
        "background":   COLORS["panel"],
        "border":       f"1px solid {COLORS['border']}",
        "borderRadius": "6px",
        "padding":      "16px 20px",
    })

def make_dotplot_figure(responder_df):
    """
    Dot plot with median bar and IQR band.
    One panel per cell population: responders (blue) vs non-responders (orange).
    """
    fig = go.Figure()

    group_map  = {"yes": ("Responder",     COLORS["responder"]),
                  "no":  ("Non-Responder", COLORS["nonresponder"])}
    x_positions = {pop: i for i, pop in enumerate(POPULATIONS)}
    n_groups    = 2
    offsets     = {"yes": -0.18, "no": 0.18}

    first = {"yes": True, "no": True}

    for resp, (label, color) in group_map.items():
        for pop in POPULATIONS:
            xi   = x_positions[pop] + offsets[resp]
            vals = responder_df[
                (responder_df["population"] == pop) &
                (responder_df["response"]   == resp)
            ]["percentage"].values

            if len(vals) == 0:
                continue

            median = float(np.median(vals))
            q1     = float(np.percentile(vals, 25))
            q3     = float(np.percentile(vals, 75))

            # IQR band
            fig.add_trace(go.Scatter(
                x=[xi, xi], y=[q1, q3],
                mode="lines",
                line=dict(color=color, width=6),
                opacity=0.35,
                showlegend=False,
                hoverinfo="skip",
            ))

            # Median tick
            hw = 0.10
            fig.add_trace(go.Scatter(
                x=[xi - hw, xi + hw], y=[median, median],
                mode="lines",
                line=dict(color=color, width=2.5),
                showlegend=False,
                hoverinfo="skip",
            ))

            # Individual data points
            rng = np.random.default_rng(seed=hash(pop + resp) & 0xFFFF)
            jitter = rng.uniform(-0.07, 0.07, size=len(vals))

            fig.add_trace(go.Scatter(
                x=xi + jitter,
                y=vals,
                mode="markers",
                name=label,
                legendgroup=label,
                showlegend=first[resp],
                marker=dict(
                    color=color,
                    size=7,
                    opacity=0.75,
                    line=dict(color="white", width=0.8),
                ),
                hovertemplate=(
                    f"<b>{label}</b><br>"
                    f"{pop.replace('_', ' ').title()}<br>"
                    "Frequency: %{y:.2f}%<extra></extra>"
                ),
            ))
            first[resp] = False

    pop_labels = [p.replace("_", " ").title() for p in POPULATIONS]

    fig.update_layout(
        **LAYOUT_BASE,
        height=440,
        xaxis=axis_style(
            tickmode="array",
            tickvals=list(range(len(POPULATIONS))),
            ticktext=pop_labels,
            title="Cell Population",
        ),
        yaxis=axis_style(title="Relative Frequency (%)"),
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
            font=dict(size=12, color=COLORS["text"]),
        ),
    )
    return fig

def make_stats_table(stats_df):
    rows = []
    for _, r in stats_df.iterrows():
        rows.append({
            "population":       r["population"],
            "n_responders":     int(r["n_responders"]),
            "n_non_responders": int(r["n_non_responders"]),
            "p_value":          r["p_value"],
            "p_value_adj":      r["p_value_adj"],
            "effect_size":      r["effect_size"],
            "significant":      "✓  Yes" if r["significant"] else "✗  No",
        })

    return dash_table.DataTable(
        columns=[
            {"name": "Population",        "id": "population"},
            {"name": "n Responders",      "id": "n_responders",     "type": "numeric"},
            {"name": "n Non-Responders",  "id": "n_non_responders", "type": "numeric"},
            {"name": "p-value",           "id": "p_value",          "type": "numeric",
             "format": dash_table.Format.Format(precision=6, scheme=dash_table.Format.Scheme.fixed)},
            {"name": "Adj p-value (FDR)", "id": "p_value_adj",      "type": "numeric",
             "format": dash_table.Format.Format(precision=6, scheme=dash_table.Format.Scheme.fixed)},
            {"name": "Effect Size (r)",   "id": "effect_size",      "type": "numeric",
             "format": dash_table.Format.Format(precision=4, scheme=dash_table.Format.Scheme.fixed)},
            {"name": "Significant?",      "id": "significant"},
        ],
        data=rows,
        style_table={"overflowX": "auto"},
        style_header={
            "backgroundColor": "#f8f9fb",
            "fontFamily":      "'IBM Plex Mono', monospace",
            "fontSize":        "11px",
            "color":           COLORS["muted"],
            "textTransform":   "uppercase",
            "letterSpacing":   "0.06em",
            "border":          f"1px solid {COLORS['border']}",
            "padding":         "10px 12px",
        },
        style_cell={
            "fontFamily":      "'IBM Plex Mono', monospace",
            "fontSize":        "12px",
            "color":           COLORS["text"],
            "border":          f"1px solid {COLORS['border']}",
            "padding":         "8px 14px",
            "backgroundColor": "white",
        },
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfc"},
            {"if": {"filter_query": '{significant} = "✓  Yes"', "column_id": "significant"},
             "color": COLORS["sig"], "fontWeight": "600"},
            {"if": {"filter_query": '{significant} = "✗  No"', "column_id": "significant"},
             "color": COLORS["muted"]},
            {"if": {"filter_query": "{p_value_adj} < 0.05", "column_id": "p_value_adj"},
             "color": COLORS["sig"], "fontWeight": "600"},
            {"if": {"column_id": "p_value"}, "color": COLORS["muted"]},
        ],
    )

def make_interpretation(stats_df):
    sig_pops = stats_df[stats_df["significant"]]["population"].tolist()
    n_r      = int(stats_df.iloc[0]["n_responders"])
    n_nr     = int(stats_df.iloc[0]["n_non_responders"])

    if sig_pops:
        return (
            f"The following populations show a statistically significant difference "
            f"(p < 0.05) between responders and non-responders: "
            f"{', '.join(sig_pops)}. These may be relevant biomarkers for predicting "
            f"Miraclib response in melanoma patients."
        )
    return (
        f"No cell population reached statistical significance (p < 0.05). "
        f"All five populations showed comparable relative frequencies between "
        f"responders (n = {n_r}) and non-responders (n = {n_nr}). "
        f"Consider longitudinal analysis across treatment timepoints or additional "
        f"biomarker modalities to identify predictors of Miraclib response."
    )

def make_freq_table_component(freq_table):
    pop_color_map = dict(zip(POPULATIONS, POP_COLORS))

    return html.Div([
        dash_table.DataTable(
            id="freq-datatable",
            columns=[
            {
                "name": "Sample",
                "id": "sample",
                "filter_options": {"placeholder_text": "e.g. sample001"},
            },
            {
                "name": "Total Count",
                "id": "total_count",
                "type": "numeric",
                "format": dash_table.Format.Format(group=","),
                "filter_options": {"placeholder_text": "e.g. >1000"},
            },
            {
                "name": "Population",
                "id": "population",
                "filter_options": {"placeholder_text": "e.g. b_cell"},
            },
            {
                "name": "Count",
                "id": "count",
                "type": "numeric",
                "format": dash_table.Format.Format(group=","),
                "filter_options": {"placeholder_text": "e.g. >200"},
            },
            {
                "name": "% Frequency",
                "id": "percentage",
                "type": "numeric",
                "format": dash_table.Format.Format(precision=2, scheme=dash_table.Format.Scheme.fixed),
                "filter_options": {"placeholder_text": "e.g. >20"},
            },
],
            data=freq_table.to_dict("records"),
            page_size=15,
            page_action="native",
            sort_action="native",
            filter_action="native",
            style_table={"overflowX": "auto"},
            style_header={
                "backgroundColor": "#f8f9fb",
                "fontFamily":      "'IBM Plex Mono', monospace",
                "fontSize":        "11px",
                "color":           COLORS["muted"],
                "textTransform":   "uppercase",
                "letterSpacing":   "0.06em",
                "border":          f"1px solid {COLORS['border']}",
                "borderBottom":    f"2px solid {COLORS['border']}",
                "padding":         "10px 12px",
            },
            style_cell={
                "fontFamily":      "'IBM Plex Mono', monospace",
                "fontSize":        "12px",
                "color":           COLORS["text"],
                "border":          f"1px solid {COLORS['border']}",
                "padding":         "8px 12px",
                "backgroundColor": "white",
                "whiteSpace":      "nowrap",
            },
            style_data_conditional=[
                {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfc"},
                *[
                    {"if": {"filter_query": f'{{population}} = "{pop}"',
                            "column_id": "population"},
                     "color": color, "fontWeight": "500"}
                    for pop, color in pop_color_map.items()
                ],
                {"if": {"column_id": "sample"}, "color": COLORS["accent"]},
            ],
        ),
    ])

def make_bar(df, x_col, y_col, title, color):
    layout = {**LAYOUT_BASE, "margin": dict(l=40, r=20, t=40, b=40)}
    fig    = go.Figure(go.Bar(
        x=df[x_col], y=df[y_col],
        marker_color=color,
        marker_line_width=0,
        text=df[y_col], textposition="outside",
        textfont=dict(family="IBM Plex Mono", size=11, color=COLORS["text"]),
    ))
    fig.update_layout(
        **layout,
        title=dict(text=title, font=dict(size=13, color=COLORS["muted"]), x=0),
        height=260,
        xaxis=axis_style(),
        yaxis=axis_style(title="Count"),
        bargap=0.4,
    )
    return fig

def make_subset_figures(subset):
    d    = subset
    figs = []
    _layout = {**LAYOUT_BASE, "margin": dict(l=40, r=20, t=40, b=40)}

    figs.append(make_bar(
        d["samples_per_project"], "project", "n_samples",
        "4a  -  Samples per Project", COLORS["accent"]
    ))

    resp = d["response_counts"].copy()
    resp["label"] = resp["response"].map({"yes": "Responder", "no": "Non-Responder"})
    fig_b = go.Figure(go.Bar(
        x=resp["label"], y=resp["n_subjects"],
        marker_color=[COLORS["responder"], COLORS["nonresponder"]],
        marker_line_width=0,
        text=resp["n_subjects"], textposition="outside",
        textfont=dict(family="IBM Plex Mono", size=11, color=COLORS["text"]),
    ))
    fig_b.update_layout(
        **_layout,
        title=dict(text="4b  -  Responders vs Non-Responders",
                   font=dict(size=13, color=COLORS["muted"]), x=0),
        height=260,
        xaxis=axis_style(), yaxis=axis_style(title="Subjects"),
        bargap=0.4,
    )
    figs.append(fig_b)

    sex = d["sex_counts"].copy()
    sex["label"] = sex["sex"].map({"M": "Male", "F": "Female"})
    fig_c = go.Figure(go.Bar(
        x=sex["label"], y=sex["n_subjects"],
        marker_color=[COLORS["accent"], COLORS["muted"]],
        marker_line_width=0,
        text=sex["n_subjects"], textposition="outside",
        textfont=dict(family="IBM Plex Mono", size=11, color=COLORS["text"]),
    ))
    fig_c.update_layout(
        **_layout,
        title=dict(text="4c  -  Sex Distribution",
                   font=dict(size=13, color=COLORS["muted"]), x=0),
        height=260,
        xaxis=axis_style(), yaxis=axis_style(title="Subjects"),
        bargap=0.4,
    )
    figs.append(fig_c)
    return figs

# App Initialization
app = dash.Dash(
    __name__,
    title="Loblaw Bio : Miraclib Trial",
    suppress_callback_exceptions=True,
    external_stylesheets=[
        "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600"
        "&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap"
    ],
)

app.index_string = """<!DOCTYPE html>
<html>
<head>
  {%metas%}
  <title>{%title%}</title>
  {%favicon%}
  {%css%}
  <style>
    * { box-sizing: border-box; }
    body { margin: 0; background: #f4f5f7; font-family: 'IBM Plex Sans', sans-serif; }
    .tab--selected {
      border-top: 2px solid #2c5282 !important;
      color: #2c5282 !important;
      font-weight: 600 !important;
    }
    .dash-table-container .dash-filter input {
      font-family: 'IBM Plex Mono', monospace;
      font-size: 11px;
    }
  </style>
</head>
<body>{%app_entry%}<footer>{%config%}{%scripts%}{%renderer%}</footer></body>
</html>"""


app.layout = html.Div([

    html.Div([
        html.Span("LOBLAW BIO", style={
            "fontFamily":    "'IBM Plex Mono', monospace",
            "fontWeight":    "600",
            "fontSize":      "14px",
            "color":         COLORS["accent"],
            "letterSpacing": "0.08em",
        }),
        html.Span(" | MIRACLIB TRIAL", style={
            "fontFamily":    "'IBM Plex Mono', monospace",
            "fontSize":      "13px",
            "color":         COLORS["muted"],
            "letterSpacing": "0.06em",
        }),
    ], style={
        "background":    "white",
        "borderBottom":  f"1px solid {COLORS['border']}",
        "padding":       "0 32px",
        "height":        "52px",
        "display":       "flex",
        "alignItems":    "center",
        "gap":           "4px",
        "boxShadow":     "0 1px 3px rgba(0,0,0,0.06)",
        "position":      "sticky",
        "top":           "0",
        "zIndex":        "100",
    }),

    dcc.Store(id="analysis-store"),

    html.Div([
        dcc.Tabs(
            id="tabs",
            value="tab-overview",
            children=[
                dcc.Tab(label="Overview",             value="tab-overview"),
                dcc.Tab(label="Cohort",               value="tab-frequency"),
                dcc.Tab(label="Statistical Analysis", value="tab-stats"),
                dcc.Tab(label="Explorer",             value="tab-subset"),
                dcc.Tab(label="Prediction",           value="tab-prediction"),
            ],
            style={"fontFamily": "'IBM Plex Mono', monospace", "fontSize": "12px"},
            colors={"border": COLORS["border"],
                    "primary": COLORS["accent"],
                    "background": "white"},
        ),
        html.Div(id="tab-content",
                 style={"padding": "28px 32px",
                        "maxWidth": "1280px",
                        "margin": "0 auto"}),
    ]),

], style={"minHeight": "100vh", "background": COLORS["bg"]})


# Store Population Callback
@app.callback(
    Output("analysis-store", "data"),
    Input("tabs", "value"),
)
def populate_store(_tab):
    """
    Load all analysis data once and cache it in dcc.Store.
    The underlying functions use lru_cache keyed on db mtime, so this
    re-fetches automatically if the database file is updated.
    """
    freq_table   = compute_frequency_table()
    responder_df = compute_responder_frequencies()
    stats_df     = run_statistical_tests(responder_df)
    subset       = compute_subset_summaries()
    model        = build_response_model()

    return {
        "freq_table":   freq_table.to_json(orient="records"),
        "responder_df": responder_df.to_json(orient="records"),
        "stats_df":     stats_df.to_json(orient="records"),
        "subset": {
            "baseline_n":          len(subset["baseline_df"]),
            "samples_per_project": subset["samples_per_project"].to_json(orient="records"),
            "response_counts":     subset["response_counts"].to_json(orient="records"),
            "sex_counts":          subset["sex_counts"].to_json(orient="records"),
        },
        "model": model,
    }

# Tab Routing Callback
@app.callback(
    Output("tab-content", "children"),
    Input("tabs",          "value"),
    Input("analysis-store", "data"),
)
def render_tab(tab, store):
    def from_store(json_str):
        return pd.DataFrame(json.loads(json_str))

    if store is None:
        return html.P("Loading…", style={"color": COLORS["muted"], "padding": "40px"})

    freq_table   = from_store(store["freq_table"])
    responder_df = from_store(store["responder_df"])
    stats_df     = from_store(store["stats_df"])

    s            = store["subset"]
    n_baseline   = s["baseline_n"]
    subset_data  = {
        "samples_per_project": from_store(s["samples_per_project"]),
        "response_counts":     from_store(s["response_counts"]),
        "sex_counts":          from_store(s["sex_counts"]),
    }

    n_samples  = freq_table["sample"].nunique()
    n_records  = len(freq_table)

    if tab == "tab-overview":
        sig_pops   = stats_df[stats_df["significant"]]["population"].tolist() \
                     if "significant" in stats_df.columns else []
        sig_text   = (", ".join(sig_pops) if sig_pops
                      else "None : no population reached α = 0.05")
        sig_color  = COLORS["sig"] if sig_pops else COLORS["nonresponder"]

        finding_items = [
            ("Study",           "Miraclib (melanoma PBMC) | Mann Whitney U | FDR Corrected"),
            ("Populations",     "B cell, CD8 T, CD4 T, NK, Monocyte"),
            ("Significant",     sig_text),
            ("Interpretation",  "Consider longitudinal or multi modal biomarker analysis"
                                if not sig_pops else
                                "Highlighted populations are candidate predictive biomarkers"),
        ]

        return html.Div([
            section_label("Overview"),
            section_heading("Miraclib Clinical Trial"),
            section_sub("High level summary: use the tabs above to explore in detail"),

            html.Div([
                stat_card("Total Samples",    f"{n_samples:,}",  "across all projects"),
                stat_card("Total Records",    f"{n_records:,}",  "sample × population rows"),
                stat_card("Baseline Subset",  n_baseline,        "melanoma | miraclib | PBMC | t=0"),
                stat_card("Populations",      5,                 "B, CD8 T, CD4 T, NK, Monocyte"),
                stat_card("Statistical Test", "M-W U",           "non-parametric | two-sided"),
            ], style={"display": "grid",
                      "gridTemplateColumns": "repeat(auto-fit, minmax(180px, 1fr))",
                      "gap": "14px", "marginBottom": "28px"}),

            section_label("Key Findings"),
            card([
                html.Table([
                    html.Tbody([
                        html.Tr([
                            html.Td(k, style={
                                "fontFamily":    "'IBM Plex Mono', monospace",
                                "fontSize":      "11px",
                                "color":         COLORS["muted"],
                                "textTransform": "uppercase",
                                "letterSpacing": "0.08em",
                                "paddingRight":  "24px",
                                "paddingBottom": "10px",
                                "whiteSpace":    "nowrap",
                                "verticalAlign": "top",
                            }),
                            html.Td(v, style={
                                "fontSize":      "13px",
                                "color":         sig_color if k == "Significant" else COLORS["text"],
                                "fontWeight":    "600" if k == "Significant" else "400",
                                "paddingBottom": "10px",
                            }),
                        ])
                        for k, v in finding_items
                    ])
                ], style={"borderCollapse": "collapse", "width": "100%"}),
            ]),

            section_label("Data Model"),
            html.P(
                "The database schema diagram is available in docs/schema.png : "
                "run python3 export_schema.py to regenerate it.",
                style={"fontSize": "13px", "color": COLORS["muted"]},
            ),
        ])

    # Cohort
    elif tab == "tab-frequency":
        return html.Div([
            section_label("Part 2 : Initial Analysis"),
            section_heading("Cell Population Frequencies"),
            section_sub("Relative frequency of each immune cell population per sample - "
                        "Click column headers to sort - Use filter row to search"),
            card(make_freq_table_component(freq_table)),
        ])

    # Statistical Analysis
    elif tab == "tab-stats":
        interp     = make_interpretation(stats_df)
        sig_pops   = stats_df[stats_df["significant"]]["population"].tolist() \
                     if "significant" in stats_df.columns else []
        interp_color = COLORS["sig"] if sig_pops else COLORS["muted"]
        summary_text = (
            f"Significant: {', '.join(sig_pops)}" if sig_pops
            else "No significant differences at α = 0.05"
        )

        return html.Div([
            section_label("Part 3 : Statistical Analysis"),
            section_heading("Responders vs Non-Responders"),
            section_sub("Melanoma patients on Miraclib (PBMC samples) - "
                        "Mann-Whitney U test - α = 0.05"),

            card([
                html.P(summary_text, style={
                    "fontSize":   "13px",
                    "fontWeight": "500",
                    "color":      interp_color,
                })
            ], style={"borderLeft": f"3px solid {interp_color}"}),

            card([
                html.P("Interpretation", style={
                    "fontFamily":    "'IBM Plex Mono', monospace",
                    "fontSize":      "11px",
                    "color":         COLORS["accent"],
                    "letterSpacing": "0.12em",
                    "textTransform": "uppercase",
                    "marginBottom":  "8px",
                }),
                html.P(interp, style={
                    "fontSize":   "14px",
                    "lineHeight": "1.75",
                    "color":      COLORS["text"],
                }),
            ], style={"borderLeft": f"3px solid {COLORS['accent']}"}),

            section_label("Distribution"),
            html.P(
                "Each point is one sample. Thick bar = median; thin band = IQR. "
                "Blue = Responder - Orange = Non-Responder.",
                style={"fontSize": "12px", "color": COLORS["muted"],
                       "marginBottom": "8px"},
            ),
            card(dcc.Graph(
                figure=make_dotplot_figure(responder_df),
                config={"displayModeBar": False},
            )),

            section_label("Longitudinal Dynamics"),
            html.H3("Immune Trajectories Over Time (t = 0, 7, 14)", style={
                "fontSize":     "14px",
                "fontWeight":   "500",
                "color":        COLORS["muted"],
                "marginBottom": "12px",
            }),
            card(dcc.Graph(
                figure=make_longitudinal_figure(),
                config={"displayModeBar": False},
            )),

            section_label("Statistical Results"),
            html.H3("Mann-Whitney U : Population Comparison", style={
                "fontSize":     "14px",
                "fontWeight":   "500",
                "color":        COLORS["muted"],
                "marginBottom": "12px",
            }),
            card(make_stats_table(stats_df)),
        ])

    # Explorer
    elif tab == "tab-subset":
        rc      = subset_data["response_counts"]
        n_resp  = int((rc.query("response == 'yes'")["n_subjects"].values or [0])[0])
        n_nonr  = int((rc.query("response == 'no'")["n_subjects"].values  or [0])[0])
        n_proj  = len(subset_data["samples_per_project"])
        sfigs   = make_subset_figures(subset_data)

        return html.Div([
            section_label("Part 4 : Data Subset Analysis"),
            section_heading("Baseline Melanoma - Miraclib - PBMC"),
            section_sub("time_from_treatment_start = 0 | condition = melanoma | "
                        "treatment = miraclib | sample_type = PBMC"),

            html.Div([
                stat_card("Baseline Samples", n_baseline,  "matching all filters"),
                stat_card("Projects",         n_proj,      "contributing samples"),
                stat_card("Responders",       n_resp,      "unique subjects"),
                stat_card("Non-Responders",   n_nonr,      "unique subjects"),
                stat_card("Avg B-cells (♂ resp.)",
                          f"{avg_b_cells_melanoma_males_responders_t0():.2f}",
                          "baseline melanoma males"),
            ], style={"display": "grid",
                      "gridTemplateColumns": "repeat(auto-fit, minmax(180px, 1fr))",
                      "gap": "14px", "marginBottom": "28px"}),

            html.Div([
                card(dcc.Graph(figure=sfigs[0], config={"displayModeBar": False})),
                card(dcc.Graph(figure=sfigs[1], config={"displayModeBar": False})),
                card(dcc.Graph(figure=sfigs[2], config={"displayModeBar": False})),
            ], style={"display": "grid",
                      "gridTemplateColumns": "repeat(auto-fit, minmax(280px, 1fr))",
                      "gap": "16px"}),
        ])

    # Prediction
    elif tab == "tab-prediction":
        model = store.get("model")

        if model is None:
            return html.Div([
                section_label("Part 5 : Predictive Modelling"),
                section_heading("Response Prediction Model"),
                card(html.P(
                    "Not enough data to fit a model. "
                    "Need at least 2 samples per class at baseline.",
                    style={"color": COLORS["muted"], "fontSize": "14px"},
                )),
            ])

        auc_mean    = model["auc_mean"]
        auc_std     = model["auc_std"]
        cv_aucs     = model["cv_aucs"]
        n_splits    = model["n_splits"]
        insample_auc = model["insample_auc"]
        fpr         = model["fpr"]
        tpr         = model["tpr"]
        n_samples   = model["n_samples"]
        n_resp      = model["n_responders"]
        n_nonr      = model["n_non_responders"]
        coef_df     = pd.DataFrame(model["feature_importance"])

        # ROC curve
        roc_fig = go.Figure()
        roc_fig.add_trace(go.Scatter(
            x=[0, 1], y=[0, 1],
            mode="lines",
            name="Random chance",
            line=dict(color=COLORS["muted"], dash="dash", width=1),
        ))
        roc_fig.add_trace(go.Scatter(
            x=fpr, y=tpr,
            mode="lines",
            name=f"In-sample ROC (AUC = {insample_auc:.3f})",
            line=dict(color=COLORS["nonresponder"], width=2),
            fill="tozeroy",
            fillcolor="rgba(224,123,57,0.08)",
        ))
        roc_fig.update_layout(
            **LAYOUT_BASE,
            height=380,
            xaxis=axis_style(title="False Positive Rate", range=[-0.02, 1.02]),
            yaxis=axis_style(title="True Positive Rate",  range=[-0.02, 1.02]),
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
                font=dict(size=11),
            ),
        )

        # AUC strip
        cv_fig = go.Figure()
        fold_labels = [f"Fold {i+1}" for i in range(len(cv_aucs))]
        cv_colors   = [COLORS["sig"] if v >= 0.7
                       else COLORS["accent"] if v >= 0.5
                       else COLORS["nonresponder"]
                       for v in cv_aucs]
        cv_fig.add_trace(go.Bar(
            x=fold_labels,
            y=cv_aucs,
            marker_color=cv_colors,
            marker_line_width=0,
            text=[f"{v:.3f}" for v in cv_aucs],
            textposition="outside",
            textfont=dict(family="IBM Plex Mono", size=11, color=COLORS["text"]),
        ))
        cv_fig.add_hline(
            y=0.5,
            line_dash="dash",
            line_color=COLORS["muted"],
            annotation_text="random",
            annotation_font_size=10,
        )
        cv_fig.add_hline(
            y=auc_mean,
            line_dash="dot",
            line_color=COLORS["accent"],
            annotation_text=f"mean = {auc_mean:.3f}",
            annotation_font_size=10,
        )
        cv_fig.update_layout(
            **LAYOUT_BASE,
            height=280,
            xaxis=axis_style(),
            yaxis=axis_style(title="AUC", range=[0, 1.1]),
            bargap=0.35,
        )

        # Feature importance
        coef_colors = [
            COLORS["responder"] if c >= 0 else COLORS["nonresponder"]
            for c in coef_df["coefficient"]
        ]
        coef_fig = go.Figure(go.Bar(
            x=coef_df["population"].str.replace("_", " ").str.title(),
            y=coef_df["coefficient"],
            marker_color=coef_colors,
            marker_line_width=0,
            text=[f"{v:+.3f}" for v in coef_df["coefficient"]],
            textposition="outside",
            textfont=dict(family="IBM Plex Mono", size=11, color=COLORS["text"]),
        ))
        coef_fig.add_hline(y=0, line_color=COLORS["border"], line_width=1)
        coef_fig.update_layout(
            **LAYOUT_BASE,
            height=300,
            xaxis=axis_style(title="Cell Population"),
            yaxis=axis_style(title="Log-odds coefficient"),
            bargap=0.35,
        )

        return html.Div([
            section_label("Part 5 : Predictive Modelling"),
            section_heading("Response Prediction Model"),
            section_sub(
                f"Logistic regression | baseline PBMC frequencies : "
                f"{n_samples} samples ({n_resp} responders, {n_nonr} non-responders)"
            ),

            # Stat cards
            html.Div([
                stat_card("CV AUC",       f"{auc_mean:.3f}",     f"±{auc_std:.3f}  -  {n_splits}-fold"),
                stat_card("In-sample AUC", f"{insample_auc:.3f}", "optimistic : see caveat below"),
                stat_card("Features",     len(POPULATIONS),      "baseline immune cell frequencies"),
                stat_card("Training n",   n_samples,             f"{n_resp} R  -  {n_nonr} NR"),
            ], style={"display": "grid",
                      "gridTemplateColumns": "repeat(auto-fit, minmax(180px, 1fr))",
                      "gap": "14px", "marginBottom": "28px"}),

            # AUC caveat banner
            card([
                html.P("⚠  Interpretation caveat", style={
                    "fontFamily":    "'IBM Plex Mono', monospace",
                    "fontSize":      "11px",
                    "color":         COLORS["nonresponder"],
                    "textTransform": "uppercase",
                    "letterSpacing": "0.1em",
                    "marginBottom":  "6px",
                }),
                html.P(
                    f"The ROC curve below is computed on the same data the model was "
                    f"trained on (in sample AUC = {insample_auc:.3f}), which is "
                    f"optimistic by construction. "
                    f"The {n_splits}-fold cross validated AUC "
                    f"({auc_mean:.3f} ± {auc_std:.3f}) is the honest estimate of "
                    f"generalisation performance. With n = {n_samples} total samples "
                    f"these estimates carry wide confidence intervals. Treat them as "
                    f"exploratory, not confirmatory.",
                    style={"fontSize": "13px", "lineHeight": "1.7",
                           "color": COLORS["text"]},
                ),
            ], style={"borderLeft": f"3px solid {COLORS['nonresponder']}"}),

            # Cross validated AUC per fold
            section_label("Cross-Validated AUC"),
            html.P(
                f"{n_splits}-fold stratified - this is the number to trust",
                style={"fontSize": "12px", "color": COLORS["muted"],
                       "marginBottom": "8px"},
            ),
            card(dcc.Graph(figure=cv_fig, config={"displayModeBar": False})),

            # In-sample ROC
            section_label("ROC Curve  (in-sample : optimistic)"),
            card(dcc.Graph(figure=roc_fig, config={"displayModeBar": False})),

            # Feature importance
            section_label("Feature Importance"),
            html.P(
                "Logistic regression log odds coefficients: "
                "blue = associated with response | orange = associated with non-response",
                style={"fontSize": "12px", "color": COLORS["muted"],
                       "marginBottom": "8px"},
            ),
            card(dcc.Graph(figure=coef_fig, config={"displayModeBar": False})),
        ])


# Run
if __name__ == "__main__":
    debug = os.environ.get("DASH_DEBUG", "false").lower() == "true"
    print("Starting Loblaw Bio Dashboard → http://127.0.0.1:8050")
    app.run(debug=debug)
