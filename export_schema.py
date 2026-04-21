"""export_schema.py
Renders the database schema as a PNG and writes it to docs/schema.png."""

import os
import pathlib

import plotly.graph_objects as go

COLORS = {
    "border":  "#dde1e7",
    "muted":   "#718096",
    "accent":  "#2c5282",
    "text":    "#1a202c",
    "grid":    "#edf2f7",
}

def build_schema_figure() -> go.Figure:
    """Return a Plotly figure of the SQLite schema."""
    tables = [
        {"name": "project",           "x": 0.02,
         "cols": ["project_id  PK", "name"]},
        {"name": "subject",           "x": 0.22,
         "cols": ["subject_id  PK", "external_id", "project_id  FK",
                  "condition", "age", "sex"]},
        {"name": "subject_treatment", "x": 0.50,
         "cols": ["subject_id  FK", "arm_id  FK", "response"]},
        {"name": "treatment_arm",     "x": 0.76,
         "cols": ["arm_id  PK", "treatment"]},
        {"name": "sample",            "x": 0.22, "y_offset": -0.55,
         "cols": ["sample_id  PK", "external_id", "subject_id  FK",
                  "sample_type", "time_from_treatment_start",
                  "b_cell … monocyte"]},
    ]

    row_h    = 0.07
    header_h = 0.10
    box_w    = 0.18
    shapes, annotations = [], []
    box_centers: dict[str, tuple[float, float]] = {}

    for t in tables:
        x0    = t["x"]
        y_off = t.get("y_offset", 0)
        n     = len(t["cols"])
        total_h = header_h + n * row_h + 0.02
        y1    = 0.95 + y_off
        y0    = y1 - total_h
        cx    = x0 + box_w / 2
        box_centers[t["name"]] = (cx, y1)

        shapes.append(dict(
            type="rect", x0=x0, x1=x0 + box_w, y0=y0, y1=y1,
            fillcolor="#f8f9fb",
            line=dict(color=COLORS["border"], width=1.5),
            xref="paper", yref="paper",
        ))
        shapes.append(dict(
            type="rect", x0=x0, x1=x0 + box_w, y0=y1 - header_h, y1=y1,
            fillcolor="#edf2f7",
            line=dict(color=COLORS["border"], width=0),
            xref="paper", yref="paper",
        ))
        annotations.append(dict(
            x=cx, y=y1 - header_h / 2, text=f"<b>{t['name']}</b>",
            xref="paper", yref="paper", showarrow=False,
            font=dict(size=11, color=COLORS["accent"], family="IBM Plex Mono"),
            xanchor="center",
        ))
        for i, col in enumerate(t["cols"]):
            color = (COLORS["accent"] if "PK" in col
                     else COLORS["muted"] if "FK" in col
                     else COLORS["text"])
            cy = y1 - header_h - (i + 0.7) * row_h
            annotations.append(dict(
                x=x0 + 0.01, y=cy, text=col,
                xref="paper", yref="paper", showarrow=False,
                font=dict(size=9, color=color, family="IBM Plex Mono"),
                xanchor="left",
            ))

    rels = [
        ("project", "subject"),
        ("subject", "subject_treatment"),
        ("subject_treatment", "treatment_arm"),
        ("subject", "sample"),
    ]
    for src, dst in rels:
        sx, sy = box_centers[src]
        dx, dy = box_centers[dst]
        ex = sx + box_w / 2 if dx > sx else sx - box_w / 2
        shapes.append(dict(
            type="line",
            x0=ex,             y0=sy - 0.05,
            x1=dx - box_w / 2 if dx > sx else dx + box_w / 2,
            y1=dy - 0.05,
            line=dict(color=COLORS["border"], width=1.5, dash="dot"),
            xref="paper", yref="paper",
        ))

    fig = go.Figure()
    fig.update_layout(
        shapes=shapes,
        annotations=annotations,
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[0, 1]),
        height=400,
        margin=dict(l=10, r=10, t=10, b=10),
        font=dict(family="IBM Plex Sans, sans-serif", size=12),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return fig

def main(out_path: str = "docs/schema.png", scale: int = 2) -> None:
    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig = build_schema_figure()
    fig.write_image(out_path, scale=scale)
    print(f"Schema diagram saved → {out_path}")


if __name__ == "__main__":
    main()