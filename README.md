Loblaw Bio: Miraclib Clinical Trial Analysis

This project analyzes immune cell populations in melanoma patients treated with Miraclib. It provides a pipeline that loads raw patient data, performs statistical comparisons between responders and non-responders, and exposes the results through an interactive dashboard.

The analysis focuses on computing immune cell frequencies, comparing groups using non-parametric statistical testing, and evaluating baseline cohort behavior at time point zero. Statistical significance is assessed using the Mann Whitney U test, with multiple comparison correction applied using the Benjamini Hochberg method. Effect sizes are also computed using rank biserial correlation to give additional context beyond p-values.

Data is taken from a single input file (cell-count.csv). The pipeline automatically detects the format, processes the data, and builds a local database (cell_counts.db) to support analysis and reproducibility.


Running the Project:

1) Install dependencies:
    `make setup`

2) Run the full pipeline (loads data, runs analysis, renders schema):
    `make pipeline`

3) Launch the dashboard:
    `make dashboard`

4) Open in browser at http://127.0.0.1:8050

Other available targets:
- `make test` — run the test suite  
- `make schema` — regenerate `docs/schema.png`  
- `make clean` — remove the generated database and cached files


Dashboard

http://127.0.0.1:8050 (start locally with make dashboard)

The dashboard has five tabs:

1) Overview: key findings summary and cohort statistics
2) Cohort: per sample cell population frequency table with sorting and filtering
3) Statistical Analysis: dot plot of responder vs non-responder distributions,
   longitudinal immune trajectories across t=0/7/14, and the full Mann Whitney
   U results table with FDR adjusted p-values
4) Explorer: subset breakdown by project, response, and sex
5) Prediction: logistic regression model with cross validated AUC, ROC
   curve, and feature importance coefficients


Database schema

The raw CSV is normalised into five tables:

    project (project_id, name)
    subject (subject_id, external_id, project_id, condition, age, sex)
    treatment_arm (arm_id, treatment)
    subject_treatment (subject_id, arm_id, response)
    sample (sample_id, external_id, subject_id, sample_type,
            time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell,
            nk_cell, monocyte)

The schema diagram is located at docs/schema.png (regenerates with: make schema).

My rationale for this design:

I split subjects and samples into separate tables because a single subject can have multiple samples collected over different times. Therefore, it felt cleaner to model it explicitly instead of flattening everything into one table.

For treatments, I did not store treatment and response directly on the subject row. Instead, I used a junction table (subject_treatment). This allows for cases where a subject is enrolled in more than one treatment arm. It also keeps the subject table from becoming over loaded with fields that do not always apply one to one.

Projects are in their own table instead of being a text column on the subject. The main reason is future flexibility. If we wanted to attach metadata at the project level later, such as ownership or timelines (etc), we can do that without changing the subject schema.

For cell counts, I chose to store them directly on the sample row as raw integers instead of using a long table. All five populations are always present and always queried together, so a long format would just add extra joins and increase the number of rows without adding real value. This keeps queries simpler and more efficient for the current use case.

Right now this project runs on SQLite, which works well into the tens of thousands of samples (+/-). The schema is designed so that adding new projects, subjects, or samples is straightforward and does not require changes. The junction table already supports multi arm treatment scenarios.

If this was needed to scale to hundreds of projects and much larger datasets, I would make a few of the following changes:

Move from SQLite to PostgreSQL for better concurrency, connection pooling, and indexing support.
Add indexes on commonly filtered columns such as subject.condition, sample.sample_type, and sample.time_from_treatment_start.
Partition the sample table by project_id if individual projects grow very large.

If the number of measured populations grows beyond the current five, I would revisit how cell counts are stored. That could mean moving them into a separate wide table or normalizing into a (sample_measurement) table with columns like sample_id, marker, and value.

On the app side, I am currently using a simple lru_cache in analysis.py based on the database path and modification time. That works fine at this scale, but for larger workloads I would replace it with a proper query cache for more complex queries.


Code structure

config.py
    Single source of truth for DB_PATH and POPULATIONS. Both analysis.py and
    load_data.py import from here to avoid duplication and ensure consistency.

load_data.py
    Reads cell-count.csv (TSV or CSV auto-detected),
    normalises into the five-table schema, and writes cell_counts.db.
    uses INSERT OR IGNORE so re-running does not duplicate data.

analysis.py
    All data access and analysis logic. Functions are grouped by concern: frequency
    table computation, responder frequency computation, statistical testing, subset
    summaries, longitudinal data, and the logistic regression model. The four main
    query functions are wrapped in functools.lru_cache keyed on (db_path, mtime)
    so results are reused across tab switches but automatically invalidated when
    the database file changes. The figure building functions are kept separate from
    the data fetching functions so data can be serialised into the dashboard store
    independently of the Plotly objects.

dashboard.py
    Layout is defined once at module level. All data is loaded
    in the populate_store callback on page load and serialised into a dcc.Store
    component. The render_tab callback reads from the store and builds the tab
    content from the cached data. There is no database access happening during tab switching.
    suppress_callback_exceptions=True is set because tab content components are
    created dynamically and are not present in the initial layout.

export_schema.py
    Renders the database schema as a static PNG using Plotly shapes and annotations. 
    Run once (or via make schema) and referenced in the README. 
    Kept out of the dashboard so potential "stakeholders" see findings, not infrastructure diagrams.

tests/test_analysis.py
    Tests covering frequency computation (percentages sum to 100) and
    statistical output shape (correct columns including p_value_adj, correct
    number of rows).


--- TODO (not implemented) ---

The statistical analysis tab triggers a longitudinal data query on every tab 3 switch, 
because make_longitudinal_figure() currently fetches directly from the
database rather than reading from the dcc.Store. The fix is to call
compute_longitudinal() in populate_store, serialise the result into the store
alongside the other dataframes, and pass the pre-loaded DataFrame into
make_longitudinal_figure(long_df) in render_tab. This brings the stats tab in
line with all other tabs, which have no database access after initial page load.
