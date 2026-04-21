Loblaw Bio: Miraclib Clinical Trial Analysis

This project analyzes immune cell populations in melanoma patients treated with Miraclib. It provides a pipeline that loads raw patient data, performs statistical comparisons between responders and non-responders, and exposes the results through an interactive dashboard.

The analysis focuses on computing immune cell frequencies, comparing groups using non-parametric statistical testing, and evaluating baseline cohort behavior at time point zero. Statistical significance is assessed using the Mann Whitney U test, with multiple comparison correction applied using the Benjamini Hochberg method. Effect sizes are also computed using rank biserial correlation to give additional context beyond p-values.

Data is taken from a single input file (cell-count.csv). The pipeline automatically detects the format, processes the data, and builds a local database (cell_counts.db) to support analysis and reproducibility.

The project is split into three main scripts:

1) load_data.py handles ingestion and database creation
2) analysis.py performs statistical comparisons and generates derived results
3) dashboard.py launches an interface for further viewing

To run the full pipeline:
1) install dependencies with:
    pip install -r requirements.txt

2) execute the scripts in order: 
    python3 load_data.py
    python3 analysis.py
    python3 dashboard.py 

3) Access locally at http://127.0.0.1:8050.

The output is designed to be fully reproducible from raw input to visualization, with all intermediate computations stored locally for traceability.