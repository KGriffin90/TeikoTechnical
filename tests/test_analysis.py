import numpy as np
import pandas as pd

from analysis import (
    compute_frequency_table,
    compute_responder_frequencies,
    run_statistical_tests,
)

# Verify that percentages for each sample sum to 100%
def test_frequencies_sum_to_100():
    df = compute_frequency_table()
    sums = df.groupby("sample")["percentage"].sum()
    assert np.allclose(sums, 100.0, atol=0.01)

# Validate that statistical output contains required columns and correct row count
def test_stats_shape():
    freq = compute_responder_frequencies()
    stats_df = run_statistical_tests(freq)

    assert set(stats_df.columns) >= {
        "population",
        "p_value",
        "effect_size",
    }

    assert len(stats_df) == 5