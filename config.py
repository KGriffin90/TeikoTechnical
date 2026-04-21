"""
config.py
   
Single source of truth for shared constants.
Both analysis.py and load_data.py import from here.
"""

DB_PATH = "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]