.PHONY: install db dashboard analysis clean help

# Config
CSV    ?= cell-count.csv
DB     ?= cell_counts.db
PYTHON  = python3

# Default target
help:
	@echo ""
	@echo "  Loblaw Bio : Miraclib Trial Analysis"
	@echo ""
	@echo "  Usage:"
	@echo "    make install     Install Python dependencies"
	@echo "    make db          Load cell-count.csv into SQLite ($(DB))"
	@echo "    make dashboard   Launch the interactive Dash dashboard"
	@echo "    make analysis    Run Parts 2–4 and print results to terminal"
	@echo "    make all         install → db → dashboard"
	@echo "    make clean       Remove the generated database and plot files"
	@echo ""
	@echo "  Options:"
	@echo "    CSV ?= data/cell-count.csv   Override input CSV  (default: $(CSV))"
	@echo "    DB ?= data/cell_counts.db      Override database path (default: $(DB))"
	@echo ""

# Install 
install:
	$(PYTHON) -m pip install -r requirements.txt

# Load data
db: $(CSV)
	$(PYTHON) load_data.py --csv $(CSV) --db $(DB)

$(CSV):
	@echo "ERROR: $(CSV) not found. Please place your cell-count CSV file in the project root."
	@exit 1

# Dashboard
dashboard: $(DB)
	$(PYTHON) dashboard.py

# CLI analysis
analysis: $(DB)
	$(PYTHON) analysis.py

# Full pipeline
all: install db dashboard

# Tests
test:
	$(PYTHON) -m pytest -v

# Clean
clean:
	rm -f $(DB) boxplot.png
	@echo "Removed $(DB) and boxplot.png"