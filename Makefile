.PHONY: setup pipeline dashboard install db schema analysis test clean help bootstrap

# Config
CSV    ?= data/cell-count.csv
DB     ?= cell_counts.db
PYTHON  = python3

# Default target
help:
	@echo ""
	@echo "  Loblaw Bio · Miraclib Trial"
	@echo ""
	@echo "  Required targets:"
	@echo "    make setup       Install all dependencies"
	@echo "    make pipeline    Run full pipeline: load data → analysis → schema"
	@echo "    make dashboard   Launch dashboard at http://127.0.0.1:8050"
	@echo ""
	@echo "  Additional targets:"
	@echo "    make test        Run pytest"
	@echo "    make clean       Remove generated files"
	@echo ""
	@echo "  Overrides:"
	@echo "    CSV=$(CSV)"
	@echo "    DB=$(DB)"
	@echo ""

# Setup
setup:
	$(PYTHON) -m pip install -r requirements.txt

# Pipeline
pipeline: $(CSV)
	$(PYTHON) load_data.py --csv $(CSV) --db $(DB)
	$(PYTHON) analysis.py
	$(PYTHON) export_schema.py

$(CSV):
	@echo "ERROR: $(CSV) not found. Place cell-count.csv in data/."
	@exit 1

# Dashbaord
dashboard: $(DB)
	$(PYTHON) dashboard.py

# Aliases
install: setup
db: pipeline
bootstrap: setup pipeline

# Schema PNG
schema: $(DB)
	$(PYTHON) export_schema.py

# CLI analysis
analysis: $(DB)
	$(PYTHON) analysis.py

# Tests
test:
	$(PYTHON) -m pytest tests/ -v

# Clean
clean:
	rm -f $(DB) docs/schema.png
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
	@echo "Clean."