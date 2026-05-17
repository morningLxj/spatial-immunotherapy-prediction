from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
GEO = ROOT / "04_figure_source_data" / "Fig6_GEO_7Cohort_C1Q_PerSD_HR_Source.csv"

with GEO.open("r", encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))

print("External survival traceability source table")
print(f"Rows: {len(rows)}")
print("Columns:", ", ".join(rows[0].keys()) if rows else "none")

for row in rows:
    joined = " | ".join(f"{k}={v}" for k, v in row.items())
    print(joined)
