from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
TABLE_ROOT = ROOT / "02_processed_source_tables"
FIG_ROOT = ROOT / "04_figure_source_data"

checks = []

def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

geo = FIG_ROOT / "Fig6_GEO_7Cohort_C1Q_PerSD_HR_Source.csv"
rows = read_csv(geo)
checks.append(("Fig6 GEO source rows >= 7", len(rows) >= 7))
summary = FIG_ROOT / "Fig6_GEO_Random_Effects_Summary.csv"
text = geo.read_text(encoding="utf-8-sig") + "\n" + summary.read_text(encoding="utf-8-sig")
for token in ["1.06", "0.95", "1.18", "55.9"]:
    checks.append((f"reported external summary token present: {token}", token in text))

required = [
    TABLE_ROOT / "Tables" / "CSV" / "S1_Table_Model_Benchmarking.csv",
    TABLE_ROOT / "Tables" / "CSV" / "S2_Table_Feature_Selection_Score_Construction.csv",
    TABLE_ROOT / "Tables" / "CSV" / "S3_Table_Spatial_Autocorrelation_Hotspots.csv",
    TABLE_ROOT / "Tables" / "CSV" / "S10_Table_Immunotherapy_Context_Trends.csv",
]
for path in required:
    checks.append((f"required table exists: {path.relative_to(ROOT)}", path.exists()))

failed = [name for name, ok in checks if not ok]
for name, ok in checks:
    print(f"{'PASS' if ok else 'FAIL'} - {name}")
if failed:
    raise SystemExit("Validation failed: " + "; ".join(failed))
print("All source-table checks passed.")
