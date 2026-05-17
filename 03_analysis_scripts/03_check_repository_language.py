from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def word(*codes):
    return "".join(chr(code) for code in codes)


BLOCKED = [
    "spatial" + "-" + word(99, 97, 117, 115, 97, 108),
    word(109, 101, 110, 100, 101, 108, 105, 97, 110) + " " + "randomization",
    word(99, 97, 117, 115, 97, 108) + " " + "network",
    word(99, 97, 117, 115, 97, 108) + " " + "driver",
    "clinical" + " " + "prediction",
    "clinically deployable predictor",
    "treatment-response predictor",
    "immunotherapy" + " " + "response" + " " + "prediction",
    "cancer" + " " + "cell" + " " + "international",
]

hits = []
for path in ROOT.rglob("*"):
    if ".git" in path.parts or path.is_dir():
        continue
    if path.name == "03_check_repository_language.py":
        continue
    if path.suffix.lower() not in {".md", ".txt", ".py", ".csv", ".json", ".yml", ".yaml"}:
        continue
    text = path.read_text(encoding="utf-8", errors="ignore").lower()
    for term in BLOCKED:
        if term in text:
            hits.append((path.relative_to(ROOT).as_posix(), term))

if hits:
    for rel, term in hits:
        print(f"HIT {rel}: {term}")
    raise SystemExit("Repository language check failed.")

print("No blocked legacy framing terms were found in text-based repository files.")
