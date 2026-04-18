from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier


def load_model_data(root: Path, top_n_features: int = 100):
    df = pd.read_csv(root / "integrated_data.csv")
    feat = pd.read_csv(root / "multi_omics_feature_importance.csv")
    df = df.dropna(subset=["OS.time", "OS", "response"]).copy()
    df = df[df["OS.time"] > 0].copy()
    df["response"] = df["response"].astype(int)
    ranked = feat["feature"].astype(str).tolist()
    valid_features = [f for f in ranked if f in df.columns][:top_n_features]
    X = df[valid_features].copy()
    y = df["response"].to_numpy(dtype=int)
    return X, y


def main():
    parser = argparse.ArgumentParser(description="Append XGBoost repeated-CV AUC distribution to the summary files.")
    parser.add_argument("--root", type=str, default=".", help="Workspace root containing integrated_data.csv and feature files.")
    parser.add_argument("--out-dir", type=str, default=None, help="Output directory containing model_auc_distribution.csv.")
    parser.add_argument("--top-n-features", type=int, default=100, help="Number of ranked features to include.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (root / "supplementary_figures_rebuild")

    X, y = load_model_data(root, args.top_n_features)
    splitter = RepeatedStratifiedKFold(n_splits=5, n_repeats=20, random_state=42)
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                XGBClassifier(
                    random_state=42,
                    n_estimators=100,
                    learning_rate=0.1,
                    max_depth=3,
                    eval_metric="logloss",
                ),
            ),
        ]
    )

    rows = []
    for split_idx, (train_idx, test_idx) in enumerate(splitter.split(X, y), start=1):
        model.fit(X.iloc[train_idx], y[train_idx])
        y_prob = model.predict_proba(X.iloc[test_idx])[:, 1]
        auc = roc_auc_score(y[test_idx], y_prob)
        rows.append(
            {
                "model": "XGBoost",
                "split_id": split_idx,
                "auc": float(auc),
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
            }
        )

    auc_path = out_dir / "model_auc_distribution.csv"
    old = pd.read_csv(auc_path)
    old = old[old["model"].astype(str) != "XGBoost"].copy()
    merged = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    merged.to_csv(auc_path, index=False)

    summary = (
        merged.groupby("model", as_index=False)["auc"]
        .agg(["mean", "std", "median", "min", "max"])
        .reset_index()
        .rename(
            columns={
                "mean": "auc_mean",
                "std": "auc_std",
                "median": "auc_median",
                "min": "auc_min",
                "max": "auc_max",
            }
        )
    )
    summary.to_csv(out_dir / "model_auc_distribution_summary.csv", index=False)
    print(f"Updated AUC distribution files in: {out_dir}")


if __name__ == "__main__":
    main()
