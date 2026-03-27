#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    roc_auc_score,
    confusion_matrix,
    classification_report,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler


def read_matrix_auto(path: str, sample_ids: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """
    Read a numeric matrix (TSV/CSV). Tries to auto-detect orientation:
    - If columns intersect sample_ids more than rows do, assume columns are samples (genes/species x samples) and transpose to samples x features.
    - Else, if rows intersect sample_ids, assume rows are samples already.
    - If sample_ids is None, assume columns are samples and transpose (common for counts).

    Returns a DataFrame indexed by samples with features as columns.
    """
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, header=0, index_col=0)

    if sample_ids is None:
        # Heuristic: typical omics matrices are features x samples
        return df.T

    sample_ids = set(map(str, sample_ids))
    cols_match = len(sample_ids.intersection(map(str, df.columns)))
    rows_match = len(sample_ids.intersection(map(str, df.index)))

    if cols_match >= rows_match:
        return df.T
    else:
        return df


essential_bracken_cols = {
    "name",
    "taxonomy_id",
    "taxonomy_lvl",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
    "fraction_total_reads",
}


def build_bracken_matrix_from_dir(
    bracken_dir: str,
    level: str = "S",
    value_col: str = "fraction_total_reads",
) -> pd.DataFrame:
    """
    Build a samples x species abundance matrix from Bracken output files.
    Accepts files ending with ".bracken.tsv" or ".tsv" in the directory.
    Filters by taxonomy level (default species 'S').
    value_col: typically 'fraction_total_reads' (relative abundance) or 'new_est_reads'.
    """
    bracken_dir = Path(bracken_dir)
    files = sorted(list(bracken_dir.glob("*.bracken.tsv")))
    if not files:
        files = sorted(list(bracken_dir.glob("*.tsv")))
    if not files:
        raise FileNotFoundError(f"No Bracken TSV files found in {bracken_dir}")

    matrices = []
    sample_names = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception as e:
            raise RuntimeError(f"Failed to read {f}: {e}")
        # Validate columns
        if not set(df.columns).issuperset(essential_bracken_cols):
            raise ValueError(
                f"Bracken file {f} missing required columns; found {list(df.columns)}"
            )
        df = df[df["taxonomy_lvl"] == level]
        # Use species names as index; in case of duplicates, aggregate by sum/mean
        agg = (
            df.groupby("name")[value_col]
            .sum()
            .to_frame(name=f.stem)
        )
        matrices.append(agg)
        sample_names.append(f.stem)

    wide = pd.concat(matrices, axis=1).fillna(0.0)
    # Now wide is species x samples; transpose to samples x species
    return wide.T


def align_inputs(
    X_a: Optional[pd.DataFrame], X_b: Optional[pd.DataFrame]
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Align two sample-indexed DataFrames on shared samples. If one is None, return as is.
    """
    if X_a is None and X_b is None:
        raise ValueError("At least one feature matrix must be provided")
    if X_b is None:
        return X_a, None
    if X_a is None:
        return X_b, None
    common = X_a.index.intersection(X_b.index)
    if len(common) == 0:
        raise ValueError("No overlapping samples between feature sets")
    return X_a.loc[common], X_b.loc[common]


def make_labels(metadata: pd.DataFrame, sample_col: str, group_col: str,
                positive: Optional[str], negative: Optional[str]) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Produce binary labels y from metadata. If positive/negative provided, filter to those values.
    Returns (y, filtered_metadata)
    """
    if sample_col not in metadata.columns:
        raise KeyError(f"Metadata missing sample column '{sample_col}'")
    if group_col not in metadata.columns:
        raise KeyError(f"Metadata missing group column '{group_col}'")

    md = metadata.copy()
    if positive is not None and negative is not None:
        md = md[md[group_col].isin([positive, negative])].copy()
        mapping = {negative: 0, positive: 1}
        y = md[group_col].map(mapping)
    else:
        # Assume numeric/binary already (0/1)
        vals = md[group_col].astype(str)
        unique_vals = sorted(vals.unique())
        if len(unique_vals) != 2:
            raise ValueError(
                f"Group column '{group_col}' must be binary or provide --positive/--negative. Found: {unique_vals}"
            )
        mapping = {unique_vals[0]: 0, unique_vals[1]: 1}
        y = vals.map(mapping)

    y.index = md[sample_col].astype(str)
    return y, md


def filter_and_scale(
    X: pd.DataFrame,
    min_prevalence: float = 0.0,
    apply_zscore: bool = False,
    min_variance: float = 0.0,
) -> pd.DataFrame:
    """
    - Drop columns with prevalence (non-zero proportion) < min_prevalence
    - Drop near-zero variance features using VarianceThreshold
    - Optionally z-score scale features (mean=0, std=1)
    """
    X_proc = X.copy()

    if min_prevalence > 0:
        nonzero_prop = (X_proc != 0).sum(axis=0) / X_proc.shape[0]
        keep = nonzero_prop[nonzero_prop >= min_prevalence].index
        X_proc = X_proc[keep]

    if min_variance > 0:
        vt = VarianceThreshold(threshold=min_variance)
        X_proc = pd.DataFrame(
            vt.fit_transform(X_proc.values),
            index=X_proc.index,
            columns=X_proc.columns[vt.get_support(indices=True)],
        )

    if apply_zscore:
        scaler = StandardScaler()
        X_proc.iloc[:, :] = scaler.fit_transform(X_proc.values)

    return X_proc


def run_rf_cv(
    X: pd.DataFrame,
    y: pd.Series,
    n_splits: int = 5,
    n_estimators: int = 1000,
    max_depth: Optional[int] = None,
    class_weight: Optional[str] = "balanced",
    random_state: int = 42,
) -> Tuple[pd.DataFrame, dict, pd.Series]:
    """
    Train RF with StratifiedKFold CV. Returns:
    - fold_metrics DataFrame
    - summary dict
    - feature_importances (mean across folds)
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    fold_rows = []
    importances = []
    cm_total = np.zeros((2, 2), dtype=int)

    for fold, (train_idx, test_idx) in enumerate(skf.split(X.values, y.values), start=1):
        X_train, X_test = X.values[train_idx], X.values[test_idx]
        y_train, y_test = y.values[train_idx], y.values[test_idx]

        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            class_weight=class_weight,
            random_state=random_state + fold,
            n_jobs=-1,
        )
        clf.fit(X_train, y_train)

        y_pred = clf.predict(X_test)
        acc = accuracy_score(y_test, y_pred)

        # ROC-AUC (binary)
        try:
            y_proba = clf.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, y_proba)
        except Exception:
            auc = np.nan

        cm = confusion_matrix(y_test, y_pred, labels=[0, 1])
        cm_total += cm

        fold_rows.append({
            "fold": fold,
            "accuracy": acc,
            "roc_auc": auc,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
        })

        imp = pd.Series(clf.feature_importances_, index=X.columns)
        importances.append(imp)

    fold_metrics = pd.DataFrame(fold_rows)
    mean_importance = pd.concat(importances, axis=1).mean(axis=1).sort_values(ascending=False)

    summary = {
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "n_splits": int(n_splits),
        "mean_accuracy": float(np.nanmean(fold_metrics["accuracy"].values)),
        "std_accuracy": float(np.nanstd(fold_metrics["accuracy"].values)),
        "mean_roc_auc": float(np.nanmean(fold_metrics["roc_auc"].values)),
        "std_roc_auc": float(np.nanstd(fold_metrics["roc_auc"].values)),
        "confusion_matrix_total": cm_total.tolist(),
    }

    return fold_metrics, summary, mean_importance


def save_results(
    outdir: str,
    fold_metrics: pd.DataFrame,
    summary: dict,
    feature_importances: pd.Series,
    top_k: int = 100,
):
    Path(outdir).mkdir(parents=True, exist_ok=True)

    fold_metrics.to_csv(os.path.join(outdir, "cv_folds_metrics.csv"), index=False)
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    feature_importances.to_frame("importance").to_csv(
        os.path.join(outdir, "feature_importances.csv")
    )
    feature_importances.head(top_k).to_frame("importance").to_csv(
        os.path.join(outdir, f"top{top_k}_feature_importances.csv")
    )

    # Save a simple text report
    with open(os.path.join(outdir, "report.txt"), "w") as f:
        f.write(
            "RF CV Classification Report\n"
            f"Samples: {summary['n_samples']}, Features: {summary['n_features']}, Folds: {summary['n_splits']}\n"
            f"Accuracy (mean±std): {summary['mean_accuracy']:.3f} ± {summary['std_accuracy']:.3f}\n"
            f"ROC-AUC (mean±std): {summary['mean_roc_auc']:.3f} ± {summary['std_roc_auc']:.3f}\n"
            f"Confusion matrix (aggregated over folds) [rows=true 0/1, cols=pred 0/1]:\n{np.array(summary['confusion_matrix_total'])}\n"
        )


def main():
    p = argparse.ArgumentParser(
        description=(
            "Random Forest classifier using microbial species abundances and/or human transcriptome profiles.\n"
            "Inputs should be sample x feature matrices (auto-detected if transposed)."
        )
    )
    p.add_argument("--metadata", required=True, help="Path to metadata file (CSV/TSV)")
    p.add_argument("--metadata_sep", default="auto", choices=["auto", "tab", "comma", "semicolon"], help="Metadata separator")
    p.add_argument("--sample_col", default="sample", help="Column in metadata with sample IDs")
    p.add_argument("--group_col", default="group", help="Column in metadata with case/control labels")
    p.add_argument("--positive", default=None, help="Label value for positive class (case)")
    p.add_argument("--negative", default=None, help="Label value for negative class (control)")

    p.add_argument("--expr_tsv", default=None, help="Path to transcriptome matrix (genes x samples or samples x genes)")
    p.add_argument("--species_tsv", default=None, help="Path to species abundance matrix (species x samples or samples x species)")
    p.add_argument("--species_dir", default=None, help="Directory with Bracken TSV files to aggregate into a species matrix")
    p.add_argument("--bracken_level", default="S", help="Bracken taxonomy level to use (e.g., S)")
    p.add_argument("--bracken_value_col", default="fraction_total_reads", choices=["fraction_total_reads", "new_est_reads"], help="Bracken value column to use for abundance")

    p.add_argument("--feature_set", default="both", choices=["species", "expression", "both"], help="Which features to use")
    p.add_argument("--log1p_expr", action="store_true", help="Apply log1p transform to expression features")
    p.add_argument("--log1p_species", action="store_true", help="Apply log1p transform to species features")
    p.add_argument("--min_prevalence", type=float, default=0.0, help="Drop features seen (non-zero) in < this proportion of samples")
    p.add_argument("--min_variance", type=float, default=0.0, help="Drop features with variance <= threshold")
    p.add_argument("--zscore", action="store_true", help="Z-score scale features")

    p.add_argument("--n_estimators", type=int, default=1000)
    p.add_argument("--max_depth", type=int, default=None)
    p.add_argument("--cv", type=int, default=5, help="Number of StratifiedKFold splits")
    p.add_argument("--random_state", type=int, default=42)

    p.add_argument("--outdir", default="results/rf_classifier", help="Output directory")

    args = p.parse_args()

    # Load metadata
    if args.metadata_sep == "auto":
        with open(args.metadata, "r", errors="ignore") as fh:
            head = fh.readline()
        if "\t" in head:
            md_sep = "\t"
        elif ";" in head and "," not in head:
            md_sep = ";"
        else:
            md_sep = ","
    else:
        md_sep = {"tab": "\t", "comma": ",", "semicolon": ";"}[args.metadata_sep]

    metadata = pd.read_csv(args.metadata, sep=md_sep, header=0)

    # Prepare labels
    y, md_filt = make_labels(
        metadata, sample_col=args.sample_col, group_col=args.group_col,
        positive=args.positive, negative=args.negative
    )

    # Load expression
    X_expr = None
    if args.expr_tsv:
        X_expr = read_matrix_auto(args.expr_tsv, sample_ids=y.index)
        X_expr = X_expr.loc[X_expr.index.intersection(y.index)]
        X_expr = X_expr.select_dtypes(include=[np.number]).fillna(0.0)
        if args.log1p_expr:
            X_expr = np.log1p(X_expr)
        X_expr = filter_and_scale(
            X_expr, min_prevalence=args.min_prevalence,
            apply_zscore=args.zscore, min_variance=args.min_variance
        )

    # Load species from matrix or build from directory
    X_sp = None
    if args.species_tsv:
        X_sp = read_matrix_auto(args.species_tsv, sample_ids=y.index)
        X_sp = X_sp.loc[X_sp.index.intersection(y.index)]
        X_sp = X_sp.select_dtypes(include=[np.number]).fillna(0.0)
    elif args.species_dir:
        X_sp = build_bracken_matrix_from_dir(
            args.species_dir, level=args.bracken_level, value_col=args.bracken_value_col
        )
        X_sp = X_sp.loc[X_sp.index.intersection(y.index)]
    if X_sp is not None:
        if args.log1p_species:
            X_sp = np.log1p(X_sp)
        X_sp = filter_and_scale(
            X_sp, min_prevalence=args.min_prevalence,
            apply_zscore=args.zscore, min_variance=args.min_variance
        )

    # Decide feature set
    feature_set = args.feature_set
    if feature_set == "species" and X_sp is None:
        raise ValueError("--feature_set species selected but no species features were provided")
    if feature_set == "expression" and X_expr is None:
        raise ValueError("--feature_set expression selected but no expression features were provided")
    if feature_set == "both" and (X_sp is None or X_expr is None):
        raise ValueError("--feature_set both selected but missing species or expression features")

    # Align samples and merge
    if feature_set == "species":
        X_all, _ = align_inputs(X_sp, None)
    elif feature_set == "expression":
        X_all, _ = align_inputs(X_expr, None)
    else:
        Xa, Xb = align_inputs(X_expr, X_sp)
        X_all = pd.concat([Xa, Xb], axis=1)

    # Align with labels
    common_samples = X_all.index.intersection(y.index)
    if len(common_samples) == 0:
        raise ValueError("No overlapping samples between features and labels")
    X_all = X_all.loc[common_samples]
    y_all = y.loc[common_samples]

    # Run RF CV
    fold_metrics, summary, feat_imp = run_rf_cv(
        X_all, y_all,
        n_splits=args.cv,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        class_weight="balanced",
        random_state=args.random_state,
    )

    # Save outputs
    save_results(args.outdir, fold_metrics, summary, feat_imp)

    # Print headline result
    print(f"Mean accuracy: {summary['mean_accuracy']:.3f} (± {summary['std_accuracy']:.3f})")
    if not np.isnan(summary["mean_roc_auc"]):
        print(f"Mean ROC-AUC: {summary['mean_roc_auc']:.3f} (± {summary['std_roc_auc']:.3f})")


if __name__ == "__main__":
    main()