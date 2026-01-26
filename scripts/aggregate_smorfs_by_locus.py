#!/usr/bin/env python3
# I run it with this command, but change the directories as needed:
# python aggregate_smorfs_by_locus.py \
#   --merged_dir /storage/scratch01/users/sbarber/Visceral/merged_per_sample \
#   --out_prefix /storage/scratch01/users/sbarber/Visceral/smorf_locus_summary

import argparse
from pathlib import Path
import pandas as pd


LOCUS_COLS = ["cds_chr", "cds_starts", "cds_ends", "cds_strand"]


def first_nonnull(series: pd.Series):
    """Return the first non-null value (as string) or NA if none."""
    s = series.dropna()
    if len(s) == 0:
        return pd.NA
    return s.iloc[0]

def unique_nonnull(series: pd.Series):
    s = [str(x) for x in series.dropna() if str(x) != ""]
    if not s:
        return pd.NA
    return ",".join(sorted(set(s)))

def most_common_nonnull(series: pd.Series):
    s = [str(x) for x in series.dropna() if str(x) != ""]
    if not s:
        return pd.NA
    counts = {}
    for val in s:
        counts[val] = counts.get(val, 0) + 1
    max_count = max(counts.values())
    top = sorted([k for k, v in counts.items() if v == max_count])
    return ",".join(top)


def main():
    ap = argparse.ArgumentParser(
        description="Aggregate per-sample merged ShortStop tables by genomic locus and find smORFs shared across patients."
    )
    ap.add_argument("--merged_dir", required=True,
                    help="Directory containing per-sample merged CSVs (e.g., merged_per_sample/)")
    ap.add_argument("--out_prefix", required=True,
                    help="Prefix for output files (e.g., smorf_locus_summary)")
    ap.add_argument("--min_patients", type=int, default=2,
                    help="Keep loci observed in at least this many patients (default: 2)")
    args = ap.parse_args()

    merged_dir = Path(args.merged_dir)
    files = sorted(merged_dir.glob("*.merged.csv"))
    if not files:
        raise SystemExit(f"No *.merged.csv files found in: {merged_dir}")

    rows = []
    for f in files:
        sample = f.name.replace(".merged.csv", "")
        df = pd.read_csv(f)

        missing = [c for c in LOCUS_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"{f}: missing required locus columns: {missing}")

        # ensure numeric starts/ends if possible
        df["cds_starts"] = pd.to_numeric(df["cds_starts"], errors="coerce")
        df["cds_ends"]   = pd.to_numeric(df["cds_ends"], errors="coerce")

        # locus key
        df["locus"] = (
            df["cds_chr"].astype(str) + ":" +
            df["cds_starts"].astype("Int64").astype(str) + "-" +
            df["cds_ends"].astype("Int64").astype(str) + ":" +
            df["cds_strand"].astype(str)
        )

        # keep key columns + some useful fields if present
        keep_cols = ["locus"] + LOCUS_COLS
        for extra in ["orf_id", "sam_probability", "classification", "aa_seq", "length", "type", "cds_seq", "smorf_type"]:
            if extra in df.columns:
                keep_cols.append(extra)

        sub = df[keep_cols].copy()
        sub["sample"] = sample
        rows.append(sub)

    all_df = pd.concat(rows, ignore_index=True)

    # Build aggregation spec (conditionally include optional columns)
    agg_spec = dict(
        cds_chr=("cds_chr", "first"),
        cds_starts=("cds_starts", "first"),
        cds_ends=("cds_ends", "first"),
        cds_strand=("cds_strand", "first"),
        n_patients=("sample", "nunique"),
        patients=("sample", lambda x: ",".join(sorted(set(x)))),
        n_rows=("sample", "size"),
    )

    if "sam_probability" in all_df.columns:
        agg_spec["max_prob"] = ("sam_probability", "max")
    else:
        agg_spec["max_prob"] = ("sample", "size")

    if "cds_seq" in all_df.columns:
        agg_spec["cds_seq"] = ("cds_seq", first_nonnull)
    else:
        agg_spec["cds_seq"] = ("sample", "size")

    # Add aa_seq to the locus-level output (first non-null across rows)
    if "aa_seq" in all_df.columns:
        agg_spec["aa_seq"] = ("aa_seq", first_nonnull)

    if "smorf_type" in all_df.columns:
        agg_spec["smorf_type"] = ("smorf_type", most_common_nonnull)
        agg_spec["smorf_types"] = ("smorf_type", unique_nonnull)

    # Aggregate across samples per locus
    agg = (
        all_df.groupby("locus", as_index=False)
        .agg(**agg_spec)
    )

    # Output full summary
    full_out = Path(f"{args.out_prefix}.all_loci.csv")
    agg.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True]
    ).to_csv(full_out, index=False)

    # Output loci shared across >= min_patients
    shared = agg[agg["n_patients"] >= args.min_patients].copy()
    shared_out = Path(f"{args.out_prefix}.shared_ge{args.min_patients}.csv")
    shared.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True, True]
    ).to_csv(shared_out, index=False)

    print(f"Total loci: {len(agg)}")

    # Print additional shared-loci counts, for reference when running manually
    for k in [2, 5, 10, 15]:
        n_k = (agg["n_patients"] >= k).sum()
        print(f"Shared loci (>= {k} patients): {n_k}")

    print(f"[OK] Wrote: {full_out}")
    print(f"[OK] Wrote: {shared_out}")


if __name__ == "__main__":
    main()
