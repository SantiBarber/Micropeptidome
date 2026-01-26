#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def load_isoform_tpms(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    # RSEM isoforms.results normally has transcript_id + TPM
    id_col = "transcript_id" if "transcript_id" in df.columns else df.columns[0]
    if "TPM" not in df.columns:
        raise SystemExit(f"No TPM column found in {path}")
    return df.set_index(id_col)["TPM"].astype(float)

def add_tpms(summary_csv: Path, rsem_dir: Path, out_csv: Path) -> None:
    summ = pd.read_csv(summary_csv)
    if "locus" not in summ.columns or "patients" not in summ.columns:
        raise SystemExit(f"{summary_csv} must contain columns: locus, patients")

    # Load all per-sample TPM vectors once
    tpm_by_sample = {}
    for sample_dir in sorted(rsem_dir.glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        iso = sample_dir / f"{sample}.isoforms.results"
        if iso.exists():
            tpm_by_sample[sample] = load_isoform_tpms(iso)

    def row_tpms(row) -> str:
        locus = str(row["locus"])
        pats = str(row["patients"]) if pd.notna(row["patients"]) else ""
        pat_list = [p for p in pats.split(",") if p]
        vals = []
        for p in pat_list:
            s = tpm_by_sample.get(p)
            v = 0.0
            if s is not None and locus in s.index:
                v = float(s.loc[locus])
            vals.append(f"{v:.6f}")
        return ",".join(vals)

    summ["tpms"] = summ.apply(row_tpms, axis=1)
    summ.to_csv(out_csv, index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all_loci_csv", required=True)
    ap.add_argument("--shared_csv", required=True)
    ap.add_argument("--rsem_dir", required=True)
    ap.add_argument("--out_all_loci_csv", required=True)
    ap.add_argument("--out_shared_csv", required=True)
    args = ap.parse_args()

    rsem_dir = Path(args.rsem_dir)
    add_tpms(Path(args.all_loci_csv), rsem_dir, Path(args.out_all_loci_csv))
    add_tpms(Path(args.shared_csv),   rsem_dir, Path(args.out_shared_csv))

if __name__ == "__main__":
    main()
