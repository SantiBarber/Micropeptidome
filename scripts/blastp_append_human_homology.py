#!/usr/bin/env python3
import argparse
import re
import subprocess
import shutil
from pathlib import Path
import pandas as pd

AA_RE = re.compile(r"^[A-Za-z\*]+$")

def clean_seq(s: str) -> str:
    s = str(s).strip()
    if s.endswith("*"):
        s = s[:-1]
    return s.upper()

def run(cmd, **kwargs):
    p = subprocess.run(cmd, text=True, capture_output=True, **kwargs)
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  {' '.join(cmd)}\n\nSTDOUT:\n{p.stdout}\n\nSTDERR:\n{p.stderr}\n"
        )
    return p

def main():
    ap = argparse.ArgumentParser(description="BLASTP aa_seq vs human proteome and append best hit to CSV.")
    ap.add_argument("--in_csv", required=True, help="Input locus summary CSV (must contain locus and aa_seq).")
    ap.add_argument("--out_csv", required=True, help="Output CSV with appended BLASTP best-hit columns.")
    ap.add_argument("--db", required=True, help="BLAST database prefix (as used with -db).")
    ap.add_argument(
        "--blastp",
        default="blastp",
        help="blastp executable (default: resolve 'blastp' from PATH). You may also pass a full path.",
    )
    ap.add_argument("--evalue", type=float, default=1e-3, help="E-value cutoff.")
    ap.add_argument("--threads", type=int, default=4, help="blastp threads.")
    ap.add_argument("--max_targets", type=int, default=25, help="How many target hits to keep per query.")
    ap.add_argument("--seg", default="no", choices=["yes", "no"], help="Low complexity filtering (default no for short peptides).")
    args = ap.parse_args()

    # Resolve blastp to an absolute path to avoid PATH issues such as ENOTDIR
    # when PATH contains non-directory entries (common in some HPC module setups).
    blastp_exe = args.blastp
    if "/" not in blastp_exe:
        resolved = shutil.which(blastp_exe)
        if resolved is None:
            raise RuntimeError(
                f"Could not find executable '{blastp_exe}' in PATH. "
                "Run this rule in the conda env that provides BLAST+, or pass --blastp /full/path/to/blastp."
            )
        blastp_exe = resolved

    in_csv = Path(args.in_csv)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # low_memory=False avoids pandas DtypeWarning for large mixed-type CSVs
    df = pd.read_csv(in_csv, low_memory=False)
    if "locus" not in df.columns or "aa_seq" not in df.columns:
        raise ValueError("Input CSV must contain columns: locus, aa_seq")

    # Build query FASTA
    tmp_dir = out_csv.parent / (out_csv.stem + ".blast_tmp")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    query_fa = tmp_dir / "queries.fa"
    blast_tsv = tmp_dir / "blast.tsv"

    qmap = {}  # qid -> locus
    records = []

    for _, row in df[["locus", "aa_seq"]].iterrows():
        locus = str(row["locus"])
        aa = row["aa_seq"]
        if pd.isna(aa):
            continue
        seq = clean_seq(aa)
        if not seq or len(seq) < 5:
            continue
        if not AA_RE.match(seq):
            continue

        # qseqid must not contain spaces for robust parsing
        qid = locus.replace(" ", "_")
        qmap[qid] = locus
        records.append((qid, seq))

    with open(query_fa, "w") as fh:
        for qid, seq in records:
            fh.write(f">{qid}\n{seq}\n")

    # If no sequences, just write NA columns and exit
    if not records:
        df["human_best_pident"] = pd.NA
        df["human_best_qcov"] = pd.NA
        df["human_best_hit_id"] = pd.NA
        df["human_best_hit_title"] = pd.NA
        df.to_csv(out_csv, index=False)
        return

    # BLASTP
    # Include qlen so we can compute query coverage; stitle gives protein description line
    outfmt = "6 qseqid sseqid pident length qlen evalue bitscore stitle"
    cmd = [
        blastp_exe,
        "-query", str(query_fa),
        "-db", args.db,
        "-evalue", str(args.evalue),
        "-max_target_seqs", str(args.max_targets),
        "-num_threads", str(args.threads),
        "-seg", args.seg,
        "-outfmt", outfmt,
        "-out", str(blast_tsv),
    ]
    run(cmd)

    # Parse BLAST output
    if blast_tsv.stat().st_size == 0:
        best = pd.DataFrame(columns=["locus","human_best_pident","human_best_qcov","human_best_hit_id","human_best_hit_title"])
    else:
        cols = ["qseqid","sseqid","pident","aln_len","qlen","evalue","bitscore","stitle"]
        b = pd.read_csv(blast_tsv, sep="\t", header=None, names=cols)

        # Choose best hit per query by bitscore (then pident, then aln_len)
        b = b.sort_values(["qseqid","bitscore","pident","aln_len"], ascending=[True, False, False, False])
        b = b.drop_duplicates("qseqid", keep="first").copy()

        b["human_best_qcov"] = (b["aln_len"] / b["qlen"]) * 100.0
        b["locus"] = b["qseqid"].map(qmap)

        best = b[["locus"]].copy()
        best["human_best_pident"] = b["pident"]
        best["human_best_qcov"] = b["human_best_qcov"]
        best["human_best_hit_id"] = b["sseqid"]
        best["human_best_hit_title"] = b["stitle"]

    # Merge back
    df2 = df.merge(best, on="locus", how="left")
    df2.to_csv(out_csv, index=False)

if __name__ == "__main__":
    main()
