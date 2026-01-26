#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci_csv", required=True) # the merged .csv file from the previous step
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tx2gene", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.loci_csv) # load the file

    required = {"locus", "cds_seq"} # these are the columns we need
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.loci_csv}: {sorted(missing)}")

    # Write FASTA
    with open(args.fasta, "w") as f_fa, open(args.tx2gene, "w") as f_map:
        for _, row in df.iterrows():
            locus = str(row["locus"])
            seq = str(row["cds_seq"]).replace(" ", "").replace("\n", "").upper()
            if seq in ("", "nan", "NA"): # safeguards
                continue
            # transcript_id must be first token on the header line
            f_fa.write(f">{locus}\n{seq}\n")
            # transcript_id \t gene_id
            f_map.write(f"{locus}\t{locus}\n")

if __name__ == "__main__":
    main()
