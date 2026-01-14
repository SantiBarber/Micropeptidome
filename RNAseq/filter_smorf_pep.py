#!/usr/bin/env python3

import argparse
from Bio import SeqIO

print(r"""
          _____           _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  almoraco 03/11/2025
""")

print(
    "filters TransDecoder peptides by length\n"
    "retaining only those within the smORF size range\n"
)

def main():
    parser = argparse.ArgumentParser(description="Filter TransDecoder peptide FASTA for small ORFs")
    parser.add_argument("pep_file", help="TransDecoder peptide FASTA file")
    parser.add_argument("--min_len", type=int, default=10, help="Minimum AA length")
    parser.add_argument("--max_len", type=int, default=150, help="Maximum AA length")
    parser.add_argument("--out_fasta", default="smorfs.fa", help="Output FASTA for smORFs")
    parser.add_argument("--out_ids", default="smorf_ids.txt", help="Output file listing retained smORF IDs")

    args = parser.parse_args()

    kept = []
    with open(args.out_fasta, "w") as fasta_out, open(args.out_ids, "w") as id_out:
        for record in SeqIO.parse(args.pep_file, "fasta"):
            length = len(record.seq)

            if args.min_len <= length <= args.max_len:
                kept.append(record.id)
                fasta_out.write(f">{record.id}\n{record.seq}\n")
                id_out.write(record.id + "\n")

    print(f"Kept {len(kept)} smORFs (AA length between {args.min_len} and {args.max_len})")

if __name__ == "__main__":
    main()
