#!/usr/bin/env python3
import argparse


def parse_gtf_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, val = part.split(" ", 1)
        attrs[key] = val.strip().strip('"')
    return attrs


def load_annotations(path: str) -> dict:
    ann = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene_id = parts[0]
            smorf_type = parts[1]
            ann[gene_id] = smorf_type
    return ann


def main():
    ap = argparse.ArgumentParser(description="Append smORF type as 10th column in a GTF.")
    ap.add_argument("--gtf", required=True, help="Input GTF file")
    ap.add_argument("--annotations", required=True, help="Annotator output file (gene_id\\tannotation\\t...)")
    ap.add_argument("--out", required=True, help="Output annotated GTF")
    args = ap.parse_args()

    ann = load_annotations(args.annotations)

    with open(args.gtf, "r") as infile, open(args.out, "w") as out:
        for line in infile:
            if not line.strip():
                continue
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                out.write(line)
                continue
            attrs = parse_gtf_attrs(parts[8])
            gene_id = attrs.get("gene_id")
            smorf_type = ann.get(gene_id, "NA")

            attr_str = parts[8].strip()
            if attr_str and not attr_str.endswith(";"):
                attr_str += ";"
            attr_str += f' smorf_type "{smorf_type}";'
            parts[8] = attr_str.strip()
            out.write("\t".join(parts) + "\n")


if __name__ == "__main__":
    main()
