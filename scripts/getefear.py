import pandas as pd

df = pd.read_csv("data/microproteinasSEER.csv", sep="\t", encoding="utf-8-sig")

# Limpiar nombres de columnas (quitar espacios extra)
df.columns = df.columns.str.strip()

gtf = pd.DataFrame({
    "seqname": df["sseqid"],
    "source": "microproteomics",
    "feature": "microprotein_peptide",
    "start": df["sstart"],
    "end": df["send"],
    "score": df["bitscore"],
    "strand": "+",
    "frame": ".",
    "attribute": (
        'gene_id "' + df["sseqid"].str.extract(r'\|([^|]+)_')[0] + '"; '
        'peptide "' + df["microprotein-derived peptide"] + '"; '
        'qseq "' + df["qseq"] + '"; '
        'sseq "' + df["sseq"] + '";'
    )
})

gtf.to_csv("data/microproteins.gtf", sep="\t", index=False, header=False)