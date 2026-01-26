#############################################
# Snakefile: per-sample ShortStop workflow
# Uses: StringTie, TransDecoder, ShortStop, 
# BlastP, Dr. Brendan Miller's smORF annotator.
#############################################

import csv
from pathlib import Path

configfile: "config.yaml"

# Absolute outdir avoids SLURM working-dir surprises
OUTDIR = str(Path(config["outdir"]).resolve())

# STAR: user may point this at a prebuilt STAR genome index directory.
STAR_INDEX_DIR = str(Path(config.get("star_index_dir", f"{OUTDIR}/STAR_index")).resolve())

# STAR output BAM path
def star_bam(sample: str) -> str:
    return f"{OUTDIR}/star/{sample}.aligned.bam"

def bam_path(wc):
    return star_bam(wc.sample)

RESULTS_SHORTSTOP_DIR = f"{OUTDIR}/results_shortstop"
MERGED_DIR = config.get("merged_dir", f"{OUTDIR}/merged_per_sample")
COHORT_PREFIX = config["cohort_prefix"]
MIN_PATIENTS = int(config.get("min_patients", 2))

UNITS_CSV = config["units_csv"]

RSEM_DIR = config.get("rsem_dir", f"{OUTDIR}/results_rsem_smorf")
RSEM_REF_DIR = f"{RSEM_DIR}/reference"
RSEM_REF_PREFIX = f"{RSEM_REF_DIR}/smorfs"
RSEM_STRANDEDNESS = config.get("rsem_strandedness", "none")
MAKE_SMORF_RSEM_REF_SCRIPT = config["make_smorf_rsem_ref_script"]
ADD_RSEM_TPMS_SCRIPT = config["add_rsem_tpms_script"]

HUMAN_PROTEOME_FA = config["human_proteome_fa"]
HUMAN_DB_PREFIX = config.get("human_blastdb_prefix", f"{OUTDIR}/blastdb/human_proteome")
BLAST_EVALUE = float(config.get("blastp_evalue", 1e-3))


# Build mapping: sample -> (r1, r2)
UNITS = {}
with open(UNITS_CSV, newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        sample = (row.get("sample") or row.get("name") or "").strip()
        r1 = (row.get("r1") or row.get("fastq_r1") or row.get("read1") or "").strip()
        r2 = (row.get("r2") or row.get("fastq_r2") or row.get("read2") or "").strip()
        if not sample:
            continue
        if not r1 or not r2:
            raise ValueError(f"Missing r1/r2 for sample '{sample}' in {UNITS_CSV}")
        if sample in UNITS:
            raise ValueError(f"Duplicate sample '{sample}' in {UNITS_CSV}")
        UNITS[sample] = (r1, r2)

SAMPLES = sorted(UNITS.keys())
if not SAMPLES:
    raise ValueError(f"No samples found in {UNITS_CSV}. Check units.csv.")

def fastq_r1(wc):
    try:
        return UNITS[wc.sample][0]
    except KeyError as e:
        raise ValueError(f"No FASTQ entry for sample '{wc.sample}' in {UNITS_CSV}") from e

def fastq_r2(wc):
    try:
        return UNITS[wc.sample][1]
    except KeyError as e:
        raise ValueError(f"No FASTQ entry for sample '{wc.sample}' in {UNITS_CSV}") from e

include: "rules/star_align.smk"
include: "rules/stringtie.smk"
include: "rules/transdecoder.smk"
include: "rules/smorfs.smk"
include: "rules/shortstop.smk"
include: "rules/annotator.smk"
include: "rules/merge.smk"
include: "rules/rsem.smk"
include: "rules/blastp.smk"


### --------------------------------------------------oOo------------------------------------------------- ###

rule all:
    input:
        expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/predict.done", sample=SAMPLES),
        expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.gtf", sample=SAMPLES),
        expand(f"{MERGED_DIR}/{{sample}}.merged.csv", sample=SAMPLES),
        f"{COHORT_PREFIX}.all_loci.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.csv",
        f"{COHORT_PREFIX}.all_loci.with_tpms.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.csv",
        f"{COHORT_PREFIX}.all_loci.with_tpms.blastp_human.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.blastp_human.csv"
