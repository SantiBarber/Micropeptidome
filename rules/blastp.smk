rule make_human_proteome_blastdb:
    input:
        fa=HUMAN_PROTEOME_FA
    output:
        done=f"{OUTDIR}/blastdb/human_proteome.db.done"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=240
    params:
        db_prefix=HUMAN_DB_PREFIX
    conda:
        "../envs/BlastP.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{params.db_prefix}")"

        makeblastdb \
          -in "{input.fa}" \
          -dbtype prot \
          -parse_seqids \
          -out "{params.db_prefix}"

        touch "{output.done}"
        """

rule blastp_human_homology_locus_summary:
    input:
        db_done=f"{OUTDIR}/blastdb/human_proteome.db.done",
        loci_csv=f"{COHORT_PREFIX}.all_loci.with_tpms.csv",
        script=config["blastp_append_script"]
    output:
        out_csv=f"{COHORT_PREFIX}.all_loci.with_tpms.blastp_human.csv"
    threads: 8
    resources:
        mem_mb=24000,
        runtime=240
    params:
        db_prefix=HUMAN_DB_PREFIX,
        evalue=BLAST_EVALUE
    conda:
        "../envs/BlastP.yaml"
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" \
          --in_csv "{input.loci_csv}" \
          --out_csv "{output.out_csv}" \
          --db "{params.db_prefix}" \
          --evalue {params.evalue} \
          --threads {threads}
        """

rule blastp_human_homology_shared_summary:
    input:
        db_done=f"{OUTDIR}/blastdb/human_proteome.db.done",
        shared_csv=f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.csv",
        script=config["blastp_append_script"]
    output:
        out_csv=f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.blastp_human.csv"
    threads: 8
    resources:
        mem_mb=24000,
        runtime=240
    params:
        db_prefix=HUMAN_DB_PREFIX,
        evalue=BLAST_EVALUE
    conda:
        "../envs/BlastP.yaml"
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" \
          --in_csv "{input.shared_csv}" \
          --out_csv "{output.out_csv}" \
          --db "{params.db_prefix}" \
          --evalue {params.evalue} \
          --threads {threads}
        """
