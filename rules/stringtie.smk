rule stringtie_assemble:
    input:
        bam=bam_path,
        ref_gtf=config["gencode_gtf"]
    output:
        gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf"
    threads: config.get("threads_stringtie", 8)
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/stringtie"

        stringtie "{input.bam}" \
          -G "{input.ref_gtf}" \
          -o "{output.gtf}" \
          -p {threads}
        """
