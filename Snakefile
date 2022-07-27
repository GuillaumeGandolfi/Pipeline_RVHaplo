configfile: "config/config.yaml"

import os

BWA_INDEX = ['amb','ann','bwt','pac','sa']
reference_file = config["FILES"]["reference_file"]
fastq_file = config["FILES"]["fastq_file"]
work_directory = config["DIRECTORIES"]["working_directory"]
seq_summary = config["FILES"]["sequencing_summary"]

ref_without_ext = os.path.splitext(os.path.basename(reference_file))[0]
fastq_without_ext = os.path.splitext(os.path.basename(fastq_file))[0]
summary_without_ext = os.path.splitext(os.path.basename(seq_summary))[0]

basename_ref = os.path.basename(reference_file)

rule all:
    input:
        expand(f"{work_directory}/sequencing/{{summary}}.html",summary=summary_without_ext),
        expand(f"{work_directory}/mapped/{{fastq}}_map{{reference}}.sam", fastq=fastq_without_ext, reference=ref_without_ext)


rule pycoQC:
    input:
        summary = seq_summary
    output:
        html = f"{work_directory}/sequencing/{{summary_without_ext}}.html"
    conda:
        "envs/pycoQC.yml"
    shell:
        "pycoQC --summary_file {input} --html_outfile {output}"


rule convert_fasta:
    input:
        fastq = fastq_file
    output:
        conv_fasta = f"{work_directory}/reads/{{fastq_without_ext}}.fasta"
    envmodules:
        "bioinfo/seqtk/1.3-r106"
    shell:
        "seqtk seq -a {input} > {output}"

rule bwa_index:
    input:
        reference = reference_file
    output:
        index_file = expand(f"{work_directory}/ref/{basename_ref}.{{suffix}}", suffix=BWA_INDEX)
    envmodules:
        "bioinfo/bwa/0.7.17"
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        reference = rules.bwa_index.input.reference,
        index_file = rules.bwa_index.output.index_file,
        fasta = rules.convert_fasta.output.conv_fasta
    output:
        sam_file = f"{work_directory}/mapped/{{fastq_without_ext}}_map{{ref_without_ext}}.sam"
    envmodules:
        "bioinfo/bwa/0.7.17"
    shell:
        "bwa mem {input.reference} {input.fasta} > {output}"