configfile: "config/config.yaml"

import os

BWA_INDEX = ['amb','ann','bwt','pac','sa']
reference_file = config["FILES"]["reference_file"]
path_fastq = config["DIRECTORIES"]["fastq_file"]
work_directory = config["DIRECTORIES"]["working_directory"]
path_summary = config["DIRECTORIES"]["sequencing_summary"]
res_directory = config["DIRECTORIES"]["rvhaplo_results"]

ref_without_ext = os.path.splitext(os.path.basename(reference_file))[0]
basename_ref = os.path.basename(reference_file)
path_ref = os.path.dirname(reference_file)

summary_files = []
for file in os.listdir(path_summary):
    if file.endswith('.txt'):
        summary_files.append(os.path.splitext(os.path.basename(file))[0])


fastq_files = []
for file in os.listdir(path_fastq):
    if file.endswith('.fastq'):
        fastq_files.append(os.path.splitext(os.path.basename(file))[0])

rule all:
    input:
        expand(f"{path_summary}/{{summary}}.html",summary=summary_files),
        expand(f"{work_directory}/mapped/{{fastq}}_map{{reference}}.sam", fastq=fastq_files, reference=ref_without_ext)


rule pycoQC:
    input:
        summary = f"{path_summary}/{{summary}}.txt"
    output:
        html = f"{path_summary}/{{summary}}.html"
    conda:
        "envs/pycoQC.yml"
    shell:
        "pycoQC --summary_file {input} --html_outfile {output}"


rule convert_fasta:
    input:
        fastq = expand(f"{path_fastq}/{{fastq}}.fastq", fastq=fastq_files)
    output:
        conv_fasta = expand(f"{path_fastq}/{{fastq}}.fasta", fastq=fastq_files)
    envmodules:
        "bioinfo/seqtk/1.3-r106"
    shell:
        "seqtk seq -a {input} > {output}"

rule bwa_index:
    input:
        reference = reference_file
    output:
        index_file = expand(f"{path_ref}/{basename_ref}.{{suffix}}", suffix=BWA_INDEX)
    envmodules:
        "bioinfo/bwa/0.7.17"
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        reference = reference_file,
        index_file = rules.bwa_index.output.index_file,
        fasta = rules.convert_fasta.output.conv_fasta
    output:
        sam_file = expand(f"{work_directory}/mapped/{{fastq}}_map{{reference}}.sam", fastq=fastq_files, reference=ref_without_ext)
    envmodules:
        "bioinfo/bwa/0.7.17"
    shell:
        "bwa mem {input.reference} {input.fasta} > {output}"