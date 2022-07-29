configfile: "config/config.yaml"

import os

BWA_INDEX = ['amb','ann','bwt','pac','sa']
reference_file = config["FILES"]["reference_file"]
path_fastq = config["DIRECTORIES"]["fastq_file"]
work_directory = config["DIRECTORIES"]["working_directory"]
path_summary = config["DIRECTORIES"]["sequencing_summary"]
res_directory = config["DIRECTORIES"]["rvhaplo_results"]

filter_seqkit = config["config_seqkit"]["cut"]

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
        expand(f"{res_directory}/{{fastq}}/{{fastq}}_sup{filter_seqkit}/rvhaplo_haplotypes.fasta", fastq=fastq_files) if filter_seqkit !=0
        else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/rvhaplo_haplotypes.fasta", fastq=fastq_files)


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

rule cut_fasta:
    input:
        fasta = rules.convert_fasta.output.conv_fasta
    output:
        fasta_cut = expand(f"{path_fastq}/{{fastq}}_sup{filter_seqkit}.fasta", fastq=fastq_files)
    envmodules:
        "bioinfo/seqkit/2.1.0"
    params:
        filter = config["config_seqkit"]["cut"]
    shell:
        "seqkit seq -g -m {params.filter} {input} > {output}"

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
        fasta = rules.cut_fasta.output.fasta_cut if filter_seqkit != 0 else rules.convert_fasta.output.conv_fasta
    output:
        sam_file = expand(f"{work_directory}/mapped/{{fastq}}_sup{filter_seqkit}_map{{reference}}.sam", fastq=fastq_files, reference=ref_without_ext)
        if filter_seqkit != 0 else expand(f"{work_directory}/mapped/{{fastq}}_map{{reference}}.sam", fastq=fastq_files, reference=ref_without_ext)
    envmodules:
        "bioinfo/bwa/0.7.17"
    shell:
        "bwa mem {input.reference} {input.fasta} > {output}"


rule RVHaplo:
    input:
        sam = rules.bwa_mem.output.sam_file,
        reference = reference_file
    output:
        haplo = expand(f"{res_directory}/{{fastq}}/{{fastq}}_sup{filter_seqkit}/rvhaplo_haplotypes.fasta", fastq=fastq_files) if filter_seqkit != 0
        else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/rvhaplo_haplotypes.fasta")
    conda:
        "envs/rvhaplo.yml"
    shell:
        "./rvhaplo.sh -i {input.sam} -r {input.reference} -o {res_directory}/{fastq_files}/{fastq_files}_sup{filter_seqkit} -t 32 || true" if filter_seqkit != 0
        else "./rvhaplo.sh -i {input.sam} -r {input.reference} -o {res_directory}/{fastq_files}/{fastq_files}_allreads -t 32 || true"