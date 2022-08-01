configfile: "config/config.yaml"

import os

filter_seqkit = config["config_seqkit"]["cut"]
summary_presence = config["summary_file"]["presence"]

BWA_INDEX = ['amb','ann','bwt','pac','sa']

work_directory = config["DIRECTORIES"]["working_directory"]
res_directory = config["DIRECTORIES"]["rvhaplo_results"]
reference_file = config["FILES"]["reference_file"]
all_refs = config["FILES"]["all_refs_file"]

path_fastq = config["DIRECTORIES"]["fastq_file"]
fastq_files = []
for file in os.listdir(path_fastq):
    if file.endswith('.fastq'):
        fastq_files.append(os.path.splitext(os.path.basename(file))[0])

path_summary = config["DIRECTORIES"]["sequencing_summary"]
if summary_presence == True:
    summary_files = []
    for file in os.listdir(path_summary):
        if file.endswith('.txt'):
            summary_files.append(os.path.splitext(os.path.basename(file))[0])

ref_without_ext = os.path.splitext(os.path.basename(reference_file))[0]
basename_ref = os.path.basename(reference_file)
path_ref = os.path.dirname(reference_file)


if summary_presence == True:
    rule all:
        input:
            expand(f"{path_summary}/{{summary}}.html",summary=summary_files),
            haplo_refs=expand(f"{res_directory}/{{fastq}}/{{fastq}}_sup{filter_seqkit}/haplo_refs.fasta",fastq=fastq_files) if filter_seqkit != 0
            else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/haplo_refs.fasta",fastq=fastq_files)
else:
    rule all:
        input:
            haplo_refs=expand(f"{res_directory}/{{fastq}}/{{fastq}}_sup{filter_seqkit}/haplo_refs.fasta",fastq=fastq_files) if filter_seqkit != 0
            else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/haplo_refs.fasta",fastq=fastq_files)


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
        else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/rvhaplo_haplotypes.fasta", fastq=fastq_files)
    conda:
        "envs/rvhaplo.yml"
    shell:
        "./rvhaplo.sh -i {input.sam} -r {input.reference} -o {res_directory}/{fastq_files}/{fastq_files}_sup{filter_seqkit} -t 32 || true" if filter_seqkit != 0
        else "./rvhaplo.sh -i {input.sam} -r {input.reference} -o {res_directory}/{fastq_files}/{fastq_files}_allreads -t 32 || true"


rule merge:
    input:
        haplotypes = rules.RVHaplo.output.haplo,
        all_refs_file = all_refs
    output:
        haplo_refs = expand(f"{res_directory}/{{fastq}}/{{fastq}}_sup{filter_seqkit}/haplo_refs.fasta", fastq=fastq_files) if filter_seqkit != 0
        else expand(f"{res_directory}/{{fastq}}/{{fastq}}_allreads/haplo_refs.fasta", fastq=fastq_files)
    run:
        with open(input.all_refs_file) as file:
            data = file.read()
        for haplo_file in input.haplotypes:
            with open(haplo_file) as file:
                data2 = file.read()
        data2 += data
        for merge_file in output.haplo_refs:
            with open(merge_file, "w") as file:
                file.write(data2)