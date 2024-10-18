
import os
import pandas as pd

sample_reads = pd.read_csv("./samples.csv")
samples = sample_reads["sample_id"]

star_index = "/home/jverwilt/resources/indexes/chm13/STAR_index"
hisat_index = "/home/jverwilt/resources/indexes/chm13/HISAT_index/genome"
bwa_index = "/home/jverwilt/resources/indexes/chm13/bwa_index/"
gtf = "/home/jverwilt/resources/gtf/chm13/chm13_MT_spikes.gtf"
bed = "/home/jverwilt/resources/gtf/chm13/chm13_MT_spikes.bed"
path = "/home/jverwilt/JV2408_tear_RNA/output/20240906_AV242402_4843/"
config_circ = "/home/jverwilt/JV2408_tear_RNA/code/config.yml"
prepare_deconv = "/home/jverwilt/JV2408_tear_RNA/code/chm13_MT_prepare_deconv.tsv"
matrix_deconv = ""

rule all: 
    input:
        os.path.join(path, "multiqc_report.html")

rule combine:
    input:
        R1L1=lambda wc: sample_reads[sample_reads["sample_id"] == wc.sample]["R1L1"].item(),
        R1L2=lambda wc: sample_reads[sample_reads["sample_id"] == wc.sample]["R1L2"].item(),
        R2L1=lambda wc: sample_reads[sample_reads["sample_id"] == wc.sample]["R2L1"].item(),
        R2L2=lambda wc: sample_reads[sample_reads["sample_id"] == wc.sample]["R2L2"].item()
    output:
        R1=os.path.join(path, "{sample}/01_combine/{sample}_R1.fastq.gz"),
        R2=os.path.join(path, "{sample}/01_combine/{sample}_R2.fastq.gz")
    run:
        shell("cat {input.R1L1} {input.R1L2} > {output.R1}")
        shell("cat {input.R2L1} {input.R2L2} > {output.R2}")

rule extract_umi:
    input:
        R1=os.path.join(path, "{sample}/01_combine/{sample}_R1.fastq.gz"),
        R2=os.path.join(path, "{sample}/01_combine/{sample}_R2.fastq.gz")
    output:
        R1=os.path.join(path, "{sample}/02_extract_umi/{sample}_umi_R1.fastq.gz"),
        R2=os.path.join(path, "{sample}/02_extract_umi/{sample}_umi_R2.fastq.gz"),
        umilog=os.path.join(path, "{sample}/02_extract_umi/UMIextract.log")
    conda:
        "base"
    shell:
        "umi_tools extract --extract-method=string -I {input.R2} -L {output.umilog} --bc-pattern=NNNNNNNN --read2-in={input.R1} -S {output.R2} --read2-out={output.R1}"


rule map:
    input:
        R1=os.path.join(path, "{sample}/02_extract_umi/{sample}_umi_R1.fastq.gz"),
        R2=os.path.join(path, "{sample}/02_extract_umi/{sample}_umi_R2.fastq.gz")
    params:
        star_index = star_index,
        path = path
    output:
        os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam")
    shell:
        "STAR --chimSegmentMin 10 --runThreadN 2 --genomeDir {params.star_index} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.path}{wildcards.sample}/03_map/{wildcards.sample}. --readFilesIn {input.R1} {input.R2}"

rule index:
    input:
        os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam.bai")
    shell:
        "samtools index {input}"

rule dedup:
    input:
        bam=os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam"),
        bami=os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam.bai")
    output:
        bam=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam"),
        log=os.path.join(path, "{sample}/04_dedup/umidedup.txt")
    params:
        path = path
    conda:
        "base"
    shell:
        "umi_tools dedup -I {input.bam} --output-stats={params.path}{wildcards.sample}/04_dedup/deduplicated --paired -S {output.bam} -L {output.log}"

rule count:
    input:
        os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam")
    params:
        gtf = gtf
    output:
        os.path.join(path, "{sample}/05_count/{sample}_htseq_dedup_counts.txt")
    shell:
        "htseq-count --format bam --order pos --nonunique none --stranded reverse {input} {params.gtf} > {output}"

    #split bam file into multiple smaller files and then process, change CPU number

rule bamtofastq:
    input:
        os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam")
    output:
        R1=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup_R1.fastq"),
        R2=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup_R2.fastq")
    shell:
        "samtools sort -n {input} | bedtools bamtofastq -i stdin -fq {output.R1} -fq2 {output.R2}"

rule ciriquant:
    input:
        R1=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup_R1.fastq"),
        R2=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup_R2.fastq")
    output:
        os.path.join(path, "{sample}/06_circRNA/{sample}_circRNA.gtf")
    params:
        config_circ = config_circ,
        path = path
    conda:
        "/home/jverwilt/resources/tools/CIRIquant/CIRIquant_env"
    shell:
        "CIRIquant --config {params.config_circ} -l 2 -t 2 -1 {input.R1} -2 {input.R2} -o {params.path}{wildcards.sample}/06_circRNA/ -p {wildcards.sample}_circRNA"

rule prepare_deconvolution:
    input:
        os.path.join(path, "{sample}/05_count/{sample}_htseq_dedup_counts.txt")
    output:
        os.path.join(path, "{sample}/07_deconv/{sample}_todeconvolve.csv")
    params:
        prepare_deconv = prepare_deconv
    shell:
        "Rscript /home/jverwilt/JV2408_tear_RNA/code/prepare_deconvolution.R {input} {params.prepare_deconv} {wildcards.sample} {output}"


rule deconvolution:
    input:
        os.path.join(path, "{sample}/07_deconv/{sample}_todeconvolve.csv")
    output:
        os.path.join(path, "{sample}/07_deconv/{sample}_deconvolutionCoefs.csv")
    params:
        path = path
    shell:
        """
        cd /home/jverwilt/resources/tools/deconvolution/deconvolve_cfrna_tutorial
        python3 -c 'import deconvolve as deconv; deconv.main(1, "{input}", ["{wildcards.sample}"] , "nuSVR", "{wildcards.sample}", "{params.path}{wildcards.sample}/07_deconv/", jackknife = False)'
        """

rule rseqc:
    input:
        dedup=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam"),
        raw=os.path.join(path, "{sample}/03_map/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        bam_stat_dedup=os.path.join(path, "{sample}/08_rseqc/{sample}_dedup_bam_stat.txt"),
        bam_stat_raw=os.path.join(path, "{sample}/08_rseqc/{sample}_raw_bam_stat.txt"),
        strandedness=os.path.join(path, "{sample}/08_rseqc/{sample}_strandedness.txt")
    params:
        bed = bed
    shell:
        """
        bam_stat.py -i {input.dedup} > {output.bam_stat_dedup}
        bam_stat.py -i {input.raw} > {output.bam_stat_raw}
        infer_experiment.py -r {params.bed} -i {input.dedup} > {output.strandedness}
        """

rule fastqc:
    input:
        R1=os.path.join(path, "{sample}/01_combine/{sample}_R1.fastq.gz"),
        R2=os.path.join(path, "{sample}/01_combine/{sample}_R2.fastq.gz"),
        bam=os.path.join(path, "{sample}/04_dedup/{sample}.Aligned.sortedByCoord.out.dedup.bam")
    output:
        os.path.join(path, "{sample}/08_fastqc/")
    shell:
        "fastqc -o {output} {input.R1} {input.R2} {input.bam}"

rule multiqc:
    input:
        path
    output:
        os.path.join(path, "multiqc_report.html")
    shell:
        "multiqc {input} -d --filename stdout > {output}"