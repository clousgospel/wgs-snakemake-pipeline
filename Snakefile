SAMPLES = ["sample1"]
REF = "ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

rule all:
    input:
        expand("results/{sample}/final.vcf.gz", sample=SAMPLES),
        "qc/multiqc/multiqc_report.html"


rule bwa_index:
    input: REF
    output:
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt",
        REF + ".pac",
        REF + ".sa"
    shell: "bwa index {input}"

rule samtools_index:
    input: REF
    output: REF + ".fai"
    shell: "samtools faidx {input}"

rule gatk_dict:
    input: REF
    output: "ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"


rule fastqc:
    input:
        R1 = "data/{sample}_R1.fastq.gz",
        R2 = "data/{sample}_R2.fastq.gz"
    output:
        "qc/fastqc/{sample}_R1_fastqc.html",
        "qc/fastqc/{sample}_R2_fastqc.html"
    shell: "fastqc {input.R1} {input.R2} -o qc/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=['R1', 'R2']),
        expand("qc/metrics/{sample}_mark_dup_metrics.txt", sample=SAMPLES)
    output: "qc/multiqc/multiqc_report.html"
    shell: "multiqc qc/fastqc/ qc/metrics/ -o qc/multiqc/"


rule align_bwa:
    input:
        R1 = "data/{sample}_R1.fastq.gz",
        R2 = "data/{sample}_R2.fastq.gz",
        ref = REF,
        idx = [REF + ".amb", REF + ".bwt", REF + ".pac", REF + ".fai", "ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict"]
    output: temp("results/{sample}/aligned.bam")
    threads: 4
    shell:
        "bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | samtools view -Sb - > {output}"

rule sort_bam:
    input: "results/{sample}/aligned.bam"
    output: temp("results/{sample}/sorted.bam")
    shell: "samtools sort -@ 4 {input} -o {output}"

rule add_read_groups:
    input: "results/{sample}/sorted.bam"
    output: temp("results/{sample}/rg_added.bam")
    params:
        rg_id = "{sample}",
        rg_sm = "{sample}",
        rg_pl = "ILLUMINA",
        rg_lb = "Lib1"
    shell:
        """
        gatk --java-options "-Xmx4g" AddOrReplaceReadGroups \
        -I {input} -O {output} \
        --RGLB {params.rg_lb} --RGPL {params.rg_pl} \
        --RGPU {params.rg_id} --RGSM {params.rg_sm}
        """

rule mark_duplicates:
    input: "results/{sample}/rg_added.bam"
    output:
        bam = "results/{sample}/dedup.bam",
        metrics = "qc/metrics/{sample}_mark_dup_metrics.txt"
    shell:
        """
        gatk --java-options "-Xmx4g" MarkDuplicates \
        -I {input} -O {output.bam} -M {output.metrics} \
        --CREATE_INDEX true
        """


rule variant_call_gatk:
    input:
        bam = "results/{sample}/dedup.bam",
        ref = REF,
        fai = REF + ".fai",
        dct = "ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
    output: "results/{sample}/raw.vcf.gz"
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R {input.ref} -I {input.bam} -O {output}
        """

rule filter_vcf:
    input: "results/{sample}/raw.vcf.gz"
    output: "results/{sample}/final.vcf.gz"
    shell:
        """
        bcftools filter -O z -o {output} {input} && \
        tabix -p vcf {output}
        """