import glob
import os

configfile: 'variables_config.yaml'

SAMPLES, = glob_wildcards(config["input_files"])
path_to = lambda file: os.path.join(config["output_dir"], file)

rule all:
    input:
        expand(config["output_dir"] + "{number}.teinsertions", number=SAMPLES)


rule sort_bam:
    input:
        config["input_files"]
    output:
        temp(path_to("{number}_sorted.bam"))
    threads: 16
    shell:
        "samtools sort -@ {threads} -T {wildcards.number}_tmp -n -o {output} {input}"

rule make_bwa_index:
    input:
        config["reference_genome"]
    output:
        expand(config["reference_genome"] + ".{ext}", ext=['sa', 'pac', 'bwt', 'ann', 'amb'])
        # touch('index.done')
    shell:
        """
        bwa index {input}
        """

rule make_fastq_files:
    input:
        rules.sort_bam.output
    output:
        fastq_1 = path_to("{number}_1.fastq.gz"),
        fastq_2 = path_to("{number}_2.fastq.gz")
    threads: 16
    shell:
        "samtools fastq -@ {threads} -n -1 {output.fastq_1} -2 {output.fastq_2} -0 /dev/null -s /dev/null {input}"

rule bwa_mem:
    input:
        idx_done = rules.make_bwa_index.output,
        fastq = path_to("{number}_{read_num}.fastq.gz")
    output:
        path_to("{number}_{read_num}.bam")
    threads: 16
    params:
        ref_genome = config['reference_genome']
    shell:
        """
        bwa mem -t {threads} {params.ref_genome} {input.fastq} | samtools view -b -o {output} -
        """

rule se2pe:
    input:
        fastq_1 = rules.make_fastq_files.output.fastq_1,
        fastq_2 = rules.make_fastq_files.output.fastq_2,
        bam_1 = path_to("{number}_1.bam"),
        bam_2 = path_to("{number}_2.bam")
    output:
        path_to("{number}_paired_aligned.bam")
    shell:
        """
        popte2 se2pe --fastq1 {input.fastq_1} --fastq2 {input.fastq_2} \
        --bam1 {input.bam_1} --bam2 {input.bam_1} --sort \
        --output {output}
        """

rule ppileup:
    input:
        rules.se2pe.output
    output:
        path_to("{number}.ppileup.gz")
    params:
        hier = config['hierarchy_file']
    shell:
        "popte2 ppileup --bam {input} --hier {params.hier} --output {output}"

rule identifySignatures:
    input:
        rules.ppileup.output
    output:
        temp(path_to("{number}.signatures"))
    params:
        mode="separate",
        min_count=3
    shell:
        """
        popte2 identifySignatures --ppileup {input} \
        --mode {params.mode} --min-count {params.min_count} \
        --output {output}
        """

rule freq:
    input:
        signature = rules.identifySignatures.output,
        ppileup = rules.ppileup.output
    output:
        path_to("{number}.freqsig")
    shell:
        """
        popte2 frequency --ppileup {input.ppileup} \
        --signature {input.signature} \
        --output {output}
        """

rule filter_freq:
    input:
        rules.freq.output
    output:
        path_to("{number}.filtered.freqsig")
    params:
        max_otherte_count=0,
        max_structvar_count=0
    shell:
        """
        popte2 filterSignatures \
        --max-otherte-count {params.max_otherte_count} \
        --max-structvar-count {params.max_structvar_count} \
        --input {input} \
        --output {output}
        """

rule pairup_filtered_signatures:
    input:
        rules.filter_freq.output
    output:
        path_to("{number}.teinsertions")
    params:
        ref_genome = config['reference_genome'],
        hier = config['hierarchy_file']
    shell:
        "popte2 pairupSignatures --signature {input} --ref-genome {params.ref_genome} --hier {params.hier} --output {output}"
