rule download_tg_list:
    output: config['ref']['tg_list']
    log: "logs/download_tg_list.log"
    params: url = config['ref']['tg_list_url']
    message: "Executing {rule}: Downloading a list of loci with disease-associated repeats to {output}."
    shell: "(wget -qO - {params.url} > {output}) > {log} 2>&1"


rule generate_tg_bed:
    input:
        tg_list = config['ref']['tg_list'],
        fai = config['ref']['index']
    output: config['ref']['tg_bed']
    log: "logs/generate_tg_bed.log"
    conda: "envs/bedtools.yaml"
    params: slop = 1000
    message: "Executing {rule}: Adding {params.slop}bp slop to {input.tg_list} to generate {output}."
    shell: 
        """
        (grep -v '^#' {input.tg_list} | sort -k1,1V -k2,2g \
        | bedtools slop -b {params.slop} -g {input.fai} -i - \
        > {output}) > {log} 2>&1
        """


rule generate_last_index:
    input: config['ref']['fasta']
    output:
        [f"{config['ref']['last_index']}.{suffix}"
         for suffix in ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']]
    log: "logs/generate_last_index.log"
    conda: "envs/last.yaml"
    params: "-uRY32 -R01"
    threads: 24
    message: "Executing {rule}: Generating last index of {input} using params: {params}"
    shell: f"(lastdb -P{{threads}} {{params}} {config['ref']['last_index']} {{input}}) > {{log}} 2>&1"


rule last_align:
    input:
        [f"{config['ref']['last_index']}.{suffix}"
         for suffix in ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']],
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bed = config['ref']['tg_bed'],
        score_matrix = config['score_matrix']
    output: temp(f"samples/{sample}/tandem-genotypes/{sample}.maf.gz")
    log: f"samples/{sample}/logs/last/align/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/last/align/{sample}.tsv"
    conda: "envs/last.yaml"
    params: 
        last_index = config['ref']['last_index'],
        extra = "-C2"
    threads: 24
    message: "Executing {rule}: Aligning {input.bed} regions of {input.bam} to {params.last_index} using lastal with {input.score_matrix} score matrix."
    shell:
        """
        (samtools view -@3 -bL {input.bed} {input.bam} | samtools fasta \
         | lastal -P20 -p {input.score_matrix} {params.extra} {params.last_index} - \
         | last-split | bgzip > {output}) > {log} 2>&1
        """


rule tandem_genotypes:
    input:
        maf = f"samples/{sample}/tandem-genotypes/{sample}.maf.gz",
        repeats = config['ref']['tg_bed']
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
    log: f"samples/{sample}/logs/tandem-genotypes/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.tsv"
    conda: "envs/tandem-genotypes.yaml"
    message: "Executing {rule}: Genotyping tandem repeats from {input.repeats} regions in {input.maf}."
    shell: "(tandem-genotypes {input.repeats} {input.maf} > {output}) > {log} 2>&1"


rule tandem_genotypes_plot:
    input: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.txt"
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.pdf"
    log: f"samples/{sample}/logs/tandem-genotypes/plot/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/plot/{sample}.tsv"
    conda: "envs/tandem-genotypes.yaml"
    params: top_N_plots = 100
    message: "Executing {rule}: Plotting tandem repeats from {input}."
    shell: "(tandem-genotypes-plot -n {params.top_N_plots} {input} {output}) > {log} 2>&1"


rule tandem_repeat_coverage_dropouts:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bed = config['ref']['tg_bed']
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.dropouts.txt"
    log: f"samples/{sample}/logs/tandem-genotypes/{sample}.dropouts.log"
    benchmark: f"samples/{sample}/benchmarks/tandem-genotypes/{sample}.dropouts.tsv"
    conda: "envs/tandem-genotypes.yaml"
    message: "Executing {rule}: Identify coverage dropouts in {input.bed} regions in {input.bam}."
    shell: "(python3 workflow/scripts/check_tandem_repeat_coverage.py {input.bed} {input.bam} > {output}) > {log} 2>&1"
