localrules: bgzip_vcf, tabix_vcf, tabix_bcf, create_ped


rule bgzip_vcf:
    input: f"cohorts/{cohort}/{{prefix}}.vcf"
    output: f"cohorts/{cohort}/{{prefix}}.vcf.gz"
    log: f"cohorts/{cohort}/logs/bgzip/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bgzip/{{prefix}}.tsv"
    threads: 2
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: f"cohorts/{cohort}/{{prefix}}.vcf.gz"
    output: f"cohorts/{cohort}/{{prefix}}.vcf.gz.tbi"
    log: f"cohorts/{cohort}/logs/tabix/index/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/index/{{prefix}}.tsv"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(tabix {params} {input}) > {log} 2>&1"


rule tabix_bcf:
    input: f"cohorts/{cohort}/{{prefix}}.bcf"
    output: temp(f"cohorts/{cohort}/{{prefix}}.bcf.csi")
    log: f"cohorts/{cohort}/logs/tabix/index/{{prefix}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/index/{{prefix}}.tsv"
    params: "-p bcf"
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(tabix {params} {input}) > {log} 2>&1"


rule create_ped:
    input: config['cohort_yaml']
    output: f"cohorts/{cohort}/{cohort}.ped"
    log: f"cohorts/{cohort}/logs/yaml2ped/{cohort}.log"
    conda: "envs/pyyaml.yaml"
    message: f"Executing {{rule}}: Creating pedigree file for {cohort}."
    shell: f"(python3 workflow/scripts/yaml2ped.py {{input}} {cohort} {{output}}) > {{log}} 2>&1"


rule calculate_phrank:
    input:
        hpoterms = config['hpo']['terms'],
        hpodag = config['hpo']['dag'],
        hpoannotations = config['hpo']['annotations'],
        ensembltohgnc = config['ensembl_to_hgnc'],
        cohort_yaml = config['cohort_yaml']
    output: f"cohorts/{cohort}/{cohort}_phrank.tsv"
    log: f"cohorts/{cohort}/logs/calculate_phrank/{cohort}.log"
    conda: "envs/pyyaml.yaml"
    message: f"Executing {{rule}}: Calculate Phrank scores for {cohort}."
    shell:
        f"""(python3 workflow/scripts/calculate_phrank.py \
        {{input.hpoterms}} {{input.hpodag}} {{input.hpoannotations}} \
        {{input.ensembltohgnc}} {{input.cohort_yaml}} {cohort} {{output}}) > {{log}} 2>&1
        """
