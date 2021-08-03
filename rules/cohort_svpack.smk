rule svpack_filter_annotated:
    input:
        pbsv_vcf = svpack_input,
        eee_vcf = config['ref']['eee_vcf'],
        gnomadsv_vcf = config['ref']['gnomadsv_vcf'],
        hprc_pbsv_vcf = config['ref']['hprc_pbsv_vcf'],
        decode_vcf = config['ref']['decode_vcf'],
        gff = config['ref']['ensembl_gff']
    output: f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf"
    log: f"cohorts/{cohort}/logs/svpack/{cohort}.{ref}.pbsv.svpack.log"
    benchmark: f"cohorts/{cohort}/benchmarks/svpack/{cohort}.{ref}.pbsv.svpack.tsv"
    params:
        min_sv_length = 50
    conda: "envs/svpack.yaml"
    shell:
        """
        (python workflow/scripts/svpack/svpack filter --pass-only \
            --min-svlen {params.min_sv_length} {input.pbsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.eee_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.gnomadsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.hprc_pbsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.decode_vcf} | \
            python workflow/scripts/svpack/svpack consequence - {input.gff} | \
            python workflow/scripts/svpack/svpack tagzygosity - > {output}) 2> {log}
        """


info_fields = [
    'SVTYPE',
    'SVLEN',
    'SVANN',
    'CIPOS',
    'MATEID',
    'END'
    ]


rule slivar_svpack_tsv:
    input:
        filt_vcf = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf.gz",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        lof_lookup = config['lof_lookup'],
        clinvar_lookup = config['clinvar_lookup'],
        phrank_lookup = f"cohorts/{cohort}/{cohort}_phrank.tsv"
    output:
        filt_tsv = f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.tsv",
    log: f"cohorts/{cohort}/logs/slivar/tsv/{cohort}.{ref}.pbsv.svpack.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/tsv/{cohort}.{ref}.pbsv.svpack.tsv"
    params: info = "".join([f"--info-field {x} " for x in info_fields])
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Converting annotated VCFs to TSVs for easier interpretation."
    shell:
        """
        (slivar tsv \
            {params.info} \
            --sample-field hetalt \
            --sample-field homalt \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out {output.filt_tsv} \
            {input.filt_vcf}) > {log} 2>&1
        """
