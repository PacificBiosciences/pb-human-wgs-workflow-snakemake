ruleorder: samtools_fasta > seqtk_fastq_to_fasta


rule samtools_fasta:
    input: lambda wildcards: ubam_dict[wildcards.movie]
    output: temp(f"samples/{sample}/fasta/{{movie}}.fasta")
    log: f"samples/{sample}/logs/samtools/fasta/{{movie}}.log"
    benchmark: f"samples/{sample}/benchmarks/samtools/fasta/{{movie}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Converting {input} to {output}."
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule seqtk_fastq_to_fasta:
    input: lambda wildcards: fastq_dict[wildcards.movie]
    output: temp(f"samples/{sample}/fasta/{{movie}}.fasta")
    log: f"samples/{sample}/logs/seqtk/seq/{{movie}}.log"
    benchmark: f"samples/{sample}/benchmarks/seqtk/seq/{{movie}}.tsv"
    conda: "envs/seqtk.yaml"
    message: "Converting {input} to {output}."
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule hifiasm_assemble:
    input: expand(f"samples/{sample}/fasta/{{movie}}.fasta", movie=movies)
    output:
        hap1_p_ctg = temp(f"samples/{sample}/hifiasm/{sample}.asm.hap1.p_ctg.gfa"),
        hap1_p_ctg_bed = f"samples/{sample}/hifiasm/{sample}.asm.hap1.p_ctg.lowQ.bed",
        hap1_p_ctg_noseq = f"samples/{sample}/hifiasm/{sample}.asm.hap1.p_ctg.noseq.gfa",
        hap2_p_ctg = temp(f"samples/{sample}/hifiasm/{sample}.asm.hap2.p_ctg.gfa"),
        hap2_p_ctg_bed = f"samples/{sample}/hifiasm/{sample}.asm.hap2.p_ctg.lowQ.bed",
        hap2_p_ctg_noseq = "samples/{sample}/hifiasm/{sample}.asm.hap2.p_ctg.noseq.gfa",
        p_ctg = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_ctg.gfa"),
        p_ctg_bed = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_ctg.lowQ.bed"),
        p_ctg_noseq = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_ctg.noseq.gfa"),
        p_utg = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_utg.gfa"),
        p_utg_bed = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_utg.lowQ.bed"),
        p_utg_noseq = temp(f"samples/{sample}/hifiasm/{sample}.asm.p_utg.noseq.gfa"),
        r_utg = temp(f"samples/{sample}/hifiasm/{sample}.asm.r_utg.gfa"),
        r_utg_bed = temp(f"samples/{sample}/hifiasm/{sample}.asm.r_utg.lowQ.bed"),
        r_utg_noseq = temp(f"samples/{sample}/hifiasm/{sample}.asm.r_utg.noseq.gfa"),
        ec_bin = temp(f"samples/{sample}/hifiasm/{sample}.asm.ec.bin"),
        ovlp_rev_bin = temp(f"samples/{sample}/hifiasm/{sample}.asm.ovlp.reverse.bin"),
        ovlp_src_bin = temp(f"samples/{sample}/hifiasm/{sample}.asm.ovlp.source.bin")
    log: f"samples/{sample}/logs/hifiasm.log"
    benchmark: f"samples/{sample}/benchmarks/hifiasm.tsv"
    conda: "envs/hifiasm.yaml"
    params: prefix = f"samples/{sample}/hifiasm/{sample}.asm"
    threads: 48
    message: f"Assembling sample {sample} from {{input}}"
    shell: "(hifiasm -o {params.prefix} -t {threads} {input}) > {log} 2>&1"


rule gfa2fa:
    input: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.gfa"
    output: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.fasta"
    log: f"samples/{sample}/logs/gfa2fa/{sample}.asm.{{infix}}.log"
    benchmark: f"samples/{sample}/benchmarks/gfa2fa/{sample}.asm.{{infix}}.tsv"
    conda: "envs/gfatools.yaml"
    message: "Extracting fasta from assembly {input}."
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.fasta"
    output: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.fasta.gz"
    log: f"samples/{sample}/logs/bgzip/{sample}.asm.{{infix}}.log"
    benchmark: f"samples/{sample}/benchmarks/bgzip/{sample}.asm.{{infix}}.tsv"
    threads: 4
    conda: "envs/htslib.yaml"
    message: "Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule asm_stats:
    input: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.fasta.gz"
    output: f"samples/{sample}/hifiasm/{sample}.asm.{{infix}}.fasta.stats.txt"
    log: f"samples/{sample}/logs/asm_stats/{sample}.asm.{{infix}}.fasta.log"
    benchmark: f"samples/{sample}/benchmarks/asm_stats/{sample}.asm.{{infix}}.fasta.tsv"
    conda: "envs/k8.yaml"
    message: "Calculating stats for {input}."
    shell: f"(k8 workflow/scripts/calN50/calN50.js -f {config['ref']['index']} {{input}} > {{output}}) > {{log}} 2>&1"


rule align_hifiasm:
    input:
        target = config['ref']['fasta'],
        query = [f"samples/{sample}/hifiasm/{sample}.asm.{infix}.fasta.gz"
                 for infix in ["hap1.p_ctg", "hap2.p_ctg"]]
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.bam"
    log: f"samples/{sample}/logs/align_hifiasm/{sample}.asm.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/align_hifiasm/{sample}.asm.{ref}.tsv"
    params:
        max_chunk = 200000,
        minimap2_args = "-L --secondary=no --eqx -ax asm5",
        minimap2_threads = 10,
        readgroup = f"@RG\\tID:{sample}_hifiasm\\tSM:{sample}",
        samtools_threads = 3
    threads: 16  # minimap2 + samtools(+1) + 2x awk + seqtk + cat
    conda: "envs/align_hifiasm.yaml"
    message: "Aligning {input.query} to {input.target}."
    shell:
        """
        (cat {input.query} \
            | seqtk seq -l {params.max_chunk} - \
            | awk '{{ if ($1 ~ />/) {{ n=$1; i=0; }} else {{ i++; print n "." i; print $0; }} }}' \
            | minimap2 -t {params.minimap2_threads} {params.minimap2_args} \
                -R '{params.readgroup}' {input.target} - \
                | awk '{{ if ($1 !~ /^@/) \
                                {{ Rct=split($1,R,"."); N=R[1]; for(i=2;i<Rct;i++) {{ N=N"."R[i]; }} print $0 "\tTG:Z:" N; }} \
                              else {{ print; }} }}' \
                | samtools sort -@ {params.samtools_threads} > {output}) > {log} 2>&1
        """


rule htsbox:
    input:
        query = [f"samples/{sample}/hifiasm/{sample}.asm.{infix}.fasta.gz"
                 for infix in ["hap1.p_ctg", "hap2.p_ctg"]],
        reference = config['ref']['fasta']
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf"
    log: f"samples/{sample}/logs/htsbox/{sample}.asm.log"
    benchmark: f"samples/{sample}/benchmarks/htsbox/{sample}.asm.tsv"
    params:
        minimap2_threads = 10,
        minimap2_args = "-a -k19 -w10 -O5,56 -E4,1 -A2 -B5 -z400,50 -r2000 \
                         -g5000 -L -Y --lj-min-ratio 0.5 --eqx --secondary=no",
        readgroup = f"@RG\\tID:{sample}_hifiasm\\tSM:{sample}",
        samtools_threads = 3,
        htsbox_args = "-q20"
    threads: 16
    conda: "envs/htsbox.yaml"
    message: "Calling variants from {input.query} using htsbox."
    shell:
        """
        (zcat {input.query} \
            | minimap2 -t {params.minimap2_threads} {params.minimap2_args} \
                -R '{params.readgroup}' {input.reference} - \
            | samtools sort -@ {params.samtools_threads} \
            | htsbox pileup {params.htsbox_args} -c -f {input.reference} - \
            > {output}) > {log} 2>&1
        """


rule htsbox_bcftools_stats:
    input: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf.gz"
    output: f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.vcf.stats.txt"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.asm.{ref}.htsbox.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/stats/{sample}.asm.{ref}.htsbox.vcf.tsv"
    params: f"--fasta-ref {config['ref']['fasta']} -s samples/{sample}/hifiasm/{sample}.asm.{ref}.bam"
    threads: 4
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Calculating VCF statistics for {input}."
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"

