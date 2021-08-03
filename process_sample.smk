import re
from pathlib import Path


shell.prefix("set -o pipefail; umask 002; ")  # set g+w
configfile: "workflow/reference.yaml"         # reference information
configfile: "workflow/config.yaml"            # general configuration


# sample will be provided at command line with `--config sample=$SAMPLE`
sample = config['sample']
ref = config['ref']['shortname']
all_chroms = config['ref']['autosomes'] + config['ref']['sex_chrom'] + config['ref']['mit_chrom']
print(f"Processing sample {sample} with reference {ref}.")

# scan samples/{sample}/aligned to generate a dict of aBAMs
pattern = re.compile(r'samples/(?P<sample>[A-Za-z0-9_-]+)/aligned/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6})\.(?P<reference>.*).bam')
abam_dict = {}
for infile in Path(f"samples/{sample}/aligned").glob('*.bam'):
    match = pattern.search(str(infile))
    if match:
        if ref == match.group('reference'):
            abam_dict[match.group('movie')] = str(infile)
movies = list(abam_dict.keys())  # list of all movie contexts for this sample
abams = list(abam_dict.values()) # list of all aBAM paths for this sample

print(f"{sample} movies available: {movies}")

# scan smrtcells/ready directory for uBAMs or FASTQs for hifiasm
# uBAMs will have priority over FASTQs in downstream processes if both are available
ubam_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam')
ubam_dict = {}
fastq_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).fastq.gz')
fastq_dict = {}
for infile in Path(f'smrtcells/ready/{sample}').glob('*.bam'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match and (ubam_match.group('movie') in movies):
        # create a dict to link movie context to uBAM filenames
        ubam_dict[ubam_match.group('movie')] = str(infile)
for infile in Path(f'smrtcells/ready/{sample}').glob('*.fastq.gz'):
    fastq_match = fastq_pattern.search(str(infile))
    if fastq_match and (fastq_match.group('movie') in movies):
        # create a dict to link movie context to FASTQ filenames
        fastq_dict[fastq_match.group('movie')] = str(infile)


# build a list of targets
targets = []
include: 'rules/sample_common.smk'

# call structural variants with pbsv
include: 'rules/sample_pbsv.smk'
if 'pbsv_vcf' in config['sample_targets']:
    # pbsv VCFs
    targets.extend([f"samples/{sample}/pbsv/{sample}.{ref}.pbsv.{suffix}"
                    for suffix in ['vcf.gz', 'vcf.gz.tbi']])

# call small variants with DeepVariant
include: 'rules/sample_deepvariant.smk'
if 'deepvariant' in config['sample_targets']:
    # deepvariant VCFs, gVCFs, reports, and stats
    targets.extend([f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.{suffix}"
                    for suffix in ['vcf.gz', 'vcf.gz.tbi', 'g.vcf.gz', 'g.vcf.gz.tbi',
                                'visual_report.html', 'vcf.stats.txt']])

# phase small variants with WhatsHap and haplotag BAM
include: 'rules/sample_whatshap.smk'
if 'whatshap' in config['sample_targets']:
    # phased VCFs, stats, phase block GTFs, and haplotagged BAMs
    targets.extend([f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.{suffix}"
                    for suffix in ['phased.vcf.gz', 'phased.vcf.gz.tbi', 'phased.gtf',
                                'phased.tsv', 'phased.blocklist',
                                'haplotagged.bam', 'haplotagged.bam.bai']])

# genotype STRs
include: 'rules/sample_tandem_genotypes.smk'
if 'tandem-genotypes' in config['sample_targets']:
    # tandem-genotypes tabular output and plots
    targets.extend([f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.{suffix}"
                   for suffix in ['txt', 'pdf', 'dropouts.txt']])

# calculate coverage of haplotagged sample aBAM with mosdepth
include: 'rules/sample_mosdepth.smk'
include: 'rules/sample_gc_coverage.smk'
if 'coverage' in config['sample_targets']:
    # coverage from merged haplotagged aBAM
    targets.extend([f"samples/{sample}/mosdepth/{sample}.{ref}.deepvariant.haplotagged.{suffix}"
                    for suffix in ['mosdepth.global.dist.txt', 'mosdepth.region.dist.txt',
                                'mosdepth.summary.txt', 'regions.bed.gz']])
    targets.extend([f"samples/{sample}/mosdepth/{sample}.{ref}.gc_coverage.summary.txt"])

# merge kmers with jellyfish
include: 'rules/sample_jellyfish.smk'
include: 'rules/sample_kmer_consistency.smk'
if 'kmers' in config['sample_targets']:
    # jellyfish kmer database
    targets.extend([f"samples/{sample}/jellyfish/{sample}.{suffix}"
                    for suffix in ['jf']])
    # measure kmer consistency of SMRT Cells
    targets.append(f"samples/{sample}/jellyfish/{sample}.kmerconsistency.txt")

# assemble with hifiasm
include: 'rules/sample_hifiasm.smk'
if 'assembly' in config['sample_targets']:
    # assembly and stats
    targets.extend([f"samples/{sample}/hifiasm/{sample}.asm.{infix}.{suffix}"
                for suffix in ['fasta.gz', 'fasta.stats.txt']
                for infix in ['a_ctg', 'p_ctg']])
    # assembly alignments
    targets.extend([f"samples/{sample}/hifiasm/{sample}.asm.{ref}.{suffix}"
                for suffix in ['bam', 'bam.bai']])
    # assembly htsbox variants
    targets.extend([f"samples/{sample}/hifiasm/{sample}.asm.{ref}.htsbox.{suffix}"
                for suffix in ['vcf.gz', 'vcf.gz.tbi', 'vcf.stats.txt']])


localrules: all


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
