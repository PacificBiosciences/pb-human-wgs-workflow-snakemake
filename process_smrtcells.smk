import re
from pathlib import Path
from collections import defaultdict


shell.prefix("set -o pipefail; umask 002; ")  # set g+w
configfile: "workflow/reference.yaml"         # reference configuration
configfile: "workflow/config.yaml"            # general configuration

ref = config['ref']['shortname']

# scan smrtcells/ready directory for uBAMs or FASTQs that are ready to process
# uBAMs have priority over FASTQs in downstream processes if both are available
ubam_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam')
ubam_dict = defaultdict(dict)
fastq_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).fastq.gz')
fastq_dict = defaultdict(dict)
for infile in Path('smrtcells/ready').glob('**/*.bam'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match:
        # create a dict-of-dict to link samples to movie context to uBAM filenames
        ubam_dict[ubam_match.group('sample')][ubam_match.group('movie')] = str(infile)
for infile in Path('smrtcells/ready').glob('**/*.fastq.gz'):
    fastq_match = fastq_pattern.search(str(infile))
    if fastq_match:
        # create a dict-of-dict to link samples to movie context to FASTQ filenames
        fastq_dict[fastq_match.group('sample')][fastq_match.group('movie')] = str(infile)

print(f"uBAMs available for samples: {list(ubam_dict.keys())}")
for s in ubam_dict.keys():
    print(f"movies available for {s}: {list(ubam_dict[s].keys())}")

print(f"FASTQs available for samples: {list(fastq_dict.keys())}")
for s in fastq_dict.keys():
    print(f"movies available for {s}: {list(fastq_dict[s].keys())}")


# build a list of targets
targets = []

# add ubams and fastqs to targets so that md5sums will be written
[targets.extend(ubam_dict[s].values()) for s in ubam_dict.keys()]
[targets.extend(fastq_dict[s].values()) for s in fastq_dict.keys()]

# align reads with pbmm2
include: 'rules/smrtcell_pbmm2.smk'
if 'alignment' in config['smrtcells_targets']:
    targets.extend([f"samples/{sample}/aligned/{movie}.{ref}.{suffix}"
                    for suffix in ['bam', 'bam.bai']
                    for sample in list(ubam_dict.keys())
                    for movie in list(ubam_dict[sample].keys())])  # aBAMs from uBAMs
    targets.extend([f"samples/{sample}/aligned/{movie}.{ref}.{suffix}"
                    for suffix in ['bam', 'bam.bai']
                    for sample in list(fastq_dict.keys())
                    for movie in list(fastq_dict[sample].keys())]) # aBAMs from FASTQs

# calculate coverage with mosdepth
include: 'rules/smrtcell_mosdepth.smk'
if 'coverage' in config['smrtcells_targets']:
    targets.extend([f"samples/{sample}/mosdepth/{movie}.{ref}.{suffix}"
                    for suffix in ['mosdepth.global.dist.txt', 'mosdepth.region.dist.txt',
                                'mosdepth.summary.txt', 'regions.bed.gz']
                    for sample in list(ubam_dict.keys())
                    for movie in list(ubam_dict[sample].keys())])  # coverage from uBAMs
    targets.extend([f"samples/{sample}/mosdepth/{movie}.{ref}.{suffix}"
                    for suffix in ['mosdepth.global.dist.txt', 'mosdepth.region.dist.txt',
                                'mosdepth.summary.txt', 'regions.bed.gz']
                    for sample in list(fastq_dict.keys())
                    for movie in list(fastq_dict[sample].keys())]) # coverage from FASTQs

# smrtcell summary statistics
include: 'rules/smrtcell_stats.smk'
if 'stats' in config['smrtcells_targets']:
    targets.extend([f"samples/{sample}/smrtcell_stats/{movie}.{suffix}"
                    for suffix in ['read_length_and_quality.tsv', 'read_length_summary.tsv', 'read_quality_summary.tsv']
                    for sample in list(ubam_dict.keys())
                    for movie in list(ubam_dict[sample].keys())])  # summary stats from uBAMs
    targets.extend([f"samples/{sample}/smrtcell_stats/{movie}.{suffix}"
                    for suffix in ['read_length_and_quality.tsv', 'read_length_summary.tsv', 'read_quality_summary.tsv']
                    for sample in list(fastq_dict.keys())
                    for movie in list(fastq_dict[sample].keys())]) # summary stats from FASTQs

# coverage-based library and sample QC metrics: mtDNA:autosome ratio and inferred chromosomal sex
include: 'rules/smrtcell_coverage_qc.smk'
if 'coverage_qc' in config['smrtcells_targets']:
    targets.extend([f"samples/{sample}/mosdepth/{movie}.{ref}.{suffix}"
                    for suffix in ['mosdepth.inferred_sex.txt', 'mosdepth.M2_ratio.txt', 'gc_coverage.summary.txt']
                    for sample in list(ubam_dict.keys())
                    for movie in list(ubam_dict[sample].keys())])  # QC from uBAMs
    targets.extend([f"samples/{sample}/mosdepth/{movie}.{ref}.{suffix}"
                    for suffix in ['mosdepth.inferred_sex.txt', 'mosdepth.M2_ratio.txt', 'gc_coverage.summary.txt']
                    for sample in list(fastq_dict.keys())
                    for movie in list(fastq_dict[sample].keys())]) # QC from FASTQs

# count kmers with jellyfish and save modimers
include: 'rules/smrtcell_jellyfish.smk'
if 'kmers' in config['smrtcells_targets']:
    targets.extend([f"samples/{sample}/jellyfish/{movie}.{suffix}"
                    for suffix in ['jf', 'modimers.tsv.gz']
                    for sample in list(ubam_dict.keys())
                    for movie in list(ubam_dict[sample].keys())])  # kmers from uBAMs
    targets.extend([f"samples/{sample}/jellyfish/{movie}.{suffix}"
                    for suffix in ['jf', 'modimers.tsv.gz']
                    for sample in list(fastq_dict.keys())
                    for movie in list(fastq_dict[sample].keys())]) # kmers from FASTQs


localrules: all


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
