#!/usr/bin/env python3
"""
Identify tandem repeat regions in the genome that might have missed
variants due to coverage dropouts perhaps induced by repeat expansions.
"""


__version__ = "0.1.0"


import argparse
import collections
import pysam


class TandemRepeatRegion:
    def __init__(self, chrom, chromStart, chromEnd, motif, label, regionType):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.motif = motif
        self.label = label
        self.regionType = regionType


def main(args):
    # Read the tandem repeat regions of interest
    tandemRepeats = list()
    f = open(args.trbed)
    for l in f:
        cols = l.rstrip().split("\t")
        tandemRepeats.append(TandemRepeatRegion(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]))
    f.close()


    bam = pysam.AlignmentFile(args.bam, "rb")
    # Infer sample sex from the chrX:chr2 ratio (MALE if chrX:chr2 < 0.7).
    # The sex helps with interpreting coverage for chrX variants.
    chromReadsPerMb = dict()
    chromLens = dict()
    for l in pysam.idxstats(args.bam).rstrip("\n").split("\n"):
        chrom,chromLen,mappedReads,unmappedReads = l.split("\t")
        if float(chromLen)>1:
            chromReadsPerMb[chrom] = 1e6*float(mappedReads)/float(chromLen)
            chromLens[chrom] = int(chromLen)
    chrX = chromReadsPerMb.get("chrX", chromReadsPerMb.get("X"))
    chrY = chromReadsPerMb.get("chrY", chromReadsPerMb.get("Y"))
    chr2 = chromReadsPerMb.get("chr2", chromReadsPerMb.get("2"))
    if chrX is not None and chrY is not None and chr2 is not None:
        sex = "MALE" if chrX/chr2 < 0.7 else "FEMALE"
    else:
        sex = "UNKNOWN"


    # Evaluate the spanning read coverage for each region,
    # considering reads by phase (HP=1, HP=2, unphased).
    PAD = 1000
    for tandemRepeat in tandemRepeats:
        hp1Spanners,hp2Spanners,unphasedSpanners = [],[],[]
        readStarts = collections.defaultdict(int)     # map from read name to lowest start position in the region
        readEnds = collections.defaultdict(int)       # map from read name to highest end position in the region
        readHaplotypes = collections.defaultdict(set) # map from read name to (phase set, haplotype) tags

        # Fetch reads +/- PAD from the region of interest.
        # Record the lowest start and highest end postion of alignments for the read.
        alignments = bam.fetch(tandemRepeat.chrom, max(0,tandemRepeat.chromStart-PAD), min(tandemRepeat.chromEnd+PAD, chromLens[tandemRepeat.chrom]))
        for al in alignments:
            if al.query_name not in readStarts:
                readStarts[al.query_name] = al.reference_start
                readEnds[al.query_name] = al.reference_end
            readStarts[al.query_name] = min(readStarts[al.query_name], al.reference_start)
            readEnds[al.query_name] = min(readEnds[al.query_name], al.reference_end)
            if al.has_tag("HP") and al.has_tag("PS"):
                readHaplotypes[al.query_name].add((al.get_tag("PS"),al.get_tag("HP")))

        # Identify reads that span the region of interest.
        for query_name in readStarts:
            if readStarts[query_name] <= tandemRepeat.chromStart and readEnds[query_name] >= tandemRepeat.chromEnd:
                if len(readHaplotypes[query_name]) == 1:
                    if list(readHaplotypes[query_name])[0][1] == 1:
                        hp1Spanners.append(query_name)
                    else:
                        hp2Spanners.append(query_name)
                else:
                    unphasedSpanners.append(query_name)
        phasedSpanners = hp1Spanners + hp2Spanners
        totalSpanners = phasedSpanners + unphasedSpanners

        # Consider a region as a full dropout if it has fewer than `coverage` spanning reads.
        if len(totalSpanners) < args.coverage:
            print("%s\t%d\t%d\t%s\t%s\t%s\tFullDropout" % (tandemRepeat.chrom, tandemRepeat.chromStart, tandemRepeat.chromEnd,
                tandemRepeat.motif, tandemRepeat.label, tandemRepeat.regionType))
        # If the region is phased (i.e. spanned by fewer than `coverage` unphased reads) then evaluate haplotype-specific coverage.
        elif len(unphasedSpanners) < args.coverage:
            # Ignore chrX variants in a male.
            if len(phasedSpanners) and not (sex == "MALE" and tandemRepeat.chrom in ("X","chrX")):
                if len(hp1Spanners) < args.coverage or len(hp2Spanners) < args.coverage:
                    print("%s\t%d\t%d\t%s\t%s\t%s\tHaplotypeDropout" % (tandemRepeat.chrom, tandemRepeat.chromStart, tandemRepeat.chromEnd,
                        tandemRepeat.motif, tandemRepeat.label, tandemRepeat.regionType))

    bam.close()


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trbed", help="Regions of interest, BED4+ (chrom<TAB>chromStart<TAB>chromEnd<TAB>motif<TAB>label<TAB>region)")
    parser.add_argument("bam", help="Alignments, BAM")
    parser.add_argument("--coverage", type=int, default=2, help="Minimum coverage to consider a region as covered [default: %(default)s]")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
