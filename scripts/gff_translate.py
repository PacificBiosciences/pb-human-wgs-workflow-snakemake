#!/usr/bin/env python3
"""Translate chromosomes in Ensembl GFF3, ignore chromosomes missing from lookup."""

lookup_table = {f'{x}': f'chr{x}' for x in [*range(1,23), 'X', 'Y']}
lookup_table['MT'] = 'chrM'

with open('/dev/stdin', 'r') as gff:
    for row in gff:
        fields = row.split()
        if fields[0] == '##sequence-region':
            if fields[1] in lookup_table:
                fields[1] = lookup_table[fields[1]]
                print('\t'.join(fields))
        elif fields[0].startswith('#'):
            print('\t'.join(fields))
        else:
            if fields[0] in lookup_table:
                fields[0] = lookup_table[fields[0]]
                print('\t'.join(fields))
