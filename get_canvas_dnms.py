#!/usr/bin/env python3
import sys
import os
import re
import logging
import argparse
from vase.ped_file import PedFile
from parse_vcf import VcfReader

def get_options():
    parser = argparse.ArgumentParser(description='''Extract de novo CNVs from a
                                     Canvas VCF file.''')
    parser.add_argument("vcf", help='''CNV VCF file produced by Canvas using
                        the small pedigree WGS option.''')
    parser.add_argument("-d", "--dq", type=float, default=0.0,
                        help='''Only output variants from the CNV.vcf.gz in the
                        results directory with a de novo quality (DQ) equal to
                        or greater than this value. Default=0.0''')
    parser.add_argument("-q", "--qs", type=float, default=0.0,
                        help='''Only output variants where the PHRED scale
                        quality score for the child is equal to or greater than
                        this value. Default=0.0''')
    parser.add_argument("-b", "--bc", type=int, default=0,
                        help='''Only output variants where the number of bins
                        for the child is equal to or greater than this value.
                        Default=0''')
    parser.add_argument("--pass_filters", action='store_true',
                        help='''When outputting variant regions only output
                        those which PASS filters.''')
    return parser

def main(vcf, dq=0.0, qs=0.0, bc=0, pass_filters=False,):
    vreader = VcfReader(vcf)
    print(str(vreader.header), end='')
    for record in vreader:
        if pass_filters and record.FILTER != 'PASS':
            continue
        gts = record.parsed_gts(fields=['DQ', 'QS', 'BC'])
        dqs = dict((k,v) for k,v in gts['DQ'].items() if v is not None and
                   v>= dq)
        if not dqs:
            continue
        if qs:
            pass_quals = [x for x in dqs if gts['QS'][x] is not None and
                          gts['QS'][x] >= qs]
            if not pass_quals:
                continue
        if bc:
            pass_bins = [x for x in dqs if gts['BC'][x] is not None and
                          gts['BC'][x] >= bc]
            if not pass_bins:
                continue
        print(record)


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))

