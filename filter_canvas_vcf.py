#!/usr/bin/env python3
import sys
import os
import re
import logging
import argparse
from vase.ped_file import PedFile
from parse_vcf import VcfReader

ref_alleles = set((0, None))

def get_options():
    parser = argparse.ArgumentParser(description='''Filter CNVs from a Canvas
                                                    VCF file.''')
    parser.add_argument("vcf", help='''CNV VCF file produced by Canvas using
                        the small pedigree WGS option.''')
    parser.add_argument("-d", "--dq", type=float, default=0.0,
                        help='''Only output variants from the CNV.vcf.gz in the
                        results directory with a de novo quality (DQ) equal to
                        or greater than this value. Default=0.0''')
    parser.add_argument("-q", "--qs", type=float, default=0.0,
                        help='''Only output variants where the PHRED scale
                        quality score for samples with a variant is equal to or
                        greater than this value. Default=0.0''')
    parser.add_argument("-b", "--bc", type=int, default=0,
                        help='''Only output variants where the number of bins
                        for samples with a variant is equal to or greater than
                        this value. Default=0''')
    parser.add_argument("--pass_filters", action='store_true',
                        help='''When outputting variant regions only output
                        those which PASS filters.''')
    parser.add_argument("-s", "--samples", nargs='+', default=[],
                        help='''Only keep variants if called (and passing
                        filters) in these samples.''')
    parser.add_argument("-c", "--controls", nargs='+', default=[],
                        help='''Ignore variants if called in these samples.''')
    return parser

def is_variant(gt):
    return not set(gt).issubset(ref_alleles)

def main(vcf, dq=0.0, qs=0.0, bc=0, pass_filters=False, samples=[],
         controls=[]):
    gt_fields = ['GT', 'QS', 'BC']
    if dq:
        gt_fields.append('DQ')
    vreader = VcfReader(vcf)
    if not samples:
        samples = vreader.header.samples
    for c in samples + controls:
        if c not in vreader.header.samples:
            sys.exit("ERROR - sample {} not found in VCF".format(c))
    print(str(vreader.header).rstrip())
    for record in vreader:
        if pass_filters and record.FILTER != 'PASS':
            continue
        gts = record.parsed_gts(fields=gt_fields)
        if dq:
            dqs = dict((k,v) for k,v in gts['DQ'].items() if v is not None and
                       v>= dq)
            if not dqs:
                continue
        var_samps = [x for x in samples if is_variant(gts['GT'][x])]
        if not var_samps:
            continue
        if qs:
            pass_quals = [x for x in var_samps if gts['QS'][x] is not None and
                          gts['QS'][x] >= qs]
            if not pass_quals:
                continue
        if bc:
            pass_bins = [x for x in var_samps if gts['BC'][x] is not None and
                          gts['BC'][x] >= bc]
            if not pass_bins:
                continue
        var_gts = set(gts['GT'][x] for x in var_samps)
        control_gts = set(gts['GT'][x] for x in controls)
        print_record = False
        for gt in var_gts:
            if gt in control_gts:
                continue
            else:
                alleles = set(x for x in gt if x not in ref_alleles)
                for al in alleles:
                    if (al, al) in control_gts:
                        continue
                    print_record = True
        if print_record:
            print(record)


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))

