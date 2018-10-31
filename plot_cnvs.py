#!/usr/bin/env python3
import sys
import os
import re
import gzip
import logging
import operator
import argparse
import matplotlib
if os.environ.get('DISPLAY','') == '':
    matplotlib.use('Agg') #for headless servers
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.collections import BrokenBarHCollection
import seaborn as sns
import pandas as pd
import het_plotter.plot_heterozygosity as het_counts
from vase.ped_file import PedFile
from parse_vcf import VcfReader

logger = logging.getLogger("CNV plotter")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
       '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

reg_matcher = re.compile(r'''^((chr)?\S):(\d+)-(\d+)$''')

def get_cytobands(build):
    refdir = os.path.join(os.path.dirname(__file__), 'data', build)
    cyto_file = os.path.join(refdir, 'cytoBandIdeo.txt.gz')
    if os.path.exists(cyto_file):
        logger.info("Reading cytoband information for {} genome".format(build))
        ideo = pd.read_table(
            cyto_file,
            names=['chrom', 'start', 'end', 'name', 'gieStain'])
        ideo['width'] = ideo.end - ideo.start
        color_lookup = {
            'gneg': (1., 1., 1.),
            'gpos25': (.6, .6, .6),
            'gpos50': (.4, .4, .4),
            'gpos75': (.2, .2, .2),
            'gpos100': (0., 0., 0.),
            'acen': (.8, .4, .4),
            'gvar': (.8, .8, .8),
            'stalk': (.9, .9, .9),
        }
        ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])
        return ideo
    else:
        logger.warn("No cytobands file for build {}".format(build))
        logger.info("To plot an ideogram for this genome download the " +
                    "relevant cytoBandIdeo.txt.gz from genome.ucsc.edu and " +
                    "place in the directory {}".format(refdir))
    return None

def read_coverage(results_dir, samples=None, chrom=None):
    df = pd.DataFrame()
    logger.info("Reading coverage files in {}".format(results_dir))
    for vis_tmp in [d for d in os.listdir(results_dir) if
                    d.startswith('VisualizationTemp')]:
        sample = vis_tmp.replace('VisualizationTemp', '')
        if samples is not None and sample not in samples:
            logger.warn("Skipping directory {} as sample ".format(vis_tmp) +
                        "{} is not in provided PED file.".format(sample))
            continue
        f = os.path.join(results_dir, vis_tmp, "coverage.bedgraph")
        if not os.path.exists(f):
            if os.path.exists(f + '.gz'):
                f = f + '.gz'
            else:
                sys.exit("ERROR: No coverage.bedgraph file in " +
                         "{}".format(vis_tmp))
        logger.info("Reading {}".format(f))
        temp_df = pd.read_table(f, names=('Chrom', 'Start', 'End', 'Ploidy'))
        if chrom:
            temp_df = temp_df[temp_df.Chrom == chrom]
        temp_df['Sample'] = sample
        df = pd.concat([df, temp_df])
    return df

def plot_chromosomes(df, outdir, ideo, samples, fig_dimensions, ymax=6,
                     het_vcf=None, window_length=1e5):
    samples = check_samples_in_df(df, samples)
    cols = sns.color_palette("colorblind", len(samples))
    for chrom in df.Chrom.unique():
        pan_per_sample = 1
        het_df = None
        if het_vcf is not None:
            logger.info("Getting heterzoygosity counts for chromosome " +
                        "{}".format(chrom))
            gt_counts = het_counts.get_gt_counts(het_vcf, samples, chrom,
                                                 logger=logger, plot=False)
            het_df = het_counts.counts_to_df(gt_counts, window_length)
            pan_per_sample = 2
        logger.info("Plotting chromosome {}".format(chrom))
        png = os.path.join(outdir, chrom + '.png')
        chr_df = df[df.Chrom == chrom]
        i = 0
        fig, axes = plt.subplots(nrows=1 + (pan_per_sample * len(samples)),
                                 ncols=1,
                                 sharex=True,
                                 gridspec_kw={'height_ratios':[1] + [8] *
                                              (pan_per_sample * len(samples)),
                                              'hspace': 0.4},
                                 figsize=fig_dimensions)
        if ideo is not None:
            chr_ideo = ideo[ideo.chrom == chrom]
            axes[0].add_collection(BrokenBarHCollection(
                chr_ideo[['start', 'width']].values,
                (0, 1),
                facecolors=chr_ideo['colors']))
            axes[0].axes.set_xlim((0, max(chr_df.End)))
            axes[0].axes.set_ylim((0.1, 1))
        for s in samples:
            row = i * pan_per_sample + 1
            tmp_df = chr_df[chr_df.Sample == s]
            alpha = min(50000.0/len(tmp_df), 0.5)
            axes[row].scatter(tmp_df.Start, tmp_df.Ploidy, color=cols[i],
                        alpha=alpha, s=[0.1] * len(tmp_df))
            axes[row].plot([chr_df.Start.min(), chr_df.Start.max()],
                             [2.0, 2.0 ],
                             '--',
                             color='black',
                             alpha=0.4)
            axes[row].axes.set_ylim((-0.5, ymax))
            axes[row].axes.set_xlim((0, max(chr_df.End)))
            axes[row].axes.set_title(s)
            axes[row].axes.set_ylabel('CN')
            if het_df is not None:
                samp_df = het_df[het_df.sample_id == s]
                axes[row + 1].plot(samp_df['pos'], samp_df['het'], '-.',
                                   color=cols[i])
                axes[row + 1].axis('tight')
                axes[row + 1].axes.set_xlim((0, max(chr_df.End)))
                axes[row + 1].axes.set_ylabel('Fraction Heterozygosity')
                axes[row].axes.set_ylim((-0.5, 1.0))
            if row == len(samples) * pan_per_sample:#add label for bottom plot
                axes[row].axes.set_xlabel('Pos')
                #make the position labels a bit prettier
                axes[row].xaxis.set_major_formatter(
                    mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
            i += 1
        axes[0].set_title(chrom)
        axes[0].axis('tight')
        axes[0].tick_params(axis='both',
                            which='both',
                            bottom=False,
                            top=False,
                            left=False,
                            right=False,
                            labelbottom=False,
                            labelleft=False,
                            labelcolor='white')
        fig.savefig(png)
        plt.close('all')

def get_region_from_record(record):
    l = record.SPAN - record.POS
    padding = int(max(5000, l))
    start = record.POS - padding
    end = record.SPAN + padding
    return (start, end)

def check_samples_in_df(df, samples):
    if samples:
        s_uniq = df.Sample.unique()
        samples = [x for x in samples if x in s_uniq]
    else:
        samples = list(df.Sample.unique())
    if not samples:
        sys.exit("No data to process! Please check that your results " +
                 "directory is not empty and if using a PED file that your " +
                 "PED file and results directory contain at least one " +
                 "overlapping sample.")
    return samples

def plot_variants(df, vcf, outdir, samples, ideo, fig_dimensions, dq=None,
                  ymax=6, filters=False, chrom=None, het_vcf=None,
                  window_length=1e5):
    samples = check_samples_in_df(df, samples)
    vreader = VcfReader(vcf)
    if chrom is not None:
        vreader.set_region(chrom=chrom)
    region_name = None
    prev_chrom = None
    het_df = None
    for record in vreader:
        if record.ALT == '.':
            continue
        if filters and record.FILTER != 'PASS':
            continue
        (start, end) = get_region_from_record(record)
        alt = record.ALT.replace("<", "").replace(">", "")
        if dq is not None:
            gts = record.parsed_gts(fields=['DQ'])
            dqs = [x for x in gts['DQ'].values() if x is not None and x >= dq]
            if not dqs:
                continue
            region_name = "DQ_{}_{}_{}_{}_{}".format("/".join((str(x) for x in
                                                               dqs)),
                                                     record.CHROM,
                                                     record.POS,
                                                     record.SPAN,
                                                     alt)
        else:
            region_name = "{}_{}_{}_{}".format(record.CHROM, record.POS,
                                            record.SPAN, alt)
        if prev_chrom is None or record.CHROM != prev_chrom:
            chr_df = df[df.Chrom == record.CHROM]
            chr_ideo = ideo[(ideo.chrom == record.CHROM)]
            if het_vcf is not None:
                logger.info("Getting heterzoygosity counts for chromosome " +
                            "{}".format(record.CHROM))
                gt_counts = het_counts.get_gt_counts(het_vcf, samples,
                                                     record.CHROM,
                                                     logger=logger, plot=False)
                het_df = het_counts.counts_to_df(gt_counts, window_length)
            prev_chrom = record.CHROM
        var_span = (record.POS, record.SPAN)
        plot_region(df=chr_df, chrom=record.CHROM, start=start, end=end,
                    outdir=outdir, ideo=chr_ideo, samples=samples,
                    name=region_name, ymax=ymax, fig_dimensions=fig_dimensions,
                    roi=var_span, het_df=het_df)

def plot_region(df, chrom, start, end, outdir, ideo, samples, fig_dimensions,
                ymax=6, name=None, roi=None, het_df=None):
    cols = sns.color_palette("colorblind", len(samples))
    region = "{}:{}-{}".format(chrom, start, end)
    name = name if name else region
    logger.info("Plotting region " + name)
    png = os.path.join(outdir, name + '.png')
    pan_per_sample = 1 if het_df is None else 2
    fig, axes = plt.subplots(nrows=1 + (pan_per_sample * len(samples)),
                             ncols=1,
                             sharex=True,
                             gridspec_kw={'height_ratios':[1] + [8] *
                                          (pan_per_sample * len(samples)),
                                          'hspace': 0.4},
                             figsize=fig_dimensions)
    reg_df = df[(df.Start >= start) & (df.Start <= end)]
    if ideo is not None:
        reg_ideo = ideo[(ideo.start >= start) & (ideo.start < end)]
        axes[0].add_collection(BrokenBarHCollection(
            reg_ideo[['start', 'width']].values,
            (0, 1),
            facecolors=reg_ideo['colors']))
        axes[0].axes.set_xlim((start, end))
        axes[0].axes.set_ylim((0.1, 1))
    i = 0
    for s in samples:
        row = i * pan_per_sample + 1
        tmp_df = reg_df[reg_df.Sample == s]
        alpha = min(50000.0/len(tmp_df), 0.5)
        size = min(10000.0/len(tmp_df), 100)
        axes[row].scatter(tmp_df.Start, tmp_df.Ploidy, color=cols[i],
                            alpha=alpha, s=[size] * len(tmp_df))
        axes[row].plot([start, end],
                         [2.0, 2.0 ],
                         '--',
                         color='black',
                         alpha=0.4)
        if roi:
            axes[row].plot([roi[0], roi[1]],
                             [ymax-0.5, ymax-0.5],
                             '-',
                             color='red',
                             alpha=0.8)
        axes[row].axes.set_xlim((start, end))
        axes[row].axes.set_ylim((-0.5, ymax))
        axes[row].axes.set_title(s)
        axes[row].axes.set_ylabel('CN', rotation='horizontal', ha='right')
        if het_df is not None:
            samp_df = het_df[(het_df.sample_id == s) & (het_df.pos >= start) &
                             (het_df.pos <= end)]
            axes[row + 1].plot(samp_df['pos'], samp_df['het'], '-.',
                               color=cols[i])
            axes[row + 1].axes.set_xlim((start, end))
            axes[row + 1].axes.set_ylabel('Het', rotation='horizontal',
                                          ha='right')
            axes[row].axes.set_ylim((-0.5, 1.0))
        i += 1
        if row == len(samples) * pan_per_sample:#add label for bottom plot
            axes[row].axes.set_xlabel('Pos')
            #make the position labels a bit prettier
            axes[row].xaxis.set_major_formatter(
                mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
    axes[0].set_title(chrom)
    axes[0].tick_params(axis='both',
                        which='both',
                        bottom=False,
                        top=False,
                        left=False,
                        right=False,
                        labelbottom=False,
                        labelleft=False,
                        labelcolor='white')
    fig.savefig(png)
    plt.close('all')



def sample_order_from_ped(pedfile):
    samples = []
    ped = PedFile(pedfile)
    for family in ped.families:
        samples.extend(indv.iid for indv in
                       ped.families[family].individuals.values() if
                       indv.is_affected())
        samples.extend(indv.iid for indv in
                       ped.families[family].individuals.values() if
                       indv.is_unknown_phenotype())
        samples.extend(indv.iid for indv in
                       ped.families[family].individuals.values() if
                       indv.is_unaffected())
    return samples

def get_options():
    parser = argparse.ArgumentParser(description='Plot Canvas CNV results.')
    reg_group = parser.add_mutually_exclusive_group()
    parser.add_argument("-r", "--results_directories", nargs='+', required=True,
                        help='''Results directory from a Canvas run. This
                        program requires the 'VisualizationTemp...' directories
                        produced by Canvas and optionally the 'CNV.vcf.gz' VCF
                        file.''')
    parser.add_argument("-o", "--output_directory", required=True,
                        help='''Directory to output plot images. Required.''')
    parser.add_argument("-p", "--ped", help='''Optional PED file for plotting
                        samples in the same family next to each other with
                        affected individuals above unaffected individuals.''')
    parser.add_argument("-d", "--dq", type=float,
                        help='''Only output variants from the CNV.vcf.gz in the
                        results directory with a de novo quality (DQ) equal to
                        or greater than this value.''')
    parser.add_argument("--pass_filters", action='store_true', default=False,
                        help='''When outputting variant regions only output
                        those which PASS filters.''')
    parser.add_argument("-v", "--variants", action='store_true',
                        help='''Output regions with a variant (as shown in the
                        CNV.vcf.gz file in the results directory) rather than
                        whole chromosomes. Use --vcf argument to use a
                        different VCF file.''')
    reg_group.add_argument("--vcfs", nargs='+', help='''Output variant regions
                           in this/these VCF file(s).''')
    reg_group.add_argument("--regions", nargs='+', help='''Output these
                           regions. Each region should be in the format
                           "chr1:1000000-2000000".''')
    parser.add_argument("-b", "--build", default='hg38',
                        help='''Genome build for plotting cytobands. A
                        cytoBandIdeo.txt.gz (as downloaded from
                        http://genome.ucsc.edu) must be present in the
                        'data/<build>/' subdirectory or no ideograms will be
                        drawn. Default=hg38.''')
    parser.add_argument("-z", "--zygosity", metavar='VCF',
                        help='''Also plot heterozygosity for each sample using
                        variants in this VCF file.''')
    parser.add_argument("--height", type=float,default=12,
                        help="Figure height in inches. Default=12.")
    parser.add_argument("--width", type=float, default=18,
                        help="Figure width in inches. Default=18.")
    parser.add_argument("--context", default='paper',
                        help='''Figure context to use. Can be set to any of the
                        seaborn plotting contexts ("notebook", "paper", "talk",
                        or "poster"). Default="paper".''')
    parser.add_argument("--ymax", type=float, default=6.0,
                        help='''Maximium limit of the y-scale for plots.
                        Default=6.0.''')
    parser.add_argument("-s", "--separate_chromosomes", action='store_true',
                        help='''If memory is limiting or processing many
                        samples at once, use this flag to only read in one
                        chromosome worth of data at a time (resulting in longer
                        runtimes).''')
    return parser

def process_data(results_directories, output_directory, cyto, samples,
                 variants, vcfs, pass_filters, dq, fig_dimensions, context,
                 ymax, chrom=None, het_vcf=None, regions=[]):
    df = pd.DataFrame()
    for results_directory in results_directories:
        df = pd.concat((df, read_coverage(results_directory, samples, chrom)))
    if variants or vcfs or dq:
        if not vcfs:
            vcfs = [os.path.join(rd, "CNV.vcf.gz") for rd in
                    results_directories]
        for vcf in vcfs:
            plot_variants(df=df, vcf=vcf, outdir=output_directory,
                          samples=samples, ideo=cyto, ymax=ymax, dq=dq,
                          fig_dimensions=fig_dimensions, chrom=chrom,
                          het_vcf=het_vcf)
    elif regions:
        process_regions(regions, df=df, outdir=output_directory, ideo=cyto,
                        samples=samples, ymax=ymax,
                        fig_dimensions=fig_dimensions, het_vcf=het_vcf)
    else:
        plot_chromosomes(df=df, outdir=output_directory, ideo=cyto,
                         samples=samples, fig_dimensions=fig_dimensions,
                         ymax=ymax, het_vcf=het_vcf)


def process_regions(regions, df, outdir, samples, ymax, fig_dimensions,
                    het_vcf, ideo=None, window_length=1e5):
    samples = check_samples_in_df(df, samples)
    intervals = []
    for region in regions:
        match = reg_matcher.match(region)
        if not match:
            logger.warn("Skipping invalid region: {}".format(region))
            continue
        chrom = match.group(1)
        start = int(match.group(3))
        end = int(match.group(4))
        intervals.append((chrom, start, end))
    intervals.sort(key=operator.itemgetter(0, 1, 2))
    prev_chrom = None
    het_df = None
    for region in intervals:
        if prev_chrom is None or region[0] != prev_chrom:
            chr_df = df[df.Chrom == region[0]]
            chr_ideo = ideo[(ideo.chrom == region[0])]
            if het_vcf is not None:
                logger.info("Getting heterzoygosity counts for chromosome " +
                            "{}".format(region[0]))
                gt_counts = het_counts.get_gt_counts(het_vcf, samples,
                                                     region[0],
                                                     logger=logger, plot=False)
                het_df = het_counts.counts_to_df(gt_counts, window_length)
        region_name = "{}_{}_{}".format(region[0], region[1], region[2])
        plot_region(df=chr_df, chrom=region[0], start=region[1], end=region[2],
                    outdir=outdir, ideo=ideo, samples=samples,
                    name=region_name, ymax=ymax, fig_dimensions=fig_dimensions,
                    het_df=het_df)

def get_chroms(results_dirs):
    logger.info("Checking chromosomes in coverage files...")
    have_chr = False
    no_chr = False
    chroms = set()
    for rdir in results_dirs:
        for vis_tmp in [d for d in os.listdir(rdir) if
                        d.startswith('VisualizationTemp')]:
            o_func = open
            f = os.path.join(rdir, vis_tmp, "coverage.bedgraph")
            if not os.path.exists(f):
                if os.path.exists(f + '.gz'):
                    f = f + '.gz'
                    o_func = gzip.open
                else:
                    sys.exit("ERROR: No coverage.bedgraph file in " +
                             "{}".format(vis_tmp))
                logger.info("Checking {}".format(f))
                with o_func(f, 'rt') as infile:
                    for line in infile:
                        c = line.split()[0]
                        if c.startswith("chr"):
                            have_chr = True
                        else:
                            no_chr = True
                        if have_chr and no_chr:
                            sys.exit("ERROR: mixed chromosome types in coverage " +
                                     "files - some chromosomes begin 'chr' while" +
                                     " others do not. Please only use with data " +
                                     "from the same reference.")
                        chroms.add(c)
    return chroms


def main(results_directories, output_directory, ped=None, variants=False,
         vcfs=[], pass_filters=False, dq=None, zygosity=None, height=18,
         width=12, context='paper', build='hg38', ymax=6.0, regions=[],
         separate_chromosomes=False):
    if os.path.exists(output_directory):
        logger.info("Using existing directory '{}'".format(output_directory))
    else:
        logger.info("Creating output directory '{}'".format(output_directory))
        os.mkdir(output_directory)
    sns.set_context(context)
    fig_dimensions = (width, height)
    samples = None
    cyto = get_cytobands(build)
    if ped is not None:
        logger.info("Getting samples from PED")
        samples = sample_order_from_ped(ped)
    if separate_chromosomes:
        for c in sorted(get_chroms(results_directories)):
                process_data(results_directories=results_directories,
                             output_directory=output_directory,
                             cyto=cyto,
                             samples=samples,
                             variants=variants,
                             vcfs=vcfs,
                             pass_filters=pass_filters,
                             dq=dq,
                             fig_dimensions=fig_dimensions,
                             context=context,
                             ymax=ymax,
                             het_vcf=zygosity,
                             chrom=c,
                             regions=regions)
    else:
        process_data(results_directories=results_directories,
                     output_directory=output_directory,
                     cyto=cyto,
                     samples=samples,
                     variants=variants,
                     vcfs=vcfs,
                     pass_filters=pass_filters,
                     dq=dq,
                     fig_dimensions=fig_dimensions,
                     context=context,
                     ymax=ymax,
                     het_vcf=zygosity,
                     regions=regions)

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
