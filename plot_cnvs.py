#!/usr/bin/env python3
import sys
import os
import re
import logging
import argparse
import matplotlib
if os.environ.get('DISPLAY','') == '':
    matplotlib.use('Agg') #for headless servers
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.collections import BrokenBarHCollection
import seaborn as sns
import pandas as pd
from vase.ped_file import PedFile

logger = logging.getLogger("CNV plotter")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
       '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

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

def read_coverage(results_dir, samples=None):
    df = pd.DataFrame()
    logger.info("Reading coverage files")
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
        temp_df['Sample'] = sample
        df = pd.concat([df, temp_df])
    return df

def plot_samples(df, outdir, ideo, samples, fig_dimensions, ymax=6):
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
    cols = sns.color_palette("colorblind", len(samples))
    for chrom in df.Chrom.unique():
        logger.info("Plotting chromosome {}".format(chrom))
        png = os.path.join(outdir, chrom + '.png')
        chr_df = df[df.Chrom == chrom]
        i = 0
        fig, axes = plt.subplots(nrows=1 + len(samples),
                                 ncols=1,
                                 sharex=True,
                                 gridspec_kw={'height_ratios':[1] + [8] *
                                                              len(samples),
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
            tmp_df = chr_df[chr_df.Sample == s]
            alpha = min(50000.0/len(tmp_df), 0.5)
            axes[i + 1].scatter(tmp_df.Start, tmp_df.Ploidy, color=cols[i],
                        alpha=alpha, s=[0.1] * len(tmp_df))
            axes[i + 1].plot([chr_df.Start.min(), chr_df.Start.max()],
                             [2.0, 2.0 ],
                             '--',
                             color='black',
                             alpha=0.4)
            axes[i + 1].axes.set_xlim((0, max(chr_df.End)))
            axes[i + 1].axes.set_ylim((-0.5, ymax))
            axes[i + 1].axes.set_title(s)
            axes[i + 1].axes.set_ylabel('CN')
            i += 1
            if i == len(samples):#add label for bottom plot
                axes[i].axes.set_xlabel('Pos')
                #make the position labels a bit prettier
                axes[i].xaxis.set_major_formatter(
                    mtick.FuncFormatter(lambda x, p: format(int(x), ',')))


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
    parser.add_argument("-r", "--results_directory", required=True,
                        help='''One or more results directories from a Canvas
                        run. This program requires the 'VisualizationTemp...'
                        directories produced by Canvas.''')
    parser.add_argument("-o", "--output_directory", required=True,
                        help='''Directory to output plot images. Required.''')
    parser.add_argument("-p", "--ped", help='''Optional PED file for plotting
                        samples in the same family next to each other with
                        affected individuals above unaffected individuals.''')
    parser.add_argument("-b", "--build", default='hg38',
                        help='''Genome build for plotting cytobands. A
                        cytoBandIdeo.txt.gz (as downloaded from
                        http://genome.ucsc.edu) must be present in the
                        'data/<build>/' subdirectory or no ideograms will be
                        drawn. Default=hg38.''')
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
    return parser

def main(results_directory, output_directory, ped=None, height=18, width=12,
         context='paper', build='hg38', ymax=6.0):
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
    df = read_coverage(results_directory, samples)
    plot_samples(df, output_directory, cyto, samples, fig_dimensions, ymax)

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
