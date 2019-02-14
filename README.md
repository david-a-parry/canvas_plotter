
# Canvas Plotter

Plots CNV data from [Canvas](https://github.com/Illumina/canvas) output.

Plots information from bedgraphs produced by Canvas and optionally can plot
heterozygosity from a VCF with SNV calls for the same samples.

## Install

### Dependencies

Python 3 is required. This requires some modules from vase
(https://github.com/gantzgraf/vase) in order to run. You can install vase via
pip3 and git:

    pip3 install git+git://github.com/gantzgraf/vase.git --user

Full installation instructions for vase are on the [vase github page](https://github.com/gantzgraf/vase.git)

### Installation

You need to clone this project and the het_plotter submodule. The easiest way to
do this is as follows:

    git clone --recursive https://git.ecdf.ed.ac.uk/dparry/canvas_plotter.git

This will clone the project into a directory called canvas_plotter. Run the main
script with the --help option to see usage information:

    canvas_plotter/plot_cnvs.py --help
    #or "python3 canvas_plotter/plot_cnvs.py --help" if you get a permissions error

## Examples

To simply plot each chromosome from a Canvas run:

    canvas_plotter/plot_cnvs.py \
        -r /path/to/canvas/results/ \
        -o output_directory

If you have affected individuals and unaffected individuals you can specify a 
PED file to ensure affecteds are plotted on top of unaffecteds:

    canvas_plotter/plot_cnvs.py \
        --ped /path/to/family.ped
        -r /path/to/canvas/results/ \
        -o output_directory

To also plot heterozygosity:

    canvas_plotter/plot_cnvs.py \
        -z /path/to/snvs.vcf.gz \
        -r /path/to/canvas/results/ \
        -o output_directory

To plot a specific region:

    canvas_plotter/plot_cnvs.py \
        --region chr1:1000000-2000000 \
        -r /path/to/canvas/results/ \
        -o output_directory

To plot variants called by Canvas instead of whole chromosomes:

    canvas_plotter/plot_cnvs.py \
        --variants \
        -r /path/to/canvas/results/ \
        -o output_directory

To plot SVs/CNVs called in another VCF:

    canvas_plotter/plot_cnvs.py \
        --vcfs /path/to/another_caller.vcf \
        -r /path/to/canvas/results/ \
        -o output_directory

See below for full usage. 

## Usage

~~~~

usage: plot_cnvs.py [-h] -r RESULTS_DIRECTORIES [RESULTS_DIRECTORIES ...] -o
                    OUTPUT_DIRECTORY [-p PED] [-d DQ] [--pass_filters] [-v]
                    [--vcfs VCFS [VCFS ...] | --regions REGIONS [REGIONS ...]]
                    [-b BUILD] [-z VCF] [--height HEIGHT] [--width WIDTH]
                    [--context CONTEXT] [--ymax YMAX] [-s]

Plot Canvas CNV results.

optional arguments:
  -h, --help            show this help message and exit
  -r RESULTS_DIRECTORIES [RESULTS_DIRECTORIES ...], --results_directories RESULTS_DIRECTORIES [RESULTS_DIRECTORIES ...]
                        Results directory from a Canvas run. This program
                        requires the 'VisualizationTemp...' directories
                        produced by Canvas and optionally the 'CNV.vcf.gz' VCF
                        file.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Directory to output plot images. Required.
  -p PED, --ped PED     Optional PED file for plotting samples in the same
                        family next to each other with affected individuals
                        above unaffected individuals.
  -d DQ, --dq DQ        Only output variants from the CNV.vcf.gz in the
                        results directory with a de novo quality (DQ) equal to
                        or greater than this value.
  --pass_filters        When outputting variant regions only output those
                        which PASS filters.
  -v, --variants        Output regions with a variant (as shown in the
                        CNV.vcf.gz file in the results directory) rather than
                        whole chromosomes. Use --vcf argument to use a
                        different VCF file.
  --vcfs VCFS [VCFS ...]
                        Output variant regions in this/these VCF file(s).
  --regions REGIONS [REGIONS ...]
                        Output these regions. Each region should be in the
                        format "chr1:1000000-2000000".
  -b BUILD, --build BUILD
                        Genome build for plotting cytobands. A
                        cytoBandIdeo.txt.gz (as downloaded from
                        http://genome.ucsc.edu) must be present in the
                        'data/<build>/' subdirectory or no ideograms will be
                        drawn. Default=hg38.
  -z VCF, --zygosity VCF
                        Also plot heterozygosity for each sample using
                        variants in this VCF file.
  --height HEIGHT       Figure height in inches. Default=12.
  --width WIDTH         Figure width in inches. Default=18.
  --context CONTEXT     Figure context to use. Can be set to any of the
                        seaborn plotting contexts ("notebook", "paper",
                        "talk", or "poster"). Default="paper".
  --ymax YMAX           Maximium limit of the y-scale for plots. Default=6.0.
  -s, --separate_chromosomes
                        If memory is limiting or processing many samples at
                        once, use this flag to only read in one chromosome
                        worth of data at a time (resulting in longer
                        runtimes).
~~~~
