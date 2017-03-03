import matplotlib
matplotlib.use('Agg')
import argparse
import pandas as pd
import numpy as np
import pybedtools
import pysam
import tabix
import os


def cuts_in_window(bam_file, intervals, n=1001):
    """
    Counts the start position of reads in a window from a bam file.
    """
    counts = np.zeros((len(intervals), len(range(n))))

    bam = pysam.AlignmentFile(bam_file)

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    for i, interval in enumerate(intervals):
        if interval.split(":")[0] in chroms:
            for aln in bam.fetch(region=interval):
                pos_in_window = aln.reference_start - int(interval.split(":")[1].split("-")[0])
                if 0 <= pos_in_window < n:
                    counts[i, pos_in_window] += 1
    bam.close()

    return pd.DataFrame(counts, index=intervals, columns=np.arange(0 - (n / 2), (n / 2) + 1))


def signal_in_window(bedgraph_file, intervals, n=1001):
    """
    Counts signal from a bedgraph (tabix indexed) file within a window.
    """
    counts = np.zeros((len(intervals), len(range(n))))

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    tb = tabix.open(bedgraph_file)

    for i, interval in enumerate(intervals):
        if interval.split(":")[0] in chroms:
            for record in tb.querys(interval):
                pos_in_window = int(record[1]) - int(interval.split(":")[1].split("-")[0])
                if 0 <= pos_in_window < n:
                    counts[i, pos_in_window] = record[3]

    return pd.DataFrame(counts, index=intervals, columns=np.arange(0 - (n / 2), (n / 2) + 1))


def fragments_in_window(bam_file, intervals, n=1001):
    """
    Counts the midpoint of paired-end fragments in a window from a bam file.
    """
    counts = np.zeros((len(intervals), len(range(n))))

    bam = pysam.AlignmentFile(bam_file)

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    for i, interval in enumerate(intervals):
        if interval.split(":")[0] in chroms:
            for aln in bam.fetch(region=interval):
                if not aln.is_reverse:
                    midpoint = aln.reference_start + (aln.template_length / 2)
                    pos_in_window = midpoint - int(interval.split(":")[1].split("-")[0])
                    if 0 <= pos_in_window < n:
                        counts[i, pos_in_window] += 1
    bam.close()

    return pd.DataFrame(counts, index=intervals, columns=np.arange(0 - (n / 2), (n / 2) + 1))


def bed_in_window(bed_file, intervals, n=1001):
    """
    Counts the number of region start point in a window from a bam file (for dyad positions).
    """
    counts = np.zeros((len(intervals), len(range(n))))

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    tb = tabix.open(bed_file)

    for i, interval in enumerate(intervals):
        if interval.split(":")[0] in chroms:
            try:
                for record in tb.querys(interval):
                    pos_in_window = int(record[1]) - int(interval.split(":")[1].split("-")[0])
                    if 0 <= pos_in_window < n:
                        counts[i, pos_in_window] += 1
            except tabix.TabixError:
                pass

    return pd.DataFrame(counts, index=intervals, columns=np.arange(0 - (n / 2), (n / 2) + 1))


parser = argparse.ArgumentParser()

parser.add_argument("--bed-file", dest="bed_file", required=True)
parser.add_argument("--bam-file", dest="bam_file", required=True)
parser.add_argument("--coverage-type", dest="coverage_type", required=True)
parser.add_argument("--name", dest="name", required=True)
parser.add_argument("--output-dir", dest="output_dir", required=True)
parser.add_argument("--plot", dest="plot", default=False, action="store_true")

args = parser.parse_args()

bedtool = pybedtools.BedTool(args.bed_file)

if args.coverage_type == "signal":
    coverage_function = cuts_in_window
elif args.coverage_type == "nucleoatac":
    coverage_function = signal_in_window
elif args.coverage_type == "dyads":
    coverage_function = bed_in_window
else:
    coverage_function = fragments_in_window

coverage_matrix = coverage_function(
    args.bam_file,
    [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in bedtool]
)
coverage_matrix.to_csv(os.path.join(args.output_dir, "%s.%s.coverage_matrix.csv" % (args.name, args.coverage_type)))

# coverage_matrix = pd.read_csv(os.path.join(args.output_dir, "%s.coverage_matrix.csv" % args.name), index_col=0)

if args.plot:
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set settings
    sns.set(context="paper", style="white", palette="pastel", color_codes=True)
    sns.set_palette(sns.color_palette("colorblind"))
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rc('text', usetex=False)

    # Plot profile
    fig, axis = plt.subplots(1)
    axis.plot(np.arange(-500, 501), coverage_matrix.mean(0))
    sns.despine(fig)
    fig.savefig(os.path.join(args.output_dir, "%s.mean.png" % args.name), bbox_inches="tight")

    # Plot image
    fig, axis = plt.subplots(1)
    axis.imshow(coverage_matrix, aspect="auto", interpolation="none")
    fig.savefig(os.path.join(args.output_dir, "%s.imshow.png" % args.name), bbox_inches="tight")
