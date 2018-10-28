#!/bin/env python

import os
from peppy import Project
import pandas as pd
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from ngs_toolkit.general import unsupervised_analysis


# Set settings
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def main():
    prj = Project(os.path.join("metadata", "project_config.yaml"))    
    prj._samples = [s for s in prj.samples if (s.library == "ATAC-seq")]
    resolutions = ["1000kb", "100kb", "20kb", "10kb"]
    data_dir = "/home/arendeiro/baf_complex/copywriter/processed/"


    # from ngs_toolkit.cnv import CNVAnalysis
    analysis = CNVAnalysis(
        name="baf_complex.cnv", prj=prj, data_dir=data_dir,
        resolutions=resolutions)

    # Get data
    c = analysis.get_cnv_data(resolutions=resolutions, save=True, assign=True, permissive=True)

    # Normalize and observe changes
    analysis.normalize(method="median", resolutions=resolutions)
    analysis.plot_all_data(sample_labels=True, rasterized=True, output_prefix=analysis.name + ".all_data.median")
    analysis.to_pickle()

    # only HAP1 now
    resolutions = ["1000kb", "100kb", "20kb", "10kb"]
    analysis = CNVAnalysis(
        name="baf_complex.hap1.cnv", prj=prj, data_dir=data_dir,
        samples=sorted([s for s in prj.samples if s.cell_line == "HAP1"], key=lambda x: x['sample_name']),
        resolutions=resolutions)
    # Get data
    c = analysis.get_cnv_data(save=True, assign=True, permissive=True)

    # Normalize and observe changes
    analysis.normalize(method="median")
    analysis.plot_all_data(sample_labels=True, robust=False, rasterized=True, output_prefix=analysis.name + ".all_data.median")
    analysis.to_pickle()

    for resolution in resolutions:
        to_igv(
            analysis.coverage_norm[resolution],
            output_file=os.path.join(analysis.results_dir, analysis.name + ".{}.igv".format(resolution)),
            save=True, viewLimits=[-2, 2])

    # Plot stats
    analysis.plot_stats_per_chromosome(sample_labels=True, resolutions=resolutions)

    # Unsupervised analysis
    for resolution in resolutions:
        attrs = ["sample_name", "knockout", "clone", "replicate"]
        analysis.cnv = analysis.coverage_norm[resolution]
        analysis.cnv = analysis.annotate_with_sample_metadata(quant_matrix="cnv", attributes=attrs)

        ## bulk
        unsupervised_analysis(
            analysis, quant_matrix="cnv", data_type="CNV",
            output_dir="{results_dir}",
            plot_prefix=analysis.name,
            attributes_to_plot=attrs[1:])

    # Not run:
    # # Segment genome
    # analysis.segment_genome()
    # # Annotate segments
    # analysis.annotate_with_chrom_bands(resolutions)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
