#!/usr/bin/env python

"""
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import cPickle as pickle
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns
from matplotlib.pyplot import cm

from peppy import Project
from ngs_toolkit.rnaseq import RNASeqAnalysis, knockout_plot
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_analysis, collect_differential_analysis,
                                 differential_enrichment, differential_overlap,
                                 plot_differential,
                                 plot_differential_enrichment,
                                 unsupervised_analysis,
                                 least_squares_fit,
                                 subtract_principal_component)


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def main():
    # Start project
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj._samples = [s for s in prj.samples if (
        (s.library == "RNA-seq") &
        (s.cell_line in ["HAP1"]) and
        (s.to_use != "0") and
        (s.main_analysis == "1"))]
    for sample in prj.samples:
        if sample.library in ["RNA-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome), "bitSeq", sample.name + ".counts")

    # RNA-seq
    data_type = "RNA-seq"
    quant_matrix = "expression_annotated"
    feature_name = "genes"
    rnaseq_analysis = RNASeqAnalysis(name="baf_complex.rnaseq", prj=prj, samples=prj.samples)


    # Get and normalize expression data
    rnaseq_analysis.get_gene_expression(sample_attributes=prj.sample_attributes)
    rnaseq_analysis.to_pickle()


    # # Unsupervised analysis
    # unsupervised_analysis(
    #     rnaseq_analysis, quant_matrix=quant_matrix, samples=None,
    #     attributes_to_plot=prj.group_attributes,
    #     plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
    #     axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
    #     output_dir="{}/unsupervised_analysis_{}".format(rnaseq_analysis.results_dir, data_type))


    # Fix batch effect
    import patsy
    import pandas as pd
    import rpy2
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    from rpy2.rinterface import RRuntimeError
    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r('require("limma")')
    _removeBatchEffect = robjects.r('removeBatchEffect')
    _as_formula = robjects.r('as.formula')
    _DESeq = robjects.r('DESeq')
    _results = robjects.r('results')
    _as_data_frame = robjects.r('as.data.frame')

    fixed = _removeBatchEffect(
        x=rnaseq_analysis.expression_annotated.values,
        batch=matrix.columns.get_level_values("batch"),
        design=patsy.dmatrix("~knockout - 1", matrix.columns.to_frame()))
    rnaseq_analysis.limma_fixed = pd.DataFrame(np.asarray(fixed), index=matrix.index, columns=matrix.columns)

    unsupervised_analysis(
        rnaseq_analysis, quant_matrix='limma_fixed', plot_prefix="limma_fixed",
        samples=None, attributes_to_plot=prj.group_attributes,
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
        axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
        output_dir="{}/unsupervised_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        test_pc_association=False, standardize_matrix=False)

    # Supervised analysis
    alpha = 0.05
    abs_fold_change = 0

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['cell_type'] == "HAP1") &
        (comparison_table['data_type'] == data_type) &
        (comparison_table['comparison_type'] == 'differential')]
    rnaseq_analysis.differential_results = differential_analysis(
        rnaseq_analysis,
        comparison_table,
        data_type=data_type,
        samples=[s for s in rnaseq_analysis.samples if s.name in comparison_table['sample_name'].tolist()],
        output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        covariates=None,
        alpha=alpha,
        overwrite=True, distributed=True)
    rnaseq_analysis.differential_results = collect_differential_analysis(
        comparison_table,
        data_type=data_type,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=False,
        overwrite=True)
    rnaseq_analysis.to_pickle()

    # Knockout plot
    knockout_plot(
        rnaseq_analysis,
        expression_matrix=rnaseq_analysis.expression_annotated,
        comparison_results=rnaseq_analysis.differential_results)

    # Overlap between differential regions
    differential_overlap(
        rnaseq_analysis.differential_results[
            (rnaseq_analysis.differential_results['padj'] < alpha) &
            (rnaseq_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        getattr(rnaseq_analysis, quant_matrix).shape[0],
        output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        data_type=data_type)

    quant_matrix = 'limma_fixed'
    plot_differential(
        rnaseq_analysis,
        rnaseq_analysis.differential_results, # rnaseq_analysis.differential_results[~rnaseq_analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")], 
        matrix=getattr(rnaseq_analysis, quant_matrix),
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        output_prefix="differential_analysis_" + quant_matrix,
        data_type=data_type,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
        group_variables=prj.group_attributes)

    # Enrichment
    differential_enrichment(
        rnaseq_analysis,
        rnaseq_analysis.differential_results[
            (rnaseq_analysis.differential_results['padj'] < alpha) &
            (rnaseq_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        genome="hg19",
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True)

    collect_differential_enrichment(
        rnaseq_analysis.differential_results[
            (rnaseq_analysis.differential_results['padj'] < alpha) &
            (rnaseq_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
        permissive=False)

    for enrichment_name, enrichment_type in [
        ('enrichr', 'enrichr')]:
        try:
            enrichment_table = pd.read_csv(
                os.path.join("{}/differential_analysis_{}".format(
                    rnaseq_analysis.results_dir, data_type), "differential_analysis" + ".{}.csv".format(enrichment_type)))
        except pd.errors.EmptyDataError:
            print("Enrichment dataframe of {} is empty.".format(enrichment_type))
            continue

        plot_differential_enrichment(
            enrichment_table,
            enrichment_name,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, data_type),
            output_prefix="differential_analysis",
            direction_dependent=True,
            barplots=False, correlation_plots=False,
            top_n=5 if enrichment_name != "motif" else 300)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
