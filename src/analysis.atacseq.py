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
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_analysis, collect_differential_analysis,
                                 differential_enrichment, differential_overlap,
                                 plot_differential,
                                 plot_differential_enrichment,
                                 least_squares_fit,
                                 subtract_principal_component,
                                 unsupervised_analysis)

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
        (s.library == "ATAC-seq") &
        (s.cell_line in ["HAP1"]) and
        (s.to_use != "0") and
        (s.main_analysis == "1"))]
    for sample in prj.samples:
        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome), "bitSeq", sample.name + ".counts")

    # ATAC-seq
    data_type = "ATAC-seq"
    quant_matrix = "accessibility"
    feature_name = "sites"
    atac_analysis = ATACSeqAnalysis(name="baf_complex.atacseq", prj=prj, samples=prj.samples)


    # Usual pipeline
    atac_analysis.get_consensus_sites(region_type="summits", extension=250, blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed")
    atac_analysis.measure_coverage()
    atac_analysis.normalize(method="gc_content")
    atac_analysis.annotate_with_sample_metadata(quant_matrix="coverage_gc_corrected", attributes=prj.sample_attributes)
    atac_analysis.to_pickle()


    # atac_analysis.coverage = atac_analysis.coverage.loc[:, [s.name for s in prj.samples]]
    # atac_analysis.coverage_qnorm = atac_analysis.coverage_qnorm.loc[:, [s.name for s in prj.samples]]
    # atac_analysis.coverage_gc_corrected = atac_analysis.coverage_gc_corrected.loc[:, [s.name for s in prj.samples]]
    # atac_analysis.annotate_with_sample_metadata(quant_matrix="coverage_gc_corrected", attributes=prj.sample_attributes)


    # Unsupervised analysis
    unsupervised_analysis(
        atac_analysis, quant_matrix=quant_matrix, samples=None,
        attributes_to_plot=prj.group_attributes,
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
        axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
        output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type))


    # Fix batch effect
    import pandas as pd
    import rpy2
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    from rpy2.rinterface import RRuntimeError
    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r('require("limma")')
    _removeBatchEffect = robjects.r('removeBatchEffect')

    fixed = _removeBatchEffect(
        x=atac_analysis.accessibility.values,
        batch=matrix.columns.get_level_values("batch"),
        design=patsy.dmatrix("~knockout - 1", matrix.columns.to_frame()))
    atac_analysis.limma_fixed = pd.DataFrame(np.asarray(fixed), index=matrix.index, columns=matrix.columns)

    unsupervised_analysis(
        atac_analysis, quant_matrix='limma_fixed', plot_prefix="limma_fixed",
        samples=None, attributes_to_plot=prj.group_attributes,
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
        axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
        output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type),
        test_pc_association=False, standardize_matrix=False)

    # matrix = atac_analysis.accessibility
    # standardize_data = True

    # fits = least_squares_fit(
    #     quant_matrix=matrix.T,
    #     design_matrix=matrix.columns.to_frame(),
    #     test_model="~ batch -1",
    #     standardize_data=standardize_data)

    # batches = matrix.columns.get_level_values("batch").unique()
    # fixed = pd.DataFrame(np.zeros(matrix.shape), index=matrix.index, columns=matrix.columns)
    # i_fixed = pd.DataFrame(np.zeros(matrix.shape), index=matrix.index, columns=matrix.columns)
    # for batch in batches:
    #     # c = fits['Intercept'] if batch == batches[0] else fits['batch[T.{}]'.format(batch)]
    #     c = fits['batch[{}]'.format(batch)]
    #     o = matrix.loc[:, matrix.columns.get_level_values("batch") == batch].T
    #     fixed.loc[:, matrix.columns.get_level_values("batch") == batch] = (o - c).T
    #     i_fixed.loc[:, matrix.columns.get_level_values("batch") == batch] = (o + (c * 2)).T
    # atac_analysis.fixed = fixed
    # atac_analysis.fixed_rescaled = ((fixed.T + matrix.mean(axis=1)) * matrix.std(axis=1)).T
    # atac_analysis.i_fixed = i_fixed
    # atac_analysis.i_fixed_rescaled = ((i_fixed.T + matrix.mean(axis=1)) * matrix.std(axis=1)).T

    # unsupervised_analysis(
    #     atac_analysis, quant_matrix='fixed', plot_prefix="fixed",
    #     samples=None, attributes_to_plot=prj.group_attributes,
    #     plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
    #     axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
    #     output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type),
    #     test_pc_association=False, standardize_matrix=True)
    # unsupervised_analysis(
    #     atac_analysis, quant_matrix='fixed_rescaled', plot_prefix="fixed_rescaled",
    #     samples=None, attributes_to_plot=prj.group_attributes,
    #     plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
    #     axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
    #     output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type),
    #     test_pc_association=False, standardize_matrix=True)
    # unsupervised_analysis(
    #     atac_analysis, quant_matrix='i_fixed_rescaled', plot_prefix="i_fixed_rescaled",
    #     samples=None, attributes_to_plot=prj.group_attributes,
    #     plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
    #     axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
    #     output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type),
    #     test_pc_association=False)

    # f = subtract_principal_component(
    #     X=atac_analysis.i_fixed_rescaled.T, pc=1, norm=False, plot=False)
    # atac_analysis.i_fixed_rescaled_pcafixed = pd.DataFrame(f.T, index=matrix.index, columns=matrix.columns)

    # # f = subtract_principal_component(
    # #     X=atac_analysis.i_fixed_rescaled_pcafixed.T, pc=1, norm=False, plot=False)
    # # atac_analysis.i_fixed_rescaled_pcafixed = pd.DataFrame(f.T, index=matrix.index, columns=matrix.columns)

    # unsupervised_analysis(
    #     atac_analysis, quant_matrix='i_fixed_rescaled_pcafixed', plot_prefix="i_fixed_rescaled_pcafixed",
    #     samples=None, attributes_to_plot=prj.group_attributes,
    #     plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True,
    #     axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
    #     output_dir="{}/unsupervised_analysis_{}".format(atac_analysis.results_dir, data_type),
    #     test_pc_association=False,
    #     standardize_matrix=True)

    # Knockout plot
    from ngs_toolkit.rnaseq import knockout_plot
    atac_g = atac_analysis.get_gene_level_accessibility()
    knockout_plot(atac_analysis, expression_matrix=atac_g, output_prefix="knockout_expression.atacseq_gene_level")


    # Supervised analysis
    alpha = 0.01
    abs_fold_change = 1

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['data_type'] == data_type) &
        (comparison_table['comparison_type'] == 'differential')]
    atac_analysis.differential_results = differential_analysis(
        atac_analysis,
        comparison_table,
        data_type=data_type,
        samples=[s for s in atac_analysis.samples if s.name in comparison_table['sample_name'].tolist()],
        output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
        covariates=None,
        alpha=alpha,
        overwrite=True, distributed=True)
    atac_analysis.differential_results = collect_differential_analysis(
        comparison_table,
        data_type=data_type,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=False,
        overwrite=True)
    atac_analysis.to_pickle()

    # Overlap between differential regions
    differential_overlap(
        atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        getattr(atac_analysis, quant_matrix).shape[0],
        output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
        data_type=data_type)

    quant_matrix = 'limma_fixed'
    plot_differential(
        atac_analysis,
        atac_analysis.differential_results, # atac_analysis.differential_results[~atac_analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")], 
        matrix=getattr(atac_analysis, quant_matrix),
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
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
    # first annotate peaks with closest gene, genomic regions, ChromHMM state from CD19+ cells
    atac_analysis.get_peak_gene_annotation()
    atac_analysis.get_peak_genomic_location()
    atac_analysis.get_peak_chromatin_state(
        chrom_state_file="data/external/HAP1_12_segments.annotated.bed")
    atac_analysis.annotate(quant_matrix='coverage_gc_corrected')
    atac_analysis.to_pickle()

    differential_enrichment(
        atac_analysis,
        atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
        genome="hg19",
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True)

    collect_differential_enrichment(
        atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
        permissive=True)

    for enrichment_name, enrichment_type in [
        ('motif', 'meme_ame'),
        ('lola', 'lola'),
        ('enrichr', 'enrichr')]:
        try:
            enrichment_table = pd.read_csv(
                os.path.join("{}/differential_analysis_{}".format(
                    atac_analysis.results_dir, data_type), "differential_analysis" + ".{}.csv".format(enrichment_type)))
        except pd.errors.EmptyDataError:
            print("Enrichment dataframe of {} is empty.".format(enrichment_type))
            continue

        plot_differential_enrichment(
            enrichment_table,
            enrichment_name,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(atac_analysis.results_dir, data_type),
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
