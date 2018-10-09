#!/usr/bin/env python

"""
This is the main script of the baf_complex project.
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

from looper.models import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.chipseq import ChIPSeqAnalysis, homer_peaks_to_bed
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_analysis,
                                 differential_enrichment, differential_overlap,
                                 plot_differential,
                                 plot_differential_enrichment,
                                 unsupervised_analysis)
from ngs_toolkit.rnaseq import RNASeqAnalysis, knockout_plot

# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def main():
    # Start project
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj._samples = [s for s in prj.samples if s.to_use == "1"]
    for sample in prj.samples:
        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome), "bitSeq", sample.name + ".counts")

    # Sample's attributes
    sample_attributes = ['sample_name', 'cell_line', 'knockout', 'treatment', 'replicate', 'clone', 'batch']
    plotting_attributes = ['knockout', 'treatment', 'replicate', 'clone', 'batch']

    # HAP1 ANALYSIS
    # ATAC-seq
    atacseq_samples = [s for s in prj.samples if (s.library == "ATAC-seq") & (s.cell_line in ["HAP1"])]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    # atacseq_samples = [s for s in atacseq_samples if s.knockout != "SMARCC2" and s.clone != "GFP"]
    atac_analysis = ATACSeqAnalysis(name="baf_complex.atacseq", prj=prj, samples=atacseq_samples)
    atac_analysis = main_analysis_pipeline(atac_analysis, data_type="ATAC-seq", cell_type="HAP1")

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["HAP1"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    # rnaseq_samples = [s for s in rnaseq_samples if s.knockout != "SMARCC2" and s.clone != "GFP"]
    rnaseq_analysis = RNASeqAnalysis(name="baf-complex.rnaseq", prj=prj, samples=rnaseq_samples)
    rnaseq_analysis = main_analysis_pipeline(rnaseq_analysis, data_type="RNA-seq", cell_type="HAP1")

    # export tables
    res = atac_analysis.differential_results
    diff = res[(res['padj'] < 0.01) & (res['log2FoldChange'].abs() > 1)]
    q = atac_analysis.accessibility.loc[diff.index.unique(), :]
    q.to_csv(os.path.join("results", 'differential_analysis_ATAC-seq',
                            "differential_accessibility.regions.csv"))
    q = q.apply(scipy.stats.zscore, axis=1)
    q.to_csv(os.path.join("results", 'differential_analysis_ATAC-seq',
                          "differential_accessibility.regions.zscore.csv"))
    baf_genes = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv"), squeeze=True)
    chrom_list = pd.read_csv(os.path.join("metadata", "Bocklab_chromatin_genes.csv"))
    for metric in ['padj', 'log2FoldChange']:
        piv = pd.pivot_table(
            data=rnaseq_analysis.differential_results.loc[baf_genes].reset_index(),
            index='index', columns='comparison_name', values=metric)
        piv.to_csv(os.path.join("results", 'differential_analysis_RNA-seq', "baf_complex_members.differential_expression.{}.csv".format(metric)))
        piv = pd.pivot_table(
            data=rnaseq_analysis.differential_results.loc[chrom_list['HGNC_symbol']].reset_index(),
            index='index', columns='comparison_name', values=metric)
        piv.to_csv(os.path.join("results", 'differential_analysis_RNA-seq', "chromatin_proteins.differential_expression.{}.csv".format(metric)))


    # OTHER CELL LINES
    # OV90 ANALYSIS
    # ATAC-seq
    atacseq_samples = [s for s in prj.samples if (s.library == "ATAC-seq") & (s.cell_line in ["OV90"])]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    ov90_atac_analysis = ATACSeqAnalysis(name="baf_complex.ov90.atacseq", prj=prj, samples=atacseq_samples, results_dir="results_ov90")
    ov90_atac_analysis = main_analysis_pipeline(ov90_atac_analysis, data_type="ATAC-seq", cell_type="OV90")

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["OV90"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    ov90_rnaseq_analysis = RNASeqAnalysis(name="baf-complex.ov90.rnaseq", prj=prj, samples=rnaseq_samples, results_dir="results_ov90")
    ov90_rnaseq_analysis = main_analysis_pipeline(ov90_rnaseq_analysis, data_type="RNA-seq", cell_type="OV90")

    # A549 ANALYSIS
    # ATAC-seq
    atacseq_samples = [s for s in prj.samples if (s.library == "ATAC-seq") & (s.cell_line in ["A549"])]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    a549_atac_analysis = ATACSeqAnalysis(name="baf_complex.a549.atacseq", prj=prj, samples=atacseq_samples, results_dir="results_a549")
    a549_atac_analysis = main_analysis_pipeline(a549_atac_analysis, data_type="ATAC-seq", cell_type="A549")

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["A549"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    a549_rnaseq_analysis = RNASeqAnalysis(name="baf-complex.a549.rnaseq", prj=prj, samples=rnaseq_samples, results_dir="results_a549")
    a549_rnaseq_analysis = main_analysis_pipeline(a549_rnaseq_analysis, data_type="RNA-seq", cell_type="A549")

    # H2122 ANALYSIS
    # ATAC-seq
    atacseq_samples = [s for s in prj.samples if (s.library == "ATAC-seq") & (s.cell_line in ["H2122"])]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    h2122_atac_analysis = ATACSeqAnalysis(name="baf_complex.h2122.atacseq", prj=prj, samples=atacseq_samples, results_dir="results_h2122")
    h2122_atac_analysis = main_analysis_pipeline(h2122_atac_analysis, data_type="ATAC-seq", cell_type="H2122")

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["H2122"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    h2122_rnaseq_analysis = RNASeqAnalysis(name="baf-complex.h2122.rnaseq", prj=prj, samples=rnaseq_samples, results_dir="results_h2122")
    h2122_rnaseq_analysis = main_analysis_pipeline(h2122_rnaseq_analysis, data_type="RNA-seq", cell_type="H2122")


    # ChIP-seq data
    # Start project and analysis objects
    chipseq_samples = [s for s in prj.samples if s.library in ["ChIP-seq", "ChIPmentation"]]
    chipseq_samples = [s for s in chipseq_samples if os.path.exists(s.filtered)]
    chipseq_analysis = ChIPSeqAnalysis(
        name="baf_complex.chipseq", prj=prj, samples=chipseq_samples)

    # read in comparison table
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[comparison_table['comparison_type'] == 'peaks']

    # call peaks
    chipseq_analysis.call_peaks_from_comparisons(
        comparison_table=comparison_table, overwrite=False)
    # summarize peaks
    chipseq_analysis.peak_summary = summarize_peaks_from_comparisons(
        chipseq_analysis,
        comparison_table=comparison_table)
    chipseq_analysis.peak_summary.to_csv(
        os.path.join(chipseq_analysis.results_dir, chipseq_analysis.name + "_peak_summary.csv"), index=False)

    # First, set peak set to ATAC-seq and quantify all ChIP-seq samples in there
    chipseq_analysis.set_consensus_sites(atac_analysis.sites.fn)
    chipseq_analysis.measure_coverage()
    # c = chipseq_analysis.coverage.iloc[:, :-3]
    # fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 1 * 4))
    # mean = np.log2(1 + c.mean(axis=1))
    # qv2 = (c.std(axis=1) / mean) ** 2
    # axis[0].scatter(mean, c.std(axis=1), s=5, alpha=0.5)
    # axis[1].scatter(mean, qv2, s=5, alpha=0.5)
    chipseq_analysis.normalize()
    chipseq_analysis.annotate(quant_matrix="coverage_qnorm")
    chipseq_analysis.annotate_with_sample_metadata(attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    chipseq_analysis.to_pickle()

    # Unsupervised analysis
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=None,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_all_peaks",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")

    # without histone marks
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=[s for s in chipseq_analysis.samples if "H3" not in s.ip],
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_all_peaks.no_histones",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")


    diff_atac_heatmap_with_chip(atac_analysis, chipseq_analysis)


    # Now let's try to use the ChIP-seq peaks on their own
    chipseq_samples = [s for s in prj.samples if s.library in ["ChIP-seq", "ChIPmentation"]]
    chipseq_samples = [s for s in chipseq_samples if os.path.exists(s.filtered)]
    chipseq_analysis = ChIPSeqAnalysis(name="baf_complex.chipseq.peaks", prj=prj, samples=chipseq_samples)
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    c = comparison_table[
        (comparison_table['comparison_type'] == 'peaks') &
        (comparison_table['comparison_name'].str.contains("ARID|SMARC|PBRM"))]
    c['comparison_genome'] = 'hg19'
    chipseq_analysis.get_consensus_sites(
        comparison_table=c,
        region_type="peaks", blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed")
    chipseq_analysis.calculate_peak_support(comparison_table=c)
    chipseq_analysis.measure_coverage()
    chipseq_analysis.normalize()
    chipseq_analysis.annotate(quant_matrix="coverage_qnorm")
    chipseq_analysis.annotate_with_sample_metadata(attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    chipseq_analysis.to_pickle()

    # Unsupervised analysis
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=None,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")
    # without histone marks
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=[s for s in chipseq_analysis.samples if "H3" not in s.ip],
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks.no_histones",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")
    # only complex members
    sel_samples = [s for s in chipseq_analysis.samples if ("ARID" in s.ip) | ("SMAR" in s.ip) | ("PBRM" in s.ip)]
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=sel_samples,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks.only_baf",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")


    # distinguish between BAF and pBAF specific
    pBAF_vs_BAF(chispeq_analysis)

    # # Investigate global changes in accessibility
    # global_changes(atacseq_samples)




    # See ATAC and RNA together
    rnaseq_analysis.accessibility_expression()


    # Deeper ATAC-seq data
    nucleosome_changes(atacseq_analysis, atacseq_samples)
    investigate_nucleosome_positions(atacseq_samples)
    phasograms(atacseq_samples)


    interaction_new_data()


def main_analysis_pipeline(analysis, data_type, cell_type):
    """
    Main analysis pipeline for ATAC-seq and RNA-seq.
    Gets quantification matrices, normalizes them,
    performes unsupervised and supervised analysis and
    gets and plots enrichments for supervised analysis.
    """
    # set thresholds
    if data_type == "ATAC-seq" and cell_type == "HAP1":
        alpha = 0.01
        abs_fold_change = 1
    else:
        alpha = 0.05
        abs_fold_change = 0


    if data_type == "ATAC-seq":

        # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT
        # Get consensus peak set from all samples
        analysis.get_consensus_sites(region_type="summits")
        analysis.calculate_peak_support(region_type="summits")

        # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
        # Get coverage values for each peak in each sample of ATAC-seq and ChIPmentation
        analysis.measure_coverage()
        # normalize coverage values
        analysis.normalize_coverage_quantiles()
        analysis.get_peak_gccontent_length()
        analysis.normalize_gc_content()

        # Annotate peaks with closest gene
        analysis.get_peak_gene_annotation()
        # Annotate peaks with genomic regions
        analysis.get_peak_genomic_location()
        # Annotate peaks with ChromHMM state from CD19+ cells
        analysis.get_peak_chromatin_state(chrom_state_file="data/external/HAP1_12_segments.annotated.bed")
        # Annotate peaks with closest gene, chromatin state,
        # genomic location, mean and variance measurements across samples
        analysis.annotate()
        analysis.annotate_with_sample_metadata(attributes=sample_attributes)
        analysis.to_pickle()

        # QC plots
        # plot general peak set features
        analysis.plot_peak_characteristics()
        # plot coverage features across peaks/samples
        analysis.plot_coverage()
        analysis.plot_variance()

        quant_matrix = "accessibility"
        feature_name = "sites"


    if data_type == "RNA-seq":
        # Get gene expression
        analysis.get_gene_expression(
            samples=analysis.samples,
            sample_attributes=["sample_name", "knockout", 'treatment', "replicate", "clone", "batch"])

        # see expression of knocked-out genes + other complex members
        baf_genes = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv"), squeeze=True)
        knockout_plot(
            analysis=analysis,
            knockout_genes=baf_genes,
            output_prefix="complex_expression")

        quant_matrix = "expression_annotated"
        feature_name = "genes"

    # Unsupervised analysis
    red_samples = [s for s in analysis.samples if (
        # ("sh" not in s.name) and
        ("dBet" not in s.name) and
        ("BRD4" not in s.name) and
        ("parental" not in s.name))]
    unsupervised_analysis(
        analysis,
        quant_matrix=quant_matrix,
        samples=red_samples,
        attributes_to_plot=plotting_attributes,
        plot_prefix="all_{}".format(feature_name),
        plot_max_attr=20,
        plot_max_pcs=6,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False,
        output_dir="{results_dir}/unsupervised.20180131")

    # fix batch effect
    analysis.matrix_batch_fix = fix_batch_effect(
        getattr(analysis, quant_matrix), analysis.samples,
        batch_variable="batch", standardize=True, intermediate_save=True)
    file = os.path.join(analysis.results_dir, analysis.name + ".{}.annotated_metadata.batch_fix.csv".format(quant_matrix))
    analysis.matrix_batch_fix.to_csv(file)

    unsupervised_analysis(
        analysis, quant_matrix="matrix_batch_fix", samples=None,
        attributes_to_plot=plotting_attributes, plot_prefix="all_{}".format(feature_name),
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False, display_corr_values=False,
        output_dir="{results_dir}/unsupervised.batch_fix")
    quant_matrix = "matrix_batch_fix"

    if data_type == "RNA-seq":
        knockout_plot(
            analysis=analysis,
            knockout_genes=baf_genes,
            expression_matrix=analysis.matrix_batch_fix,
            output_prefix="complex_expression.batch_fix")


    # only with certain subunits
    if cell_type == "HAP1":
        to_exclude = ['SMARCA4', "SMARCC1", "ARID1A", "ARID1B", "BCL7B"]
        # no strong KOs
        unsupervised_analysis(
            analysis, quant_matrix=quant_matrix,
            samples=[s for s in analysis.samples if s.knockout not in to_exclude],
            attributes_to_plot=plotting_attributes, plot_prefix="all_{}_no_strong_knockouts".format(feature_name),
            plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
            output_dir="{results_dir}/unsupervised.batch_fix")

        unsupervised_analysis(
            analysis, quant_matrix=quant_matrix,
            samples=[s for s in analysis.samples if (("BRD4" not in s.name) and ("sh" not in s.name) and ("dBet" not in s.name) and ("parental" not in s.name))],
            attributes_to_plot=plotting_attributes, plot_prefix="all_{}_no_kd_dBet_BRD4".format(feature_name),
            plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
            output_dir="{results_dir}/unsupervised.batch_fix")

        unsupervised_analysis(
            analysis, quant_matrix=quant_matrix,
            samples=[s for s in analysis.samples if (("sh" in s.name) or ("dBet" in s.name) or ("parental" in s.name))],
            attributes_to_plot=plotting_attributes, plot_prefix="all_{}_only_kd_dBet".format(feature_name),
            plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
            output_dir="{results_dir}/unsupervised.batch_fix")

    # Supervised analysis
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['data_type'] == data_type) &
        (comparison_table['cell_type'] == cell_type) &
        (comparison_table['comparison_type'] == 'differential')]
    analysis.differential_results = differential_analysis(
        analysis,
        comparison_table,
        data_type=data_type,
        samples=[s for s in analysis.samples if s.name in comparison_table['sample_name'].tolist()],
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        covariates=None,
        alpha=0.05,
        overwrite=True)
    analysis.differential_results = analysis.differential_results.set_index("index")
    analysis.to_pickle()

    if cell_type == "HAP1":
        differential_overlap(
            analysis.differential_results[
                (analysis.differential_results['padj'] < alpha) &
                (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
            getattr(analysis, quant_matrix).shape[0],
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            data_type=data_type)

    plot_differential(
        analysis,
        analysis.differential_results, # analysis.differential_results[~analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")], 
        matrix=getattr(analysis, quant_matrix),
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        output_prefix="differential_analysis.batch_fix.top1000",
        data_type=data_type,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
        group_variables=plotting_attributes)

    # repeat without SMARCA4, ARID1A, SMARCC1
    if cell_type == "HAP1":
        plot_differential(
            analysis,
            analysis.differential_results[
                ~analysis.differential_results['comparison_name'].isin(
                    ['SMARCA4', 'ARID1A', 'SMARCC1'])],
            comparison_table=comparison_table,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_very_strong",
            alpha=alpha,
            corrected_p_value=True,
            fold_change=abs_fold_change,
            rasterized=True, robust=True)
        # repeat without the strongest guys
        plot_differential(
            analysis,
            analysis.differential_results[
                ~analysis.differential_results['comparison_name'].isin(
                    ['SMARCA4', 'ARID1A', 'SMARCC1', "BCL7B", "ARID1B"])],
            matrix=getattr(analysis, quant_matrix),
            comparison_table=comparison_table,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_strong",
            alpha=alpha,
            corrected_p_value=True,
            fold_change=abs_fold_change,
            rasterized=True, robust=True)

    # repeat for specific subsets
    if cell_type == "HAP1":
        plot_differential(
            analysis,
            analysis.differential_results[
                ~analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")],
            matrix=getattr(analysis, quant_matrix),
            samples=[s for s in analysis.samples if not (("sh" in s.name) or ("dBet" in s.name) or ("BRD4" in s.name) or ("parental" in s.name))],
            comparison_table=comparison_table[~comparison_table['comparison_name'].str.contains("sh|dBet|BRD4")],
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_knockdown_dBet_BRD4",
            alpha=alpha,
            corrected_p_value=True,
            fold_change=abs_fold_change,
            rasterized=True, robust=True)

        differential_overlap(
            analysis.differential_results[
                (~analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")) &
                (analysis.differential_results['padj'] < alpha) &
                (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
            getattr(analysis, quant_matrix).shape[0],
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_knockdown_dBet_BRD4",
            data_type=data_type)

        plot_differential(
            analysis,
            analysis.differential_results[
                analysis.differential_results['comparison_name'].str.contains("sh|dBet")],
            matrix=getattr(analysis, quant_matrix),
            samples=[s for s in analysis.samples if (("sh" in s.name) or ("dBet" in s.name))],
            comparison_table=comparison_table[comparison_table['comparison_name'].str.contains("sh|dBet")],
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.only_kd_dBet",
            alpha=alpha,
            corrected_p_value=True,
            fold_change=abs_fold_change,
            rasterized=True, robust=True)

    differential_enrichment(
        analysis,
        analysis.differential_results[
            (analysis.differential_results['padj'] < alpha) &
            (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        genome="hg19",
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True)

    collect_differential_enrichment(
        analysis.differential_results[
            (analysis.differential_results['padj'] < alpha) &
            (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        permissive=False)

    if data_type == "RNA-seq":
        enrichment_table = pd.read_csv(
            os.path.join("{}/differential_analysis_{}".format(analysis.results_dir, data_type), "differential_analysis.enrichr.csv"))
        enrichment_table = enrichment_table[~enrichment_table['comparison_name'].str.contains("sh|dBet|BRD4")]

        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis",
            direction_dependent=True,
            top_n=5)
        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_knockdown_dBet_BRD4",
            direction_dependent=True,
            top_n=5)
        plot_differential_enrichment(
            enrichment_table[enrichment_table['gene_set_library'].isin(
                ['WikiPathways_2016'])],
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.top3",
            direction_dependent=True,
            top_n=10)
        plot_differential_enrichment(
            enrichment_table[
                (~enrichment_table['comparison_name'].str.contains("sh|dBet|BRD4")) &
                (enrichment_table['gene_set_library'].str.contains("GO_Biological_Process"))],
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_knockdown_dBet_BRD4.only_top",
            direction_dependent=True,
            top_n=1)
    elif data_type == "ATAC-seq":
        for enrichment_name, enrichment_type in [('motif', 'meme_ame'), ('lola', 'lola'), ('enrichr', 'enrichr')]:
            try:
                enrichment_table = pd.read_csv(
                    os.path.join("{}/differential_analysis_{}".format(analysis.results_dir, data_type), "differential_analysis" + ".{}.csv".format(enrichment_type)))
            except pd.errors.EmptyDataError:
                print("Enrichment dataframe of {} is empty.".format(enrichment_type))
                continue

            # plot_differential_enrichment(
            #     enrichment_table,
            #     enrichment_name,
            #     data_type=data_type,
            #     output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            #     direction_dependent=True,
            #     top_n=5 if enrichment_name != "motif" else 300)


            plot_differential_enrichment(
                enrichment_table[~enrichment_table['comparison_name'].str.contains("sh|dBet|BRD4")],
                enrichment_name,
                data_type=data_type,
                output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
                output_prefix="differential_analysis.no_knockdown_dBet_BRD4",
                direction_dependent=True,
                top_n=5 if enrichment_name != "motif" else 300)

    return analysis


def fix_batch_effect(
        matrix, samples,
        batch_variable="pool_id", standardize=True, intermediate_save=True):
    """
    """
    import pandas as pd
    from statsmodels.formula.api import ols
    from statsmodels.graphics.factorplots import interaction_plot
    import matplotlib.pyplot as plt
    from scipy import stats
    from tqdm import tqdm

    m = matrix.loc[:, [s.name for s in samples]]
    m.index = ['feature'] * len(m.index)

    # standardize variables
    if standardize:
        m = m.apply(stats.zscore)

    # get coefficients for batches
    levels = m.columns.get_level_values(batch_variable).unique()
    n = len(levels)
    res = pd.DataFrame(np.zeros((m.shape[0], n)))
    for i in tqdm(range(m.shape[0]), total=m.shape[0]):
        data = m.T.iloc[:, i].reset_index()

        formula = 'feature ~ {} - 1'.format(batch_variable)
        model = ols(formula, data).fit()
        res.iloc[i] = model.params.values

        if i % 10000 == 0:
            if intermediate_save:
                res.iloc[:i, :].to_csv(os.path.join(
                    "batch.{}.fit.csv".format(batch_variable)), index=False)
    res.columns = model.params.index
    res.index = matrix.index
    res.to_csv(os.path.join(
        "batch.{}.fit.csv".format(batch_variable)), index=True)

    # fix
    m.index = matrix.index
    for batch in levels:
        s = m.columns[m.columns.get_level_values(batch_variable) == batch]
        for s_ in s:
            print(batch, s_)
            m.loc[:, s_] = m.loc[:, s_] - res.loc[:, res.columns.str.contains(batch)].squeeze()

    return m


def diff_atac_heatmap_with_chip(
        atac_analysis, chipseq_analysis,
        output_dir="results"):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    output_prefix = "atac-chip_comparison"
    quantity = "accesibility"
    var_name = "regions"
    robust = True
    rasterized = True

    results = atac_analysis.differential_results
    results = results[~results['comparison_name'].str.contains("sh|dBet|parental|BRD4")]
    results['diff'] = (results["padj"] < 0.01) & (results["log2FoldChange"].abs() > 1)

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[(
        comparison_table['data_type'] == "ATAC-seq") &
        (comparison_table['cell_type'] == "HAP1") &
        (comparison_table['comparison_type'] == 'differential')]

    atac_matrix = atac_analysis.accessibility_batch_fix.copy()
    atac_matrix = atac_matrix.loc[:, ~atac_matrix.columns.get_level_values("sample_name").str.contains("dBet|sh|parental|BRD4")]

    # PLOTS
    comparisons = sorted(results["comparison_name"].drop_duplicates())
    n_side = int(np.ceil(np.sqrt(len(comparisons))))

    # Observe values of variables across all comparisons
    all_diff = results[results["diff"] == True].index.drop_duplicates()
    sample_cols = atac_matrix.columns.get_level_values("sample_name").tolist()

    if comparison_table is not None:
        if results["comparison_name"].drop_duplicates().shape[0] > 1:
            groups = pd.DataFrame()
            for sample_group in results["comparison_name"].drop_duplicates():
                c = comparison_table.loc[comparison_table["sample_group"] == sample_group, "sample_name"].drop_duplicates()
                if c.shape[0] > 0:
                    groups.loc[:, sample_group] = atac_matrix[[d for d in c if d in sample_cols]].mean(axis=1)

            if groups.empty:
                # It seems comparisons were not done in a all-versus-all fashion
                for group in comparison_table["sample_group"].drop_duplicates():
                    c = comparison_table.loc[comparison_table["sample_group"] == group, "sample_name"].drop_duplicates()
                    if c.shape[0] > 0:
                        groups.loc[:, group] = atac_matrix[c].mean(axis=1)

            # Select only differential regions from groups
            groups = groups.loc[all_diff, :]

            figsize = (max(5, 0.12 * groups.shape[1]), 5)
            # Heatmaps
            # Comparison level
            g = sns.clustermap(
                groups,
                xticklabels=True, yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
                metric="correlation", robust=robust, rasterized=rasterized, figsize=figsize)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(
                groups,
                xticklabels=True, yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
                cmap="RdBu_r", metric="correlation", robust=robust, rasterized=rasterized, figsize=figsize)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Fold-changes and P-values
    # pivot table of genes vs comparisons
    fold_changes = pd.pivot_table(results.loc[all_diff, :].reset_index(), index=results.index.name, columns="comparison_name", values="log2FoldChange")
    p_values = -np.log10(pd.pivot_table(results.loc[all_diff, :].reset_index(), index=results.index.name, columns="comparison_name", values="padj"))

    # fold
    if fold_changes.shape[1] > 1:
        figsize = (max(5, 0.12 * fold_changes.shape[1]), 5)

        g = sns.clustermap(fold_changes.corr(),
            xticklabels=False, yticklabels=True, cbar_kws={"label": "Pearson correlation\non fold-changes"},
            cmap="BuGn", vmin=0, vmax=1, metric="correlation", rasterized=rasterized, figsize=(figsize[0], figsize[0]))
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=300, metric="correlation")

        g = sns.clustermap(fold_changes.loc[all_diff, :],
            xticklabels=True, yticklabels=False, cbar_kws={"label": "Fold-change of\ndifferential {}".format(var_name)},
            cmap="RdBu_r", robust=True, metric="correlation", rasterized=rasterized, figsize=figsize)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Sample level
    if type(atac_matrix.columns) is pd.core.indexes.multi.MultiIndex:
        atac_matrix.columns = atac_matrix.columns.get_level_values("sample_name")

    atac_matrix = atac_matrix.loc[all_diff, :]
    figsize = (max(5, 0.12 * atac_matrix.shape[1]), 5)

    col = sns.clustermap(
        atac_matrix,
        xticklabels=False, yticklabels=False, metric="correlation", figsize=(1, 1), rasterized=True, robust=robust)

    g = sns.clustermap(
        atac_matrix,
        col_linkage=col.dendrogram_col.linkage,
        xticklabels=True, yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
        cmap="RdBu_r", metric="euclidean", figsize=figsize, rasterized=rasterized, robust=robust)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.new.euclidean.svg".format(var_name)), bbox_inches="tight", dpi=100)

    g2 = sns.clustermap(
        atac_matrix,
        xticklabels=True, yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
        cmap="RdBu_r", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust)
    g2.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)


    # add ChIP-seq data with same clustering
    # but z-scored internally
    regs = atac_matrix.iloc[g.dendrogram_row.reordered_ind, :].index

    background_mean = chipseq_analysis.accessibility.loc[:, chipseq_analysis.accessibility.columns.get_level_values("ip").str.contains("IgG|Input")].mean(axis=1)
    chip_over_background = chipseq_analysis.accessibility.copy()
    for s in chipseq_analysis.accessibility.columns:
        chip_over_background.loc[:, s] = chipseq_analysis.accessibility.loc[:, s] - background_mean

    chip_over_background = chip_over_background.drop(
        ['ChIP-seq_HAP1_WT_H3K27ac_r1'], axis=1)

    c = chip_over_background.loc[regs, ~chip_over_background.columns.get_level_values("ip").str.contains("IgG|Input")]

    figsize = (0.12 * c.shape[1], 5)

    from scipy.stats import zscore
    # chip_z = pd.DataFrame(zscore(c, axis=0), index=regs, columns=chipseq_analysis.accessibility.columns)
    for label, i in [
            ("all", "H3|CTCF|Pol|ARID|BRD|PBRM|SMARC"),
            ("histones", "H3|CTCF|Pol"),
            ("baf", "ARID|PBRM|SMARC")
        ]:
        c2 = c.loc[:, c.columns.get_level_values(0).str.contains(i)]

        fig, axis = plt.subplots(1, figsize=figsize)
        sns.heatmap(c2, xticklabels=True, yticklabels=False, ax=axis, vmin=-3, vmax=3, rasterized=rasterized)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.chip_over_input.only_{}.heatmap.svg".format(var_name, label)), bbox_inches="tight", dpi=300)

        g3 = sns.clustermap(c2, xticklabels=True, yticklabels=False, row_cluster=False, vmin=-3, vmax=3, rasterized=rasterized, figsize=figsize)
        g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g3.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.chip_over_input.only_{}.clustermap.svg".format(var_name, label)), bbox_inches="tight", dpi=300)

        c_mean = c2.T.groupby('ip').mean().T
        figsize = (0.12 * c_mean.shape[1], 5)

        fig, axis = plt.subplots(1, figsize=figsize)
        sns.heatmap(c_mean, xticklabels=True, yticklabels=False, ax=axis, vmin=-3, vmax=3, rasterized=rasterized)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.chip_over_input.only_{}.mean.heatmap.svg".format(var_name, label)), bbox_inches="tight", dpi=300)

        g4 = sns.clustermap(c_mean, xticklabels=True, yticklabels=False, row_cluster=False, vmin=-3, vmax=3, rasterized=rasterized, figsize=figsize)
        g4.ax_heatmap.set_xticklabels(g4.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g4.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.chip_over_input.only_{}.mean.clustermap.svg".format(var_name, label)), bbox_inches="tight", dpi=300)

        g4 = sns.clustermap(c_mean, xticklabels=True, yticklabels=False, row_cluster=False, vmin=-3, vmax=3, rasterized=rasterized, figsize=figsize, cmap="PuOr_r")
        g4.ax_heatmap.set_xticklabels(g4.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g4.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.chip_over_input.only_{}.mean.clustermap.PuOr_r.svg".format(var_name, label)), bbox_inches="tight", dpi=300)


def pBAF_vs_BAF(chispeq_analysis, output_dir="{results_dir}/pBAF_vs_BAF"):
    from scipy.stats import mannwhitneyu
    from scipy.stats import gaussian_kde
    from statsmodels.nonparametric.smoothers_lowess import lowess

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=chipseq_analysis.results_dir)

    chip_mean = chipseq_analysis.accessibility.T.groupby("ip").mean().T
    chip_std = chipseq_analysis.accessibility.T.groupby("ip").std().dropna().T
    background_mean = chip_mean[["IgG", "Input"]].mean(1)
    stats = np.log2(chip_mean["ARID2"] /
                    chip_mean["ARID1A"]).to_frame(name="fold_change")
    stats = stats.join(((chip_mean["ARID2"] / background_mean) - (
        chip_mean["ARID1A"] / background_mean)).to_frame(name="fold_change_over_background"))
    stats = stats.join(chip_mean[["ARID2", "ARID1A"]].mean(
        1).to_frame(name="comparison_mean"))
    stats = stats.join(chip_mean.mean(1).to_frame(name="global_mean"))
    stats = stats.join(chip_std['ARID2'].to_frame(name="ARID2_std"))
    stats = stats.join(np.log2(
        (stats['comparison_mean'] / background_mean)).to_frame(name="comparison_mean_over_background"))
    stats = stats.join(np.log2(
        (stats['global_mean'] / background_mean)).to_frame(name="global_mean_over_background"))

    # standardize fold change
    bounds = np.linspace(0, stats['comparison_mean'].max(), 250)
    for i, (start, end) in enumerate(zip(bounds[:-2], bounds[1: -1])):
        for metric in ['fold_change', 'fold_change_over_background']:
            r = stats.loc[(stats['comparison_mean'] > start) &
                          (stats['comparison_mean'] < end)].index
            v = stats.loc[r, metric]
            stats.loc[r, 'norm_{}'.format(metric)] = (v - v.mean()) / v.std()

    # let's try a bivariate gaussian kernel
    # separately for positive and negative to avoid biases in center of mass
    kernel = gaussian_kde(stats.loc[stats['norm_fold_change'] > 0, [
                          "comparison_mean", "norm_fold_change"]].T.values)
    stats.loc[stats['norm_fold_change'] > 0, "density"] = kernel(
        stats.loc[stats['norm_fold_change'] > 0, ["comparison_mean", "norm_fold_change"]].T.values)
    kernel = gaussian_kde(stats.loc[stats['norm_fold_change'] <= 0, [
                          "comparison_mean", "norm_fold_change"]].T.values)
    stats.loc[stats['norm_fold_change'] <= 0, "density"] = kernel(
        stats.loc[stats['norm_fold_change'] <= 0, ["comparison_mean", "norm_fold_change"]].T.values)

    diff = stats.loc[
        (stats['comparison_mean_over_background'] > 0) &
        (stats['ARID2_std'] < 2) &
        (
            ((stats['norm_fold_change'] > 1.5)) |
            ((stats['norm_fold_change'] < -1.5))) &
        (stats['density'] < np.percentile(stats['density'].dropna(), 5))
    ]
    diff['direction'] = (diff['norm_fold_change'] >= 0).astype(int).replace(0, -1)
    # for ip in chip_mean.columns:
    #     stats = stats.join(np.log2(chip_mean[ip] / background_mean).to_frame(name=ip))
    diff.sort_values(['direction', 'norm_fold_change'], ascending=False).to_csv(
        os.path.join(output_dir, chipseq_analysis.name + ".baf_peaks.pBAF_vs_BAF.csv"), index=True)

    chipseq_analysis.pbaf_sites = diff
    chipseq_analysis.to_pickle()

    cmap = plt.get_cmap("RdBu_r")
    cmap2 = plt.get_cmap("inferno_r")
    norm = matplotlib.colors.Normalize(vmin=-3, vmax=3)
    norm2 = matplotlib.colors.Normalize(
        vmin=0, vmax=stats.loc[:, 'density'].max())
    fig, axis = plt.subplots(1, 4, figsize=(4 * 4, 4))
    axis[0].scatter(stats.loc[:, 'comparison_mean'], stats.loc[:, 'fold_change'], alpha=0.2,
                    s=3, color=cmap(norm(stats.loc[:, 'fold_change'].values)), rasterized=True)
    axis[1].scatter(stats.loc[:, 'comparison_mean'], stats.loc[:, 'norm_fold_change'], alpha=0.2,
                    s=3, color=cmap(norm(stats.loc[:, 'norm_fold_change'].values)), rasterized=True)
    axis[2].scatter(stats.loc[:, 'comparison_mean'], stats.loc[:, 'norm_fold_change'],
                    alpha=0.2, s=3, color=cmap2(norm2(stats.loc[:, 'density'])), rasterized=True)
    axis[3].scatter(stats.loc[diff.index, 'comparison_mean'], stats.loc[diff.index,
                                                                        'norm_fold_change'], alpha=0.2, s=3, color='red', rasterized=True)
    for ax in axis[1:]:
        ax.set_xlim(0, stats.loc[:, 'comparison_mean'].max())
        ax.set_ylim(-6, 6)
    for ax in axis:
        ax.axhline(0, color="black", linestyle="--")
        ax.set_xlabel("Mean (ARID2, ARID1A)")
    axis[0].set_ylabel("log2(fold-change)")
    axis[1].set_ylabel("normalized log2(fold-change)")
    axis[2].set_ylabel("normalized log2(fold-change)")
    axis[3].set_ylabel("normalized log2(fold-change)")
    axis[0].set_title("MA plot")
    axis[1].set_title("MA plot with standardized fold-change in percentiles")
    axis[2].set_title(
        "MA plot with bivariate gaussian kernel density estimator")
    axis[3].set_title("Final set of differential pBAF/BAF peaks")
    axis[1].axhline(-1, color="red", linestyle="--")
    axis[1].axhline(1, color="red", linestyle="--")
    axis[1].axvline(3, color="red", linestyle="--")
    fig.savefig(os.path.join(output_dir, chipseq_analysis.name +
                             ".baf_peaks.pBAF_vs_BAF.finding_differential.svg"), bbox_inches="tight", dpi=300)

    # get values over input
    chip_over_background = chipseq_analysis.accessibility.copy()
    for s in chipseq_analysis.accessibility.columns:
        chip_over_background.loc[:, s] = chipseq_analysis.accessibility.loc[:, s] - background_mean

    chip_over_background = chip_over_background.drop(
        ['ChIP-seq_HAP1_WT_H3K27ac_r1'], axis=1)

    m = chip_over_background.loc[diff.index, ~chip_over_background.columns.get_level_values(
        "ip").str.contains("IgG|Input")]
    g = sns.clustermap(m, yticklabels=False, rasterized=True, robust=True,
                       col_cluster=True, row_cluster=True, metric="euclidean")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_dir, chipseq_analysis.name +
                           ".baf_peaks.pBAF_vs_BAF.differential.heatmap.over_background.svg"), bbox_inches="tight", dpi=300)

    g2 = sns.clustermap(
        m.T.groupby("ip").mean().T,
        row_linkage=g.dendrogram_row.linkage,
        yticklabels=False, rasterized=True, robust=True, col_cluster=True)
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.savefig(os.path.join(output_dir, chipseq_analysis.name +
                            ".baf_peaks.pBAF_vs_BAF.differential.heatmap.over_background.groups.svg"), bbox_inches="tight", dpi=300)


    # Get histone mark enrichment in pBAF vs BAF complex
    ips = chip_over_background.columns.get_level_values("ip").unique()
    s = chip_over_background.loc[diff.index].T.groupby("ip").mean().T.join(diff)

    ss = pd.melt(s, id_vars=[c for c in s.columns if c not in ips], var_name="ip")
    fig, axis = plt.subplots(figsize=(4 * 2, 4))
    sns.barplot(data=ss, x="direction", y="value", hue="ip", ax=axis)
    fig.savefig(os.path.join(output_dir, chipseq_analysis.name +
                            ".baf_peaks.pBAF_vs_BAF.differential.histone_enrichment.barplot.svg"), bbox_inches="tight", dpi=300)

    # Get ATAC-seq enrichment in pBAF vs BAF complex
    atac_analysis.gene_annotation = atac_analysis.gene_annotation.reset_index()
    atac_analysis.gene_annotation.index = (
        atac_analysis.gene_annotation['chrom'] + ":" +
        atac_analysis.gene_annotation['start'].astype(str) + "-" +
        atac_analysis.gene_annotation['end'].astype(str))

    enrichment = pd.DataFrame()
    for direction, label in [(1, "pBAF"), (-1, "BAF")]:
        # get ATAC-seq peaks that overlap with BAF/pBAF peaks
        i = diff.loc[diff['direction'] == direction].index
        sites_bed = os.path.join(output_dir, chipseq_analysis.name + ".baf_peaks.pBAF_vs_BAF.{}_specific.bed".format(label))
        chipseq_analysis.coverage.loc[i, ['chrom', 'start', 'end']].to_csv(sites_bed, sep="\t", header=False, index=False)

        # Get total accesibility in those sites per knockout
        intersect = atac_analysis.sites.intersect(sites_bed, wa=True).to_dataframe()
        j = intersect.iloc[:, 0] + ':' + intersect.iloc[:, 1].astype(str) + "-" + intersect.iloc[:, 2].astype(str)

        m = pd.melt(atac_analysis.accessibility.loc[j].T.groupby("knockout").mean().T.reset_index(), id_vars=["index"])
        m['direction'] = label
        m['value_type'] = "absolute"
        m['data_type'] = "accesibility"
        enrichment = enrichment.append(m, ignore_index=True)

        # Get fold-change compared to control
        m = (
            atac_analysis.differential_results[~atac_analysis.differential_results['comparison_name'].str.contains("sh|parental|dBet")]
            .loc[j, ["log2FoldChange", "comparison_name"]]
            .rename(columns={"comparison_name": "knockout", "log2FoldChange": "value"}))
        m['direction'] = label
        m['value_type'] = "fold_change"
        m['data_type'] = "accesibility"
        enrichment = enrichment.append(m, ignore_index=True)

        # Get gene expression of genes associated with these regions
        g = atac_analysis.gene_annotation.loc[j]
        g = g.loc[g['distance'] <= 10000, "gene_name"]
        m = pd.melt(rnaseq_analysis.expression_annotated.loc[g].T.groupby("knockout").mean().T.reset_index(), id_vars=["gene_name"]).rename(columns={"gene_name": "index"})
        m['direction'] = label
        m['value_type'] = "absolute"
        m['data_type'] = "expression"
        enrichment = enrichment.append(m, ignore_index=True)

        # Get fold-change compared to control
        m = (
            rnaseq_analysis.differential_results[~rnaseq_analysis.differential_results['comparison_name'].str.contains("sh|parental|dBet")]
            .loc[g, ["log2FoldChange", "comparison_name"]]
            .rename(columns={"comparison_name": "knockout", "log2FoldChange": "value"}))
        m['direction'] = label
        m['value_type'] = "fold_change"
        m['data_type'] = "expression"
        enrichment = enrichment.append(m, ignore_index=True)

    enrichment.to_csv(os.path.join(output_dir, chipseq_analysis.name +
                            ".baf_peaks.pBAF_vs_BAF.differential.enrichments.csv"), index=False)

    fig, axis = plt.subplots(2, 2, figsize=(4 * 4, 2 * 4))
    axis[0, 0].set_title("Accessibility")
    axis[0, 1].set_title("Expression")
    sns.barplot(data=enrichment[
        (enrichment['data_type'] == "accesibility") &
        (enrichment['value_type'] == "absolute")
    ], x="knockout", y="value", hue="direction", ax=axis[0, 0])
    sns.barplot(data=enrichment[
        (enrichment['data_type'] == "accesibility") &
        (enrichment['value_type'] == "fold_change")
    ], x="knockout", y="value", hue="direction", ax=axis[1, 0])
    sns.barplot(data=enrichment[
        (enrichment['data_type'] == "expression") &
        (enrichment['value_type'] == "absolute")
    ], x="knockout", y="value", hue="direction", ax=axis[0, 1])
    sns.barplot(data=enrichment[
        (enrichment['data_type'] == "expression") &
        (enrichment['value_type'] == "fold_change")
    ], x="knockout", y="value", hue="direction", ax=axis[1, 1])
    axis[0, 0].set_ylabel("Accessiblity")
    axis[1, 0].set_ylabel("Log2 fold-change")
    axis[0, 1].set_ylabel("Expression")
    axis[1, 1].set_ylabel("Log2 fold-change")
    for ax in axis.flatten():
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, chipseq_analysis.name +
                            ".baf_peaks.pBAF_vs_BAF.differential.all_enrichments.<10kb.barplot.svg"), bbox_inches="tight", dpi=300)


def synthetic_letality(
        output_dir="{results_dir}/synthetic_lethality"):

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format("results")
    analysis_name = "baf_complex"

    sel_vars = [
        "CellLine", "GeneSymbol", "Replicate", "Hit",
        "TreatmentMean_all", "ControlMean_all", "TreatmentSD_all", "ControlSD_all",
        "pValue_all", "pLogP_all", "DeltaPOC"]
    num_vars = sel_vars[4:]

    df = pd.read_csv(os.path.join("metadata", "original", "20171120_inclOutliers_inclHits_3.txt"), sep="\t", decimal=",")
    df = df[sel_vars]
    df['cell_line'] = "HAP1"

    # Curate/sanitize
    df.loc[df['CellLine'].str.contains("A549"), "cell_line"] = "A549"
    df['knockout'] = df['CellLine']
    df['knockout'] = df['knockout'].str.replace("A549 ", "")
    df.loc[(df['knockout'] == "WT") & (df['cell_line'] == "A549"), 'knockout'] = "SMARCA4"
    df.loc[df['knockout'] == "SMARCA4 REC", 'knockout'] = "WT"
    df.loc[df['knockout'].str.startswith("ARID1A"), 'knockout'] = "ARID1A"
    df.loc[df['knockout'].str.startswith("SMARCA4"), 'knockout'] = "SMARCA4"
    df["clone"] = df['CellLine'].str.split(" ").apply(lambda x: x[-1])
    df['knockdown'] = df['GeneSymbol']
    for var_ in num_vars:
        df[var_] = df[var_].astype(float)

    # only single knockdowns
    df2 = df.loc[~df["knockdown"].astype(str).str.contains("_| "), :]

    for d, label in [(df, "multiple_KDs"), (df2, "single_KDs")]:

        # heatmaps
        fig, axis = plt.subplots(2, 2, figsize=(8, 8))
        for i, cell_line in enumerate(d['cell_line'].unique()):
            piv = pd.pivot_table(data=d.loc[d['cell_line'] == cell_line],
                index="knockdown", columns="knockout", values="pLogP_all")

            sns.heatmap(
                piv, cbar_kws={"label": "-log10(synthetic interaction)"},
                ax=axis[i, 0], xticklabels=True, yticklabels=True)
            sns.heatmap(
                pd.DataFrame(zscore(piv, 0), index=piv.index, columns=piv.columns), xticklabels=True, yticklabels=True,
                cbar_kws={"label": "Synthetic interaction\n(column Z-score)"}, ax=axis[i, 1])

            axis[i, 0].set_title(cell_line)
            axis[i, 1].set_title(cell_line)

        for ax in axis.flatten():
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize="xx-small")
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
        fig.savefig(os.path.join(output_dir, analysis_name +
                                ".{}.p_value.heatmaps.svg".format(label)), bbox_inches="tight", dpi=300)

        cell_line = "HAP1"
        piv = pd.pivot_table(data=d.loc[d['cell_line'] == cell_line],
            index="knockdown", columns="knockout", values="pLogP_all")
        g = sns.clustermap(
            piv.corr(), metric="euclidean",
            cbar_kws={"label": "Distance in synthetic interactions\n(Euclidean distance)"},
            cmap="RdBu_r", vmin=-1, vmax=1,
            xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, analysis_name + ".{}.correlation.euclidean.svg".format(label)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(
            piv.corr(), metric="correlation",
            cbar_kws={"label": "Distance in synthetic interactions\n(Pearson correlation)"},
            cmap="RdBu_r", vmin=-1, vmax=1,
            xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, analysis_name + ".{}.correlation.pearson.svg".format(label)), bbox_inches="tight", dpi=300)


    # Compare with similarity in ATAC-seq and RNA-seq
    screen = pd.pivot_table(data=df2.loc[df2['cell_line'] == "HAP1"],
                index="knockdown", columns="knockout", values="pLogP_all")
    screen_corr = screen.corr()
    screen_corr.index.name = "k1"
    screen_corr.columns.name = "k2"
    screen_m = pd.melt(screen_corr.reset_index(), id_vars='k1', value_name='screen')
    screen_m['k1'] = screen_m['k1'].str.replace("ACTIN", "ACTB")
    screen_m['k2'] = screen_m['k2'].str.replace("ACTIN", "ACTB")
    atac_corr = atac_analysis.accessibility.T.groupby("knockout").mean().T.corr()
    atac_corr.index.name = "k1"
    atac_corr.columns.name = "k2"
    atac_corr_m = pd.melt(atac_corr.reset_index(), id_vars='k1', value_name='atac')
    rna_corr = rnaseq_analysis.expression_annotated.T.groupby("knockout").mean().T.corr()
    rna_corr.index.name = "k1"
    rna_corr.columns.name = "k2"
    rna_corr_m = pd.melt(rna_corr.reset_index(), id_vars='k1', value_name='rna')

    corrs = (screen_m.set_index(['k1', 'k2'])
        .join(atac_corr_m.set_index(['k1', 'k2']))
        .join(rna_corr_m.set_index(['k1', 'k2'])))

    baf_members = ["ARID1A", "ARID1B"]
    pbaf_members = ["ARID2", "PBRM1", "BRD9", "PHF10", "DPF1"]

    corrs['baf'] = (corrs.index.get_level_values("k1").isin(baf_members).astype(int) + corrs.index.get_level_values("k2").isin(baf_members).astype(int))
    corrs['pbaf'] = -(corrs.index.get_level_values("k1").isin(pbaf_members).astype(int) + corrs.index.get_level_values("k2").isin(pbaf_members).astype(int))
    corrs['color'] = corrs['baf'] + corrs['pbaf']
    corrs.loc[(corrs['baf'] == 1) & (corrs['pbaf'] == -1), 'color'] = 3

    fig, axis = plt.subplots(1, 2, figsize=(8, 4))
    cmap = plt.get_cmap("coolwarm")
    cmap2 = plt.get_cmap("PuOr")
    norm = matplotlib.colors.Normalize(vmin=-2, vmax=2)
    for i, data_type in enumerate(['atac', 'rna']):

        axis[i].scatter(
            corrs[data_type], corrs['screen'], alpha=0.5, s=25,
            color=cmap(norm(corrs['color'])), rasterized=True)
        axis[i].scatter(
            corrs.loc[(corrs.color == 3), data_type], corrs.loc[(corrs.color == 3), 'screen'], alpha=0.5, s=25,
            color="black", rasterized=True)
        axis[i].set_xlabel(data_type.swapcase())
        axis[i].set_ylabel("Screen")
    fig.savefig(os.path.join(output_dir, analysis_name + ".correlation_with_ATAC_RNA.svg"), bbox_inches="tight", dpi=300)


    # Correlation of fold-changes
    atac_corr = pd.pivot_table(atac_analysis.differential_results.reset_index(), index='region', columns='comparison_name', values='log2FoldChange').corr()
    atac_corr.index.name = "k1"
    atac_corr.columns.name = "k2"
    atac_corr_m = pd.melt(atac_corr.reset_index(), id_vars='k1', value_name='atac')
    rna_corr = pd.pivot_table(rnaseq_analysis.differential_results.reset_index(), index='gene_name', columns='comparison_name', values='log2FoldChange').corr()
    rna_corr.index.name = "k1"
    rna_corr.columns.name = "k2"
    rna_corr_m = pd.melt(rna_corr.reset_index(), id_vars='k1', value_name='rna')


    corrs = (screen_m.set_index(['k1', 'k2'])
        .join(atac_corr_m.set_index(['k1', 'k2']))
        .join(rna_corr_m.set_index(['k1', 'k2'])))

    baf_members = ["ARID1A", "ARID1B"]
    pbaf_members = ["ARID2", "PBRM1", "BRD9", "PHF10", "DPF1"]

    corrs['baf'] = (corrs.index.get_level_values("k1").isin(baf_members).astype(int) + corrs.index.get_level_values("k2").isin(baf_members).astype(int))
    corrs['pbaf'] = -(corrs.index.get_level_values("k1").isin(pbaf_members).astype(int) + corrs.index.get_level_values("k2").isin(pbaf_members).astype(int))
    corrs['color'] = corrs['baf'] + corrs['pbaf']
    corrs.loc[(corrs['baf'] == 1) & (corrs['pbaf'] == -1), 'color'] = 3

    fig, axis = plt.subplots(1, 2, figsize=(8, 4))
    cmap = plt.get_cmap("coolwarm")
    cmap2 = plt.get_cmap("PuOr")
    norm = matplotlib.colors.Normalize(vmin=-2, vmax=2)
    for i, data_type in enumerate(['atac', 'rna']):

        axis[i].scatter(
            corrs[data_type], corrs['screen'], alpha=0.5, s=25,
            color=cmap(norm(corrs['color'])), rasterized=True)
        axis[i].scatter(
            corrs.loc[(corrs.color == 3), data_type], corrs.loc[(corrs.color == 3), 'screen'], alpha=0.5, s=25,
            color="black", rasterized=True)
        axis[i].set_xlabel(data_type.swapcase())
        axis[i].set_ylabel("Screen")
    fig.savefig(os.path.join(output_dir, analysis_name + ".correlation_with_ATAC_RNA.fold_changes.svg"), bbox_inches="tight", dpi=300)


def cross_cell_type_comparison(
        cell_types=["A549", "H2122", "HAP1", "OV90"]):

    for cell_type in cell_types:
        diff = pd.read_csv(
            os.path.join(
                "results{}".format("_" + cell_type if cell_type != "HAP1" else ""),
                "differential_analysis_ATAC-seq",
                "differential_analysis.deseq_result.ARID1A.csv"),
            index_col=0)
    
        if cell_type == "HAP1":
            up = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] > 1) ].index
            down = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] < -1) ].index
        else:
            up = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] > 0) ].index
            down = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] < 0) ].index

        # to bed
        up_bed = pd.DataFrame([
            map(lambda x: x[0], up.str.split(":")),
            map(lambda x: x[1].split("-")[0], up.str.split(":")),
            map(lambda x: x[1].split("-")[1], up.str.split(":"))], index=['chrom', 'start', 'end']).T
        up_bed['start'] = up_bed['start'].astype(int)
        up_bed['end'] = up_bed['end'].astype(int)
        up_bed = up_bed.sort_values(['chrom', 'start', 'end'])
        down_bed = pd.DataFrame([
            map(lambda x: x[0], down.str.split(":")),
            map(lambda x: x[1].split("-")[0], down.str.split(":")),
            map(lambda x: x[1].split("-")[1], down.str.split(":"))], index=['chrom', 'start', 'end']).T
        down_bed['start'] = down_bed['start'].astype(int)
        down_bed['end'] = down_bed['end'].astype(int)
        down_bed = down_bed.sort_values(['chrom', 'start', 'end'])


    cell_type = "A549"
    diff = pd.read_csv(
        os.path.join(
            "results{}".format("_" + cell_type.lower() if cell_type != "HAP1" else ""),
            "differential_analysis_ATAC-seq",
            "differential_analysis.deseq_result.REC.csv"),
        index_col=0)

    if cell_type == "HAP1":
        up = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] > 1) ].index
        down = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] < -1) ].index
    else:
        up = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] > 0) ].index
        down = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] < 0) ].index

    # to bed
    a549_up_bed = pd.DataFrame([
        map(lambda x: x[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[1], up.str.split(":"))], index=['chrom', 'start', 'end']).T
    a549_up_bed['start'] = a549_up_bed['start'].astype(int)
    a549_up_bed['end'] = a549_up_bed['end'].astype(int)
    a549_up_bed = a549_up_bed.sort_values(['chrom', 'start', 'end'])
    a549_down_bed = pd.DataFrame([
        map(lambda x: x[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[1], down.str.split(":"))], index=['chrom', 'start', 'end']).T
    a549_down_bed['start'] = a549_down_bed['start'].astype(int)
    a549_down_bed['end'] = a549_down_bed['end'].astype(int)
    a549_down_bed = a549_down_bed.sort_values(['chrom', 'start', 'end'])


    cell_type = "H2122"
    diff = pd.read_csv(
        os.path.join(
            "results{}".format("_" + cell_type.lower() if cell_type != "HAP1" else ""),
            "differential_analysis_ATAC-seq",
            "differential_analysis.deseq_result.ARID1A.csv"),
        index_col=0)

    if cell_type == "HAP1":
        up = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] > 1) ].index
        down = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] < -1) ].index
    else:
        up = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] > 0) ].index
        down = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] < 0) ].index

    # to bed
    h2122_up_bed = pd.DataFrame([
        map(lambda x: x[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[1], up.str.split(":"))], index=['chrom', 'start', 'end']).T
    h2122_up_bed['start'] = h2122_up_bed['start'].astype(int)
    h2122_up_bed['end'] = h2122_up_bed['end'].astype(int)
    h2122_up_bed = h2122_up_bed.sort_values(['chrom', 'start', 'end'])
    h2122_down_bed = pd.DataFrame([
        map(lambda x: x[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[1], down.str.split(":"))], index=['chrom', 'start', 'end']).T
    h2122_down_bed['start'] = h2122_down_bed['start'].astype(int)
    h2122_down_bed['end'] = h2122_down_bed['end'].astype(int)
    h2122_down_bed = h2122_down_bed.sort_values(['chrom', 'start', 'end'])


    cell_type = "HAP1"
    diff = pd.read_csv(
        os.path.join(
            "results{}".format("_" + cell_type.lower() if cell_type != "HAP1" else ""),
            "differential_analysis_ATAC-seq",
            "differential_analysis.deseq_result.ARID1A.csv"),
        index_col=0)

    if cell_type == "HAP1":
        up = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] > 1) ].index
        down = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] < -1) ].index
    else:
        up = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] > 0) ].index
        down = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] < 0) ].index

    # to bed
    hap1_up_bed = pd.DataFrame([
        map(lambda x: x[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[1], up.str.split(":"))], index=['chrom', 'start', 'end']).T
    hap1_up_bed['start'] = hap1_up_bed['start'].astype(int)
    hap1_up_bed['end'] = hap1_up_bed['end'].astype(int)
    hap1_up_bed = hap1_up_bed.sort_values(['chrom', 'start', 'end'])
    hap1_down_bed = pd.DataFrame([
        map(lambda x: x[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[1], down.str.split(":"))], index=['chrom', 'start', 'end']).T
    hap1_down_bed['start'] = hap1_down_bed['start'].astype(int)
    hap1_down_bed['end'] = hap1_down_bed['end'].astype(int)
    hap1_down_bed = hap1_down_bed.sort_values(['chrom', 'start', 'end'])

    cell_type = "OV90"
    diff = pd.read_csv(
        os.path.join(
            "results{}".format("_" + cell_type.lower() if cell_type != "HAP1" else ""),
            "differential_analysis_ATAC-seq",
            "differential_analysis.deseq_result.ARID1A.csv"),
        index_col=0)

    if cell_type == "HAP1":
        up = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] > 1) ].index
        down = diff[(diff['padj'] < 0.01) & (diff['log2FoldChange'] < -1) ].index
    else:
        up = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] > 0) ].index
        down = diff[(diff['padj'] < 0.05) & (diff['log2FoldChange'] < 0) ].index

    # to bed
    ov90_up_bed = pd.DataFrame([
        map(lambda x: x[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[0], up.str.split(":")),
        map(lambda x: x[1].split("-")[1], up.str.split(":"))], index=['chrom', 'start', 'end']).T
    ov90_up_bed['start'] = ov90_up_bed['start'].astype(int)
    ov90_up_bed['end'] = ov90_up_bed['end'].astype(int)
    ov90_up_bed = ov90_up_bed.sort_values(['chrom', 'start', 'end'])
    ov90_down_bed = pd.DataFrame([
        map(lambda x: x[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[0], down.str.split(":")),
        map(lambda x: x[1].split("-")[1], down.str.split(":"))], index=['chrom', 'start', 'end']).T
    ov90_down_bed['start'] = ov90_down_bed['start'].astype(int)
    ov90_down_bed['end'] = ov90_down_bed['end'].astype(int)
    ov90_down_bed = ov90_down_bed.sort_values(['chrom', 'start', 'end'])


    a549_up_bed.to_csv("a549_up.bed", sep="\t", header=None, index=False)
    a549_down_bed.to_csv("a549_down.bed", sep="\t", header=None, index=False)
    h2122_up_bed.to_csv("h2122_up.bed", sep="\t", header=None, index=False)
    h2122_down_bed.to_csv("h2122_down.bed", sep="\t", header=None, index=False)
    hap1_up_bed.to_csv("hap1_up.bed", sep="\t", header=None, index=False)
    hap1_down_bed.to_csv("hap1_down.bed", sep="\t", header=None, index=False)
    ov90_up_bed.to_csv("ov90_up.bed", sep="\t", header=None, index=False)
    ov90_down_bed.to_csv("ov90_down.bed", sep="\t", header=None, index=False)


    """bedtools multiinter -i \
        a549_up.bed a549_down.bed \
        h2122_up.bed h2122_down.bed \
        hap1_up.bed hap1_down.bed \
        ov90_up.bed ov90_down.bed \
        -names \
        a549_up a549_down \
        h2122_up h2122_down \
        hap1_up hap1_down \
        ov90_up ov90_down \
        > cross_cell_type_overlap.bed"""

    overlap = pd.read_csv(os.path.join("cross_cell_type_overlap.bed"), sep="\t", header=None)

    # plot overlap
    sns.distplot(overlap[3])

    # get peaks with two cell types overlaping
    o = overlap.loc[overlap[3] == 2]

    # from those, count
    o['down'] = o[4].str.count("down")
    o['up'] = o[4].str.count("down")

    o['disagreement'] = (o[4].str.contains("down") & o[4].str.contains("up")).astype(int)

    sns.distplot(o['disagreement'])


def enrichment_network():
    # network of enrichment terms
    corr = enrichr_pivot.corr()
    corr['from'] = corr.columns
    net = pd.melt(corr.reset_index(), id_vars=["from"])

    net2 = net[((6 ** net['value']) > 4) & (net['value'] != 1)]

    net3 = pd.merge(net2, enrichr_pivot.T.reset_index(), left_on="from", right_on="description").drop("description_y", axis=1).rename(columns={"description_x": "description"})

    import networkx as nx
    g = nx.from_pandas_dataframe(net3, "from", "description", ["value"])
    net3 = net3.set_index("from")
    for attr in net3.columns[2:]:
        nx.set_node_attributes(g, attr, {k: float(v) for k, v in net3[attr].to_dict().items()})

    nx.write_graphml(g, 'net_go.graphml')


def accessibility_expression(
        atac_analysis,
        rnaseq_analysis,
        acce_dir="deseq_knockout",
        expr_dir="deseq_expression_knockout",
        trait="knockout",
        output_dir="results/cross_datatype_comparison",
        output_prefix="cross_datatype_comparison"
):
    """
    """
    from scipy.stats import pearsonr
    from statsmodels.nonparametric.smoothers_lowess import lowess

    def signed_max(x, f=0.66):
        """
        Return maximum or minimum depending on the sign of the majority of values.
        If there isn't a clear majority (at least `f` fraction in one side), return mean of values.
        """
        obs_f = max(sum(x < 0) / float(len(x)), sum(x > 0) / float(len(x)))
        if obs_f >= f:
            if sum(x < 0) > sum(x > 0):
                return min(x)
            else:
                return max(x)
        else:
            return np.mean(x)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Accessibility fold-changes
    acce = atac_analysis.differential_results

    # annotate reg. elements with genes
    g = atac_analysis.gene_annotation
    acce = acce.join(
        g['gene_name'].str.split(",")
        .apply(pd.Series).stack()
        .reset_index(drop=True, level=1)
        .reset_index().set_index("index")
    ).rename(columns={0: "gene_name"})
    # get signed max for each gene and comparison
    acce_fc = pd.pivot_table(data=acce, values="log2FoldChange", index="gene_name", columns="comparison_name", aggfunc=signed_max)
    acce_fc.to_csv(os.path.join(output_dir, "accessibility.fold_changes.signed_max.per_gene.csv"), index=True)
    # acce_fc = pd.read_csv(os.path.join(output_dir, "accessibility.fold_changes.signed_max.per_gene.csv"), index_col=0)

    # Expression fold-changes
    expr = rnaseq_analysis.differential_results
    expr.index.name = "gene_name"
    expr_fc = pd.pivot_table(data=expr.reset_index(), values="log2FoldChange", index="gene_name", columns="comparison_name", aggfunc=np.mean)
    expr_fc.to_csv(os.path.join(output_dir, "expression.fold_changes.signed_max.per_gene.csv"), index=True)
    # expr_fc = pd.read_csv(os.path.join(output_dir, "expression.fold_changes.signed_max.per_gene.csv"), index_col=0)

    # match indexes (genes quantified)
    expr_fc = expr_fc.ix[acce_fc.index].dropna()
    acce_fc = acce_fc.ix[expr_fc.index].dropna()

    # Plot correlation of fold_changes
    joint = expr_fc.join(acce_fc, lsuffix="_RNA-seq", rsuffix="_ATAC-seq")
    c = joint.corr()

    # all vs all
    g = sns.clustermap(
        c,
        xticklabels=True, yticklabels=True,
        robust=True, col_cluster=False, row_cluster=False, cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".fold_change.signed_max.correlation.all.ordered.svg"))

    g = sns.clustermap(
        c,
        xticklabels=True, yticklabels=True,
        robust=True, col_cluster=True, row_cluster=True, cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".fold_change.signed_max.correlation.all.clustered.svg"))

    # one vs other
    q = c.loc[
            c.index.str.contains("RNA-seq"),
            c.columns.str.contains("ATAC-seq")].sort_index(0).sort_index(1)
    v = q.abs().values.flatten().max()
    v += (v / 5.)
    g = sns.clustermap(
        q,
        xticklabels=True, yticklabels=True,
        cmap="RdBu_r", vmin=-v, vmax=v, col_cluster=False, row_cluster=False, robust=False, cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".fold_change.signed_max.correlation.cross-data_types.ordered.svg"))
    g = sns.clustermap(
        q,
        xticklabels=True, yticklabels=True,
        cmap="RdBu_r", vmin=-v, vmax=v, col_cluster=True, row_cluster=True, robust=False, cbar_kws={"label": "Pearson correlation"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".fold_change.signed_max.correlation.cross-data_types.clustered.svg"))


    # Plot fold-changes in accessibility vs expression
    # Using gene-level measurements for ATAC-seq
    n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    for i, ko in enumerate(expr_fc.columns):
        print(ko)
        # get values
        a = expr_fc[ko]
        b = acce_fc[ko]

        # correlate
        r, p = pearsonr(a, b)

        # # fit lowess
        # if i == 0:
        #     l = lowess(a, b, return_sorted=False)
        #     # pass
        # dist = pd.DataFrame([abs(a - l), abs(b - l)]).mean()

        # plot scatter
        axis.flat[i].scatter(a, b, s=0.5, rasterized=True, color="gray", alpha=0.1)  #, color=plt.cm.GnBu(dist))

        # add title and correlation values
        axis.flat[i].set_title(ko)
        axis.flat[i].text(1, -5, s="r = {:.3f}".format(r))  # \np = {:.3f}

        # Color significant differently
        sig = expr_fc[(abs(a) > 1) & (abs(b) > 1)].index
        # sig = dist[dist > np.percentile(dist, 99)].index

        axis.flat[i].scatter(a.ix[sig], b.ix[sig], s=1)  # , color=sns.color_palette("Set3")[3])

    axis[2, 0].set_ylabel("log2 fold-change (ATAC-seq)")
    axis[4, 2].set_xlabel("log2 fold-change (RNA-seq)")
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.99_perc.no_lowess_color.svg"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.99_perc.png"), bbox_inches="tight", dpi=300)

    # By fixed bins of fold-change
    step = 0.25
    range_max = 3.0
    phase = 0.1
    # bins = zip(-np.arange(step, range_max + step, step)[::-1], -np.arange(0, range_max, step)[::-1])
    bins = zip(np.arange(0, range_max, step), np.arange(step, range_max + step, step))
    bins += zip(np.arange(phase, range_max, step + phase), np.arange(step + phase, range_max + step + phase, step + phase))
    bins += zip(np.arange(phase * 1.25, range_max, step + phase * 1.25), np.arange(step + phase * 1.25, range_max + step + phase * 1.25, step + phase * 1.25))
    bins += zip(np.arange(phase * 1.5, range_max, step + phase * 1.5), np.arange(step + phase * 1.5, range_max + step + phase * 1.5, step + phase * 1.5))
    bins += zip(np.arange(phase * 1.75, range_max, step + phase * 1.75), np.arange(step + phase * 1.75, range_max + step + phase * 1.75, step + phase * 1.75))
    bins += zip(np.arange(phase * 2.0, range_max, step + phase * 2.0), np.arange(step + phase * 2.0, range_max + step + phase * 2.0, step + phase * 2.0))
    bins += zip(np.arange(phase * 2.25, range_max, step + phase * 2.25), np.arange(step + phase * 2.25, range_max + step + phase * 2.25, step + phase * 2.25))
    bins += zip(np.arange(phase * 2.5, range_max, step + phase * 2.5), np.arange(step + phase * 2.5, range_max + step + phase * 2.5, step + phase * 2.5))
    bins += zip(np.arange(phase * 2.75, range_max, step + phase * 2.75), np.arange(step + phase * 2.75, range_max + step + phase * 2.75, step + phase * 2.75))
    bins += zip(np.arange(phase * 3.0, range_max, step + phase * 3.0), np.arange(step + phase * 3.0, range_max + step + phase * 3.0, step + phase * 3.0))
    # bins = bins[:5] + bins[-5:]

    dists = pd.DataFrame()
    means = pd.DataFrame()
    for i, ko in enumerate(acce_fc.columns):
        # get values
        try:
            a = abs(expr_fc[ko])
            b = abs(acce_fc[ko])
        except KeyError:
            continue
        print(ko)
        for j, (l1, l2) in enumerate(bins):
            sel_genes = b[(b > l1) & (b < l2)].index
            means = means.append(pd.Series([a[sel_genes].dropna().mean(), ko, l1, l2], index=["mean", "knockout", "min", "max"]), ignore_index=True)
            d = a[sel_genes].dropna().to_frame()
            d.columns = ['change']
            d["knockout"] = ko
            d['min'] = l1
            d['max'] = l2
            dists = dists.append(d)

    means2 = means.groupby(['max', 'min'])['mean'].mean().reset_index()
    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))
    axis.scatter(means2.loc[:, "min"], means2.loc[:, "mean"])
    axis.set_ylabel("abs log2 fold-change (RNA-seq)")
    axis.set_xlabel("abs log2 fold-change (ATAC-seq)")
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.absolute.cross_knockouts.svg"), bbox_inches="tight", dpi=300)

    means2 = dists[~dists['knockout'].str.contains("HIRA")].groupby(['max', 'min'])['change'].mean().reset_index()
    means2 = means2[means2['max'] < 2.5] 
    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))
    axis.scatter(means2.loc[:, "min"], means2.loc[:, "change"])
    axis.axhline(0, alpha=0.2, color="black", linestyle="--")
    axis.set_ylabel("abs log2 fold-change (RNA-seq)")
    axis.set_xlabel("abs log2 fold-change (ATAC-seq)")
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.absolute.cross_knockouts.svg"), bbox_inches="tight", dpi=300)

    dists2 = dists[dists['min'] <= 3.0]
    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))
    sns.violinplot(dists2.loc[:, "min"], dists2.loc[:, "change"], trim=True, cut=0, ax=axis)
    axis.set_ylabel("abs log2 fold-change (RNA-seq)")
    axis.set_xlabel("abs log2 fold-change (ATAC-seq)")
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.absolute.cross_knockouts.violinplot.svg"), bbox_inches="tight", dpi=300)

    dists2 = dists[dists['min'] <= 3.0]
    dists2 = dists2.groupby(['knockout', 'min', 'max']).mean().reset_index()
    n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    for i, ko in enumerate(acce_fc.columns):
        try:
            axis.flat[i].scatter(dists2.loc[dists2['knockout'] == ko, "min"], dists2.loc[dists2['knockout'] == ko, "change"])
        except IndexError:
            continue
        axis.flat[i].set_title(ko)
    fig.savefig(os.path.join(output_dir, output_prefix + ".fold_changes.signed_max.absolute.each_knockout.svg"), bbox_inches="tight", dpi=300)

    # Correlate gene-level enrichments from ATAC-seq and RNA-seq
    atac_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(atac_analysis.results_dir, "ATAC-seq"), "differential_analysis.enrichr.csv"))
    rnaseq_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, "RNA-seq"), "differential_analysis.enrichr.csv"))

    for gene_set_library in atac_enrichment_table['gene_set_library'].unique():
        fig, axis = plt.subplots(21, 2, figsize=(2 * 4, 21 * 4))
        for i, comparison_name in enumerate(atac_enrichment_table["comparison_name"].unique()):
            for j, direction in enumerate(atac_enrichment_table["direction"].unique()):
                print(gene_set_library, comparison_name, direction)
                a = atac_enrichment_table[
                    (atac_enrichment_table['gene_set_library'] == gene_set_library) &
                    (atac_enrichment_table['comparison_name'] == comparison_name) &
                    (atac_enrichment_table['direction'] == direction)]
                r = rnaseq_enrichment_table[
                    (rnaseq_enrichment_table['gene_set_library'] == gene_set_library) &
                    (rnaseq_enrichment_table['comparison_name'] == comparison_name) &
                    (rnaseq_enrichment_table['direction'] == direction)]
                joint = a.set_index("description")['p_value'].to_frame(name="ATAC-seq")
                joint = joint.join(r.set_index("description")['p_value'].to_frame(name="RNA-seq")).dropna().drop_duplicates()

                joint = -np.log10(joint)
                axis[i, j].scatter(joint['ATAC-seq'], joint['RNA-seq'], alpha=0.5, s=3, rasterized=True)
                axis[i, j].set_title("{} - {}".format(comparison_name, direction))

                for term in joint.mean(1).sort_values().tail(top_n).index:
                    axis[i, j].text(x=joint.loc[term, 'ATAC-seq'], y=joint.loc[term, 'RNA-seq'], s=term)
        fig.savefig(os.path.join(output_dir, output_prefix + ".enrichments.{}.scatter+text.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

    # Now directions combined
    for gene_set_library in atac_enrichment_table['gene_set_library'].unique():
        fig, axis = plt.subplots(21, 2, figsize=(2 * 4, 21 * 4))
        for i, comparison_name in enumerate(atac_enrichment_table["comparison_name"].unique()):
            for j, direction in enumerate(atac_enrichment_table["direction"].unique()):
                print(gene_set_library, comparison_name, direction)
                au = -np.log10(atac_enrichment_table[
                    (atac_enrichment_table['gene_set_library'] == gene_set_library) &
                    (atac_enrichment_table['comparison_name'] == comparison_name) &
                    (atac_enrichment_table['direction'] == "up")].set_index("description")['p_value'])
                ad = -np.log10(atac_enrichment_table[
                    (atac_enrichment_table['gene_set_library'] == gene_set_library) &
                    (atac_enrichment_table['comparison_name'] == comparison_name) &
                    (atac_enrichment_table['direction'] == "down")].set_index("description")['p_value'])

                ru = rnaseq_enrichment_table[
                    (rnaseq_enrichment_table['gene_set_library'] == gene_set_library) &
                    (rnaseq_enrichment_table['comparison_name'] == comparison_name) &
                    (rnaseq_enrichment_table['direction'] == "up")]
                rd = rnaseq_enrichment_table[
                    (rnaseq_enrichment_table['gene_set_library'] == gene_set_library) &
                    (rnaseq_enrichment_table['comparison_name'] == comparison_name) &
                    (rnaseq_enrichment_table['direction'] == "up")]
                
                
                joint = a.set_index("description")['p_value'].to_frame(name="ATAC-seq")
                joint = joint.join(r.set_index("description")['p_value'].to_frame(name="RNA-seq")).dropna().drop_duplicates()

                joint = -np.log10(joint)
                axis[i, j].scatter(joint['ATAC-seq'], joint['RNA-seq'], alpha=0.5, s=3, rasterized=True)
                axis[i, j].set_title("{} - {}".format(comparison_name, direction))

                for term in joint.mean(1).sort_values().tail(top_n).index:
                    axis[i, j].text(x=joint.loc[term, 'ATAC-seq'], y=joint.loc[term, 'RNA-seq'], s=term)
        fig.savefig(os.path.join(output_dir, output_prefix + ".enrichments.{}.scatter+text.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)




    # ARID2 & SMARCA4
    sig = rnaseq_analysis.differential_results[
        (rnaseq_analysis.differential_results['padj'] < alpha) &
        (rnaseq_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)]

    a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "up")].index
    b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "down")].index
    a[a.isin(b.tolist())].to_series().to_clipboard(index=False)
    a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "down")].index
    b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "up")].index
    a[a.isin(b.tolist())].to_series().to_clipboard(index=False)

    a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "up")].index
    b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "up")].index
    a[a.isin(b.tolist())].to_series().to_clipboard(index=False)

    a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "down")].index
    b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "down")].index
    a[a.isin(b.tolist())].to_series().to_clipboard(index=False)


    # Find the interaction
    rnaseq_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, "RNA-seq"), "differential_analysis.enrichr.csv"))

    q = rnaseq_enrichment_table.loc[
        # (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
        (rnaseq_enrichment_table['gene_set_library'].isin(["NCI-Nature_2016"])) &
        (rnaseq_enrichment_table['direction'] == 'down') &
        (rnaseq_enrichment_table['p_value'] < 0.05) &
        (
            rnaseq_enrichment_table['description'].str.contains("E2F") |
            rnaseq_enrichment_table['description'].str.contains("cell cycle", case=False)), :]

    genes = q.loc[:, "genes"]

    cc_genes = list(set(genes.str.replace("[", " ").str.replace(
        ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # clustermap
    g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature&WikiPathways.clustermap.svg"), bbox_inches="tight", dpi=300)


    g = sns.clustermap(rnaseq_analysis.expression_annotated.loc[cc_genes, :].dropna().T.groupby("knockout").mean().T, z_score=0, rasterized=True, square=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature.group.clustermap.svg"), bbox_inches="tight", dpi=300)



    # cc_genes = ["CCND3", "RBL1", "CCND2", "CDK2", "CDC25A"]
    q = rnaseq_enrichment_table.loc[
        (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
        # (rnaseq_enrichment_table['gene_set_library'] == 'NCI-Nature_2016') &
        (rnaseq_enrichment_table['direction'] == 'down') &
        (rnaseq_enrichment_table['p_value'] < 0.05), :]

    genes = q.loc[:, "genes"]

    cc_genes = list(set(genes.str.replace("[", " ").str.replace(
        ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # clustermap
    g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.svg"), bbox_inches="tight", dpi=300)

    # investigate genes in second cluster (ARID2-specific)
    clusters = scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=2, criterion='maxclust')
    # plot again just to confirm clusters
    g2 = sns.clustermap(
        rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True,
        row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage, row_colors=plt.get_cmap("Paired")(clusters))
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.clusters_labeled.svg"), bbox_inches="tight", dpi=300)

    pbaf_genes = pd.Series(g.data.index).iloc[clusters == pd.Series(clusters).value_counts().argmin()].sort_values()

    g3 = sns.clustermap(rnaseq_analysis.expression.loc[pbaf_genes, :], z_score=0, rasterized=True, metric="correlation")
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90)
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_yticklabels(), rotation=0)
    g3.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.pbaf_genes.svg"), bbox_inches="tight", dpi=300)


    # Cell cycle signature
    joint = rnaseq_analysis.expression_annotated.loc[cc_genes, :].T.groupby('knockout').mean().mean(1)
    joint_z = pd.Series(scipy.stats.zscore(joint), index=joint.index)

    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4))
    axis[0].scatter(joint.rank(), joint, cmap="RdBu", vmin=joint.min(), vmax=joint.max())
    for ko in joint.index:
        axis[0].text(joint.rank()[ko], joint[ko], ko)
    axis[0].axhline(joint.mean(), linestyle='--', color="black")
    axis[0].set_xlabel("Cell cycle signature (Rank)")
    axis[0].set_ylabel("Cell cycle signature score")
    axis[1].scatter(joint_z.rank(), joint_z, cmap="RdBu", vmin=joint_z.min(), vmax=joint_z.max())
    for ko in joint_z.index:
        axis[1].text(joint_z.rank()[ko], joint_z[ko], ko)
    axis[1].axhline(0, linestyle='--', color="black")
    axis[1].set_xlabel("Cell cycle signature (Rank)")
    axis[1].set_ylabel("Cell cycle signature (Z-score)")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.mean_zscore.rank_vs_zscore.svg"), bbox_inches="tight")


    # Classic cell cycle signature
    cc_genes = pd.read_table("regev_lab_cell_cycle_genes.txt", header=None, squeeze=True)

    # clustermap
    g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_regev_signature.clustermap.svg"), bbox_inches="tight", dpi=300)


def transcription_factor_accessibility():
    from glob import glob
    import pybedtools

    bed_dir = "/home/arendeiro/resources/regions/LOLACore/hg19/encode_tfbs/regions/"
    output_dir = "results"

    tfs = [
        "E2F1", "E2F4", "E2F6",
        "POL2", "CTCF", "TAF1", "SP1", "SP2", "ELK1", "ELK4", 
    ]

    all_res = pd.DataFrame()
    for factor_name in tfs:
        print(factor_name)
        # get consensus TF regions from LOLA database
        transcription_factor_regions_sets = glob(
            os.path.join(bed_dir, "*{}*".format(factor_name.capitalize())))[:5]
        bt = pybedtools.BedTool(transcription_factor_regions_sets[0]).sort()
        for fn in transcription_factor_regions_sets[1:]:
            bt = bt.cat(pybedtools.BedTool(fn).sort().merge())
        bt = bt.merge()
        bt.saveas(os.path.join("data", "external", "TF_binding_sites" + factor_name + ".bed"))

        # get regions overlapping with TF sites
        transcription_factor_r = analysis.sites.intersect(bt, wa=True).to_dataframe(names=['chrom', 'start', 'end'])
        transcription_factor_r.index = transcription_factor_r['chrom'] + ":" + transcription_factor_r['start'].astype(str) + "-" + transcription_factor_r['end'].astype(str)
        transcription_factor_a = analysis.accessibility.loc[transcription_factor_r.index].dropna()

        # group regions by quantile of accessibility across all experiments
        lower = 0.0
        upper = 1.0
        n_groups = 10
        r = np.arange(lower, upper + (upper / n_groups), upper / n_groups)
        mean = transcription_factor_a.mean(axis=1)

        res = pd.DataFrame()
        for l_quantile, u_quantile in zip(r, r[1:]):
            i = mean[(mean.quantile(l_quantile) > mean) & (mean < mean.quantile(u_quantile))].index

            m = transcription_factor_a.loc[i, :].mean()
            m.index = m.index.get_level_values("sample_name")
            m['upper_quantile'] = u_quantile
            res = res.append(m, ignore_index=True)
        res = res.set_index('upper_quantile')
        i = pd.DataFrame(map(pd.Series, res.columns.str.split("_")))
        res.columns = pd.MultiIndex.from_arrays(i[[2, 3]].values.T, names=['knockout', 'replicate'])

        res = res.sort_index(axis=1, level=['knockout', 'replicate'], sort_remaining=False)

        d = res.dropna().T.groupby(level=['knockout', 'replicate']).mean().mean(1).reset_index()
        d["transcription_factor"] = factor_name
        d["binding_sites"] = transcription_factor_a.shape[0]
        d = d.rename(columns={0: "accessibility"})
        all_res = all_res.append(d, ignore_index=True)


        g = sns.clustermap(res.dropna().T, col_cluster=False, z_score=1, rasterized=True, figsize=(res.shape[0] * 0.12, res.shape[1] * 0.12), row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, "{}_binding.per_quantile.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")

        res_mean = res.dropna().T.groupby(level=['knockout', 'replicate']).mean()
        g = sns.clustermap(res_mean.dropna(), col_cluster=False, z_score=1, rasterized=False, square=True, row_cluster=False)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.savefig(os.path.join(output_dir, "{}_binding.per_quantile.mean_replicates.sorted.svg".format(factor_name)), dpi=300, bbox_inches="tight")


    # Get "background" accessibility per cell type per replicate
    # group regions by quantile of accessibility across all experiments
    lower = 0.0
    upper = 1.0
    n_groups = 10
    r = np.arange(lower, upper + (upper / n_groups), upper / n_groups)
    mean = analysis.accessibility.mean(axis=1)

    res = pd.DataFrame()
    for l_quantile, u_quantile in zip(r, r[1:]):
        i = mean[(mean.quantile(l_quantile) > mean) & (mean < mean.quantile(u_quantile))].index

        m = analysis.accessibility.loc[i, :].mean()
        m.index = m.index.get_level_values("sample_name")
        m['upper_quantile'] = u_quantile
        res = res.append(m, ignore_index=True)
    res = res.set_index('upper_quantile')
    i = pd.DataFrame(map(pd.Series, res.columns.str.split("_")))
    res.columns = pd.MultiIndex.from_arrays(i[[2, 3]].values.T, names=['knockout', 'replicate'])

    res = res.sort_index(axis=1, level=['knockout', 'replicate'], sort_remaining=False)
    d = res.dropna().T.groupby(level=['knockout', 'replicate']).mean().mean(1).reset_index()
    d["transcription_factor"] = "background"
    d["binding_sites"] = analysis.accessibility.shape[0]
    d = d.rename(columns={0: "accessibility"})
    all_res = all_res.append(d, ignore_index=True)

    # Save
    all_res.to_csv(os.path.join(output_dir, "all_factor_binding.csv"), index=False)
    all_res = pd.read_csv(os.path.join(output_dir, "all_factor_binding.csv"))


    # Normalize to background
    for knockout in all_res['knockout'].drop_duplicates():
        for tf in all_res['transcription_factor'].drop_duplicates():
            s = all_res.loc[(all_res['knockout'] == knockout) & (all_res['transcription_factor'] == tf)].set_index(['knockout', 'replicate'])['accessibility']
            b = all_res.loc[(all_res['knockout'] == knockout) & (all_res['transcription_factor'] == "background")].set_index(['knockout', 'replicate'])['accessibility']
            all_res.loc[(all_res['knockout'] == knockout) & (all_res['transcription_factor'] == tf), 'norm_accessibility'] = (s - b).values  # ((s - b) / b.std()).values

    all_res.to_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"), index=False)
    all_res = pd.read_csv(os.path.join(output_dir, "all_factor_binding.normalized.csv"))

    # Plot
    piv = pd.pivot_table(data=all_res, index="knockout", columns='transcription_factor', values='accessibility').drop('background', axis=1)

    g = sns.clustermap(pd.DataFrame(scipy.stats.zscore(piv, axis=1), index=piv.index, columns=piv.columns), z_score=1)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "all_factor_binding.mean_replicates.heatmap.zscore.svg"), dpi=300, bbox_inches="tight")



def characterize_regions_structure(df, prefix, output_dir, universe_df=None, plot=True):
    # use all sites as universe
    if universe_df is None:
        universe_df = analysis.coverage_annotated

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # compare genomic regions and chromatin_states
    enrichments = pd.DataFrame()
    for i, var in enumerate(['genomic_region', 'chromatin_state']):
        # prepare:
        # separate comma-delimited fields:
        df_var = df[var].str.split(',').apply(pd.Series).stack().reset_index(drop=True, level=1)
        df_var.name = "test_set"
        universe_var = universe_df[var].str.split(',').apply(pd.Series).stack().reset_index(drop=True, level=1)
        universe_var.name = "universe_set"

        counts = df_var.value_counts().to_frame().join(universe_var.value_counts())
        counts.index.name = "region"

        # get also % of genome space "used"
        df["length"] = df["end"] - df["start"]
        universe_df["length"] = universe_df["end"] - universe_df["start"]

        df_size = df.join(df_var).groupby("test_set")["length"].sum()
        df_size = (df_size / pd.Series(analysis.genome_space)).dropna() * 100
        df_size.name = "df_space"
        counts = counts.join(df_size)
        universe_size = universe_df.join(universe_var).groupby("universe_set")["length"].sum()
        universe_size = (universe_size / pd.Series(analysis.genome_space)).dropna() * 100
        universe_size.name = "universe_space"
        counts = counts.join(universe_size)

        # transform to percentages
        counts["%_test_set"] = (counts["test_set"] / counts["test_set"].sum()) * 100
        counts["%_universe_set"] = (counts["universe_set"] / counts["universe_set"].sum()) * 100
        counts["%_df_space"] = (counts["df_space"] / counts["df_space"].sum()) * 100
        counts["%_universe_space"] = (counts["universe_space"] / counts["universe_space"].sum()) * 100

        # calculate change
        counts["fold_change"] = np.log2(counts["%_test_set"] / counts["%_universe_set"])
        counts["fold_size_change"] = np.log2(counts["%_df_space"] / counts["%_universe_space"])

        # sort for same order
        counts.sort_values('fold_size_change', inplace=True)

        if plot:
            g = sns.FacetGrid(data=pd.melt(counts.reset_index(), id_vars=["region"]), col="variable", col_wrap=3, sharex=True, sharey=False)
            g.map(sns.barplot, "region", "value")
            plt.savefig(os.path.join(output_dir, "{}_regions.{}.svg".format(prefix, var)), bbox_inches="tight")

        # append
        enrichments = enrichments.append(counts)

    # save
    enrichments.to_csv(os.path.join(output_dir, "%s_regions.region_enrichment.csv" % prefix), index=True)


def nucleosome_changes(analysis, samples):
    # select only ATAC-seq samples
    df = analysis.prj.sheet.df[analysis.prj.sheet.df["library"] == "ATAC-seq"]

    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    # nucleosomes per sample
    nucpos_metrics = {
        "z-score": 4,
        "nucleosome_occupancy_estimate": 5,
        "lower_bound_for_nucleosome_occupancy_estimate": 6,
        "upper_bound_for_nucleosome_occupancy_estimat": 7,
        "log_likelihood_ratio": 8,
        "normalized_nucleoatac_signal": 9,
        "cross-correlation_signal_value_before_normalization": 10,
        "number_of_potentially_nucleosome-sized_fragments": 11,
        "number_of_fragments_smaller_than_nucleosome_sized": 12,
        "fuzziness": 13
    }
    # nucleosome-free regions per sample
    nfrpos_metrics = {
        "mean_occupancy": 4,
        "minimum_upper_bound_occupancy": 5,
        "insertion_density": 6,
        "bias_density": 7,
    }
    for data_type, metrics in [("nucpos", nucpos_metrics), ("nfrpos", nfrpos_metrics)]:
        figs = list()
        axizes = list()
        for metric in metrics:
            fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharex=True, sharey=True)
            figs.append(fig)
            axizes.append(axis.flatten())

        counts = pd.Series()
        for i, group in enumerate(groups):
            print(data_type, group)
            s = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}.bed.gz".format(data_type)),
                sep="\t", header=None)
            counts[group] = s.shape[0]

            for j, (metric, col) in enumerate(metrics.items()):
                print(data_type, group, metric, col)
                sns.distplot(s[col - 1].dropna(), hist=False, ax=axizes[j][i])
                axizes[j][i].set_title(group)

        for i, metric in enumerate(metrics):
            sns.despine(figs[i])
            figs[i].savefig(os.path.join("results", "nucleoatac", "plots", "{}.{}.svg".format(data_type, metric)), bbox_inches="tight")

        fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 1 * 4))
        sns.barplot(counts, counts.index, orient="horizontal", ax=axis, color=sns.color_palette("colorblind")[0])
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
        axis.set_xlabel("Calls")
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.count_per_sample.svg".format(data_type)), bbox_inches="tight")

    # fragment distribution
    for data_type in ["InsertionProfile", "InsertSizes", "fragmentsizes"]:
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4))
        axis = axis.flatten()

        data = pd.DataFrame()
        for i, group in enumerate(groups):
            print(data_type, group)
            s = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}.txt".format(data_type)),
                sep="\t", header=None, squeeze=True, skiprows=5 if data_type == "fragmentsizes" else 0)
            if data_type == "InsertionProfile":
                a = (len(s.index) - 1) / 2.
                s.index = np.arange(-a, a + 1)
            if data_type == "fragmentsizes":
                s = s.squeeze()
            data[group] = s

            axis[i].plot(s)
            axis[i].set_title(group)
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.svg".format(data_type)), bbox_inches="tight")

        norm_data = data.apply(lambda x: x / data['ATAC-seq_HAP1_WT_C631'])
        if data_type == "fragmentsizes":
            norm_data = norm_data.loc[50:, :]
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharey=True)
        axis = axis.flatten()
        for i, group in enumerate(groups):
            print(data_type, group)
            axis[i].plot(norm_data[group])
            axis[i].set_title(group)
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.WT_norm.svg".format(data_type)), bbox_inches="tight")

    # Vplots and
    v_min, v_max = (105, 251)
    # Vplots over WT
    for data_type in ["VMat"]:
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4))
        fig2, axis2 = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharey=True)
        axis = axis.flatten()
        axis2 = axis2.flatten()

        group = 'ATAC-seq_HAP1_WT_C631'
        wt = pd.read_csv(
            os.path.join("results", "nucleoatac", group, group + ".{}".format(data_type)),
            sep="\t", header=None, skiprows=7)
        wt.index = np.arange(v_min, v_max)
        a = (len(wt.columns) - 1) / 2.
        wt.columns = np.arange(-a, a + 1)
        wt = wt.loc[0:300, :]

        for i, group in enumerate(groups):
            print(data_type, group)
            m = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}".format(data_type)),
                sep="\t", header=None, skiprows=7)
            m.index = np.arange(v_min, v_max)
            a = (len(m.columns) - 1) / 2.
            m.columns = np.arange(-a, a + 1)
            m = m.loc[0:300, :]

            n = m / wt

            axis[i].imshow(m.sort_index(ascending=False))
            axis[i].set_title(group)
            axis2[i].imshow(n.sort_index(ascending=False))
            axis2[i].set_title(group)
        sns.despine(fig)
        sns.despine(fig2)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.svg".format(data_type)), bbox_inches="tight")
        fig2.savefig(os.path.join("results", "nucleoatac", "plots", "{}.WT_norm.svg".format(data_type)), bbox_inches="tight")


def investigate_nucleosome_positions(self, samples, cluster=True):
    df = pd.DataFrame([s.as_series() for s in samples])
    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    def get_differential(diff, conditions_for, directions_for):
        # filter conditions
        d = diff[diff['comparison'].isin(conditions_for)]

        # filter directions
        d = pd.concat([d[f(d['log2FoldChange'], x)] for f, x in directions_for])

        # make bed format
        df = pd.Series(d.index.str.split(":")).apply(lambda x: pd.Series([x[0]] + x[1].split("-"))).drop_duplicates().reset_index(drop=True)
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(np.int64)
        df[2] = df[2].astype(np.int64)
        return df

    def center_window(bedtool, width=1000):
        chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

        sites = list()
        for site in bedtool:
            if site.chrom in chroms:
                mid = site.start + ((site.end - site.start) / 2)
                sites.append(pybedtools.Interval(site.chrom, mid - (width / 2), (mid + 1) + (width / 2)))
        return pybedtools.BedTool(sites)

    def center_series(series, width=1000):
        mid = series[1] + ((series[2] - series[1]) / 2)
        return pd.Series([series[0], mid - (width / 2), (mid + 1) + (width / 2)])

    # Regions to look at
    regions_pickle = os.path.join(self.results_dir, "nucleoatac", "all_types_of_regions.pickle")
    if os.path.exists():
        regions = pickle.load(open(regions_pickle, "rb"))
    else:
        regions = dict()
        # All accessible sites
        out = os.path.join(self.results_dir, "nucleoatac", "all_sites.bed")
        center_window(self.sites).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
        regions["all_sites"] = out

        # get bed file of promoter proximal sites
        promoter_sites = os.path.join(self.results_dir, "nucleoatac", "promoter_proximal.bed")
        self.coverage_annotated[self.coverage_annotated['distance'].astype(int) < 5000][['chrom', 'start', 'end']].to_csv(
            promoter_sites, index=False, header=None, sep="\t")
        regions["promoter"] = promoter_sites

        # All accessible sites - that are distal
        out = os.path.join(self.results_dir, "nucleoatac", "distal_sites.bed")
        center_window(self.sites.intersect(pybedtools.BedTool(promoter_sites), v=True, wa=True)).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
        regions["distal_sites"] = out

        # Differential sites
        diff = pd.read_csv(os.path.join("results", "deseq_knockout", "deseq_knockout.knockout.csv"), index_col=0)
        diff = diff[(diff["padj"] < 0.01) & (abs(diff["log2FoldChange"]) > 1.)]

        # Loosing accessibility with ARID1A/SMARCA4 KO
        less_ARID1ASMARCA4 = get_differential(
            diff,
            ['ARID1A-WT', 'SMARCA4-WT'],
            [(np.less_equal, -1), (np.less_equal, -1)]).apply(center_series, axis=1)
        out = os.path.join(self.results_dir, "nucleoatac", "diff_sites.less_ARID1ASMARCA4.bed")
        less_ARID1ASMARCA4.to_csv(out, index=False, header=None, sep="\t")
        regions["diff_sites.less_ARID1ASMARCA4"] = out

        # Gaining accessibility with ARID1A/SMARCA4 KO
        more_ARID1ASMARCA4 = get_differential(
            diff,
            ['ARID1A-WT', 'SMARCA4-WT'],
            [(np.greater_equal, 1), (np.greater_equal, 1)]).apply(center_series, axis=1)
        out = os.path.join(self.results_dir, "nucleoatac", "diff_sites.more_ARID1ASMARCA4.bed")
        more_ARID1ASMARCA4.to_csv(out, index=False, header=None, sep="\t")
        regions["diff_sites.more_ARID1ASMARCA4"] = out

        # Gaining accessibility with BCL7B/ARID1B KO
        more_BCL7BARID1B = get_differential(
            diff,
            ['BCL7B-WT', 'ARID1B-WT'],
            [(np.greater_equal, 1), (np.greater_equal, 1)]).apply(center_series, axis=1)
        out = os.path.join(self.results_dir, "nucleoatac", "diff_sites.more_BCL7BARID1B.bed")
        more_BCL7BARID1B.to_csv(out, index=False, header=None, sep="\t")
        regions["diff_sites.more_BCL7BARID1B"] = out

        # TFBSs
        tfs = ["CTCF", "BCL", "SMARC", "POU5F1", "SOX2", "NANOG", "TEAD4"]
        for tf in tfs:
            tf_bed = pybedtools.BedTool("/home/arendeiro/resources/genomes/hg19/motifs/TFs/{}.true.bed".format(tf))
            out = os.path.join(self.results_dir, "nucleoatac", "tfbs.%s.bed" % tf)
            center_window(tf_bed.intersect(self.sites, wa=True)).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
            regions["tfbs.%s" % tf] = out

        pickle.dump(regions, open(regions_pickle, "wb"))

    # Launch jobs
    for group in groups:
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)

        # Signals to measure in regions
        signal_files = [
            ("signal", os.path.join(self.data_dir, "merged", group + ".merged.sorted.bam")),
            ("nucleosome", os.path.join(self.data_dir, "merged", group + ".nucleosome_reads.bam")),
            ("nucleosome_free", os.path.join(self.data_dir, "merged", group + ".nucleosome_free_reads.bam")),
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz")),
            ("dyads", os.path.join("results", "nucleoatac", group, group + ".nucpos.bed.gz"))
        ]
        for region_name, bed_file in regions.items():
            for label, signal_file in signal_files:
                print(group, region_name, label)
                # run job
                run_coverage_job(bed_file, signal_file, label, ".".join([group, region_name, label]), output_dir, window_size=2001)
                # run vplot
                if label == "signal":
                    run_vplot_job(bed_file, signal_file, ".".join([group, region_name, label]), output_dir)

    # Collect signals
    signals = pd.DataFrame(columns=['group', 'region', 'label'])
    # signals = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"))

    # import re
    for group in groups:  # [g for g in groups if any([re.match(".*%s.*" % x, g) for x in ["C631", "HAP1_ARID1", "HAP1_SMARCA"]])]
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)
        signal_files = [
            ("signal", os.path.join(self.data_dir, "merged", group + ".merged.sorted.bam")),
            ("nucleosome", os.path.join(self.data_dir, "merged", group + ".nucleosome_reads.bam")),
            ("nucleosome_free", os.path.join(self.data_dir, "merged", group + ".nucleosome_free_reads.bam")),
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz")),
            ("dyads", os.path.join("results", "nucleoatac", group, group + ".nucpos.bed.gz"))
        ]
        for region_name, bed_file in regions.items():
            for label, signal_file in signal_files:  # [signal_files.items()[-1]]
                # Skip already done
                if len(signals[
                        (signals["group"] == group) &
                        (signals["region"] == region_name) &
                        (signals["label"] == label)
                ]) > 0:
                    print("Continuing", group, region_name, label)
                    continue

                print(group, region_name, label)
                df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join([group, region_name, label, label]))), index_col=0)
                df = df.mean(0).reset_index(name="value").rename(columns={"index": "distance"})
                df["group"] = group
                df["region"] = region_name
                df["label"] = label
                signals = signals.append(df, ignore_index=True)
    signals.to_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"), index=False)

    signals = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"))

    signals = signals[(signals["distance"] > -400) & (signals["distance"] < 400)]

    # plot together
    region_order = sorted(signals["region"].unique(), reverse=True)
    group_order = sorted(signals["group"].unique(), reverse=True)
    label_order = sorted(signals["label"].unique(), reverse=True)

    # raw
    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "value")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.raw_mean_coverage.svg"), bbox_inches="tight")

    # normalized
    signals["norm_values"] = signals.groupby(["region", "label", "group"])["value"].apply(lambda x: (x - x.mean()) / x.std())

    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.norm_mean_coverage.svg"), bbox_inches="tight")

    # normalized smoothed
    signals["norm_smooth_values"] = signals.groupby(["region", "label", "group"])["value"].apply(lambda x: pd.rolling_window(((x - x.mean()) / x.std()), 10))

    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_smooth_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.norm_mean_coverage.smooth.svg"), bbox_inches="tight")

    #
    # specific regions/samples
    specific_signals = signals[
        (signals["group"].str.contains("ARID1|BCL7B|SMARCA4|C631")) &
        (~signals["group"].str.contains("OV90|GFP")) &

        # (signals["region"].str.contains("diff_sites")) &

        (signals["label"] == "nucleoatac")
    ]

    region_order = sorted(specific_signals["region"].unique(), reverse=True)
    group_order = sorted(specific_signals["group"].unique(), reverse=True)
    label_order = sorted(specific_signals["label"].unique(), reverse=True)

    g = sns.FacetGrid(specific_signals, hue="group", col="region", col_wrap=4, hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "value")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.extended.svg"), bbox_inches="tight")

    # zoom in center
    g = sns.FacetGrid(specific_signals[
        (specific_signals["distance"] < 250) &
        (specific_signals["distance"] > -250)],
        hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "value")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.extended.zoom.svg"), bbox_inches="tight")

    # normalized (centered)
    specific_signals["norm_values"] = specific_signals.groupby(["region", "label", "group"])["value"].apply(lambda x: (x - x.mean()) / x.std())
    g = sns.FacetGrid(specific_signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.norm.svg"), bbox_inches="tight")

    # zoom in center
    g = sns.FacetGrid(specific_signals[
        (specific_signals["distance"] < 250) &
        (specific_signals["distance"] > -250)], hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.norm.zoom.svg"), bbox_inches="tight")

    #

    # Violinplots of nucleosome occupancies

    # Heatmap of nucleosome occupancies
    # Collect signals
    sel_groups = [x for x in groups if "OV90" not in x and "GFP" not in x and ("ARID1" in x or "SMARCA" in x or "WT" in x)]
    regions = pickle.load(open(regions_pickle, "rb"))
    sel_regions = {k: v for k, v in regions.items() if "diff" in k}

    # get parameters based on WT accessibility
    region_order = dict()
    region_vmax = dict()
    region_norm_vmax = dict()
    output_dir = os.path.join(self.results_dir, "nucleoatac", "ATAC-seq_HAP1_WT_C631")
    for region_name in sel_regions.keys():
        df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join(["ATAC-seq_HAP1_WT_C631", region_name, "nucleoatac", "nucleoatac"]))), index_col=0)
        # vmax
        region_vmax[region_name] = np.percentile(df, 95)
        region_norm_vmax[region_name] = np.percentile((df - df.mean(0)) / df.std(0), 95)
        # region order
        region_order[region_name] = df.sum(axis=1).sort_values().index  # sorted by mean
        # region_order[region_name] = g.dendrogram_row.dendrogram

    # plot all
    fig, axis = plt.subplots(len(sel_regions), len(sel_groups), figsize=(len(sel_groups) * 4, len(sel_regions) * 4))
    fig2, axis2 = plt.subplots(len(sel_regions), len(sel_groups), figsize=(len(sel_groups) * 4, len(sel_regions) * 4))
    for j, group in enumerate(sorted(sel_groups, reverse=True)):
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)
        signal_files = [
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz"))
        ]
        for i, (region_name, bed_file) in enumerate(sel_regions.items()):
            for label, signal_file in signal_files:
                print(group, region_name, label)
                df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join([group, region_name, label, label]))), index_col=0)

                d = df.ix[region_order[region_name]]
                axis[i, j].imshow(
                    d,
                    norm=None, cmap="inferno", vmax=region_vmax[region_name], extent=[-500, 500, 0, 10], aspect="auto")  # aspect=100
                d_norm = (d - d.mean(0)) / d.std(0)
                axis2[i, j].imshow(
                    d_norm,
                    norm=None, cmap="inferno", vmax=region_norm_vmax[region_name], extent=[-500, 500, 0, 10], aspect="auto")  # aspect=100
                for ax in [axis, axis2]:
                    ax[i, j].set_title(group)
                    ax[i, j].set_xlabel("distance")
                    ax[i, j].set_ylabel(region_name)
    sns.despine(fig, top=True, right=True, left=True, bottom=True)
    sns.despine(fig2, top=True, right=True, left=True, bottom=True)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.svg"), bbox_inches="tight")
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.centered.png"), bbox_inches="tight", dpi=300)
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.centered.svg"), bbox_inches="tight")


def phasograms(self, samples, max_dist=10000, rolling_window=50, plotting_window=(0, 500)):
    df = self.prj.sheet.df[self.prj.sheet.df["library"] == "ATAC-seq"]
    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    def difference_matrix(a):
        x = np.reshape(a, (len(a), 1))
        return x - x.transpose()

    distances = dict()

    for group in groups:
        print(group)
        # Get dyad calls from nucleoatac
        df = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", group, group + ".nucpos.bed.gz"), sep="\t", header=None)

        # separate by chromosome (groupby?)
        # count pairwise distance
        dists = list()
        for chrom in df[0].unique():
            d = abs(difference_matrix(df[df[0] == chrom][1]))
            dd = d[(d < max_dist) & (d != 0)]
            dists += dd.tolist()
        distances[group] = dists

    pickle.dump(distances, open(os.path.join(self.results_dir, "nucleoatac", "phasogram.distances.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    distances = pickle.load(open(os.path.join(self.results_dir, "nucleoatac", "phasogram.distances.pickle"), "rb"))

    # Plot distances between dyads
    from scipy.ndimage.filters import gaussian_filter1d
    n_rows = n_cols = int(np.ceil(np.sqrt(len(groups))))
    n_rows -= 1
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    fig2, axis2 = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    axis = axis.flatten()
    axis2 = axis2.flatten()
    for i, group in enumerate(groups):
        # Count frequency of dyad distances
        x = pd.Series(distances[group])
        y = x.value_counts().sort_index()
        y = y.ix[range(plotting_window[0], plotting_window[1])]
        y /= y.sum()

        # Find peaks
        y2 = pd.Series(gaussian_filter1d(y, 5), index=y.index)
        peak_indices = detect_peaks(y2.values, mpd=73.5)[:3]
        print(group, y2.iloc[peak_indices].index)

        # Plot distribution and peaks
        axis[i].plot(y.index, y, color="black", alpha=0.6, linewidth=0.5)
        axis[i].plot(y2.index, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis[i].scatter(y2.iloc[peak_indices].index, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis[i].axvline(peak, color="black", linestyle="--")
        axis[i].set_title(group)

        # Transform into distances between nucleosomes
        # Plot distribution and peaks
        axis2[i].plot(y.index - 147, y, color="black", alpha=0.6, linewidth=0.5)
        axis2[i].plot(y2.index - 147, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis2[i].scatter(y2.iloc[peak_indices].index - 147, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis2[i].axvline(peak - 147, color="black", linestyle="--")
        axis2[i].set_title(group)
    sns.despine(fig)
    sns.despine(fig2)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.dyad_distances.peaks.svg"), bbox_inches="tight")
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.nucleosome_distances.peaks.svg"), bbox_inches="tight")

    # Get NFR per knockout
    lengths = dict()

    for group in groups:
        print(group)
        # Get NFR calls from nucleoatac
        df = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", group, group + ".nfrpos.bed.gz"), sep="\t", header=None)
        # Get lengths
        lengths[group] = (df[2] - df[1]).tolist()

    pickle.dump(lengths, open(os.path.join(self.results_dir, "nucleoatac", "nfr.lengths.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    lengths = pickle.load(open(os.path.join(self.results_dir, "nucleoatac", "nfr.lengths.pickle"), "rb"))

    # plot NFR lengths
    from scipy.ndimage.filters import gaussian_filter1d
    n_rows = n_cols = int(np.ceil(np.sqrt(len(groups))))
    n_rows -= 1
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    axis = axis.flatten()
    for i, group in enumerate(groups):
        # Count NFR lengths
        x = pd.Series(lengths[group])
        y = x.value_counts().sort_index()
        y = y.ix[range(plotting_window[0], 300)]
        y /= y.sum()

        # Find peaks
        y2 = pd.Series(gaussian_filter1d(y, 5), index=y.index)
        peak_indices = [detect_peaks(y2.values, mpd=73.5)[0]]
        print(group, y2.iloc[peak_indices].index)

        # Plot distribution and peaks
        axis[i].plot(y.index, y, color="black", alpha=0.6, linewidth=0.5)
        axis[i].plot(y2.index, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis[i].scatter(y2.iloc[peak_indices].index, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis[i].axvline(peak, color="black", linestyle="--")
        axis[i].set_title(group)

    sns.despine(fig)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.nfr_lengths.peaks.svg"), bbox_inches="tight")


def run_coverage_job(bed_file, bam_file, coverage_type, name, output_dir, window_size=1001):
    from pypiper import NGSTk
    tk = NGSTk()
    job_file = os.path.join(output_dir, "%s.run_enrichment.sh" % name)
    log_file = os.path.join(output_dir, "%s.run_enrichment.log" % name)

    cmd = """#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --nodes=1

#SBATCH --job-name=baf-kubicek-run_enrichment_{}
#SBATCH --output={}

#SBATCH --mail-type=end
#SBATCH --mail-user=

# Start running the job
hostname
date

cd /home/arendeiro/baf-kubicek/

python /home/arendeiro/jobs/run_profiles.py \
--bed-file {} \
--bam-file {} \
--coverage-type {} \
--window-size {} \
--name {} \
--output-dir {}

date
""".format(
        name,
        log_file,
        bed_file,
        bam_file,
        coverage_type,
        window_size,
        name,
        output_dir
    )

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(cmd)

    tk.slurm_submit_job(job_file)


def run_vplot_job(bed_file, bam_file, name, output_dir):
    from pypiper import NGSTk
    tk = NGSTk()
    job_file = os.path.join(output_dir, "%s.run_enrichment.sh" % name)
    log_file = os.path.join(output_dir, "%s.run_enrichment.log" % name)

    cmd = """#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=1-12:00:00

#SBATCH --cpus-per-task=8
#SBATCH --mem=24000
#SBATCH --nodes=1

#SBATCH --job-name=baf-kubicek-run_enrichment_{name}
#SBATCH --output={log}

#SBATCH --mail-type=end
#SBATCH --mail-user=

# Start running the job
hostname
date

\t\t# Remove everything to do with your python and env.  Even reset your home dir
\t\tunset PYTHONPATH
\t\tunset PYTHON_HOME
\t\tmodule purge
\t\tmodule load python/2.7.6
\t\tmodule load slurm
\t\tmodule load gcc/4.8.2

\t\tENV_DIR=/scratch/users/arendeiro/nucleoenv
\t\texport HOME=$ENV_DIR/home

\t\t# Activate your virtual env
\t\texport VIRTUALENVWRAPPER_PYTHON=/cm/shared/apps/python/2.7.6/bin/python
\t\tsource $ENV_DIR/bin/activate

\t\t# Prepare to install new python packages
\t\texport PATH=$ENV_DIR/install/bin:$PATH
\t\texport PYTHONPATH=$ENV_DIR/install/lib/python2.7/site-packages


cd /home/arendeiro/baf-kubicek/

pyatac vplot \
--bed {peaks} \
--bam {bam} \
--out {output_prefix}.vplot \
--cores 8 \
--lower 30 \
--upper 1000 \
--flank 500 \
--scale \
--plot_extra

pyatac sizes \
--bam {bam} \
--bed {peaks} \
--out {output_prefix}.sizes  \
--lower 30 \
--upper 1000

pyatac bias \
--fasta ~/resources/genomes/hg19/hg19.fa \
--bed {peaks} \
--out {output_prefix}.bias \
--cores 8

pyatac bias_vplot \
--bed {peaks} \
--bg {output_prefix}.bias.Scores.bedgraph.gz \
--fasta ~/resources/genomes/hg19/hg19.fa \
--sizes {output_prefix}.sizes.fragmentsizes.txt \
--out {output_prefix}.bias_vplot \
--cores 8 \
--lower 30 \
--upper 1000 \
--flank 500 \
--scale \
--plot_extra

date
""".format(
        name=name,
        log=log_file,
        peaks=bed_file,
        bam=bam_file,
        output_prefix=os.path.join(output_dir, name))

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(cmd)

    tk.slurm_submit_job(job_file)


def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`

    The function can handle NaN's

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan - 1, indnan + 1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size - 1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    return ind


def footprint(samples, attributes=["cell_line", "knockout"], output_dir="footprinting"):

    def chunks(l, n):
        n = max(1, n)
        return (l[i:i + n] for i in range(0, len(l), n))

    df = pd.DataFrame([s.as_series() for s in samples])
    df.to_csv(os.path.join("metadata", "to_footprint.annotation.csv"), index=False)

    # Prepare motifs
    piq_prepare_motifs()

    # Prepare BAM files
    # for each knockout
    for attrs, index in df.groupby(attributes).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        bams = [s.filtered for s in samples if s.name in df.loc[index, "sample_name"].tolist()]

        piq_prepare_bams(
            bams,
            name,
            output_dir=output_dir,
            piq_source_dir="/home/arendeiro/workspace/piq-single/")
    # for all samples
    piq_prepare_bams([s.filtered for s in samples], "all_samples")

    # Footprint
    for motif_numbers in range(1, 367):
        for attrs, index in df.groupby(attributes).groups.items():
            name = "_".join([a for a in attrs if not pd.isnull(a)])
            print(name)
            footprint(name, motif_numbers=[motif_numbers])

    footprint("all_samples", motif_numbers=range(1, 367))

    # Collect footprints, assemble TF-gene network
    for attrs, index in df.groupby(attributes).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        piq_to_network(
            group_name=name, motif_numbers=range(1, 367),
            peak_universe_file=os.path.join("results", "baf-kubicek_peak_set.slop_b500.bed"))
    piq_to_network(
        group_name="all_samples", motif_numbers=range(1, 367),
        peak_universe_file=os.path.join("results", "baf-kubicek_peak_set.slop_b500.bed"))

    # OR in parallel:
    # launch
    batch_size = 10
    for attrs, index in df.groupby(attributes).groups.items():
        for chunk in [",".join([str(j) for j in i]) for i in chunks(range(1, 367), batch_size)]:
            name = "_".join([a for a in attrs if not pd.isnull(a)])
            cmd = """\
sbatch -c 1 --mem 20000 -p shortq \
-J PIQ.gather_interactions.{g}.{m} \
-o ~/projects/baf-kubicek/footprinting/footprint_calls/jobs/PIQ.gather_interactions.{g}.{m}.log \
-D ~/projects/baf-kubicek \
--wrap "python ~/projects/baf-kubicek/gather_tfnetwork.job.py -g {g} -m {m}" """.format(g=name, m=chunk)
            os.system(cmd)

#             cmd = """\
# sbatch -c 1 --mem 20000 -p shortq \
# -J PIQ.gather_interactions.{g}.{m} \
# -o ~/projects/baf-kubicek/footprinting/footprint_calls/jobs/PIQ.gather_interactions.{g}.{m}.log \
# -D ~/projects/baf-kubicek \
# --wrap "python ~/projects/baf-kubicek/gather_tfnetwork.job.py -g {g} -m {m}" """.format(g="all_samples", m=chunk)
#             os.system(cmd)

    # OR in parallel:
    # gather
    batch_size = 10
    for attrs, index in df.groupby(attributes).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        cmd = """\
sbatch -c 1 --mem 20000 -p shortq \
-J PIQ.gather_interactions.reduce.{g}.all \
-o ~/projects/baf-kubicek/footprinting/footprint_calls/jobs/PIQ.gather_interactions.reduce.{g}.all.log \
-D ~/projects/baf-kubicek \
--wrap "python ~/projects/baf-kubicek/gather_tfnetwork.job.py -g {g} " """.format(g=name)
        os.system(cmd)

#      cmd = """\
# sbatch -c 1 --mem 20000 -p shortq \
# -J PIQ.gather_interactions.{g}.all \
# -o ~/projects/baf-kubicek/footprinting/footprint_calls/jobs/PIQ.gather_interactions.{g}.all.log \
# -D ~/projects/baf-kubicek \
# --wrap "python ~/projects/baf-kubicek/gather_tfnetwork.job.py -g {g} -m {m}" """.format(g="all_samples", m=chunk)
#      os.system(cmd)
    groups = map(
        lambda x: "_".join([a for a in x[0] if not pd.isnull(a)]),
        df.groupby(attributes).groups.items())
    for group in groups:
        differential_interactions(group, "HAP1_WT")

    # gather all
    all_stats = pd.DataFrame()
    for group in groups:
        if group == "HAP1_WT":
            continue
        comparison_name = "{}-{}".format(group, "HAP1_WT")
        tf_stats = pd.read_csv(os.path.join(diff_dir, "tf_differential_binding.{}.csv".format(comparison_name)))
        tf_stats.loc[:, "knockout"] = group.replace("HAP1_", "")

        all_stats = all_stats.append(tf_stats, ignore_index=True)

    all_stats.loc[:, "sign_log_q_value"] = -np.log10(all_stats.loc[:, "q_value"] + 1e-80) * (all_stats.loc[:, "log_fold_change"] > 0).astype(int).replace(0, -1)

    tf_pivot = pd.pivot_table(all_stats, index="tf_name", columns="knockout", values="sign_log_q_value")

    # Knockout correlation of changes in binding
    g = sns.clustermap(tf_pivot.corr(), cbar_kws={"label": "Pearson correlation of KOs"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("footprinting", "signed_qvalue.KO_change.correlation.svg"), bbox_inches="tight")

    # TF correlation of changes in binding
    c = tf_pivot.dropna().T.dropna().corr().drop_duplicates().T.drop_duplicates()
    c = c.loc[~c.isnull().all(axis=0), ~c.T.isnull().all(axis=1)]

    g = sns.clustermap(c, cbar_kws={"label": "Pearson correlation of TFs"}, rasterized=True, figsize=(20, 20))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("footprinting", "signed_qvalue.TF_change.correlation.svg"), bbox_inches="tight", dpi=300)

    # Full matrix of changes
    g = sns.clustermap(tf_pivot.dropna().T.dropna(), cbar_kws={"label": "Sign * -log10(Q value) of change"}, rasterized=True, figsize=(22, 8))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("footprinting", "signed_qvalue.KO-TF_change.svg"), bbox_inches="tight", dpi=300)


    tf_pivot = pd.pivot_table(all_stats, index="tf_name", columns="knockout", values="log_fold_change")

    # Knockout correlation of changes in binding
    g = sns.clustermap(tf_pivot.corr(), cbar_kws={"label": "Pearson correlation of KOs"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("footprinting", "log_fold_change.KO_change.correlation.svg"), bbox_inches="tight")

    # TF correlation of changes in binding
    c = tf_pivot.dropna().T.dropna().corr().drop_duplicates().T.drop_duplicates()
    c = c.loc[~c.isnull().all(axis=0), ~c.T.isnull().all(axis=1)]

    g = sns.clustermap(c, cbar_kws={"label": "Pearson correlation of TFs"}, rasterized=True, figsize=(20, 20))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.savefig(os.path.join("footprinting", "log_fold_change.TF_change.correlation.svg"), bbox_inches="tight", dpi=300)

    # Full matrix of changes
    g = sns.clustermap(tf_pivot.dropna().T.dropna(), cbar_kws={"label": "Sign * -log10(Q value) of change"}, rasterized=True, figsize=(22, 8))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("footprinting", "log_fold_change.KO-TF_change.svg"), bbox_inches="tight", dpi=300)


def chromatin_protein_expression():
    """
    Explore the expression of chromatin-related genes across the various knockouts.
    """
        
    trait = "knockout"
    output_suffix = "deseq_expression_knockout"
    results_dir = "results"
    output_dir = os.path.join(results_dir, output_suffix)
    df = pd.read_csv(os.path.join(output_dir, output_suffix) + ".%s.annotated.csv" % trait, index_col=0)

    chrom_list = pd.read_csv(os.path.join("metadata", "Bocklab_chromatin_genes.csv"))

    comps = [x for x in df['comparison'].unique().tolist() + ['WT'] if x != "SMARCC2"]
    # Pure expression
    exp = np.log2(1 + df.loc[:, comps].drop_duplicates())
    e = exp.loc[chrom_list['HGNC_symbol'], :].dropna()
    g = sns.clustermap(e.T, z_score=1, xticklabels=False, figsize=(e.shape[0] * 0.01, e.shape[1] * 0.12), rasterized=True)
    g.ax_heatmap.set_xlabel("Chromatin genes", ha="center")
    g.ax_heatmap.set_ylabel("Knockouts", ha="center")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", "baf-kubicek.chromatin_gene_expression.svg"), dpi=300, bbox_inches="tight")

    e.to_csv(os.path.join("results", "baf-kubicek.chromatin_gene_expression.csv"), index=True)


    # Only significant in any comparison
    all_diff = df[(df['padj'] < 0.01)].index.unique().tolist()

    diff = chrom_list['HGNC_symbol'][chrom_list['HGNC_symbol'].isin(all_diff)]
    e = exp.loc[diff, :].dropna()

    g = sns.clustermap(e.T, z_score=1, yticklabels=True, figsize=(e.shape[0] * 0.12, e.shape[1] * 0.12), rasterized=True)
    g.ax_heatmap.set_xlabel("Chromatin genes (all differential)", ha="center")
    g.ax_heatmap.set_ylabel("Knockouts", ha="center")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", "baf-kubicek.chromatin_gene_expression.differential_only.svg"), dpi=300, bbox_inches="tight")


    # Only significant in any comparison
    all_diff = df[(df['padj'] < 0.01) & (abs(df['log2FoldChange']) > 1)].index.unique().tolist()

    diff = chrom_list['HGNC_symbol'][chrom_list['HGNC_symbol'].isin(all_diff)]
    e = exp.loc[diff, :].dropna()

    g = sns.clustermap(e.T, z_score=1, yticklabels=True, figsize=(e.shape[0] * 0.12, e.shape[1] * 0.12), rasterized=True)
    g.ax_heatmap.set_xlabel("Chromatin genes (all differential with FC > 1)", ha="center")
    g.ax_heatmap.set_ylabel("Knockouts", ha="center")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="left", fontsize="xx-small")
    g.savefig(os.path.join("results", "baf-kubicek.chromatin_gene_expression.differential_fc1.svg"), dpi=300, bbox_inches="tight")


def discordance_analysis(atac_analysis, chipseq_analysis, rnaseq_analysis):
    """
    """
    def index_to_bed_file(index, bed_file):
        d = pd.DataFrame(
            zip(
                map(lambda x: x[0], index.str.split(":")),
                map(lambda x: x[1].split("-")[0], index.str.split(":")),
                map(lambda x: x[1].split("-")[1], index.str.split(":"))))
        d.to_csv(bed_file, index=False, header=False, sep="\t")

    
    # Get ChIP-seq data
    chipseq_analysis.sites = pybedtools.BedTool("results/baf_complex.chipseq.peaks_peak_set.bed")
    chipseq_analysis.support = pd.read_csv("results/baf_complex.chipseq.peaks_peaks.support.csv", index_col=0, header=range(2))

    # filter to get sites bound by ARID1A/B
    comparison_table = pd.read_csv(
        os.path.join("metadata", "comparison_table.csv"))
    c = comparison_table[
        (comparison_table['comparison_type'] == 'peaks') &
        (comparison_table['comparison_name'].str.contains("ARID|SMARC|PBRM"))]
    chipseq_analysis.calculate_peak_support(comparison_table=c)

    support = chipseq_analysis.support.drop('support', axis=1)
    support_sum = (
        support
        .loc[:, support.columns.get_level_values("comparison").str.contains("ARID1")]
        .astype(bool).astype(int).sum(axis=1))

    # get sites with at least two peaks called
    bound_sites = support_sum[support_sum >= 2].index
    # make bed file
    index_to_bed_file(
        bound_sites,
        os.path.join("results", "ARID1A-ARID1B_bound_sites.bed"))

    # get only relevant comparisons (KOs only)
    diff = atac_analysis.differential_results
    diff = diff.loc[~diff['comparison_name'].str.contains("sh|dBet6")]
    diff = diff.loc[(diff['padj'] < 0.01) & (diff['log2FoldChange'].abs() > 1)]

    # intersect differential ATAC-seq peaks with bound sites
    index_to_bed_file(
        diff.index.unique(),
        os.path.join("results", "differential_sites.bed"))

    diff_bound = (
        pybedtools.BedTool(os.path.join("results", "differential_sites.bed"))
        .intersect(b=os.path.join("results", "ARID1A-ARID1B_bound_sites.bed"), wa=True))
    diff_bound = diff_bound.to_dataframe()
    diff_bound = diff_bound['chrom'] + ":" + diff_bound['start'].astype(str) + '-' + diff_bound['end'].astype(str)

    diff_unbound = (
        pybedtools.BedTool(os.path.join("results", "differential_sites.bed"))
        .intersect(b=os.path.join("results", "ARID1A-ARID1B_bound_sites.bed"), wa=True, v=True))
    diff_unbound = diff_unbound.to_dataframe()
    diff_unbound = diff_unbound['chrom'] + ":" + diff_unbound['start'].astype(str) + '-' + diff_unbound['end'].astype(str)

    # compare proportion of discordant changes vs concordant in bound sites vs all differntial for these subunits
    differential_overlap(
        diff.loc[diff_bound],
        total=atac_analysis.accessibility.shape[0],
        output_dir="{}/differential_analysis_ATAC-seq".format(atac_analysis.results_dir),
        data_type="ATAC-seq", output_prefix='differential_analysis.only_bound_sites')
    differential_overlap(
        diff.loc[diff_unbound],
        total=len(bound_sites),
        output_dir="{}/differential_analysis_ATAC-seq".format(atac_analysis.results_dir),
        data_type="ATAC-seq", output_prefix='differential_analysis.only_unbound_sites')

    differential_overlap(
        diff,
        total=len(bound_sites),
        output_dir="{}/differential_analysis_ATAC-seq".format(atac_analysis.results_dir),
        data_type="ATAC-seq", output_prefix='differential_analysis.all_sites_sites')

    # plot gene expression fold-change for genes with discordant vs concordant changes between the two subunits

    rnaseq_analysis.differential_results = pd.read_csv(
        os.path.join(
            "results", "differential_analysis_RNA-seq",
            "differential_analysis" + ".deseq_result.all_comparisons.csv"), index_col=0)
    rna = rnaseq_analysis.differential_results
    rna = rna.loc[rna['comparison_name'].isin(['ARID1A', "ARID1B"])]
    # annotate genes with regions
    g = atac_analysis.gene_annotation
    rna = rna.join(
        g['gene_name'].str.split(",")
        .apply(pd.Series).stack()
        .reset_index(drop=True, level=1)
        .reset_index().set_index(0)
    ).rename(columns={"index": "region"})


    d = diff.loc[diff['comparison_name'].isin(['ARID1A', "ARID1B"])]
    if "direction" not in d.columns:
        d["direction"] = d["log2FoldChange"].apply(lambda x: "up" if x > 0 else "down")
    intersections = pd.DataFrame(columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
    perms = list(itertools.permutations(d.groupby(['comparison_name', 'direction']).groups.items(), 2))
    res = pd.DataFrame()
    for ((k1, dir1), i1), ((k2, dir2), i2) in tqdm(perms, total=len(perms)):
        i1 = set(i1)
        i2 = set(i2)
        inter = i1.intersection(i2)
        union = i1.union(i2)

        k1_inter_mean = rna.loc[
            (rna['comparison_name'] == k1) &
            (rna['region'].isin(inter)), "log2FoldChange"].mean()
        k1_union_mean = rna.loc[
            (rna['comparison_name'] == k1) &
            (rna['region'].isin(union)), "log2FoldChange"].mean()
        k2_inter_mean = rna.loc[
            (rna['comparison_name'] == k2) &
            (rna['region'].isin(inter)), "log2FoldChange"].mean()
        k2_union_mean = rna.loc[
            (rna['comparison_name'] == k2) &
            (rna['region'].isin(union)), "log2FoldChange"].mean()
        res = res.append(
            pd.Series([k1, dir1, k2, dir2,
                       k1_inter_mean, k1_union_mean,
                       k2_inter_mean, k2_union_mean],
                      index=["k1", "dir1", "k2", "dir2",
                             "k1_inter_mean", "k1_union_mean",
                             "k2_inter_mean", "k2_union_mean"]),
            ignore_index=True)
    
    res['label'] = res['k1'] + " " + res['dir1'] + "\n" + res['k2'] + " " + res['dir2']
    res2 = pd.melt(res.dropna()[['label', 'k1_inter_mean', 'k2_inter_mean']], id_vars=['label'])

    fig, axis = plt.subplots(1, figsize=(1 * 4, 1 * 4))
    sns.barplot(data=res2, x="label", y="value", hue="variable", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join("results", "differential_analysis.differential_overlap.bound_sites.disagreement.expression.barplot.svg"), bbox_inches="tight")


    # Investigate discordant changes within the same gene for one specific subunit
    # plot expression fold change for genes depending on the number of discordant/concordant reg.elements for each subunit


def chips_across_knockouts():

    # ChIP-seq on ARID1A across all knockouts
    chipseq_samples = [s for s in prj.samples if s.library in ["ChIP-seq", "ChIPmentation"]]
    # chipseq_samples = [s for s in chipseq_samples if os.path.exists(s.filtered)]
    chipseq_analysis = ChIPSeqAnalysis(name="baf_complex.chipseq.arid1a", prj=prj, samples=chipseq_samples)
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    c = comparison_table[
        (comparison_table['comparison_type'] == 'peaks') &
        (comparison_table['comparison_name'].str.contains("ARID|SMARC|PBRM")) &
        (comparison_table['toggle'] == 1)]
    c['comparison_genome'] = 'hg19'
    chipseq_analysis.call_peaks_from_comparisons(comparison_table=c)
    chipseq_analysis.summarize_peaks_from_comparisons(comparison_table=c, output_dir="{results_dir}/chipseq_peaks", permissive=False)
    chipseq_analysis = get_consensus_sites(comparison_table=c, region_type="peaks", blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed")
    chipseq_analysis.calculate_peak_support(comparison_table=c)
    chipseq_analysis.measure_coverage()
    chipseq_analysis.normalize()
    chipseq_analysis.annotate(quant_matrix="coverage_qnorm")
    chipseq_analysis.annotate_with_sample_metadata(attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    chipseq_analysis.to_pickle()

    # Unsupervised analysis
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=None,
        attributes_to_plot=['ip', 'knockout', 'replicate'], plot_prefix="chipseq_baf_peaks.arid1a",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")
    # without histone marks
    unsupervised_analysis(
        chipseq_analysis, data_type="ATAC-seq", quant_matrix=None, samples=[s for s in chipseq_analysis.samples if ("ARID" in s.name) or ("SMARC" in s.name)],
        attributes_to_plot=['ip', 'knockout', 'replicate'], plot_prefix="chipseq_baf_peaks.arid1a.arid1a_only",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")

    # Supervised analysis
    data_type = "ChIP-seq"
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['data_type'] == data_type) &
        (comparison_table['comparison_type'] == 'differential')]
    # chipseq_analysis.differential_results = differential_analysis(
    #     chipseq_analysis,
    #     comparison_table,
    #     data_type="ATAC-seq",
    #     samples=[s for s in chipseq_analysis.samples if s.name in comparison_table['sample_name'].tolist()],
    #     output_dir="{}/differential_analysis_{}".format(chipseq_analysis.results_dir, data_type),
    #     covariates=None,
    #     alpha=0.05,
    #     overwrite=True,
    #     distributed=True)
    # chipseq_analysis.differential_results = chipseq_analysis.differential_results.set_index("index")
    from ngs_toolkit.general import differential_from_bivariate_fit

    chipseq_analysis.differential_results = differential_from_bivariate_fit(
        comparison_table,
        chipseq_analysis.accessibility,
        output_dir=os.path.join("results", "differential_analysis_ChIP-seq"),
        output_prefix="differential_analysis_ChIP-seq.bivariate_fit",
        n_bins=250, multiple_correction_method="fdr_bh",
        plot=True, palette="colorblind", make_values_positive=True)
    chipseq_analysis.to_pickle()

    plot_differential(
        chipseq_analysis,
        chipseq_analysis.differential_results.rename(columns={"global_mean": "baseMean", "pval": "pvalue"}), # chipseq_analysis.differential_results[~chipseq_analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")], 
        # matrix=getattr(chipseq_analysis, "accessibility"),
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(chipseq_analysis.results_dir, data_type),
        output_prefix="differential_analysis_ChIP-seq.bivariate_fit",
        data_type="ATAC-seq",
        alpha=0.005,
        corrected_p_value=False,
        fold_change=None,
        rasterized=True,
        robust=True,
        group_wise_colours=False,
        group_variables=plotting_attributes[:3])



def interaction_new_data():

    # RNA-seq on subunit combination perturbation

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["HAP1"])]#  & (s.flowcell in ["BSF_0426_HTHTJBBXX"])]
    # rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    kd_rnaseq_analysis = RNASeqAnalysis(name="baf-complex.hap1.rnaseq.KD", prj=prj, samples=rnaseq_samples, results_dir="results")
    kd_rnaseq_analysis = main_analysis_pipeline(kd_rnaseq_analysis, data_type="RNA-seq", cell_type="HAP1")
    # kd_rnaseq_analysis.samples = kd_rnaseq_analysis.samples + rnaseq_analysis.samples
    # kd_rnaseq_analysis.expression_annotated_all = kd_rnaseq_analysis.expression_annotated.join(rnaseq_analysis.expression_annotated)

    # Unsupervised analysis
    red_samples = [s for s in kd_rnaseq_analysis.samples if (
        ("dBet" not in s.name) and
        ("BRD4" not in s.name) and
        ("parental" not in s.name) and
        (s.name in kd_rnaseq_analysis.expression_annotated_all.columns))]
    unsupervised_analysis(
        kd_rnaseq_analysis,
        quant_matrix="expression_annotated",
        samples=red_samples,
        attributes_to_plot=plotting_attributes,
        plot_prefix="all_{}".format(feature_name),
        plot_max_attr=20,
        plot_max_pcs=4,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False,
        test_pc_association=False,
        output_dir="{results_dir}/unsupervised.20180131")

    # fix batch effect
    kd_rnaseq_analysis.matrix_batch_fix = fix_batch_effect(
        getattr(kd_rnaseq_analysis, "expression_annotated"), kd_rnaseq_analysis.samples,
        batch_variable="batch", standardize=True, intermediate_save=True)
    file = os.path.join(kd_rnaseq_analysis.results_dir, kd_rnaseq_analysis.name + ".{}.annotated_metadata.batch_fix.csv".format("expression_annotated"))
    kd_rnaseq_analysis.matrix_batch_fix.to_csv(file)

    # rescale variables
    m = kd_rnaseq_analysis.expression_annotated
    mean = m.mean(axis=1)
    std = m.std(axis=1)
    kd_rnaseq_analysis.matrix_batch_fix_rescalled = pd.DataFrame(
        np.multiply(np.add(kd_rnaseq_analysis.matrix_batch_fix.values.T, mean.values), std.values).T,
        index=m.index, columns=m.columns)

    unsupervised_analysis(
        kd_rnaseq_analysis,
        quant_matrix="matrix_batch_fix_rescalled",
        samples=red_samples,
        attributes_to_plot=plotting_attributes,
        plot_prefix="all_{}-matrix_batch_fix_rescalled".format(feature_name),
        plot_max_attr=20,
        plot_max_pcs=4,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False,
        test_pc_association=False,
        output_dir="{results_dir}/unsupervised.20180131")


    # Find the interaction
    rnaseq_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format("results", "RNA-seq"), "differential_analysis.enrichr.csv"))

    q = rnaseq_enrichment_table.loc[
        (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
        (rnaseq_enrichment_table['gene_set_library'].isin(["NCI-Nature_2016", "WikiPathways_2016"])) &
        (rnaseq_enrichment_table['direction'] == 'down') &
        # (rnaseq_enrichment_table['p_value'] < 0.05) &
        (
            rnaseq_enrichment_table['description'].str.contains("E2F") |
            rnaseq_enrichment_table['description'].str.contains("cell cycle", case=False)), :]

    genes = q.loc[:, "genes"]

    cc_genes = list(set(genes.str.replace("[", " ").str.replace(
        ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # clustermap
    g = sns.clustermap(kd_rnaseq_analysis.matrix_batch_fix.loc[cc_genes, :].dropna(), z_score=0, rasterized=True, xticklabels=True, cmap="RdBu_r")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature&WikiPathways.clustermap.svg"), bbox_inches="tight", dpi=300)


    g = sns.clustermap(kd_rnaseq_analysis.matrix_batch_fix.loc[cc_genes, :].dropna().T.groupby(["knockout", "treatment"]).mean().T, z_score=0, rasterized=True, xticklabels=True, cmap="RdBu_r", square=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature&WikiPathways.group.clustermap.svg"), bbox_inches="tight", dpi=300)



    # cc_genes = ["CCND3", "RBL1", "CCND2", "CDK2", "CDC25A"]
    q = rnaseq_enrichment_table.loc[
        (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
        # (rnaseq_enrichment_table['gene_set_library'] == 'NCI-Nature_2016') &
        (rnaseq_enrichment_table['direction'] == 'down') &
        (rnaseq_enrichment_table['p_value'] < 0.05), :]

    genes = q.loc[:, "genes"]

    cc_genes = list(set(genes.str.replace("[", " ").str.replace(
        ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # clustermap
    g = sns.clustermap(kd_rnaseq_analysis.expression_annotated.loc[cc_genes, :].dropna(), z_score=0, rasterized=True, xticklabels=True, yticklabels=False, cmap="RdBu_r")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=4)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.svg"), bbox_inches="tight", dpi=300)

    # investigate genes in second cluster (ARID2-specific)
    clusters = scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=2, criterion='maxclust')
    # plot again just to confirm clusters
    g2 = sns.clustermap(
        kd_rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True,
        row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage, row_colors=plt.get_cmap("Paired")(clusters))
    g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    g2.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.clusters_labeled.svg"), bbox_inches="tight", dpi=300)

    pbaf_genes = pd.Series(g.data.index).iloc[clusters == pd.Series(clusters).value_counts().argmin()].sort_values()

    g3 = sns.clustermap(kd_rnaseq_analysis.expression.loc[pbaf_genes, :], z_score=0, rasterized=True, metric="correlation")
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90)
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_yticklabels(), rotation=0)
    g3.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.pbaf_genes.svg"), bbox_inches="tight", dpi=300)



def interaction_new_data_independent():

    # Supervised analysis
    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["HAP1"]) & (s.flowcell in ["BSF_0426_HTHTJBBXX"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
   # Get gene expression
    rnaseq_analysis = RNASeqAnalysis(name="baf-complex.rnaseq.only_kd_samples", prj=prj, samples=rnaseq_samples)
    rnaseq_analysis.get_gene_expression(
        samples=rnaseq_analysis.samples,
        sample_attributes=["sample_name", "knockout", 'treatment', "replicate", "clone", "batch"])

    quant_matrix = "expression_annotated"
    feature_name = "genes"

    # Unsupervised analysis
    red_samples = [s for s in rnaseq_analysis.samples]
    unsupervised_analysis(
        rnaseq_analysis,
        quant_matrix=quant_matrix,
        samples=red_samples,
        attributes_to_plot=plotting_attributes,
        plot_prefix="all_{}".format(feature_name),
        plot_max_attr=20,
        plot_max_pcs=6,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False,
        output_dir="{results_dir}/unsupervised.20180320")

    data_type = "RNA-seq"
    cell_type = "HAP1"
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['data_type'] == data_type) &
        (comparison_table['cell_type'] == cell_type) &
        (comparison_table['comparison_type'] == 'differential')]
    rnaseq_analysis.differential_results = differential_analysis(
        rnaseq_analysis,
        comparison_table,
        data_type=data_type,
        samples=[s for s in rnaseq_analysis.samples if s.name in comparison_table['sample_name'].tolist()],
        output_dir="{}/differential_analysis_{}.only_kd_samples".format(rnaseq_analysis.results_dir, data_type),
        covariates=None,
        alpha=0.05,
        overwrite=True)
    rnaseq_analysis.differential_results = rnaseq_analysis.differential_results.set_index("index")
    rnaseq_analysis.to_pickle()

    alpha = 0.05
    abs_fold_change = 0

    differential_overlap(
        rnaseq_analysis.differential_results[
            (rnaseq_analysis.differential_results['padj'] < alpha) &
            (rnaseq_analysis.differential_results['log2FoldChange'].abs() >= abs_fold_change)],
        getattr(rnaseq_analysis, quant_matrix).shape[0],
        output_dir="{}/differential_analysis_{}.only_kd_samples".format(rnaseq_analysis.results_dir, data_type),
        data_type=data_type)

    for (alpha, label) in [(0.05, ""), (0.1, ".p0.1")]:
        plot_differential(
            rnaseq_analysis,
            rnaseq_analysis.differential_results, # rnaseq_analysis.differential_results[~rnaseq_analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")], 
            matrix=getattr(rnaseq_analysis, quant_matrix),
            comparison_table=comparison_table,
            output_dir="{}/differential_analysis_{}.only_kd_samples{}".format(rnaseq_analysis.results_dir, data_type, label),
            output_prefix="differential_analysis",
            data_type=data_type,
            alpha=0.1,
            corrected_p_value=True,
            fold_change=None,
            rasterized=True,
            robust=True,
            group_wise_colours=True,
            group_variables=plotting_attributes)

        differential_enrichment(
            rnaseq_analysis,
            rnaseq_analysis.differential_results[
                (rnaseq_analysis.differential_results['padj'] < alpha) &
                (rnaseq_analysis.differential_results['log2FoldChange'].abs() >= abs_fold_change)],
            data_type=data_type,
            output_dir="{}/differential_analysis_{}.only_kd_samples{}".format(rnaseq_analysis.results_dir, data_type, label),
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
            output_dir="{}/differential_analysis_{}.only_kd_samples{}".format(rnaseq_analysis.results_dir, data_type, label),
            permissive=False)

        enrichment_table = pd.read_csv(
            os.path.join("{}/differential_analysis_{}.only_kd_samples{}".format(rnaseq_analysis.results_dir, data_type, label), "differential_analysis.enrichr.csv"))
        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}.only_kd_samples{}".format(rnaseq_analysis.results_dir, data_type, label),
            output_prefix="differential_analysis",
            direction_dependent=True,
            top_n=5)

    # Directed analysis
    # Miracle plot
    # shSMARCA4_over_ARID2KO
    # vs
    # shSMARCA4_over_WT
    d = rnaseq_analysis.differential_results
    double_changes = d[
        (d['comparison_name'] == "shSMARCA4shARID2_over_WT") &
        (d['padj'] < 0.05)]
    a = d.loc[d['comparison_name'] == "shSMARCA4_over_ARID2KO", "log2FoldChange"]
    b = d.loc[d['comparison_name'] == "shSMARCA4_over_WT", "log2FoldChange"]

    fig, axis = plt.subplots(1, figsize=(1 * 4, 1 * 4))
    axis.scatter(a, b, rasterized=True, alpha=0.2, s=2)
    axis.scatter(a.loc[double_changes.index], b.loc[double_changes.index], color="red", rasterized=True, alpha=0.2, s=2)
    done = list()
    for s in a.sort_values().head(25).index:
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().head(25).index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().tail(25).index:
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().tail(25).index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in double_changes.index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
    axis.plot((-3, 3), (-3, 3), color="black", alpha=0.4)
    axis.set_xlabel("shSMARCA4_over_ARID2KO")
    axis.set_ylabel("shSMARCA4_over_WT")
    sns.despine(fig)
    fig.savefig(
        os.path.join(
            "{}/differential_analysis_{}.only_kd_samples".format(rnaseq_analysis.results_dir, data_type),
            "ARID2_SMARCA4_interaction.miracle_plot.svg"), bbox_inches="tight", dpi=300)

    # miracle plot2
    d = rnaseq_analysis.differential_results
    double_changes = d[
        (d['comparison_name'] == "shSMARCA4shARID2_over_WT") &
        (d['padj'] < 0.05)]

    a = d.loc[d['comparison_name'] == "shARID2_over_SMARCA4KO", "log2FoldChange"]
    b = d.loc[d['comparison_name'] == "shARID2_over_WT", "log2FoldChange"]

    fig, axis = plt.subplots(1, figsize=(1 * 4, 1 * 4))
    axis.scatter(a, b, rasterized=True, alpha=0.2, s=2)
    axis.scatter(a.loc[double_changes.index], b.loc[double_changes.index], color="red", rasterized=True, alpha=0.2, s=2)
    done = list()
    for s in a.sort_values().head(25).index:
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().head(25).index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().tail(25).index:
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in a.sort_values().tail(25).index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
        done.append(s)
    for s in double_changes.index:
        if s in done: continue
        axis.text(a.loc[s], b.loc[s], s, color="black", rasterized=True, alpha=0.5, fontsize=3)
    axis.plot((-3, 3), (-3, 3), color="black", alpha=0.4)
    axis.set_xlabel("shARID2_over_SMARCA4KO")
    axis.set_ylabel("shARID2_over_WT")
    sns.despine(fig)
    fig.savefig(
        os.path.join(
            "{}/differential_analysis_{}.only_kd_samples".format(rnaseq_analysis.results_dir, data_type),
            "ARID2_SMARCA4_interaction.miracle_plot2.svg"), bbox_inches="tight", dpi=300)

    # heatmap with top differential from either comparison
    d = rnaseq_analysis.differential_results
    double_changes = d[
        (d['comparison_name'] == "shSMARCA4shARID2_over_WT") &
        (d['padj'] < 0.05)]
    g = d.loc[d['comparison_name'] == "shSMARCA4_over_ARID2KO", "log2FoldChange"].abs().sort_values().tail(25).index.tolist()
    g += d.loc[d['comparison_name'] == "shSMARCA4_over_WT", "log2FoldChange"].abs().sort_values().tail(25).index.tolist()
    d = rnaseq_analysis.differential_results
    double_changes = d[
        (d['comparison_name'] == "shSMARCA4shARID2_over_WT") &
        (d['padj'] < 0.05)]

    g += d.loc[d['comparison_name'] == "shARID2_over_SMARCA4KO", "log2FoldChange"].abs().sort_values().tail(25).index.tolist()
    g += d.loc[d['comparison_name'] == "shARID2_over_WT", "log2FoldChange"].abs().sort_values().tail(25).index.tolist()

    # clustermap
    grid = sns.clustermap(
        rnaseq_analysis.expression_annotated.loc[set(g), :].dropna(),
        z_score=0, rasterized=True, xticklabels=rnaseq_analysis.expression_annotated.columns.get_level_values("sample_name"),
        yticklabels=True, cmap="RdBu_r", cbar_kws={"label": "Expression (Z-score)\non extreme genes in miracle plot"}, robust=True)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    grid.ax_heatmap.set_xlabel("Samples")
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=6)
    grid.ax_heatmap.set_ylabel("Genes")
    grid.savefig(
        os.path.join(
            "{}/differential_analysis_{}.only_kd_samples".format(rnaseq_analysis.results_dir, data_type),
            "ARID2_SMARCA4_interaction.miracle_plot.extreme_genes.svg"), bbox_inches="tight", dpi=300)



def joint_heatmap():
    import itertools

    def z_score(x, axis=1):
        """
        Compute a Z-score, defined as (x - mean(x)) / std(x).

        :param numpy.array x: Numeric array.
        """
        return (x - np.nanmean(x, axis=axis)) / np.nanstd(x, axis=axis)


    knockout_genes = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv"), squeeze=True)

    output_dir = "results"
    output_prefix = "all_data_types_combined"

    # Start project
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj._samples = [s for s in prj.samples if s.to_use == "1"]
    for sample in prj.samples:
        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome), "bitSeq", sample.name + ".counts")

    # Sample's attributes
    sample_attributes = ['sample_name', 'cell_line', 'knockout', 'treatment', 'replicate', 'clone', 'batch']
    plotting_attributes = ['knockout', 'treatment', 'replicate', 'clone', 'batch']

    # RNA ANALYSIS
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["HAP1"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]
    rnaseq_analysis = RNASeqAnalysis(name="baf-complex.rnaseq", prj=prj, samples=rnaseq_samples)
    rnaseq_analysis = rnaseq_analysis.from_pickle()

    # RNA-seq logfoldchanges
    fc_table = pd.pivot_table(rnaseq_analysis.differential_results.reset_index(), index="comparison_name", columns="index", values="log2FoldChange")
    fc_table.index.name = "Knockout gene"
    fc_table.columns.name = "Gene"
    fc_table = fc_table.loc[knockout_genes, knockout_genes].dropna()

    # IP-MS data
    # load the data
    arid = pd.read_csv(os.path.join("metadata", "original", "ARID1A-data.csv"), index_col=0).T
    arid = arid.loc[~arid.index.str.contains("ARID1A"), ~arid.columns.str.contains("ARID1A")]
    smarc = pd.read_csv(os.path.join("metadata", "original", "BRG1-data.csv"), index_col=0).T
    smarc = smarc.loc[~smarc.index.str.contains("SMARCA4"), ~smarc.columns.str.contains("SMARCA4")]
    arid.index.name = "knockout"
    arid.columns.name = "subunit"
    smarc.index.name = "knockout"
    smarc.columns.name = "subunit"

    arid.index = arid.index.str.replace("\..*", "")
    arid = arid.groupby(level=0).mean()
    smarc.index = smarc.index.str.replace("\..*", "")
    smarc = smarc.groupby(level=0).mean()

    # arid = np.log2(arid)
    # smarc = np.log2(smarc)

    grid = sns.clustermap(arid.T, cmap="RdBu_r", center=1, metric="correlation", vmin=0, vmax=2, square=True)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.log2FoldChange.ARID1A.svg"), bbox_inches="tight")
    grid = sns.clustermap(smarc.T, cmap="RdBu_r", center=1, metric="correlation", vmin=0, vmax=2, square=True)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.log2FoldChange.SMARCA4.svg"), bbox_inches="tight")

    disc = list()
    all_ = list()
    joint = pd.DataFrame()
    for ko in set(smarc.index.tolist() + arid.index.tolist()):
        for sub in set(smarc.columns.tolist() + arid.columns.tolist()):

            try:
                s = smarc.loc[ko, sub]
            except KeyError:
                s = np.nan
            try:
                a = arid.loc[ko, sub]
            except KeyError:
                a = np.nan

            all_.append([ko, sub, a, s])

            if pd.isnull(s) & pd.isnull(a):
                joint.loc[ko, sub] = np.nan
                del a, s
                continue

            if pd.isnull(s):
                joint.loc[ko, sub] = a
                del a, s
                continue
            if pd.isnull(a):
                joint.loc[ko, sub] = s
                del a, s
                continue

            if (s < 1) and (a < 1):
                joint.loc[ko, sub] = min(a, s)
            elif (s > 1) and (a > 1):
                joint.loc[ko, sub] = max(a, s)
            else:
                disc.append([ko, sub, a, s])
                joint.loc[ko, sub] = np.mean([a, s])

            del a, s

    joint = joint.sort_index(0).sort_index(1)
    joint.index.name = "knockout"
    joint.columns.name = "subunit"

    grid = sns.clustermap(
        joint.fillna(1).T,
        cmap="RdBu_r", cbar_kws={"label": "fold-change"},
        center=1, metric="correlation", vmin=0, vmax=2, square=True, row_cluster=False, col_cluster=False)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.FoldChange.joint.sorted.svg"), bbox_inches="tight")

    grid = sns.clustermap(
        np.log2(joint.fillna(1).T),
        cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        center=0, metric="correlation", vmin=-1, vmax=1, square=True, row_cluster=False, col_cluster=False)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.log2FoldChange.joint.sorted.svg"), bbox_inches="tight")

    grid = sns.clustermap(
        joint.fillna(1).T,
        cmap="RdBu_r", cbar_kws={"label": "fold-change"},
        center=1, metric="correlation", vmin=0, vmax=2, square=True)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.FoldChange.joint.svg"), bbox_inches="tight")

    grid = sns.clustermap(
        np.log2(joint.fillna(1).T),
        cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        center=0, metric="correlation", vmin=-1, vmax=1, square=True)
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.savefig(os.path.join("ip-ms.log2FoldChange.joint.svg"), bbox_inches="tight")

    # 

    # Western
    # get western data
    west = pd.read_csv(os.path.join("metadata", "westernblot.subunit_abundance.csv"), index_col=[0, 1])
    west = west.groupby("gene").mean()
    # west = west.replace(-10, -5)
    west[(west < 0)] = -3
    west[(west > 0)] = 3
    west = west.reindex(index=knockout_genes, columns=knockout_genes)

    # interaction screen
    sel_vars = [
        "CellLine", "GeneSymbol", "Replicate", "Hit",
        "POC", "POC_noOutliers",
        "Z", "Z_noOutliers",
        "TreatmentMean_all", "ControlMean_all", "TreatmentSD_all", "ControlSD_all",
        "pValue_all", "pValue_noOutliers", "pLogP_all", "pLogP_NoOutliers", "DeltaPOC", "deltaPOC_noOutliers"]
    num_vars = sel_vars[4:]

    df = pd.read_csv(os.path.join("metadata", "original", "20171120_inclOutliers_inclHits_3.txt"), sep="\t", decimal=",")
    df = df[sel_vars]
    df['cell_line'] = "HAP1"

    # remove clone 4
    df = df[df['CellLine'] != 'SMARCA4 4']

    # Curate/sanitize
    df.loc[df['CellLine'].str.contains("A549"), "cell_line"] = "A549"
    df['knockout'] = df['CellLine']
    df['knockout'] = df['knockout'].str.replace("A549 ", "")
    df.loc[(df['knockout'] == "WT") & (df['cell_line'] == "A549"), 'knockout'] = "SMARCA4"
    df.loc[df['knockout'] == "SMARCA4 REC", 'knockout'] = "WT"
    df.loc[df['knockout'].str.startswith("ARID1A"), 'knockout'] = "ARID1A"
    df.loc[df['knockout'].str.startswith("SMARCA4"), 'knockout'] = "SMARCA4"
    df["clone"] = df['CellLine'].str.split(" ").apply(lambda x: x[-1])
    df['knockdown'] = df['GeneSymbol']
    for var_ in num_vars:
        df[var_] = df[var_].astype(float)

    # only single knockdowns
    screen = df.loc[
        (~df["knockdown"].astype(str).str.contains("_| ")) &
        (df['cell_line'] == "HAP1"), :]
    screen['knockout'] = screen['knockout'].replace("ACTIN", "ACTB")

    # correct p-value
    from statsmodels.sandbox.stats.multicomp import multipletests
    screen.loc[~screen['pValue_noOutliers'].isnull(), "pValue_noOutliers_FDR"] = multipletests(
        screen.loc[~screen['pValue_noOutliers'].isnull(), "pValue_noOutliers"])[1]

    screen2 = screen.copy()
    screen2.loc[screen2['pValue_noOutliers'] > 0.05, 'deltaPOC_noOutliers'] = 0

    fig, axis = plt.subplots(1, figsize=(4, 4))
    sns.distplot(screen['pLogP_all'].dropna(), ax=axis)
    axis.set_xlabel("-log10(p-value)")
    axis.set_ylabel("Density")
    fig.savefig(os.path.join("screen_table.pLogP_all.svg"), bbox_inches="tight")


    screen_table = pd.pivot_table(
        screen2, index="knockout", columns="knockdown",
        values="deltaPOC_noOutliers", aggfunc=np.mean)  # , aggfunc=lambda x: scipy.stats.combine_pvalues(x)[1])

    # normalize to WT
    # sign = (screen_table >= 0).astype(int).replace(0, -1)
    # screen_table_fc = sign * np.log2(screen_table.abs() / screen_table.loc['WT', :])
    # screen_table_diff = screen_table - screen_table.loc['WT', :]
    # screen_table_fc = screen_table_fc.loc[knockout_genes, knockout_genes]
    # screen_table_diff = screen_table_diff.loc[knockout_genes, knockout_genes]
    # screen_table_fc = screen_table_fc.loc[~screen_table_fc.isnull().all(axis=1), ~screen_table_fc.isnull().all(axis=0)]
    # screen_table_diff = screen_table_diff.loc[~screen_table_diff.isnull().all(axis=1), ~screen_table_diff.isnull().all(axis=0)]

    screen_table = screen_table.loc[knockout_genes, knockout_genes]
    screen_table = screen_table.loc[~screen_table.isnull().all(axis=1), ~screen_table.isnull().all(axis=0)]

    # scale
    # screen_table_z = (((screen_table - 0) / (100)) - 1) * 2
    screen_table_z = screen_table / 10.

    # label and join
    fc_table['data_type'] = "RNA-seq"
    ip_ms_log = np.log2(joint) * 3
    ip_ms_log['data_type'] = "MS"
    west['data_type'] = "Prot"
    screen_table_z['data_type'] = "Synth_screen"
    # screen_table = -screen_table
    # screen_table['data_type'] = "Synth_screen"
    joint_table = ip_ms_log.append(fc_table).append(west).append(screen_table_z)
    joint_table = joint_table.reset_index().set_index(['index', 'data_type']).sort_index()

    # plot
    m = (joint_table
            .drop(['ACTL6A', 'ACTL6B', 'SMARCE1', "SS18", "SS18L1"], axis=0)
            .drop(['ACTL6A', 'ACTL6B', 'SMARCE1', "SS18", "SS18L1"], axis=1))


    # reshape matrix manually
    # df = joint_table
    m3 = pd.melt(joint_table.reset_index(), id_vars=['index', 'data_type'], var_name="gene").rename(columns={"index": "knockout"})

    # m2 = m[m['data_type'] != "Prot"]
    # west.index.name = "knockout"
    # west_m = pd.melt(west.reset_index(), id_vars=['knockout'], var_name="gene").rename(columns={"index": "knockout"})
    # west_m["data_type"] = "Prot"
    # m3 = m2.append(west_m)

    techs = ["Prot", "MS", "RNA-seq", "Synth_screen"]
    q = np.empty((29 * 2, 29 * 2))
    for i, knockout in enumerate(sorted(m3['knockout'].unique())):
        for j, gene in enumerate(sorted(m3['gene'].unique())):

            for k, tech in enumerate(techs):
                if k >= 2: a = 1
                else: a = 0
                if k in [1, 3]: b = 1
                else: b = 0
                v = m3.loc[
                    (m3['knockout'] == knockout) &
                    (m3['gene'] == gene) &
                    (m3['data_type'] == tech), "value"].squeeze()
                q[(i * 2) + a, (j * 2) + b] = v if type(v) is np.float64 else np.nan

    q2 = pd.DataFrame(
        q,
        index=list(itertools.chain.from_iterable(itertools.repeat(x, 2) for x in sorted(m3['knockout'].unique()))),
        columns=list(itertools.chain.from_iterable(itertools.repeat(x, 2) for x in sorted(m3['gene'].unique()))))

    # plot
    # q3 = q2
    q3 = (q2
            .drop(['ACTL6A', 'ACTL6B', 'BCL7C', 'SMARCB1', 'SMARCE1', "SS18", "SS18L1"], axis=0))
            # .drop(['ACTL6A', 'ACTL6B', 'BCL7C', 'SMARCB1', 'SMARCE1', "SS18", "SS18L1"], axis=1))

    sns.set_style("darkgrid")
    matplotlib.rc('font', family='sans-serif') 
    matplotlib.rc('font', serif='Arial') 
    matplotlib.rc('text', usetex='false') 
    matplotlib.rcParams.update({'font.size': 22})

    fig, axis = plt.subplots(1, figsize=(4, 4))
    sns.heatmap(
        q3,
        center=0, vmin=-3, vmax=3, cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        xticklabels=True, yticklabels=True, square=True, ax=axis,
        linewidths=0, 
        linecolor="black")
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize=4)
    axis.set_yticklabels(axis.get_yticklabels(), rotation=0, fontsize=4)
    axis.set_xlabel("Gene")
    axis.set_ylabel("Knockout")
    fig.savefig(os.path.join("all_data_types_combined.all.square.20180411.svg"), bbox_inches="tight")


def subunit_network():
    from scipy.cluster.hierarchy import fcluster
    import networkx as nx

    results_dir = "results"

    # ATAC-seq
    data_type = "ATAC-seq"
    # atac_enr = pd.read_csv(
    #     os.path.join("{}/differential_analysis_{}".format(results_dir, data_type), "differential_analysis.enrichr.csv"))
    # atac_enr.loc[atac_enr['direction'] == "down", 'combined_score'] = -atac_enr.loc[atac_enr['direction'] == "down", 'combined_score']
    # atac_enr_comb = atac_enr.groupby(['gene_set_library', 'description', 'comparison_name'])['combined_score'].mean().reset_index()
    # atac_enr_comb = atac_enr_comb[~atac_enr_comb['comparison_name'].str.contains(r"_|sh|BRD4|SMARCC2")]
    atac_enr = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(results_dir, data_type), "differential_analysis.lola.csv"))
    atac_enr.loc[atac_enr['direction'] == "down", 'pValueLog'] = -atac_enr.loc[atac_enr['direction'] == "down", 'pValueLog']
    atac_enr_comb = atac_enr.groupby(['collection', 'description', 'comparison_name'])['pValueLog'].mean().reset_index()
    atac_enr_comb = atac_enr_comb[~atac_enr_comb['comparison_name'].str.contains(r"_|sh|BRD4|SMARCC2")]

    atac = pd.pivot_table(atac_enr_comb, index="description", columns="comparison_name", values="pValueLog")
    atac_c = atac.corr(method="pearson")
    grid = sns.clustermap(atac_c, cmap="inferno", metric="correlation")
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.ax_heatmap.set_xlabel("Gene")
    grid.ax_heatmap.set_ylabel("Gene")
    grid.savefig(os.path.join("subunit_interplay.ATAC-seq.all_enrichment.svg"), bbox_inches="tight")

    # Make network
    atac_c.index.name = "to"
    net = pd.melt(atac_c.reset_index(), id_vars=["to"], var_name="from")
    net['value'] = 6 ** (1 - net['value'])
    net2 = net[(net['value'] < 2) & (net['value'] != 1)]
    net2['description'] = "ATAC"
    g = nx.from_pandas_dataframe(net2, "from", "to", ['description', "value"])
    nx.write_graphml(g, 'subunit_interplay.network.atac-seq.graphml')

    # cluster subunits
    clusters = pd.Series(
        dict(zip(
            atac_c.index,
            fcluster(grid.dendrogram_col.linkage, t=7, criterion="maxclust"))))
    clusters.loc['BCL7A'] = 5
    clusters.loc['PBRM1'] = 6

    # get terms most specific to each cluster
    diff_enr = pd.DataFrame()
    for cluster in sorted(clusters.unique()):
        own = atac_enr_comb[atac_enr_comb['comparison_name'].isin(clusters[clusters == cluster].index.tolist())]
        rest = atac_enr_comb[~atac_enr_comb['comparison_name'].isin(clusters[clusters == cluster].index.tolist())]
        # get diference of mean across subunits of each cluster
        diff_enr["cluster" + str(cluster) + " " + ",".join(clusters[clusters == cluster].index.tolist())] = (
            own.groupby(['collection', 'description'])['pValueLog'].mean() -
            rest.groupby(['collection', 'description'])['pValueLog'].mean())
    diff_enr.to_csv(os.path.join("subunit_interplay.ATAC-seq.enrichment.20180403.csv"))


    # RNA-seq
    data_type = "RNA-seq"
    rna_enr = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(results_dir, data_type), "differential_analysis.enrichr.csv"))
    rna_enr['p_value'] = -np.log10(rna_enr['p_value'])
    rna_enr.loc[rna_enr['direction'] == "down", 'p_value'] = -rna_enr.loc[rna_enr['direction'] == "down", 'p_value']
    rna_enr_comb = rna_enr.groupby(['gene_set_library', 'description', 'comparison_name'])['p_value'].mean().reset_index()
    rna_enr_comb = rna_enr_comb[~rna_enr_comb['comparison_name'].str.contains(r"_|sh|BRD4|SMARCC2")]

    rna = pd.pivot_table(rna_enr_comb, index="description", columns="comparison_name", values="p_value")
    # rna = pd.pivot_table(rna_enr_comb[rna_enr_comb['gene_set_library'] == gene_set_library], index="description", columns="comparison_name")
    rna_c = rna.corr()
    grid = sns.clustermap(rna_c, cmap="inferno", metric="correlation")
    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.ax_heatmap.set_xlabel("Gene")
    grid.ax_heatmap.set_ylabel("Gene")
    grid.savefig(os.path.join("subunit_interplay.RNA-seq.all_enrichment.svg"), bbox_inches="tight")

    # Make network
    rna_c.index.name = "to"
    net = pd.melt(rna_c.reset_index(), id_vars=["to"], var_name="from")
    net['value'] = 6 ** net['value']
    net2 = net[(net['value'] != 6)].sort_values("value").tail(200)
    net2['description'] = "RNA"
    g = nx.from_pandas_dataframe(net2, "from", "to", ['description', "value"])
    nx.write_graphml(g, 'subunit_interplay.network.rna-seq.graphml')

    # cluster subunits
    clusters = pd.Series(
        dict(zip(
            rna_c.index.get_level_values("to"),
            fcluster(grid.dendrogram_col.linkage, t=10, criterion="maxclust"))))
    clusters.loc['DPF2'] = 9

    # get terms most specific to each cluster
    diff_enr = pd.DataFrame()
    for cluster in sorted(clusters.unique()):
        own = rna_enr_comb[rna_enr_comb['comparison_name'].isin(clusters[clusters == cluster].index.tolist())]
        rest = rna_enr_comb[~rna_enr_comb['comparison_name'].isin(clusters[clusters == cluster].index.tolist())]
        # get diference of mean across subunits of each cluster
        diff_enr["cluster" + str(cluster) + " " + ",".join(clusters[clusters == cluster].index.tolist())] = (
            own.groupby(['gene_set_library', 'description'])['p_value'].mean() -
            rest.groupby(['gene_set_library', 'description'])['p_value'].mean())
    diff_enr.to_csv(os.path.join("subunit_interplay.RNA-seq.enrichment.csv"))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
