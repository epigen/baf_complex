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
from ngs_toolkit.chipseq import ChIPSeqAnalysis, homer_peaks_to_bed
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_analysis, collect_differential_analysis,
                                 differential_enrichment, differential_overlap,
                                 plot_differential,
                                 plot_differential_enrichment,
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
        (s.library in ["ChIP-seq", "ChIPmentation"]) &
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

    # ChIP-seq
    data_type = "ChIP-seq"
    quant_matrix = "binding"
    feature_name = "sites"
    chipseq_analysis = ChIPSeqAnalysis(name="baf_complex.chipseq", prj=prj, samples=prj.samples)

    # read in comparison table
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[comparison_table['comparison_type'] == 'peaks']

    # call peaks
    call_peaks_from_comparisons(chipseq_analysis,
        comparison_table=comparison_table, overwrite=False)
    # summarize peaks
    chipseq_analysis.peak_summary = summarize_peaks_from_comparisons(
        chipseq_analysis,
        comparison_table=comparison_table)
    chipseq_analysis.peak_summary.to_csv(
        os.path.join(chipseq_analysis.results_dir, chipseq_analysis.name + "_peak_summary.csv"), index=False)

    # First, set peak set to ATAC-seq and quantify all ChIP-seq samples in there
    chipseq_analysis.set_consensus_sites(os.path.join("results", "baf_complex.atacseq_peak_set.bed"))
    chipseq_analysis.measure_coverage()
    chipseq_analysis.normalize()
    # chipseq_analysis.annotate(quant_matrix="coverage_qnorm")
    chipseq_analysis.annotate_with_sample_metadata(
        attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    chipseq_analysis.to_pickle()

    # Unsupervised analysis
    unsupervised_analysis(
        chipseq_analysis, data_type=data_type, samples=None,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_all_peaks",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{}/unsupervised_analysis_{}".format(chipseq_analysis.results_dir, data_type), test_pc_association=False)

    # without histone marks
    unsupervised_analysis(
        chipseq_analysis, data_type=data_type, samples=[s for s in chipseq_analysis.samples if "H3" not in s.ip],
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_all_peaks.no_histones",
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{}/unsupervised_analysis_{}".format(chipseq_analysis.results_dir, data_type), test_pc_association=False)


    # Add ChIP signal to ATAC-seq differential regions
    from ngs_toolkit.atacseq import ATACSeqAnalysis
    atac_analysis = ATACSeqAnalysis(name="baf_complex.atacseq").from_pickle()
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


def diff_atac_heatmap_with_chip(
        atac_analysis, chipseq_analysis,
        output_dir="results/differential_analysis_ATAC-seq",
        output_prefix="atac-chip_comparison"):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

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

    atac_matrix = atac_analysis.limma_fixed.copy()
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

    background_mean = chipseq_analysis.binding.loc[:, chipseq_analysis.binding.columns.get_level_values("ip").str.contains("IgG|Input")].mean(axis=1)
    chip_over_background = chipseq_analysis.binding.copy()
    for s in chipseq_analysis.binding.columns:
        chip_over_background.loc[:, s] = chipseq_analysis.binding.loc[:, s] - background_mean

    chip_over_background = chip_over_background.drop(
        ['ChIP-seq_HAP1_WT_H3K27ac_r1'], axis=1)

    c = chip_over_background.loc[regs, ~chip_over_background.columns.get_level_values("ip").str.contains("IgG|Input")]

    figsize = (0.12 * c.shape[1], 5)

    from scipy.stats import zscore
    # chip_z = pd.DataFrame(zscore(c, axis=0), index=regs, columns=chipseq_analysis.binding.columns)
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

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    chip_mean = chipseq_analysis.binding.T.groupby("ip").mean().T
    chip_std = chipseq_analysis.binding.T.groupby("ip").std().dropna().T
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

    stats = stats.dropna()
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
                    alpha=0.2, s=3, color=cmap2(norm2(stats.loc[:, 'density'].values)), rasterized=True)
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
    chip_over_background = chipseq_analysis.binding.copy()
    for s in chipseq_analysis.binding.columns:
        chip_over_background.loc[:, s] = chipseq_analysis.binding.loc[:, s] - background_mean

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
        atac_analysis.gene_annotation['distance'] = atac_analysis.gene_annotation['distance'].astype(int)
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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
