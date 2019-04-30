#!/usr/bin/env python

"""
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seaborn as sns

from ngs_toolkit.chipseq import ChIPSeqAnalysis
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.rnaseq import RNASeqAnalysis


# Set settings
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def main():
    # ChIP-seq
    a = ChIPSeqAnalysis(
        name="baf_complex.chipseq",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    a.samples = [s for s in a.samples
                 if s.protocol in ["ChIP-seq", "ChIPmentation"]]
    baf_samples = [s for s in a.samples
                   if ("H3" not in s.ip) and
                      ("CTCF" not in s.ip) and
                      ("IgG" not in s.ip) and
                      ("PolII" not in s.ip)]
    arid1a_ip_samples = [s for s in a.samples
                         if ("ARID1A" in s.ip) and
                            ("r1" in s.replicate)]

    # First, set peak set to ATAC-seq and quantify all ChIP-seq samples in there
    a.set_consensus_sites(os.path.join("results", "baf_complex.atacseq_peak_set.bed"))
    a.measure_coverage()
    a.normalize(method="total")
    # a.annotate(quant_matrix="coverage_rpm")
    a.annotate_with_sample_metadata(
        attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    a.to_pickle()

    # Unsupervised analysis
    a.unsupervised_analysis(
        steps=['correlation'],
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.all_samples",
        standardize_matrix=False)
    a.unsupervised_analysis(
        steps=['pca', 'manifold'],
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.all_samples",
        standardize_matrix=True)

    # Only BAF samples
    a.unsupervised_analysis(
        steps=['correlation'],
        samples=baf_samples,
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.only_baf",
        standardize_matrix=False)
    a.unsupervised_analysis(
        steps=['pca', 'manifold'],
        samples=baf_samples,
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.only_baf",
        standardize_matrix=True)

    # Only ARID1A IP samples
    a.unsupervised_analysis(
        steps=['correlation'],
        samples=arid1a_ip_samples,
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.only_ARID1A_ips",
        standardize_matrix=False)
    a.unsupervised_analysis(
        steps=['pca', 'manifold'],
        samples=arid1a_ip_samples,
        attributes_to_plot=['ip', 'knockout', 'replicate'],
        plot_prefix="chipseq_all_peaks.only_ARID1A_ips",
        standardize_matrix=True)

    # Add ChIP signal to ATAC-seq differential regions
    atac_analysis = ATACSeqAnalysis(
        name="baf_complex.atacseq",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    atac_analysis.samples = [s for s in atac_analysis.samples if s.protocol in ["ATAC-seq"]]
    atac_analysis.load_data(only_these_keys=['differential_results'])
    diff_atac_heatmap_with_chip(atac_analysis, a)

    # Now let's try to use the ChIP-seq peaks on their own
    # read in comparison table
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['cell_type'] == 'HAP1') &
        (comparison_table['comparison_type'] == 'peaks') &
        (~comparison_table['comparison_name'].str.contains("KO"))]

    # call peaks
    a.call_peaks_from_comparisons(comparison_table=comparison_table, overwrite=False)

    # filter peaks
    a.filter_peaks(
        comparison_table,
        filter_bed=os.path.join("data", "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))

    # summarize peaks
    a.peak_summary = a.summarize_peaks_from_comparisons(
        comparison_table=comparison_table)
    a.peak_summary.to_csv(
        os.path.join(a.results_dir, a.name + "_peak_summary.csv"), index=False)

    comparison_table['comparison_genome'] = 'hg19'
    a. get_consensus_sites(
        a,
        comparison_table=comparison_table,
        region_type="peaks",
        blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed")
    # a.calculate_peak_support(comparison_table=c)
    a.measure_coverage()
    a.normalize()
    a.annotate(quant_matrix="coverage_qnorm")
    a.annotate_with_sample_metadata(
        attributes=['sample_name', 'ip', 'cell_line', 'knockout', 'replicate', 'clone'])
    a.to_pickle()

    # Now let's work only with BAF ChIPs
    a = ChIPSeqAnalysis(
        name="baf_complex.chipseq.peaks",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    a.samples = [s for s in a.samples
                 if s.protocol in ["ChIP-seq", "ChIPmentation"]]
    baf_samples = [s for s in a.samples
                   if ("H3" not in s.ip) and
                      ("CTCF" not in s.ip) and
                      ("IgG" not in s.ip) and
                      ("PolII" not in s.ip)]
    arid1a_ip_samples = [s for s in a.samples
                         if ("ARID1A" in s.ip) and
                            ("r1" in s.replicate)]
    a.comparison_table = a.comparison_table[
        (a.comparison_table['cell_type'] == 'HAP1') &
        (a.comparison_table['comparison_type'] == 'peaks') &
        (~a.comparison_table['comparison_name'].str.contains("KO")) &
        (a.comparison_table['comparison_name'].str.contains("ARID|SMARC|PBRM"))]
    a.comparison_table.loc[:, 'comparison_genome'] = 'hg19'
    a.get_consensus_sites(region_type="peaks")
    # a.calculate_peak_support(comparison_table=c)
    a.measure_coverage()
    a.normalize(method="quantile")
    a.annotate(quant_matrix="coverage_qnorm")
    a.annotate_with_sample_metadata()
    a.to_pickle()

    # Unsupervised analysis
    a.unsupervised_analysis(
        a,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks")
    # without histone marks
    a.unsupervised_analysis(
        a, samples=[s for s in a.samples if "H3" not in s.ip],
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks.no_histones")
    # only complex members
    sel_samples = [s for s in a.samples if ("ARID" in s.ip) | ("SMAR" in s.ip) | ("PBRM" in s.ip)]
    a.unsupervised_analysis(
        a, samples=sel_samples,
        attributes_to_plot=['ip', 'replicate'], plot_prefix="chipseq_baf_peaks.only_baf")

    # distinguish between BAF and pBAF specific
    pbaf_vs_baf(a)


def diff_atac_heatmap_with_chip(
        atac_analysis, chipseq_analysis,
        output_dir="results/differential_analysis_ATAC-seq",
        output_prefix="atac-chip_comparison"
        ):
    import pandas as pd
    # import numpy as np
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
    if type(atac_matrix.columns) is pd.core.indexes.multi.MultiIndex:
        atac_matrix = atac_matrix.loc[:, ~atac_matrix.columns.get_level_values("sample_name").str.contains("dBet|sh|parental|BRD4")]
    else:
        atac_matrix = atac_matrix.loc[:, ~atac_matrix.columns.str.contains("dBet|sh|parental|BRD4")]

    # PLOTS
    # Observe values of variables across all comparisons
    all_diff = results[results["diff"].isin([True])].index.drop_duplicates()

    # Sample level
    if type(atac_matrix.columns) is pd.core.indexes.multi.MultiIndex:
        atac_matrix.columns = atac_matrix.columns.get_level_values("sample_name")

    atac_matrix = atac_matrix.loc[all_diff, :]
    n = atac_matrix.isnull().sum().sum()
    if n > 0:
        print("WARNING! {} {} (across all comparisons) were not found in quantification matrix!".format(n, var_name))
        print("Proceeding without those.")
        atac_matrix = atac_matrix.dropna()

    figsize = (max(5, 0.12 * atac_matrix.shape[1]), 5)

    col = sns.clustermap(
        atac_matrix,
        z_score=0, center=0, cmap="RdBu_r",
        xticklabels=False, yticklabels=False, metric="correlation", figsize=(1, 1), rasterized=True, robust=robust)

    g = sns.clustermap(
        atac_matrix,
        col_linkage=col.dendrogram_col.linkage,
        z_score=0,
        xticklabels=True, yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
        cmap="RdBu_r", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(
        output_dir, output_prefix +
        "20190221.diff_{}.samples.clustermap.new.correlation.svg"
        .format(var_name)),
        bbox_inches="tight", dpi=300)

    g = sns.clustermap(
        atac_matrix,
        col_linkage=col.dendrogram_col.linkage,
        z_score=0,
        xticklabels=True, yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
        cmap="RdBu_r", metric="euclidean", figsize=figsize, rasterized=rasterized, robust=robust)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.fig.savefig(os.path.join(
        output_dir, output_prefix +
        "20190221.diff_{}.samples.clustermap.new.euclidean.svg".format(var_name)),
        bbox_inches="tight", dpi=300)

    g2 = sns.clustermap(
        atac_matrix,
        xticklabels=True, yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
        cmap="RdBu_r", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust)
    g2.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.fig.savefig(os.path.join(
        output_dir, output_prefix +
        "20190221.diff_{}.samples.clustermap.z0.svg".format(var_name)),
        bbox_inches="tight", dpi=300)
    plt.close('all')

    # add ChIP-seq data with same clustering
    # but z-scored internally
    regs = atac_matrix.iloc[g.dendrogram_row.reordered_ind, :].index

    background_mean = chipseq_analysis.binding.loc[:, chipseq_analysis.binding.columns.get_level_values("ip").str.contains("IgG|Input")].mean(axis=1)
    chip_over_background = chipseq_analysis.binding.copy()
    for s in chipseq_analysis.binding.columns:
        chip_over_background.loc[:, s] = chipseq_analysis.binding.loc[:, s] - background_mean

    chip_over_background = chip_over_background.drop(
        ['ChIP-seq_HAP1_WT_H3K27ac_r1'], axis=1)

    chip_over_background = chip_over_background.loc[regs, ~chip_over_background.columns.get_level_values("ip").str.contains("IgG|Input")]
    # chip_over_background = chip_over_background.loc[regs, :]

    for label1, c in [
            # ("raw", chipseq_analysis.binding),
            ("chip_over_input", chip_over_background)
            ]:
        for label2, i in [
                ("all", "H3|CTCF|Pol|ARID|BRD|PBRM|SMARC"),
                ("histones", "BRD4|H3|CTCF|Pol"),
                ("baf", "ARID|PBRM|SMARC")
                ]:
            for param, value in [("col_cluster", True), ("col_cluster", False)]:
                c2 = c.loc[:, c.columns.get_level_values(0).str.contains(i)]
                figsize = (0.12 * c2.shape[1], 5)
                kwargs = {
                    "xticklabels": True, "yticklabels": False, "robust": True,
                    "rasterized": rasterized, "row_cluster": False,
                    "metric": "correlation", "figsize": figsize}
                ax_kwargs = {
                    "rotation": 90, "fontsize": "xx-small"}

                g3 = sns.clustermap(c2, **{param: value}, **kwargs)
                g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), **ax_kwargs)
                g3.savefig(os.path.join(
                    output_dir, output_prefix + "20190221.diff_{}.samples+chip.{}.only_{}.clustermap.{}_{}.svg"
                    .format(var_name, label1, label2, param, value)), bbox_inches="tight", dpi=300)

                c_mean = c2.T.groupby(['ip', 'knockout']).mean().T
                figsize = (0.12 * c_mean.shape[1], 5)
                kwargs.update({"figsize": figsize})

                g4 = sns.clustermap(c_mean, **{param: value}, **kwargs)
                g4.ax_heatmap.set_xticklabels(g4.ax_heatmap.get_xticklabels(), **ax_kwargs)
                g4.savefig(os.path.join(
                    output_dir, output_prefix + "20190221.diff_{}.samples+chip.{}.only_{}.mean.clustermap.{}_{}.svg"
                    .format(var_name, label1, label2, param, value)), bbox_inches="tight", dpi=300)
                plt.close('all')

                # Z-score
                c2 = (c2 - c2.mean()) / c2.std()
                figsize = (0.12 * c2.shape[1], 5)
                kwargs.update({"cmap": "PuOr_r", "center": 0, "figsize": figsize})

                g3 = sns.clustermap(c2, **{param: value}, **kwargs)
                g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), **ax_kwargs)
                g3.savefig(os.path.join(
                    output_dir, output_prefix + "20190221.diff_{}.samples+chip.{}.only_{}.z_score.clustermap.{}_{}.svg"
                    .format(var_name, label1, label2, param, value)), bbox_inches="tight", dpi=300)

                c_mean = c2.T.groupby(['ip', 'knockout']).mean().T
                figsize = (0.12 * c_mean.shape[1], 5)
                kwargs.update({"figsize": figsize})

                g4 = sns.clustermap(c_mean, **{param: value}, **kwargs)
                g4.ax_heatmap.set_xticklabels(g4.ax_heatmap.get_xticklabels(), **ax_kwargs)
                g4.savefig(os.path.join(
                    output_dir, output_prefix + "20190221.diff_{}.samples+chip.{}.only_{}.z_score.mean.clustermap.{}_{}.svg"
                    .format(var_name, label1, label2, param, value)), bbox_inches="tight", dpi=300)
                plt.close('all')


def pbaf_vs_baf(chipseq_analysis, output_dir="{results_dir}/pBAF_vs_BAF"):
    from scipy.stats import gaussian_kde

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
    atac_analysis = ATACSeqAnalysis(
        name="baf_complex.atacseq",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    atac_analysis = atac_analysis.from_pickle()
    rnaseq_analysis = RNASeqAnalysis(
        name="baf_complex.atacseq",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    rnaseq_analysis = rnaseq_analysis.from_pickle()

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
            atac_analysis.differential_results[
                ~atac_analysis.differential_results['comparison_name'].str.contains("sh|parental|dBet")]
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
        m = pd.melt(
            rnaseq_analysis.expression_annotated.loc[g]
            .T
            .groupby("knockout").mean()
            .T.reset_index(),
            id_vars=["gene_name"]).rename(columns={"gene_name": "index"})
        m['direction'] = label
        m['value_type'] = "absolute"
        m['data_type'] = "expression"
        enrichment = enrichment.append(m, ignore_index=True)

        # Get fold-change compared to control
        m = (
            rnaseq_analysis.differential_results[
                ~rnaseq_analysis.differential_results['comparison_name'].str.contains("sh|parental|dBet")]
            .loc[g, ["log2FoldChange", "comparison_name"]]
            .rename(columns={"comparison_name": "knockout", "log2FoldChange": "value"}))
        m['direction'] = label
        m['value_type'] = "fold_change"
        m['data_type'] = "expression"
        enrichment = enrichment.append(m, ignore_index=True)

    enrichment.to_csv(os.path.join(output_dir, chipseq_analysis.name +
                      ".baf_peaks.pBAF_vs_BAF.differential.enrichments.csv"), index=False)
    chipseq_analysis.pbaf_enrichment = enrichment
    chipseq_analysis.to_pickle()

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
    fig.savefig(
        os.path.join(output_dir, chipseq_analysis.name +
                     ".baf_peaks.pBAF_vs_BAF.differential.all_enrichments.<10kb.barplot.svg"),
        bbox_inches="tight", dpi=300)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
