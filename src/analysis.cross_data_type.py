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
from ngs_toolkit.rnaseq import RNASeqAnalysis


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

    # ATAC-seq
    atac_analysis = ATACSeqAnalysis(
        name="baf_complex.atacseq", prj=prj,
        samples=[s for s in prj.samples if s.library == "ATAC-seq"])
    atac_analysis = atac_analysis.from_pickle()

    # RNA-seq
    rnaseq_analysis = RNASeqAnalysis(
        name="baf_complex.rnaseq", prj=prj,
        samples=[s for s in prj.samples if s.library == "RNA-seq"])
    rnaseq_analysis = rnaseq_analysis.from_pickle()


    accessibility_expression(atac_analysis, rnaseq_analysis)

    joint_heatmap(rnaseq_analysis)


def accessibility_expression(
        atac_analysis,
        rnaseq_analysis,
        acce_dir="deseq_knockout",
        expr_dir="deseq_expression_knockout",
        trait="knockout",
        output_dir="results/cross_datatype_comparison",
        output_prefix="cross_datatype_comparison",
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
    top_n = 10
    atac_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(atac_analysis.results_dir, "ATAC-seq"), "differential_analysis.enrichr.csv"))
    rnaseq_enrichment_table = pd.read_csv(
        os.path.join("{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, "RNA-seq"), "differential_analysis.enrichr.csv"))

    for gene_set_library in atac_enrichment_table['gene_set_library'].unique():
        i = len(atac_enrichment_table['comparison_name'].unique())
        fig, axis = plt.subplots(i, 2, figsize=(2 * 4, i * 4))
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
        i = len(atac_enrichment_table['comparison_name'].unique())
        fig, axis = plt.subplots(i, 2, figsize=(2 * 4, i * 4))
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




    # # ARID2 & SMARCA4
    # sig = rnaseq_analysis.differential_results[
    #     (rnaseq_analysis.differential_results['padj'] < alpha) &
    #     (rnaseq_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)]

    # a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "up")].index
    # b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "down")].index
    # a[a.isin(b.tolist())].to_series().to_clipboard(index=False)
    # a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "down")].index
    # b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "up")].index
    # a[a.isin(b.tolist())].to_series().to_clipboard(index=False)

    # a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "up")].index
    # b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "up")].index
    # a[a.isin(b.tolist())].to_series().to_clipboard(index=False)

    # a = sig[(sig['comparison_name'] == "ARID2") & (sig['direction'] == "down")].index
    # b = sig[(sig['comparison_name'] == "SMARCA4") & (sig['direction'] == "down")].index
    # a[a.isin(b.tolist())].to_series().to_clipboard(index=False)


    # # Find the interaction
    # rnaseq_enrichment_table = pd.read_csv(
    #     os.path.join("{}/differential_analysis_{}".format(rnaseq_analysis.results_dir, "RNA-seq"), "differential_analysis.enrichr.csv"))

    # q = rnaseq_enrichment_table.loc[
    #     # (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
    #     (rnaseq_enrichment_table['gene_set_library'].isin(["NCI-Nature_2016"])) &
    #     (rnaseq_enrichment_table['direction'] == 'down') &
    #     (rnaseq_enrichment_table['p_value'] < 0.05) &
    #     (
    #         rnaseq_enrichment_table['description'].str.contains("E2F") |
    #         rnaseq_enrichment_table['description'].str.contains("cell cycle", case=False)), :]

    # genes = q.loc[:, "genes"]

    # cc_genes = list(set(genes.str.replace("[", " ").str.replace(
    #     ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # # clustermap
    # g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    # g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature&WikiPathways.clustermap.svg"), bbox_inches="tight", dpi=300)


    # g = sns.clustermap(rnaseq_analysis.expression_annotated.loc[cc_genes, :].dropna().T.groupby("knockout").mean().T, z_score=0, rasterized=True, square=True)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    # g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.E2F_in_NCI-Nature.group.clustermap.svg"), bbox_inches="tight", dpi=300)



    # # cc_genes = ["CCND3", "RBL1", "CCND2", "CDK2", "CDC25A"]
    # q = rnaseq_enrichment_table.loc[
    #     (rnaseq_enrichment_table['comparison_name'] == 'ARID2') &
    #     # (rnaseq_enrichment_table['gene_set_library'] == 'NCI-Nature_2016') &
    #     (rnaseq_enrichment_table['direction'] == 'down') &
    #     (rnaseq_enrichment_table['p_value'] < 0.05), :]

    # genes = q.loc[:, "genes"]

    # cc_genes = list(set(genes.str.replace("[", " ").str.replace(
    #     ']', ' ').str.replace(', ', ' ').sum().split(' ')))

    # # clustermap
    # g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    # g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.svg"), bbox_inches="tight", dpi=300)

    # # investigate genes in second cluster (ARID2-specific)
    # clusters = scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=2, criterion='maxclust')
    # # plot again just to confirm clusters
    # g2 = sns.clustermap(
    #     rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True,
    #     row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage, row_colors=plt.get_cmap("Paired")(clusters))
    # g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
    # g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), rotation=0)
    # g2.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.clusters_labeled.svg"), bbox_inches="tight", dpi=300)

    # pbaf_genes = pd.Series(g.data.index).iloc[clusters == pd.Series(clusters).value_counts().argmin()].sort_values()

    # g3 = sns.clustermap(rnaseq_analysis.expression.loc[pbaf_genes, :], z_score=0, rasterized=True, metric="correlation")
    # g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xticklabels(), rotation=90)
    # g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_yticklabels(), rotation=0)
    # g3.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.clustermap.pbaf_genes.svg"), bbox_inches="tight", dpi=300)


    # # Cell cycle signature
    # joint = rnaseq_analysis.expression_annotated.loc[cc_genes, :].T.groupby('knockout').mean().mean(1)
    # joint_z = pd.Series(scipy.stats.zscore(joint), index=joint.index)

    # fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4))
    # axis[0].scatter(joint.rank(), joint, cmap="RdBu", vmin=joint.min(), vmax=joint.max())
    # for ko in joint.index:
    #     axis[0].text(joint.rank()[ko], joint[ko], ko)
    # axis[0].axhline(joint.mean(), linestyle='--', color="black")
    # axis[0].set_xlabel("Cell cycle signature (Rank)")
    # axis[0].set_ylabel("Cell cycle signature score")
    # axis[1].scatter(joint_z.rank(), joint_z, cmap="RdBu", vmin=joint_z.min(), vmax=joint_z.max())
    # for ko in joint_z.index:
    #     axis[1].text(joint_z.rank()[ko], joint_z[ko], ko)
    # axis[1].axhline(0, linestyle='--', color="black")
    # axis[1].set_xlabel("Cell cycle signature (Rank)")
    # axis[1].set_ylabel("Cell cycle signature (Z-score)")
    # sns.despine(fig)
    # fig.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_signature.mean_zscore.rank_vs_zscore.svg"), bbox_inches="tight")


    # # Classic cell cycle signature
    # cc_genes = pd.read_table("regev_lab_cell_cycle_genes.txt", header=None, squeeze=True)

    # # clustermap
    # g = sns.clustermap(rnaseq_analysis.expression.loc[cc_genes, :].dropna(), z_score=0, rasterized=True)
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    # g.savefig(os.path.join(output_dir, "ARID2_SMARCA4_interaction.cell_cycle_regev_signature.clustermap.svg"), bbox_inches="tight", dpi=300)



def joint_heatmap(
        rnaseq_analysis,
        ):
    import itertools
    from statsmodels.sandbox.stats.multicomp import multipletests

    def z_score(x, axis=1):
        """
        Compute a Z-score, defined as (x - mean(x)) / std(x).

        :param numpy.array x: Numeric array.
        """
        return (x - np.nanmean(x, axis=axis)) / np.nanstd(x, axis=axis)


    knockout_genes = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv"), squeeze=True)

    output_dir = "results"
    output_prefix = "all_data_types_combined"

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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
