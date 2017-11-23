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
    for sample in prj.samples:
        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome), "bitSeq", sample.name + ".counts")

    # Sample's attributes
    sample_attributes = ['sample_name', 'cell_line', 'knockout', 'replicate', 'clone']
    plotting_attributes = ['knockout', 'replicate', 'clone']

    # HAP1 ANALYSIS
    # ATAC-seq
    atacseq_samples = [s for s in prj.samples if (s.library == "ATAC-seq") & (s.cell_line in ["HAP1"])]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    atacseq_samples = [s for s in atacseq_samples if s.knockout != "SMARCC2" and s.clone != "GFP"]
    atac_analysis = ATACSeqAnalysis(name="baf_complex.atacseq", prj=prj, samples=atacseq_samples)
    atac_analysis = main_analysis_pipeline(atac_analysis, data_type="ATAC-seq", cell_type="HAP1")

    # RNA-seq
    rnaseq_samples = [s for s in prj.samples if (s.library == "RNA-seq") & (s.cell_line in ["HAP1"])]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(s.bitseq_counts)]  # and s.pass_qc == 1
    rnaseq_samples = [s for s in rnaseq_samples if s.knockout != "SMARCC2" and s.clone != "GFP"]
    rnaseq_analysis = RNASeqAnalysis(name="baf-complex.rnaseq", prj=prj, samples=rnaseq_samples)
    rnaseq_analysis = main_analysis_pipeline(rnaseq_analysis, data_type="RNA-seq", cell_type="HAP1")


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


    # ChIP-seq data
    # Start project and analysis objects
    chipseq_samples = [s for s in prj.samples if s.library in ["ChIP-seq", "ChIPmentation"]]
    chipseq_samples = [s for s in chipseq_samples if os.path.exists(s.filtered)]
    chipseq_analysis = ChIPSeqAnalysis(
        name="baf_complex.chipseq", prj=prj, samples=chipseq_samples)

    # read in comparison table
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    # call peaks
    chipseq_analysis.call_peaks_from_comparisons(
        comparison_table=comparison_table[comparison_table['comparison_type'] == 'peaks'], overwrite=False)
    # summarize peaks
    chipseq_analysis.summarize_peaks_from_comparisons(comparison_table=comparison_table[comparison_table['comparison_type'] == 'peaks'])

    # set peak set to ATAC-seq and quantify all ChIP-seq samples in there
    chipseq_analysis.set_consensus_sites(atac_analysis.sites.fn)
    chipseq_analysis.measure_coverage()
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

    # Compare with ATAC-seq
    diff_atac_heatmap_with_chip(
        atac_analysis, chipseq_analysis,


    # # Investigate global changes in accessibility
    # global_changes(atacseq_samples)




    # See ATAC and RNA together
    rnaseq_analysis.accessibility_expression()


    # Deeper ATAC-seq data
    nucleosome_changes(atacseq_analysis, atacseq_samples)
    investigate_nucleosome_positions(atacseq_samples)
    phasograms(atacseq_samples)



def main_analysis_pipeline(analysis, data_type, cell_type):
    """
    Main analysis pipeline for ATAC-seq and RNA-seq.
    Gets quantification matrices, normalizes them,
    performes unsupervised and supervised analysis and
    gets and plots enrichments for supervised analysis.
    """

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
        analysis.get_gene_expression(samples=analysis.samples, attributes=["sample_name", "knockout", "replicate", "clone"])

        # see expression of knocked-out genes + other complex members
        baf_genes = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv"), squeeze=True)
        knockout_plot(
            analysis=analysis,
            knockout_genes=baf_genes,
            output_prefix="complex_expression")
        quant_matrix = "expression_annotated"
        feature_name = "genes"


    # Unsupervised analysis
    unsupervised_analysis(
        analysis, quant_matrix=quant_matrix, samples=None,
        attributes_to_plot=plotting_attributes, plot_prefix="all_{}".format(feature_name),
        plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
        output_dir="{results_dir}/unsupervised")
    # only with certain subunits
    if not 'ov90' in analysis.name:
        to_exclude = ['SMARCA4', "SMARCC1", "ARID1A", "ARID1B", "BCL7B"]
        unsupervised_analysis(
            analysis, quant_matrix=quant_matrix,
            samples=[s for s in analysis.samples if s.knockout not in to_exclude],
            attributes_to_plot=plotting_attributes, plot_prefix="all_{}_no_strong_knockouts".format(feature_name),
            plot_max_attr=20, plot_max_pcs=8, plot_group_centroids=True, axis_ticklabels=False, axis_lines=True, always_legend=False,
            output_dir="{results_dir}/unsupervised")


    # Supervised analysis
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[(
        comparison_table['data_type'] == data_type) &
        (comparison_table['cell_type'] == cell_type) &
        (comparison_table['comparison_type'] == 'differential')]
    analysis.differential_results = differential_analysis(
        analysis,
        comparison_table,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        covariates=None,
        alpha=0.05,
        overwrite=True)
    analysis.differential_results = analysis.differential_results.set_index("index")
    analysis.to_pickle()

    if not 'ov90' in analysis.name:
        differential_overlap(
            analysis.differential_results,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            data_type=data_type)

    plot_differential(
        analysis,
        analysis.differential_results,
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        data_type=data_type,
        alpha=0.01,
        corrected_p_value=True,
        fold_change=1,
        rasterized=True, robust=True)

    # repeat without SMARCA4, ARID1A, SMARCC1
    if not 'ov90' in analysis.name:
        plot_differential(
            analysis,
            analysis.differential_results[
                ~analysis.differential_results['comparison_name'].isin(
                    ['SMARCA4', 'ARID1A', 'SMARCC1'])],
            comparison_table=comparison_table,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_very_strong",
            alpha=0.01,
            corrected_p_value=True,
            fold_change=1,
            rasterized=True, robust=True)
        # repeat without the strongest guys
        plot_differential(
            analysis,
            analysis.differential_results[
                ~analysis.differential_results['comparison_name'].isin(
                    ['SMARCA4', 'ARID1A', 'SMARCC1', "BCL7B", "ARID1B"])],
            comparison_table=comparison_table,
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            output_prefix="differential_analysis.no_strong",
            alpha=0.01,
            corrected_p_value=True,
            fold_change=1,
            rasterized=True, robust=True)

    differential_enrichment(
        analysis,
        analysis.differential_results,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        genome="hg19",
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True)
    collect_differential_enrichment(
        analysis.differential_results,
        directional=True,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        permissive=True)

    if data_type == "RNA-seq":
        enrichment_table = pd.read_csv(
            os.path.join("{}/differential_analysis_{}".format(analysis.results_dir, data_type), "differential_analysis.enrichr.csv"))

        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type=data_type,
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            direction_dependent=True,
            top_n=5)
    elif data_type == "ATAC-seq":
        for enrichment_name, enrichment_type in [('motif', 'meme_ame'), ('lola', 'lola'), ('enrichr', 'enrichr')]:
            enrichment_table = pd.read_csv(
                os.path.join("{}/differential_analysis_{}".format(analysis.results_dir, data_type), "differential_analysis" + ".{}.csv".format(enrichment_type)))

            plot_differential_enrichment(
                enrichment_table,
                enrichment_name,
                data_type=data_type,
                output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
                direction_dependent=True,
                top_n=5 if enrichment_name != "motif" else 300)

    return analysis


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
    results['diff']

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[(
        comparison_table['data_type'] == "ATAC-seq") &
        (comparison_table['cell_type'] == "HAP1") &
        (comparison_table['comparison_type'] == 'differential')]

    atac_matrix = atac_analysis.accessibility

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
                yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
                metric="correlation", rasterized=True, figsize=figsize, vmin=0)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(
                groups,
                yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
                metric="correlation", rasterized=True, figsize=figsize)
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
            xticklabels=False, cbar_kws={"label": "Pearson correlation\non fold-changes"},
            cmap="BuGn", vmin=0, vmax=1, metric="correlation", rasterized=True, figsize=(figsize[0], figsize[0]))
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=300, metric="correlation")

        g = sns.clustermap(fold_changes.loc[all_diff, :],
            yticklabels=False, cbar_kws={"label": "Fold-change of\ndifferential {}".format(var_name)},
            robust=True, metric="correlation", rasterized=True, figsize=figsize)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Sample level
    if type(atac_matrix.columns) is pd.core.indexes.multi.MultiIndex:
        atac_matrix.columns = atac_matrix.columns.get_level_values("sample_name")

    atac_matrix = atac_matrix.loc[all_diff, :]
    figsize = (max(5, 0.12 * atac_matrix.shape[1]), 5)

    g2 = sns.clustermap(atac_matrix,
        yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
        xticklabels=True, row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
        vmin=0, metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust)
    g2.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g2.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    g3 = sns.clustermap(atac_matrix,
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True, row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
        metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust)
    g3.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g3.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)


    # add ChIP-seq data with same clustering
    # but z-scored internally
    regs = atac_matrix.iloc[g.dendrogram_row.reordered_ind, :].index
    c = chipseq_analysis.accessibility.loc[regs, :]

    figsize = (0.12 * c.shape[1], 5)

    from scipy.stats import zscore
    chip_z = pd.DataFrame(zscore(c, axis=0), index=regs, columns=chipseq_analysis.accessibility.columns)

    fig, axis = plt.subplots(1, figsize=figsize)
    sns.heatmap(chip_z, yticklabels=False, ax=axis, vmin=-3, vmax=3, rasterized=rasterized)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize="xx-small")
    fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    chip_z = pd.DataFrame(zscore(c, axis=1), index=regs, columns=chipseq_analysis.accessibility.columns)
    fig, axis = plt.subplots(1, figsize=figsize)
    sns.heatmap(chip_z, yticklabels=False, ax=axis, vmin=-3, vmax=3, rasterized=rasterized)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90, fontsize="xx-small")
    fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples+chip.clustermap.z1.svg".format(var_name)), bbox_inches="tight", dpi=300)


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
        self,
        acce_dir="deseq_knockout",
        expr_dir="deseq_expression_knockout",
        trait="knockout"
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

    # read in accessibility
    acce = pd.read_csv(os.path.join("results", acce_dir, acce_dir) + ".%s.annotated.csv" % trait)
    acce['comparison'] = acce['comparison'].str.replace("-WT", "")

    # split accessiblity by gene (shared regulatory elements will count twice)
    acce_split = pd.DataFrame([pd.Series([row['log2FoldChange'], row['comparison'], g])
                                for _, row in acce.iterrows() for g in row['gene_name'].split(',')])
    acce_split.columns = ['log2FoldChange', 'comparison', 'gene_name']
    acce_split.to_csv(os.path.join("results", "accessibility.fold_changes.long_format.gene_split.csv"), index=False)

    # pivot tables
    acce_fc = pd.pivot_table(data=acce_split, values="log2FoldChange", index="gene_name", columns="comparison", aggfunc=signed_max)

    expr = pd.read_csv(os.path.join("results", expr_dir, expr_dir) + ".%s.annotated.csv" % trait)
    expr_fc = pd.pivot_table(data=expr, values="log2FoldChange", index="gene_name", columns="comparison", aggfunc=np.mean)

    # save
    acce_fc.to_csv(os.path.join("results", "accessibility.fold_changes.pivot.gene_split.signed_max.csv"), index=True)
    expr_fc.to_csv(os.path.join("results", "expression.fold_changes.pivot.signed_max.csv"), index=True)

    # read
    acce_fc = pd.read_csv(os.path.join("results", "accessibility.fold_changes.pivot.gene_split.signed_max.csv"))
    expr_fc = pd.read_csv(os.path.join("results", "expression.fold_changes.pivot.signed_max.csv"))

    # match indexes (genes quantified)
    expr_fc = expr_fc.ix[acce_fc.index].dropna()
    acce_fc = acce_fc.ix[expr_fc.index].dropna().drop("WT_GFP", axis=1)

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
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.99_perc.no_lowess_color.svg"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.99_perc.png"), bbox_inches="tight", dpi=300)

    # Plot per percentiles
    # all genes
    n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    dists = pd.DataFrame()
    means = pd.DataFrame()
    for i, ko in enumerate(acce_fc.columns):
        print(ko)
        last_p = 0
        for j, p in enumerate(np.arange(0, 101, 1)):
            if p == 0:
                continue
            b = acce_fc[ko]
            sel_genes = b[(b > np.percentile(b, last_p)) & (b < np.percentile(b, p))].index
            a = expr_fc.loc[sel_genes, ko].to_frame()
            a.columns = ["gene"]
            a["knockout"] = ko
            a["percentile"] = p
            dists = dists.append(a)
            means = means.append(pd.Series([a['gene'].mean(), ko, p], index=["mean", "knockout", "percentile"]), ignore_index=True)
            last_p = p

    g = sns.FacetGrid(data=dists, col="knockout", col_wrap=4)
    g.map(sns.violinplot, data=dists, x="percentile", y="gene")

    g = sns.FacetGrid(data=means, col="knockout", col_wrap=4)
    g.map(plt.scatter, data=means, x="percentile", y="mean", alpha=0.2)

    # Plot per percentiles - same values for every knockout
    # significant in both
    percentils = np.concatenate([np.logspace(-2, np.log10(50), 50, base=10), (100 - np.logspace(-2, np.log10(50), 50, base=10))[::-1]])
    percentil_values = pd.Series({p: np.percentile(expr_fc, p) for p in percentils})

    dists = pd.DataFrame()
    means = pd.DataFrame()
    for i, ko in enumerate(acce_fc.columns):
        print(ko)

        # get values
        a = expr_fc[ko]
        b = acce_fc[ko]
        sig = expr_fc[(abs(a) > 1) & (abs(b) > 1)].index

        last_v = 0
        for p, v in enumerate(percentil_values):
            # if p == 0:
            #     continue
            # b = b[sig]
            sel_genes = b[(b > last_v) & (b < v)].index
            means = means.append(pd.Series([expr_fc.loc[sel_genes, ko].mean(), ko, p], index=["mean", "knockout", "percentile"]), ignore_index=True)
            # a = expr_fc.loc[sel_genes, ko].to_frame()
            # a.columns = ["gene"]
            # a["knockout"] = ko
            # a["percentile"] = p
            # dists = dists.append(a)
            last_v = v

    n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    for i, ko in enumerate(acce_fc.columns):
        axis.flat[i].scatter(means.loc[means['knockout'] == ko, "percentile"], means.loc[means['knockout'] == ko, "mean"])
        axis.flat[i].set_title(ko)

    fig, axis = plt.subplots(1, figsize=(8, 6))
    sns.violinplot(data=dists.groupby(['knockout', 'percentile']).mean().reset_index(), x="percentile", y="gene", ax=axis)
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.all_percentiles.cross_knockouts.svg"), bbox_inches="tight", dpi=300)

    g = sns.FacetGrid(data=dists, col="knockout", col_wrap=4)
    g.map(sns.violinplot, data=dists, x="percentile", y="gene", col="knockout")

    g = sns.FacetGrid(data=means, col="knockout", col_wrap=4)
    g.map(plt.scatter, data=means, x="percentile", y="mean", alpha=0.2)

    # By fixed bins of fold-change
    step = 0.25
    range_max = 5.0
    # bins = zip(-np.arange(step, range_max + step, step)[::-1], -np.arange(0, range_max, step)[::-1])
    bins = zip(np.arange(0, range_max, step), np.arange(step, range_max + step, step))
    # bins = bins[:5] + bins[-5:]

    dists = pd.DataFrame()
    means = pd.DataFrame()
    for i, ko in enumerate(acce_fc.columns):
        print(ko)

        # get values
        a = abs(expr_fc[ko])
        b = abs(acce_fc[ko])
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
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.absolute.cross_knockouts.svg"), bbox_inches="tight", dpi=300)

    means2 = dists.groupby(['max', 'min'])['change'].mean().reset_index()
    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))
    axis.scatter(means2.loc[:, "min"], means2.loc[:, "change"])
    axis.axhline(0, alpha=0.2, color="black", linestyle="--")
    axis.set_ylabel("abs log2 fold-change (RNA-seq)")
    axis.set_xlabel("abs log2 fold-change (ATAC-seq)")
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.absolute.cross_knockouts.svg"), bbox_inches="tight", dpi=300)

    dists2 = dists[dists['min'] <= 3.0]
    fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))
    sns.violinplot(dists2.loc[:, "min"], dists2.loc[:, "change"], trim=True, cut=0, ax=axis)
    axis.set_ylabel("abs log2 fold-change (RNA-seq)")
    axis.set_xlabel("abs log2 fold-change (ATAC-seq)")
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.absolute.cross_knockouts.violinplot.svg"), bbox_inches="tight", dpi=300)

    dists2 = dists[dists['min'] <= 3.0]
    dists2 = dists2.groupby(['knockout', 'min', 'max']).mean().reset_index()
    n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    for i, ko in enumerate(acce_fc.columns):
        axis.flat[i].scatter(dists2.loc[dists2['knockout'] == ko, "min"], dists2.loc[dists2['knockout'] == ko, "change"])
        axis.flat[i].set_title(ko)
    fig.savefig(os.path.join("results", "accessibility-expression.fold_changes.signed_max.absolute.each_knockout.svg"), bbox_inches="tight", dpi=300)

    # n_rows = n_cols = int(np.ceil(np.sqrt(expr_fc.shape[1])))
    # fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(2 * n_rows, 2 * n_cols))
    # for i, ko in enumerate(acce_fc.columns):
    #     axis.flat[i].scatter(means.loc[means['knockout'] == ko, "min"], means.loc[means['knockout'] == ko, "mean"])
    #     axis.flat[i].axhline(0, alpha=0.2, color="black", linestyle="--")
    #     # axis.flat[i].set_xticklabels([x[0] for x in bins])
    #     axis.flat[i].set_title(ko)

    # # Check if genes with high agreement and increase are just already more expressed or accessible to begin with
    # g_up = sig[(pd.DataFrame([a.ix[sig], b.ix[sig]]).T > 0).all(1)]
    # g_down = sig[(pd.DataFrame([a.ix[sig], b.ix[sig]]).T < 0).all(1)]
    # self.expression.ix[g_up]
    # self.expression.ix[g_down]
    # self.accessibility.ix[sig]
    # self.expression.ix[sig]

    # # Check if genes with increasing expression but no increasing accessibility are already more accessible to begin with

    # # Using the native many-to-one relationships for ATAC-seq


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
        # for region_name, bed_file in regions.items():
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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
