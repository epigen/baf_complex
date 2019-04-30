#!/usr/bin/env python

"""
"""

import os
import sys

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

from ngs_toolkit.atacseq import ATACSeqAnalysis


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def main():
    atac_analysis = ATACSeqAnalysis(
        name="baf_complex.atacseq",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    atac_analysis._overwride_sample_representation()
    atac_analysis.samples = [s for s in atac_analysis.samples
                             if (s.protocol == "ATAC-seq") and
                                (s.cell_line == "HAP1")]

    # Usual pipeline
    atac_analysis.get_consensus_sites()
    atac_analysis.measure_coverage()
    atac_analysis.normalize(method="gc_content")
    atac_analysis.annotate_with_sample_metadata(
        quant_matrix="coverage_gc_corrected")
    atac_analysis.to_pickle()

    # Unsupervised analysis
    atac_analysis.unsupervised_analysis()

    # Fix batch effect
    import patsy
    import pandas as pd
    import rpy2
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r('require("limma")')
    _removeBatchEffect = robjects.r('removeBatchEffect')

    fixed = _removeBatchEffect(
        x=atac_analysis.accessibility.values,
        batch=atac_analysis.accessibility.columns.get_level_values("batch"),
        design=patsy.dmatrix("~knockout - 1", atac_analysis.accessibility.columns.to_frame()))
    atac_analysis.limma_fixed = pd.DataFrame(
        np.asarray(fixed),
        index=atac_analysis.accessibility.index,
        columns=atac_analysis.accessibility.columns)
    atac_analysis.limma_fixed.to_csv(os.path.join(atac_analysis.results_dir, atac_analysis.name + "_peaks.limma_fixed.csv"))
    atac_analysis.to_pickle()

    atac_analysis.unsupervised_analysis(
        quant_matrix='limma_fixed', plot_prefix="limma_fixed")

    # Supervised analysis
    alpha = 0.01
    abs_fold_change = 1
    quant_matrix = 'limma_fixed'

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['toggle'] == 1) &
        (comparison_table['data_type'] == atac_analysis.data_type) &
        (comparison_table['comparison_type'] == 'differential')]

    atac_analysis.differential_analysis(
        overwrite=True, distributed=True)

    atac_analysis.collect_differential_analysis(
        permissive=False,
        overwrite=True)
    atac_analysis.to_pickle()

    # Overlap between differential regions
    atac_analysis.differential_overlap(
        atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        getattr(atac_analysis, quant_matrix).shape[0])

    atac_analysis.plot_differential(
        results=atac_analysis.differential_results,
        # results=atac_analysis.differential_results[~atac_analysis.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")],
        output_prefix="differential_analysis_" + quant_matrix,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
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

    atac_analysis.differential_enrichment(
        atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True)

    atac_analysis.collect_differential_enrichment(
        differential=atac_analysis.differential_results[
            (atac_analysis.differential_results['padj'] < alpha) &
            (atac_analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        permissive=True)

    atac_analysis.plot_differential_enrichment(
        direction_dependent=True)

    # Knockout plot
    from ngs_toolkit.rnaseq import knockout_plot
    g = pd.read_csv("metadata/baf_complex_subunits.csv", squeeze=True)
    atac_g = atac_analysis.get_gene_level_accessibility()
    knockout_plot(
        atac_analysis, knockout_genes=g, expression_matrix=atac_g,
        output_prefix="knockout_expression.atacseq_gene_level")

    synth_lethal()


def synth_lethal():
    from ngs_toolkit.utils import bed_to_index
    import matplotlib.pyplot as plt

    # Load PBAF or BAF-specific regions
    baf = os.path.join("results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.BAF_specific.bed")
    pbaf = os.path.join("results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.pBAF_specific.bed")
    output = os.path.join("results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.both_sites.bed")

    os.system("cat {} {} > {}".format(baf, pbaf, output))

    a = ATACSeqAnalysis(
        name="baf_complex.atacseq-" + "BAF-PBAF_specific",
        from_pep=os.path.join("metadata", "project_config.yaml"))
    a.samples = [
        s for s in a.samples if
        # (("ARID2" in s.name) or ("SMARCA4" in s.name)) and
        (s.knockout in ["SMARCA4", "ARID2", "WT"]) and
        ("parental" not in s.name) and
        ("dBet" not in s.name) and
        (s.protocol == "ATAC-seq") and
        (s.to_use == "1") and
        (s.cell_line == "HAP1")]

    # Get accessibility at ChIP-seq sites
    a.set_consensus_sites(output)
    a.measure_coverage(distributed=True, peak_set_name="BAF-PBAF")
    a.collect_coverage(peak_set_name="BAF-PBAF")
    a.normalize(method="quantile")

    # Get a color vector depending on site belonging to BAF or PBAF
    baf = bed_to_index(pd.read_csv(baf, sep="\t", header=None, names=['chrom', 'start', 'end']))
    pbaf = bed_to_index(pd.read_csv(pbaf, sep="\t", header=None, names=['chrom', 'start', 'end']))

    g = sns.clustermap(
        a.coverage_qnorm.T,
        col_colors=plt.get_cmap("Paired")(a.coverage_qnorm.index.isin(baf).astype(int) * 256),
        robust=True,
        xticklabels=False, yticklabels=True, rasterized=True,
        metric="correlation")
    g.savefig(os.path.join(
        "results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.both.ATAC-seq_signal.svg"),
        bbox_inches="tight", dpi=300)
    g = sns.clustermap(
        a.coverage_qnorm.T,
        col_colors=plt.get_cmap("Paired")(a.coverage_qnorm.index.isin(baf).astype(int) * 256),
        robust=True,
        xticklabels=False, yticklabels=True, rasterized=True,
        metric="correlation", z_score=1, cmap="RdBu_r", center=0)
    g.savefig(os.path.join(
        "results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.both.ATAC-seq_signal.z_score.svg"),
        bbox_inches="tight", dpi=300)

    baf_g = sns.clustermap(
        a.coverage_qnorm.loc[baf].T,
        robust=True,
        xticklabels=False, yticklabels=True, rasterized=True,
        metric="correlation", z_score=1, cmap="RdBu_r", center=0)
    baf_g.savefig(os.path.join(
        "results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.BAF_specific.ATAC-seq_signal.svg"),
        bbox_inches="tight", dpi=300)
    pbaf_g = sns.clustermap(
        a.coverage_qnorm.loc[pbaf].T,
        robust=True,
        xticklabels=False, yticklabels=True, rasterized=True,
        metric="correlation", z_score=1, cmap="RdBu_r", center=0)
    pbaf_g.savefig(os.path.join(
        "results", "pBAF_vs_BAF", "baf_complex.chipseq.baf_peaks.pBAF_vs_BAF.pBAF_specific.ATAC-seq_signal.svg"),
        bbox_inches="tight", dpi=300)

    # Get differences in ATAC-seq regions
    suffix = "synth_lethal"
    a = ATACSeqAnalysis(
        name="baf_complex.atacseq-" + suffix,
        from_pep=os.path.join("metadata", "project_config.yaml"))
    a._overwride_sample_representation()
    a.samples = [s for s in a.samples if
                 (s.protocol == "ATAC-seq") and
                 (s.cell_line == "HAP1") and
                 ("sh" in s.name)]
    a.comparison_table = a.comparison_table.loc[a.comparison_table['data_type'] == "ATAC-seq"]
    a.comparison_table = a.comparison_table.loc[a.comparison_table['comparison_name'].str.contains("sh")]

    a.set_consensus_sites(os.path.join("results", "baf_complex.atacseq_peak_set.bed"))
    a.measure_coverage(distributed=True, peak_set_name=suffix)
    a.collect_coverage(peak_set_name=suffix)
    a.normalize(method="quantile")
    a.annotate(matrix="coverage_qnorm")
    a.annotate_with_sample_metadata(matrix="coverage_qnorm")

    a.calculate_peak_support()
    a.differential_analysis(distributed=True)
    a.collect_differential_analysis(output_prefix=suffix)

    a.plot_differential(output_prefix=suffix)

    # # take top N regions and do enrichment
    a.differential_enrichment(max_diff=100, output_prefix=suffix, distributed=True)
    a.collect_differential_enrichment(input_prefix=suffix)
    a.plot_differential_enrichment(output_prefix=suffix)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
