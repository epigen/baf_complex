#!/usr/bin/env python

"""
"""

import os
import sys

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

from ngs_toolkit.rnaseq import RNASeqAnalysis, knockout_plot


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc("text", usetex=False)


def main():
    a = RNASeqAnalysis(
        name="baf_complex.rnaseq",
        from_pep=os.path.join("metadata", "project_config.yaml"),
    )
    a._overwride_sample_representation()
    a.samples = [
        s for s in a.samples if (s.protocol == "RNA-seq") and (s.cell_line == "HAP1")
    ]
    for sample in a.samples:
        sample.bitseq_counts = os.path.join(
            sample.paths.sample_root,
            "bowtie1_{}".format(sample.transcriptome),
            "bitSeq",
            sample.name + ".counts",
        )

    # Get and normalize expression data
    a.get_gene_expression(sample_attributes=a.sample_attributes)
    a.to_pickle()

    # # Unsupervised analysis
    # a.unsupervised_analysis()

    # Fix batch effect
    import patsy
    import pandas as pd
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects

    numpy2ri.activate()
    pandas2ri.activate()
    matrix = a.expression_annotated

    robjects.r('require("limma")')
    _removeBatchEffect = robjects.r("removeBatchEffect")
    fixed = _removeBatchEffect(
        x=a.expression_annotated.values,
        batch=matrix.columns.get_level_values("batch"),
        design=patsy.dmatrix("~knockout - 1", matrix.columns.to_frame()),
    )
    a.limma_fixed = pd.DataFrame(
        np.asarray(fixed), index=matrix.index, columns=matrix.columns
    )

    a.unsupervised_analysis(
        quant_matrix="limma_fixed", plot_prefix="limma_fixed", standardize_matrix=False
    )

    # Supervised analysis
    alpha = 0.05
    abs_fold_change = 0

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table["toggle"] == 1)
        & (comparison_table["cell_type"] == "HAP1")
        & (comparison_table["data_type"] == a.data_type)
        & (comparison_table["comparison_type"] == "differential")
    ]
    a.differential_analysis(
        comparison_table,
        samples=[
            s for s in a.samples if s.name in comparison_table["sample_name"].tolist()
        ],
        overwrite=True,
        distributed=True,
    )
    a.collect_differential_analysis(comparison_table, permissive=False, overwrite=True)
    a.to_pickle()

    # Knockout plot
    g = pd.read_csv("metadata/baf_complex_subunits.csv", squeeze=True)
    knockout_plot(
        a,
        knockout_genes=g,
        expression_matrix=a.expression_annotated,
        comparison_results=a.differential_results,
    )

    # Overlap between differential regions
    a.differential_overlap(
        a.differential_results[
            (a.differential_results["padj"] < alpha)
            & (a.differential_results["log2FoldChange"].abs() > abs_fold_change)
        ],
        getattr(a, "limma_fixed").shape[0],
    )

    quant_matrix = "limma_fixed"
    a.plot_differential(
        results=a.differential_results,  # a.differential_results[~a.differential_results['comparison_name'].str.contains("sh|dBet|BRD4")],
        matrix=getattr(a, quant_matrix),
        comparison_table=comparison_table,
        output_prefix="differential_analysis_" + quant_matrix,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
        group_variables=a.group_attributes,
    )

    # Enrichment
    a.differential_enrichment(
        a.differential_results[
            (a.differential_results["padj"] < alpha)
            & (a.differential_results["log2FoldChange"].abs() > abs_fold_change)
        ],
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        distributed=True,
    )

    a.collect_differential_enrichment(
        a.differential_results[
            (a.differential_results["padj"] < alpha)
            & (a.differential_results["log2FoldChange"].abs() > abs_fold_change)
        ],
        directional=True,
        permissive=False,
    )

    a.plot_differential_enrichment(direction_dependent=True)

    synth_lethal()


def synth_lethal():
    import matplotlib.pyplot as plt

    suffix = "synth_lethal"
    a = RNASeqAnalysis(
        name="baf_complex.rnaseq-" + suffix,
        from_pep=os.path.join("metadata", "project_config.yaml"),
    )
    a._overwride_sample_representation()
    a.samples = [
        s
        for s in a.samples
        if (s.protocol == "RNA-seq")
        and (s.cell_line == "HAP1")
        and (
            ("sh" in s.name)
            or ("SMARCA4_ARID2" in s.name)
            or ("ARID2_SMARCA4" in s.name)
        )
    ]
    for sample in a.samples:
        sample.bitseq_counts = os.path.join(
            sample.paths.sample_root,
            "bowtie1_{}".format(sample.transcriptome),
            "bitSeq",
            sample.name + ".counts",
        )
    a.comparison_table = a.comparison_table.loc[
        a.comparison_table["data_type"] == "RNA-seq"
    ]
    a.comparison_table = a.comparison_table.loc[
        a.comparison_table["comparison_name"].str.contains(
            "sh|SMARCA4_ARID2|ARID2_SMARCA4"
        )
    ]

    # prefix = os.path.join(a.results_dir, a.name)
    # kwargs = {"index_col": 0}
    # output_map = {
    #     "matrix_raw": (prefix + ".expression_counts.gene_level.csv", kwargs),
    #     "matrix_norm": (
    #         prefix + ".expression_counts.gene_level.quantile_normalized.log2_tpm.annotated_metadata.csv",
    #         {"index_col": 0, "header": None},
    #     )
    # }
    # a.load_data(output_map=output_map)

    a.get_gene_expression()

    a.differential_analysis(distributed=True)
    a.collect_differential_analysis(output_prefix=suffix)

    a.plot_differential(output_prefix=suffix)

    # take significant genes and do enrichment
    a.differential_enrichment(
        differential=a.differential_results.loc[a.differential_results["padj"] < 0.05],
        max_diff=100,
        output_prefix=suffix,
        distributed=True,
    )
    a.collect_differential_enrichment(input_prefix=suffix)
    a.plot_differential_enrichment(output_prefix=suffix)

    # # have a deeper look at the genes in some of the terms
    enr = a.enrichment_results["enrichr"]
    enr = enr.query("p_value < 0.05")

    genes = (
        enr.loc[
            (enr["gene_set_library"].str.contains("BioCarta")) &
            (enr["description"].str.contains("cell cycle|cyclins|p53|BTG", case=False)),
            "genes"]
        .str.replace(r"[\[\]]", "")
        .str.split(", ")
        .apply(pd.Series)
        .stack()
        .drop_duplicates()
    )

    grid = sns.clustermap(
        a.matrix_norm.reindex(genes).dropna().T,
        robust=True,
        xticklabels=True, yticklabels=a.matrix_norm.columns.get_level_values(0),
        cbar_kws={"label": "Expression (log)"},
        figsize=(25, 4), rasterized=True)
    grid.ax_heatmap.set_xlabel("Cell cycle genes")
    grid.ax_heatmap.set_ylabel("Sample")
    grid.savefig(os.path.join(a.results_dir, a.name + ".enriched_cell_cycle_terms.genes.svg"), dpi=300, bbox_inches="tight")

    grid = sns.clustermap(
        a.matrix_norm.reindex(genes).dropna().T,
        z_score=1, cmap="RdBu_r", center=0, robust=True,
        xticklabels=True, yticklabels=a.matrix_norm.columns.get_level_values(0),
        cbar_kws={"label": "Expression (Z-score)"},
        figsize=(25, 4), rasterized=True)
    grid.ax_heatmap.set_xlabel("Cell cycle genes")
    grid.ax_heatmap.set_ylabel("Sample")
    grid.savefig(os.path.join(a.results_dir, a.name + ".enriched_cell_cycle_terms.genes.z_score.svg"), dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(3, 3))
    egr1 = a.matrix_norm.loc["EGR1"].reset_index(level=list(range(1, 5)), drop=1)
    sns.barplot(egr1, egr1.index, orient="horiz", ax=axis)
    fig.savefig(os.path.join(a.results_dir, a.name + ".enriched_cell_cycle_terms.genes.EGR1.svg"), dpi=300, bbox_inches="tight")

    sb = a.matrix_norm.columns.levels[0][a.matrix_norm.columns.levels[0].str.contains("ARID2_SMARCA4|SMARCA4_ARID2")]
    sa = a.matrix_norm.columns.levels[0][a.matrix_norm.columns.levels[0].str.contains("SMARCA4_shControl|ARID2_shControl")]
    d = a.matrix_norm[sb].mean(1) - a.matrix_norm[sa].mean(1)

    genes = d.nlargest(50).index.tolist() + d.nsmallest(50).index.tolist()

    grid = sns.clustermap(
        a.matrix_norm.reindex(genes).dropna().T,
        metric="euclidean",
        robust=True,
        xticklabels=True, yticklabels=a.matrix_norm.columns.get_level_values(0),
        cbar_kws={"label": "Expression (log)"},
        figsize=(16, 4), rasterized=True)
    grid.ax_heatmap.set_xlabel("Cell cycle genes")
    grid.ax_heatmap.set_ylabel("Sample")
    grid.savefig(os.path.join(a.results_dir, a.name + ".top_100genes.combi_vs_rest.svg"), dpi=300, bbox_inches="tight")

    grid = sns.clustermap(
        a.matrix_norm.reindex(genes).dropna().T,
        metric="correlation",
        z_score=1, cmap="RdBu_r", center=0, robust=True,
        xticklabels=True, yticklabels=a.matrix_norm.columns.get_level_values(0),
        cbar_kws={"label": "Expression (Z-score)"},
        figsize=(16, 4), rasterized=True)
    grid.ax_heatmap.set_xlabel("Cell cycle genes")
    grid.ax_heatmap.set_ylabel("Sample")
    grid.savefig(os.path.join(a.results_dir, a.name + ".top_100genes.combi_vs_rest.z_score.svg"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
