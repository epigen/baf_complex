
import os
from looper.models import Project
import pandas as pd
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


prj = Project(os.path.join("metadata", "project_config.yaml"))
prj.add_sample_sheet()

control_name = "ATAC-seq_HAP1_WT"

output_folder = "/home/arendeiro/baf-kubicek/results/copywriter/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Make normalized visualizations
for resolution in ["100kb", "20kb", "10kb"]:
    igv2 = pd.DataFrame()

    names = list()
    for attrs, index in prj.sheet.df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])

        # Read igv file
        try:
            igv = pd.read_csv(os.path.join("results", "copywriter", resolution + "_" + name, "CNAprofiles", "log2_read_counts.igv"), sep="\t", comment="#").set_index("Feature")
        except:
            print("Sample {} does not exist".format(name))
            continue
        igv.columns = igv.columns.str.replace("log2.", "").str.replace(".merged.sorted.bam", "")

        for sample in igv.drop(['Chromosome', 'Start', 'End'], axis=1).columns:
            igv2[sample] = pd.np.log2(((0.1 + (2 ** igv[sample])) / (0.1 + (2 ** igv[control_name]))))
        names.append(name)

    # save
    igv2 = igv2.join(igv[['Chromosome', 'Start', 'End']]).reset_index()
    igv2 = igv2[['Chromosome', 'Start', 'End', 'Feature'] + sorted(names)]
    tracks = os.path.join(output_folder, "copywriter.{}.normalized.igv".format(resolution))
    open(tracks, 'w')
    output_handle = open(tracks, 'a')
    output_handle.write('#track viewLimits=-2:2 graphType=heatmap color=255,0,0\n')
    igv2.to_csv(output_handle, sep="\t", index=False)
    output_handle.close()

    # Plot sample correlation
    g = sns.clustermap(igv2[names].corr(), cbar_kws={"label": "Pearson correlation"}, cmap="Spectral_r")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_folder, "copywriter.{}.normalized.sample_correlation.svg".format(resolution)), bbox_inches="tight")

    g = sns.clustermap(igv2[names].loc[~igv2['Chromosome'].str.contains("X|Y"), names].corr(), cbar_kws={"label": "Pearson correlation"}, cmap="Spectral_r")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_folder, "copywriter.{}.normalized.sample_correlation.no_sex_chroms.svg".format(resolution)), bbox_inches="tight")

    # Plot variation per chromosome
    var = igv2[names + ["Chromosome"]].groupby("Chromosome").std()
    g = sns.clustermap(var, cbar_kws={"label": "Variation"}, cmap="Greens")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_folder, "copywriter.{}.normalized.chrom_variation.svg".format(resolution)), bbox_inches="tight")

    var = igv2.loc[~igv2['Chromosome'].str.contains("X|Y"), names + ["Chromosome"]].groupby("Chromosome").std()
    g = sns.clustermap(var, cbar_kws={"label": "Variation"}, cmap="Greens")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.savefig(os.path.join(output_folder, "copywriter.{}.normalized.chrom_variation.no_sex_chroms.svg".format(resolution)), bbox_inches="tight")

min_bins = 20
alpha = 0.05


# get annotation
"""
http://grch37.ensembl.org/biomart/martview/3ee51bd6fe6cad2220e6b26d21547a1e?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.band&FILTERS=&VISIBLEPANEL=resultspanel
"""


for min_bins in [1, 2, 4, 8, 16, 32, 64]:
    results = pd.DataFrame()

    for sample in prj.samples:
        print(sample.name)
        try:
            segments = pd.read_csv(os.path.join(output_folder, sample.name, "output.csv"))
        except IOError:
            continue
        segments.columns = ["id", "chrom", "start", "end", "bins", "mean"]

        # Fix coordinate types
        segments[["start", "end"]] = segments[["start", "end"]].astype(int)
        segments["chrom"] = segments["chrom"].astype(str)

        # Fix chromosome names
        segments["chrom"] = "chr" + segments["chrom"]
        segments["chrom"] = segments["chrom"].replace({"chr23": "chrX", "chr24": "chrY"})

        # Fix sample names
        segments[["sample", "background"]] = segments["id"].str.split(".vs.").apply(pd.Series)
        segments["sample"] = segments["sample"].str.replace("log2.", "").str.replace(".trimmed.bowtie2.bam", "").str.replace(".", "-")
        segments["background"] = segments["background"].str.replace("log2.", "").str.replace(".trimmed.bowtie2.bam", "")

        # Get CNV length
        segments["length"] = segments["end"] - segments["start"]

        # Filter by confidence
        filtered = segments[
            (segments["bins"] >= min_bins) &
            (~segments["chrom"].str.contains("chrX|chrY"))]

        # Get calls for sample
        calls = filtered[
            (~filtered["sample"].str.contains("healthy")) &
            (~filtered["background"].str.contains("none"))]

        # Get distribution of variability in healthy donors and compute p-values
        background = filtered[
            (filtered["sample"].str.contains("healthy")) &
            (filtered["background"].str.contains("none"))]["mean"]

        param = scipy.stats.norm.fit(background)
        calls["p_value"] = scipy.stats.norm.sf(abs(calls["mean"]), *param[:-2], loc=param[-2], scale=param[-1]) * 2

        # Correct p-values
        calls["q_value"] = multipletests(calls["p_value"], method="fdr_bh")[1]

        # Annotate
        annotation = pybedtools.BedTool("/home/arendeiro/cll-patients/hg19_annotation.bed")
        calls[["chrom", "start", "end", "id"]].to_csv("tmp.bed", sep="\t", index=False, header=False)
        p = pybedtools.BedTool("tmp.bed").intersect(annotation, wa=True, wb=True).to_dataframe().icol([0, 1, 2, 3, -2, -1])
        p.columns = [["chrom", "start", "end", "id", "gene", "chrom_band"]]

        g = p.groupby(["chrom", "start", "end", "id"])["gene"].apply(lambda x: ",".join(sorted(list(set(x))))).reset_index()
        m = pd.merge(calls, g, on=["chrom", "start", "end", "id"])
        g = p.groupby(["chrom", "start", "end", "id"])["chrom_band"].apply(lambda x: ",".join(sorted(list(set(x))))).reset_index()
        calls = pd.merge(m, g, on=["chrom", "start", "end", "id"])

        # Save
        calls.to_csv(os.path.join(output_folder, sample.name, "calls.csv"), index=False)

        significant = calls[calls["q_value"] < alpha].drop(["id", "background"], axis=1).sort("q_value")
        significant.to_csv(os.path.join(output_folder, sample.name, "calls.significant.csv"), index=False)

        # Append
        results = results.append(significant)

        # Check distributions
        fig, axis = plt.subplots(1, sharex=True, sharey=True)
        sns.distplot(background, label="healthy background", ax=axis)
        sns.distplot(calls["mean"], label=sample.name, ax=axis)
        axis.legend("top left")
        fig.savefig(os.path.join(output_folder, sample.name, "calls.distributions.%i_bins.%f_alpha.png" % (min_bins, alpha)), dpi=300)

        # Volcano plot
        g = sns.jointplot(calls["mean"], -pd.np.log10(calls["q_value"]))
        g.savefig(os.path.join(output_folder, sample.name, "calls.volcano_plot.%i_bins.%f_alpha.png" % (min_bins, alpha)), dpi=300)

    results.to_csv(os.path.join(output_folder, "CNA_calls.all_samples.%i_bins.%f_alpha.csv" % (min_bins, alpha)), index=False)
