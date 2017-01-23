#!/usr/bin/env python

"""
This is the main script of the baf-kubicek project.
"""
import matplotlib
matplotlib.use('Agg')

# %logstart  # log ipython session
from argparse import ArgumentParser
import os
import sys
from looper.models import Project
import pybedtools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
import multiprocessing
import parmap
import pysam
import numpy as np
import pandas as pd
import cPickle as pickle
from collections import Counter


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def pickle_me(function):
    """
    Decorator for some methods of Analysis class.
    """
    def wrapper(obj, *args):
        function(obj, *args)
        pickle.dump(obj, open(obj.pickle_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    return wrapper


class Analysis(object):
    """
    Class to hold functions and data from analysis.
    """
    def __init__(
            self,
            name="analysis",
            data_dir=os.path.join(".", "data"),
            results_dir=os.path.join(".", "results"),
            pickle_file=None,
            samples=None,
            prj=None,
            from_pickle=False,
            **kwargs):
        # parse kwargs with default
        self.name = name
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.samples = samples
        if pickle_file is None:
            pickle_file = os.path.join(results_dir, "analysis.{}.pickle".format(name))
        self.pickle_file = pickle_file

        for directory in [self.data_dir, self.results_dir]:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # parse remaining kwargs
        self.__dict__.update(kwargs)

        # reload itself if required
        if from_pickle:
            self.__dict__.update(self.from_pickle().__dict__)

    @pickle_me
    def to_pickle(self):
        pass

    def from_pickle(self):
        return pickle.load(open(self.pickle_file, 'rb'))

    @pickle_me
    def get_consensus_sites(self, samples, region_type="peaks", extension=250):
        """Get consensus (union) sites across samples"""
        import re

        for i, sample in enumerate(samples):
            print(sample.name)
            # Get peaks
            try:
                if region_type == "summits":
                    peaks = pybedtools.BedTool(re.sub("_peaks.narrowPeak", "_summits.bed", sample.peaks)).slop(b=extension, genome=sample.genome)
                else:
                    peaks = pybedtools.BedTool(sample.peaks)
            except ValueError:
                print("Peaks for sample {} not found!".format(sample))
                continue
            # Merge overlaping peaks within a sample
            peaks = peaks.merge()
            if i == 0:
                sites = peaks
            else:
                # Concatenate all peaks
                sites = sites.cat(peaks)

        # Merge overlaping peaks across samples
        sites = sites.merge()

        # Filter
        # remove blacklist regions
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))
        # remove chrM peaks and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

        # Read up again
        self.sites = pybedtools.BedTool(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    @pickle_me
    def set_consensus_sites(self, bed_file, overwrite=True):
        """Get consensus (union) sites across samples"""
        self.sites = pybedtools.BedTool(bed_file)
        if overwrite:
            self.sites.saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    @pickle_me
    def calculate_peak_support(self, samples, region_type="peaks"):
        import re

        # calculate support (number of samples overlaping each merged peak)
        for i, sample in enumerate(samples):
            print(sample.name)
            if region_type == "summits":
                peaks = re.sub("_peaks.narrowPeak", "_summits.bed", sample.peaks)
            else:
                peaks = sample.peaks

            if i == 0:
                support = self.sites.intersect(peaks, wa=True, c=True)
            else:
                support = support.intersect(peaks, wa=True, c=True)

        try:
            support = support.to_dataframe()
        except:
            support.saveas("_tmp.peaks.bed")
            support = pd.read_csv("_tmp.peaks.bed", sep="\t", header=None)

        support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.binary_overlap_support.csv"), index=False)

        # get % of total consensus regions found per sample
        m = pd.melt(support, ["chrom", "start", "end"], var_name="sample_name")
        # groupby
        n = m.groupby("sample_name").apply(lambda x: len(x[x["value"] == 1]))

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support[range(len(samples))].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(self.samples)), axis=1)
        # save
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=False)

        self.support = support

    def get_supported_peaks(self, samples):
        """
        Mask peaks with 0 support in the given samples.
        Returns boolean pd.Series of length `peaks`.
        """
        # calculate support (number of samples overlaping each merged peak)
        return self.support[[s.name for s in samples]].sum(1) != 0

    @pickle_me
    def measure_coverage(self, samples):
        def count_reads_in_intervals(bam, intervals):
            """
            Counts reads in a iterable holding strings
            representing genomic intervals of the type chrom:start-end.
            """
            counts = dict()
            bam = pysam.Samfile(bam, 'rb')
            chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

            for interval in intervals:
                if interval.split(":")[0] not in chroms:
                    continue
                counts[interval] = bam.count(region=interval)
            bam.close()

            return counts

        missing = [s for s in samples if not os.path.exists(s.filtered)]
        if len(missing) > 0:
            print("Samples have missing BAM file: {}".format(missing))
            samples = [s for s in samples if s not in missing]

        # Count reads with pysam
        # make strings with intervals
        sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in self.sites]
        # count, create dataframe
        self.coverage = pd.DataFrame(
            map(
                lambda x:
                    pd.Series(x),
                    parmap.map(
                        count_reads_in_intervals,
                        [sample.filtered for sample in samples],
                        sites_str,
                        parallel=True
                    )
            ),
            index=[sample.name for sample in samples]
        ).T

        # Add interval description to df
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            self.coverage.index
        )
        self.coverage["chrom"] = [x[0] for x in ints]
        self.coverage["start"] = [int(x[1]) for x in ints]
        self.coverage["end"] = [int(x[2]) for x in ints]

        # save to disk
        self.coverage.to_csv(os.path.join(self.results_dir, self.name + "_peaks.raw_coverage.csv"), index=True)

    @pickle_me
    def normalize_coverage_quantiles(self, samples):
        def normalize_quantiles_r(array):
            # install R package
            # source('http://bioconductor.org/biocLite.R')
            # biocLite('preprocessCore')

            import rpy2.robjects as robjects
            import rpy2.robjects.numpy2ri
            rpy2.robjects.numpy2ri.activate()

            robjects.r('require("preprocessCore")')
            normq = robjects.r('normalize.quantiles')
            return np.array(normq(array))

        # Quantifle normalization
        to_norm = self.coverage[[s.name for s in samples]]

        self.coverage_qnorm = pd.DataFrame(
            normalize_quantiles_r(to_norm.values),
            index=to_norm.index,
            columns=to_norm.columns
        )
        self.coverage_qnorm = self.coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        self.coverage_qnorm.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_qnorm.csv"), index=True)

    @pickle_me
    def get_peak_gccontent_length(self, bed_file=None, genome="hg19", fasta_file="/home/arendeiro/resources/genomes/{g}/{g}.fa"):
        """
        Bed file must be a 3-column BED!
        """
        if bed_file is None:
            sites = self.sites
        else:
            sites = pybedtools.BedTool(bed_file)

        nuc = sites.nucleotide_content(fi=fasta_file.format(g=genome)).to_dataframe(comment="#")[["score", "blockStarts"]]
        nuc.columns = ["gc_content", "length"]
        nuc.index = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites]

        # get only the sites matching the coverage (not overlapping blacklist)
        self.nuc = nuc.ix[self.coverage.index]

        self.nuc.to_csv(os.path.join(self.results_dir, self.name + "_peaks.gccontent_length.csv"), index=True)

    @pickle_me
    def normalize_gc_content(self, samples):
        def cqn(cov, gc_content, lengths):
            # install R package
            # source('http://bioconductor.org/biocLite.R')
            # biocLite('cqn')
            import rpy2
            rpy2.robjects.numpy2ri.deactivate()

            import rpy2.robjects as robjects
            import rpy2.robjects.pandas2ri
            rpy2.robjects.pandas2ri.activate()

            robjects.r('require("cqn")')
            cqn = robjects.r('cqn')

            cqn_out = cqn(cov, x=gc_content, lengths=lengths)

            y_r = cqn_out[list(cqn_out.names).index('y')]
            y = pd.DataFrame(
                np.array(y_r),
                index=cov.index,
                columns=cov.columns)
            offset_r = cqn_out[list(cqn_out.names).index('offset')]
            offset = pd.DataFrame(
                np.array(offset_r),
                index=cov.index,
                columns=cov.columns)

            return y + offset

        if not hasattr(self, "nuc"):
            self.normalize_coverage_quantiles(samples)

        if not hasattr(self, "nuc"):
            self.get_peak_gccontent_length()

        to_norm = self.coverage_qnorm[[s.name for s in samples]]

        self.coverage_gc_corrected = (
            cqn(cov=to_norm, gc_content=self.nuc["gc_content"], lengths=self.nuc["length"])
            .join(self.coverage[['chrom', 'start', 'end']]))

        self.coverage_gc_corrected.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_gc_corrected.csv"), index=True)

    def get_peak_gene_annotation(self):
        """
        Annotates peaks with closest gene.
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with hg19 TSS positions
        hg19_refseq_tss = pybedtools.BedTool(os.path.join(self.data_dir, "external", "refseq.refflat.tss.bed"))
        # get closest TSS of each cll peak
        gene_annotation = self.sites.closest(hg19_refseq_tss, d=True).to_dataframe()
        gene_annotation = gene_annotation[['chrom', 'start', 'end'] + gene_annotation.columns[-3:].tolist()]  # TODO: check this
        gene_annotation.columns = ['chrom', 'start', 'end', 'gene_name', "strand", 'distance']

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        self.gene_annotation = gene_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.gene_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.gene_annotation.csv"), index=False)

        # save distances to all TSSs (for plotting)
        self.closest_tss_distances = gene_annotation['distance'].tolist()
        pickle.dump(self.closest_tss_distances, open(os.path.join(self.results_dir, self.name + "_peaks.closest_tss_distances.pickle"), 'wb'))

    def get_peak_genomic_location(self):
        """
        Annotates peaks with its type of genomic location.
        Needs files downloaded by prepare_external_files.py
        """
        regions = [
            "ensembl_genes.bed", "ensembl_tss2kb.bed",
            "ensembl_utr5.bed", "ensembl_exons.bed", "ensembl_introns.bed", "ensembl_utr3.bed"]

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        for i, region in enumerate(regions):
            region_name = region.replace(".bed", "").replace("ensembl_", "")
            r = pybedtools.BedTool(os.path.join(self.data_dir, "external", region))
            if region_name == "genes":
                region_name = "intergenic"
                df = self.sites.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
                dfb = background.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
            else:
                df = self.sites.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
                dfb = background.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
            df['genomic_region'] = region_name
            dfb['genomic_region'] = region_name
            if i == 0:
                region_annotation = df
                region_annotation_b = dfb
            else:
                region_annotation = pd.concat([region_annotation, df])
                region_annotation_b = pd.concat([region_annotation_b, dfb])

        # sort
        region_annotation.sort(['chrom', 'start', 'end'], inplace=True)
        region_annotation_b.sort(['chrom', 'start', 'end'], inplace=True)
        # remove duplicates (there shouldn't be anyway)
        region_annotation = region_annotation.reset_index(drop=True).drop_duplicates()
        region_annotation_b = region_annotation_b.reset_index(drop=True).drop_duplicates()
        # join various regions per peak
        self.region_annotation = region_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        self.region_annotation_b = region_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.region_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.region_annotation.csv"), index=False)
        self.region_annotation_b.to_csv(os.path.join(self.results_dir, self.name + "_peaks.region_annotation_background.csv"), index=False)

    def get_peak_chromatin_state(self):
        """
        Annotates peaks with chromatin states.
        (For now states are from CD19+ cells).
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with CD19 chromatin states
        states_cd19 = pybedtools.BedTool(os.path.join(self.data_dir, "external", "HAP1_12_segments.annotated.bed"))

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        # intersect with cll peaks, to create annotation, get original peaks
        chrom_state_annotation = self.sites.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        chrom_state_annotation_b = background.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

        # aggregate annotation per peak, concatenate various annotations (comma-separated)
        self.chrom_state_annotation = chrom_state_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation.columns = ['chrom', 'start', 'end', 'chromatin_state']

        self.chrom_state_annotation_b = chrom_state_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation_b.columns = ['chrom', 'start', 'end', 'chromatin_state']

        # save to disk
        self.chrom_state_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.chromatin_state.csv"), index=False)
        self.chrom_state_annotation_b.to_csv(os.path.join(self.results_dir, self.name + "_peaks.chromatin_state_background.csv"), index=False)

    @pickle_me
    def annotate(self, samples):
        # add closest gene
        self.coverage_annotated = pd.merge(
            self.coverage_gc_corrected,
            self.gene_annotation, on=['chrom', 'start', 'end'], how="left")
        # add genomic location
        self.coverage_annotated = pd.merge(
            self.coverage_annotated,
            self.region_annotation[['chrom', 'start', 'end', 'genomic_region']], on=['chrom', 'start', 'end'], how="left")
        # add chromatin state
        self.coverage_annotated = pd.merge(
            self.coverage_annotated,
            self.chrom_state_annotation[['chrom', 'start', 'end', 'chromatin_state']], on=['chrom', 'start', 'end'], how="left")

        # add support
        self.coverage_annotated = pd.merge(
            self.coverage_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'], how="left")

        # calculate mean coverage
        self.coverage_annotated['mean'] = self.coverage_annotated[[sample.name for sample in samples]].mean(axis=1)
        # calculate coverage variance
        self.coverage_annotated['variance'] = self.coverage_annotated[[sample.name for sample in samples]].var(axis=1)
        # calculate std deviation (sqrt(variance))
        self.coverage_annotated['std_deviation'] = np.sqrt(self.coverage_annotated['variance'])
        # calculate dispersion (variance / mean)
        self.coverage_annotated['dispersion'] = self.coverage_annotated['variance'] / self.coverage_annotated['mean']
        # calculate qv2 (std / mean) ** 2
        self.coverage_annotated['qv2'] = (self.coverage_annotated['std_deviation'] / self.coverage_annotated['mean']) ** 2

        # calculate "amplitude" (max - min)
        self.coverage_annotated['amplitude'] = (
            self.coverage_annotated[[sample.name for sample in samples]].max(axis=1) -
            self.coverage_annotated[[sample.name for sample in samples]].min(axis=1)
        )

        # Pair indexes
        assert self.coverage.shape[0] == self.coverage_annotated.shape[0]
        self.coverage_annotated.index = self.coverage.index

        # Save
        self.coverage_annotated.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_qnorm.log2.annotated.csv"), index=True)

    @pickle_me
    def annotate_with_sample_metadata(
            self,
            attributes=["cell_line", "knockout", "replicate", "clone"]):

        samples = [s for s in self.samples if s.name in self.coverage_annotated.columns]

        attrs = list()
        for attr in attributes:
            l = list()
            for sample in samples:  # keep order of samples in matrix
                try:
                    l.append(getattr(sample, attr))
                except AttributeError:
                    l.append(np.nan)
            attrs.append(l)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=attributes)
        self.accessibility = self.coverage_annotated[[s.name for s in samples]]
        self.accessibility.columns = index

        # Save
        self.accessibility.to_csv(os.path.join(self.results_dir, self.name + ".accessibility.annotated_metadata.csv"), index=True)

    def get_level_colors(self, index=None, levels=None, pallete="Paired", cmap="RdBu_r", nan_color=(0.662745, 0.662745, 0.662745, 0.5)):
        if index is None:
            index = self.accessibility.columns

        if levels is not None:
            index = index.droplevel([l.name for l in index.levels if l.name not in levels])

        _cmap = plt.get_cmap(cmap)
        _pallete = plt.get_cmap(pallete)

        colors = list()
        for level in index.levels:
            # determine the type of data in each level
            most_common = Counter([type(x) for x in level]).most_common()[0][0]
            print(level.name, most_common)

            # Add either colors based on categories or numerical scale
            if most_common in [int, float, np.float32, np.float64, np.int32, np.int64]:
                values = index.get_level_values(level.name)
                # Create a range of either 0-100 if only positive values are found
                # or symmetrically from the maximum absolute value found
                if not any(values.dropna() < 0):
                    norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
                else:
                    r = max(abs(values.min()), abs(values.max()))
                    norm = matplotlib.colors.Normalize(vmin=-r, vmax=r)

                col = _cmap(norm(values))
                # replace color for nan cases
                col[np.where(index.get_level_values(level.name).to_series().isnull().tolist())] = nan_color
                colors.append(col.tolist())
            else:
                n = len(set(index.get_level_values(level.name)))
                # get n equidistant colors
                p = [_pallete(1. * i / n) for i in range(n)]
                color_dict = dict(zip(list(set(index.get_level_values(level.name))), p))
                # color for nan cases
                color_dict[np.nan] = nan_color
                col = [color_dict[x] for x in index.get_level_values(level.name)]
                colors.append(col)

        return colors

    def plot_peak_characteristics(self):
        # Loop at summary statistics:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites], bins=300, kde=False, ax=axis)
        axis.set_xlim(0, 2000)  # cut it at 2kb
        axis.set_xlabel("peak width (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.lengths.svg"), bbox_inches="tight")
        plt.close("all")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis)
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.support.svg"), bbox_inches="tight")
        plt.close("all")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots()
        sns.distplot(self.closest_tss_distances, bins=200, ax=axis)
        axis.set_xlim(0, 100000)  # cut it at 100kb
        axis.set_xlabel("distance to nearest TSS (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.tss_distance.svg"), bbox_inches="tight")
        plt.close("all")

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = [item for sublist in self.region_annotation['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]
        all_region_annotation_b = [item for sublist in self.region_annotation_b['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_region_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_region_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        # plot individually
        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((data[1] / background[1]).astype(float))]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("genomic region")
        axis[2].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.genomic_regions.svg"), bbox_inches="tight")
        plt.close("all")

        # # Plot chromatin states
        # get long list of chromatin states (for plotting)
        all_chrom_state_annotation = [item for sublist in self.chrom_state_annotation['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]
        all_chrom_state_annotation_b = [item for sublist in self.chrom_state_annotation_b['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_chrom_state_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_chrom_state_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((data[1] / background[1]).astype(float))]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("chromatin state")
        axis[2].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.chromatin_states.svg"), bbox_inches="tight")

        # distribution of count attributes
        data = self.coverage_annotated.copy()

        fig, axis = plt.subplots(1)
        sns.distplot(data["mean"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.mean.distplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["qv2"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.qv2.distplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["dispersion"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.dispersion.distplot.svg"), bbox_inches="tight")

        # this is loaded now
        df = pd.read_csv(os.path.join(self.data_dir, self.name + "_peaks.support.csv"))
        fig, axis = plt.subplots(1)
        sns.distplot(df["support"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "baf-kubicek.support.distplot.svg"), bbox_inches="tight")

        plt.close("all")

    def plot_coverage(self):
        data = self.coverage_annotated.copy()
        # (rewrite to avoid putting them there in the first place)
        variables = ['gene_name', 'genomic_region', 'chromatin_state']

        for variable in variables:
            d = data[variable].str.split(',').apply(pd.Series).stack()  # separate comma-delimited fields
            d.index = d.index.droplevel(1)  # returned a multiindex Series, so get rid of second index level (first is from original row)
            data = data.drop([variable], axis=1)  # drop original column so there are no conflicts
            d.name = variable
            data = data.join(d)  # joins on index

        variables = [
            'chrom', 'start', 'end',
            'ensembl_transcript_id', 'distance', 'support',
            'mean', 'variance', 'std_deviation', 'dispersion', 'qv2',
            'amplitude', 'gene_name', 'genomic_region', 'chromatin_state']
        # Plot
        data_melted = pd.melt(
            data,
            id_vars=variables, var_name="sample", value_name="norm_counts")

        # transform dispersion
        data_melted['dispersion'] = np.log2(1 + data_melted['dispersion'])

        # Together in same violin plot
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.mean.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.support.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.mean.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.support.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)
        plt.close("all")

    def plot_variance(self, samples):

        g = sns.jointplot('mean', "dispersion", data=self.coverage_annotated, kind="kde")
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.dispersion.png"), bbox_inches="tight", dpi=300)

        g = sns.jointplot('mean', "qv2", data=self.coverage_annotated)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.qv2_vs_mean.png"), bbox_inches="tight", dpi=300)

        g = sns.jointplot('support', "qv2", data=self.coverage_annotated)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.support_vs_qv2.png"), bbox_inches="tight", dpi=300)

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.coverage_annotated[self.coverage_annotated[[sample.name for sample in samples]].max(axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.dispersion.filtered.png"), bbox_inches="tight", dpi=300)
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.support_vs_qv2.filtered.png"), bbox_inches="tight", dpi=300)

    def unsupervised(
            self, samples, attributes=["cell_line", "knockout", "replicate", "clone"],
            exclude=[]):
        """
        """
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from collections import OrderedDict
        import re
        import itertools
        from scipy.stats import kruskal
        from scipy.stats import pearsonr

        color_dataframe = pd.DataFrame(self.get_level_colors(levels=attributes), index=attributes, columns=[s.name for s in self.samples])

        # exclude samples if needed
        samples = [s for s in samples if s.name not in exclude]
        color_dataframe = color_dataframe[[s.name for s in samples]]
        sample_display_names = color_dataframe.columns.str.replace("_ATAC-seq", "").str.replace("_hg19", "")

        # exclude attributes if needed
        to_plot = attributes[:]
        to_exclude = ["patient_age_at_collection", "sample_name"]
        for attr in to_exclude:
            try:
                to_plot.pop(to_plot.index(attr))
            except:
                continue

        color_dataframe = color_dataframe.ix[to_plot]

        # All regions
        X = self.accessibility[[s.name for s in samples if s.name not in exclude]]

        # Pairwise correlations
        g = sns.clustermap(
            X.corr(), xticklabels=False, yticklabels=sample_display_names, annot=True, aspect=2,
            cmap="Spectral_r", figsize=(15, 15), cbar_kws={"label": "Pearson correlation"}, row_colors=color_dataframe.values.tolist())
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        g.ax_heatmap.set_xlabel(None)
        g.ax_heatmap.set_ylabel(None)
        g.fig.savefig(os.path.join(self.results_dir, "{}.all_sites.corr.clustermap.svg".format(self.name)), bbox_inches='tight')

        # MDS
        mds = MDS(n_jobs=-1)
        x_new = mds.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new)
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        fig, axis = plt.subplots(1, len(to_plot), figsize=(4 * len(to_plot), 4 * 1))
        axis = axis.flatten()
        for i, attr in enumerate(to_plot):
            for j in range(len(xx)):
                try:
                    label = getattr(samples[j], to_plot[i])
                except AttributeError:
                    label = np.nan
                axis[i].scatter(xx.ix[j][0], xx.ix[j][1], s=50, color=color_dataframe.ix[attr][j], label=label)
            axis[i].set_title(to_plot[i])
            axis[i].set_xlabel("MDS 1")
            axis[i].set_ylabel("MDS 2")
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
            axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)

            # Unique legend labels
            handles, labels = axis[i].get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 20:
                if not any([re.match("^\d", c) for c in by_label.keys()]):
                    axis[i].legend(by_label.values(), by_label.keys())
        fig.savefig(os.path.join(self.results_dir, "{}.all_sites.mds.svg".format(self.name)), bbox_inches="tight")

        # PCA
        pca = PCA()
        x_new = pca.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new)
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        # plot % explained variance per PC
        fig, axis = plt.subplots(1)
        axis.plot(
            range(1, len(pca.explained_variance_) + 1),  # all PCs
            (pca.explained_variance_ / pca.explained_variance_.sum()) * 100, 'o-')  # % of total variance
        axis.axvline(len(to_plot), linestyle='--')
        axis.set_xlabel("PC")
        axis.set_ylabel("% variance")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.all_sites.pca.explained_variance.svg".format(self.name)), bbox_inches='tight')

        # plot
        pcs = min(xx.shape[0] - 1, 10)
        fig, axis = plt.subplots(pcs, len(to_plot), figsize=(4 * len(to_plot), 4 * pcs))
        for pc in range(pcs):
            for i, attr in enumerate(to_plot):
                for j in range(len(xx)):
                    try:
                        label = getattr(samples[j], to_plot[i])
                    except AttributeError:
                        label = np.nan
                    axis[pc, i].scatter(xx.ix[j][pc], xx.ix[j][pc + 1], s=50, color=color_dataframe.ix[attr][j], label=label)
                axis[pc, i].set_title(to_plot[i])
                axis[pc, i].set_xlabel("PC {}".format(pc + 1))
                axis[pc, i].set_ylabel("PC {}".format(pc + 2))
                axis[pc, i].set_xticklabels(axis[pc, i].get_xticklabels(), visible=False)
                axis[pc, i].set_yticklabels(axis[pc, i].get_yticklabels(), visible=False)

                # Unique legend labels
                handles, labels = axis[pc, i].get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 20:
                    if not any([re.match("^\d", c) for c in by_label.keys()]):
                        axis[pc, i].legend(by_label.values(), by_label.keys())
        fig.savefig(os.path.join(self.results_dir, "{}.all_sites.pca.svg".format(self.name)), bbox_inches="tight")

        #

        # # Test association of PCs with attributes
        associations = list()
        for pc in range(pcs):
            for attr in attributes[1:]:
                print("PC {}; Attribute {}.".format(pc + 1, attr))
                sel_samples = [s for s in samples if hasattr(s, attr)]
                sel_samples = [s for s in sel_samples if not pd.isnull(getattr(s, attr))]

                # Get all values of samples for this attr
                groups = set([getattr(s, attr) for s in sel_samples])

                # Determine if attr is categorical or continuous
                if all([type(i) in [str, bool] for i in groups]) or len(groups) == 2:
                    variable_type = "categorical"
                elif all([type(i) in [int, float, np.int64, np.float64] for i in groups]):
                    variable_type = "numerical"
                else:
                    print("attr %s cannot be tested." % attr)
                    associations.append([pc + 1, attr, variable_type, np.nan, np.nan, np.nan])
                    continue

                if variable_type == "categorical":
                    # It categorical, test pairwise combinations of attributes
                    for group1, group2 in itertools.combinations(groups, 2):
                        g1_indexes = [i for i, s in enumerate(sel_samples) if getattr(s, attr) == group1]
                        g2_indexes = [i for i, s in enumerate(sel_samples) if getattr(s, attr) == group2]

                        g1_values = xx.loc[g1_indexes, pc]
                        g2_values = xx.loc[g2_indexes, pc]

                        # Test ANOVA (or Kruskal-Wallis H-test)
                        p = kruskal(g1_values, g2_values)[1]

                        # Append
                        associations.append([pc + 1, attr, variable_type, group1, group2, p])

                elif variable_type == "numerical":
                    # It numerical, calculate pearson correlation
                    indexes = [i for i, s in enumerate(samples) if s in sel_samples]
                    pc_values = xx.loc[indexes, pc]
                    trait_values = [getattr(s, attr) for s in sel_samples]
                    p = pearsonr(pc_values, trait_values)[1]

                    associations.append([pc + 1, attr, variable_type, np.nan, np.nan, p])

        associations = pd.DataFrame(associations, columns=["pc", "attribute", "variable_type", "group_1", "group_2", "p_value"])

        # write
        associations.to_csv(os.path.join(self.results_dir, "{}.all_sites.pca.variable_principle_components_association.csv".format(self.name)), index=False)

        # Plot
        # associations[associations['p_value'] < 0.05].drop(['group_1', 'group_2'], axis=1).drop_duplicates()
        # associations.drop(['group_1', 'group_2'], axis=1).drop_duplicates().pivot(index="pc", columns="attribute", values="p_value")
        pivot = associations.groupby(["pc", "attribute"]).min()['p_value'].reset_index().pivot(index="pc", columns="attribute", values="p_value").dropna(axis=1)

        # heatmap of -log p-values
        g = sns.clustermap(-np.log10(pivot), row_cluster=False, annot=True, cbar_kws={"label": "-log10(p_value) of association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.all_sites.pca.variable_principle_components_association.svg".format(self.name)), bbox_inches="tight")

        # heatmap of masked significant
        g = sns.clustermap((pivot < 0.05).astype(int), row_cluster=False, cbar_kws={"label": "significant association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.all_sites.pca.variable_principle_components_association.masked.svg".format(self.name)), bbox_inches="tight")

    def unsupervised_enrichment(self, samples):
        """
        """
        from sklearn.decomposition import PCA
        import itertools
        from statsmodels.sandbox.stats.multicomp import multipletests

        def jackstraw(data, pcs, n_vars, n_iter=100):
            """
            """
            import rpy2.robjects as robj
            from rpy2.robjects.vectors import IntVector
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()

            run = robj.r("""
                run = function(data, pcs, n_vars, B) {
                    library(jackstraw)
                    out <- jackstraw.PCA(data, r1=pcs, r=n_vars, B=B)
                    return(out$p.value)
                }
            """)
            # save to disk just in case
            data.to_csv("_tmp_matrix.jackstraw.csv", index=True)

            if type(pcs) is not int:
                pcs = IntVector(pcs)
            return run(data.values, pcs, n_vars, n_iter)

        def lola(bed_files, universe_file, output_folder):
            """
            Performs location overlap analysis (LOLA) on bedfiles with regions sets.
            """
            import rpy2.robjects as robj

            run = robj.r("""
                function(bedFiles, universeFile, outputFolder) {
                    library("LOLA")

                    # universeFile = "~/baf-kubicek/results/baf-kubicek_peak_set.bed"
                    # bedFiles = "~/baf-kubicek/results/baf-kubicek.ibrutinib_treatment/baf-kubicek.ibrutinib_treatment.timepoint_name.diff_regions.comparison_after_Ibrutinib-before_Ibrutinib.up/baf-kubicek.ibrutinib_treatment.timepoint_name.diff_regions.comparison_after_Ibrutinib-before_Ibrutinib.up_regions.bed"
                    # outputFolder = "~/baf-kubicek/results/baf-kubicek.ibrutinib_treatment/baf-kubicek.ibrutinib_treatment.timepoint_name.diff_regions.comparison_after_Ibrutinib-before_Ibrutinib.up/"

                    userUniverse  <- LOLA::readBed(universeFile)

                    dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/"
                    dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/hg19/"
                    regionDB = loadRegionDB(c(dbPath1, dbPath2))

                    if (typeof(bedFiles) == "character") {
                        userSet <- LOLA::readBed(bedFiles)
                        lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                        writeCombinedEnrichment(lolaResults, outFolder=outputFolder, includeSplits=FALSE)
                    } else if (typeof(bedFiles) == "double") {
                        for (bedFile in bedFiles) {
                            userSet <- LOLA::readBed(bedFile)
                            lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                            writeCombinedEnrichment(lolaResults, outFolder=outputFolder, includeSplits=FALSE)
                        }
                    }
                }
            """)

            # convert the pandas dataframe to an R dataframe
            run(bed_files, universe_file, output_folder)

        def vertical_line(x, **kwargs):
            plt.axvline(x.mean(), **kwargs)

        # Get accessibility matrix excluding sex chroms
        X = self.accessibility[[s.name for s in samples]]
        X = X.ix[X.index[~X.index.str.contains("chrX|chrY")]]

        # Now perform association analysis (jackstraw)
        n_vars = 20
        max_sig = 1000
        alpha = 0.01
        pcs = range(1, n_vars + 1)
        pcs += list(itertools.combinations(pcs, 2))

        p_values = pd.DataFrame(index=X.index)
        for pc in pcs:
            print(pc)
            out = jackstraw(X, pc, n_vars, 10).flatten()
            if type(pc) is int:
                p_values[pc] = out
            else:
                p_values["+".join([str(x) for x in pc])] = out
        q_values = p_values.apply(lambda x: multipletests(x, method="fdr_bh")[1])
        p_values.to_csv(os.path.join("results", "{}.PCA.PC_pvalues.csv".format(self.name)), index=True)

        # Get enrichments of each PC-regions
        lola_enrichments = pd.DataFrame()
        enrichr_enrichments = pd.DataFrame()
        for pc in p_values.columns:
            p = p_values[pc].sort_values()
            sig = p[p < alpha].index

            # Cap to a maximum number of regions
            if len(sig) > max_sig:
                sig = p.head(max_sig).index

            # Run LOLA
            # save to bed
            output_folder = os.path.join("results", "{}.PCA.PC{}_regions".format(self.name, pc))
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            universe_file = os.path.join(output_folder, "universe_sites.bed")
            self.sites.saveas(universe_file)
            bed_file = os.path.join(output_folder, "PCA.PC{}_regions.bed".format(pc))
            self.coverage[['chrom', 'start', 'end']].ix[sig].to_csv(bed_file, sep="\t", header=False, index=False)
            # lola(bed_file, universe_file, output_folder)

            # read lola
            lol = pd.read_csv(os.path.join(output_folder, "allEnrichments.txt"), sep="\t")
            lol["PC"] = pc
            lola_enrichments = lola_enrichments.append(lol)

            # Get genes, run enrichr
            sig_genes = self.coverage_annotated['gene_name'].ix[sig]
            sig_genes = [x for g in sig_genes.dropna().astype(str).tolist() for x in g.split(',')]

            enr = enrichr(pd.DataFrame(sig_genes, columns=["gene_name"]))
            enr["PC"] = pc
            enrichr_enrichments = enrichr_enrichments.append(enr)
        enrichr_enrichments.to_csv(os.path.join("results", "{}.PCA.enrichr.csv".format(self.name)), index=False, encoding="utf-8")
        enrichr_enrichments = pd.read_csv(os.path.join("results", "{}.PCA.enrichr.csv".format(self.name)))
        lola_enrichments.to_csv(os.path.join("results", "{}.PCA.lola.csv".format(self.name)), index=False, encoding="utf-8")

        # Plots

        # p-value distributions
        g = sns.FacetGrid(data=pd.melt(-np.log10(p_values), var_name="PC", value_name="-log10(p-value)"), col="PC", col_wrap=5)
        g.map(sns.distplot, "-log10(p-value)", kde=False)
        g.map(plt.axvline, x=-np.log10(alpha), linestyle="--")
        g.add_legend()
        g.fig.savefig(os.path.join("results", "{}.PCA.PC_pvalues.distplot.svg".format(self.name)), bbox_inches="tight")

        # Volcano plots (loading vs p-value)
        # get PCA loadings
        pca = PCA()
        pca.fit(X.T)
        loadings = pd.DataFrame(pca.components_.T, index=X.index, columns=range(1, X.shape[1] + 1))
        loadings.to_csv(os.path.join("results", "{}.PCA.loadings.csv".format(self.name)), index=True, encoding="utf-8")

        melted_loadings = pd.melt(loadings.reset_index(), var_name="PC", value_name="loading", id_vars=["index"]).set_index(["index", "PC"])
        melted_p_values = pd.melt((-np.log10(p_values)).reset_index(), var_name="PC", value_name="-log10(p-value)", id_vars=["index"]).set_index(["index", "PC"])
        melted = melted_loadings.join(melted_p_values)

        g = sns.FacetGrid(data=melted.dropna().reset_index(), col="PC", col_wrap=5, sharey=False, sharex=False)
        g.map(plt.scatter, "loading", "-log10(p-value)", s=2, alpha=0.5)
        g.map(plt.axhline, y=-np.log10(alpha), linestyle="--")
        g.add_legend()
        g.fig.savefig(os.path.join("results", "{}.PCA.PC_pvalues_vs_loading.scatter.png".format(self.name)), bbox_inches="tight", dpi=300)

        # Plot enrichments
        # LOLA
        # take top n per PC
        import string
        lola_enrichments["set_id"] = lola_enrichments[
            ["collection", "description", "cellType", "tissue", "antibody", "treatment"]].astype(str).apply(string.join, axis=1)

        top = lola_enrichments.set_index('set_id').groupby("PC")['pValueLog'].nlargest(25)
        top_ids = top.index.get_level_values('set_id').unique()

        pivot = pd.pivot_table(
            lola_enrichments,
            index="set_id", columns="PC", values="pValueLog").fillna(0)
        pivot.index = pivot.index.str.replace(" nan", "").str.replace("blueprint blueprint", "blueprint").str.replace("None", "")
        top_ids = top_ids.str.replace(" nan", "").str.replace("blueprint blueprint", "blueprint").str.replace("None", "")

        g = sns.clustermap(
            pivot.ix[top_ids],
            cbar_kws={"label": "p-value z-score"},
            col_cluster=True, z_score=0)
        for tick in g.ax_heatmap.get_xticklabels():
            tick.set_rotation(90)
        for tick in g.ax_heatmap.get_yticklabels():
            tick.set_rotation(0)
        g.fig.savefig(os.path.join("results", "{}.PCA.PC_pvalues.lola_enrichments.svg".format(self.name)), bbox_inches="tight", dpi=300)

        # Enrichr
        for gene_set_library in enrichr_enrichments["gene_set_library"].drop_duplicates():
            enr = enrichr_enrichments[enrichr_enrichments["gene_set_library"] == gene_set_library]

            top = enr.set_index('description').groupby("PC")['p_value'].nsmallest(20)
            top_ids = top.index.get_level_values('description').unique()

            pivot = pd.pivot_table(enr, index="description", columns="PC", values="p_value").fillna(1)
            pivot.index = pivot.index.str.extract("(.*)[,\_\(].*").str.replace("_Homo sapiens", "")
            top_ids = top_ids.str.extract("(.*)[,\_\(].*").str.replace("_Homo sapiens", "")

            g = sns.clustermap(
                -np.log10(pivot.ix[top_ids]), cmap='BuGn',
                cbar_kws={"label": "-log10(p-value)"}, col_cluster=True, figsize=(6, 15))
            for tick in g.ax_heatmap.get_xticklabels():
                tick.set_rotation(90)
            for tick in g.ax_heatmap.get_yticklabels():
                tick.set_rotation(0)
            g.fig.savefig(os.path.join("results", "{}.PCA.PC_pvalues.enrichr_enrichments.{}.svg".format(self.name, gene_set_library)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(
                -np.log10(pivot.ix[top_ids]),
                cbar_kws={"label": "-log10(p-value) z-score"}, col_cluster=True, z_score=0, figsize=(6, 15))
            for tick in g.ax_heatmap.get_xticklabels():
                tick.set_rotation(90)
            for tick in g.ax_heatmap.get_yticklabels():
                tick.set_rotation(0)
            g.fig.savefig(os.path.join("results", "{}.PCA.PC_pvalues.enrichr_enrichments.{}.z_score.svg".format(self.name, gene_set_library)), bbox_inches="tight", dpi=300)


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    """
    import json
    import requests

    ENRICHR_ADD = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_RETRIEVE = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    if gene_set_libraries is None:
        gene_set_libraries = [
            "GO_Biological_Process_2015",
            "GO_Molecular_Function_2015",
            "GO_Cellular_Component_2015",
            # "ChEA_2015",
            "KEGG_2016",
            # "ESCAPE",
            # "Epigenomics_Roadmap_HM_ChIP-seq",
            # "ENCODE_TF_ChIP-seq_2015",
            # "ENCODE_Histone_Modifications_2015",
            "OMIM_Expanded",
            "TF-LOF_Expression_from_GEO",
            "Single_Gene_Perturbations_from_GEO_down",
            "Single_Gene_Perturbations_from_GEO_up",
            "Disease_Perturbations_from_GEO_down",
            "Disease_Perturbations_from_GEO_up",
            "Drug_Perturbations_from_GEO_down",
            "Drug_Perturbations_from_GEO_up",
            "WikiPathways_2016",
            "Reactome_2016",
            "BioCarta_2016",
            "NCI-Nature_2016"
        ]
    results = pd.DataFrame()
    for gene_set_library in gene_set_libraries:
        print("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(list(set(dataframe["gene_name"].dropna().tolist())))
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(list(set(dataframe[['chrom', 'start', 'end']].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())))

        payload = {
            'list': (None, attr),
            'description': (None, gene_set_library)
        }
        # Request adding gene set
        response = requests.post(ENRICHR_ADD, files=payload)
        if not response.ok:
            raise Exception('Error adding gene list')

        # Track gene set ID
        user_list_id = json.loads(response.text)['userListId']

        # Request enriched sets in gene set
        response = requests.get(
            ENRICHR_RETRIEVE + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if len(res.columns) == 7:
            res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value"]
        elif len(res.columns) == 9:
            res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value", "old_p_value", "old_adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results

    # for F in `find . -iname *_genes.symbols.txt`
    # do
    #     if  [ ! -f ${F/symbols.txt/enrichr.csv} ]; then
    #         echo $F
    #         sbatch -J ENRICHR_${F} -o ${F/symbols.txt/}enrichr.log ~/run_Enrichr.sh $F
    #     fi
    # done


def add_args(parser):
    """
    Options for project and pipelines.
    """
    # Behaviour
    parser.add_argument("-g", "--generate", dest="generate", action="store_true",
                        help="Should we generate data and plots? Default=False")

    return parser


def main():
    # Parse arguments
    parser = ArgumentParser()
    parser = add_args(parser)
    args = parser.parse_args()

    #
    # First an analysis of the Ibrutinib samples in light of the previous CLL and other hematopoietic cells
    analysis = Analysis(name="baf-kubicek", from_pickle=args.generate)

    # Start project
    prj = Project("metadata/project_config.yaml")
    prj.add_sample_sheet()
    for sample in prj.samples:
        sample.subproject = analysis.name

    # temporary:
    for sample in prj.samples:
        sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        sample.summits = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_summits.bed")
        sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
        sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
        sample.coverage = os.path.join(sample.paths.sample_root, "coverage", sample.name + ".cov")

    # work only with ATAC-seq samples
    prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq" and sample.cell_line in ["HAP1"]]
    prj.samples = [s for s in prj.samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1

    # pair analysis and Project
    analysis.prj = prj
    analysis.samples = prj.samples

    # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT
    # Get consensus peak set from all samples
    analysis.get_consensus_sites(analysis.samples, "summits")
    analysis.calculate_peak_support(analysis.samples, "summits")

    # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
    # Get coverage values for each peak in each sample of ATAC-seq and ChIPmentation
    analysis.measure_coverage(analysis.samples)
    # normalize coverage values
    analysis.normalize_coverage_quantiles(analysis.samples)
    analysis.get_peak_gccontent_length()
    analysis.normalize_gc_content(analysis.samples)

    # Annotate peaks with closest gene
    analysis.get_peak_gene_annotation()
    # Annotate peaks with genomic regions
    analysis.get_peak_genomic_location()
    # Annotate peaks with ChromHMM state from CD19+ cells
    analysis.get_peak_chromatin_state()
    # Annotate peaks with closest gene, chromatin state,
    # genomic location, mean and variance measurements across samples
    analysis.annotate(analysis.samples)
    analysis.annotate_with_sample_metadata()

    # Plots
    # plot general peak set features
    analysis.plot_peak_characteristics()
    # Plot coverage features across peaks/samples
    analysis.plot_coverage()
    analysis.plot_variance()

    # Unsupervised analysis
    analysis.unsupervised(attributes=["cell_line", "knockout", "replicate", "clone"])

    # analysis.sites = pybedtools.BedTool(os.path.join(analysis.results_dir, analysis.name + "_peak_set.bed"))
    # analysis.support = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.support.csv"))
    # analysis.coverage = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.raw_coverage.csv"), index_col=0)
    # analysis.gene_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.gene_annotation.csv"))
    # analysis.closest_tss_distances = pickle.load(open(os.path.join(analysis.results_dir, analysis.name + "_peaks.closest_tss_distances.pickle"), "rb"))
    # analysis.region_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.region_annotation.csv"))
    # analysis.region_annotation_b = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.region_annotation_background.csv"))
    # analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.chromatin_state.csv"))
    # analysis.chrom_state_annotation_b = pd.read_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.chromatin_state_background.csv"))

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
