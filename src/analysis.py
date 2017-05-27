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
        try:
            pickle.dump(obj, open(obj.pickle_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        except:
            pass
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

    def get_peak_genomic_location(self, genome="hg19"):
        """
        Annotates peaks with its type of genomic location.
        Needs files downloaded by prepare_external_files.py
        """
        regions = [
            "ensembl_genes.bed", "ensembl_tss2kb.bed",
            "ensembl_utr5.bed", "ensembl_exons.bed", "ensembl_introns.bed", "ensembl_utr3.bed"]

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=genome, chrom=True)

        if not hasattr(self, "genome_space"):
            self.genome_space = dict()

        for i, region in enumerate(regions):
            region_name = region.replace(".bed", "").replace("ensembl_", "")
            r = pybedtools.BedTool(os.path.join(self.data_dir, "external", region))

            # sum up the total space of region locations in the genome
            self.genome_space[region_name] = sum([i.length for i in r])

            if region_name == "genes":
                region_name = "intergenic"

                # estimate intergenic space as remainder of non-gene space
                self.genome_space[region_name] = sum([x[1][1] for x in pybedtools.chromsizes(genome).items()]) - self.genome_space["genes"]

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
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with chromatin states
        chrom_states = pybedtools.BedTool(os.path.join(self.data_dir, "external", "HAP1_12_segments.annotated.bed"))

        # sum up the total space of region locations in the genome
        if not hasattr(self, "genome_space"):
            self.genome_space = dict()
        c = chrom_states.to_dataframe()
        self.genome_space.update(c.groupby(['name']).apply(lambda x: x['end'] - x['start']).groupby(level=0).sum().to_dict())

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        # intersect with cll peaks, to create annotation, get original peaks
        chrom_state_annotation = self.sites.intersect(chrom_states, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        chrom_state_annotation_b = background.intersect(chrom_states, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

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
            attributes=["sample_name", "cell_line", "knockout", "replicate", "clone"]):

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

    def collect_bitseq_output(self, samples):
        first = True
        for i, sample in enumerate(samples):
            if first:
                try:
                    # read the "tr" file of one sample to get indexes
                    tr = pd.read_csv(
                        os.path.join(
                            sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                            "bitSeq",
                            sample.name + ".tr"),
                        sep=" ", header=None, skiprows=1,
                        names=["ensembl_gene_id", "ensembl_transcript_id", "v1", "v2"])
                except IOError:
                    print("Sample {} is missing.".format(sample.name))
                    continue
                # add id index
                tr.set_index("ensembl_gene_id", append=False, inplace=True)
                tr.set_index("ensembl_transcript_id", append=True, inplace=True)
                # create dataframe
                expr = pd.DataFrame(index=tr.index)
                first = False

            # load counts
            try:
                counts = pd.read_csv(os.path.join(
                    sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                    "bitSeq",
                    sample.name + ".counts"), sep=" ")
            except IOError:
                print("Sample {} is missing.".format(sample.name))
                continue
            counts.index = expr.index

            # Append
            expr[sample.name] = counts

        return expr

    def get_gene_expression(self, samples=None, attributes=["sample_name", "cell_line", "knockout", "replicate", "clone", "complex_subunit"]):
        """
        """
        def quantile_normalize(df_input):
            df = df_input.copy()
            # compute rank
            dic = {}
            for col in df:
                dic.update({col: sorted(df[col])})
            sorted_df = pd.DataFrame(dic)
            rank = sorted_df.mean(axis=1).tolist()
            # sort
            for col in df:
                t = np.searchsorted(np.sort(df[col]), df[col])
                df[col] = [rank[i] for i in t]
            return df

        if samples is None:
            samples = [s for s in self.samples if s.library == "RNA-seq"]

        self.count_matrix = self.collect_bitseq_output(samples)

        # Map ensembl gene IDs to gene names
        import requests
        url_query = "".join([
            """http://grch37.ensembl.org/biomart/martservice?query=""",
            """<?xml version="1.0" encoding="UTF-8"?>""",
            """<!DOCTYPE Query>""",
            """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
            """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
            """<Attribute name = "ensembl_gene_id" />""",
            """<Attribute name = "external_gene_name" />""",
            """</Dataset>""",
            """</Query>"""])
        req = requests.get(url_query, stream=True)
        mapping = pd.DataFrame((x.strip().split(",") for x in list(req.iter_lines())), columns=["ensembl_gene_id", "gene_name"])
        self.count_matrix = self.count_matrix.reset_index()
        self.count_matrix['ensembl_gene_id'] = self.count_matrix['ensembl_gene_id'].str.replace("\..*", "")

        self.count_matrix = pd.merge(self.count_matrix, mapping, how="outer").dropna()
        self.count_matrix.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.csv"), index=False)

        # Quantile normalize
        self.matrix_qnorm = quantile_normalize(self.count_matrix.drop(["ensembl_transcript_id", "ensembl_gene_id", "gene_name"], axis=1))

        # Log2 TPM
        self.matrix_qnorm_log = np.log2((1 + self.matrix_qnorm) / self.matrix_qnorm.sum(axis=0) * 1e6)
        self.matrix_qnorm_log = self.count_matrix[["gene_name", "ensembl_gene_id", "ensembl_transcript_id"]].join(self.matrix_qnorm_log)
        self.matrix_qnorm_log.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.transcript_level.quantile_normalized.log2_tpm.csv"), index=False)

        # Reduce to gene-level measurements by max of transcripts
        self.expression_matrix_counts = self.count_matrix.drop(['ensembl_transcript_id'], axis=1).dropna().set_index(['ensembl_gene_id', 'gene_name']).groupby(level=["gene_name"]).max()

        # Quantile normalize
        self.expression_matrix_qnorm = quantile_normalize(self.expression_matrix_counts)

        # Log2 TPM
        self.expression = np.log2((1 + self.expression_matrix_qnorm) / self.expression_matrix_qnorm.sum(axis=0) * 1e6)

        # Annotate with sample metadata
        _samples = [s for s in samples if s.name in self.expression.columns]
        attrs = list()
        for attr in attributes:
            l = list()
            for sample in _samples:  # keep order of samples in matrix
                try:
                    l.append(getattr(sample, attr))
                except AttributeError:
                    l.append(np.nan)
            attrs.append(l)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=attributes)
        self.expression_annotated = self.expression[[s.name for s in _samples]]
        self.expression_annotated.columns = index

        # Save
        self.expression_matrix_counts.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.csv"), index=True)
        self.expression_matrix_qnorm.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.csv"), index=True)
        self.expression.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.log2_tpm.csv"), index=True)
        self.expression_annotated.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.log2_tpm.annotated_metadata.csv"), index=True)

    def plot_peak_characteristics(self, samples=None):
        def get_sample_reads(bam_file):
            return pysam.AlignmentFile(bam_file).count()

        def get_peak_number(bed_file):
            return len(open(bed_file, "r").read().split("\n"))

        def get_total_open_chromatin(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return (peaks.iloc[:, 2] - peaks.iloc[:, 1]).sum()

        def get_peak_lengths(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return (peaks.iloc[:, 2] - peaks.iloc[:, 1])

        def get_peak_chroms(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return peaks.iloc[:, 0].value_counts()

        # Peaks per sample:
        if samples is None:
            samples = self.samples
        stats = pd.DataFrame([
            map(get_sample_reads, [s.filtered for s in samples]),
            map(get_peak_number, [s.peaks for s in samples]),
            map(get_total_open_chromatin, [s.peaks for s in samples])],
            index=["reads_used", "peak_number", "open_chromatin"],
            columns=[s.name for s in samples]).T

        stats["peaks_norm"] = (stats["peak_number"] / stats["reads_used"]) * 1e3
        stats["open_chromatin_norm"] = (stats["open_chromatin"] / stats["reads_used"])
        stats.to_csv(os.path.join(self.results_dir, "{}.open_chromatin_space.csv".format(self.name)), index=True)
        stats = pd.read_csv(os.path.join(self.results_dir, "{}.open_chromatin_space.csv".format(self.name)), index_col=0)

        stats['knockout'] = stats.index.str.extract("HAP1_(.*?)_")
        stats = stats.ix[stats.index[~stats.index.str.contains("GFP|HAP1_ATAC-seq_WT_Bulk_")]]
        # median lengths per sample (split-apply-combine)
        stats = pd.merge(stats, stats.groupby('knockout')['open_chromatin'].median().to_frame(name='group_open_chromatin').reset_index())
        stats = pd.merge(stats, stats.groupby('knockout')['open_chromatin_norm'].median().to_frame(name='group_open_chromatin_norm').reset_index())

        # plot
        stats = stats.sort_values("open_chromatin_norm")
        fig, axis = plt.subplots(2, 1, figsize=(4 * 2, 6 * 1))
        sns.barplot(x="index", y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
        sns.barplot(x="index", y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
        axis[0].set_ylabel("Total open chromatin space (bp)")
        axis[1].set_ylabel("Total open chromatin space (normalized)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=45, ha="right")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.total_open_chromatin_space.per_sample.svg".format(self.name)), bbox_inches="tight")

        # per knockout group
        fig, axis = plt.subplots(2, 1, figsize=(4 * 2, 6 * 1))
        stats = stats.sort_values("group_open_chromatin")
        sns.barplot(x="knockout", y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
        sns.stripplot(x="knockout", y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
        stats = stats.sort_values("group_open_chromatin_norm")
        sns.barplot(x="knockout", y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
        sns.stripplot(x="knockout", y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
        axis[0].axhline(stats.groupby('knockout')['open_chromatin'].median()["WT"], color="black", linestyle="--")
        axis[1].axhline(stats.groupby('knockout')['open_chromatin_norm'].median()["WT"], color="black", linestyle="--")
        axis[0].set_ylabel("Total open chromatin space (bp)")
        axis[1].set_ylabel("Total open chromatin space (normalized)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=45, ha="right")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.total_open_chromatin_space.per_knockout.svg".format(self.name)), bbox_inches="tight")

        # plot distribution of peak lengths
        sample_peak_lengths = map(get_peak_lengths, [s.peaks for s in samples])
        lengths = pd.melt(pd.DataFrame(sample_peak_lengths, index=[s.name for s in samples]).T, value_name='peak_length', var_name="sample_name").dropna()
        lengths['knockout'] = lengths['sample_name'].str.extract("HAP1_(.*?)_")
        lengths = lengths[~lengths['sample_name'].str.contains("GFP|HAP1_ATAC-seq_WT_Bulk_")]
        # median lengths per sample (split-apply-combine)
        lengths = pd.merge(lengths, lengths.groupby('sample_name')['peak_length'].median().to_frame(name='mean_peak_length').reset_index())
        # median lengths per group (split-apply-combine)
        lengths = pd.merge(lengths, lengths.groupby('knockout')['peak_length'].median().to_frame(name='group_mean_peak_length').reset_index())

        lengths = lengths.sort_values("mean_peak_length")
        fig, axis = plt.subplots(2, 1, figsize=(8 * 1, 4 * 2))
        sns.boxplot(x="sample_name", y="peak_length", data=lengths, palette="summer", ax=axis[0], showfliers=False)
        axis[0].set_ylabel("Peak length (bp)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
        sns.boxplot(x="sample_name", y="peak_length", data=lengths, palette="summer", ax=axis[1], showfliers=False)
        axis[1].set_yscale("log")
        axis[1].set_ylabel("Peak length (bp)")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.peak_lengths.per_sample.svg".format(self.name)), bbox_inches="tight")

        lengths = lengths.sort_values("group_mean_peak_length")
        fig, axis = plt.subplots(2, 1, figsize=(8 * 1, 4 * 2))
        sns.boxplot(x="knockout", y="peak_length", data=lengths, palette="summer", ax=axis[0], showfliers=False)
        axis[0].set_ylabel("Peak length (bp)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
        sns.boxplot(x="knockout", y="peak_length", data=lengths, palette="summer", ax=axis[1], showfliers=False)
        axis[1].set_yscale("log")
        axis[1].set_ylabel("Peak length (bp)")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.peak_lengths.per_knockout.svg".format(self.name)), bbox_inches="tight")

        # peaks per chromosome per sample
        chroms = pd.DataFrame(map(get_peak_chroms, [s.peaks for s in samples]), index=[s.name for s in samples]).fillna(0).T
        chroms_norm = (chroms / chroms.sum(axis=0)) * 100
        chroms_norm = chroms_norm.ix[["chr{}".format(i) for i in range(1, 23) + ['X', 'Y', 'M']]]

        fig, axis = plt.subplots(1, 1, figsize=(8 * 1, 8 * 1))
        sns.heatmap(chroms_norm, square=True, cmap="summer", ax=axis)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(self.results_dir, "{}.peak_location.per_sample.svg".format(self.name)), bbox_inches="tight")

        # Peak set across samples:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites if interval.length < 2000], bins=300, kde=False, ax=axis)
        axis.set_xlabel("peak width (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.lengths.svg".format(self.name)), bbox_inches="tight")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis)
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.support.svg".format(self.name)), bbox_inches="tight")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2), sharex=False, sharey=False)
        sns.distplot([x for x in self.closest_tss_distances if x < 1e5], bins=1000, kde=False, hist=True, ax=axis[0][0])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e5], bins=1000, kde=False, hist=True, ax=axis[0][1])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e6], bins=1000, kde=True, hist=False, ax=axis[1][0])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e6], bins=1000, kde=True, hist=False, ax=axis[1][1])
        for ax in axis.flat:
            ax.set_xlabel("distance to nearest TSS (bp)")
            ax.set_ylabel("frequency")
        axis[0][1].set_yscale("log")
        axis[1][1].set_yscale("log")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.tss_distance.svg".format(self.name)), bbox_inches="tight")

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = self.region_annotation['genomic_region'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_region_annotation.name = "foreground"
        all_region_annotation_b = self.region_annotation_b['genomic_region'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_region_annotation_b.name = "background"

        # count region frequency
        data = all_region_annotation.value_counts().sort_values(ascending=False)
        background = all_region_annotation_b.value_counts().sort_values(ascending=False)
        data = data.to_frame().join(background)
        data["fold_change"] = np.log2(data['foreground'] / data['background'])
        data.index.name = "region"

        # plot also % of genome space "used"
        self.region_annotation["length"] = self.region_annotation["end"] - self.region_annotation["start"]

        s = self.region_annotation.join(all_region_annotation).groupby("foreground")["length"].sum()
        s.name = "size"
        s = (s / pd.Series(self.genome_space)).dropna() * 100
        s.name = "percent_space"
        data = data.join(s)

        # plot together
        g = sns.FacetGrid(data=pd.melt(data.reset_index(), id_vars='region'), col="variable", col_wrap=2, sharex=True, sharey=False, size=4)
        g.map(sns.barplot, "region", "value")
        sns.despine(fig)
        g.savefig(os.path.join(self.results_dir, "{}.genomic_regions.svg".format(self.name)), bbox_inches="tight")

        # # Plot chromatin states
        all_chrom_state_annotation = self.chrom_state_annotation['chromatin_state'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_chrom_state_annotation.name = "foreground"
        all_chrom_state_annotation_b = self.chrom_state_annotation_b['chromatin_state'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_chrom_state_annotation_b.name = "background"

        # count region frequency
        data = all_chrom_state_annotation.value_counts().sort_values(ascending=False)
        background = all_chrom_state_annotation_b.value_counts().sort_values(ascending=False)
        data = data.to_frame().join(background)
        data["fold_change"] = np.log2(data['foreground'] / data['background'])
        data.index.name = "region"

        # plot also % of genome space "used"
        self.chrom_state_annotation["length"] = self.chrom_state_annotation["end"] - self.chrom_state_annotation["start"]

        s = self.chrom_state_annotation.join(all_chrom_state_annotation).groupby("foreground")["length"].sum()
        s.name = "size"
        s = (s / pd.Series(self.genome_space)).dropna() * 100
        s.name = "percent_space"
        data = data.join(s)

        # plot together
        g = sns.FacetGrid(data=pd.melt(data.reset_index(), id_vars='region'), col="variable", col_wrap=2, sharex=True, sharey=False, size=4)
        g.map(sns.barplot, "region", "value")
        sns.despine(fig)
        g.savefig(os.path.join(self.results_dir, "{}.chromatin_states.svg".format(self.name)), bbox_inches="tight")

        # distribution of count attributes
        data = self.coverage_annotated.copy()

        fig, axis = plt.subplots(1)
        sns.distplot(data["mean"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.mean.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["qv2"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.qv2.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["dispersion"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.dispersion.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["support"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.support.distplot.svg".format(self.name)), bbox_inches="tight")

        # joint
        for metric in ["support", "variance", "std_deviation", "dispersion", "qv2", "amplitude"]:
            p = data[(data["mean"] > 0) & (data[metric] < np.percentile(data[metric], 99) * 3)]
            g = sns.jointplot(p["mean"], p[metric], s=2, alpha=0.2, rasterized=True)
            sns.despine(g.fig)
            g.fig.savefig(os.path.join(self.results_dir, "{}.mean_{}.svg".format(self.name, metric)), bbox_inches="tight", dpi=300)

    def plot_raw_coverage(self, samples=None, by_attribute=None):
        if samples is None:
            samples = self.samples

        if by_attribute is None:
            cov = pd.melt(
                np.log2(1 + self.coverage[[s.name for s in samples]]),
                var_name="sample_name", value_name="counts"
            )
            fig, axis = plt.subplots(1, 1, figsize=(6, 1 * 4))
            sns.violinplot("sample_name", "counts", data=cov, ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.svg"), bbox_inches="tight")

        else:
            attrs = set([getattr(s, by_attribute) for s in samples])
            fig, axis = plt.subplots(len(attrs), 1, figsize=(8, len(attrs) * 6))
            for i, attr in enumerate(attrs):
                print(attr)
                cov = pd.melt(
                    np.log2(1 + self.coverage[[s.name for s in samples if getattr(s, by_attribute) == attr]]),
                    var_name="sample_name", value_name="counts"
                )
                sns.violinplot("sample_name", "counts", data=cov, ax=axis[i])
                axis[i].set_title(attr)
                axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.by_{}.svg".format(by_attribute)), bbox_inches="tight")

    def plot_coverage(self, samples=None, by_attribute=None):
        if samples is None:
            samples = self.samples

        if by_attribute is None:
            cov = pd.melt(
                np.log2(1 + self.coverage_qnorm[[s.name for s in samples]]),
                var_name="sample_name", value_name="counts"
            )
            fig, axis = plt.subplots(1, 1, figsize=(6, 1 * 4))
            sns.violinplot("sample_name", "counts", data=cov, ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".qnorm_counts.violinplot.svg"), bbox_inches="tight")

        else:
            attrs = set([getattr(s, by_attribute) for s in samples])
            fig, axis = plt.subplots(len(attrs), 1, figsize=(8, len(attrs) * 6))
            for i, attr in enumerate(attrs):
                print(attr)
                cov = pd.melt(
                    np.log2(1 + self.coverage_qnorm[[s.name for s in samples if getattr(s, by_attribute) == attr]]),
                    var_name="sample_name", value_name="counts"
                )
                sns.violinplot("sample_name", "counts", data=cov, ax=axis[i])
                axis[i].set_title(attr)
                axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".qnorm_counts.violinplot.by_{}.svg".format(by_attribute)), bbox_inches="tight")

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

    def plot_expression_characteristics(self):
        fig, axis = plt.subplots(figsize=(4, 4 * np.log10(len(self.expression_matrix_counts.columns))))
        sns.barplot(data=self.expression_matrix_counts.sum().sort_values().reset_index(), y="index", x=0, orient="horiz", color=sns.color_palette("colorblind")[0], ax=axis)
        axis.set_xlabel("Reads")
        axis.set_ylabel("Samples")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "expression.reads_per_sample.svg"), bbox_inches="tight")

        cov = pd.DataFrame()
        for i in [1, 2, 3, 6, 12, 24, 48, 96, 200, 300, 400, 500, 1000]:
            cov[i] = self.expression_matrix_counts.apply(lambda x: sum(x >= i))

        fig, axis = plt.subplots(1, 2, figsize=(6 * 2, 6))
        sns.heatmap(cov.drop([1, 2], axis=1), ax=axis[0], cmap="GnBu", cbar_kws={"label": "Genes covered"})
        sns.heatmap(cov.drop([1, 2], axis=1).apply(lambda x: (x - x.mean()) / x.std(), axis=0), ax=axis[1], cbar_kws={"label": "Z-score"})
        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_ylabel("Samples")
            ax.set_xlabel("Genes with #reads")
        axis[1].set_yticklabels(ax.get_yticklabels(), visible=False)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "expression.genes_with_reads.svg"), bbox_inches="tight")

        for name, matrix in [("counts", self.expression_matrix_counts), ("qnorm_TPM", self.expression)]:
            # Boxplot with values per sample
            to_plot = pd.melt(matrix, var_name="sample")

            fig, axis = plt.subplots()
            sns.boxplot(data=to_plot, y="sample", x="value", orient="horiz", ax=axis)
            axis.set_xlabel(name)
            axis.set_ylabel("Samples")
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, "expression.boxplot_per_sample.{}.svg".format(name)), bbox_inches="tight")

        # Plot gene expression along chromossomes
        import requests
        url_query = "".join([
            """http://grch37.ensembl.org/biomart/martservice?query=""",
            """<?xml version="1.0" encoding="UTF-8"?>""",
            """<!DOCTYPE Query>""",
            """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
            """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
            """<Attribute name = "chromosome_name" />""",
            """<Attribute name = "start_position" />""",
            """<Attribute name = "end_position" />""",
            """<Attribute name = "external_gene_name" />""",
            """</Dataset>""",
            """</Query>"""])
        req = requests.get(url_query, stream=True)
        annot = pd.DataFrame((x.strip().split(",") for x in list(req.iter_lines())), columns=["chr", "start", "end", "gene_name"])
        gene_order = annot[annot['chr'].isin([str(x) for x in range(22)] + ['X', 'Y'])].sort_values(["chr", "start", "end"])["gene_name"]

        exp = self.expression.ix[gene_order].dropna()
        data = exp.rolling(int(3e4), axis=0).mean().dropna()
        plt.pcolor(data)

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

        color_dataframe = pd.DataFrame(
            self.get_level_colors(index=self.accessibility[[s.name for s in samples]].columns, levels=attributes, pallete="Vega20b"),
            index=attributes, columns=[s.name for s in samples])

        # exclude samples if needed
        samples = [s for s in samples if s.name not in exclude and s.library == "ATAC-seq"]
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
            X.corr(), xticklabels=False, yticklabels=sample_display_names, annot=True,
            cmap="YlGnBu", figsize=(15, 15), cbar_kws={"label": "Pearson correlation"}, row_colors=color_dataframe.values.tolist())
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xlabel(None)
        g.ax_heatmap.set_ylabel(None)
        g.fig.savefig(os.path.join(self.results_dir, "{}.all_sites.corr.clustermap.svg".format(self.name)), bbox_inches='tight')

        # MDS
        mds = MDS(n_jobs=-1)
        x_new = mds.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new, index=X.columns, columns=range(x_new.shape[1]))
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        fig, axis = plt.subplots(1, len(to_plot), figsize=(4 * len(to_plot), 4 * 1))
        axis = axis.flatten()
        for i, attr in enumerate(to_plot):
            for j, sample in enumerate(xx.index):
                try:
                    label = getattr(samples[j], to_plot[i])
                except AttributeError:
                    label = np.nan
                axis[i].scatter(xx.loc[sample, 0], xx.loc[sample, 1], s=50, color=color_dataframe.iloc[i, j], label=label)
            axis[i].set_title(to_plot[i])
            axis[i].set_xlabel("MDS 1")
            axis[i].set_ylabel("MDS 2")
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
            axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)

            # Unique legend labels
            handles, labels = axis[i].get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 30:
                # if not any([re.match("^\d", c) for c in by_label.keys()]):
                axis[i].legend(by_label.values(), by_label.keys())

        # Add group text labels to middle
        for i, attr in enumerate(to_plot):
            group_xx = xx.groupby(level=[attr]).mean()
            # average colour across samples of same group
            color_dataframe.columns = xx.index
            group_color_dataframe = color_dataframe.T.groupby(level=[attr])[attr].apply(lambda x: np.array(list(x)).mean(0))
            for j, group in enumerate(group_xx.index):
                axis[i].text(group_xx.loc[group, 0], group_xx.loc[group, 1], s=group, color=group_color_dataframe.ix[j], ha="center")

        fig.savefig(os.path.join(self.results_dir, "{}.all_sites.mds.svg".format(self.name)), bbox_inches="tight")

        # PCA
        pca = PCA()
        x_new = pca.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new, index=X.columns, columns=range(x_new.shape[1]))
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
        pcs = min(xx.shape[0] - 1, 4)
        fig, axis = plt.subplots(pcs, len(to_plot), figsize=(4 * len(to_plot), 4 * pcs))
        for pc in range(pcs):
            for i, attr in enumerate(to_plot):
                for j, sample in enumerate(xx.index):
                    try:
                        label = getattr(samples[j], to_plot[i])
                    except AttributeError:
                        label = np.nan
                    axis[pc, i].scatter(xx.loc[sample, pc], xx.loc[sample, pc + 1], s=50, color=color_dataframe.iloc[i, j], label=label)
                axis[pc, i].set_title(to_plot[i])
                axis[pc, i].set_xlabel("PC {}".format(pc + 1))
                axis[pc, i].set_ylabel("PC {}".format(pc + 2))
                axis[pc, i].set_xticklabels(axis[pc, i].get_xticklabels(), visible=False)
                axis[pc, i].set_yticklabels(axis[pc, i].get_yticklabels(), visible=False)

                # Unique legend labels
                handles, labels = axis[pc, i].get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 30:
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    axis[pc, i].legend(by_label.values(), by_label.keys())

        # Add group text labels to middle
        for pc in range(pcs):
            for i, attr in enumerate(to_plot):
                group_xx = xx.groupby(level=[attr]).mean()
                # average colour across samples of same group
                color_dataframe.columns = xx.index
                group_color_dataframe = color_dataframe.T.groupby(level=[attr])[attr].apply(lambda x: np.array(list(x)).mean(0))
                for j, group in enumerate(group_xx.index):
                    axis[pc, i].text(group_xx.loc[group, pc], group_xx.loc[group, pc + 1], s=group, color=group_color_dataframe.ix[group], ha="center")

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
                    from scipy.stats import pearsonr
                    pearsonr

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
                            writeCombinedEnrichment(lolaResults, outFolder=dirname(bedFile), includeSplits=FALSE)
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
        n_vars = 10
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
        p_values = pd.read_csv(os.path.join("results", "{}.PCA.PC_pvalues.csv".format(self.name)), index_col=0)

        # Get enrichments of each PC-regions
        lola_enrichments = pd.DataFrame()
        enrichr_enrichments = pd.DataFrame()
        for pc in p_values.columns[16:]:
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

            # # read lola
            # lol = pd.read_csv(os.path.join(output_folder, "allEnrichments.txt"), sep="\t")
            # lol["PC"] = pc
            # lola_enrichments = lola_enrichments.append(lol)

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

    def unsupervised_expression(
            self, samples, attributes=["knockout", "replicate", "clone", "complex_subunit"],
            exclude=["RNA-seq_HAP1_WT_r1_GFP", "RNA-seq_HAP1_WT_r2_GFP", "RNA-seq_HAP1_WT_r3_GFP", "RNA-seq_HAP1_WT_r4_GFP"]):
        """
        """
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from collections import OrderedDict
        import re
        import itertools
        from scipy.stats import kruskal
        from scipy.stats import pearsonr
        pearsonr

        samples = [s for s in samples if s.name in self.expression_annotated.columns.get_level_values("sample_name") and s.library == "RNA-seq" and s.cell_line == "HAP1" and s.name not in exclude]

        color_dataframe = pd.DataFrame(
            self.get_level_colors(index=self.expression_annotated[[s.name for s in samples]].columns, levels=attributes, pallete="Vega20b"),
            index=attributes, columns=[s.name for s in samples])

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
        X = self.expression_annotated[[s.name for s in samples]]

        # Pairwise correlations
        g = sns.clustermap(
            X.corr(), xticklabels=False, yticklabels=sample_display_names, annot=True,
            cmap="Spectral_r", figsize=(15, 15), cbar_kws={"label": "Pearson correlation"}, row_colors=color_dataframe.values.tolist())
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        g.ax_heatmap.set_xlabel(None)
        g.ax_heatmap.set_ylabel(None)
        g.fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.corr.clustermap.svg".format(self.name)), bbox_inches='tight')

        # MDS
        mds = MDS(n_jobs=-1)
        x_new = mds.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new, index=X.columns, columns=range(x_new.shape[1]))
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        fig, axis = plt.subplots(1, len(to_plot), figsize=(4 * len(to_plot), 4 * 1))
        axis = axis.flatten()
        for i, attr in enumerate(to_plot):
            for j, sample in enumerate(xx.index):
                try:
                    label = getattr(samples[j], to_plot[i])
                except AttributeError:
                    label = np.nan
                axis[i].scatter(xx.loc[sample, 0], xx.loc[sample, 1], s=50, color=color_dataframe.iloc[i, j], label=label)
            axis[i].set_title(to_plot[i])
            axis[i].set_xlabel("MDS 1")
            axis[i].set_ylabel("MDS 2")
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
            axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)

            # Unique legend labels
            handles, labels = axis[i].get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 30:
                # if not any([re.match("^\d", c) for c in by_label.keys()]):
                axis[i].legend(by_label.values(), by_label.keys())

        # Add group text labels to middle
        for i, attr in enumerate(to_plot):
            group_xx = xx.groupby(level=[attr]).mean()
            # average colour across samples of same group
            color_dataframe.columns = xx.index
            group_color_dataframe = color_dataframe.T.groupby(level=[attr])[attr].apply(lambda x: np.array(list(x)).mean(0))
            for j, group in enumerate(group_xx.index):
                axis[i].text(group_xx.loc[group, 0], group_xx.loc[group, 1], s=group, color=group_color_dataframe.ix[j], ha="center")

        fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.mds.svg".format(self.name)), bbox_inches="tight")

        # PCA
        pca = PCA()
        x_new = pca.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new, index=X.columns, columns=range(x_new.shape[1]))
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
        fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.pca.explained_variance.svg".format(self.name)), bbox_inches='tight')

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
                if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 30:
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    axis[pc, i].legend(by_label.values(), by_label.keys())

        # Add group text labels to middle
        for pc in range(pcs):
            for i, attr in enumerate(to_plot):
                group_xx = xx.groupby(level=[attr]).mean()
                # average colour across samples of same group
                color_dataframe.columns = xx.index
                group_color_dataframe = color_dataframe.T.groupby(level=[attr])[attr].apply(lambda x: np.array(list(x)).mean(0))
                for j, group in enumerate(group_xx.index):
                    axis[pc, i].text(group_xx.loc[group, pc], group_xx.loc[group, pc + 1], s=group, color=group_color_dataframe.ix[group], ha="center")

        fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.pca.svg".format(self.name)), bbox_inches="tight")

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
                    from scipy.stats import pearsonr
                    pearsonr

                    associations.append([pc + 1, attr, variable_type, np.nan, np.nan, p])

        associations = pd.DataFrame(associations, columns=["pc", "attribute", "variable_type", "group_1", "group_2", "p_value"])

        # write
        associations.to_csv(os.path.join(self.results_dir, "{}.expression.all_sites.pca.variable_principle_components_association.csv".format(self.name)), index=False)

        # Plot
        # associations[associations['p_value'] < 0.05].drop(['group_1', 'group_2'], axis=1).drop_duplicates()
        # associations.drop(['group_1', 'group_2'], axis=1).drop_duplicates().pivot(index="pc", columns="attribute", values="p_value")
        pivot = associations.groupby(["pc", "attribute"]).min()['p_value'].reset_index().pivot(index="pc", columns="attribute", values="p_value").dropna(axis=1)

        # heatmap of -log p-values
        g = sns.clustermap(-np.log10(pivot), row_cluster=False, annot=True, cbar_kws={"label": "-log10(p_value) of association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.pca.variable_principle_components_association.svg".format(self.name)), bbox_inches="tight")

        # heatmap of masked significant
        g = sns.clustermap((pivot < 0.05).astype(int), row_cluster=False, cbar_kws={"label": "significant association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.expression.all_sites.pca.variable_principle_components_association.masked.svg".format(self.name)), bbox_inches="tight")

    def differential_region_analysis(
            self, samples, trait="knockout",
            variables=["knockout", "replicate"],
            output_suffix="deseq_knockout"):
        """
        Discover differential regions across samples that are associated with a certain trait.
        """
        sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]

        # Get matrix of counts
        counts_matrix = self.coverage[[s.name for s in sel_samples]]

        # Get experiment matrix
        experiment_matrix = pd.DataFrame([sample.as_series() for sample in sel_samples], index=[sample.name for sample in sel_samples])
        # keep only variables
        experiment_matrix = experiment_matrix[["sample_name"] + variables].fillna("Unknown")

        # Make output dir
        output_dir = os.path.join(self.results_dir, output_suffix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Run DESeq2 analysis
        deseq_table = DESeq_analysis(
            counts_matrix, experiment_matrix, trait, covariates=[x for x in variables if x != trait], output_prefix=os.path.join(output_dir, output_suffix), alpha=0.05)

        # to just read in
        # deseq_table = pd.read_csv(os.path.join(output_dir, output_suffix) + ".%s.csv" % trait, index_col=0)
        # self.coverage_annotated = pd.read_csv(os.path.join(self.results_dir, "breg_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index_col=0)

        df = self.coverage_annotated.join(deseq_table)
        df['comparison'] = df['comparison'].str.replace("-WT", "")
        df.to_csv(os.path.join(output_dir, output_suffix) + ".%s.annotated.csv" % trait)
        df = pd.read_csv(os.path.join(output_dir, output_suffix) + ".%s.annotated.csv" % trait, index_col=0)
        df['comparison'] = df['comparison'].str.replace("-WT", "")
        df = df[~df['comparison'].str.contains("-")]
        df = df[~df['comparison'].str.contains("SMARCC2|GFP")]

        # Extract significant based on p-value and fold-change
        diff = df[(df["padj"] < 0.01) & (abs(df["log2FoldChange"]) > 1.)]

        if diff.shape[0] < 1:
            print("No significantly different regions found.")
            return

        # groups = list(set([getattr(s, trait) for s in sel_samples]))
        comparisons = pd.Series(df['comparison'].unique())
        groups = sorted(comparisons.tolist())

        # Statistics of differential regions
        import string
        total_sites = float(len(self.sites))

        total_diff = diff.groupby(["comparison"])['stat'].count().sort_values(ascending=False)
        fig, axis = plt.subplots(1)
        sns.barplot(total_diff.values, total_diff.index, orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.total.svg" % (output_suffix, trait)), bbox_inches="tight")
        # percentage of total
        fig, axis = plt.subplots(1)
        sns.barplot(
            (total_diff.values / total_sites) * 100,
            total_diff.index,
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.total_percentage.svg" % (output_suffix, trait)), bbox_inches="tight")

        # direction-dependent
        diff["direction"] = diff["log2FoldChange"].apply(lambda x: "up" if x >= 0 else "down")

        split_diff = diff.groupby(["comparison", "direction"])['stat'].count().sort_values(ascending=False)

        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.barplot(
            split_diff.values,
            split_diff.reset_index()[['comparison', 'direction']].apply(string.join, axis=1),
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split.svg" % (output_suffix, trait)), bbox_inches="tight")
        # percentage of total
        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.barplot(
            (split_diff.values / total_sites) * 100,
            split_diff.reset_index()[['comparison', 'direction']].apply(string.join, axis=1),
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split_percentage.svg" % (output_suffix, trait)), bbox_inches="tight")

        # in same
        split_diff = diff.groupby(["comparison", "direction"])['stat'].count().sort_values(ascending=False)

        split_diff = split_diff.reset_index()
        split_diff.loc[split_diff['direction'] == "down", "stat"] = -split_diff.loc[split_diff['direction'] == "down", "stat"]
        fig, axis = plt.subplots(1, figsize=(6, 3))
        sns.barplot(
            x="comparison", y="stat", data=split_diff, hue="direction",  # split=False,
            order=split_diff.groupby(['comparison'])['stat'].apply(lambda x: sum(abs(x))).sort_values(ascending=False).index,
            orient="v", ax=axis)
        axis.axhline(0, color="black", linestyle="--")
        axis.set_xticklabels(axis.get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split.pos_neg.svg" % (output_suffix, trait)), bbox_inches="tight")

        # percentage of total
        split_diff['stat'] = (split_diff['stat'].values / total_sites) * 100
        fig, axis = plt.subplots(1, figsize=(6, 3))
        sns.barplot(
            x="comparison", y="stat", data=split_diff, hue="direction",  # split=False,
            order=split_diff.groupby(['comparison'])['stat'].apply(lambda x: sum(abs(x))).sort_values(ascending=False).index,
            orient="v", ax=axis)
        axis.axhline(0, color="black", linestyle="--")
        axis.set_xticklabels(axis.get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split_percentage.pos_neg.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Overlap between sets
        # # Pyupset
        # import pyupset as pyu
        # # Build dict
        # diff["comparison_direction"] = diff[["comparison", "direction"]].apply(string.join, axis=1)
        # df_dict = {group: diff[diff["comparison_direction"] == group].reset_index()[['index']] for group in set(diff["comparison_direction"])}
        # # Plot
        # plot = pyu.plot(df_dict, unique_keys=['index'], inters_size_bounds=(10, np.inf))
        # plot['figure'].set_size_inches(20, 8)
        # plot['figure'].savefig(os.path.join(output_dir, "%s.%s.number_differential.upset.svg" % (output_suffix, trait)), bbox_inched="tight")
        diff["intersect"] = 1
        diff["direction"] = diff["log2FoldChange"].apply(lambda x: "up" if x >= 0 else "down")
        piv = pd.pivot_table(diff.reset_index(), index='index', columns=['comparison', 'direction'], values='intersect', fill_value=0)

        import itertools
        intersections = pd.DataFrame(columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
        for ((k1, dir1), i1), ((k2, dir2), i2) in itertools.permutations(piv.T.groupby(level=['comparison', 'direction']).groups.items(), 2):
            print(k1, k2)
            i1 = set(piv[i1][piv[i1] == 1].dropna().index)
            i2 = set(piv[i2][piv[i2] == 1].dropna().index)
            intersections = intersections.append(
                pd.Series(
                    [k1, k2, dir1, dir2, len(i1), len(i2), len(i1.intersection(i2)), len(i1.union(i2))],
                    index=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"]
                ),
                ignore_index=True
            )
        # convert to %
        intersections['intersection_perc'] = (intersections['intersection'] / intersections['size2']) * 100.

        # make pivot tables
        piv_up = pd.pivot_table(intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "up")], index="group1", columns="group2", values="intersection_perc").fillna(0)
        piv_down = pd.pivot_table(intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "down")], index="group1", columns="group2", values="intersection_perc").fillna(0)
        np.fill_diagonal(piv_up.values, np.nan)
        np.fill_diagonal(piv_down.values, np.nan)

        # heatmaps
        fig, axis = plt.subplots(1, 2, figsize=(8 * 2, 8))
        sns.heatmap(piv_down, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis[0])
        sns.heatmap(piv_up, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis[1])
        axis[0].set_title("Downregulated regions")
        axis[1].set_title("Upregulated regions")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="right")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="right")
        axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
        axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.up_down_split.svg" % (output_suffix, trait)), bbox_inches="tight")

        # combined heatmap
        # with upregulated regions in upper square matrix and downredulated in down square
        piv_combined = pd.DataFrame(np.triu(piv_up), index=piv_up.index, columns=piv_up.columns).replace(0, np.nan)
        piv_combined.update(pd.DataFrame(np.tril(piv_down), index=piv_down.index, columns=piv_down.columns).replace(0, np.nan))
        piv_combined = piv_combined.fillna(0)
        np.fill_diagonal(piv_combined.values, np.nan)

        fig, axis = plt.subplots(1, figsize=(8, 8))
        sns.heatmap(piv_combined, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.up_down_together.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Observe disagreement between knockouts
        # (overlap of down-regulated with up-regulated and vice-versa)
        piv_up = pd.pivot_table(intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "down")], index="group1", columns="group2", values="intersection")
        piv_down = pd.pivot_table(intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "up")], index="group1", columns="group2", values="intersection")

        piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()

        fig, axis = plt.subplots(1, 2, figsize=(16, 8))
        sns.heatmap(piv_disagree, square=True, cmap="Greens", cbar_kws={"label": "Discordant regions"}, ax=axis[0])
        sns.heatmap(np.log2(1 + piv_disagree), square=True, cmap="Greens", cbar_kws={"label": "Discordant regions (log2)"}, ax=axis[1])

        norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
        cmap = plt.get_cmap("Greens")
        log_norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log2(1 + piv_disagree).max().max())
        for i, g1 in enumerate(piv_disagree.columns):
            for j, g2 in enumerate(piv_disagree.index):
                axis[0].scatter(
                    len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                    s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1, color=cmap(norm(piv_disagree.loc[g1, g2])), marker="o")
                axis[1].scatter(
                    len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                    s=(100 ** (log_norm(np.log2(1 + piv_disagree).loc[g1, g2]))) - 1, color=cmap(log_norm(np.log2(1 + piv_disagree).loc[g1, g2])), marker="o")

        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.disagreement.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Changes per chromosome per group
        diff['chrom'] = [x[0] for x in diff.index.str.split(":")]
        chroms = diff.groupby(["comparison"])['chrom'].value_counts()
        chroms.name = "count"
        pivot_chroms = pd.pivot_table(chroms.reset_index(), index="chrom", columns='comparison', values="count", fill_value=0)
        chroms_norm = (pivot_chroms / pivot_chroms.sum(axis=0)) * 100
        chroms_norm = chroms_norm.ix[["chr{}".format(i) for i in range(1, 23) + ['X']]]

        fig, axis = plt.subplots(1, 1, figsize=(8 * 1, 8 * 1))
        sns.heatmap(chroms_norm, square=True, cmap="summer", ax=axis)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.chrom_location.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Pairwise scatter plots
        n_rows = n_cols = int(np.ceil(np.sqrt(len(comparisons))))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            if "-" not in cond1:
                cond2 = "WT"
            else:
                cond1, cond2 = cond1.split("-")

            # Hexbin plot
            axis.hexbin(np.log2(1 + df2[cond1]), np.log2(1 + df2[cond2]), alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)
            axis.set_xlabel(cond1)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # Scatter plot
                axis.scatter(np.log2(1 + diff2[cond1]), np.log2(1 + diff2[cond2]), alpha=0.1, color="red", s=2, rasterized=True)
            m = max(np.log2(1 + df2[cond1]).max(), np.log2(1 + df2[cond2]).max())
            axis.plot([0, m], [0, m], color="black", alpha=0.8, linestyle="--")
            axis.set_ylabel(cond1)
            axis.set_ylabel(cond2)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.scatter_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Volcano plots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            if "-" not in cond1:
                cond2 = "WT"
            else:
                cond1, cond2 = cond1.split("-")
            # hexbin
            axis.hexbin(df2["log2FoldChange"], -np.log10(df2['pvalue']), alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # significant scatter
                axis.scatter(diff2["log2FoldChange"], -np.log10(diff2['pvalue']), alpha=0.2, color="red", s=2, rasterized=True)
            axis.axvline(0, linestyle="--", color="black", alpha=0.8)
            axis.set_title(comparison)
            axis.set_xlabel("log2(fold change)")
            axis.set_ylabel("-log10(p-value)")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.volcano_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # MA plots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            if "-" not in cond1:
                cond2 = "WT"
            else:
                cond1, cond2 = cond1.split("-")
            # hexbin
            axis.hexbin(np.log2(df2["baseMean"]), df2["log2FoldChange"], alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # significant scatter
                axis.scatter(np.log2(diff2["baseMean"]), diff2["log2FoldChange"], alpha=0.2, color="red", s=2, rasterized=True)
            axis.axhline(0, linestyle="--", color="black", alpha=0.8)
            axis.set_title(comparison)
            axis.set_xlabel("Intensity")
            axis.set_ylabel("log2(fold change)")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.ma_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Pairwise scatter plots of differential genes
        import itertools
        new_groups = [g for g in groups if "ARID" in g or "SMARC" in g]
        n_rows = n_cols = int(np.ceil(np.sqrt(len(list(itertools.combinations(new_groups, 2))))))

        # with accessibility values
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        axes = iter(axes.flatten())
        for cond1, cond2 in sorted(itertools.combinations(new_groups, 2)):
            diff1 = set(diff[diff["comparison"] == cond1].index)
            diff2 = set(diff[diff["comparison"] == cond2].index)

            if (len(diff1) == 0) or (len(diff2) == 0):
                continue
            print(cond1, cond2)
            axis = axes.next()

            # get union of differential regions
            inter = diff1.intersection(diff2)
            exc1 = diff1.difference(diff2)
            exc2 = diff2.difference(diff1)

            dots = list(exc1) + list(exc2)
            colors = ["red" for _ in exc1] + ["green" for _ in exc2]
            d = pd.DataFrame([dots, colors], index=["region", "color"]).T.sample(frac=1)
            df1 = df[df["comparison"] == cond1].ix[d['region']]
            df2 = df[df["comparison"] == cond2].ix[d['region']]
            axis.scatter(np.log2(1 + df1[cond1]), np.log2(1 + df2[cond2]), alpha=0.5, color=d['color'], s=5, rasterized=True)
            df1 = df[df["comparison"] == cond1].ix[inter]
            df2 = df[df["comparison"] == cond2].ix[inter]
            axis.scatter(np.log2(1 + df1[cond1]), np.log2(1 + df2[cond2]), alpha=0.85, color="orange", s=15, rasterized=True)
            axis.set_xlabel(cond1)
            axis.set_ylabel(cond2)

            # add X=Y line
            m = max(np.log2(1 + df2[cond1]).max(), np.log2(1 + df2[cond2]).max())
            axis.plot([0, m], [0, m], color="black", alpha=0.8, linestyle="--")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.group_combinations.scatter_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # with fold-change values
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        axes = iter(axes.flatten())
        for cond1, cond2 in sorted(itertools.combinations(new_groups, 2)):
            diff1 = set(diff[diff["comparison"] == cond1].index)
            diff2 = set(diff[diff["comparison"] == cond2].index)

            if (len(diff1) == 0) or (len(diff2) == 0):
                continue
            print(cond1, cond2)
            axis = axes.next()

            # get union of differential regions
            inter = diff1.intersection(diff2)
            exc1 = diff1.difference(diff2)
            exc2 = diff2.difference(diff1)

            dots = list(exc1) + list(exc2)
            colors = ["red" for _ in exc1] + ["green" for _ in exc2]
            d = pd.DataFrame([dots, colors], index=["region", "color"]).T.sample(frac=1)
            df1 = df[df["comparison"] == cond1].ix[d['region']]
            df2 = df[df["comparison"] == cond2].ix[d['region']]
            axis.scatter(df1["log2FoldChange"], df2["log2FoldChange"], alpha=0.5, color=d['color'], s=5, rasterized=True)
            df1 = df[df["comparison"] == cond1].ix[inter]
            df2 = df[df["comparison"] == cond2].ix[inter]
            axis.scatter(df1["log2FoldChange"], df2["log2FoldChange"], alpha=0.85, color="orange", s=15, rasterized=True)
            axis.set_xlabel(cond1)
            axis.set_ylabel(cond2)

            # add 0 lines
            axis.axhline(0, color="black", alpha=0.8, linestyle="--")
            axis.axvline(0, color="black", alpha=0.8, linestyle="--")

            # add text with correlation
            _all = d['region'].tolist() + list(inter)
            from scipy.stats import pearsonr
            r, p = pearsonr(df[df["comparison"] == cond1].ix[_all]["log2FoldChange"], df[df["comparison"] == cond2].ix[_all]["log2FoldChange"])
            axis.text(0.1, -3, "Pearson's r:{}".format(r))

        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.group_combinations.scatter_plots.fold_change.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # save unique differential regions
        groups = diff2.columns.tolist()
        groups = [g for g in groups if "-" not in g]
        groups = [g for g in groups if "SMARCC2" not in g]
        groups = [g for g in groups if "GFP" not in g]
        diff2 = diff[groups].ix[diff.index.unique()].drop_duplicates()
        diff2.to_csv(os.path.join(output_dir, "%s.%s.differential_regions.csv" % (output_suffix, trait)))
        diff2 = pd.read_csv(os.path.join(output_dir, "%s.%s.differential_regions.csv" % (output_suffix, trait)), index_col=0)

        # Exploration of differential regions
        # get unique differential regions
        df2 = diff2.join(self.coverage_annotated)

        # Characterize regions
        prefix = "%s.%s.diff_regions" % (output_suffix, trait)
        # region's structure
        # characterize_regions_structure(df=df2, prefix=prefix, output_dir=output_dir)
        # region's function
        # characterize_regions_function(df=df2, prefix=prefix, output_dir=output_dir)

        # Heatmaps
        # Comparison level
        g = sns.clustermap(np.log2(1 + df2[groups]).corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non differential regions"}, cmap="Spectral_r", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(np.log2(1 + df2[groups]), yticklabels=False, cbar_kws={"label": "Accessibility of\ndifferential regions"}, cmap="BuGn", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(
            np.log2(1 + df2[groups]), yticklabels=False, z_score=0, metric="correlation",
            cbar_kws={"label": "Z-score of accessibility\non differential regions"}, vmin=-3, vmax=3, rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.clustermap.z0.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # # "Cluster" regions based on patttern
        # ll = g.z_score(np.log2(1 + df2[groups]), axis=0)
        # arid1b = ll.loc[:, "ARID1B"]
        # arid1a = ll.loc[:, ["ARID1A", "SMARCA4", "SMARCC1"]].mean(1)
        # smarca4 = ll.loc[:, ["SMARCA4"]].mean(1)
        # left = ll.loc[:, ["BCL11A", "BCL11B", "SMARCD3"]].mean(1)

        # ll['cluster'] = 3
        # ll.loc[(left > 0) & (arid1a < -1), "cluster"] = 1
        # ll.loc[(arid1b > 2) & (smarca4 < 0.5) & (left < 0.5), "cluster"] = 0
        # ll.loc[(left < 0) & (smarca4 > 1.5), "cluster"] = 2

        # g2 = sns.clustermap(
        #     ll, yticklabels=False, metric="correlation",
        #     row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
        #     cbar_kws={"label": "Z-score of accessibility\non differential regions"}, vmin=-3, vmax=3,
        #     row_colors=[(sns.color_palette("colorblind")[:3] + [(1., 1., 1.)])[i] for i in ll['cluster']])
        # g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), rotation=90)
        # g2.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.clustermap.z0.kmeans_labeled.png" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Fold-changes and P-values
        # pivot table of regions vs comparisons
        sig = diff.index.drop_duplicates()
        fold_changes = pd.pivot_table(df.reset_index(), index="index", columns="comparison", values="log2FoldChange")
        p_values = -np.log10(pd.pivot_table(df.reset_index(), index="index", columns="comparison", values="padj"))

        # fold
        g = sns.clustermap(fold_changes.corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non fold-changes"}, cmap="Spectral_r", vmin=0, vmax=1, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.fold_changes.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(fold_changes.ix[sig], metric="correlation", yticklabels=False, cbar_kws={"label": "Fold-changes of\ndifferential regions"}, vmin=-3, vmax=3, rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.fold_changes.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # p-values
        g = sns.clustermap(p_values.corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non p-values"}, cmap="Spectral_r", vmin=0, vmax=1, rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.p_values.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(p_values.ix[sig], yticklabels=False, cbar_kws={"label": "-log10(p-value) of\ndifferential regions"}, vmin=0, vmax=20, cmap="Spectral_r", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.groups.p_values.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Sample level
        g = sns.clustermap(df2[[s.name for s in sel_samples]].corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non differential regions"}, cmap="Spectral_r", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.samples.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(df2[[s.name for s in sel_samples]], yticklabels=False, cbar_kws={"label": "Accessibility of\ndifferential regions"}, cmap="BuGn", vmin=0, rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.samples.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        sel_samples = [s for s in sel_samples if s.knockout not in ["SMARCC2"]]
        sel_samples = [s for s in sel_samples if "GFP" not in s.name]
        g = sns.clustermap(df2[[s.name for s in sel_samples]], yticklabels=False, z_score=0, metric="correlation", cbar_kws={"label": "Z-score of accessibility\non differential regions"}, rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.samples.clustermap.z0.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        from scipy.cluster.hierarchy import fcluster
        clusters = fcluster(g.dendrogram_row.linkage, 6, criterion="maxclust")
        colors = sns.color_palette("colorblind") + sns.color_palette("Set1")
        g2 = sns.clustermap(
            df2[[s.name for s in sel_samples]], yticklabels=False,
            row_linkage=g.dendrogram_row.linkage, col_linkage=g.dendrogram_col.linkage,
            z_score=0, metric="correlation", cbar_kws={"label": "Z-score of accessibility\non differential regions"},
            row_colors=[colors[i] for i in clusters], rasterized=True)
        # g2.ax_col_dendrogram.set_rasterized(True)
        g2.ax_row_dendrogram.set_rasterized(True)
        g2.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g2.fig.savefig(os.path.join(output_dir, "%s.%s.diff_regions.samples.clustermap.z0.clusters_labeled.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Examine each region cluster
        region_enr = pd.DataFrame()
        lola_enr = pd.DataFrame()
        motif_enr = pd.DataFrame()
        homer_enr = pd.DataFrame()
        pathway_enr = pd.DataFrame()
        for cluster, cluster_name in [(1, "ARID_SMARC_down"), (5, "SMARCA4_up"), (6, "ARID2_up")]:
            # Separate in up/down-regulated regions
            comparison_df = self.coverage_annotated.ix[df2.index[clusters == cluster]]

            # Characterize regions
            prefix = "%s.%s.diff_regions.%s" % (output_suffix, trait, cluster_name)

            comparison_dir = os.path.join(output_dir, prefix)

            print("Doing regions of cluster %s, with prefix %s" % (cluster_name, prefix))

            # region's structure
            if not os.path.exists(os.path.join(comparison_dir, prefix + "_regions.region_enrichment.csv")):
                print(prefix)
                characterize_regions_structure(df=comparison_df, prefix=prefix, output_dir=comparison_dir)
            # region's function
            if not os.path.exists(os.path.join(comparison_dir, prefix + "_genes.enrichr.csv")):
                print(prefix)
                characterize_regions_function(df=comparison_df, prefix=prefix, output_dir=comparison_dir)

            # Read/parse enrichment outputs and add to DFs
            enr = pd.read_csv(os.path.join(comparison_dir, prefix + "_regions.region_enrichment.csv"))
            enr["comparison"] = prefix
            region_enr = region_enr.append(enr, ignore_index=True)

            enr = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.txt"), sep="\t")
            enr["comparison"] = prefix
            lola_enr = lola_enr.append(enr, ignore_index=True)

            enr = parse_ame(comparison_dir).reset_index()
            enr["comparison"] = prefix
            motif_enr = motif_enr.append(enr, ignore_index=True)

            enr = parse_homer(os.path.join(comparison_dir, "homerResults.html"))
            enr["comparison"] = prefix
            homer_enr = homer_enr.append(enr, ignore_index=True)

            enr = pd.read_csv(os.path.join(comparison_dir, prefix + "_genes.enrichr.csv"))
            enr["comparison"] = prefix
            pathway_enr = pathway_enr.append(enr, ignore_index=True)

        # Examine each region cluster
        max_diff = 1000
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()

            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            # Separate in up/down-regulated regions
            for f, direction, top in [(np.less, "down", "head"), (np.greater, "up", "tail")]:
                comparison_df = self.coverage_annotated.ix[diff[
                    (diff["comparison"] == comparison) &
                    (f(diff["log2FoldChange"], 0))
                ].index]

                # Handle extremes of regions
                if comparison_df.shape[0] < 1:
                    continue

                if comparison_df.shape[0] > max_diff:
                    comparison_df = comparison_df.ix[
                        getattr(
                            diff[
                                (diff["comparison"] == comparison) &
                                (f(diff["log2FoldChange"], 0))]
                            ["log2FoldChange"].sort_values(), top)
                        (max_diff).index]

                # Characterize regions
                prefix = "%s.%s.diff_regions.comparison_%s.%s" % (output_suffix, trait, comparison, direction)

                comparison_dir = os.path.join(output_dir, prefix)

                print("Doing regions of comparison %s, with prefix %s" % (comparison, prefix))

                # region's structure
                if not os.path.exists(os.path.join(comparison_dir, prefix + "_regions.region_enrichment.csv")):
                    print(prefix)
                    characterize_regions_structure(df=comparison_df, prefix=prefix, output_dir=comparison_dir)
                # region's function
                if not os.path.exists(os.path.join(comparison_dir, prefix + "_genes.enrichr.csv")):
                    print(prefix)
                    characterize_regions_function(df=comparison_df, prefix=prefix, output_dir=comparison_dir)

                # Read/parse enrichment outputs and add to DFs
                enr = pd.read_csv(os.path.join(comparison_dir, prefix + "_regions.region_enrichment.csv"))
                enr["comparison"] = prefix
                region_enr = region_enr.append(enr, ignore_index=True)

                enr = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.txt"), sep="\t")
                enr["comparison"] = prefix
                lola_enr = lola_enr.append(enr, ignore_index=True)

                enr = parse_ame(comparison_dir).reset_index()
                enr["comparison"] = prefix
                motif_enr = motif_enr.append(enr, ignore_index=True)

                enr = parse_homer(os.path.join(comparison_dir, "homerResults.html"))
                enr["comparison"] = prefix
                homer_enr = homer_enr.append(enr, ignore_index=True)

                enr = pd.read_csv(os.path.join(comparison_dir, prefix + "_genes.enrichr.csv"))
                enr["comparison"] = prefix
                pathway_enr = pathway_enr.append(enr, ignore_index=True)

        # write combined enrichments
        region_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_regions.regions.csv" % (output_suffix, trait)), index=False)
        lola_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_regions.lola.csv" % (output_suffix, trait)), index=False)
        pathway_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_regions.enrichr.csv" % (output_suffix, trait)), index=False)
        motif_enr.columns = ["id", "p_value", "comparison"]
        motif_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_regions.motifs.csv" % (output_suffix, trait)), index=False)
        homer_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_regions.homer.csv" % (output_suffix, trait)), index=False)

    def investigate_differential_regions(
            self, samples, trait="knockout",
            variables=["knockout", "replicate"],
            output_suffix="deseq_knockout", n=20, method="groups"):
        import string
        from scipy.cluster.hierarchy import fcluster

        output_dir = os.path.join(self.results_dir, output_suffix)

        # REGION TYPES
        # read in
        regions = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_regions.regions.csv" % (output_suffix, trait)))
        # pretty names
        regions["comparison"] = regions["comparison"].str.replace("%s.%s.diff_regions.comparison_" % (output_suffix, trait), "")
        regions["comparison"] = regions["comparison"].str.replace("%s.%s.diff_regions." % (output_suffix, trait), "")
        regions = regions[~regions["comparison"].str.contains("_|SMARCC2")]

        for metric in ["test_set", "%_df_space", "fold_change", "fold_size_change"]:
            # pivot table
            regions_pivot = pd.pivot_table(regions, values=metric, columns="region", index="comparison", fill_value=0)

            # plot correlation
            fig = sns.clustermap(regions_pivot)
            for tick in fig.ax_heatmap.get_xticklabels():
                tick.set_rotation(90)
            for tick in fig.ax_heatmap.get_yticklabels():
                tick.set_rotation(0)
            fig.savefig(os.path.join(output_dir, "region_type_enrichment.{}.svg".format(metric)), bbox_inches="tight")

        #

        # LOLA
        # read in
        lola = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_regions.lola.csv" % (output_suffix, trait)))
        # pretty names
        lola["comparison"] = lola["comparison"].str.replace("%s.%s.diff_regions.comparison_" % (output_suffix, trait), "")
        lola["comparison"] = lola["comparison"].str.replace("%s.%s.diff_regions." % (output_suffix, trait), "")
        lola = lola[~lola["comparison"].str.contains("_|SMARCC2")]

        # unique ids for lola sets
        cols = ['description', u'cellType', u'tissue', u'antibody', u'treatment', u'dataSource', u'filename']
        lola['label'] = lola[cols].astype(str).apply(string.join, axis=1)
        lola["label"] = lola["label"].str.decode("utf-8")
        lola["label"] = (
            lola["label"]
            .str.replace("bock_regions_mm10", "")
            .str.replace(" nan", "")
            .str.replace(".bed", ""))

        # pivot table
        lola_pivot = pd.pivot_table(lola, values="pValueLog", columns="label", index="comparison")
        # lola_pivot.columns = lola_pivot.columns.str.decode("utf-8")

        # plot correlation
        g = sns.clustermap(lola_pivot.T.corr(), cbar_kws={"label": "Correlation on region enrichemnts\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "lola.correlation.svg"), bbox_inches="tight")
        g.fig.savefig(os.path.join(output_dir, "lola.correlation.png"), bbox_inches="tight", dpi=300)

        #

        if method == "clustering":
            # Get top n terms which are more in each cluster compared with all others
            cluster_assignment = fcluster(fig.dendrogram_col.linkage, 3, criterion="maxclust")

            top_terms = list()
            cluster_means = pd.DataFrame()
            for cluster in set(cluster_assignment):
                cluster_comparisons = lola_pivot.index[cluster_assignment == cluster].tolist()
                other_comparisons = lola_pivot.index[cluster_assignment != cluster].tolist()

                terms = (lola_pivot.ix[cluster_comparisons].mean() - lola_pivot.ix[other_comparisons].mean()).sort_values()

                top_terms += terms.dropna().head(n).index.tolist()

                # additionallly, get mean of cluster
                cluster_means[cluster] = lola_pivot.ix[cluster_comparisons].mean()
        elif method == "groups":
            top = lola.set_index('label').groupby("comparison")['pValueLog'].nlargest(n)
            top_terms = top.index.get_level_values('label').unique()
            top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

        # plot clustered heatmap
        g = sns.clustermap(lola_pivot[list(set(top_terms))], square=True,
            figsize=(20, 12), cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation", vmax=100)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "lola.cluster_specific.not_rast.svg"), bbox_inches="tight")

        # plot clustered heatmap
        g = sns.clustermap(lola_pivot[list(set(top_terms))],
            figsize=(20, 12), z_score=1, cbar_kws={"label": "Z-score of p-values of enrichment\nof differential regions"}, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "lola.cluster_specific.z_score.svg"), bbox_inches="tight")

        #

        # MOTIFS
        # read in
        motifs = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_regions.motifs.csv" % (output_suffix, trait)))
        motifs.columns = ["motif", "p_value", "comparison"]
        # pretty names
        motifs["comparison"] = motifs["comparison"].str.replace("%s.%s.diff_regions.comparison_" % (output_suffix, trait), "")
        motifs["comparison"] = motifs["comparison"].str.replace("%s.%s.diff_regions." % (output_suffix, trait), "")
        motifs = motifs[~motifs["comparison"].str.contains("_|SMARCC2")]

        # pivot table
        motifs_pivot = pd.pivot_table(motifs, values="p_value", columns="motif", index="comparison")

        # transform p-values
        motifs_pivot = -np.log10(motifs_pivot.fillna(1))
        motifs_pivot = motifs_pivot.replace({np.inf: 300})

        if method == "clustering":
            # Get top n terms which are more in each cluster compared with all others
            cluster_assignment = fcluster(fig.dendrogram_col.linkage, 3, criterion="maxclust")

            top_terms = list()
            cluster_means = pd.DataFrame()
            for cluster in set(cluster_assignment):
                cluster_comparisons = motifs_pivot.index[cluster_assignment == cluster].tolist()
                other_comparisons = motifs_pivot.index[cluster_assignment != cluster].tolist()

                terms = (motifs_pivot.ix[cluster_comparisons].mean() - motifs_pivot.ix[other_comparisons].mean()).sort_values()

                top_terms += terms.dropna().head(n).index.tolist()

                # additionallly, get mean of cluster
                cluster_means[cluster] = motifs_pivot.ix[cluster_comparisons].mean()
        elif method == "groups":
            top = motifs.set_index('motif').groupby("comparison")['p_value'].nlargest(n)
            top_terms = top.index.get_level_values('motif').unique()
            top_terms = top_terms[top_terms.isin(motifs_pivot.columns[motifs_pivot.sum() > 5])]

        # plot correlation
        g = sns.clustermap(motifs_pivot.T.corr(), cbar_kws={"label": "Correlation of motif enrichemnt\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.corr.svg"), bbox_inches="tight")
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.corr.png"), bbox_inches="tight", dpi=300)

        # plot sorted heatmap
        p = motifs_pivot[motifs_pivot.index.str.contains("ACTB|ARID1|SMARCA4|SMARCC1|BCL7B")]
        g = sns.clustermap(p[top_terms], figsize=(32, 10), cbar_kws={"label": "-log10(p-value) of motif enrichment\nof differential regions"}, square=True, row_cluster=False)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.sorted.svg"), bbox_inches="tight")
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.png"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        g = sns.clustermap(motifs_pivot[top_terms], figsize=(32, 10), cbar_kws={"label": "-log10(p-value) of motif enrichment\nof differential regions"}, square=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.svg"), bbox_inches="tight")
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.png"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        g = sns.clustermap(motifs_pivot[top_terms], figsize=(32, 10), z_score=1, cbar_kws={"label": "Z-score of p-values of motif enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.z_score.svg"), bbox_inches="tight")
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.z_score.png"), bbox_inches="tight", dpi=300)

        # HOMER
        # read in
        motifs = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_regions.homer.csv" % (output_suffix, trait)))
        # pretty names
        motifs["comparison"] = motifs["comparison"].str.replace("%s.%s.diff_regions.comparison_" % (output_suffix, trait), "")
        motifs["comparison"] = motifs["comparison"].str.replace("%s.%s.diff_regions." % (output_suffix, trait), "")
        motifs = motifs[~motifs["comparison"].str.contains("_|SMARCC2")]

        # try joining them by family
        motifs['family'] = motifs['Best Match/Details'].str.extract("^(.*?)[:/\(\.].*$", expand=False)

        # pivot table
        motifs_pivot = pd.pivot_table(motifs, values="-log_p_value", columns="motif_name", index="comparison", fill_value=1)

        if method == "clustering":
            # Get top n terms which are more in each cluster compared with all others
            cluster_assignment = fcluster(fig.dendrogram_col.linkage, 3, criterion="maxclust")

            top_terms = list()
            cluster_means = pd.DataFrame()
            for cluster in set(cluster_assignment):
                cluster_comparisons = motifs_pivot.index[cluster_assignment == cluster].tolist()
                other_comparisons = motifs_pivot.index[cluster_assignment != cluster].tolist()

                terms = (motifs_pivot.ix[cluster_comparisons].mean() - motifs_pivot.ix[other_comparisons].mean()).sort_values()

                top_terms += terms.dropna().head(n).index.tolist()

                # additionallly, get mean of cluster
                cluster_means[cluster] = motifs_pivot.ix[cluster_comparisons].mean()
        elif method == "groups":
            top = motifs.set_index('motif_name').groupby("comparison")['-log_p_value'].nlargest(n)
            top_terms = top.index.get_level_values('motif_name').unique()
            top_terms = top_terms[top_terms.isin(motifs_pivot.columns[motifs_pivot.sum() > 5])]

        # plot correlation
        g = sns.clustermap(motifs_pivot.T.corr(), cbar_kws={"label": "Correlation of motif enrichemnt\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.corr.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        g = sns.clustermap(motifs_pivot[top_terms], figsize=(32, 10), cbar_kws={"label": "-log10(p-value) of motif enrichment\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        g = sns.clustermap(motifs_pivot[top_terms], figsize=(32, 10), z_score=1, cbar_kws={"label": "Z-score of p-values of motif enrichment\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "motifs.cluster_specific.z_score.svg"), bbox_inches="tight", dpi=300)

        #

        # ENRICHR
        # read in
        enrichr = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_regions.enrichr.csv" % (output_suffix, trait)))
        # pretty names
        enrichr["comparison"] = enrichr["comparison"].str.replace("%s.%s.diff_regions.comparison_" % (output_suffix, trait), "")
        enrichr["comparison"] = enrichr["comparison"].str.replace("%s.%s.diff_regions." % (output_suffix, trait), "")
        enrichr = enrichr[~enrichr["comparison"].str.contains("_|SMARCC2")]
        enrichr["description"] = enrichr["description"].str.decode("utf-8")

        for gene_set_library in enrichr["gene_set_library"].unique():
            print(gene_set_library)
            if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
                continue

            # pivot table
            enrichr_pivot = pd.pivot_table(
                enrichr[enrichr["gene_set_library"] == gene_set_library],
                values="p_value", columns="description", index="comparison", fill_value=1)

            # transform p-values
            enrichr_pivot = -np.log10(enrichr_pivot)
            enrichr_pivot = enrichr_pivot.replace({np.inf: 300})

            # plot correlation
            g = sns.clustermap(enrichr_pivot.T.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential regions"}, rasterized=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.correlation.svg" % gene_set_library), bbox_inches="tight", dpi=300)

            if method == "clustering":
                # Get top n terms which are more in each cluster compared with all others
                cluster_assignment = fcluster(fig.dendrogram_col.linkage, 3, criterion="maxclust")

                top_terms = list()
                cluster_means = pd.DataFrame()
                for cluster in set(cluster_assignment):
                    cluster_comparisons = lola_pivot.index[cluster_assignment == cluster].tolist()
                    other_comparisons = lola_pivot.index[cluster_assignment != cluster].tolist()

                    terms = (lola_pivot.ix[cluster_comparisons].mean() - lola_pivot.ix[other_comparisons].mean()).sort_values()

                    top_terms += terms.dropna().head(n).index.tolist()

                    # additionallly, get mean of cluster
                    cluster_means[cluster] = lola_pivot.ix[cluster_comparisons].mean()
            elif method == "groups":
                top = enrichr[enrichr["gene_set_library"] == gene_set_library].set_index('description').groupby("comparison")['p_value'].nlargest(n)
                top_terms = top.index.get_level_values('description').unique()
                # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))], figsize=(20, 12), cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, rasterized=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.svg" % gene_set_library), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))], figsize=(20, 12), z_score=1, cbar_kws={"label": "Z-score of enrichment\nof differential regions"}, rasterized=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.z_score.svg" % gene_set_library), bbox_inches="tight", dpi=300)

    def differential_expression_analysis(
            self, samples, trait="knockout",
            variables=["knockout", "replicate"],
            output_suffix="deseq_expression_knockout"):
        """
        Discover differential regions across samples that are associated with a certain trait.
        """
        sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]

        # Get matrix of counts
        counts_matrix = self.expression_matrix_counts[[s.name for s in sel_samples]]

        # Get experiment matrix
        experiment_matrix = pd.DataFrame([sample.as_series() for sample in sel_samples], index=[sample.name for sample in sel_samples])
        # keep only variables
        experiment_matrix = experiment_matrix[["sample_name"] + variables].fillna("Unknown")

        # Make output dir
        output_dir = os.path.join(self.results_dir, output_suffix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Run DESeq2 analysis
        deseq_table = DESeq_analysis(
            counts_matrix, experiment_matrix, trait, covariates=[x for x in variables if x != trait], output_prefix=os.path.join(output_dir, output_suffix), alpha=0.05)

        # to just read in
        # deseq_table = pd.read_csv(os.path.join(output_dir, output_suffix) + ".%s.csv" % trait, index_col=0)
        # self.coverage_annotated = pd.read_csv(os.path.join(self.results_dir, "breg_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index_col=0)

        df = self.expression.join(deseq_table)
        df['comparison'] = df['comparison'].str.replace("-WT", "")
        df.to_csv(os.path.join(output_dir, output_suffix) + ".%s.annotated.csv" % trait)
        df = pd.read_csv(os.path.join(output_dir, output_suffix) + ".%s.annotated.csv" % trait, index_col=0)
        df['comparison'] = df['comparison'].str.replace("-WT", "")

        # Extract significant based on p-value and fold-change
        diff = df[(df["padj"] < 0.01) & (abs(df["log2FoldChange"]) > 1.)]

        if diff.shape[0] < 1:
            print("No significantly different regions found.")
            return

        # groups = list(set([getattr(s, trait) for s in sel_samples]))
        comparisons = pd.Series(df['comparison'].unique())
        groups = sorted(comparisons.tolist())

        # Expression of the knocked-out genes
        genes = set(s.knockout for s in sel_samples).discard("WT")

        # ordered
        x = self.expression.loc[genes, [s.name for s in sel_samples]].dropna().sort_index()

        fig, axis = plt.subplots(1, figsize=(12, 4))
        sns.heatmap(
            x,
            # xticklabels=self.expression.columns.get_level_values("sample_name"),
            cbar_kws={"label": "log2(1 + TPM)"}, cmap="Spectral_r", vmin=-2, ax=axis, square=True)
        axis.set_xlabel("Sample")
        axis.set_ylabel("Gene")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
        fig.savefig(os.path.join(self.results_dir, "expression.knocked_out_genes.heatmap.svg"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1, figsize=(12, 4))
        sns.heatmap(
            x.apply(lambda x: (x - x.mean()) / x.std(), axis=1),
            # xticklabels=self.expression.columns.get_level_values("sample_name"),
            cbar_kws={"label": "Z-score log2(TPM)"}, ax=axis, square=True)
        axis.set_xlabel("Sample")
        axis.set_ylabel("Gene")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
        fig.savefig(os.path.join(self.results_dir, "expression.knocked_out_genes.heatmap.z_score.svg"), bbox_inches="tight", dpi=300)

        # clustered
        g = sns.clustermap(
            x,
            cbar_kws={"label": "log2(1 + TPM)"}, cmap="Spectral_r", vmin=-2)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(self.results_dir, "expression.knocked_out_genes.clustermap.svg"), bbox_inches="tight", dpi=300)

        g = sns.clustermap(
            x,
            cbar_kws={"label": "Z-score log2(TPM)"}, z_score=0)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(self.results_dir, "expression.knocked_out_genes.clustermap.z_score.svg"), bbox_inches="tight", dpi=300)

        # Statistics of differential regions
        import string
        total_sites = float(self.expression.shape[0])

        total_diff = diff.groupby(["comparison"])['stat'].count().sort_values(ascending=False)
        fig, axis = plt.subplots(1)
        sns.barplot(total_diff.values, total_diff.index, orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.total.svg" % (output_suffix, trait)), bbox_inches="tight")
        # percentage of total
        fig, axis = plt.subplots(1)
        sns.barplot(
            (total_diff.values / total_sites) * 100,
            total_diff.index,
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.total_percentage.svg" % (output_suffix, trait)), bbox_inches="tight")

        # direction-dependent
        diff["direction"] = diff["log2FoldChange"].apply(lambda x: "up" if x >= 0 else "down")

        split_diff = diff.groupby(["comparison", "direction"])['stat'].count().sort_values(ascending=False)
        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.barplot(
            split_diff.values,
            split_diff.reset_index()[['comparison', 'direction']].apply(string.join, axis=1),
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split.svg" % (output_suffix, trait)), bbox_inches="tight")
        # percentage of total
        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.barplot(
            (split_diff.values / total_sites) * 100,
            split_diff.reset_index()[['comparison', 'direction']].apply(string.join, axis=1),
            orient="h", ax=axis)
        for t in axis.get_xticklabels():
            t.set_rotation(0)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split_percentage.svg" % (output_suffix, trait)), bbox_inches="tight")

        # in same
        split_diff = diff.groupby(["comparison", "direction"])['stat'].count().sort_values(ascending=False)

        split_diff = split_diff.reset_index()
        split_diff.loc[split_diff['direction'] == "down", "stat"] = -split_diff.loc[split_diff['direction'] == "down", "stat"]
        fig, axis = plt.subplots(1, figsize=(6, 3))
        sns.barplot(
            x="comparison", y="stat", data=split_diff, hue="direction",  # split=False,
            order=split_diff.groupby(['comparison'])['stat'].apply(lambda x: sum(abs(x))).sort_values(ascending=False).index,
            orient="v", ax=axis)
        axis.axhline(0, color="black", linestyle="--")
        axis.set_xticklabels(axis.get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split.pos_neg.svg" % (output_suffix, trait)), bbox_inches="tight")

        # percentage of total
        split_diff['stat'] = (split_diff['stat'].values / total_sites) * 100
        fig, axis = plt.subplots(1, figsize=(6, 3))
        sns.barplot(
            x="comparison", y="stat", data=split_diff, hue="direction",  # split=False,
            order=split_diff.groupby(['comparison'])['stat'].apply(lambda x: sum(abs(x))).sort_values(ascending=False).index,
            orient="v", ax=axis)
        axis.axhline(0, color="black", linestyle="--")
        axis.set_xticklabels(axis.get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.split_percentage.pos_neg.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Overlap between sets
        # # Pyupset
        # import pyupset as pyu
        # # Build dict
        # diff["comparison_direction"] = diff[["comparison", "direction"]].apply(string.join, axis=1)
        # df_dict = {group: diff[diff["comparison_direction"] == group].reset_index()[['index']] for group in set(diff["comparison_direction"])}
        # # Plot
        # plot = pyu.plot(df_dict, unique_keys=['index'], inters_size_bounds=(10, np.inf))
        # plot['figure'].set_size_inches(20, 8)
        # plot['figure'].savefig(os.path.join(output_dir, "%s.%s.number_differential.upset.svg" % (output_suffix, trait)), bbox_inched="tight")
        diff["intersect"] = 1
        diff["direction"] = diff["log2FoldChange"].apply(lambda x: "up" if x >= 0 else "down")
        piv = pd.pivot_table(diff.reset_index(), index='gene_name', columns=['comparison', 'direction'], values='intersect', fill_value=0)

        import itertools
        intersections = pd.DataFrame(columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
        for ((k1, dir1), i1), ((k2, dir2), i2) in itertools.permutations(piv.T.groupby(level=['comparison', 'direction']).groups.items(), 2):
            print(k1, k2)
            i1 = set(piv[i1][piv[i1] == 1].dropna().index)
            i2 = set(piv[i2][piv[i2] == 1].dropna().index)
            intersections = intersections.append(
                pd.Series(
                    [k1, k2, dir1, dir2, len(i1), len(i2), len(i1.intersection(i2)), len(i1.union(i2))],
                    index=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"]
                ),
                ignore_index=True
            )
        # convert to %
        intersections['intersection_perc'] = (intersections['intersection'] / intersections['size2']) * 100.

        # make pivot tables
        piv_up = pd.pivot_table(intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "up")], index="group1", columns="group2", values="intersection_perc").fillna(0)
        piv_down = pd.pivot_table(intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "down")], index="group1", columns="group2", values="intersection_perc").fillna(0)
        np.fill_diagonal(piv_up.values, np.nan)
        np.fill_diagonal(piv_down.values, np.nan)

        # heatmaps
        fig, axis = plt.subplots(1, 2, figsize=(8 * 2, 8))
        sns.heatmap(piv_down, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis[0])
        sns.heatmap(piv_up, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis[1])
        axis[0].set_title("Downregulated regions")
        axis[1].set_title("Upregulated regions")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="right")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="right")
        axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
        axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.up_down_split.svg" % (output_suffix, trait)), bbox_inches="tight")

        # combined heatmap
        # with upregulated regions in upper square matrix and downredulated in down square
        piv_combined = pd.DataFrame(np.triu(piv_up), index=piv_up.index, columns=piv_up.columns).replace(0, np.nan)
        piv_combined.update(pd.DataFrame(np.tril(piv_down), index=piv_down.index, columns=piv_down.columns).replace(0, np.nan))
        piv_combined = piv_combined.fillna(0)
        np.fill_diagonal(piv_combined.values, np.nan)

        fig, axis = plt.subplots(1, figsize=(8, 8))
        sns.heatmap(piv_combined, square=True, cmap="summer", cbar_kws={"label": "Concordant regions (% of group 2)"}, ax=axis)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.up_down_together.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Observe disagreement between knockouts
        # (overlap of down-regulated with up-regulated and vice-versa)
        piv_up = pd.pivot_table(intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "down")], index="group1", columns="group2", values="intersection")
        piv_down = pd.pivot_table(intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "up")], index="group1", columns="group2", values="intersection")

        piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()

        fig, axis = plt.subplots(1, 2, figsize=(16, 8))
        sns.heatmap(piv_disagree, square=True, cmap="Greens", cbar_kws={"label": "Discordant regions"}, ax=axis[0])
        sns.heatmap(np.log2(1 + piv_disagree), square=True, cmap="Greens", cbar_kws={"label": "Discordant regions (log2)"}, ax=axis[1])

        norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
        cmap = plt.get_cmap("Greens")
        log_norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log2(1 + piv_disagree).max().max())
        for i, g1 in enumerate(piv_disagree.columns):
            for j, g2 in enumerate(piv_disagree.index):
                axis[0].scatter(
                    len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                    s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1, color=cmap(norm(piv_disagree.loc[g1, g2])), marker="o")
                axis[1].scatter(
                    len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                    s=(100 ** (log_norm(np.log2(1 + piv_disagree).loc[g1, g2]))) - 1, color=cmap(log_norm(np.log2(1 + piv_disagree).loc[g1, g2])), marker="o")

        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, "%s.%s.number_differential.overlap.disagreement.svg" % (output_suffix, trait)), bbox_inches="tight")

        # Pairwise scatter plots
        cond2 = "WT"
        n_rows = n_cols = int(np.ceil(np.sqrt(len(comparisons))))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            # Hexbin plot
            axis.hexbin(np.log2(1 + df2[cond1]), np.log2(1 + df2[cond2]), alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)
            axis.set_xlabel(cond1)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # Scatter plot
                axis.scatter(np.log2(1 + diff2[cond1]), np.log2(1 + diff2[cond2]), alpha=0.1, color="red", s=2)
            m = max(np.log2(1 + df2[cond1]).max(), np.log2(1 + df2[cond2]).max())
            axis.plot([0, m], [0, m], color="black", alpha=0.8, linestyle="--")
            axis.set_ylabel(cond1)
            axis.set_ylabel(cond2)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.scatter_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Volcano plots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            # hexbin
            axis.hexbin(df2["log2FoldChange"], -np.log10(df2['pvalue']), alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # significant scatter
                axis.scatter(diff2["log2FoldChange"], -np.log10(diff2['pvalue']), alpha=0.2, color="red", s=2)
            axis.axvline(0, linestyle="--", color="black", alpha=0.8)
            axis.set_title(comparison)
            axis.set_xlabel("log2(fold change)")
            axis.set_ylabel("-log10(p-value)")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.volcano_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # MA plots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), sharex=True, sharey=True)
        if n_rows > 1 or n_cols > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()
            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            df2 = df[df["comparison"] == comparison]
            if df2.shape[0] == 0:
                continue
            axis = axes.next()

            # hexbin
            axis.hexbin(np.log2(df2["baseMean"]), df2["log2FoldChange"], alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            diff2 = diff[diff["comparison"] == comparison]
            if diff2.shape[0] > 0:
                # significant scatter
                axis.scatter(np.log2(diff2["baseMean"]), diff2["log2FoldChange"], alpha=0.2, color="red", s=2)
            axis.axhline(0, linestyle="--", color="black", alpha=0.8)
            axis.set_title(comparison)
            axis.set_xlabel("Intensity")
            axis.set_ylabel("log2(fold change)")
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "%s.%s.ma_plots.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # save unique differential regions
        diff2 = diff[groups].ix[diff.index.unique()].drop_duplicates()
        diff2.to_csv(os.path.join(output_dir, "%s.%s.differential_regions.csv" % (output_suffix, trait)))

        # Exploration of differential regions
        # get unique differential regions
        df2 = diff2.join(self.expression)

        # Characterize regions
        prefix = "%s.%s.diff_regions" % (output_suffix, trait)
        # region's structure
        # characterize_regions_structure(df=df2, prefix=prefix, output_dir=output_dir)
        # region's function
        # characterize_regions_function(df=df2, prefix=prefix, output_dir=output_dir)

        # Heatmaps
        # Comparison level
        g = sns.clustermap(np.log2(1 + df2[groups]).corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non differential genes"}, cmap="Spectral_r", metric="correlation", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.groups.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300, metric="correlation")

        g = sns.clustermap(np.log2(1 + df2[groups]), yticklabels=False, cbar_kws={"label": "Accessibility of\ndifferential genes"}, cmap="BuGn", metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.groups.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(np.log2(1 + df2[groups]), yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of expression\non differential genes"}, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.groups.clustermap.z0.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Fold-changes and P-values
        # pivot table of genes vs comparisons
        sig = diff.index.drop_duplicates()
        fold_changes = pd.pivot_table(df.reset_index(), index="gene_name", columns="comparison", values="log2FoldChange")
        p_values = -np.log10(pd.pivot_table(df.reset_index(), index="gene_name", columns="comparison", values="padj"))

        # fold
        g = sns.clustermap(fold_changes.corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non fold-changes"}, cmap="Spectral_r", vmin=0, vmax=1, metric="correlation", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.groups.fold_changes.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300, metric="correlation")

        g = sns.clustermap(fold_changes.ix[sig], yticklabels=False, cbar_kws={"label": "Fold-changes of\ndifferential genes"}, vmin=-4, vmax=4, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.groups.fold_changes.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        # Sample level
        g = sns.clustermap(df2[[s.name for s in sel_samples]].corr(), xticklabels=False, cbar_kws={"label": "Pearson correlation\non differential genes"}, cmap="Spectral_r", metric="correlation", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.samples.clustermap.corr.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(df2[[s.name for s in sel_samples]], yticklabels=False, cbar_kws={"label": "Expression of\ndifferential genes"}, cmap="BuGn", vmin=0, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.samples.clustermap.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        g = sns.clustermap(df2[[s.name for s in sel_samples]], yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of expression\non differential genes"}, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, "%s.%s.diff_genes.samples.clustermap.z0.svg" % (output_suffix, trait)), bbox_inches="tight", dpi=300)

        #
        # import scipy
        # D = scipy.spatial.distance.pdist(fold_changes.ix[sig], "euclidean")
        # Z = scipy.cluster.hierarchy.linkage(D, 'complete')
        # dn = scipy.cluster.hierarchy.dendrogram(Z, truncate_mode='lastp', p=2000)

        # Examine each region cluster
        max_diff = 1000
        pathway_enr = pd.DataFrame()
        for cond1 in sorted(groups):
            # get comparison
            comparison = comparisons[comparisons.str.contains(cond1)].squeeze()

            if type(comparison) is pd.Series:
                if len(comparison) > 1:
                    comparison = comparison.iloc[0]

            # Separate in up/down-regulated genes
            for f, direction, top in [(np.less, "down", "head"), (np.greater, "up", "tail")]:
                comparison_df = self.expression.ix[diff[
                    (diff["comparison"] == comparison) &
                    (f(diff["log2FoldChange"], 0))
                ].index]

                # Handle extremes of regions
                if comparison_df.shape[0] < 1:
                    continue

                if comparison_df.shape[0] > max_diff:
                    comparison_df = comparison_df.ix[
                        getattr(
                            diff[
                                (diff["comparison"] == comparison) &
                                (f(diff["log2FoldChange"], 0))]
                            ["log2FoldChange"].sort_values(), top)
                        (max_diff).index]

                # Characterize regions
                prefix = "%s.%s.diff_regions.comparison_%s.%s" % (output_suffix, trait, comparison, direction)

                comparison_dir = os.path.join(output_dir, prefix)
                if not os.path.exists(comparison_dir):
                    os.makedirs(comparison_dir)

                print("Doing regions of comparison %s, with prefix %s" % (comparison, prefix))
                if not os.path.exists(os.path.join(comparison_dir, "enrichr.csv")):
                    comparison_df.reset_index()[['gene_name']].drop_duplicates().to_csv(os.path.join(comparison_dir, "genes.txt"), index=False)
                    enr = enrichr(comparison_df.reset_index())
                    enr.to_csv(os.path.join(comparison_dir, "enrichr.csv"), index=False)
                else:
                    enr = pd.read_csv(os.path.join(comparison_dir, "enrichr.csv"))
                    enr["comparison"] = prefix
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)

        # write combined enrichments
        pathway_enr.to_csv(
            os.path.join(output_dir, "%s.%s.diff_genes.enrichr.csv" % (output_suffix, trait)), index=False)

    def investigate_differential_genes(
            self, samples, trait="knockout",
            variables=["knockout", "replicate"],
            output_suffix="deseq_expression_knockout", n=5, method="groups"):
        output_dir = os.path.join(self.results_dir, output_suffix)

        # ENRICHR
        # read in
        pathway_enr = pd.read_csv(os.path.join(output_dir, "%s.%s.diff_genes.enrichr.csv" % (output_suffix, trait)))

        pathway_enr = pathway_enr[~pathway_enr["comparison"].str.contains("SMARCC2")]

        # pretty names
        pathway_enr["comparison"] = pathway_enr["comparison"].str.extract("%s.%s.diff_regions.comparison_(.*)" % (output_suffix, trait), expand=True)
        pathway_enr["description"] = pathway_enr["description"].str.decode("utf-8")

        for gene_set_library in pathway_enr["gene_set_library"].unique()[9:]:
            print(gene_set_library)
            if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
                continue

            # pivot table
            enrichr_pivot = pd.pivot_table(
                pathway_enr[pathway_enr["gene_set_library"] == gene_set_library],
                values="p_value", columns="description", index="comparison").fillna(1)

            # transform p-values
            enrichr_pivot = -np.log10(enrichr_pivot.fillna(1))
            enrichr_pivot = enrichr_pivot.replace({np.inf: 300})

            # plot correlation
            g = sns.clustermap(enrichr_pivot.T.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential genes"})
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            # g.fig.savefig(os.path.join(output_dir, "enrichr.%s.correlation.svg" % gene_set_library), bbox_inches="tight")
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.correlation.svg" % gene_set_library), bbox_inches="tight", dpi=300)

            top = pathway_enr[pathway_enr["gene_set_library"] == gene_set_library].set_index('description').groupby("comparison")['p_value'].nsmallest(n)
            top_terms = top.index.get_level_values('description').unique()
            # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))], figsize=(20, 12), cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            # g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.svg" % gene_set_library), bbox_inches="tight")
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.svg" % gene_set_library), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))], figsize=(20, 12), z_score=1, cbar_kws={"label": "Z-score of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            # g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.z_score.svg" % gene_set_library), bbox_inches="tight")
            g.fig.savefig(os.path.join(output_dir, "enrichr.%s.cluster_specific.z_score.svg" % gene_set_library), bbox_inches="tight", dpi=300)

            for direction in ["up", "down"]:
                # pivot table
                enrichr_pivot = pd.pivot_table(
                    (pathway_enr[(pathway_enr["gene_set_library"] == gene_set_library) & (pathway_enr["comparison"].str.contains(direction))]),
                    values="p_value", columns="description", index="comparison").fillna(1)

                # transform p-values
                enrichr_pivot = -np.log10(enrichr_pivot.fillna(1))
                enrichr_pivot = enrichr_pivot.replace({np.inf: 300})

                top = pathway_enr[pathway_enr["gene_set_library"] == gene_set_library].set_index('description').groupby("comparison")['p_value'].nsmallest(n)
                top_terms = [t for t in top.index.get_level_values('description').unique().astype(str).astype(set).tolist() if t in enrichr_pivot.columns]
                # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

                # plot clustered heatmap
                g = sns.clustermap(enrichr_pivot[top_terms], figsize=(20, 12), cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                # g.fig.savefig(os.path.join(output_dir, "enrichr.%s.only_%s.cluster_specific.svg" % (gene_set_library, direction)), bbox_inches="tight")
                g.fig.savefig(os.path.join(output_dir, "enrichr.%s.only_%s.cluster_specific.svg" % (gene_set_library, direction)), bbox_inches="tight", dpi=300)

                # plot clustered heatmap
                g = sns.clustermap(enrichr_pivot[top_terms], figsize=(20, 12), z_score=1, cbar_kws={"label": "Z-score of enrichment\nof differential genes"}, metric="correlation")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize=4)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                # g.fig.savefig(os.path.join(output_dir, "enrichr.%s.only_%s.cluster_specific.z_score.svg" % (gene_set_library, direction)), bbox_inches="tight")
                g.fig.savefig(os.path.join(output_dir, "enrichr.%s.only_%s.cluster_specific.z_score.svg" % (gene_set_library, direction)), bbox_inches="tight", dpi=300)

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

            nx.write_graphml(g, '/home/arendeiro/net_go.graphml')

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


def DESeq_analysis(counts_matrix, experiment_matrix, variable, covariates, output_prefix, alpha=0.05):
    """
    """
    import rpy2.robjects as robj
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    run = robj.r("""
        run = function(countData, colData, variable, covariates, output_prefix, alpha) {
            library(DESeq2)

            alpha = 0.05
            output_prefix = "results/deseq_expression_knockout/deseq_expression_knockout"
            countData = read.csv("results/deseq_expression_knockout/counts_matrix.csv", sep=",", row.names=1)
            colData = read.csv("results/deseq_expression_knockout/experiment_matrix.csv", sep=",")
            variable = "knockout"
            covariates = "replicate + "

            colData$knockout = as.character(colData$knockout)
            colData[colData$clone == "GFP", "knockout"] = "WT_GFP"
            colData$knockout = as.factor(colData$knockout)

            colData$replicate = as.factor(colData$replicate)
            colData$clone = as.factor(colData$clone)

            design = as.formula((paste("~", covariates, variable)))
            print(design)
            dds <- DESeqDataSetFromMatrix(
                countData = countData, colData = colData,
                design)

            dds <- DESeq(dds, parallel=TRUE)
            save(dds, file=paste0(output_prefix, ".deseq_dds_object.Rdata"))
            # load(paste0(output_prefix, ".deseq_dds_object.Rdata"))

            # Get group means
            # get groups with one sample and set mean to the value of that sample
            single_levels = names(table(colData[, variable])[table(colData[, variable]) == 1])
            single_levels_values = sapply(
                single_levels,
                function(lvl) counts(dds, normalized=TRUE)[, dds[, variable] == lvl]
            )
            # all others, get sample means
            multiple_levels = names(table(colData[, variable])[table(colData[, variable]) > 1])
            multiple_levels_values = sapply(
                multiple_levels,
                function(lvl) rowMeans(counts(dds, normalized=TRUE)[, colData[, variable] == lvl])
            )
            group_means = cbind(single_levels_values, multiple_levels_values)
            rownames(group_means) = rownames(countData)
            write.table(group_means, paste0(output_prefix, ".", variable, ".group_means.csv"), sep=",")

            # pairwise combinations
            knockouts = sort(unique(colData[, variable]), descending=FALSE)

            # keep track of output files
            result_files = list()

            for (i in 1:length(knockouts)) {

                cond1 = as.character(knockouts[i])
                cond2 = "WT"
                if (cond1 == cond2){
                    next
                }
                contrast = c(variable, cond1, cond2)
                print(contrast)

                # get results
                res <- results(dds, contrast=contrast, alpha=alpha, independentFiltering=TRUE, parallel=TRUE)
                res <- as.data.frame(res)

                # append group means
                res <- cbind(group_means, res)

                # append to results
                comparison_name = paste(cond1, cond2, sep="-")
                output_name = paste0(output_prefix, ".", variable, ".", comparison_name, ".csv")
                res["comparison"] = comparison_name

                # coherce to character
                res = data.frame(lapply(res, as.character), stringsAsFactors=FALSE)

                # add index
                rownames(res) = rownames(countData)

                write.table(res, output_name, sep=",")
                result_files[i] = output_name
            }

            # ARID1A vs ARID1B
            i = length(knockouts) + 1
            cond1 = "ARID1A"
            cond2 = "ARID1B"
            contrast = c(variable, cond1, cond2)

            # get results
            res <- results(dds, contrast=contrast, alpha=alpha, independentFiltering=TRUE, parallel=TRUE)
            res <- as.data.frame(res)

            # append group means
            res <- cbind(group_means, res)

            # append to results
            comparison_name = paste(cond1, cond2, sep="-")
            output_name = paste0(output_prefix, ".", variable, ".", comparison_name, ".csv")
            res["comparison"] = comparison_name

            # coherce to character
            res = data.frame(lapply(res, as.character), stringsAsFactors=FALSE)

            # add index
            rownames(res) = rownames(countData)

            write.table(res, output_name, sep=",")
            result_files[i] = output_name

            # SMARCA4 vs SMARCA2
            i = length(knockouts) + 1
            cond1 = "SMARCA4"
            cond2 = "SMARCA2"
            contrast = c(variable, cond1, cond2)

            # get results
            res <- results(dds, contrast=contrast, alpha=alpha, independentFiltering=TRUE, parallel=TRUE)
            res <- as.data.frame(res)

            # append group means
            res <- cbind(group_means, res)

            # append to results
            comparison_name = paste(cond1, cond2, sep="-")
            output_name = paste0(output_prefix, ".", variable, ".", comparison_name, ".csv")
            res["comparison"] = comparison_name

            # coherce to character
            res = data.frame(lapply(res, as.character), stringsAsFactors=FALSE)

            # add index
            rownames(res) = rownames(countData)

            write.table(res, output_name, sep=",")
            result_files[i] = output_name


            return(result_files)
        }

    """)

    # replace names
    counts_matrix.columns = ["S" + str(i) for i in range(len(counts_matrix.columns))]
    experiment_matrix.index = ["S" + str(i) for i in range(len(experiment_matrix.index))]
    experiment_matrix.index.name = "sample"

    # save to disk just in case
    counts_matrix.to_csv(os.path.join(os.path.dirname(output_prefix), "counts_matrix.csv"), index=True)
    experiment_matrix.to_csv(os.path.join(os.path.dirname(output_prefix), "experiment_matrix.csv"), index=True)

    result_files = run(counts_matrix, experiment_matrix, variable, " + ".join(covariates) + " + " if len(covariates) > 0 else "", output_prefix, alpha)

    # concatenate all files
    import glob
    # result_files = glob.glob(output_prefix + ".*-WT.csv")
    result_files = glob.glob(output_prefix + "*-*.csv")
    results = pd.DataFrame()
    for result_file in result_files:
        df = pd.read_csv(result_file)
        # df.index = counts_matrix.index

        results = results.append(df)

    # save all
    results.to_csv(os.path.join(output_prefix + ".%s.csv" % variable), index=True)

    # return
    return results


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

    # for F in `find . -name "*symbols.txt"`
    # do
    #     if  [ ! -f ${F/symbols.txt/enrichr.csv} ]; then
    #         echo $F
    #         wc -l $F
    #         sbatch -J ENRICHR_${F} -o ${F/symbols.txt/enrichr.log} ~/run_Enrichr.sh $F
    #     fi
    # done


def lola(bed_files, universe_file, output_folder):
    """
    Performs location overlap analysis (LOLA) on bedfiles with regions sets.
    """
    import rpy2.robjects as robj

    run = robj.r("""
        function(bedFiles, universeFile, outputFolder) {
            library("LOLA")

            userUniverse  <- LOLA::readBed(universeFile)

            dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/"
            dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/hg19/"
            regionDB = loadRegionDB(c(dbPath1, dbPath2))

            if (typeof(bedFiles) == "character") {
                userSet <- LOLA::readBed(bedFiles)
                lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                lolaResults[order(support, decreasing=TRUE), ]
                writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
            } else if (typeof(bedFiles) == "double") {
                for (bedFile in bedFiles) {
                    userSet <- LOLA::readBed(bedFile)
                    lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                    lolaResults[order(support, decreasing=TRUE), ]
                    writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
                }
            }
        }
    """)

    # convert the pandas dataframe to an R dataframe
    run(bed_files, universe_file, output_folder)

    # for F in `find . -iname *_regions.bed`
    # do
    #     if  [ ! -f `dirname $F`/allEnrichments.txt ]; then
    #         echo $F
    #         sbatch -J LOLA_${F} -o ${F/_regions.bed/}_lola.log ~/run_LOLA.sh $F ~/projects/baf-kubicek/results/baf-kubicek_peak_set.bed hg19 `dirname $F`
    #     fi
    # done


def bed_to_fasta(bed_file, fasta_file):
    # write name column
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed['name'] = bed[0] + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
    bed[1] = bed[1].astype(int)
    bed[2] = bed[2].astype(int)
    bed.to_csv(bed_file + ".tmp.bed", sep='\t', header=None, index=False)

    # do enrichment
    cmd = "twoBitToFa ~/resources/genomes/hg19/hg19.2bit -bed={0} {1}".format(bed_file + ".tmp.bed", fasta_file)

    os.system(cmd)
    # os.system("rm %s" % bed_file + ".tmp.bed")


def meme_ame(input_fasta, output_dir, background_fasta=None):
    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 1 -f {0} > {1}
        """.format(input_fasta, shuffled)
        os.system(cmd)

    cmd = """
    ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \\
    --control {0} -o {1} {2} ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv9.meme
    """.format(background_fasta if background_fasta is not None else shuffled, output_dir, input_fasta)
    os.system(cmd)

    os.system("rm %s" % shuffled)

    # for F in `find . -iname *fa`
    # do
    #     if  [ ! -f `dirname $F`/ame.txt ]; then
    #         echo $F
    #         sbatch -J MEME-AME_${F} -o ${F/fa/}ame.log ~/run_AME.sh $F human
    #     fi
    # done


def homer_motifs(bed_file, output_dir, genome="hg19"):
    cmd = "findMotifsGenome.pl {bed} {genome}r {out_dir}  -size 1000 -h -p 2 -len 8,10,12,14 -noknown".format(
        bed=bed_file, genome=genome, out_dir=output_dir
    )
    os.system(cmd)

    # # for missing
    # for F in `find . -iname *regions.bed`
    # do
    #     if  [ ! -f `dirname $F`/homerResults.html ]; then
    #         echo $F
    #         DIRNAME=`dirname $F`
    #         sbatch -J HOMER_MOTIF_${F} -o ${F/bed/homer_motif.log} ~/run_homer_motifs.sh ${F} hg19 ${DIRNAME}
    #     fi
    # done


def parse_ame(ame_dir):

    with open(os.path.join(ame_dir, "ame.txt"), 'r') as handle:
        lines = handle.readlines()

    output = list()
    for line in lines:
        # skip header lines
        if line[0] not in [str(i) for i in range(10)]:
            continue

        # get motif string and the first half of it (simple name)
        motif = line.strip().split(" ")[5].split("_")[0]
        # get corrected p-value
        q_value = float(line.strip().split(" ")[-2])
        # append
        output.append((motif, q_value))

    return pd.Series(dict(output))


def parse_homer(homer_html):
    import re

    with open(homer_html, 'r') as handle:
        content = handle.read().strip()

    content2 = re.sub("(\d)\n", "\g<1>", content)
    rows = [x for x in content2.split("\n") if x[:4] == "<TR>"]

    additional_header = ["STD(Bg STD)", "Best Match/Details", "Motif File"]
    df = pd.DataFrame(
        [re.sub("</TD>", "", re.sub("</TR>", "", re.sub("<TR>", "", x))).split("<TD>") for x in rows[1:]],
        columns=re.sub("</TD>", "", re.sub("</TR>", "", re.sub("<TR>", "", rows[0]))).split("<TD>") + additional_header)
    df = pd.concat([df, df["Best Match/Details"].str.extract("^(?P<motif_name>.*)\((?P<motif_match_score>.*)\)<BR/")], axis=1)
    df["-log_p_value"] = -df["log P-pvalue"].astype(float)
    df["% of Targets"] = df["% of Targets"].str.extract("(.*)%").astype(float)
    df["% of Background"] = df["% of Background"].str.extract("(.*)%").astype(float)
    df["fold_enrichment_over_background"] = df["% of Targets"] / df["% of Background"]
    df["log2_fold_enrichment_over_background"] = np.log2(df["fold_enrichment_over_background"])

    return df


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


def characterize_regions_function(df, output_dir, prefix, results_dir="results", universe_file=None):
    def scale_score(x):
        return (x - x.min()) / (x.max() - x.min())

    # use all sites as universe
    if universe_file is None:
        universe_file = os.path.join(analysis.results_dir, analysis.name + "_peak_set.bed")

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # save to bed
    bed_file = os.path.join(output_dir, "%s_regions.bed" % prefix)
    df[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)
    # save as tsv
    tsv_file = os.path.join(output_dir, "%s_regions.tsv" % prefix)
    df[['chrom', 'start', 'end']].reset_index().to_csv(tsv_file, sep="\t", header=False, index=False)

    # export ensembl gene names
    df['gene_name'].str.split(",").apply(pd.Series, 1).stack().drop_duplicates().to_csv(os.path.join(output_dir, "%s_genes.symbols.txt" % prefix), index=False)

    # export gene symbols with scaled absolute fold change
    if "log2FoldChange" in df.columns:
        df["score"] = scale_score(abs(df["log2FoldChange"]))
        df["abs_fc"] = abs(df["log2FoldChange"])

        d = df[['gene_name', 'score']].sort_values('score', ascending=False)

        # split gene names from score if a reg.element was assigned to more than one gene
        a = (
            d['gene_name']
            .str.split(",")
            .apply(pd.Series, 1)
            .stack()
        )
        a.index = a.index.droplevel(1)
        a.name = 'gene_name'
        d = d[['score']].join(a)
        # reduce various ranks to mean per gene
        d = d.groupby('gene_name').mean().reset_index()
        d.to_csv(os.path.join(output_dir, "%s_genes.symbols.score.csv" % prefix), index=False)

    # Motifs
    fasta_file = os.path.join(output_dir, "%s_regions.fa" % prefix)
    bed_to_fasta(bed_file, fasta_file)

    # motif enrichment
    meme_ame(fasta_file, output_dir)

    # de novo motif finding
    homer_motifs(bed_file, output_dir, genome="hg19")

    # Lola
    try:
        lola(bed_file, universe_file, output_dir)
    except:
        print("LOLA analysis for %s failed!" % prefix)

    # Enrichr
    results = enrichr(df[['chrom', 'start', 'end', "gene_name"]])

    # Save
    results.to_csv(os.path.join(output_dir, "%s_regions.enrichr.csv" % prefix), index=False, encoding='utf-8')


def metagene_plot(bams, labels, output_prefix, region="genebody", genome="hg19"):
    from pypiper import NGSTk
    import textwrap
    import os
    tk = NGSTk()

    job_name = output_prefix
    job_file = output_prefix + ".sh"
    job_log = output_prefix + ".log"

    # write ngsplot config file to disk
    config_file = os.path.join(os.environ['TMPDIR'], "ngsplot_config.txt")
    with open(config_file, "w") as handle:
        for i in range(len(bams)):
            handle.write("\t".join([bams[i], "-1", labels[i]]) + "\n")

    cmd = tk.slurm_header(job_name, job_log, queue="mediumq", time="1-10:00:00", mem_per_cpu=8000, cpus_per_task=8)

    # build plot command
    if region == "genebody":
        cmd += """xvfb-run ngs.plot.r -G {0} -R {1} -C {2} -O {3} -L 3000 -GO km\n""".format(genome, region, config_file, output_prefix)
    elif region == "tss":
        cmd += """xvfb-run ngs.plot.r -G {0} -R {1} -C {2} -O {3} -L 3000 -FL 300\n""".format(genome, region, config_file, output_prefix)

    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def global_changes(samples, trait="knockout"):
    import glob
    import re

    output_dir = os.path.join("data", "merged")
    # sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]
    # groups = sorted(list(set([getattr(s, trait) for s in sel_samples])))
    groups = [os.path.basename(re.sub(".merged.sorted.bam", "", x)) for x in glob.glob(output_dir + "/*.merged.sorted.bam")]

    for region in ["genebody", "tss"]:
        print metagene_plot(
            [os.path.abspath(os.path.join(output_dir, group + ".merged.sorted.bam")) for group in groups],
            groups,
            os.path.abspath(os.path.join(output_dir, "%s.metaplot" % region)),
            region=region
        )


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

    # pair analysis and Project
    analysis.prj = prj
    analysis.samples = prj.samples

    analysis = analysis.from_pickle()

    # work only with ATAC-seq samples
    atacseq_samples = [sample for sample in analysis.samples if sample.library == "ATAC-seq" and sample.cell_line in ["HAP1"]]
    atacseq_samples = [s for s in atacseq_samples if os.path.exists(s.filtered)]  # and s.pass_qc == 1
    rnaseq_samples = [sample for sample in analysis.samples if sample.library == "RNA-seq" and sample.cell_line in ["HAP1"]]
    rnaseq_samples = [s for s in rnaseq_samples if os.path.exists(os.path.join(
                      sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                      "bitSeq",
                      sample.name + ".counts"))]  # and s.pass_qc == 1

    # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT
    # Get consensus peak set from all samples
    if not hasattr(analysis, "sites"):
        analysis.get_consensus_sites(atacseq_samples, "summits")
    if not hasattr(analysis, "support"):
        analysis.calculate_peak_support(atacseq_samples, "summits")

    # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
    # Get coverage values for each peak in each sample of ATAC-seq and ChIPmentation
    analysis.measure_coverage(atacseq_samples)
    # normalize coverage values
    analysis.normalize_coverage_quantiles(atacseq_samples)
    analysis.get_peak_gccontent_length()
    analysis.normalize_gc_content(atacseq_samples)

    # Annotate peaks with closest gene
    analysis.get_peak_gene_annotation()
    # Annotate peaks with genomic regions
    analysis.get_peak_genomic_location()
    # Annotate peaks with ChromHMM state from CD19+ cells
    analysis.get_peak_chromatin_state()
    # Annotate peaks with closest gene, chromatin state,
    # genomic location, mean and variance measurements across samples
    analysis.annotate(atacseq_samples)
    analysis.annotate_with_sample_metadata()
    analysis.to_pickle()

    # Unsupervised analysis
    analysis.unsupervised(atacseq_samples, attributes=["cell_line", "knockout", "replicate", "clone"])

    # Plots
    # plot general peak set features
    analysis.plot_peak_characteristics()
    # Plot coverage features across peaks/samples
    analysis.plot_coverage()
    analysis.plot_variance()

    #

    # Supervised analysis
    analysis.differential_region_analysis(
        atacseq_samples,
        trait="knockout",
        variables=["knockout", "replicate"],
        output_suffix="deseq_knockout")

    # Investigate global changes in accessibility
    global_changes(atacseq_samples)

    #

    rnaseq_samples = [s for s in rnaseq_samples if s.knockout != "SMARCC2" and "GFP" not in s.name]

    # RNA-seq
    analysis.get_gene_expression(samples=rnaseq_samples)
    # Unsupervised analysis
    analysis.unsupervised_expression(rnaseq_samples, attributes=["knockout", "replicate", "clone"])
    # Supervised analysis
    analysis.differential_expression_analysis()

    #

    # See ATAC and RNA together
    analysis.accessibility_expression()

    #

    # Deeper ATAC-seq data
    nucleosome_changes(analysis, atacseq_samples)
    investigate_nucleosome_positions(atacseq_samples)
    phasograms(atacseq_samples)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
