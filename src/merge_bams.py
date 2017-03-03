#!/usr/bin/env python

"""
This script merges bam files from groups of samples
and generated bigWig files from these.
"""

import os
import sys
import pandas as pd
from looper.models import Project
from pypiper import NGSTk
import textwrap
import re

global tk
tk = NGSTk()


def merge_bams(bams, output_bam):
    """
    Decorator for some methods of Analysis class.
    """
    job_file = os.path.join(os.path.dirname(output_bam), os.path.basename(output_bam).split(".")[0] + ".sh")
    log_file = os.path.join(os.path.dirname(output_bam), os.path.basename(output_bam).split(".")[0] + ".log")

    cmd = tk.slurm_header("merge_bams", log_file, cpus_per_task=8, time='10:00:00', queue="shortq", mem_per_cpu=8000)

    cmd += """
\t\tsamtools merge {0} {1}

\t\tsambamba sort -t 8 {0}
""".format(output_bam, " ".join(bams))
    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def bamToBigWig(inputBam, outputBigWig, tagmented=False, normalize=True):
    import os
    import re

    genomeSizes = "/data/groups/lab_bock/shared/resources/genomes/hg19/hg19.chromSizes"
    genome = "hg19"

    cmd = tk.slurm_header("bam_to_bigwig", os.path.join("/scratch/users/arendeiro/", "merge_bams.slurm.log"), cpus_per_task=8, time='6-10:00:00', queue="longq", mem_per_cpu=8000)

    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))

    cmd1 = """
    bedtools bamtobed -i {0} |""".format(inputBam)
    if not tagmented:
        cmd1 += " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
    cmd1 += " genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
            "-5 " if tagmented else "",
            genomeSizes,
            transientFile
    )
    cmd += cmd1

    if normalize:
        cmd += """
    awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov
    """.format(transientFile)

    cmd += """
    bedGraphToBigWig {0}{1}.cov {2} {3}
    """.format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig)

    # remove tmp files
    cmd += """
    if [[ -s {0}.cov ]]; then rm {0}.cov; fi
    """.format(transientFile)
    if normalize:
        cmd += """
    if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi
    """.format(transientFile)

    cmd += """
    chmod 755 {0}
    """.format(outputBigWig)

    cmd += tk.slurm_footer()

    # write job to file
    job_file = "/scratch/users/arendeiro/tmp.sh"
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def get_nucleosome_free_reads(bam, output_bam):
    output_prefix = re.sub(".bam", "", bam)
    job_file = output_prefix + ".get_nucfreads.sh"
    log_file = output_prefix + ".get_nucfreads.log"

    cmd = tk.slurm_header(
        "get_nucleosome_free_reads." + output_prefix.split("/")[-1],
        log_file,
        cpus_per_task=2, time='10:00:00', queue="shortq", mem_per_cpu=8000)

    cmd += """
\t\tsambamba view -f bam -t 2 -o {0} -F "(template_length < 100) and (template_length > -100)" {1}
""".format(output_bam, bam)

    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def get_nucleosome_reads(bam, output_bam):
    output_prefix = re.sub(".bam", "", bam)
    job_file = output_prefix + ".get_nucreads.sh"
    log_file = output_prefix + ".get_nucreads.log"

    cmd = tk.slurm_header(
        "get_nucleosome_reads." + output_prefix.split("/")[-1],
        log_file,
        cpus_per_task=2, time='10:00:00', queue="shortq", mem_per_cpu=8000)

    cmd += """
\t\tsambamba view -f bam -t 2 -o {0} -F "((template_length > 180) and (template_length < 247)) or ((template_length < -180) and (template_length > -247))" {1}
""".format(output_bam, bam)

    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def main():
    # Start project
    prj = Project("metadata/project_config.yaml")

    output_dir = os.path.join(prj.metadata.results_subdir, "merged")
    try:
        os.makedirs(output_dir)
    except:
        pass

    for sample in prj.samples:
        sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")

    for attrs, index in prj.sheet.df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        bams = [s.filtered for s in prj.samples if s.name in prj.sheet.df.loc[index, "sample_name"].tolist()]

        merged_bam = os.path.join(output_dir, name + ".merged.bam")

        if not os.path.exists(merged_bam):
            merge_bams(bams, merged_bam)

        # Split reads
        if not os.path.exists(os.path.join(output_dir, name + ".nucleosome_reads.bam")):
            get_nucleosome_free_reads(
                os.path.join(output_dir, name + ".merged.sorted.bam"),
                os.path.join(output_dir, name + ".nucleosome_free_reads.bam"))
            get_nucleosome_reads(
                os.path.join(output_dir, name + ".merged.sorted.bam"),
                os.path.join(output_dir, name + ".nucleosome_reads.bam"))

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
