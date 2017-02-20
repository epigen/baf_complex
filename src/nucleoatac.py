#!/usr/bin/env python

"""
This script submits jobs running NucleoATAC for merged sample groups.
"""

import os
import sys
import pandas as pd
from looper.models import Project
from pypiper import NGSTk
import textwrap

global tk
tk = NGSTk()


def nucleoatac(bed_file, bam_file, fasta_file, output_basename, cpus=8):
    """
    Decorator for some methods of Analysis class.
    """
    job_file = os.path.join(output_basename + ".sh")
    log_file = os.path.join(output_basename + ".log")

    cmd = tk.slurm_header(
        "".join(["nucleoatac-", os.path.basename(output_basename)]), log_file,
        cpus_per_task=cpus, time='30-00:00:00', queue="longq", mem_per_cpu=24000)

    cmd += """

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


\t\t# To run the whole nucleoatac procedure:

\t\tnucleoatac run \\
\t\t--write_all \\
\t\t--cores {0} \\
\t\t--bed {1} \\
\t\t--bam {2} \\
\t\t--fasta {3} \\
\t\t--out {4}


\t\t# To run each step separately:

\t\t#nucleoatac occ \\
\t\t#--cores {0} \\
\t\t#--bed {1} \\
\t\t#--bam {2} \\
\t\t#--fasta {3} \\
\t\t#--out {4}

\t\t#nucleoatac nuc \\
\t\t#--write_all \\
\t\t#--cores {0} \\
\t\t#--bed {1} \\
\t\t#--vmat {4}.VMat \\
\t\t#--occ_track {4}.occ.bedgraph.gz \\
\t\t#--bam {2} \\
\t\t#--fasta {3} \\
\t\t#--out {4}

\t\t#nucleoatac merge \\
\t\t#--occpeaks {4}.occpeaks.bed.gz \\
\t\t#--nucpos {4}.nucpos.bed.gz \\
\t\t#--out {4}

\t\t#nucleoatac nfr \\
\t\t#--cores {0} \\
\t\t#--bed {1} \\
\t\t#--occ_track {4}.occ.bedgraph.gz \\
\t\t#--calls {4}.nucmap_combined.bed.gz \\
\t\t#--bam {2} \\
\t\t#--fasta {3} \\
\t\t#--out {4}

\t\t#pyatac vplot \\
\t\t#--out {4} \\
\t\t#--cores {0} \\
\t\t#--bed {1} \\
\t\t#--bam {2} \\
\t\t#--plot_extra \\
\t\t#--lower 20 \\
\t\t#--upper 750 \\
\t\t#--flank 147

\t\tdeactivate
""".format(cpus, bed_file, bam_file, fasta_file, output_basename)
    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def main():
    # Start project
    prj = Project("metadata/project_config.yaml")

    # merged replicates/clones
    merged_dir = os.path.join(prj.metadata.results_subdir, "merged")

    output_dir = os.path.join(os.path.abspath("results"), "nucleoatac")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # static files
    fasta_file = "/data/groups/lab_bock/shared/resources/genomes/hg19/hg19.fa"
    bed_file = "/scratch/lab_bock/shared/projects/baf-kubicek/results/baf-kubicek_peak_set.slop_b500.bed"

    # select only ATAC-seq samples
    df = prj.sheet.df[prj.sheet.df["library"] == "ATAC-seq"]

    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])

        merged_bam = os.path.join(merged_dir, name + ".merged.sorted.bam")
        output_basename = os.path.join(output_dir, name, name)
        if not os.path.exists(os.path.dirname(output_basename)):
            os.mkdir(os.path.dirname(output_basename))

        nucleoatac(bed_file, merged_bam, fasta_file, output_basename)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
