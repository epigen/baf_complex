#!/usr/bin/env python

"""
Call STARR-seq peaks.
"""

from looper.models import Project
from pypiper import NGSTk
import pandas as pd
import os
import textwrap
import pybedtools


tk = NGSTk()


def call_peaks(samples, controls, name, base_path="/scratch/lab_bock/shared/projects/baf-kubicek/data/"):

    for parameter in [("", ""), ("--nomodel --extsize 160 ", "_nomodel")]:

        output_path = os.path.join(base_path, "peaks", name + parameter[1])

        if not os.path.exists(output_path):
            os.mkdir(output_path)

        job_name = "macs2_%s" % name + parameter[1]

        # Build job script
        # slurm header
        cmd = tk.slurm_header(
            job_name + parameter[1],
            os.path.join(output_path, name + parameter[1] + ".log"),
            cpus_per_task=4)

        # load macs2
        cmd += """
\t\t/home/arendeiro/.local/bin/macs2 callpeak {0}-t {1} -c {2} -n {3}{4} --outdir {5}
""".format(parameter[0], " ".join([s.mapped for s in samples]), " ".join([s.mapped for s in controls]), name, parameter[1], output_path)

        # Slurm footer
        cmd += "\t\t" + tk.slurm_footer() + "\n"

        # Write job to file
        job_file = os.path.join(output_path, name + parameter[1] + ".sh")
        with open(job_file, "w") as handle:
            handle.write(textwrap.dedent(cmd))

        # Submit
        tk.slurm_submit_job(job_file)

        print(job_file)


def call_peaks_no_background(samples, name, base_path="/scratch/lab_bock/shared/projects/baf-kubicek/data/"):

    ext = "_no_background"
    output_path = os.path.join(base_path, "peaks", name + ext)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    job_name = "macs2_%s" % name + ext

    # Build job script
    # slurm header
    cmd = tk.slurm_header(
        job_name,
        os.path.join(output_path, name + ext + ".log"),
        cpus_per_task=4)

    # load macs2
    cmd += """
\t\t/home/arendeiro/.local/bin/macs2 callpeak -t {0} -n {1}{2} --outdir {3}
""".format(" ".join([s.mapped for s in samples]), name, ext, output_path)

    # Slurm footer
    cmd += "\t\t" + tk.slurm_footer() + "\n"

    # Write job to file
    job_file = os.path.join(output_path, name + ext + ".sh")
    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)

    print(job_file)


def call_peaks_extra_params(samples, controls, name, base_path="/scratch/lab_bock/shared/projects/baf-kubicek/data/"):

    lengths = [75, 100, 125, 130, 150, 175]

    for parameter in [("--nomodel --extsize %i " % x, "_nomodel-%i" % x) for x in lengths]:

        output_path = os.path.join(base_path, "peaks", name + parameter[1])

        if not os.path.exists(output_path):
            os.mkdir(output_path)

        job_name = "macs2_%s" % name + parameter[1]

        # Build job script
        # slurm header
        cmd = tk.slurm_header(
            job_name,
            os.path.join(output_path, name + parameter[1] + ".log"),
            cpus_per_task=4)

        # load macs2
        cmd += """
\t\t/home/arendeiro/.local/bin/macs2 callpeak {0}-t {1} -c {2} -n {3}{4} --outdir {5}
""".format(parameter[0], " ".join([s.mapped for s in samples]), " ".join([s.mapped for s in controls]), name, parameter[1], output_path)

        # Slurm footer
        cmd += "\t\t" + tk.slurm_footer() + "\n"

        # Write job to file
        job_file = os.path.join(output_path, name + parameter[1] + ".sh")
        with open(job_file, "w") as handle:
            handle.write(textwrap.dedent(cmd))

        # Submit
        tk.slurm_submit_job(job_file)

        print(job_file)


def call_peaks_no_background_single(sample):

    sample.paths.peaks = os.path.join(sample.paths.sample_root, "peaks")
    sample.peaks = os.path.join(sample.paths.peaks, sample.name + "_peaks.narrowPeak")

    if os.path.exists(sample.peaks):
        return

    if not os.path.exists(sample.paths.peaks):
        os.mkdir(sample.paths.peaks)

    job_name = "macs2_%s" % sample.name
    log_file = os.path.join(sample.paths.peaks, sample.name + ".log")
    job_file = os.path.join(sample.paths.peaks, sample.name + ".sh")

    print(log_file)

    # Build job script
    # slurm header
    cmd = tk.slurm_header(job_name, log_file, cpus_per_task=4)

    # load macs2
    cmd += """
\t\t/home/arendeiro/.local/bin/macs2 callpeak --nomodel --extsize 160 -t {0} -n {1} --outdir {2}
""".format(sample.mapped, sample.name, sample.paths.peaks)

    # Slurm footer
    cmd += "\t\t" + tk.slurm_footer() + "\n"

    # Write job to file
    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)


# Start project
prj = Project("metadata/project_config.yaml")
prj.samples = [s for s in prj.samples if s.library == "ChIP-seq"]
for s in prj.samples:
    s.mapped = os.path.join(s.paths.sample_root, "mapped", s.name + ".trimmed.bowtie2.bam")
    s.filtered = os.path.join(s.paths.sample_root, "mapped", s.name + ".trimmed.bowtie2.filtered.bam")

# Read comparison table
comparisons = pd.read_csv("metadata/comparison_table.csv")

ext = "_no_background"
base_path = "/scratch/lab_bock/shared/projects/baf-kubicek/data/"


for comparison in pd.np.unique(comparisons["comparison_name"].dropna()):

    # If there aren't two sides to each comparison, skip it for now
    if len(set(comparisons[
        (comparisons["comparison_name"] == comparison)
    ]["comparison_side"].tolist())) != 2:
        continue

    pos_names = comparisons[
        (comparisons["comparison_name"] == comparison) &
        (comparisons["comparison_side"] == 1)
    ]["sample_name"].tolist()
    neg_names = comparisons[
        (comparisons["comparison_name"] == comparison) &
        (comparisons["comparison_side"] == -1)
    ]["sample_name"].tolist()

    signal_samples = [s for s in prj.samples if s.name in pos_names]
    control_samples = [s for s in prj.samples if s.name in neg_names]

    print([s.name for s in signal_samples], [s.name for s in control_samples], comparison)
    # Call peaks
    call_peaks(signal_samples, control_samples, name=comparison, base_path=prj.metadata.results_subdir)

    # Call peaks without background
    call_peaks_no_background(signal_samples, name=comparison, base_path=prj.metadata.results_subdir)


for sample in prj.samples:
    call_peaks_no_background_single(sample)
