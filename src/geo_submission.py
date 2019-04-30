#!/usr/bin/env python

"""
"""


import os
import sys

from peppy import Project
from ngs_toolkit.general import project_to_geo


def collect_md5_sums(df):
    """
    Given a dataframe with columns with paths to md5sum files ending in '_md5sum',
    replace the paths to the md5sum files with the actual checksums.

    Useful to use in combination with ``project_to_geo``.

    :param df: A dataframe with columns ending in '_md5sum'.
    :type df: pandas.DataFrame
    :returns: pandas.DataFrame with md5sum columns replaced with the actual md5sums.
    :rtype: pandas.DataFrame
    """
    import pandas as pd
    cols = df.columns[df.columns.str.endswith("_md5sum")]
    for col in cols:
        for i, path in df.loc[:, col].iteritems():
            if not pd.isnull(path):
                cont = open(path, 'r').read().strip()
                if any([x.isspace() for x in cont]):
                    cont = cont.split(" ")[0]
                df.loc[i, col] = cont
    return df


def main():
    prj = Project(os.path.join("metadata", "project_config.yaml"))

    for sample in prj.samples:
        if sample.protocol == "ATAC-seq":
            sample.peaks = os.path.join(
                prj.metadata.results_subdir, "{0}/peaks/{0}_peaks.narrowPeak".format(sample.name))

    samples = [s for s in prj.samples if s.geo == '1']

    for sample in samples:
        cov = os.path.join(
            sample.paths.sample_root,
            "coverage",
            sample.name + ".bigWig"
        )
        if os.path.exists(cov):
            sample.bigwig = cov
        else:
            cov = os.path.join(
                os.path.expanduser('~'),
                "public_html", "baf_complex",
                "hg19",
                sample.name + ".bigWig"
            )
            if os.path.exists(cov):
                sample.bigwig = cov
            else:
                print(f"FAIL {sample.name} ")
    try:
        annot = project_to_geo(prj, dry_run=True, samples=samples, distributed=True)
        annot.to_csv(os.path.join("geo_submission", "metadata.csv"))
    except Exception:
        raise

    annot = project_to_geo(prj, dry_run=False, samples=samples, distributed=True)

    annot = collect_md5_sums(annot)
    annot.to_csv(os.path.join("geo_submission", "metadata.csv"))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
