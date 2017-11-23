import os

import matplotlib
import pandas as pd
import seaborn as sns

import maplotlib.pyplot as plt

# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


## MASSSPEC
df = pd.read_csv(os.path.join("metadata", 'original', 'ARID1A_A.tab'), sep="\t", index_col=[0, 1])
df = df[df.columns[~df.columns.str.contains("Abund")]]
df.columns = df.columns.str.replace(r"\..*", "")

df2 = df.T.groupby(level=0).mean().T
df3 = df2[df2.index.get_level_values(1).str.contains('GN=')]
df3.index = list(map(lambda x: x[1].split(' ')[0], df3.index.get_level_values(1).str.split("GN=")))
df3 = df3.sort_index()
df4 = df3.dropna()
df5 = np.log2(1 + (df4 / df4.sum(axis=0)) * 1e3)

# Plot all prots
g = sns.clustermap(df5.dropna(), cmap="RdBu_r", robust=True, rasterized=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

# Get variable
qv2 = (df5.std(axis=1) / df5.mean(axis=1)) ** 2
fig, axis = plt.subplots(1)
axis.scatter(df5.mean(axis=1), qv2, s=5, alpha=0.5)
axis.axhline(0.25, color="black", alpha=0.2, linestyle="--", zorder=0)
axis.axvline(df5.mean(axis=1).quantile(0.1), color="black", alpha=0.2, linestyle="--", zorder=0)
axis.set_xlabel("Mean (log PPM)")
axis.set_ylabel("QV2")
sns.despine(fig)
fig.savefig(os.path.join("results", "ip-lcmsms.ppm.log2.variable.scatter.svg"), dpi=300)

variable_prots = qv2[(qv2 > 0.25) & (df5.mean(axis=1) > df5.mean(axis=1).quantile(0.1))].index.unique()
g = sns.clustermap(df5.loc[variable_prots], cmap="RdBu_r", robust=True, rasterized=True, metric="correlation", cbar_kws={"label": "Expression (log2)"})
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=2, family="arial")
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
g.savefig(os.path.join("results", "ip-lcmsms.ppm.log2.variable.clustermap.svg"), dpi=300)

g = sns.clustermap(df5.loc[variable_prots], cmap="RdBu_r", robust=0, rasterized=True, z_score=0, metric="correlation", cbar_kws={"label": "Expression (Z-score)"})
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=2, family="arial")
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
g.savefig(os.path.join("results", "ip-lcmsms.ppm.log2.variable.clustermap.zscore.svg"), dpi=300)
