import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster
import networkx as nx

# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


## MASSSPEC
smarc = pd.DataFrame()
for replicate in ["A", "B"]:
    smarc_ip = pd.read_csv(os.path.join("metadata", 'original', 'BRG1_{}.tab'.format(replicate)), sep="\t", index_col=[0, 1])
    smarc_ip = smarc_ip[smarc_ip.columns[~smarc_ip.columns.str.contains("Abund")]]
    smarc_ip.columns = smarc_ip.columns.str.replace(r"\..*", "")

    smarc_ip2 = smarc_ip.T.groupby(level=0).mean().T
    smarc_ip3 = smarc_ip2[smarc_ip2.index.get_level_values(1).str.contains('GN=')]
    smarc_ip3.index = list(map(lambda x: x[1].split(' ')[0], smarc_ip3.index.get_level_values(1).str.split("GN=")))
    smarc_ip3 = smarc_ip3.sort_index()
    # smarc_ip4 = smarc_ip3.dropna()
    # smarc_ip5 = np.log2(1 + (smarc_ip4 / smarc_ip4.sum(axis=0)) * 1e3)

    smarc_ip3['replicate'] = replicate
    smarc = smarc.append(smarc_ip3)

smarc = smarc.groupby(level=0).mean()
smarc.index.name = 'protein'
smarc.columns.name = 'knockout'

arid = pd.DataFrame()
for replicate in ["A", "B"]:
    arid_ip = pd.read_csv(os.path.join("metadata", 'original', 'ARID1A_{}.tab'.format(replicate)), sep="\t", index_col=[0, 1])
    arid_ip = arid_ip[arid_ip.columns[~arid_ip.columns.str.contains("Abund")]]
    arid_ip.columns = arid_ip.columns.str.replace(r"\..*", "")

    arid_ip2 = arid_ip.T.groupby(level=0).mean().T
    arid_ip3 = arid_ip2[arid_ip2.index.get_level_values(1).str.contains('GN=')]
    arid_ip3.index = list(map(lambda x: x[1].split(' ')[0], arid_ip3.index.get_level_values(1).str.split("GN=")))
    arid_ip3 = arid_ip3.sort_index()
    # arid_ip4 = arid_ip3.dropna()
    # arid_ip5 = np.log2(1 + (arid_ip4 / arid_ip4.sum(axis=0)) * 1e3)

    arid_ip3['replicate'] = replicate
    arid = arid.append(arid_ip3)

arid = arid.groupby(level=0).mean()
arid.index.name = 'protein'
arid.columns.name = 'knockout'

# smarc.index = "SMARCA4, " + smarc.index
# arid.index = "ARID1A, " + arid.index
ip_ms = arid.T.append(smarc.T)


# Plot all prots
g = sns.clustermap(arid_ip5.dropna(), cmap="RdBu_r", robust=True, rasterized=True)
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


# Get only complex members
baf_members = pd.read_csv(os.path.join("metadata", "baf_complex_subunits.csv")).squeeze()
baf_members = baf_members.drop(0)

# b = ip_ms.loc[:, pd.Series(list(map(lambda x: x[1], ip_ms.columns.str.split(", "))), index=ip_ms.columns).isin(baf_members)]
b = ip_ms.loc[:, baf_members]
b = b.loc[~b.isnull().all(axis=1), :]
b = b.loc[~(b.T.std() == 0), :]
b.columns.name = "protein"
bb = b.groupby(level=0).mean()

for m, label in [(b, "separate"), (bb, "mean")]:
    with sns.axes_style("darkgrid"):
        g = sns.clustermap(
            m.fillna(-1), mask=m.isnull(),
            robust=True, rasterized=False, metric="braycurtis", cbar_kws={"label": "Expression(log2)"}, vmin=0)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
        g.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.braycurtis.svg".format(label)), dpi=300, bbox_inches="tight")

    for metric in ['braycurtis', 'correlation']:
        with sns.axes_style("darkgrid"):
            g = sns.clustermap(
                m.T.dropna().T,
                robust=True, rasterized=False, metric=metric, cbar_kws={"label": "Expression(log2)"}, vmin=0)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
            g.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.no_nan.{}.svg".format(label, metric)), dpi=300, bbox_inches="tight")

# Make pan-BAF network

g = sns.clustermap(
    b.T.dropna().T.corr(),
    robust=True, rasterized=False, metric=metric, cbar_kws={"label": "Expression(log2)"}, vmin=0)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
g.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.BAF_all.svg".format(label)), dpi=300, bbox_inches="tight")


fig, axis = plt.subplots(1, figsize=(6, 6))
G = nx.from_pandas_adjacency((b.T.dropna().T.corr()) ** 2)
nx.draw_spectral(
    G, with_labels=True,
    node_size=10 ** (3 * b.T.dropna().T.mean()),
    node_color=((b.T.dropna().T.var() / b.T.dropna().T.mean()) ** 2), vmin=0, vmax=1, cmap="RdBu_r",
    ax=axis)
fig.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.BAF_all.networkx.svg".format("mean")), dpi=300, bbox_inches="tight")

# # q = (1 - b.T.dropna().T.corr()) ** 6
# q = (2 * b.T.dropna().T.corr()) ** 2
# q.index.name = 'knockout'
# q['mean'] = b.T.dropna().T.mean()
# nx.set_node_attributes(G, b.T.dropna().T.mean().to_dict(), name="mean")
# q['std'] = b.T.dropna().T.std()
# nx.set_node_attributes(G, b.T.dropna().T.std().to_dict(), name="std")
# q['var'] = b.T.dropna().T.var()
# nx.set_node_attributes(G, b.T.dropna().T.var().to_dict(), name="var")
# q['qv2'] = ((b.T.dropna().T.var() / b.T.dropna().T.mean()) ** 2)
# nx.set_node_attributes(G, ((b.T.dropna().T.var() / b.T.dropna().T.mean()) ** 2).to_dict(), name="qv2")
# q = pd.melt(q.reset_index(), id_vars=["knockout", 'mean', 'std', 'var', 'qv2'])
# q.to_csv(os.path.join("results", "ip-lcmsms.net.csv"), index=False)
# nx.write_gexf(G, os.path.join("results", "ip-lcmsms.net.gexf"))


# Separate the two complexes
g = sns.clustermap(
    bb.T.dropna().T,
    robust=True, rasterized=False, metric="braycurtis", cbar_kws={"label": "Expression(log2)"}, vmin=0)
labels = pd.Series(fcluster(g.dendrogram_col.linkage, t=3, criterion="maxclust"), index=bb.T.dropna().T.columns)
gg = sns.clustermap(
    bb.T.dropna().T,
    col_colors=plt.get_cmap("Paired")(labels),
    robust=True, rasterized=False, metric="braycurtis", cbar_kws={"label": "Expression(log2)"}, vmin=0)
gg.ax_heatmap.set_yticklabels(gg.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
gg.ax_heatmap.set_xticklabels(gg.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
gg.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.braycurtis.labeled.svg".format("mean")), dpi=300, bbox_inches="tight")

sub_baf = bb.loc[:, labels[labels == 1].index].columns
sub_baf = sub_baf[sub_baf.isin(bb.index)].tolist() + ['SMARCA2']
sub_pbaf = bb.loc[:, labels[labels == 2].index].columns
sub_pbaf = sub_pbaf[sub_pbaf.isin(bb.index)].tolist() + ['SMARCA2']

# BAF
baf_c = bb.loc[sub_baf, :].T.dropna().T
gg = sns.clustermap(
    baf_c.corr(),
    robust=True, rasterized=False, metric="correlation", cbar_kws={"label": "Expression(log2)"}, vmin=0)
gg.ax_heatmap.set_yticklabels(gg.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
gg.ax_heatmap.set_xticklabels(gg.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
gg.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.BAF_only.svg".format("mean")), dpi=300, bbox_inches="tight")

fig, axis = plt.subplots(1, figsize=(6, 6))
G = nx.from_pandas_adjacency((baf_c.corr() ** (1 / 4.)))
nx.draw_spectral(
    G,
    alpha=0.5, edge_cmap="greys",
    with_labels=True,
    cmap="RdBu_r", node_size=10 ** (3 * bb.loc["WT", baf_c.columns]),
    node_color=((baf_c.var() / baf_c.mean()) ** (1/4.)), vmin=0, vmax=1,
    ax=axis)
fig.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.BAF_only.networkx.svg".format("mean")), dpi=300, bbox_inches="tight")


# pBAF
pbaf_c = bb.loc[sub_pbaf, :].T.dropna().T
gg = sns.clustermap(
    pbaf_c.corr(),
    col_colors=plt.get_cmap("Paired")(labels),
    robust=True, rasterized=False, metric="correlation", cbar_kws={"label": "Expression(log2)"}, vmin=0)
gg.ax_heatmap.set_yticklabels(gg.ax_heatmap.get_yticklabels(), rotation=0, family="arial")
gg.ax_heatmap.set_xticklabels(gg.ax_heatmap.get_xticklabels(), rotation=90, family="arial")
gg.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.pBAF_only.svg".format("mean")), dpi=300, bbox_inches="tight")

fig, axis = plt.subplots(1, figsize=(6, 6))
G = nx.from_pandas_adjacency((pbaf_c.corr() ** (1 / 4.)))
nx.draw_spectral(
    G,
    alpha=0.5, edge_cmap="greys",
    with_labels=True,
    cmap="RdBu_r", node_size=10 ** (3 * bb.loc["WT", pbaf_c.columns]),
    node_color=((pbaf_c.var() / pbaf_c.mean()) ** (1/4.)), vmin=0, vmax=1,
    ax=axis)
fig.savefig(os.path.join("results", "ip-lcmsms.{}.clustermap.correlation.pBAF_only.networkx.svg".format("mean")), dpi=300, bbox_inches="tight")
