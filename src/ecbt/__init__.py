from pandas import DataFrame
from numpy import NaN
from scipy.cluster import hierarchy
from scipy.spatial import distance
from matplotlib import pyplot as plt
from matplotlib import gridspec


def eclevel(level, coln='ec', dismiss_undefined=True):
    def _(row):
        parts = row[coln].split('.')
        levelec = ".".join(parts[:level])
        if dismiss_undefined and levelec.endswith('-'):
            return NaN
        return levelec
    return _


def aggdf(btdf, level, rcoln='reaction', btcoln='bt'):
    levels = ['ec_1st', 'ec_2nd', 'ec_3rd', 'ec_4th']
    l = levels.pop(level-1)
    return btdf.groupby([btcoln,l]).nunique()\
        .loc[:,[rcoln] + levels[level-1:]]


def makeheatdf(grouped, coln='reaction'):
    
    heatdf = DataFrame()
    
    btids, ecids = grouped.index.levels
    for btn, ecn in zip(grouped.index.codes[0], grouped.index.codes[1]):
        heatdf.loc[btids[btn],ecids[ecn]] = grouped.loc[(btids[btn],ecids[ecn]), coln]
    
    heatdf = heatdf.rename(columns=str)
    heatdf = heatdf[sorted(heatdf.columns)]
    
    return heatdf


def penalize(data, penalty):
    df = data.copy()
    for i in df.index.values:
        for c in df.columns.values:
            if df.loc[i,c] != df.loc[i,c]:
                df.loc[i,c] = -penalty
    return df


def blank(data, erasable):
    df = data.copy()
    for i in df.index.values:
        for c in df.columns.values:
            if df.loc[i,c] == erasable:
                df.loc[i,c] = NaN
    return df


def clusterindex(data):
    pdist = distance.pdist(data)
    linkage = hierarchy.linkage(pdist)
    leaves = hierarchy.leaves_list(linkage)
    return [data.index.values[int(l)] for l in leaves]


def normalize(df):
    return df.div(df.sum(axis=1), axis=0)


def makeclusteredheatmap(heatdf, plot,
                         nan_penalty=100,
                         num_row_clusters=0,
                         num_col_clusters=0,
                         excluded=['-']):
    from heatmapcluster import heatmapcluster
    import matplotlib.pyplot as plt

    # convert NaN to its penalty
    datamatrix = penalize(heatdf, nan_penalty)

    # remove excluded columns
    # by default don't take missing EC numbers into account for clustering
    for coln in excluded:
        try:
            datamatrix = datamatrix.drop(coln, 1)
        except ValueError:
            pass
    
    row_labels = datamatrix.index.values
    col_labels = datamatrix.columns.values

    h = heatmapcluster(datamatrix.as_matrix(),
                       row_labels, col_labels,
                       num_row_clusters=num_row_clusters, num_col_clusters=num_col_clusters,
                       label_fontsize=9,
                       xlabel_rotation=-90,
                       cmap=plt.cm.Purples,
                       show_colorbar=False,
                       top_dendrogram=False,
                       figsize=[25,25])
    plt.savefig(plot)
    
    return [yl.get_text() for yl in h.left_dendrogram_axis.get_yticklabels()]


def addEClevels(df, coln='ec', dismiss_undefined=True):
    df['ec_1st'] = df.apply(eclevel(1, coln=coln, dismiss_undefined=dismiss_undefined), axis=1)
    df['ec_2nd'] = df.apply(eclevel(2, coln=coln, dismiss_undefined=dismiss_undefined), axis=1)
    df['ec_3rd'] = df.apply(eclevel(3, coln=coln, dismiss_undefined=dismiss_undefined), axis=1)
    df['ec_4th'] = df.apply(eclevel(4, coln=coln, dismiss_undefined=dismiss_undefined), axis=1)


def generate_heatmap(btdf, level=3, rcoln='reaction', btcoln='bt'):
    grouped = aggdf(btdf, level, rcoln, btcoln)
    return makeheatdf(grouped, coln=rcoln)


def generate_penalized_heatmap(data, newindex):
    penalty = max(data.max())
    reindexed = data.reindex(newindex)
    return penalize(reindexed, penalty=penalty)


def generate_clustered_heatmap(data):
    penalty = max(data.max())
    penalized = penalize(data, penalty=penalty)
    newindex = clusterindex(penalized)
    return newindex, penalized.reindex(newindex)


def align(df1, df2):
    missing1 = [x for x in df2.columns if x not in df1.columns]
    missing2 = [x for x in df1.columns if x not in df2.columns]
    penalty1 = max(df1.max())
    penalty2 = max(df2.max())
    df1r = df1.copy()
    for bt in df1.index:
        for ec in missing1:
            df1r.loc[bt,ec] = -penalty1
    df2r = df2.copy()
    for bt in df2.index:
        for ec in missing2:
            df2r.loc[bt,ec] = -penalty2
    return df1r[sorted(df1r.columns)], df2r[sorted(df1r.columns)]


def blank_align(df1, df2):
    missing1 = [x for x in df2.columns if x not in df1.columns]
    missing2 = [x for x in df1.columns if x not in df2.columns]
    df1r = df1.copy()
    for bt in df1.index:
        for ec in missing1:
            df1r.loc[bt,ec] = NaN
    df2r = df2.copy()
    for bt in df2.index:
        for ec in missing2:
            df2r.loc[bt,ec] = NaN
    return df1r[sorted(df1r.columns)], df2r[sorted(df1r.columns)]


def make_comparable(hm1, hm2):
    cl_ind, cdf = generate_clustered_heatmap(
        hm1.add(hm2,fill_value=0))
    cl_hm1 = generate_penalized_heatmap(hm1, cl_ind)
    cl_hm2 = generate_penalized_heatmap(hm2, cl_ind)
    return align(cl_hm1, cl_hm2)


def overplot(plt, ground, overlay, 
             labelsize=None,
             ground_cmap='autumn', overlay_cmap='Blues_r', alpha=0.5,
             title='BT rule - EC number annotation',
             ylabel='BT rules',
             xlabel='EC numbers',
             dumpfile=None):
    plt.set_title(title)
    plt.set_ylabel(ylabel)
    plt.set_xlabel(xlabel)
    plt.set_yticks(range(ground.shape[0]))
    plt.set_yticklabels(labels=ground.index.values, size=labelsize)
    plt.set_xticks(range(ground.shape[1]))
    plt.set_xticklabels(ground.columns.values, rotation=90, size=labelsize)
    
    im1 = plt.imshow(ground, cmap=ground_cmap)
    im2 = plt.imshow(overlay, cmap=overlay_cmap, alpha=alpha)


def colorscheme(ax, ground, overlay, ground_label, overlay_label, ground_cmap='Reds', overlay_cmap='Greens', alpha=0.5):
    maxx=int(ground.max().max())
    maxy=int(overlay.max().max())
    stepx = 5*(maxx//45+1)
    stepy = 5*(maxy//45+1)
    bsteps = [-maxy, 1] + list(range(stepy, maxy, stepy)) + [maxy]
    ksteps = [0, 1] + list(range(stepx, maxx, stepx)) + [maxx]
    kground = DataFrame(dict([
        (0, -maxx) ]+[
        (i, [i] * len(bsteps)) for i in ksteps[1:]
    ]))
    kground.index = [0] + bsteps[1:]
    boverlay = DataFrame(dict([
        (i, bsteps[:]) for i in ksteps
    ]))
    boverlay.index =  [0] + bsteps[1:]
    overplot(ax, kground, boverlay, ground_cmap=ground_cmap, overlay_cmap=overlay_cmap, alpha=alpha,
            title='Color Scheme', xlabel=ground_label, ylabel=overlay_label)


def combined_histogram(ground, overlay, dumpfile=None,
                       hist_h=40, hist_w=30, scheme_h=3, 
                       ground_cmap='Reds', overlay_cmap='Greens', alpha=0.5):
    fig = plt.figure(figsize=(hist_w, hist_h + scheme_h)) 
    gs = gridspec.GridSpec(2, 1, height_ratios=[hist_h, scheme_h]) 

    overplot(plt.subplot(gs[0]), ground, overlay,
             ground_cmap=ground_cmap, overlay_cmap=overlay_cmap, alpha=alpha)
    colorscheme(plt.subplot(gs[1]), ground, overlay,
                ground_label='# KEGG evidence',
                overlay_label='# EAWAG-BBD evidence',
                ground_cmap=ground_cmap, overlay_cmap=overlay_cmap, alpha=alpha)
    if dumpfile:
        plt.savefig(dumpfile)
