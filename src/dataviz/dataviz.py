#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib as mpl
from matplotlib import gridspec 

def barplot(data, x, y, xlabel, ylabel, title='', orient='h', figsize=(10,7)):

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if orient == 'h':
        x, y = y, x # switch
        xlabel, ylabel = ylabel, xlabel

    sns.barplot(
        data=data, x=x, y=y, orient=orient, ax=ax
    )

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)

    fig.tight_layout()

    plt.close(fig)

    return fig

def plot_choropleth(df, column, left_on='country_name', cbar_label='', figsize=(10, 7), title='', how='inner', ax=None, vmin=None, vmax=None, na_val=None, **args):

    df = df.copy()

    # replace usaworld
    df.replace('USA', 'United States of America', inplace=True)

    world = gpd.read_file('/Users/hanmar/repos/db_work/data/ne_10m_admin_0_map_units/ne_10m_admin_0_map_units.shp')#gpd.datasets.get_path('naturalearth_lowres'))
    #world = world.loc[world.TYPE.isin(['Sovereign country', 'Country'])]
    gpd_df = gpd.GeoDataFrame(
        world.merge(df, left_on='NAME', right_on=left_on, how=how),
        geometry='geometry',
    )            

    # fill nan
    if na_val is not None:
        gpd_df[column].fillna(na_val, inplace=True)

    return_fig = False
    if ax is None:
        return_fig = True
        fig, ax = plt.subplots(figsize=figsize)

    # borders
    world.plot(
        color='white',
        edgecolor='black',
        lw=1,
        ax=ax,
        alpha=.2
    )

    # colors
    gpd_df.plot(
        column=column, 
        ax=ax, 
        legend_kwds={
            'orientation': 'horizontal', 
            'label': cbar_label,
            'fraction': 0.046,
        },
        vmin=vmin,
        vmax=vmax,
        **args
    )

    ax.set_title(title)

    if return_fig: 
        fig.tight_layout()
        plt.close(fig)
        return fig
    

def swarm_plot(df, x, y, hue, xlabel, ylabel, alpha=.5, figsize=(15, 10), title='', rotation=90, ax=None, orientation='v'):

    hue_order = sorted(df[hue].unique().tolist())

    return_fig = False
    if ax is None:
        return_fig = True
        fig, ax = plt.subplots(1,1, figsize=figsize)
    
    if orientation == 'h':
        x, y = y, x # switch
        xlabel, ylabel = ylabel, xlabel

    sns.catplot(
        data=df.sort_values(by=[x, hue]), x=x, y=y, 
        hue=hue, hue_order=hue_order,
        kind='swarm', alpha=alpha, ax=ax,
        orient=orientation
    )

    # rotate xaxis labels
    for label in ax.get_xticklabels():
        label.set_rotation(rotation)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if return_fig:
        fig.tight_layout()
        plt.close(fig)
        return fig


def closure_barplot(df, stacked=True, figsize=(10, 7), ax=None, title='', xlabel='', ylabel='', add_sum=False, add_sum_name='Total', replace_zero=True):

    df = df.copy()

    if add_sum:
        df.loc[add_sum_name, :] = df.sum()
    
    if replace_zero:
        df = df.coda.zero_replacement()
        
    df = df.coda.closure(100)
    
    return_fig = False
    if ax is None:
        return_fig = True
        fig, ax = plt.subplots(1,1, figsize=figsize)

    df.plot.bar(
        ax=ax,
        stacked=stacked,
        legend=False
    )
    ax.legend(bbox_to_anchor=(1, 0.5))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    if return_fig:
        fig.tight_layout()
        plt.close(fig)
        return fig

# def year_bar_plot(df, year_col, pivot_cols, ncols=2, subtitle='', as_int=True, xlabel='', ylabel=''):
#     years = sorted(df[year_col].dropna().unique().tolist())

#     if as_int:
#         years = [int(x) for x in years]
    
#     nrows = len(years) // ncols

#     fig, axes = plt.subplots(nrows, ncols, figsize=(10*ncols, 7*nrows))

#     for year, ax in zip(years, axes.flatten()):
#         sub_df = df.loc[df[year_col] == year].copy()

#         sub_cmat = count_matrix(
#             sub_df, 
#             groupby_col=pivot_cols,
#             index_pivot=pivot_cols[0],
#             column_pivot=pivot_cols[1]
#         )

#         closure_barplot(
#             sub_cmat, 
#             ax=ax,
#             title=f'Distribution of {subtitle} in {year}',
#             xlabel=xlabel,
#             ylabel=ylabel
#         )

#     fig.tight_layout()
#     plt.close(fig)
    
#     return fig

def gene_formatter(g, split='-'):
    s = r"${}$" + split + "{}"
    return s.format(*g.split(split))

def make_timeline(ax, dates, labels):
    
    levels = np.repeat(1, len(dates))
    
    markerline, _, _ = ax.stem(
        dates, levels,
        #linefmt="C3-",
        basefmt='k-',
        use_line_collection=True
    )
    
    markerline.set_ydata(np.zeros(len(dates)))

    vert = np.array(['top', 'bottom'])[(levels >0).astype(int)]
    for d, l, r, va in zip(dates, levels, labels, vert):
        ax.annotate(
            r, xy=(d, l), xytext=(-3, np.sign(l)*2),
            textcoords="offset points",
            va=va,
            ha='center'
        )

    ax.get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.margins(y=0.1)

def make_stacked_bars(ax, data, labels, index, width):
    
    ax.bar(index, data[labels[0]], width, label=labels[0])
    bottom = np.array(data[labels[0]])
    
    i = 1
    while i < len(labels):
        ax.bar(index, data[labels[i]], width, label=labels[i], bottom=bottom)
        bottom += np.array(data[labels[i]])
        i += 1

def make_bars(ax, data, x, ys, width, colors, alpha=.5):
    shift=0
    for y, c in zip(ys, colors):
        ax.bar(
            data[x] + shift,
            data[y],
            width = width,
            color=c,
            label=y,
            alpha=alpha
        )
        
        shift += width

def plot_timeline(cm, cm_closed, gene_timeline, timesamples, samples_ys, figsize=(12, 7), legend_kwargs={}, ylabel_kwargs={}):
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=1, nrows=3, height_ratios=[.15, 1, .5])

    ax_timeline = fig.add_subplot(gs[0])
    ax_alr = fig.add_subplot(gs[1], sharex=ax_timeline)
    ax_samples = ax_alr.twinx()
    ax_bar = fig.add_subplot(gs[2], sharex=ax_timeline)
    
    # plot the timeline according to the litterature
    timeline = gene_timeline.groupby('Year discovered')['Gene'].apply(
        lambda x: "\n".join([gene_formatter(g) for g in x])
    )
    years = timeline.index.astype('int').tolist()
    labels = timeline.values.tolist()

    make_timeline(
        ax = ax_timeline,
        dates = years,
        labels = labels
    )

    # plot the count of samples with and without gene hits    
    make_bars(
        ax = ax_samples,
        data = timesamples,
        x = 'year',
        ys = samples_ys,
        width = .2,
        colors=['dimgray', 'lightgray']
    )
    
    # plot ALR data
    for gene in cm:
        ax_alr.plot(
            cm.index,
            cm[gene],
            '-o',
            label=gene_formatter(gene),
            lw = 1.5,
            alpha = .8
        )
        
    # plot relative abundance(s) in the composition
    make_stacked_bars(
        ax = ax_bar,
        data = cm_closed,
        labels = cm_closed.columns,
        index = cm_closed.index,
        width = .5
    )
    
    # show all xticks
    xticks = sorted(list(set(cm_closed.index.tolist() + years)))
    xticks = np.arange(min(xticks), max(xticks) + 1, 1) # no gaps
    for ax in [ax_timeline, ax_alr, ax_bar]:
        ax.set_xticks(xticks)
    
    # hide axis on middle plot
    ax_alr.get_xaxis().set_visible(False)
    
    # create figure legends
    ax_alr.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., **legend_kwargs)
    ax_samples.legend(bbox_to_anchor=(1.05, 0.35), loc=2, borderaxespad=0., **legend_kwargs)
    
    # create labels for the y-axes
    ax_alr.set_ylabel('ALR', **ylabel_kwargs)
    ax_bar.set_ylabel('Relative abundance', **ylabel_kwargs)
    ax_samples.set_ylabel(r"$\log_{10}$(Count) (bars)", **ylabel_kwargs)
    
    # add title for discovery
    ax_timeline.text(xticks[0], .1, "Timeline of reporting", ylabel_kwargs)

    plt.close(fig)
    return fig

def plot_hosts(cm_tot, cm_closed, host_counts, samples_ys, sample_colors, sort_by='Samples\nwithout $mcr$', figsize=(12, 7), legend_kwargs={}, ylabel_kwargs={}):
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[.5, 1])

    ax_tot = fig.add_subplot(gs[0])
    ax_ab = fig.add_subplot(gs[1], sharex=ax_tot)
    
    # merge
    bar_data = host_counts.merge(
        cm_tot, left_index=True, right_index=True
    ).sort_values(by=sort_by, ascending=False)
    bar_data['host_name'] = bar_data.index
    host_order = bar_data.index.tolist()
    
    # plot sample counts and total ALR levels for each host
    bar_data.plot.bar(
        x = 'host_name',
        y = samples_ys,
        color = [sample_colors[y] for y in samples_ys],
        ax = ax_tot,
        width=.7
    )
    
    # plot relative levels in the composition for each host
    cm_closed.loc[host_order, :].plot.bar(
        stacked=True,
        ax = ax_ab
    )
    
    # make figure legends (look nice)
    ax_tot.legend(bbox_to_anchor=(1,1), **legend_kwargs)
    
    bar_handles, bar_labels = ax_ab.get_legend_handles_labels()
    bar_labels = [gene_formatter(t) for t in bar_labels]
    ax_ab.legend(
        handles=bar_handles,
        labels=bar_labels,
        bbox_to_anchor=(1, 1),
        **legend_kwargs
    )
    
    # set ylabels
    ax_tot.set_ylabel('$\log_{10}$(count)', **ylabel_kwargs)
    ax_ab.set_ylabel('Relative Abundance', **ylabel_kwargs)
    
    # remove xlabel
    ax_ab.set_xlabel('')
    
    plt.close(fig)
    return fig

def plot_maps(cm, cm_tot, left_on='country', figsize=(20, 14), ncols=4, nrows=5, cbar_height=.2, subtitles_kwargs={}, ytot='Total ALR', titel_kwargs={}):
    
    vmax = max(cm.max().max().item(), cm_tot[ytot].max().max().item())
    vmin = max(cm.min().min().item(), cm_tot[ytot].min().min().item())
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=ncols, nrows=nrows, height_ratios=np.repeat(1, nrows-1).tolist() + [cbar_height])

    ax_main = fig.add_subplot(gs[:-2, :-1]) 
    ax_bar = fig.add_subplot(gs[-1, :])
    
    axes = []
    # add subplots in second to last row
    for i in range(nrows-1):
        ax = fig.add_subplot(gs[-2, i])
        axes.append(ax)
    
    # add remaining subplots in the last column
    for i in range(ncols-1):
        ax = fig.add_subplot(gs[i, -1])
        axes.append(ax)
    
    # start plotting the small maps
    genes = sorted(cm.columns.tolist())
    for gene, ax in zip(genes, axes):
        plot_choropleth(
            cm[gene],
            column=gene,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            legend=False,
            cmap='Blues',
            left_on=left_on
        )
        
        ax.axis('off')
        ax.set_title(gene, **subtitles_kwargs)
        
        
    # plot total levels in the whole wide world
    plot_choropleth(
        cm_tot, 
        column=ytot, 
        cbar_label='ALR', 
        ax=ax_main,
        vmin=vmin,
        vmax=vmax,
        legend=True,
        cax=ax_bar,
        cmap='Blues',
        left_on=left_on

    )
    
    # set title for big map
    ax_main.set_title(ytot, **titel_kwargs)
    ax_main.axis('off')
    
    plt.close(fig)
    return fig
    
def norm_cmap(values, cmap, vmin=None, vmax=None):
    """
    Normalize and set colormap
    
    Parameters
    ----------
    values : Series or array to be normalized
    cmap : matplotlib Colormap
    normalize : matplotlib.colors.Normalize
    cm : matplotlib.cm
    vmin : Minimum value of colormap. If None, uses min(values).
    vmax : Maximum value of colormap. If None, uses max(values).
    
    Returns
    -------
    n_cmap : mapping of normalized values to colormap (cmap)
    
    """
    mn = vmin or min(values)
    mx = vmax or max(values)
    norm = Normalize(vmin=mn, vmax=mx)
    n_cmap = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    return n_cmap, norm

def make_colorbar(ax, cmap, norm, orientation='horizontal', label=''):
    cb1 = mpl.colorbar.ColorbarBase(
        ax=ax,
        cmap=cmap,
        norm=norm,
        orientation=orientation
    )
    
    cb1.set_label(label)