#!/usr/bin/env python3
import sys
import pandas as pd
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
import matplotlib as mpl

sys.path.append(
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), os.path.pardir
        )
    )
)
from country_lookup import CountryLookup
countryLocator = CountryLookup()

def plot_choropleth(df, column, left_on='country_name', cbar_label='', figsize=(10, 7), title='', how='inner', ax=None, vmin=None, vmax=None, na_val=None, **args):
    """Not in use anymore, should use the function `cartography_map` in the CountryLookup class"""

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


def plot_maps(cm, cm_tot, left_on='country', figsize=(20, 14), ncols=4, nrows=5, cbar_height=.2, subtitles_kwargs={}, ytot='Total ALR', title_kwargs={}, plot_args={}):
    
    vmax = max(cm.max().max().item(), cm_tot.max().max().item())
    vmin = max(cm.min().min().item(), cm_tot.min().min().item())
    
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
        # plot_choropleth(
        #     cm[gene],
        #     column=gene,
        #     ax=ax,
        #     vmin=vmin,
        #     vmax=vmax,
        #     legend=False,
        #     cmap='Blues',
        #     left_on=left_on
        # )

        countryLocator.cartopy_map(
            df = cm,
            valcol = gene,
            geocol = 'geo',
            ax_map = ax,
            plot_ocean=plot_ocean,
            **plot_args
        )
        
        ax.axis('off')
        ax.set_title(gene, **subtitles_kwargs)
        
        
    # plot total levels in the whole wide world
    # plot_choropleth(
    #     cm_tot, 
    #     column=ytot, 
    #     cbar_label='ALR', 
    #     ax=ax_main,
    #     vmin=vmin,
    #     vmax=vmax,
    #     legend=True,
    #     cax=ax_bar,
    #     cmap='Blues',
    #     left_on=left_on
    # )

    countryLocator.cartopy_map(
        df = cm_tot,
        geocol = 'geo',
        valcol = ytot,
        ax_map = ax_main,
        ax_cbar = ax_bar,
        plot_ocean=plot_ocean,
        cbar_label=ytot
        **plot_args

    )
    
    # set title for big map
    ax_main.set_title(ytot, **title_kwargs)
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

import cartopy.crs as ccrs
from cartopy.feature import BORDERS, ShapelyFeature

def make_colorbar(ax, cmap, norm, orientation='horizontal', label=''):
    cb1 = mpl.colorbar.ColorbarBase(
        ax=ax,
        cmap=cmap,
        norm=norm,
        orientation=orientation
    )
    
    cb1.set_label(label)

def cartopy_map(df, geocol, valcol, ax_map, ax_cbar=None, ncmap='Blues', cbar_label=''):
    
    df = df.copy()
    
    # setup colors
    cmap, norm = norm_cmap(df[valcol], cmap=ncmap)
    df['color'] = df[valcol].apply(cmap.to_rgba)
    
    # plot map
    ax_map.coastlines()
    ax_map.add_feature(BORDERS, linestyle=':')
    
    for _, row in df.iterrows():
        ax_map.add_feature(
            ShapelyFeature(
                [row[geocol]],
                ccrs.PlateCarree(),
                facecolor=row['color']
            )
        )
    # pretty map axis
    ax_map.background_patch.set_visible(False)
    ax_map.outline_patch.set_visible(False)

    if ax_cbar is not None:
        make_colorbar(ax_cbar, cmap.cmap, norm, label=cbar_label)