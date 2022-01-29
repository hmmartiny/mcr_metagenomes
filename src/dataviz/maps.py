#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
import matplotlib as mpl
import cartopy.crs as ccrs

import argparse
import gc

parpath = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), os.path.pardir
        )
    )
sys.path.append(
    parpath
)
from country_lookup import CountryLookup
from src.dataviz.dataviz import gene_formatter
from funcs import alr, closure,count_matrix, tot_alr

# countryLocator = CountryLookup(data_dir=os.path.join(parpath, os.path.pardir, 'data', 'country_shapes'))

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

def get_geo(cm, countryLocator, col='country'):

    cm.reset_index(inplace=True)
    cm[['geo', 'geotype']] = cm[col].apply(
        lambda x: countryLocator.name2geo(x) if isinstance(x, str) else np.nan
    ).apply(pd.Series)

    return cm

def plot_maps(cm, cm_tot, countryLocator, left_on='country', figsize=(20, 14), ncols=4, nrows=5, cbar_height=.2, subtitles_kwargs={}, ytot='Total ALR', small_scale=10, title_kwargs={}, plot_args={}, plot_water=True):
    
    vmax = max(cm.max().max().item(), cm_tot.max().max().item())
    vmin = max(cm.min().min().item(), cm_tot.min().min().item())

    # get genes
    genes = sorted(cm.columns.tolist())
    genes[ncols:] = genes[ncols:][::-1]

    # get geo
    cm = get_geo(cm, countryLocator, col = left_on)
    cm_tot = get_geo(cm_tot, countryLocator, col=left_on)
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=ncols, nrows=nrows, height_ratios=np.repeat(1, nrows-1).tolist() + [cbar_height])

    ax_main = fig.add_subplot(gs[:-2, :-1], projection=ccrs.PlateCarree())
    # ax_main.set_aspect('auto') 
    ax_bar = fig.add_subplot(gs[-1, :])
    
    axes = []
    # add subplots in second to last row
    for i in range(nrows-1):
        ax = fig.add_subplot(gs[-2, i], projection=ccrs.PlateCarree())
        # ax.set_aspect('auto')
        axes.append(ax)
    
    # add remaining subplots in the last column
    for i in range(ncols-1):
        ax = fig.add_subplot(gs[i, -1], projection=ccrs.PlateCarree())
        ax.set_aspect('auto')
        axes.append(ax)
    
    # start plotting the small maps
    for gene, ax in zip(genes, axes):
        countryLocator.cartopy_map(
            df = cm,
            valcol = gene,
            geocol = 'geo',
            ax_map = ax,
            plot_water=plot_water,
            vmax=vmax,
            vmin=vmin,
            scale=small_scale,
            **plot_args
        )
        
        ax.axis('off')
        ax.set_title(gene_formatter(gene), **subtitles_kwargs)
        
        
    # plot total levels in the whole wide world
    countryLocator.cartopy_map(
        df = cm_tot,
        geocol = 'geo',
        valcol = ytot,
        ax_map = ax_main,
        ax_cbar = ax_bar,
        plot_water=plot_water,
        cbar_label='ALR',
        vmax=vmax,
        vmin=vmin,
        **plot_args

    )
    
    # set title for big map
    ax_main.set_title(ytot, **title_kwargs)
    ax_main.axis('off')
    
    plt.close(fig)
    return fig

def plot_single_map(cm, valcol, countryLocator, countrycol='country', figsize=(20, 12), cbar_height=.2, plot_water=True, vmin=None, vmax=None):

    cm = cm.copy()
    if cm[valcol].isna().any():
        cm.dropna(subset=[valcol])

    # get geo
    cm = get_geo(cm, countryLocator, col=countrycol)

    # define fig
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[1, cbar_height])

    ax_main = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    # ax_main.set_aspect('auto') 
    ax_bar = fig.add_subplot(gs[1])

    # plot total levels in the whole wide world
    countryLocator.cartopy_map(
        df = cm,
        geocol = 'geo',
        valcol = valcol,
        ax_map = ax_main,
        ax_cbar = ax_bar,
        plot_water=plot_water,
        cbar_label='ALR',
        vmax=vmax,
        vmin=vmin,
    )

    ax_main.set_title(valcol)

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

def make_maps(df, dest, countryLocator, context='paper'):
    mdf = df.melt(
        ['run_accession', 'host', 'country','collection_date', 'collection_year', 'bacterial_fragment', 'tax_id'], 
        var_name='gene', 
        value_name='fragmentCountAln'
    )
    gc.collect()
    del df

    cm_country = count_matrix(
        df = mdf.query("collection_year >= 2000").replace(['NULL', 'not applicable', 'not available'], pd.NA),
        groupby_col=['country', 'gene'],
        index_pivot = 'country',
        column_pivot = 'gene'
        ).merge(
            mdf.groupby('country').agg(
                {'bacterial_fragment': 'sum'}
            ) / 1e6,
            right_index=True, left_index=True
    )

    cm_country.rename(index={
        'USA':'United States of America',
        'North Macedonia': 'Macedonia',
        'Tanzania': 'United Republic of Tanzania',
        'Czech Republic': 'Czechia'
    }, inplace=True)


    cm_country_alr = cm_country.apply(alr, axis=1)#.replace([np.inf, -np.inf])

    gc.collect()

    cm_country_alr_total = pd.DataFrame(cm_country.apply(tot_alr, axis=1))
    cm_country_alr_total.rename(columns={0: 'Total ALR'}, inplace=True)

    cm_country_alr_total.reset_index(inplace=True)
    cm_country_alr_total[['geo', 'geotype']] = cm_country_alr_total['country'].apply(
        lambda x: countryLocator.name2geo(x) if not isinstance(x, float) else pd.NA
    ).apply(pd.Series)
    cm_country_alr_total.set_index('country', inplace=True)

    gc.collect()
    
    cm_country_alr.reset_index(inplace=True)
    cm_country_alr[['geo', 'geotype']] = cm_country_alr['country'].apply(
        lambda x: countryLocator.name2geo(x) if not isinstance(x, float) else pd.NA
    ).apply(pd.Series)
    cm_country_alr.set_index('country', inplace=True)

    gc.collect()
    

    del cm_country

    maps = plot_maps(
        cm = cm_country_alr.select_dtypes('float'),
        cm_tot = cm_country_alr_total.select_dtypes('float'),
        countryLocator = countryLocator,
        ncols=5,nrows=6,
        left_on='country'
    )

    gc.collect()
    del countryLocator
    del cm_country_alr
    del cm_country_alr_total

    # maps.savefig(os.path.join(dest, f'maps_{context.upper()}.pdf'))
    maps.savefig(os.path.join(dest, f'maps_{context.upper()}.png'))

    gc.collect()
    del maps

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-d', '--dest',
        type=str,
        required=True,
        help='destionation folder to save figures in',
        dest='dest'
    )

    parser.add_argument(
        '-f', '--file',
        type=str,
        required=True,
        help='mcr file',
        dest='file'
    )

    return parser.parse_args()
if __name__ == '__main__':

    gc.collect()
    args = parse_args()

    df = pd.read_csv(args.file, parse_dates=['collection_date'])
    gc.collect()
    if not 'collection_year' in df.columns:
        df['collection_year'] = df['collection_date'].dt.year
    gc.collect()

    countryLocator = CountryLookup(data_dir=os.path.join(parpath, os.path.pardir, 'data', 'country_shapes'))    
    make_maps(df, dest=args.dest, countryLocator=countryLocator)

    gc.collect()
    del countryLocator