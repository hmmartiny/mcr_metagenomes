import pandas as pd
import numpy as np
import sys
import re
import pycountry
import geopy
import os
import reverse_geocoder as rg

import cartopy.io.shapereader as shpreader
from cartopy.feature import BORDERS, ShapelyFeature
import cartopy.crs as ccrs
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

parpath = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), os.path.pardir
        )
    )

sys.path.append(
    parpath
)


from src.dataviz.dataviz import norm_cmap, make_colorbar


class CountryLookup:
    """
    A class that tries to correct country labels found in the original metadata set from ENA

    Examples
    ----------
        >>> country_lookup = CountryLookup(no_fix='NOFIX') 
        >>> df['country_name'] = df['country'].apply(country_lookup.search)
    """
    
    def __init__(self, na_val='Not available', no_fix=np.nan, data_dir='data/'):
        """Initialize the lookup system

        Parameters
        ----------
        na_val : str, optional
            Label for missing information, by default 'Not available'
        no_fix : str, optional
            Label for search that did not yield a result, by default None
        """
        self.na_val = na_val
        self.no_fix = no_fix
        self.data_dir = data_dir
        
        if self.no_fix is None:
            self.no_fix = self.na_val
        
        self.chars = re.compile(r'[\.\,\:\;\&\'\"\?\{\}\\\*\=]') # that looks weird..
        self.nums = re.compile(r'^[\d-]')
        self.coord_1 = re.compile(r'(\d+\.\d+)\s([NSEW])')
        self.coord_2 = re.compile(r'(\d+\.\d+)')
        
        # Load geometry shapes
        self.load_world()

        # Load coordinates for countries
        self.load_coords_country()

        # alternative spellings for countries and oceans
        self.altGeoNames = {
            'USA': 'United States Of America',
            'North Macedonia': 'Macedonia',
            'Tanzania': 'United Republic Of Tanzania',
            'Czech Republic': 'Czechia',
            'Italia': 'Italy',
            'Not available': self.na_val,
            'NULL': self.na_val,
            None: self.na_val
        }
        
    def search(self, country):
        """
        Wrapper function that combines all search steps for a country string.

        The search is broken down into three steps:
            1. First the raw string is tested to see if a match can be found.
            2. Then all weird characters (e.g. '?' or ':') are removed and step 1 search is repeated.
            3. Loop through the indiviual words in the input string with weird characters removed.
        
        If nothing is found, then the `no_fix` label is found. 
        If the value is nan, the `na_val` label is returned.

        Parameters
        ----------
        country : str
            Input country string to be fixed

        Returns
        -------
        country_name : str
            Either the correct label for the country or one of the labels for missing results
        """
        
        if self.isnan(country):
            return self.na_val
        
        # broken down into multiple step functions that essentially just calls step_1 over and over again.
        steps = ['step_1', 'step_2', 'step_3']
        for step in steps:
            func = getattr(self, step)
            
            fixed, fixed_val = func(country)
        
            if fixed:
                return fixed_val
                
        return self.no_fix
        
    
    def step_1(self, country):
        """No string editing applied"""
        
        fixed, fixed_val = False, self.no_fix        

        # check format of string
        # if coordinate, try that.
        if self.nums.search(country):
            fixed, fixed_val = self.coordinate_lookup(country)            
        
        # if string, look up.
        elif isinstance(country, str):
            fixed, fixed_val = self.pycountry_lookup(country)
            
            if not fixed:
                fixed, fixed_val = self.pycountry_lookup(country, funcname='search_fuzzy')
        
        return fixed, fixed_val
        
        
    def step_2(self, country):
        """Replace weird characters"""        
        country = re.sub(self.chars, ' ', country)
        if self.isnan(country):
            return True, self.na_val
        
        fixed, fixed_val = self.step_1(country)
        
        return fixed, fixed_val
    
    def step_3(self, country):
        """Replace weird characters and search on each individual word in string"""  
        cs = re.split(self.chars, country)
        
        for c in cs:
            fixed, fixed_val = self.step_1(c.strip())
            
            if fixed:
                return fixed, fixed_val
        
        return False, country
    
    def isnan(self, string):
        """Check if the string is empty"""
        return pd.isna(string) or string.strip() == ''
    
    def pycountry_lookup(self, country, funcname='lookup'):
        """Utilizes pycountry to lookup if the country string is correct.

        Parameters
        ----------
        country : str
            input country string
        funcname : str, optional
            Either search exact ('lookup') or do a fuzzy search ('search_fuzzy'), by default 'lookup'.
            If fuzzy search is applied, the first match is returned.

        """
        func = getattr(pycountry.countries, funcname)
        try:
            m = func(country)
            if isinstance(m, list):
                m = m[0]
            return True, m.name
        except:
            return False, country
        
    def to_coordinate(self, coord, verbose=False):
        """Convert a string to DDM coordinate"""
        
        m = self.coord_1.findall(coord)
        coords = []
        if len(m) == 2:           
            for mm in m:
                dd, direction = mm
                dd = float(dd)
                if direction in ('S', 'W'):
                    dd *= -1
                coords.append(dd)
        else:
            m = self.coord_2.findall(coord)
            coords = [float(mm) for mm in m]
            
        if verbose:
            print("Input:", coord)
            print("re findall:", m)
            print("length of m:", len(m))
            print("Returned:", coords)
        
        if len(coords) < 2:
            coords = [np.nan, np.nan]
        
        return coords
    
    def load_geo(self, shp_fname=None, name=None, resolution='10m', category='physical'):
        """Use cartopy reader to load NaturalEarthData files"""
        if name is not None:
            shp_fname = shpreader.natural_earth(
                resolution=resolution,
                category=category, 
                name=name
            )

        reader = shpreader.Reader(shp_fname)
            
        geom = unary_union(list(reader.geometries()))
        area = prep(geom)
        
        recres = {}
        for rec in reader.records():
            try:
                ckey = rec.attributes['name']
            except:
                ckey = rec.attributes['ADMIN']
                
            if ckey is None:
                ckey = 'None'
            recres[ckey.title()] = {
                'prep': prep(rec.geometry),
                'geo': rec.geometry
            }
        return area, recres
    
    def load_world(self):
        """Load both land and marine labels and geometry coordinates"""
        self.land, self.land_recs = self.load_geo(
            #shp_fname = #'data/ne_10m_admin_0_map_units/ne_10m_admin_0_map_units.shp'
            shp_fname = os.path.join(self.data_dir, 'ne_10m_admin_0_countries.shp')
        )
        self.ocean, self.ocean_recs = self.load_geo(
            shp_fname = os.path.join(self.data_dir, 'ne_10m_geography_marine_polys.shp')
        )

        self.lakes, self.lake_recs = self.load_geo(
            shp_fname = os.path.join(self.data_dir, 'ne_10m_lakes.shp')
        )

        self.area2geo = {**self.land_recs, **self.ocean_recs, **self.lake_recs}

    def is_land(self, x, y):
        """Check if a given coordinate is on land"""
        return self.land.contains(sgeom.Point(x, y))
    def is_ocean(self, x, y):
        """Check if a given coordinate is in the ocean"""
        return self.ocean.contains(sgeom.Point(x, y))
    
    def is_lake(self, x, y):
        """Check if a given coordinate is in a lake"""
        return self.lakes.contains(sgeom.Point(x, y))
    
    def geo_match(self, x, y, records):
        """Loop through defined areas and see if a given coordinate falls within"""
        for geoName, geometry in records.items():
            match = geometry['prep'].contains(sgeom.Point(x, y))
            if match:
                return match, geoName
        
        return False, None

    def name2geo(self, geoName, return_category=True):
        """
        Convert a country or ocean name into a geographical shape.
        The geographical shape can then be used to create maps.
        
        """
        geo = np.nan
        geotype = np.nan

        if pd.isna(geoName) or not isinstance(geoName, str):
            return [geo, geotype]
        
        oceans = ['North', 'South']
        if geoName == 'Pacific Ocean' or geoName == 'Atlantic Ocean':
            geo = [
                self.area2geo[d + ' ' + geoName]['geo'] for d in oceans
            ]
            geotype='Water'
        else:

            try:
                geo = [self.area2geo[geoName.title()]['geo']]
            except: 
                pass
        
        if geoName.title() in self.land_recs: 
            geotype='Land'
        elif geoName.title() in self.ocean_recs or geoName.title() in self.lake_recs: 
            geotype='Water'
        
        if return_category:
            return [geo, geotype]
        else:
            return geo
    
    def reverse_geo(self, x, y):
        fixed, fixed_val = False, self.no_fix

        hits = rg.search((x, y))
        if len(hits) == 1:
            fixed_val = hits[0]['name']
            fixed = True
        
        return fixed, fixed_val

    def coordinate_lookup(self, country):
        fixed, fixed_val = False, self.no_fix
        
        if country.count('.') == 2:
            x, y = self.to_coordinate(country)

            fixed, fixed_val = self.reverse_geo(x, y)

            if fixed:
                return fixed, fixed_val
            
            if self.is_lake(x, y): # lakes go first, because they can be within country borders.
                fixed = fixed_val = self.geo_match(x,y, self.lake_recs)
                            
            elif self.is_land(x, y):
                fixed, fixed_val = self.geo_match(x, y, self.land_recs)
                
            elif self.is_ocean(x, y):
                fixed, fixed_val = self.geo_match(x, y, self.ocean_recs)
            
        return fixed, fixed_val
    
    def add_water_borders(self, ax):
        
        for k, v in {**self.ocean_recs, **self.lake_recs}.items():
            if len(v) > 1: continue

            geo = v['geo']
            ax.add_feature(
                ShapelyFeature(
                    [geo],
                    ccrs.PlateCarree(),
                    facecolor='white',
                    edgecolor='grey',
                    linestyle=':'
                )
            )
        
    def gps_scatter(self, df, lon, lat, valcol, ax, scale=50, fmt=1):
        ax.scatter(
                df.longitude,
                df.latitude,
                s=df[valcol].round(fmt).values*scale,
                color=df['color'],
                edgecolors='black',
                alpha=.5,
                transform=ccrs.PlateCarree(),
                zorder=2
            )

    def cartopy_map(self, df, geocol, valcol, ax_map, ax_cbar=None, ncmap='Blues', cbar_label='', plot_water=True, vmin=None, vmax=None, scale=50, **plot_args):
        
        # make sure that we are not modifying the original dataframe
        df = df.copy()
        
        df.dropna(subset=[valcol], inplace=True)

        # setup colors
        cmap, norm = norm_cmap(df[valcol], cmap=ncmap, vmin=vmin, vmax=vmax)
        df['color'] = df[valcol].apply(cmap.to_rgba)

        # style shapes
        shape_style = {
            'Water': {
                'edgecolor': 'lightgray',
                'linestyle': ':',
                'lw': 1
            },
            'Land': {
                'edgecolor': 'black',
                'linestyle': '-',
                'lw': .5
            }
        }
        
        # create map        
        for _, row in df.loc[df['geotype'] != 'Water'].iterrows():    
            geos = row[geocol]
            if not isinstance(geos, list):
                geos = [geos]

            for geo in geos: # to handle oceans
                try:
                    ax_map.add_feature(
                        ShapelyFeature(
                            [geo],
                            ccrs.PlateCarree(),
                            facecolor=row['color'],
                            **shape_style[row['geotype']],
                            **plot_args
                        )
                    )
                except KeyError:
                    pass
        
        if plot_water:
            ax_map.set_global()
            self.add_water_borders(ax=ax_map)
            
            water_df = df.loc[df['geotype'] == 'Water'].merge(
                self.country_coords, left_on='country', right_on='name'
            )

            self.gps_scatter(water_df, lat='lat', lon='lon', valcol=valcol, ax=ax_map, scale=scale)

        # pretty map axis
        ax_map.background_patch.set_visible(False)
        ax_map.outline_patch.set_visible(False)
        
        # add borders and coastlines
        ax_map.coastlines(lw=.5)
        ax_map.add_feature(BORDERS, linestyle='-', edgecolor='black', lw=.5)

        # make axis global

        if ax_cbar is not None:
            make_colorbar(ax_cbar, cmap.cmap, norm, label=cbar_label)


    def load_coords_country(self):
        self.country_coords = pd.read_csv(
            os.path.join(
                self.data_dir, 'countries.csv'
            )
        )

    def name2coordinate(self, geoName):
        """Return coordinates for a given country

        Parameters
        ----------
        geoName : str
            country name

        Returns
        -------
        float, float
            latitude, longitude
        """

        if geoName in ['USA', 'United States Of America']:
            geoName = 'United States'

        m = self.country_coords.loc[self.country_coords['name'] == geoName]
        try:
            lat = m.latitude.item()
            lon = m.longitude.item()
        except ValueError:
            lat, lon = np.nan, np.nan

        return lat, lon

    def name2region(self, geoName):
        if geoName in ['USA', 'United States of America']:
            geoName = 'United States'
        
        m = self.country_coords.loc[self.country_coords['name'] == geoName]
        if m.shape[0] > 0:
            return m.country.values[0]
        else:
            return np.nan
    
    def coordinator(self, row):
        country = row.country
        location = row.location
        
        lat, lon = np.nan, np.nan
        # FIRST CASE: location is given
        if isinstance(location, str):
            lat, lon = self.to_coordinate(location)
        # SECOND CASE: Country is given, but no coordinates
        elif np.isnan(location) and isinstance(country, str):
            lat, lon = self.name2coordinate(country)
        
        return lat, lon        

if __name__ == "__main__":
    locator = CountryLookup()