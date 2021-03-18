import os
import argparse
from dataviz.dataviz import plot_timeline, plot_hosts, plot_maps
from funcs import *
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from country_lookup import CountryLookup
from dataviz.maps import plot_maps

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Directory to store pdf plots in [required]',
        dest='output'
    )

    parser.add_argument(
        '-m', '--maps',
        action='store_true',
        help='If given, create maps of gene variant distribution. Default: False',
        dest='create_maps'
    )

    parser.add_argument(
        '-t', '--timeline',
        action='store_true',
        help='If given, create figure showing the gene variant distribution over time. Default: False',
        dest='create_timeline'
    )

    parser.add_argument(
        '-ho', '--host',
        action='store_true',
        help='If given, create figure showing the gene variant distribution over hosts. Default: False',
        dest='create_host'
    )

    parser.add_argument(
        '--db_args',
        type=str,
        default='/Users/hanmar/Documents/repos/mcr_metagenomes/data/db_args.json',
        help='Path to JSON file with arguments for connecting to db. Default: /Users/hanmar/Documents/repos/mcr_metagenomes/data/db_args.json',
        dest='db_args'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='If given, print statements while running code. Default: False',
        dest='verbose'
    )

    parser.add_argument(
        '--oldest_year',
        type=int,
        help='Filter results to only have samples collected after given year. Default: Discard samples where collection_year is before 1900.',
        dest='oldest_year'
    )

    return parser.parse_args()

def get_data(query, db_args, fix_gene=True, fix_date=True, oldest_year=None):
    
    df = query_db(query=query,db_args=args.db_args)

    # format genes
    if fix_gene:
        df['gene'] = df['refSequence'].apply(lambda x: x.split('.')[0].split('_')[0])
        df['gene_variant'] = df['refSequence'].apply(lambda x: x.split('_')[0] if x[5] == '.' else ".".join(x.split('_')[:2]))

    # format dates
    if fix_date:
        df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
        df.loc[df['collection_date'].dt.year > 2019, 'collection_date'] = pd.NaT
        df['collection_year'] = df['collection_date'].dt.year

        if isinstance(oldest_year, int):
            df = df.loc[df['collection_year'] >= oldest_year, ]
    
    return df


def create_cms(df, group_col, gene_col='gene', bac_col='bacterial_fragment'):

    df = df.copy()

    cm = count_matrix(
        df = df,
        groupby_col=[group_col, gene_col],
        index_pivot=group_col,
        column_pivot=gene_col
    ).merge(
        df.groupby(group_col).agg(
            {bac_col: 'sum'}
        ) / 1e6,
        right_index=True, left_index=True
    )

    # close the count matrix
    cm_closed = cm.drop(columns=bac_col).apply(closure, axis=1)

    # create the ALR matrix
    cm_alr = cm.apply(alr, axis=1).replace([-np.inf, np.inf], np.nan)

    # create total ALR matrix
    cm_alr_total = pd.DataFrame(cm.apply(tot_alr, axis=1)).replace([-np.inf, np.inf], np.nan)
    cm_alr_total.rename(columns={0: 'Total ALR'}, inplace=True)

    return cm, cm_closed, cm_alr, cm_alr_total

if __name__ == "__main__":
    args = parse_args()

    # Pull data
    if args.verbose:
        print("Retrieving data..")
    q_amr = "select * from get_amr2 where refSequence like 'mcr-%%' and bacterial_fragment > 0"
    df_main = get_data(
        query=q_amr,
        db_args=args.db_args,
        oldest_year=args.oldest_year
    )
    if args.verbose:
        print(" * Pulled {} rows matching the query:\n\t{}".format(df_main.shape[0], q_amr))
       

    if args.create_host or args.create_timeline:
        q_count = "select collection_date, host, country, count(distinct(run_accession)) as n_runs, count(distinct(sample_accession)) as n_samples, count(distinct(project_accession)) as n_projects from metadata inner join AMR_analysis using(run_accession) group by collection_date, host, country"
        sample_counts = get_data(
            query=q_count,
            db_args=args.db_args,
            fix_gene=False,
            oldest_year=args.oldest_year
        )

        if args.verbose:
            print(" * Pulled {} rows matching the query:\n\t{}".format(sample_counts.shape[0], q_count))
    

    if args.create_maps:
        countryLocator = CountryLookup(data_dir='../data/country_shapes/')

        df = df_main.copy()
        
        df.replace(countryLocator.altGeoNames, inplace=True)

        cm_country, cm_country_closed, cm_country_alr, cm_country_alr_total = create_cms(
            df = df,
            group_col='country',
        )
        
        map_plots = plot_maps(
            cm = cm_country_alr,
            cm_tot = cm_country_alr_total,
            ncols=5,
            nrows=6
        ) 
        map_plots.savefig(
            os.path.join(
                args.output, 'mcr_maps_cartopy.pdf'
            )
        )  
        map_plots.savefig(
            os.path.join(
                args.output, 'mcr_maps_cartopy.png'
            )
        )
