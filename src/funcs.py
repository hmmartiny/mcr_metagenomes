#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import json
import pycoda

def closure(row, k=100):
    """Closure a composition to k

    Parameters
    ----------
    row : pd.Series
        Pandas DataFrame row
    k : int, optional
        Value to close to, by default 100

    Returns
    -------
    pd.Series
        Closed compositional row
    """
    return k * row / row.sum()

def alr(row):
    """Calculate ALR values for a row by using the last component as the denominator

    Parameters
    ----------
    row : pd.Series
        Pandas DataFrame row

    Returns
    -------
    pd.Series
        ALR transformed row
    """
    return np.log(row[:-1] / row[-1])

def tot_alr(row):
    """Sum all row contents except for the last component and calculate ALR"""
    return np.log(sum(row[:-1]) / row[-1])

def count_matrix(df, groupby_col, index_pivot, column_pivot, values='fragmentCountAln'):
    """Create count matrix by pivoting a long dataframe into wide format

    Parameters
    ----------
    df : pd.DataFrame
        A compositional DataFrame
    groupby_col : str or list
        Column(s) to group data by
    index_pivot : str or list
        Column(s) to use as index in pivot operation
    column_pivot : str or list
        Column(s) to use as column in pivot operation
    values : str or list, optional
        Column(s) to use as values in pivot operation, by default 'fragmentCountAln'

    Returns
    -------
    pd.DataFrame
        Compositional dataframe transformed to be of correct format
    """
    
    agg_df = df.groupby(groupby_col).agg({values: 'sum'}).reset_index()

    pivoted = agg_df.pivot(index=index_pivot, columns=column_pivot, values=values)
    pivoted.fillna(0, inplace=True)
    
    return pivoted

def count_samples(gene_df, bac_df, gene='mcr', col='collection_year', val_col='runs', log=np.log10):
    """Count samples with and without gene hits

    Parameters
    ----------
    gene_df : pd.DataFrame
        Compositional dataframe with gene counts
    bac_df : pd.DataFrame
        DataFrame with bacterial gene counts
    
    
    
    """
    
    timesamples = pd.DataFrame(
        gene_df.drop_duplicates('run_accession')[col].value_counts()
    )
    
    timesamples = timesamples.merge(
        bac_df.groupby(col).agg({val_col: 'sum'}),
        right_index=True, left_index=True
    )
    
    timesamples[f'Samples\nwithout ${gene}$'] = log(
        timesamples[val_col] - timesamples[col]
    )
    timesamples[f'Samples\nwith ${gene}$'] = log(timesamples[col])
    
    return timesamples

def load_json(handle):
    """Load a json file"""
    with open(handle, 'r') as js:
        data = json.load(js)
    return data

def query_db(query, db_args='db_args.json'):
    """Pull data with a query from a MySQL db"""
    if isinstance(db_args, str):
        db_args = load_json(db_args)
    
    engine = create_engine(
        'mysql+pymysql://{user}:{pw}@{host}:{port}/{db}'.format(
            user=db_args['user'],
            pw=db_args['passwd'],
            host=db_args['host'],
            port=db_args['port'],
            db=db_args['db']
            )
        )    
    data = pd.read_sql_query(query, engine)
    return data
