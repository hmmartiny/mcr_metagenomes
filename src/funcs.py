#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import json
import pycoda

def closure(row, k=100):
    return k * row / row.sum()

def alr(row):
    return np.log(row[:-1] / row[-1])

def tot_alr(row):
    return np.log(sum(row[:-1]) / row[-1])

def count_matrix(df, groupby_col, index_pivot, column_pivot, values='fragmentCountAln'):
    
    agg_df = df.groupby(groupby_col).agg({values: 'sum'}).reset_index()

    pivoted = agg_df.pivot(index=index_pivot, columns=column_pivot, values=values)
    pivoted.fillna(0, inplace=True)
    
    return pivoted

def count_samples(gene_df, bac_df, gene='mcr', col='collection_year', val_col='runs', log=np.log10):
    """Count samples with and without gene hits"""
    
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

def pca(data, parts, metadf=None, left_on=None, left_index=False):
    
    # replace zeroes
    df = data[parts].coda.zero_replacement()
    
    # close
    df = df.coda.closure(100)
      
    # clr transform
    df_clr = df.coda.center().coda.clr()
    
    # run SVD
    scores, eig_val, loadings = np.linalg.svd(df_clr)
    n_eig = len(eig_val)
    
    # explained variance by each pc
    explained_variance = ( eig_val ** 2 / np.sum(eig_val ** 2)) * 100
    explained_variance = pd.DataFrame(
        explained_variance,
        index=[f'PC{i+1}' for i in range(n_eig)],
        columns=['Explained Variance']
    )
    #explained_variance.columns.name = 'Variance'
    explained_variance.index.name = 'Principal Component'

    # make scores df
    scores = pd.DataFrame(
        scores.T,
        columns=df_clr.index
    )

    # loadings df
    # take inner product of eig_val * identity matrix of number of eigenvalues
    loadings = pd.DataFrame(
        np.inner(
            eig_val * np.identity(n_eig),
            loadings.T[0:n_eig, 0:n_eig]
        ),
        columns=df_clr.columns[0:n_eig],
        index=[f'PC{i+1}' for i in range(n_eig)]
    )

    scales = [np.max(np.abs(loadings.values)),
                  [np.max(np.abs(scores.loc[idx].values)) for idx in [1, 2]]]


    scores = scores[0:2]
    scores = (scales[0] * (scores.T / scales[1])).T
    
    if meta_df is not None:
        scores_ext = scores.T.merge(metadf, right_index=True, left_index=left_index, left_on=left_on)
        return explained_variance, scores, loadings, scales, scores_ext
    
    return explained_variance, scores, loadings, scales
    
