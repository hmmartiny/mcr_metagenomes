import pandas as pd
import os
import re

def split_name(x):
    
    p = re.compile(r"((E|D|S)RR[0-9]{6,})")
    runid, template = x.split(':')
    
    m = p.findall(runid)
    runid = m[0][0]
    return pd.Series([runid, template])
    
def split_template(x):
    p = re.compile(r'^(mcr-\d(_|\.)\d+)((_1)|(_NZ))?_(\S+)$')
    m = p.findall(x)[0]
    return m[0], m[-1]

def load_coverage(dataPath, reload=False):

    columns = [
        'Template',
        'Score',
        'Expected',
        'Template_length',
        'Template_Identity',
        'Template_Coverage',
        'Query_Identity',
        'Query_Coverage',
        'Depth',
        'q_value',
        'p_value'
    ]

    mcrFile = os.path.join(dataPath, 'mcr_res_fixed.tsv')


    if reload:
        df = pd.read_csv(os.path.join(dataPath, 'mcr_res.tsv'), sep='\t')
        df.columns = columns

        mcrDF = df.loc[df.Template.str.contains(':mcr')].copy()
        mcrDF[['runid', 'Template']] = mcrDF.apply(lambda x: split_name(x['Template']), axis=1, result_type='expand')
        mcrDF['gene'] = mcrDF.Template.str.extract(r'(mcr-\d)')

        mcrDF[['Allele', 'RefSeq']] = mcrDF.apply(lambda x: split_template(x['Template']), axis=1, result_type='expand')
        
        metadata = pd.read_csv(os.path.join(dataPath, 'mcr_df.csv'), infer_datetime_format=True)[['run_accession', 'host', 'country', 'collection_date']]
        metadata['collection_date'] = pd.to_datetime(metadata['collection_date'], errors='coerce')
        metadata['year'] = metadata['collection_date'].dt.year
        mcrDF = mcrDF.merge(metadata, left_on='runid', right_on='run_accession', how='left')


        mcrDF.to_csv(mcrFile)
    else:
        mcrDF = pd.read_csv(mcrFile, index_col=0)

    return mcrDF


def filter_coverage(covDF, min_cov, min_depth, min_pval, min_qident, cov_col='Template_Coverage', depth_col='Depth', pval_col='p_value', qident_col='Query_Identity'):

    filtDF = covDF.loc[
        (covDF[cov_col] >= min_cov) & 
        (covDF[depth_col] >= min_depth) & 
        (covDF[pval_col] <= min_pval ) &
        (covDF[qident_col] >= min_qident)
    ]

    return filtDF

def get_coverage_stats(row):
    
    runCount = row['runid'].nunique()
    
    # coverage
    meanCoverage = row['Template_Coverage'].mean()
    maxCoverage = row['Template_Coverage'].max()
    medianCoverage = row['Template_Coverage'].median()
    minCoverage = row['Template_Coverage'].min()
    
    # depths
    meanDepth = row['Depth'].mean()
    maxDepth = row['Depth'].max()
    medianDepth = row['Depth'].median()
    minDepth = row['Depth'].min()
    
    res = pd.DataFrame(
        [
            runCount, 
            meanCoverage, 
            minCoverage, 
            medianCoverage, 
            maxCoverage, 
            meanDepth, 
            minDepth,
            medianDepth,
            maxDepth
        ]
    ).T

    res.rename(
        columns={
            0: 'No. of metagenomes',
            1: 'Mean template coverage',
            2: 'Minimum template coverage',
            3: 'Median template coverage',
            4: 'Maximum template cpverage',
            5: 'Mean depth of coverage',
            6: 'Minimum depth of coverage',
            7: 'Median depth of coverage',
            8: 'Maximum depth of coverage'
        }, 
        inplace=True
    )

    return res
