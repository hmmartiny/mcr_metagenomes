import argparse
import os
import re
import sys
import glob
import subprocess
import pandas as pd
from io import StringIO

import matplotlib.pyplot as plt
import seaborn as sns

from vcfs import snp_filtering, load_vcf
from coverages import load_coverage, filter_coverage
from seqs import aligner, make_anno, snipit, parse_fasta, write_multi_fasta, sequence_cleaner, fasttree, snp_dists, viz_trees

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Be verbose',
        dest='verbose'
    )

    dirArgs = parser.add_argument_group('Directory arguments')
    dirArgs.add_argument(
        '-dp', '--data-path',
        type=str,
        required=True,
        help='Data dir',
        dest='data_path'
    )

    dirArgs.add_argument(
        '-fp', '--fasta-path',
        type=str,
        required=True,
        help='Fasta dir',
        dest='fasta_path'
    )

    dirArgs.add_argument(
        '-op', '--output-path',
        type=str,
        required=True,
        help='Output directory',
        dest='output_path'
    )

    covArgs = parser.add_argument_group('Coverage filters')
    covArgs.add_argument(
        '--minimum-coverage',
        type=float,
        default=98,
        help='Minimum percentage coverage of template',
        dest='minimum_coverage'
    )

    covArgs.add_argument(
        '--minimum-depth-of-coverage',
        type=float,
        default=.2,
        help='Minimum depth of coverage',
        dest='minimum_depth_of_cov'
    )

    covArgs.add_argument(
        '--minimum-p-value',
        type=float,
        default=0.05,
        help='Minimum p-value',
        dest='min_pvalue'
    )

    covArgs.add_argument(
        '--minimum-query-identity',
        type=float,
        default=50,
        help='Minimum query identity',
        dest='min_qidentity'

    )

    snpArgs = parser.add_argument_group('SNP filters')
    snpArgs.add_argument(
        '-AD' ,'--allele-depth',
        type=float,
        default=10,
        help='Minimum depth of SNP allele.',
        dest='AD'
    )

    snpArgs.add_argument(
        '-AF', '--allele-frequency',
        type=float,
        default=0.9,
        help='Minimum allele frequency of ALT',
        dest='AF'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # load coverage statistics
    covDF = load_coverage(dataPath=args.data_path)
    
    # filter
    filtCovDF = filter_coverage(
        covDF, 
        min_cov=args.minimum_coverage, 
        min_depth=args.minimum_depth_of_cov,
        min_pval=args.min_pvalue,
        min_qident=args.min_qidentity
    )

    if args.verbose:
        print(f"Went from {covDF.shape[0]} to {filtCovDF.shape[0]} after applying coverage filters.")

    # parse reference sequences
    references = parse_fasta(os.path.join(args.fasta_path, 'ResFinder.fasta'))


    # stats
    if args.verbose:
        print("Number of KMA produced consensus sequences found per gene that passed filters:")
        print(filtCovDF['gene'].value_counts())
        print()

    # make subset directory
    destDir = os.path.join(args.fasta_path, f"chosen_cov{args.minimum_coverage}_depth{args.minimum_depth_of_cov}_AD{args.AD}_AF{args.AF}")
    os.makedirs(destDir, exist_ok=True)

    vcfFiles = os.path.join(args.fasta_path, 'all_vcfs', f'*.vcf.gz')

    occurenceDFs = []
    allConsensus = {}

    noVCFs = []
    
    for allele, alleledata in filtCovDF.dropna(subset=['run_accession']).groupby('Template'):
        if args.verbose:
            print(allele, alleledata.shape[0])
        
        destSubDir = os.path.join(destDir, allele)
        os.makedirs(destSubDir, exist_ok=True)
        
        seqs = {allele: references[allele]}

        # find vcfs
        for runid in alleledata['run_accession']:
            vcfFile = glob.glob(os.path.join(args.fasta_path, 'all_vcfs', f'ResFinder_20200125__{runid}*.vcf.gz'))
            if len(vcfFile) != 1 :
                # if args.verbose:
                    # print(f'Found either none or multiple matches for {runid}: {len(vcfFile)} .. skipping')
                noVCFs.append(runid)
                continue
            vcfFile = vcfFile[0]
            consensus = snp_filtering(
                vcfFile, 
                allele=allele, 
                targetDir=os.path.join(destSubDir, 'vcf_consensus'), 
                alleleSeq=references[allele],
                prefix=runid + ':'
            )

            seqs[runid] = consensus

        uniqueSeqs, occurenceDF = sequence_cleaner(seqs, allele)
        
        allConsensus[allele] = uniqueSeqs
        occurenceDFs.append(occurenceDF)
        
        # see if we found other than just the reference

        if len(uniqueSeqs) > 1:
            multiConsensusFile = os.path.join(destSubDir, allele + '.consensus.fasta')
            write_multi_fasta(seqs, fastaFile=multiConsensusFile)

            # align consensus sequences
            aln_file=aligner(
                fasta_file=multiConsensusFile
            )
            
            snipitFile = snipit(aln_file, verbose=args.verbose)

    if args.verbose:
        print("Number of samples where a vcf file is not found:", len(set(noVCFs)))

    with open(os.path.join(destDir, 'missing_vcfs.list'), 'w') as handle:
        print("\n".join(list(set(noVCFs))), file=handle)

    occurenceDFs = pd.concat(occurenceDFs).fillna(0)
    occurenceDFs[
        ['minCov', 'minCovDepth', 'minPval', 'minQIdent', 'minAD', 'minAF']
    ] = [
        args.minimum_coverage,
        args.minimum_depth_of_cov,
        args.min_pvalue,
        args.min_qidentity,
        args.AD,
        args.AF
        ]  
    
    occurenceDFs.to_csv(os.path.join(destDir, 'new_variant_overview.csv'))

    # flatten
    allConsensus = {vk: vv for v in allConsensus.values() for vk, vv in v.items()}
    
    # write combined consensus fasta file
    allConsensusFile = os.path.join(destDir, 'all_consensus_found.fasta')
    write_multi_fasta(allConsensus, fastaFile=allConsensusFile)

    # align consensus sequences
    aln_file=aligner(
        fasta_file=allConsensusFile
    )

    # make tree
    fasttree(aln_file)

    # calculate snp distances
    snp_dists(aln_file)

    # viz trees
    viz_trees(path=destDir, dest=destDir)