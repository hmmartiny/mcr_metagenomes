import argparse
import os
import re
import sys
import glob
import subprocess
import pandas as pd
from io import StringIO
import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from vcfs import snp_filtering, load_vcf
from coverages import load_coverage, filter_coverage, get_coverage_stats
from seqs import aligner, make_anno, snipit, parse_fasta, write_multi_fasta, sequence_cleaner, fasttree, snp_dists, viz_trees, summarize_occurences, summarize_seqCounts

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
        '-MC' ,'--minimum-coverage',
        type=float,
        default=98,
        help='Minimum percentage coverage of template',
        dest='minimum_coverage'
    )

    covArgs.add_argument(
        '-MCD','--minimum-depth-of-coverage',
        type=float,
        default=5,
        help='Minimum depth of coverage',
        dest='minimum_depth_of_cov'
    )

    covArgs.add_argument(
        '-MPv','--minimum-p-value',
        type=float,
        default=0.05,
        help='Minimum p-value',
        dest='min_pvalue'
    )

    covArgs.add_argument(
        '-MQI','--minimum-query-identity',
        type=float,
        default=90,
        help='Minimum query identity',
        dest='min_qidentity'
    )

    snpArgs = parser.add_argument_group('SNP filters')
    snpArgs.add_argument(
        '-AD' ,'--allele-depth',
        type=float,
        default=5,
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
    if args.verbose:
        s = "Input arguments:"
        s += f"\n * Minimum template coverage {args.minimum_coverage}"
        s += f"\n * Minimum average depth: {args.minimum_depth_of_cov}"
        s += f"\n * Minimum query identity: {args.min_qidentity}"
        s += f"\n * Max p-value: {args.min_pvalue}"
        s += f"\n * SNP allele depth: {args.AD}"
        s += f"\n * SNP allele frequency: {args.AF}"
        s += f"\n"

        print(s)

    # make subset directory
    destDir = os.path.join(args.fasta_path, f"chosen_cov{args.minimum_coverage}_depth{args.minimum_depth_of_cov}_qident{args.min_qidentity}_pvalue{args.min_pvalue}_AD{args.AD}_AF{args.AF}")
    os.makedirs(destDir, exist_ok=True)
    if args.verbose:
        print("Files written to:", os.path.realpath(destDir), "\n")

    # load coverage statistics
    covDF = load_coverage(dataPath=args.data_path)

    # make summary table of coverages
    if args.verbose:
        covGeneSummary = covDF.groupby(['gene']).apply(get_coverage_stats)
        covAlleleSummary = covDF.groupby(['Allele', 'RefSeq']).apply(get_coverage_stats)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None): 
            print("Coverage statistics of gene groups:")
            print(covGeneSummary)
            print()

            print("Coverage statistics of templates:")
            print(covAlleleSummary)
            print()
        sys.exit()
    
    # filter
    filtCovDF = filter_coverage(
        covDF, 
        min_cov=args.minimum_coverage, 
        min_depth=args.minimum_depth_of_cov,
        min_pval=args.min_pvalue,
        min_qident=args.min_qidentity
    )


    # parse reference sequences
    references = parse_fasta(os.path.join(args.fasta_path, 'ResFinder.fasta'))


    vcfFiles = os.path.join(args.fasta_path, 'all_vcfs', f'*.vcf.gz')

    occurenceDFs = []
    allUniqueConsensus = {}
    allConsensus = {}
    allSeqCounts = {}

    noVCFs = []

    alleles = sorted(filtCovDF.dropna(subset=['run_accession'])['Template'].unique())
    pbar = tqdm.tqdm(alleles, disable=not args.verbose)
    
    # for allele, alleledata in filtCovDF.dropna(subset=['run_accession']).groupby('Template'):
    for allele in pbar:

        pbar.set_description(allele)
        alleledata = filtCovDF.loc[
            (filtCovDF['Template'] == allele) & 
            (~filtCovDF['run_accession'].isna())
        ]

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

        seqs = {s: seq for s, seq in seqs.items() if not s[-3:-2] == '.v'}
        uniqueSeqs, occurenceDF, seqCounts = sequence_cleaner(seqs, allele)
        allUniqueConsensus[allele] = uniqueSeqs
        allConsensus[allele] = seqs
        allSeqCounts[allele] = seqCounts
        occurenceDFs.append(occurenceDF)

        occFile = os.path.join(destSubDir, 'new_variant_overview.csv')
        summarize_seqCounts(seqCounts, flatten=False).to_csv(occFile)
        # occurenceDF.to_csv(occFile)
        # see if we found other than just the reference
        
        if len(uniqueSeqs) > 1:
            multiConsensusFile = os.path.join(destSubDir, 'all_consensus_found.fasta')
            write_multi_fasta(seqs, fastaFile=multiConsensusFile)

            # align consensus sequences
            aln_file=aligner(
                fasta_file=multiConsensusFile
            )
            
            snipitFile = snipit(aln_file, verbose=args.verbose)

            # only align unique
            uniqueConsensusFile = os.path.join(destSubDir, 'all_unique_consensus_found.fasta')
            write_multi_fasta(uniqueSeqs, fastaFile=uniqueConsensusFile)
            aln_file=aligner(fasta_file=uniqueConsensusFile)
            snpdistFile=snp_dists(aln_file)
            treeFile=fasttree(aln_file)
            snipitFile = snipit(aln_file, verbose=args.verbose)

            
            viz_trees(
                treefile=treeFile, 
                snpdistfile=snpdistFile,
                occfile=occFile,
                dest=destSubDir)
   
    
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
    
    occFile = os.path.join(destDir, 'new_variant_overview.csv')
    sumOcc = summarize_seqCounts(allSeqCounts)
    sumOcc.to_csv(occFile)

    if args.verbose:
        print(sumOcc.sum(1))
    
    
    # flatten
    allConsensus = {vk: vv for v in allConsensus.values() for vk, vv in v.items()}
    allUniqueConsensus = {vk: vv for v in allUniqueConsensus.values() for vk, vv in v.items()}
    
    # write combined consensus fasta file
    allConsensusFile = os.path.join(destDir, 'all_consensus_found.fasta')
    write_multi_fasta(allConsensus, fastaFile=allConsensusFile)

    allUniqueConsensusFile = os.path.join(destDir, 'all_unique_consensus_found.fasta')
    write_multi_fasta(allUniqueConsensus, fastaFile=allUniqueConsensusFile)


    # align consensus sequences
    aln_file=aligner(
        fasta_file=allUniqueConsensusFile
    )

    # make tree
    treeFile=fasttree(aln_file)

    # calculate snp distances
    snpdistFile=snp_dists(aln_file)

    # viz trees
    viz_trees(
        treefile=treeFile,
        snpdistfile=snpdistFile,
        occfile=occFile,
        dest=destDir)