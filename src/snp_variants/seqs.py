import glob
import os
from Bio import SeqIO
from Bio import Align
from Bio.Align.Applications import MafftCommandline
import subprocess
import pandas as pd

def parse_fasta(fastaFile):
    record_dict = SeqIO.to_dict(
        SeqIO.parse(fastaFile, "fasta")
    )
    return record_dict

def write_multi_fasta(recordDict, fastaFile):

    with open(fastaFile, "w") as output_handle:
        SeqIO.write(recordDict.values(), output_handle, "fasta")

def write_fasta(record, fastaFile):
    with open(fastaFile, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

def load_allele_fasta(fastaDir, allele, subset_list=[]):

    fastaFile = glob.glob(os.path.join(fastaDir, '*' + allele + '.fa*'))
    assert len(fastaFile) == 1, f'There can only be one fasta file matching the allele {allele} (found {len(fastaFile)})'
    # record_dict = SeqIO.to_dict(SeqIO.parse(fastaFile[0], "fasta"))
    record_dict = parse_fasta(fastaFile[0])

    # update keys
    new_dict = {allele : record_dict[allele]}
    for k in record_dict.keys():
        nk = k
        if '|' in nk and (nk.endswith('_1') or nk.endswith('_2')):
            nk = nk[:-2]

        new_dict[nk] = record_dict[k]
    
    if len(subset_list) > 0:
        return {h: new_dict[h] for h in [allele] + subset_list}

    return new_dict
    
def aligner(fasta_file):
    # run mafft
    mafft_cline = MafftCommandline(input=fasta_file)
    
    # extract alignment
    stdout, _ = mafft_cline()
    aln_file = fasta_file.replace('.fa', '.aln.fa')
    with open(aln_file, "w") as handle:
        handle.write(stdout)
        
    return aln_file

def get_snp_sites(aln_file):

    output_file = aln_file.replace('.fa', '.snps.vcf')
    cmd = f"snp-sites -v -o {output_file} {aln_file}"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if not p.returncode == 0:
        print("Failed to create snp vcf file")
        print(p.stderr.decode())
        return None
    
    # return matrix
    snpDF = pd.read_csv(output_file, sep='\t', skiprows=3)

    snpDF.rename(
        columns={
            c: c[:-2] for c in snpDF.columns if c.startswith('mcr') and c.endswith('_1')
        },
        inplace=True
    )
    return snpDF


def snipit(aln_file, label_file=None, labels='header,runid', verbose=False):
    
    output_file = aln_file.replace('.fasta', '.plot')
    cmd=f"snipit {aln_file} -o {output_file}"
    if label_file is not None:
        cmd += " -l {label_file} --l-header '{labels}'"
    
    p = subprocess.run(
        cmd, 
        shell=True, 
        stderr=subprocess.PIPE, 
        stdout=subprocess.PIPE)
    
    output_file =  output_file + '.png'
    if p.returncode != 0:
        if verbose:
            print(f"Failed to run snipit for {aln_file}")
            print(p.stderr.decode())
        output_file = None
    
    return output_file
    
def make_anno(df, chosenSeqs, output_file, allele):
    
    
    # make dataframe
    meta = pd.DataFrame([], columns=['header', 'runid', 'year', 'country', 'host'])
    for i, (k,v) in enumerate(chosenSeqs.items()):
        
        meta.loc[i, 'header'] = v.id
        meta.loc[i, 'runid'] = k
        
        if k != 'ref':
            try:
                m = df.loc[df['run_accession'] == k, ['year', 'country', 'host']].copy().drop_duplicates().values[0]
            except IndexError:
                continue
            meta.loc[i, 'year'] = m[0]
            meta.loc[i, 'country'] = m[1]
            meta.loc[i, 'host'] = m[2]
    
    
    meta.fillna('', inplace=True)
    meta['label'] = meta.apply(lambda x: f"{x['runid']}|{x['year']}|{x['country']}|{x['host']}", axis=1)
    meta.loc[0, 'label'] = f'Reference {allele}'
            
    
    # save to csv
    meta.to_csv(output_file, index=None)
    
    return meta

def sequence_cleaner(seqRecords, allele):

    uniqueRecords = {}
    uniqueSeqs2name = {}
    occurenceDF = pd.DataFrame(columns=list(seqRecords.keys()))

    variantCount = 1
    for recid, record in seqRecords.items():
        sequence = record.seq

        if sequence not in uniqueSeqs2name.keys():
            if recid != allele:
                seqName = f"{allele}.v{variantCount}"
                variantCount += 1
            else:
                seqName = allele
            
            uniqueSeqs2name[sequence] = seqName
        
            record.id = seqName
            record.name = seqName
            record.description = seqName
            uniqueRecords[seqName] = record
        
        seqName = uniqueSeqs2name[sequence]
        occurenceDF.loc[seqName, recid] = 1
        
    occurenceDF.fillna(0, inplace=True)
    occurenceDF.drop(columns=[allele], inplace=True)
    
    return uniqueRecords,  occurenceDF
        
def fasttree(aln_file):
    outputFile = aln_file.replace('.fasta', '.tree')

    cmd = f"FastTree -q -gtr -nt {aln_file} > {outputFile}"
    subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

def snp_dists(aln_file):
    outputFile = aln_file.replace('.fasta', '.snp_dists.tsv')
    cmd=f"snp-dists -q {aln_file} > {outputFile}"
    subprocess.run(cmd, shell=True)

def viz_trees(path, dest, width=10, height=7):

    cmd = f"./plot_trees.R -width {width} -height {height} {path} {dest}"
    subprocess.run(cmd, shell=True)