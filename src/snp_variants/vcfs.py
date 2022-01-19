import os
import sys
import subprocess
import pandas as pd
from io import StringIO
import glob
import subprocess
from Bio import SeqIO

from seqs import write_fasta, parse_fasta

def load_vcf(vcfFile):
    o = subprocess.check_output(f"bcftools view {vcfFile} | grep -v '^##'", shell=True).decode()
    vcfDF = pd.read_csv(StringIO(o), sep='\t')
    
    info_cols = ['DP', 'AD', 'AF', 'RAF', 'DEL', 'AD6']
    p_info = r"DP=(\d+);AD=(\d+);AF=(\d+\.?\d*);RAF=(\d+\.?\d*);DEL=(\d+);AD6=(\S+)"
    vcfDF[info_cols] = vcfDF['INFO'].str.extractall(p_info).values

    return vcfDF

def snp_filtering_old(vcfFile, allele, destDir, alleleRecord, AD=5, AF=.9, verbose=False):

    empty_return = (None, None)

    # find file
    if os.path.getsize(vcfFile) == 0:
        if verbose:
            print(f"Empty vcf file: {vcfFile}")
        return empty_return
    
    fasta_file = os.path.join(destDir, 'REFERENCE.fasta')
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(alleleRecord, output_handle, "fasta")

    # check if allele in vcf before doing anything
    wc = subprocess.check_output(f"bcftools view -t {allele} {vcfFile}  2> /dev/null | grep -v '^#' | wc -l", shell=True)
    if wc.decode().strip() == '0':
        if verbose: print(f"{allele} not in vcf")
        return empty_return
    
    # subset vcf file to only contain target allele
    vcf_prefix = os.path.join(destDir, os.path.basename(vcfFile.replace('.vcf.gz', '.' + allele)))
    vcf_cmd = f"bcftools view -t {allele} -Oz -o {vcf_prefix + '.vcf.gz'} {vcfFile} 2>/dev/null"
    subprocess.run(vcf_cmd, shell=True)

    # index vcf file
    subprocess.run(f"bcftools index {vcf_prefix + '.vcf.gz'}", shell=True)
    
    # filter snps
    output_vcf = vcf_prefix + '.filtered.vcf.gz'
    subprocess.run(f"bcftools filter -i '( TYPE=\"snp\" & AD >= {AD} & AF >= {AF} )' -Oz -o {output_vcf} {vcf_prefix + '.vcf.gz'}", shell=True)
    
    # return output_vcf
    vcfDF = load_vcf(output_vcf)
    if vcfDF.shape[0] > 0:
        # make new consensus 
        consensusFastaFile = vcf_prefix + '.consensus.fa'
        cmd=f"bcftools consensus -f {fasta_file} -o {consensusFastaFile} {output_vcf}"
        subprocess.run(cmd, shell=True)

        # load sequence
        seqRecord = list(SeqIO.parse(consensusFastaFile, 'fasta'))[0]

        return vcfDF, seqRecord
    else:
        return empty_return



    # norm vcf
    # norm_cmd = f"bcftools norm --check-ref s -t {allele} --fasta-ref {fasta_file} -Oz -o {norm_vcf} {vcfFile}"
    # p = subprocess.run(norm_cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    # if p.returncode != 0:
    #     if verbose:
    #         print("Failed to normalize vcf:", p.stderr.decode())
    #     return None
    
    # # filter vcf
    # filter_cmd = f"bcftools filter -i '(ALT=\".\" || (TYPE=\"snp\" && AD >= {AD} && AF >= {AF}))' -Oz -o {output_vcf} {norm_vcf}"
    # print(filter_cmd)
    # p = subprocess.run(filter_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if p.returncode != 0:
    #     print(f"Failed to filter {vcfFile}:", p.stderr.decode())
    #     return None

    # p = subprocess.run(f"bcftools index {output_vcf}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def snp_filtering(vcfFile, allele, targetDir, alleleSeq, AD=5, AF=.9, prefix=''):

    os.makedirs(targetDir, exist_ok=True)

    # filter vcf
    filterOut = vcfFile.replace(".vcf.gz", f".filtered.vcf.gz")
    filterOut = os.path.join(targetDir, os.path.basename(filterOut))

    cmd = f"bcftools filter -i '( TYPE=\"snp\" & AD >= {AD} & AF >= {AF} )' -t {allele} -Oz -o {filterOut} {vcfFile}"
    subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    # index vcf
    subprocess.run(f"bcftools index {filterOut}", shell=True, stderr=subprocess.PIPE)

    # make allele reference file
    refFile = os.path.join(targetDir, f'{allele}.fasta')
    write_fasta(record=alleleSeq, fastaFile=refFile)

    # make consensus
    consensusFile = filterOut.replace(".filtered.vcf.gz", ".consensus.fasta")
    cmd = f"bcftools consensus -f {refFile} -o {consensusFile} -p {prefix} {filterOut}"
    subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    # read consensus file into python and return it
    consensus = parse_fasta(consensusFile)[prefix + allele]
    return consensus