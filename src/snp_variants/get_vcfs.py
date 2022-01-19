import pandas as pd
from coverages import load_coverage, filter_coverage
import subprocess
import pexpect
import tempfile
import tqdm 

data_path = 'data'
targetDir = 'data/snp_data/all_vcfs/'

# load coverage statistics
covDF = load_coverage(dataPath=data_path)


# filter
filtCovDF = filter_coverage(
    covDF, 
    min_cov=98, 
    min_depth=0.2,
    min_pval=0.05,
    min_qident=50
)

def run_scp(cmd, timeout=30):

    ''
    fname = tempfile.mktemp()                                                                                                                                                  
    fout = open(fname, 'w') 

    child = pexpect.spawn(cmd, timeout=timeout)
    child.expect(['[pP]assword: '])  
    child.sendline('mRm74zEr')                                                                                                                                                   
    child.logfile = fout                                                                                                                                                       
    child.expect(pexpect.EOF)                                                                                                                                                  
    child.close()            
    fout.close()   

    fin = open(fname, 'r')                                                                                                                                                     
    stdout = fin.read()                                                                                                                                                        
    fin.close()
    print(stdout)

run_ids = filtCovDF['run_accession'].unique().tolist()
c2paths = []
for runid in tqdm.tqdm(run_ids):
    # get projid
    o = subprocess.check_output(
        f"mysql -N -e \"select project_accession from run_results where run_accession='{runid}'\"",
        shell=True
    ).decode()
    

    c2path=f"/home/projects/cge/scratch/ena/{o.strip()}/kma/ResFinder_20200125/ResFinder_20200125__{runid}*.vcf.gz"
    c2paths.append(c2path)


    # cmd = f"scp computerome2:{c2path} {targetDir}"
    # print(cmd)
    # run_scp(cmd)

print(c2paths)
print(" ".join(c2paths))

cmd = f"scp -T computerome2:\"{' '.join(c2paths)}\" {targetDir}"
print(cmd)
subprocess.run(cmd, shell=True)