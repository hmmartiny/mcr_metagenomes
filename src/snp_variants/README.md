## Detecting SNP variants 
Takes .res files from KMA and a reference fasta files to construct consensus sequences with filters: Minimum template coverage (`MINIMU_COVERAGE`), Minimum depth of coverage (`MINIMUM_DEPTH_OF_COV`), Minimum p-value (`MAX_PVALUE`), Minimum query identity (`MIN_QIDENTITY`), and for SNPs minium allele depth (`AD`) and minimum allele frequency (`AF`).

```
> python snp_variants.py -h
usage: snp_variants.py [-h] [-v] -dp DATA_PATH -fp FASTA_PATH -op OUTPUT_PATH [-MC MINIMUM_COVERAGE] [-MCD MINIMUM_DEPTH_OF_COV] [-MPv MIN_PVALUE] [-MQI MIN_QIDENTITY] [-AD AD] [-AF AF]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Be verbose (default: False)

Directory arguments:
  -dp DATA_PATH, --data-path DATA_PATH
                        Data dir (default: None)
  -fp FASTA_PATH, --fasta-path FASTA_PATH
                        Fasta dir (default: None)
  -op OUTPUT_PATH, --output-path OUTPUT_PATH
                        Output directory (default: None)

Coverage filters:
  -MC MINIMUM_COVERAGE, --minimum-coverage MINIMUM_COVERAGE
                        Minimum percentage coverage of template (default: 98)
  -MCD MINIMUM_DEPTH_OF_COV, --minimum-depth-of-coverage MINIMUM_DEPTH_OF_COV
                        Minimum depth of coverage (default: 5)
  -MPv MAX_PVALUE, --maximum-p-value MIN_PVALUE
                        Maximum p-value (default: 0.05)
  -MQI MIN_QIDENTITY, --minimum-query-identity MIN_QIDENTITY
                        Minimum query identity (default: 90)

SNP filters:
  -AD AD, --allele-depth AD
                        Minimum depth of SNP allele. (default: 5)
  -AF AF, --allele-frequency AF
                        Minimum allele frequency of ALT (default: 0.9)
```