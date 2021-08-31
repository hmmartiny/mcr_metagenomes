# Large-scale metagenomic analysis of mobilized colistin resistance genes
Code and notebooks to rerun analysis and regenerate figures in paper.

## Overview
`notebooks/` contains:
* one Jupyter Notebook named [`Plots.ipynb`](notebooks/Plots.ipynb) that contains code for creating figures 1-3, S1-S2 and runs PCA on CLR values.
* one Rmd notebook named [`mcr_plots.Rmd`](notebooks/mcr_plots.Rmd) that contains the code for running [ALDEx2](https://github.com/ggloor/ALDEx_bioc) analysis and creating figures 4-5, S3-S4.

The notebooks imports functions from files in the [`src/`](src/) folder:
* [`src/funcs.py`](src/funcs.py): functions for retrieving data from database and compositional functions.
* [`src/dataviz/dataviz.py`](src/dataviz/dataviz.py):
* [`src/dataviz/maps.py`](src/dataviz/maps.py):