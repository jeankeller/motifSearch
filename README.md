This is the README file of motifSearch

This script aims to retrieve motifs in nucleotide or protein sequences.
motifSearch is written in Python3.6.8 and should work with all Python3.x distributions. A R script is called for the graph draw.
Requirements:
    - Python3.x
    - R vXXX
    - Basic Python packages (included in basic Python installation): sys, os, re
    - Additional packages that are not included in Python basic installation: subprocess, configparser, biopython
    - R packages: ggplot2, tidyr, dplyr, gridExtra

To install any package additional Python packages, the easiest way is to use the pip command: pip install package_name
To install R packages, use the followaing command: DESCRIPTION OF INSTALLATION (CYRIL)
Input parameters of motifSearch are given in the configfile.
Motif to be searched haeve to be provided using the syntax of regular expression (see details below).

