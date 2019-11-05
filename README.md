# **motifSearch**

This is the README file of motifSearch

This script aims to retrieve motifs in nucleotide or protein sequences.
motifSearch is written in Python3.6.8 and should work with all Python3.x distributions. A R script is called for the graph draw.
Requirements:  
- Python3.x  
- R v3.6.0  
- Basic Python packages (included in basic Python installation): sys, os, re  
- Additional packages that are not included in Python basic installation: subprocess, configparser, biopython  
- R packages: data.table, ggplot2, tidyr, dplyr, gridExtra  

To install any package additional Python packages, the easiest way is to use the pip command: `pip install package_name`  
To install R packages, use the followaing command in R:  
`install.packages("package_name", dependencies=TRUE, INSTALL_opts = c('--no-lock'))`

motifSearch usage:  
`python path/to/motifSearch_main.py -f path/to/sequence_file -o path/to/output_directory --mod {nucl or prot} -m motifs_to_search` (space separated)

Motif to be searched have to be provided using the syntax of regular expression (see details below).
Input motif format:
- ".": any character (e.g T.C means search for T + any letter + C)  
- "[]": one of the character contains between brackets (e.g T[AG]C means search for T + A or G + C); [^A] means any character except A  
- "*": any or many occurences of the previous character (e.g AT* means search for A + 0 or many T)  
- "{}": precise number of occurences for the previous character (e.g AT{2,5}C means search for A + 2 to 5 T + C; whereas AT{2,}C means search for A + at least 2 T + C)  

Useful link to use regular expression: LINK2ADD

