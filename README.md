# **motifSearch**

This is the README file of motifSearch

This script aims to retrieve motifs in nucleotide or protein sequences.
motifSearch is written in Python3.6.8 and should work with all Python3.6 distributions and higher. A R script is called for the graph draw.
Requirements:  
- Python3.6  
- R v3.6.0  
- Basic Python packages (included in basic Python installation): *sys, os, re, subprocess*  
- Additional packages that are not included in Python basic installation: *pandas, biopython*  
To install any additional Python packages, the easiest way is to use the pip command:  
`pip install package_name`  
- R packages: *data.table, ggplot2, tidyr, dplyr, gridExtra*  
To install R packages, use the following command in R:  
`install.packages("package_name", dependencies=TRUE, INSTALL_opts = c('--no-lock'))`  

motifSearch compatibility:  
- Linux (developped and tested in Ubuntu 18.04 LTS)  
- Windows (tested using the PowerShell 6: https://github.com/PowerShell/PowerShell)  

## Install ##  
```
git clone https://github.com/jeankeller/motifSearch.git  
cd motifSearch/  
python setup.py install
```

## Usage ##  
`python3 path/to/motifSearch_main.py -f path/to/sequence_file -o path/to/output_directory --mod {nucl or prot} -m motifs_to_search` (space separated, between quotes)  
Typing `python3 motifSearch_main.py -h` will bring up help and usage information.

Motif to be searched have to be provided using the syntax of regular expression (see details below).
Input motif format:
- ".": any character (e.g T.C means search for T + any letter + C)  
- "[]": one of the character contains between brackets (e.g T[AG]C means search for T + A or G + C); [^A] means any character except A  
- "*": any or many occurences of the previous character (e.g AT* means search for A + 0 or many T)  
- "{}": precise number of occurences for the previous character (e.g AT{2,5}C means search for A + 2 to 5 T + C; whereas AT{2,}C means search for A + at least 2 T + C)  

Useful link to use regular expression: https://docs.python.org/3/library/re.html#regular-expression-syntax  

## Output ##

motifSearch creates an output folder:  
- motifSearch_results_YYYYMMDD-HHMMSS  
  - data  
    - sequences.length => input sequences length (used for graphical output)  
  - res
    - motifs_count.txt => table containing motifs count (in forward, reverse and total) per sequence per motif  
    - motifs_positions_fwd.txt => table containing motifs positions for each motif identified in forward (exact motif indicated)  
    - motifs_positions_rev.txt => table containing motifs positions for each motif identified in reverse (exact motif indicated)  
    - motifs_sequences_fwd.fa => fasta file containing exact motif sequences in forward  
    - motifs_sequences_rev.fa => fasta file containing exact motif sequences in reverse  
  - stats
    - MOTIF_motifs_count_per_sequence.pdf => barplot representing number of motif per sequence  
    - MOTIF_motifs_density_plots_windows_XXbp.pdf => density plot representing the number of motifs in the promoter region (0 correspond to the start of the input fasta file sequence)  
    - Hist_MOTIFmotif.pdf => histogramms representing the number of sequences per number of motifs  (forward and reverse, only forward and only reverse)  
  - desc_results.txt => short description of results  
  - motifs_count_nonull.txt => correspond to /res/motifs_count.txt without sequences with any motifs  

## **Authors** ##

- Jean Keller (kellerjeanphd@gmail.com)
- Cyril Libourel (cyril.libourel@gmail.com)
