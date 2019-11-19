#!/usr/bin/env python
#-*- coding: utf-8 -*-
#
#Copyright 2019 Jean Keller & Cyril Libourel
#This program (motifSearch) is distributed under the terms of the GNU General Public License v3+

from motifSearch import getMotifs as gm
from motifSearch import general_functions as gf
import os, sys, subprocess
import argparse
import pkg_resources

def main():
    parser = argparse.ArgumentParser(usage="motifsearch [options] ....",
                                     description="""A Python package to search motifs in nucleotide/
                                     protein sequences""",
                                     epilog="""Written by Jean Keller (kellerjeanphd@gmail.com) and
                                     Cyril Libourel (cyril.libourel@gmail.com), LRSV University of
                                     Toulouse III/CNRS. (c) 2019- Released under the terms of the
                                     GNU General Public License v3 or later (LGPLv3+)""")
    #Mandatory arguments
    parser.add_argument("-f", "--file", required = True, type=str, help= "path to sequence file [mandatory]")
    parser.add_argument("-m", "--motif", required = True, type=str, nargs = "+", help = "motifs to search, space separated [mandatory]")
    parser.add_argument("--mod", type=str, required = True, help = "nucleotide (nucl) or protein (prot) mode [mandatory]")

    #Optional arguments
    parser.add_argument("-o", "--output", default="./", type=str, help= """path to output directory,
                        by default, use the current directory [optional]""")
    
    args = parser.parse_args()

    gen_out = args.output
    seq = args.file
    mod = args.mod
    lmot = args.motif
    path2rscript = pkg_resources.resource_filename("motifSearch", "motifSearch2_plot.R")
    
    path_out = gf.create_gen_out(gen_out)
    os.chdir(path_out)
    os.makedirs("res")
    os.makedirs("data")
    os.makedirs("stats")

    fcount = open("res" + os.path.sep + "motifs_count.txt", "w")
    fpos_fw = open("res" + os.path.sep + "motifs_positions_fwd.txt", "w")
    fpos_rv = open("res" + os.path.sep + "motifs_positions_rev.txt", "w")
    fseq_fw = open("res" + os.path.sep + "motifs_sequences_fwd.fa", "w")
    fseq_rv = open("res" + os.path.sep + "motifs_sequences_rev.fa", "w")
    flenseq = open("data" + os.path.sep + "sequences.length", "w")

    fcount.write("sequence_index\tsequence_id\tmotif\tcount_fwd\tcount_rev\ttotal\n")
    fpos_fw.write("sequence_id\tmotif\tsequence\tstart\tstop\n")
    fpos_rv.write("sequence_id\tmotif\tsequence\tstart\tstop\n")

    sys.stdout.write("\n\n################################\n\n")
    for mot in lmot:
        sys.stdout.write("Starting search for {}\n".format(mot))
        c_mot = gm.GetMotifs(mot, seq, fcount, fseq_fw, fseq_rv, fpos_fw, fpos_rv, flenseq)
        gm.GetMotifs.checkinput(c_mot, mod)
        gm.GetMotifs.writemotif(c_mot, mod)
        sys.stdout.write("done!")
        sys.stdout.write("\n\n################################\n\n")


    fcount.close()
    fpos_fw.close()
    fpos_rv.close()
    fseq_fw.close()
    fpos_rv.close()
    flenseq.close()
    gf.getstats(fcount, path_out, lmot)

    sys.stdout.write("Generating the graphs, this may take a while!\n")

    r_cmd = ["Rscript", path2rscript]
    subprocess.run(r_cmd)
    sys.stdout.write("\n\n################################\n\nJob done!\nResults are located at: {}\n".format(path_out))
    
