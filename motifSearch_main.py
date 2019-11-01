# Developed by J. Keller & C. Libourel


from scripts import getMotifs as gm
from scripts import general_functions as gf
import os, sys, subprocess
import configparser
import argparse

if __name__ == '__main__':
    #config = configparser.ConfigParser()
    #config.read(sys.argv[1])
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help= "path to sequence file")
    parser.add_argument("-o", "--output", help= "path to output directory")
    parser.add_argument("-m", "--motif", nargs = "+", help = "motifs to search (semi-colon separated)")
    parser.add_argument("--mod", help = "nucleotide (nucl) or protein (prot) mode")
    args = parser.parse_args()

    gen_out = args.output
    seq = args.file
    mod = args.mod
    lmot = args.motif
    path2rscript = os.path.abspath(os.path.dirname(sys.argv[0])) + os.path.sep + "scripts" + os.path.sep + "motifSearch2_plot.R"

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

    gf.setrpath(path2rscript, os.path.abspath("res"), os.path.abspath("stats"))
    r_cmd = ["Rscript", path2rscript]
    subprocess.run(r_cmd)
    sys.stdout.write("Job done!\nResults are located at: {}\n".format(path_out))
