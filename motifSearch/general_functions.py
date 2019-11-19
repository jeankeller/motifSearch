import time, os, sys
import pandas as pds


def getcurrenttime():
    return time.strftime("%Y%m%d-%H%M%S")


def create_gen_out(gen_out):
    os.makedirs(gen_out + os.path.sep + "motifSearch_results_{}".format(getcurrenttime()))
    return gen_out + os.path.sep + "motifSearch_results_{}".format(getcurrenttime())


def setrpath(rscript, path_in, path_out):
    with open(rscript, "r") as fin:
        lines = fin.readlines()
        print(lines)
    fin.close()
    lines[5] = "PATH = " + '"' +str(path_in + "/").replace(os.path.sep, "/") + '"' + "\n"
    lines[6] = "OUTPATH = " + '"' + str(path_out + "/").replace(os.path.sep, "/") + '"' + "\n"
    fout = open(rscript, "w")
    fout.write("".join(lines))
    fout.close()


def getstats(fcount, gen_out, lmot):
    sys.stdout.write("Generating summary files...    ")
    fc = os.path.join(gen_out, "res") + os.path.sep + "motifs_count.txt"
    df = pds.read_csv(fc, sep = "\t")
    fstats = open(gen_out + os.path.sep + "desc_results.txt", "w")
    fnonull = open(gen_out + os.path.sep + "motifs_count_nonull.txt", "w")
    dfuniq = df.groupby("sequence_index")["motif"].nunique()
    nseq, nmot = len(dfuniq), len(lmot)
    fstats.write("Number of sequences analyzed: {0}\nNumber of motifs searched: {1}\n".format(nseq, nmot))
    dfnonulltot = df.loc[df["total"] != 0]
    dfnonulltot.to_csv(fnonull, sep = "\t", index = False)
    for mot in lmot:
        fstats.write("\n#####################\n\nMotif: {0}\n".format(mot))
        dfmot = df.loc[df["motif"] == mot]
        df_fwd, df_rev, df_tot = dfmot.loc[dfmot["count_fwd"] != 0], dfmot.loc[dfmot["count_rev"] != 0], dfmot.loc[dfmot["total"] != 0]
        nseqfwd, nseqrev, nseqtot = len(df_fwd), len(df_rev), len(df_tot)
        perc_seqfwd, perc_seqrev, perc_seqtot = round(nseqfwd/nseq*100, 2), round(nseqrev/nseq*100, 2), round(nseqtot/nseq*100, 2)
        nmotfwd, nmotrev, nmottot = df_fwd["count_fwd"].sum(), df_rev["count_rev"].sum(), df_tot["total"].sum()
        min_fwd, max_fwd = dfmot["count_fwd"].min(), dfmot["count_fwd"].max()
        min_rev, max_rev = dfmot["count_rev"].min(), dfmot["count_rev"].max()
        min_tot, max_tot = dfmot["total"].min(), dfmot["total"].max()
        fstats.write("Forward:\n\tNumber of sequences with at least 1 motif: {0} ({1}% of sequences)\n\tNumber of motifs found: {2}\n\tmin: {3}\n\tmax: {4}\n".format(nseqfwd, perc_seqfwd, nmotfwd, min_fwd, max_fwd))
        fstats.write("Reverse:\n\tNumber of sequences with at least 1 motif: {0} ({1}% of sequences)\n\tNumber of motifs found: {2}\n\tmin: {3}\n\tmax: {4}\n".format(nseqrev, perc_seqrev, nmotrev, min_rev, max_rev))
        fstats.write("Total:\n\tNumber of sequences with at least 1 motif: {0} ({1}% of sequences)\n\tNumber of motifs found: {2}\n\tmin: {3}\n\tmax: {4}\n".format(nseqtot, perc_seqtot, nmottot, min_tot, max_tot))
    sys.stdout.write("done!\n\n################################\n\n")




