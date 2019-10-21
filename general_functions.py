import time
import os


def getcurrenttime():
    return time.strftime("%Y%m%d-%H%M%S")


def create_gen_out(gen_out):
    os.makedirs(gen_out + os.path.sep + "motifSearch_results_{}".format(getcurrenttime()))
    return gen_out + os.path.sep + "motifSearch_results_{}".format(getcurrenttime())

def setrpath(rscript, path_in, path_out):
    with open(rscript, "r") as fin:
        lines = fin.readlines()
    fin.close()
    lines[5] = "PATH = " + '"' +str(path_in + "/").replace(os.path.sep, "/") + '"' + "\n"
    lines[6] = "OUTPATH = " + '"' + str(path_out + "/").replace(os.path.sep, "/") + '"' + "\n"
    fout = open(rscript, "w")
    fout.write("".join(lines))
    fout.close()


