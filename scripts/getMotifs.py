import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq


class GetMotifs:
    def __init__(self, motif, file, path_count, fseqfw, fseqrv, fposfw, fposrv, fseqlen):
        self.mot = motif
        self.fseq = file
        self.fcount = path_count
        self.outseqfw = fseqfw
        self.outseqrv = fseqrv
        self.outposfw = fposfw
        self.outposrv = fposrv
        self.fseqlen = fseqlen

    def checkinput(self, mod):
        nucl_ok = ["A", "T", "C", "G", "N"]
        aa_ok = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "X"]
        if mod == "nucl":
            GetMotifs.completter(self, nucl_ok)
        if mod == "prot":
            GetMotifs.completter(self, aa_ok)

    def completter(self, lref):
        valid_mot = True
        for let in str(self.mot):
            if str(let).upper() not in lref and str(let).isalpha():
                valid_mot = False
                raise ValueError("Invalid residue {0} in motif {1}".format(let, self.mot))
        if valid_mot:
            sys.stdout.write("{} is a valid input motif\n".format(self.mot))

    def replaceN(self, mod):
        if mod == "nucl":
            mot4search = self.mot.replace("N", ".")
        else:
            mot4search = self.mot.replace("X", ".")
        return mot4search

    def searchmotif(self, mot2search, seq):
        hits = mot2search.finditer(str(seq))
        n = 0
        mot_spec = []
        for hsp in hits:
            n += 1
            mot_spec.append((hsp.start(), hsp.end(), hsp.group()))
        mot_spec.append(n)
        return mot_spec

    def writemotif(self, mod):
        mot2search = re.compile(GetMotifs.replaceN(self, mod))
        seqidx = 0
        for rec in SeqIO.parse(self.fseq, "fasta"):
            res = GetMotifs.searchmotif(self, mot2search, str(rec.seq))
            n_fwd, n_rev, c_seqlen = 0, 0, len(str(rec.seq))
            self.fseqlen.write("{0}\t{1}\n".format(rec.id, c_seqlen))
            seqidx += 1
            if len(res) > 1:
                n_fwd = res[-1]
                for i in res[:-1]:
                    pstart, pend, seqf = i[0]+1, i[1], i[2]
                    GetMotifs.writepos(self.outposfw, rec.id, self.mot, seqf, pstart, pend)
                    GetMotifs.writeseq(self.outseqfw, rec.id, self.mot, pstart, pend, seqf)

            if mod == "nucl":
                revseq = Seq(str(rec.seq)).reverse_complement()
                res_rv = GetMotifs.searchmotif(self, mot2search, str(revseq))
                if len(res_rv) > 1:
                    n_rev = res_rv[-1]
                    for i in res_rv[:-1]:
                        pstart, pend, seqf = c_seqlen-i[1]+1, c_seqlen-i[0], i[2]
                        GetMotifs.writepos(self.outposrv, rec.id, self.mot, seqf, pstart, pend)
                        GetMotifs.writeseq(self.outseqrv, rec.id, self.mot, pstart, pend, seqf)
            GetMotifs.writecount(self.fcount, rec.id, self.mot, str(n_fwd), str(n_rev), "seq"+str(seqidx))

    @staticmethod
    def writeseq(fseq, idseq, regexpr, pos_start, pos_end, seq2write):
        fseq.write(">{0}_{1}:{2} {3}\n{4}\n".format(idseq, pos_start, pos_end, regexpr, seq2write))

    @staticmethod
    def writepos(fpos, idseq, regexpr, found_comb, pos_start, pos_end):
        fpos.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(idseq, regexpr, found_comb, pos_start, pos_end))

    @staticmethod
    def writecount(fcount, idseq, regexpr, count_fwd, count_rev, seqidx):
        fcount.write("{5}\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(idseq, regexpr, count_fwd, count_rev, int(count_fwd)+int(count_rev), seqidx))
