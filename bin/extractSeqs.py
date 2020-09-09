#!/usr/bin/env python

import sys
import os.path

def readFitnesses(directory, generation, batches, nsubsamples, label):
    """Read all the *.fitness.csv files for the indicated `generation'
from `directory', and return the top `nsubsamples' ones."""
    allfitnesses = []

    for b in range(batches):
        fitfile = "{}/{}.{}.{}.indeces.fitness.csv".format(directory, label, generation, b+1)
        subfile = "{}/{}.{}.{}.indeces.subsamples.csv".format(directory, label, generation, b+1)
        print(fitfile)
        if os.path.isfile(fitfile) and os.path.isfile(subfile):
            with open(fitfile, "r") as f, open(subfile, "r") as s:
                ids = f.readline().rstrip("\r\n").split(",")
                f.readline()
                f.readline()
                ftns = f.readline().rstrip("\r\n").split(",")
                for i in range(len(ids)):
                    samples = s.readline().rstrip("\r\n").split(",")
                    allfitnesses.append([ ftns[i], samples ])

    allfitnesses.sort(key=lambda e:e[0], reverse=True)
    # print(allfitnesses)
    
    return [ f[1] for f in allfitnesses[:nsubsamples] ]

def readIDs(idsfile):
    ids = []
    with open(idsfile, "r") as f:
        for line in f:
            ids.append(line.strip())
    return ids

def getId(line):
    sp = line.find(" ")
    if sp > 0:
        return line[1:sp]
    else:
        return line[1:]

def extractSeqsOld(fastafile, idsfile):
    ids = readIDs(idsfile)
    good = False
    with open(fastafile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line) == 0:
                continue
            if line[0] == ">":
                seqid = getId(line)
                if seqid in ids:
                    sys.stdout.write(line + "\n")
                    good = True
                else:
                    good = False
            else:
                if good:
                    sys.stdout.write(line + "\n")

def extractSeqs(fastafile, outfile, indeces):
    n = 0
    good = False
    indeces = [ int(i) for i in indeces ]
    with open(outfile, "w") as out, open(fastafile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line) == 0:
                continue
            if line[0] == ">":
                n += 1
                if n in indeces:
                    out.write(line + "\n")
                    good = True
                else:
                    good = False
            else:
                if good:
                    out.write(line + "\n")

def main(args):
    outdir = args[0]
    fasta = args[1]
    generation = int(args[2]) - 1
    batches = int(args[3])
    nsubsamples = int(args[4])
    label = args[5]

    indeces = readFitnesses(outdir, generation, batches, nsubsamples, label)
    outfiles = []
    print(indeces)
    print(nsubsamples)
    for i in range(nsubsamples):
        outfile = "{}/subsample.{}.{}.fa".format(outdir, label, i+1)
        outfiles.append(outfile)
        extractSeqs(fasta, outfile, indeces[i])
    sys.stdout.write(" ".join(outfiles) + "\n")

if __name__ == "__main__":
    args = sys.argv[1:]

    if len(args) == 6:
        main(args)
    else:
        sys.stdout.write("Usage: extractSeqs.py outdir fastafile generation nbatches nsubsamples\n")
        sys.exit(1)
