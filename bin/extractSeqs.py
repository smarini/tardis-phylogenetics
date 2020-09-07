#!/usr/bin/env python

import sys

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

def extractSeqs(fastafile, idsfile):
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

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 2:
        extractSeqs(args[0], args[1])
    else:
        sys.stdout.write("Usage: extractSeqs.py fastafile idsfile\n")
        sys.exit(1)
