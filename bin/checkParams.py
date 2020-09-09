#!/usr/bin/env python

import sys

DISTOPTS = ["mean", "median", "max"]

def main(fn, fev, fel, n, opt):
    if fn + fev + fel != 1.0:
        sys.stderr.write("Error: the sum of the three fractions is not equal to 1.\n")
        sys.exit(2)
    s = round(fn * n) + round(fev * n) + round(fel * n)
    if s != n:
        sys.stderr.write("Error: rounding changes totals, adjust fractions.\n")
        sys.exit(3)
    if opt not in DISTOPTS:
        sys.stderr.write("Error: dist.opts should be one of {}.\n".format(", ".join(DISTOPTS)))
        sys.exit(4)

if __name__ == "__main__":
    args = sys.argv[1:]
    main(float(args[0]),
         float(args[1]),
         float(args[2]),
         int(args[3]),
         args[4])
