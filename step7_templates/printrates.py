from helperstuff import config
from helperstuff.combinehelpers import getrates
from helperstuff.enums import DataTree
from helperstuff.filemanager import tfiles
import os
from projections import Projections
import ROOT
import sys

if __name__ == "__main__":
    totalexp = totalobs = 0
    for production in config.productionsforcombine:
        print production
        for flavor in "2e2mu", "4e", "4mu":
            print flavor
            print
            rates = getrates(flavor, "fordata", production, *sys.argv[1:])
            print "{:.2f} {:.2f} {:.2f} {:.2f}".format(*(float(_) for _ in rates.split()[1:]))
            print "Total bkg: {:.2f}".format(sum(float(_) for _ in rates.split()[2:]))
            print "Total expected: {:.2f}".format(sum(float(_) for _ in rates.split()[1:]))
            totalexp += sum(float(_) for _ in rates.split()[1:])
            print "Observed:", tfiles[DataTree(production, flavor).treefile].candTree.GetEntries()
            totalobs += tfiles[DataTree(production, flavor).treefile].candTree.GetEntries()
            print
        print

    print
    print
    print "Total:"
    print "Expected:", totalexp
    print "Observed:", totalobs
