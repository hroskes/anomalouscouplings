from collections import Counter
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
    exp = Counter()
    for production in config.productionsforcombine:
        print production
        productionmodes = ("signal", "qqZZ", "ggZZ", "ZX")
        for flavor in "2e2mu", "4e", "4mu":
            print flavor
            print
            rates = getrates(flavor, "fordata", production, *sys.argv[1:])
            for p, rate in zip(productionmodes, rates.split()[1:]):
                exp[p,flavor] += float(rate)
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
    print "Table:"
    print
    flavors = ["4e", "4mu", "2e2mu"]
    fmt = "{:<10}"*4
    print fmt.format(*([""]+flavors))
    totals = Counter()
    for p in "qqZZ", "ZX", "ggZZ":
        print fmt.format(*([p] + ["{:.2f}".format(exp[p,f]) for f in flavors]))
        for f in flavors:
            totals[f] += exp[p,f]
    print fmt.format(*(["bkg"] + ["{:.2f}".format(totals[f]) for f in flavors]))
    print fmt.format(*(["signal"] + ["{:.2f}".format(exp["signal",f]) for f in flavors]))
    print fmt.format(*(["observed"] + [sum(tfiles[DataTree(production, flavor).treefile].candTree.GetEntries() for production in config.productionsforcombine) for flavor in flavors]))
    print
    print
    print "Total:"
    print "Expected:", totalexp
    print "Observed:", totalobs
