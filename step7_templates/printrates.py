from helperstuff import config
from helperstuff.combinehelpers import getrates
from helperstuff.enums import Analysis
import os
from projections import Projections
import ROOT
import sys

if __name__ == "__main__":
    for production in config.productionsforcombine:
        print production
        for flavor in "2e2mu", "4e", "4mu":
            print flavor
            print
            rates = getrates(flavor, "fordata", production, *sys.argv[1:], format="{:.2f} {:.2f} {:.2f} {:.2f}")
            print rates
            print "Total bkg: {:.2f}".format(sum(float(_) for _ in rates.split()[2:]))
            print "Total expected: {:.2f}".format(sum(float(_) for _ in rates.split()[1:]))
            otherargv = [_ for _ in sys.argv[1:] if _ != "fordata"]
            f = ROOT.TFile(os.path.join(Projections("fa3", flavor, "noenrich", production).saveasdir, Analysis("fa3").purediscriminant()+".root"))
            print "Observed:", f.c1.GetListOfPrimitives()[1].GetHists()[-1].Integral()
            print
        print
