from helperstuff import config
from helperstuff.combinehelpers import getrates
from helperstuff.enums import Analysis
import os
from projections import Projections
import ROOT
import sys

if __name__ == "__main__":
    for flavor in "2e2mu", "4e", "4mu":
        print flavor
        print
        print getrates(flavor, "fordata", config.productionforcombine, *sys.argv[1:])
        if "fordata" in sys.argv[1:]:
            otherargv = [_ for _ in sys.argv[1:] if _ != "fordata"]
            f = ROOT.TFile(os.path.join(Projections("fa3", flavor, "noenrich", config.productionforcombine).saveasdir, Analysis("fa3").purediscriminant()+".root"))
            print "Observed:", f.c1.GetListOfPrimitives()[1].GetHists()[-1].Integral()
        print
