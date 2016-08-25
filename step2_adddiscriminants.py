from array import array
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import flavors, decayonlyhypotheses, prodonlyhypotheses, productions
from helperstuff.samples import Sample
from helperstuff.treewrapper import TreeWrapper
import os
import ROOT
import sys

definitelyexists = Sample("ggH", "0+", "160225")
print definitelyexists.CJLSTfile()
assert xrd.exists(definitelyexists.CJLSTfile())

def adddiscriminants(*args):

    sample = Sample(*args)
    reweightingsamples = sample.reweightingsamples()

    filename = sample.CJLSTfile()
    newfilename = sample.withdiscriminantsfile()
    print newfilename
    if os.path.exists(newfilename):
        return

    isdummy = False
    if not xrd.exists(filename):
        print filename
        isdummy = True
        #give it a tree so that it can get the format, but not fill any entries
        filename = definitelyexists.CJLSTfile()

    f = ROOT.TFile.Open(filename)
    Counters = f.Get("{}/Counters".format(sample.TDirectoryname()))
    Counters_reweighted = f.Get("{}/Counters_reweighted".format(sample.TDirectoryname()))
    if not Counters_reweighted:
        Counters_reweighted = None
    if not Counters:
        raise ValueError("No Counters in file "+filename)

    t = ROOT.TChain("{}/candTree".format(sample.TDirectoryname()))
    t.Add(filename)

    treewrapper = TreeWrapper(t, sample, Counters=Counters, Counters_reweighted=Counters_reweighted, isdummy=isdummy)

    if os.path.exists(newfilename):
        return

    newf = ROOT.TFile.Open(newfilename, "recreate")
    newt = t.CloneTree(0)

    discriminants = {}
    for discriminant in treewrapper.toaddtotree:
        discriminants[discriminant] = array('d', [0])
        newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")

    for entry in treewrapper:
        for discriminant in discriminants:
            discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
        newt.Fill()

    newf.Write()
    newf.Close()

if __name__ == '__main__':
    for production in productions:
        for hypothesis in decayonlyhypotheses:
            adddiscriminants("ggH", hypothesis, production)
        if config.analysistype == "prod+dec":
            for hypothesis in prodonlyhypotheses:
                adddiscriminants("VBF", hypothesis, production)
        for flavor in flavors:
            adddiscriminants("ggZZ", flavor, production)
        adddiscriminants("qqZZ", production)
        adddiscriminants("ZX", production)
        if config.analysistype == "ICHEP16":
            adddiscriminants("VBF", "0+", production)
        adddiscriminants("ZH", "0+", production)
        adddiscriminants("WplusH", "0+", production)
        adddiscriminants("WminusH", "0+", production)
        adddiscriminants("ttH", "0+", production)
        adddiscriminants("data", production, "unblind")
        adddiscriminants("data", production, "blind")
