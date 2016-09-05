from array import array
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import flavors, decayonlyhypotheses, prodonlyhypotheses, productions
from helperstuff.samples import Sample
from helperstuff.treewrapper import TreeWrapper
import os
import ROOT
import sys

definitelyexists = Sample("ggH", "0+", "160901")
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
    couplings = f.Get("{}/couplings".format(sample.TDirectoryname()))
    if not Counters:
        raise ValueError("No Counters in file "+filename)
    if not Counters_reweighted:
        Counters_reweighted = None
    if not couplings:
        couplings = None

    t = ROOT.TChain("{}/candTree".format(sample.TDirectoryname()))
    t.Add(filename)

    treewrapper = TreeWrapper(t, sample, Counters=Counters, Counters_reweighted=Counters_reweighted, couplings=couplings, isdummy=isdummy)

    if os.path.exists(newfilename):
        return

    newf = ROOT.TFile.Open(newfilename, "recreate")
    newt = t.CloneTree(0)

    discriminants = {}
    for discriminant in treewrapper.toaddtotree:
        discriminants[discriminant] = array('d', [0])
        newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")
    for discriminant in treewrapper.toaddtotree_int:
        discriminants[discriminant] = array('i', [0])
        newt.Branch(discriminant, discriminants[discriminant], discriminant + "/I")

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
        for hypothesis in prodonlyhypotheses:
            adddiscriminants("VBF", hypothesis, production)
        for flavor in flavors:
            adddiscriminants("ggZZ", flavor, production)
            if not flavor.hastaus:
                adddiscriminants("VBF bkg", flavor, production)
        adddiscriminants("qqZZ", production)
        adddiscriminants("ZX", production)
        adddiscriminants("ZH", "0+", production)
        adddiscriminants("WplusH", "0+", production)
        adddiscriminants("WminusH", "0+", production)
        adddiscriminants("ttH", "0+", production)
        adddiscriminants("data", production, "unblind")
        adddiscriminants("data", production, "blind")
