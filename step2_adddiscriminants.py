from array import array
from helperstuff.enums import flavors, hypotheses
from helperstuff.samples import Sample
from helperstuff.treewrapper import TreeWrapper
import os
import ROOT
import sys

def adddiscriminants(*args):

    sample = Sample(*args)
    reweightingsamples = sample.reweightingsamples()

    filename = sample.CJLSTfile()
    f = ROOT.TFile.Open(filename)
    Counters = f.Get("{}/Counters".format(sample.TDirectoryname()))
    Counters_reweighted = f.Get("{}/Counters_reweighted".format(sample.TDirectoryname()))
    if not Counters_reweighted:
        Counters_reweighted = None
    if not Counters:
        raise ValueError("No Counters in file "+filename)

    newfilename = sample.withdiscriminantsfile()
    print newfilename
    if os.path.exists(newfilename):
        return

    t = ROOT.TChain("{}/candTree".format(sample.TDirectoryname()))
    t.Add(filename)

    newf = ROOT.TFile.Open(newfilename, "recreate")
    newt = t.CloneTree(0)

    treewrapper = TreeWrapper(t, sample, Counters=Counters, Counters_reweighted=Counters_reweighted)

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
    for hypothesis in hypotheses:
        adddiscriminants("ggH", hypothesis)
    for flavor in flavors:
        adddiscriminants("ggZZ", flavor)
    adddiscriminants("qqZZ")
    adddiscriminants("ZX")
    adddiscriminants("data")
