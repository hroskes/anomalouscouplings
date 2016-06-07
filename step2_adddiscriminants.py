from array import array
from helperstuff.enums import hypotheses
from helperstuff.samples import Sample
from helperstuff.treewrapper import TreeWrapper
import os
import ROOT
import sys

def adddiscriminants(*args):

    sample = Sample(*args)
    reweightingsamples = sample.reweightingsamples()

    filename = sample.CJLSTfile()
    f = ROOT.TFile(filename)
    Counters_reweighted = f.Get("ZZTree/Counters_reweighted")

    newfilename = sample.withdiscriminantsfile()
    print newfilename
    if os.path.exists(newfilename):
        return

    t = ROOT.TChain("ZZTree/candTree")
    t.Add(filename)

    newf = ROOT.TFile.Open(newfilename, "recreate")
    newt = t.CloneTree(0)

    discriminants = {}
    for discriminant in TreeWrapper.toaddtotree:
        discriminants[discriminant] = array('d', [0])
        newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")

    treewrapper = TreeWrapper(t, sample, Counters_reweighted)
    for entry in treewrapper:
        for discriminant in discriminants:
            discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
        newt.Fill()

    newf.Write()
    newf.Close()

if __name__ == '__main__':
    for hypothesis in hypotheses:
        adddiscriminants("ggH", hypothesis)


