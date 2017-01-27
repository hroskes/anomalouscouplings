#!/usr/bin/env python
from array import array
from collections import OrderedDict
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import flavors, hffhypotheses, ProductionMode, productions
from helperstuff.samples import Sample
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeepWhileOpenFile
import os
import ROOT
import sys

definitelyexists = Sample("VBF", "0+", config.productionsforcombine[0])
if not xrd.exists(definitelyexists.CJLSTfile()):
    raise ValueError("{} does not exist!".format(definitelyexists.CJLSTfile()))

def adddiscriminants(*args):

  sample = Sample(*args)
  reweightingsamples = sample.reweightingsamples()

  filename = sample.CJLSTfile()
  newfilename = sample.withdiscriminantsfile()
  print newfilename
  with KeepWhileOpenFile(newfilename+".tmp") as kwof:
    if not kwof:
        return
    if os.path.exists(newfilename):
        return

    isdummy = False
    if not xrd.exists(filename):
        isdummy = True
    else:
        f = ROOT.TFile.Open(filename)
        if not f.Get("{}/candTree".format(sample.TDirectoryname())):
            isdummy = True

    if isdummy:
        print "{} does not exist or is bad, using {}".format(filename, definitelyexists.CJLSTfile())
        #give it a tree so that it can get the format, but not fill any entries
        filename = definitelyexists.CJLSTfile()
        f = ROOT.TFile.Open(filename)

    Counters = f.Get("{}/Counters".format(sample.TDirectoryname()))
    if not Counters:
        raise ValueError("No Counters in file "+filename)

    t = ROOT.TChain("{}/candTree".format(sample.TDirectoryname()))
    t.Add(filename)

    if f.Get("{}/candTree_failed".format(sample.TDirectoryname())):
        t_failed = ROOT.TChain("{}/candTree_failed".format(sample.TDirectoryname()))
        t_failed.Add(filename)
    else:
        t_failed = None

    treewrapper = TreeWrapper(t, sample, Counters=Counters, failedtree=t_failed, isdummy=isdummy)

    if os.path.exists(newfilename):
        return

    failed = False
    try:
        newf = ROOT.TFile.Open(newfilename, "recreate")
        newt = t.CloneTree(0)
        if treewrapper.effectiveentriestree is not None:
            treewrapper.effectiveentriestree.SetDirectory(newf)

        discriminants = OrderedDict()
        for discriminant in treewrapper.toaddtotree:
            discriminants[discriminant] = array('d', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")
        for discriminant in treewrapper.toaddtotree_int:
            discriminants[discriminant] = array('i', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/I")

        try:
            for entry in treewrapper:
                for discriminant in discriminants:
                    discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
                newt.Fill()
        except:
            treewrapper.Show()
            raise
    except:
        failed = True
        raise
    finally:
        try:
            newf.Write()
            newf.Close()
        except:
            failed = True
        if failed:
            try:
                os.remove(newfilename)
            except:
                pass

if __name__ == '__main__':
    for production in productions:
        for productionmode in "ggH", "VBF", "ZH", "WH":
            for hypothesis in ProductionMode(productionmode).generatedhypotheses:
                adddiscriminants(productionmode, hypothesis, production)
        for hypothesis in hffhypotheses:
            adddiscriminants("HJJ", hypothesis, "0+", production)
            adddiscriminants("ttH", hypothesis, "0+", production)
        for flavor in flavors:
            adddiscriminants("ggZZ", flavor, production)
            if not flavor.hastaus:
                adddiscriminants("VBF bkg", flavor, production)
        for productionmode in "VBF", "ZH", "WplusH", "WminusH":
            adddiscriminants(productionmode, "0+", "POWHEG", production)
        adddiscriminants("ttH", "Hff0+", "0+", "POWHEG", production)
        adddiscriminants("qqZZ", production)
        adddiscriminants("ZX", production)
        adddiscriminants("data", production)
