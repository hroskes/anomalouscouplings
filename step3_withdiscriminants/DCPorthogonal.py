#!/usr/bin/env python

import argparse, itertools

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("VBForVH", choices="VBF VH".split())
  p.add_argument("ggHorttH", choices="ggH ttH".split())
  p.add_argument("disc", choices="full D_0minus D_0hplus D_L1 D_L1Zg D_int".split())
  args = p.parse_args()

import ROOT, uncertainties

from helperstuff.discriminants import discriminant
from helperstuff.samples import Sample
from helperstuff.utilities import TFile

def DCPorthogonal(VBForVH, ggHorttH, disc):
   category, DCP = {
     "VBF": (2, "D_CP_VBF"),
     "VH": (4, "D_CP_HadVH"),
   }[VBForVH]
   if disc == "full":
     disc = discriminant("D_4couplings_" + {"VBF": "VBF", "VH": "HadVH"}[VBForVH] + "decay")
   else:
     disc = discriminant(disc + "_" + {"VBF": "VBF", "VH": "HadVH"}[VBForVH] + ("decay_3bins" if disc != "D_int" else "_2bins"))
   discname = disc.name
   axisrange = "{0.bins},{0.min},{0.max}".format(disc)

   args = {
     "ggH": ("ggH", "MC@NLO"),
     "ttH": ("ttH",),
   }[ggHorttH]
   with TFile(Sample("fCP0.5", "0+", "190703_2017", *args).withdiscriminantsfile()) as f:
     c = ROOT.TCanvas()
     t = f.candTree
     t.Draw(DCP+":"+discname+">>"+VBForVH+"("+axisrange+",2,-1,1)", "category_0P_or_0M_or_a2_or_L1_or_L1Zg == {} && 105 < ZZMass && ZZMass < 140".format(category))
     h = getattr(ROOT, VBForVH)

     allpos = []
     allneg = []

     for x in xrange(1, h.GetNbinsX()+1):
       neg = uncertainties.ufloat(h.GetBinContent(x, 1), h.GetBinError(x, 1))
       pos = uncertainties.ufloat(h.GetBinContent(x, 2), h.GetBinError(x, 2))
       allpos.append(pos)
       allneg.append(neg)

     totalpos = sum(allpos)
     totalneg = sum(allneg)       
     print (totalpos - totalneg) / (totalpos + totalneg)

     for pos, neg in itertools.izip(allpos, allneg):
       if not pos+neg: continue
       print "{:10} {:10}".format(((pos-neg) / (pos+neg)), ((pos-neg) / (pos+neg)) / ((totalpos - totalneg) / (totalpos + totalneg)) - 1)

if __name__ == "__main__":
  DCPorthogonal(**args.__dict__)
