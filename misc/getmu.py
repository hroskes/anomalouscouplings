#!/usr/bin/env python

import glob
import sys
filename = glob.glob("cards_{}_March27_submitjobs/higgsCombine_obs_lumi35.8671.MultiDimFit.mH125.root".format(sys.argv[1]))
assert len(filename) == 1
filename = filename[0]

import ROOT

f = ROOT.TFile(filename)
t = f.limit

def absfai(i):
    t.GetEntry(i)
    return abs(t.CMS_zz4l_fai1)

entry = min((i for i, _ in enumerate(t)), key=absfai)

t.GetEntry(entry)

try:
    for a in "CMS_zz4l_fai1", "deltaNLL", "muV_scaled", "muf_scaled":
        print "{:>10}: {:.2f}".format(a, getattr(t, a))
except:
    t.Show()
    raise
