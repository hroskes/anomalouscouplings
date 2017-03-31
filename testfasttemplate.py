#!/usr/bin/env python

"""
Run this script to test the FastHisto3D classes.
It should print the same number 6 times.
"""

import itertools

import ROOT

ROOT.gSystem.Load("libRooFit")
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

D1 = ROOT.RooRealVar("D1", "D1", 0, 1)
D2 = ROOT.RooRealVar("D2", "D2", 0, 1)
D3 = ROOT.RooRealVar("D3", "D3", 0, 1)
D1.setBins(50)
D2.setBins(50)
D3.setBins(50)

t = ROOT.TH3F("h", "h", 50, 0, 1, 50, 0, 1, 50, 0, 1)
for i, j, k in itertools.product(range(50), range(50), range(50)):
    t.Fill(i/50., j/50., k/50., 2500*i + 50*j + k)  #every bin has a different content

fastdatahist_f = ROOT.FastHisto3D_f(t)
fasthistfunc_f = ROOT.FastHisto3DFunc_f("fasthistfunc_f", "", ROOT.RooArgList(D1,D2,D3), fastdatahist_f)

fastdatahist_d = ROOT.FastHisto3D_d(t)
fasthistfunc_d = ROOT.FastHisto3DFunc_d("fasthistfunc_d", "", ROOT.RooArgList(D1,D2,D3), fastdatahist_d)

slowdatahist = ROOT.RooDataHist("slowdatahist", "", ROOT.RooArgList(D1,D2,D3), t)
slowhistfunc = ROOT.RooHistFunc("slowhistfunc", "", ROOT.RooArgSet(D1,D2,D3), slowdatahist)


def indicesfromxyz(h, x, y, z):
    binglobal = h.Fill(x, y, z, 0)
    from array import array
    x = array('i', [0])
    y = array('i', [0])
    z = array('i', [0])
    h.GetBinXYZ(binglobal, x, y, z)
    return x[0], y[0], z[0]

def test(d1val, d2val, d3val):
    D1.setVal(d1val)
    D2.setVal(d2val)
    D3.setVal(d3val)

    dvals = [D1.getVal(), D2.getVal(), D3.getVal()]
    assert dvals == [d1val, d2val, d3val]

    print indicesfromxyz(t, *dvals)

    allthesame = [
                  fastdatahist_f.GetAt(*dvals),
                  fasthistfunc_f.getVal(),
                  fastdatahist_d.GetAt(*dvals),
                  fasthistfunc_d.getVal(),
                  slowhistfunc.getVal(),
                  t.GetBinContent(*indicesfromxyz(t, *dvals)),
                 ]

    for _ in allthesame:
       print _
    assert len(set(allthesame)) == 1

test(.24, .36, .84) #all on bin boundary
test(.25, .37, .85) #none on bin boundary
test(.1, .4, .85)   #mix
test(0, .1, .4)     #middle of an edge, on a boundary
test(0, .1, .45)    #middle of an edge, not on a boundary
test(0, .1, .999)   #in a 2D corner (I forget the word for that)
test(0, 0, .999)    #in a 3D corner
test(.5, .5, 999)     #outside the range
