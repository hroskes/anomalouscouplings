#!/usr/bin/env python

import os

import ROOT

from helperstuff import config
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import Sample
from helperstuff.utilities import mkdir_p

c = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

for disc in "D_0minus_HadVHdecay", "D_CP_HadVH":
    h = plotfromtree(
        reweightfrom=Sample("VBF", "0+", production),
        disc=disc,
        normalizeto1=True,
        category="VHHadrtagged",
        categorization="0P_or_0M",
        bins=50,
    )
    h.SetMinimum(0)
    h.Draw("hist")

    outdir = os.path.join(config.plotsbasedir, "offdiagonalLHEstudy", "VBFinVH")
    mkdir_p(outdir)
    for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(outdir, "{}.{}".format(disc, ext)))

for disc in "D_0minus_VBFdecay", "D_CP_VBF":
    h = plotfromtree(
        reweightfrom=Sample("VBF", "0+", production),
        disc=disc,
        normalizeto1=True,
        category="VBFtagged",
        categorization="0P_or_0M",
        bins=40,
    )
    h.SetMinimum(0)
    h.Draw("hist")

    outdir = os.path.join(config.plotsbasedir, "offdiagonalLHEstudy", "VBFinVBF")
    mkdir_p(outdir)
    for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(outdir, "{}.{}".format(disc, ext)))

for p in "ZH", "WH":
    for disc in "D_0minus_VBFdecay", "D_CP_VBF":
        h = plotfromtree(
            reweightfrom=Sample(p, "0+", production),
            disc=disc,
            normalizeto1=True,
            category="VBFtagged",
            categorization="0P_or_0M",
            bins=50,
        )
        h.SetMinimum(0)
        h.Draw("hist")

        outdir = os.path.join(config.plotsbasedir, "offdiagonalLHEstudy", "{}inVBF".format(p))
        mkdir_p(outdir)
        for ext in "png eps root pdf".split():
            c.SaveAs(os.path.join(outdir, "{}.{}".format(disc, ext)))

for p in "ZH", "WH":
    for disc in "D_0minus_HadVHdecay", "D_CP_HadVH":
        h = plotfromtree(
            reweightfrom=Sample(p, "0+", production),
            disc=disc,
            normalizeto1=True,
            category="VHHadrtagged",
            categorization="0P_or_0M",
            bins=50,
        )
        h.SetMinimum(0)
        h.Draw("hist")

        outdir = os.path.join(config.plotsbasedir, "offdiagonalLHEstudy", "{}inVH".format(p))
        mkdir_p(outdir)
        for ext in "png eps root pdf".split():
            c.SaveAs(os.path.join(outdir, "{}.{}".format(disc, ext)))
