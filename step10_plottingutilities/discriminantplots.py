#!/usr/bin/env python

import itertools, os, pprint, re

import numpy as np
import uncertainties

import ROOT

from TemplateBuilder.TemplateBuilder.moremath import weightedaverage

from helperstuff import config, style
from helperstuff.enums import Category
from helperstuff.samples import ReweightingSampleWithPdf, Sample
from helperstuff.templates import Template
from helperstuff.utilities import cache, cache_file, mkdir_p, PlotCopier, TFile

class Tree(object):
  def __init__(self, filename, treename):
    self.__filename = filename
    self.__treename = treename
    self.__histogramcomponentpieces = []
    self.__filled = False

  def __repr__(self):
    return "gettree({filename!r}, {treename!r})".format(filename=self.__filename, treename=self.__treename)
  def __str__(self):
    return self.__filename+":"+self.__treename

  def registerhistogramcomponentpiece(self, piece):
    self.__histogramcomponentpieces.append(piece)

  def fillall(self):
    if self.__filled: return
    print "Iterating through", self
    print "Filling", len(self.__histogramcomponentpieces), "histograms"
    print
    with TFile(self.__filename) as f:
      t = getattr(f, self.__treename)
      t.SetBranchStatus("*", 0)

      for hcp in self.__histogramcomponentpieces:
        hcp.setformulas(t)

      nentries = t.GetEntries()

      for i, entry in enumerate(t, start=1):
        if i % 10000 == 0 or i == nentries: print i, "/", nentries
        for hcp in self.__histogramcomponentpieces:
          hcp.fill()

      for hcp in self.__histogramcomponentpieces:
        hcp.finalize()
    self.__filled = True
    print

class HistogramComponentPiece(object):
  def __init__(self, name, tree, xformula, weightformula, cutformula, binning):
    tree.registerhistogramcomponentpiece(self)
    self.__name = name
    self.__tree = tree
    self.__xformula = xformula
    self.__weightformula = weightformula
    self.__cutformula = cutformula

    self.__histogram = ROOT.TH1F(name, "", len(binning)-1, binning)

    self.__finalized = False

  def setformulas(self, tree):
    for formula in self.__xformula, self.__weightformula, self.__cutformula:
      for branch in re.findall(r"\b[a-zA-Z_][a-zA-Z0-9_]+\b", formula):
        tree.SetBranchStatus(branch, 1)

        try:
          tree.GetEntry(0)
          getattr(tree, branch)
        except AttributeError:
          tree.Show()
          raise ValueError("Bad branch "+branch)

    self.__xtreeformula = ROOT.TTreeFormula(self.__name+"_x", self.__xformula, tree)
    self.__weighttreeformula = ROOT.TTreeFormula(self.__name+"_wt", self.__weightformula, tree)
    self.__cuttreeformula = ROOT.TTreeFormula(self.__name+"_wt", self.__cutformula, tree)

  def fill(self):
    assert not self.__finalized
    if self.__cuttreeformula.EvalInstance(): self.histogram.Fill(self.__xtreeformula.EvalInstance(), self.__weighttreeformula.EvalInstance())

  @property
  def histogram(self): return self.__histogram

  def GetBinContentError(self, x):
    return uncertainties.ufloat(self.histogram.GetBinContent(x), self.histogram.GetBinError(x))

  def finalize(self):
    del self.__xtreeformula, self.__weighttreeformula, self.__cuttreeformula
    self.__finalized = True

  def makehistogram(self):
    self.__tree.fillall()

class HistogramComponent(object):
  def __init__(self, name, trees, xformula, weightformula, cutformula, binning):
    self.__pieces = [
      HistogramComponentPiece(name+"_"+str(i), tree, xformula, weightformula, cutformula, binning) for i, tree in enumerate(trees)
    ]
    self.__histogram = self.__pieces[0].histogram.Clone(name)
    self.__finalized = False

  @property
  def histogram(self): return self.__histogram

  def makehistogram(self):
    if self.__finalized: return
    for piece in self.__pieces:
      piece.makehistogram()
    for x in xrange(1, self.__histogram.GetNbinsX()+1):
      piececontents = [piece.GetBinContentError(x) for piece in self.__pieces]
      piececontents = [piececontent for piececontent in piececontents if piececontent.n or piececontent.s]
      if piececontents:
        self.SetBinContentError(x, weightedaverage(piececontents))
    self.__finalized = True

  def SetBinContentError(self, x, value):
    assert not self.__finalized
    self.__histogram.SetBinContent(x, value.n)
    self.__histogram.SetBinError(x, value.s)

class Histogram(object):
  def __init__(self, name, trees, xformula, weightformulas, cutformula, binning, linecolor, linestyle, linewidth, fillcolor, fillstyle, legendtitle, legendlpf, addonbottom=[]):
    self.__components = [
      HistogramComponent(
        name+"_"+str(i), componenttrees, xformula, weightformula, cutformula, binning
      ) for i, (componenttrees, weightformula) in enumerate(itertools.izip_longest(trees, weightformulas))
    ]
    self.__histogram = self.__components[0].histogram.Clone(name)
    self.__finalized = False

    self.__histogram.SetLineColor(linecolor)
    self.__histogram.SetLineStyle(linestyle)
    self.__histogram.SetLineWidth(linewidth)
    self.__histogram.SetFillColor(fillcolor)
    self.__histogram.SetFillStyle(fillstyle)

    self.__addonbottom = addonbottom

    self.__legendtitle = legendtitle
    self.__legendlpf = legendlpf

  def makefinalhistogram(self):
    if self.__finalized: return
    for component in self.__components:
      component.makehistogram()
      self.__histogram.Add(component.histogram)

    for _ in self.__addonbottom:
      _.makefinalhistogram()
      self.__histogram.Add(_.histogram)

    self.__finalized = True

  @property
  def histogram(self): return self.__histogram

  def addtolegend(self, legend):
    return legend.AddEntry(self.histogram, self.__legendtitle, self.__legendlpf)

@cache
def gettree(*args, **kwargs):
  return Tree(*args, **kwargs)

@cache_file("trees_tmp.pkl")
def gettrees(*productionmodesandhypotheses):
  print "Finding trees for", productionmodesandhypotheses
  return [
    [
      gettree(
        sample.withdiscriminantsfile(),
        "candTree",
      ) for sample in sorted(st, key=lambda x: x.withdiscriminantsfile())
    ]
    for production in config.productionsforcombine
    for otherargs in productionmodesandhypotheses
    for st in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg", *otherargs).reweightfrom()
  ]


def getweights(*productionmodesandhypotheses, **kwargs):
  scaleby = kwargs.pop("scaleby", 1)
  assert not kwargs
  return [
    "({}) * ({}) * ({})".format(production.dataluminosity if "ZX" not in otherargs else 1, scaleby, weight)
    for production in config.productionsforcombine
    for otherargs in productionmodesandhypotheses
    for weight in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg", *otherargs).weightname()
  ]

class HypothesisLine(object):
  def __init__(self, hypothesis, ffHlinecolor, ffHlinestyle, VVHlinecolor, VVHlinestyle, legendname):
    self.hypothesis = hypothesis
    self.ffHlinecolor = ffHlinecolor
    self.ffHlinestyle = ffHlinestyle
    self.VVHlinecolor = VVHlinecolor
    self.VVHlinestyle = VVHlinestyle
    self.legendname = legendname

  def getweights(self, *otherargs):
    SMargs = [tuple(args) + ("0+",) for args in otherargs]
    otherargs = [tuple(args) + (self.hypothesis,) for args in otherargs]
    scaleby = uncertainties.nominal_value(
      (
        sum(ReweightingSampleWithPdf(production, *args).xsec for production in config.productionsforcombine for args in SMargs)
      ) / (
        sum(ReweightingSampleWithPdf(production, *args).xsec for production in config.productionsforcombine for args in otherargs)
      )
    )
    return getweights(*otherargs, scaleby=scaleby)

  def gettrees(self, *otherargs):
    otherargs = (tuple(args) + (self.hypothesis,) for args in otherargs)
    return gettrees(*otherargs)

  @property
  def ffHweights(self):
    return self.getweights(
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )

  @property
  def VVHweights(self):
    return self.getweights(
      ("VBF",),
      ("VH",),
    )
  @property
  def VVHtrees(self):
    return self.gettrees(
      ("VBF",),
      ("VH",),
    )

class Plot(object):
  def __init__(self, name, xtitle, ytitle, hypothesislines, xformula, cutformula, binning, legendargs, legendcolumns, saveasdir, ymax, plotcopier=ROOT):
    self.__name = name
    self.__xtitle = xtitle
    self.__ytitle = ytitle
    self.__ymax = ymax
    self.__legendargs = legendargs
    self.__legendcolumns = legendcolumns
    self.__saveasdir = saveasdir
    self.__plotcopier = plotcopier

    ZXtrees = gettrees(
      ("ZX",),
    )
    ZZtrees = gettrees(
      ("qqZZ",),
      ("ggZZ",),
    )
    ffHtrees = gettrees(
      ("ggH", "0+"),
      ("bbH", "0+"),
      ("ttH", "0+", "Hff0+"),
    )

    ZXhistogram = Histogram(
      name+"_ZX",
      ZXtrees,
      xformula,
      getweights(("ZX",),),
      cutformula,
      binning,
      linecolor=1,
      linestyle=1,
      linewidth=2,
      fillcolor=ROOT.TColor.GetColor("#669966"),
      fillstyle=1001,
      legendtitle="Z+X",
      legendlpf="f",
      addonbottom=[],
    )

    ZZhistogram = Histogram(
      name+"_ZZ",
      ZZtrees,
      xformula,
      getweights(("qqZZ",), ("ggZZ",)),
      cutformula,
      binning,
      linecolor=1,
      linestyle=1,
      linewidth=2,
      fillcolor=ROOT.kAzure-9,
      fillstyle=1001,
      legendtitle="ZZ",
      legendlpf="f",
      addonbottom=[ZXhistogram],
    )

    histograms = [ZZhistogram, ZXhistogram]

    for hypothesis in hypothesislines:
      VVH = Histogram(
        name+"_VVH_"+str(hypothesis.hypothesis),
        hypothesis.VVHtrees,
        xformula,
        hypothesis.VVHweights,
        cutformula,
        binning,
        linecolor=hypothesis.VVHlinecolor,
        linestyle=hypothesis.VVHlinestyle,
        linewidth=2,
        fillcolor=0,
        fillstyle=0,
        legendtitle="VBF+VH "+hypothesis.legendname,
        legendlpf="f",
        addonbottom=[ZZhistogram]
      )

      ffH = Histogram(
        name+"_ffH_"+str(hypothesis.hypothesis),
        ffHtrees,
        xformula,
        hypothesis.ffHweights,
        cutformula,
        binning,
        linecolor=hypothesis.ffHlinecolor,
        linestyle=hypothesis.ffHlinestyle,
        linewidth=2,
        fillcolor=0,
        fillstyle=0,
        legendtitle="Total "+hypothesis.legendname,
        legendlpf="l",
        addonbottom=[VVH]
      )

      histograms += [ffH, VVH]

    self.histograms = histograms

  def makeplot(self):
    c = self.__plotcopier.TCanvas()

    for h in self.histograms:
      h.makefinalhistogram()

    hstack = ROOT.THStack(self.__name, "")
    l = ROOT.TLegend(*self.__legendargs)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetNColumns(self.__legendcolumns)

    for h in self.histograms:
      hstack.Add(h.histogram)
      h.addtolegend(l)

    hstack.Draw("hist nostack")
    hstack.GetXaxis().SetTitle(self.__xtitle)
    hstack.GetYaxis().SetTitle(self.__ytitle)
    hstack.SetMaximum(self.__ymax)
    l.Draw()

    mkdir_p(self.__saveasdir)
    for ext in "png pdf root C".split():
      c.SaveAs(os.path.join(self.__saveasdir, self.__name+"."+ext))

def makeplots():
  with PlotCopier() as pc:
    masscut = "ZZMass>{} && ZZMass<{}".format(config.m4lmin, config.m4lmax)
    categoryname = "category_0P_or_0M_or_a2_or_L1_or_L1Zg"

    untaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("Untagged").idnumbers) + ")"
    VBFtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VBFtagged").idnumbers) + ")"
    HadVHtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VHHadrtagged").idnumbers) + ")"

    untaggedenrichcut = masscut + " && D_bkg > 0.5"
    VBFtaggedenrichcut = masscut + " && D_bkg_VBFdecay > 0.5"
    HadVHtaggedenrichcut = masscut + " && D_bkg_HadVHdecay > 0.5"

    purehypothesislines = [
      HypothesisLine("0+",   2,               1, 2,               2, "SM"),
      HypothesisLine("0-",   4,               1, 4,               2, "f_{a3}=1"),
      HypothesisLine("a2",   ROOT.kGreen+3,   1, ROOT.kGreen+3,   2, "f_{a2}=1"),
      HypothesisLine("L1",   ROOT.kMagenta+2, 1, ROOT.kMagenta+3, 2, "f_{#Lambda1}=1"),
      HypothesisLine("L1Zg", ROOT.kOrange+2,  1, ROOT.kOrange+2,  2, "f_{#Lambda1}^{Z#gamma}=1"),
    ]
    plots = [
      Plot(
        name="D_0minus_decay",
        xtitle="D_{0-}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_0hplus_decay",
        xtitle="D_{0h+}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .5, .7, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1_decay",
        xtitle="D_{#Lambda1}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .55, .8, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1Zg_decay",
        xtitle="D_{#Lambda1}^{Z#gamma,dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .4, .55, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_int_decay",
        xtitle="D_{int}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_int_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([-1, .8, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_0minus_VBFdecay",
        xtitle="D_{0-}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_0hplus_VBFdecay",
        xtitle="D_{0h+}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1_VBFdecay",
        xtitle="D_{#Lambda1}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1Zg_VBFdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .8, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_int_VBFdecay",
        xtitle="D_{int}^{VBF}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_int_VBF",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([-1., 0, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_0minus_HadVHdecay",
        xtitle="D_{0-}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_0hplus_HadVHdecay",
        xtitle="D_{0h+}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1_HadVHdecay",
        xtitle="D_{#Lambda1}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_L1Zg_HadVHdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_int_HadVH",
        xtitle="D_{int}^{VH}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_int_HadVH",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([-1, -.6, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),

      Plot(
        name="D_bkg",
        xtitle="D_{bkg}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=untaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
      ),
      Plot(
        name="D_bkg_VBFdecay",
        xtitle="D_{bkg}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_VBFdecay",
        cutformula=VBFtaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
      Plot(
        name="D_bkg_HadVHdecay",
        xtitle="D_{bkg}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_HadVHdecay",
        cutformula=untaggedcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .5, .8, .9),
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
      ),
    ]

    for plot in plots:
      plot.makeplot()

if __name__ == "__main__":
  makeplots()
