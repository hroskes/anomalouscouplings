#!/usr/bin/env python

debug = False

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  p.add_argument("--filter", type=eval, default="lambda kwargs: True")
  p.add_argument("--debug", action="store_true")
  args = p.parse_args()
  debug = args.debug
  del args.debug

import functools, itertools, os, pprint, re

import numpy as np
import uncertainties

import ROOT

from TemplateBuilder.TemplateBuilder.moremath import weightedaverage

from helperstuff import config, stylefunctions as style
from helperstuff.combinehelpers import getrate
from helperstuff.enums import categories, Category, channels, Hypothesis
from helperstuff.samples import ReweightingSample, ReweightingSampleWithPdf, Sample
from helperstuff.templates import IntTemplate, Template
from helperstuff.utilities import cache, cache_file, cache_keys, mkdir_p, PlotCopier, TFile

masscut = "ZZMass>{} && ZZMass<{}".format(config.m4lmin, config.m4lmax)

class TreeFormula(object):
  def __init__(self, formula, ttree):
    name = re.sub(r"\W", "", formula)
    self.__ttreeformula = ROOT.TTreeFormula(name, formula, ttree)
    self.__value = None

  def EvalInstance(self):
    if self.__value is None: self.__value = self.__ttreeformula.EvalInstance()
    return self.__value

  def nextevent(self):
    self.__value = None

class TTree(object):
  def __init__(self, ttree):
    self.__ttree = ttree
    self.__treeformulas = {}

  def __getattr__(self, attr):
    return getattr(self.__ttree, attr)

  def treeformula(self, formula):
    if formula not in self.__treeformulas: self.__treeformulas[formula] = TreeFormula(formula, self.__ttree)
    return self.__treeformulas[formula]

  def __iter__(self):
    nentries = self.GetEntries()
    for i, entry in enumerate(self.__ttree, start=1):
      if i % 10000 == 0 or i == nentries:
        print i, "/", nentries
        if debug: break
      self.nextevent()
      yield entry

  def nextevent(self):
    for treeformula in self.__treeformulas.itervalues():
      treeformula.nextevent()

class Tree(object):
  def __init__(self, filename, treename):
    self.__filename = filename
    self.__treename = treename
    self.__histogramcomponentpieces = []
    self.__filled = False

  @property
  def isVVH(self):
    if "VBF" in self.__filename or "ZH" in self.__filename or "WH" in self.__filename: return True
    if "ggH" in self.__filename or "ttH" in self.__filename or "bbH" in self.__filename: return False
    assert False, self.__filename

  def __repr__(self):
    return "gettree({filename!r}, {treename!r})".format(filename=self.__filename, treename=self.__treename)
  def __str__(self):
    return self.__filename+":"+self.__treename

  def registerhistogramcomponentpiece(self, piece):
    self.__histogramcomponentpieces.append(piece)

  def fillall(self):
    if self.__filled: return
    print "Iterating through", self
    print "Filling", len(self.__histogramcomponentpieces), "histograms:"
    print ", ".join(_._HistogramComponentPiece__name for _ in self.__histogramcomponentpieces)
    print
    with TFile(self.__filename) as f:
      t = TTree(getattr(f, self.__treename))
      t.SetBranchStatus("*", 0)

      for hcp in self.__histogramcomponentpieces:
        hcp.setformulas(t)

      nentries = t.GetEntries()

      for entry in t:
        for hcp in self.__histogramcomponentpieces:
          hcp.fill()

      for hcp in self.__histogramcomponentpieces:
        hcp.finalize()
    self.__filled = True
    print

class HistogramComponentPiece(object):
  def __init__(self, name, tree, xformula, weightformula, cutformula, binning, mirror=False):
    print "   ", name
    tree.registerhistogramcomponentpiece(self)
    self.__name = name
    self.__tree = tree
    self.__xformula = xformula
    self.__weightformula = weightformula
    self.__cutformula = cutformula
    self.__mirror = mirror
    if mirror: assert len(binning) == 3

    self.__histogram = ROOT.TH1F(name, "", len(binning)-1, binning)

    self.__finalized = False

  def setformulas(self, tree):
    for formula in self.__xformula, self.__weightformula, self.__cutformula:
      for branch in re.findall(r"\b[a-zA-Z_][a-zA-Z0-9_]+\b", formula):
        if branch in ("max", "min"): continue
        tree.SetBranchStatus(branch, 1)

        try:
          tree.GetEntry(0)
          getattr(tree, branch)
        except AttributeError:
          tree.Show()
          raise ValueError("Bad branch "+branch)

    self.__xtreeformula = tree.treeformula(self.__xformula)
    self.__weighttreeformula = tree.treeformula(self.__weightformula)
    self.__cuttreeformula = tree.treeformula(self.__cutformula)

  def fill(self):
    assert not self.__finalized
    if self.__cuttreeformula.EvalInstance():
      x = self.__xtreeformula.EvalInstance()
      wt = self.__weighttreeformula.EvalInstance()

      if self.__mirror:
        wt /= 2
        if x == 0: x += 1e-6
        self.histogram.Fill(-x, wt)

      self.histogram.Fill(x, wt)

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
  def __init__(self, name, trees, xformula, weightformula, cutformula, binning, scaleto, mirror=False, normalizationtrees=None, normalizationweightformula=None, rescaling=None):
    print " ", name
    self.__pieces = [
      HistogramComponentPiece(name+"_"+str(i), tree, xformula, weightformula, cutformula, binning, mirror=mirror) for i, tree in enumerate(trees)
    ]
    self.__histogram = self.__pieces[0].histogram.Clone(name)
    self.__finalized = False
    if normalizationweightformula is None: normalizationweightformula = weightformula
    if normalizationtrees is None: normalizationtrees = trees
    if scaleto is None:
      self.__normalization = None
    else:
      self.__normalization = makehistogramcomponentnormalization(name=name+"_normalization", trees=normalizationtrees, weightformula=normalizationweightformula, scaleto=scaleto)
    self.__rescaling = rescaling

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

    toprint = [self.__histogram.GetName(), self.__histogram.Integral()]
    if self.__normalization is not None:
      self.__histogram.Scale(self.__normalization.scaleby)
    toprint += [self.__histogram.Integral()]
    if self.__rescaling is not None:
      self.__histogram.Scale(self.__rescaling.scaleby)
    toprint += [self.__histogram.Integral()]
    for _ in toprint: print _,
    print

  def SetBinContentError(self, x, value):
    assert not self.__finalized
    self.__histogram.SetBinContent(x, value.n)
    self.__histogram.SetBinError(x, value.s)

class HistogramComponentNormalization(HistogramComponent):
  def __init__(self, name, trees, weightformula, scaleto):
    super(HistogramComponentNormalization, self).__init__(
      name=name, trees=trees, weightformula=weightformula,
      xformula="1", cutformula=masscut, binning=np.array([0., 2.]), mirror=False, scaleto=None
    )

    self.__scaleto = scaleto

  @property
  def scaleby(self):
    self.makehistogram()
    return self.__scaleto / self.histogram.Integral()

@cache_keys(
  name=lambda name: None,
  trees=tuple,
)
def makehistogramcomponentnormalization(**kwargs): return HistogramComponentNormalization(**kwargs)

class Histogram(object):
  def __init__(self, name, trees, xformula, weightformulas, scaletos, cutformula, binning, linecolor, linestyle, linewidth, fillcolor, fillstyle, legendname, legendlpf, addonbottom, mirror, normalizationtrees=None, normalizationweightformulas=None, makegraph=False, markercolor=None, markerstyle=None, markersize=None, rescalings=None, inmainlegend=False):
    print "initing", name
    if normalizationweightformulas is None: normalizationweightformulas = weightformulas
    if normalizationtrees is None: normalizationtrees = trees
    if rescalings is None: rescalings = [None for _ in trees]
    if scaletos is None: scaletos = [None for _ in trees]
    assert len(scaletos) == len(weightformulas) == len(trees) == len(normalizationtrees) == len(normalizationweightformulas) == len(rescalings), (name, len(scaletos), len(weightformulas), len(trees), len(normalizationtrees), len(normalizationweightformulas), len(rescalings))
    self.__components = [
      HistogramComponent(
        name=name+"_"+str(i), trees=componenttrees, xformula=xformula, weightformula=weightformula,
        cutformula=cutformula, binning=binning, mirror=mirror, scaleto=scaleto,
        normalizationtrees=componentnormalizationtrees, normalizationweightformula=normalizationweightformula,
        rescaling=rescaling,
      ) for i, (componenttrees, weightformula, scaleto, componentnormalizationtrees, normalizationweightformula, rescaling) in enumerate(itertools.izip_longest(trees, weightformulas, scaletos, normalizationtrees, normalizationweightformulas, rescalings))
    ]
    self.__histogram = self.__components[0].histogram.Clone(name)
    self.__finalized = False

    self.__histogram.SetLineColor(linecolor)
    self.__histogram.SetLineStyle(linestyle)
    self.__histogram.SetLineWidth(linewidth)
    self.__histogram.SetFillColor(fillcolor)
    self.__histogram.SetFillStyle(fillstyle)

    self.__addonbottom = addonbottom

    self.__legendname = legendname
    self.__legendlpf = legendlpf

    self.__makegraph = makegraph
    if self.__makegraph:
      self.__markerstyle = {"color": markercolor, "style": markerstyle, "size": markersize}
      assert all(_ is not None for _ in self.__markerstyle.itervalues()), self.__markerstyle

    self.__inmainlegend = inmainlegend

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
  @cache
  def graph(self):
    if not self.__makegraph: raise ValueError("Call histogram, not graph")
    assert self.__finalized
    g = style.asymmerrorsfromhistogram(self.__histogram, showemptyerrors=False)
    g.SetMarkerColor(self.__markerstyle["color"])
    g.SetMarkerStyle(self.__markerstyle["style"])
    g.SetMarkerSize(self.__markerstyle["size"])
    g.SetLineColor(self.__histogram.GetLineColor())
    g.SetLineStyle(self.__histogram.GetLineStyle())
    g.SetLineWidth(self.__histogram.GetLineWidth())
    return g

  @property
  def histogram(self):
    if self.__makegraph: raise ValueError("Call graph, not histogram")
    return self.__histogram

  def addtolegend(self, legend):
    return legend.AddEntry(self.graph if self.__makegraph else self.histogram, self.__legendname, self.__legendlpf)

  @property
  def inmainlegend(self):
    return self.__inmainlegend

@cache_keys(
  name=lambda name: None,
  trees=lambda trees: tuple(tuple(_) for _ in trees),
  addonbottom=tuple,
  weightformulas=tuple,
  binning=tuple,
  normalizationtrees=lambda normalizationtrees: tuple(tuple(_) for _ in normalizationtrees),
  normalizationweightformulas=tuple,
  scaletos=lambda scaletos: tuple(scaletos) if scaletos is not None else scaletos,
  rescalings=tuple,
)
def makehistogram(**kwargs): return Histogram(**kwargs)

class HistogramNormalization(Histogram):
  def __init__(self, name, trees, weightformulas, scaletos, normalizationtrees=None, normalizationweightformulas=None):
    super(HistogramNormalization, self).__init__(
      name=name, trees=trees, xformula="1", weightformulas=weightformulas,
      scaletos=scaletos, cutformula=masscut, binning=np.array([0., 2.]),
      linecolor=0, linestyle=0, linewidth=0, fillcolor=0, fillstyle=0, legendname="", legendlpf="", addonbottom=[], mirror=None,
      normalizationtrees=normalizationtrees, normalizationweightformulas=normalizationweightformulas,
    )

@cache_keys(
  name=lambda name: None,
  trees=lambda trees: tuple(tuple(_) for _ in trees),
  weightformulas=tuple,
  normalizationtrees=lambda normalizationtrees: tuple(tuple(_) for _ in normalizationtrees),
  normalizationweightformulas=tuple,
  scaletos=lambda scaletos: tuple(scaletos) if scaletos is not None else scaletos,
)
def makehistogramnormalization(**kwargs): return HistogramNormalization(**kwargs)

class HistogramRescaling(object):
  def __init__(self, numerator, denominator):
    self.__numerator = numerator
    self.__denominator = denominator
  @property
  def scaleby(self):
    self.__numerator.makefinalhistogram()
    self.__denominator.makefinalhistogram()
    return self.__numerator.histogram.Integral() / self.__denominator.histogram.Integral()

class MuRescaling(object):
  def __init__(self, mu):
    self.__mu = mu
  @property
  def scaleby(self):
    return self.__mu

@cache
def makehistogramrescaling(**kwargs): return HistogramRescaling(**kwargs)

@cache
def gettree(*args, **kwargs):
  return Tree(*args, **kwargs)

@cache_file(os.path.join(config.repositorydir, "data", "trees.pkl"))
def getfilenames(productions, *productionmodesandhypotheses):
  print "Finding trees for", productionmodesandhypotheses
  return [
    [
      sample.withdiscriminantsfile()
      for sample in sorted(st, key=lambda x: x.withdiscriminantsfile())
    ]
    for production in productions
    for otherargs in productionmodesandhypotheses
    for st in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg_morecategories", *otherargs).reweightfrom()
  ]

def gettrees(*productionmodesandhypotheses, **kwargs):
  result = [
    [
      gettree(
        filename,
        "candTree",
      ) for filename in lst
    ] for lst in getfilenames(tuple(config.productionsforcombine), *productionmodesandhypotheses, **kwargs)
  ]
  if debug: result = [[_[0]] for _ in result]
  return result

@cache_file(os.path.join(config.repositorydir, "data", "weights.pkl"))
def getweights(*productionmodesandhypotheses, **kwargs):
  print "Finding weights for", productionmodesandhypotheses
  scalebys = kwargs.pop("scalebys", [1 for _ in productionmodesandhypotheses])
  uselumi = kwargs.pop("uselumi", False)
  assert not kwargs
  return [
    "({}) * ({}) * ({})".format(weight, scaleby, production.dataluminosity if uselumi else 1)
    for production in config.productionsforcombine
    for otherargs, scaleby in itertools.izip_longest(productionmodesandhypotheses, scalebys)
    for weight in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg_morecategories", *otherargs).weightname()
  ]

def getscaletos(*productionmodes, **kwargs):
  isL1Zg = kwargs.pop("isL1Zg", False)
  assert not kwargs
  if isL1Zg: productionmodes = ["ZH" if p == "VH" else p for p in productionmodes]
  return [
    sum(
      getrate(
        ch,
        ca,
        "fordata",
        production,
        "fa3fa2fL1fL1Zg_morecategories",
        (
          productionmode
          if productionmode != "VH"
          else weight.split("p_Gen_")[1].split("_")[0]
        )
      )
      for ch in channels
      for ca in categories
    )
    for production in config.productionsforcombine
    for productionmode in productionmodes
    for weight in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg_morecategories", productionmode, "0+", "Hff0+" if productionmode == "ttH" else None).weightname()
  ]

class HypothesisLine(object):
  def __init__(self, hypothesis, ffHlinecolor, ffHlinestyle, VVHlinecolor, VVHlinestyle, legendname, muV=None, muf=None):
    self.hypothesis = hypothesis
    self.ffHlinecolor = ffHlinecolor
    self.ffHlinestyle = ffHlinestyle
    self.VVHlinecolor = VVHlinecolor
    self.VVHlinestyle = VVHlinestyle
    self.legendname = legendname
    self.muV = muV
    self.muf = muf
    if not (muV is muf is None) and (muV is None or muf is None):
      raise ValueError("Have to provide muV and muf, or neither")

  @property
  def ispure(self):
    return Hypothesis(self.hypothesis).ispure

  def getweights(self, *otherargs):
    scalebys = []
    for args in otherargs:
      productionmode = args[0]
      if self.muV is self.muf is None:
        productionmodes = {
          "ggH": ("ggH", "bbH", "ttH"),
          "bbH": ("ggH", "bbH", "ttH"),
          "ttH": ("ggH", "bbH", "ttH"),
          "VBF": ("VBF", "VH"),
          "VH": ("VBF", "VH"),
          "ZH": ("VBF", "VH"),
          "WH": ("VBF", "VH"),
        }[productionmode]
        numerator = sum(
          Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "0+", pm).gettemplate().Integral() * production.dataluminosity
          for ca in categories
          for ch in channels
          for pm in productionmodes
          for production in config.productionsforcombine
        )
        if Hypothesis(self.hypothesis).ispure:
          denominator = sum(
            Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, self.hypothesis, pm).gettemplate().Integral() * production.dataluminosity
            for ca in categories
            for ch in channels
            for pm in productionmodes
            for production in config.productionsforcombine
          )
        else:
          s = ReweightingSample("ggH", self.hypothesis)
          g1 = s.g1
          g2 = s.g2
          g4 = s.g4
          gL1 = s.g1prime2 / 1e4
          gL1Zg = s.ghzgs1prime2 / 1e4

          s = ReweightingSample("VBF", self.hypothesis)
          assert g1 == s.g1
          assert g2 == s.g2
          assert g4 == s.g4
          assert gL1 == s.g1prime2 / 1e4
          assert gL1Zg == s.ghzgs1prime2 / 1e4

          if productionmode in ("ggH", "bbH", "ttH"):
            maxpower = 2
          elif productionmode in ("VBF", "VH", "ZH", "WH"):
            maxpower = 4
          def inttype(power1, poweri, powerj, powerk, powerl):
            assert power1 + poweri + powerj + powerk + powerl == maxpower, (power1, poweri, powerj, powerk, powerl)
            assert min(power1, poweri, powerj, powerk, powerl) >= 0, (power1, poweri, powerj, powerk, powerl)
            result = ""
            if power1: result += "g1{}".format(power1)
            if poweri: result += "gi{}".format(poweri)
            if powerj: result += "gj{}".format(powerj)
            if powerk: result += "gk{}".format(powerk)
            if powerl: result += "gl{}".format(powerl)
            return result
          denominator = sum(
            (
                 Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "0+",         pm).gettemplate().Integral() * g1   **maxpower
            +    Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "0-",         pm).gettemplate().Integral() * g4   **maxpower
            +    Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "a2",         pm).gettemplate().Integral() * g2   **maxpower
            +    Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "L1",         pm).gettemplate().Integral() * gL1  **maxpower
            +    Template(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, "L1Zg",       pm).gettemplate().Integral() * gL1Zg**maxpower
            + sum(
              IntTemplate(production, "fa3fa2fL1fL1Zg_morecategories", ca, ch, inttype(power1, poweri, powerj, powerk, powerl), pm).gettemplate().Integral()
              * g1**power1
              * g4**poweri
              * g2**powerj
              * gL1**powerk
              * gL1Zg**powerl
              for power1 in range(maxpower+1)
              for poweri in range(maxpower+1-power1)
              for powerj in range(maxpower+1-power1-poweri)
              for powerk in range(maxpower+1-power1-poweri-powerj)
              for powerl in (maxpower-power1-poweri-powerj-powerk,)
              if power1 != maxpower and poweri != maxpower and powerj != maxpower and powerk != maxpower and powerl != maxpower
              and (
                g1**power1
              * g4**poweri
              * g2**powerj
              * gL1**powerk
              * gL1Zg**powerl
              )
            )
            ) * production.dataluminosity
            for ca in categories
            for ch in channels
            for pm in productionmodes
            for production in config.productionsforcombine
          )

        scalebys.append(numerator/denominator)

      else:
        assert self.muf is not None is not self.muV
        scalebys.append({
          "ggH": self.muf,
          "ttH": self.muf,
          "bbH": self.muf,
          "VH": self.muV,
          "VBF": self.muV,
        }[productionmode])

    otherargs = [tuple(args) + (self.hypothesis,) for args in otherargs]
    assert len(otherargs) == len(scalebys)
    return getweights(*otherargs, scalebys=tuple(scalebys))

  def gettrees(self, *otherargs, **kwargs):
    hypothesis = self.hypothesis
    otherargs = (tuple(args) + (hypothesis,) for args in otherargs)
    return gettrees(*otherargs, **kwargs)

  def ffHweights(self, isL1Zg=False):
    return self.getweights(
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )

  def VVHweights(self, isL1Zg=False):
    return self.getweights(
      ("VBF",),
      ("VH" if not isL1Zg else "ZH",),
    )
  def VVHtrees(self, isL1Zg=False, **kwargs):
    return self.gettrees(
      ("VBF",),
      ("VH" if not isL1Zg else "ZH",),
      **kwargs
    )
  def VBFweights(self, isL1Zg=False, **kwargs):
    return self.getweights(
      ("VBF",),
    )
  def VBFtrees(self, isL1Zg=False, **kwargs):
    return self.gettrees(
      ("VBF",),
      **kwargs
    )
  def notVBFweights(self, isL1Zg=False, **kwargs):
    return self.getweights(
      ("VH" if not isL1Zg else "ZH",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
      **kwargs
    )
  def notVBFtrees(self, isL1Zg=False, **kwargs):
    return self.gettrees(
      ("VH" if not isL1Zg else "ZH",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
      **kwargs
    )
  def VHweights(self, isL1Zg=False, **kwargs):
    return self.getweights(
      ("VH" if not isL1Zg else "ZH",),
      **kwargs
    )
  def VHtrees(self, isL1Zg=False, **kwargs):
    return self.gettrees(
      ("VH" if not isL1Zg else "ZH",),
      **kwargs
    )
  def notVHweights(self, isL1Zg=False, **kwargs):
    return self.getweights(
      ("VBF",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
      **kwargs
    )
  def notVHtrees(self, isL1Zg=False, **kwargs):
    return self.gettrees(
      ("VBF",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )

class Plot(object):
  mainlegendargs = .2, .1, .9, .9
  def __init__(self, **kwargs):
    name = kwargs.pop("name")
    xtitle = kwargs.pop("xtitle")
    ytitle = kwargs.pop("ytitle")
    hypothesislines = kwargs.pop("hypothesislines")
    xformula = kwargs.pop("xformula")
    cutformula = kwargs.pop("cutformula")
    binning = kwargs.pop("binning")
    legendargs = kwargs.pop("legendargs")
    legendcolumns = kwargs.pop("legendcolumns")
    categorylabel = kwargs.pop("categorylabel")
    Dbkglabel = kwargs.pop("Dbkglabel")
    saveasdir = kwargs.pop("saveasdir")
    ymax = kwargs.pop("ymax")
    plotcopier = kwargs.pop("plotcopier")
    CMStext = kwargs.pop("CMStext")
    isDCP = kwargs.pop("isDCP", None)
    iscategorydiscriminant = kwargs.pop("iscategorydiscriminant", None)
    assert not kwargs, kwargs

    self.__name = name
    self.__xtitle = xtitle
    self.__ytitle = ytitle
    self.__ymax = ymax
    self.__legendargs = legendargs
    self.__legendcolumns = legendcolumns
    self.__categorylabel = categorylabel
    self.__Dbkglabel = Dbkglabel
    self.__saveasdir = saveasdir
    self.__plotcopier = plotcopier
    self.__iscategorydiscriminant = iscategorydiscriminant

    assert isDCP in (None, "prod", "dec")
    assert iscategorydiscriminant in (None, "VBF", "VH")

    datatrees = gettrees(
      ("data",),
    )
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

    self.histograms = histograms = []
    self.graphs = graphs = []

    datahistogram = makehistogram(
      name=name+"_data",
      trees=datatrees,
      xformula=xformula,
      weightformulas=["1" for production in config.productionsforcombine],
      cutformula=cutformula,
      binning=binning,
      linecolor=1,
      linestyle=1,
      linewidth=2,
      fillcolor=ROOT.TColor.GetColor("#669966"),
      fillstyle=1001,
      markercolor=1,
      markerstyle=20,
      markersize=1.2,
      legendname="data",
      legendlpf="lp",
      addonbottom=[],
      mirror=False,
      scaletos=None,
      makegraph=True,
      inmainlegend=True,
    )

    ZXhistogram = makehistogram(
      name=name+"_ZX",
      trees=ZXtrees,
      xformula=xformula,
      weightformulas=getweights(("ZX",),),
      cutformula=cutformula,
      binning=binning,
      linecolor=1,
      linestyle=1,
      linewidth=2,
      fillcolor=ROOT.TColor.GetColor("#669966"),
      fillstyle=1001,
      legendname="Z+X",
      legendlpf="f",
      addonbottom=[],
      mirror=isDCP is not None,
      scaletos=None,
      inmainlegend=True,
    )

    ZZhistogram = makehistogram(
      name=name+"_ZZ",
      trees=ZZtrees,
      xformula=xformula,
      weightformulas=getweights(("qqZZ",), ("ggZZ",), uselumi=True),
      cutformula=cutformula,
      binning=binning,
      linecolor=1,
      linestyle=1,
      linewidth=2,
      fillcolor=ROOT.kAzure-9,
      fillstyle=1001,
      legendname="ZZ/Z#gamma*",
      legendlpf="f",
      addonbottom=[ZXhistogram],
      mirror=isDCP is not None,
      scaletos=None,
      inmainlegend=True,
    )

    SMhypothesis = [_ for _ in hypothesislines if Hypothesis(_.hypothesis) == "0+"][0]

    SMVVHweights = SMhypothesis.VVHweights
    SMffHweights = SMhypothesis.ffHweights
    SMVBFweights = SMhypothesis.VBFweights
    SMnotVBFweights = SMhypothesis.notVBFweights
    SMVHweights = SMhypothesis.VHweights
    SMnotVHweights = SMhypothesis.notVHweights

    SMVVHtrees = SMhypothesis.VVHtrees
    SMVBFtrees = SMhypothesis.VBFtrees
    SMnotVBFtrees = SMhypothesis.notVBFtrees
    SMVHtrees = SMhypothesis.VHtrees
    SMnotVHtrees = SMhypothesis.notVHtrees

    VVHscaletos = functools.partial(getscaletos, "VBF", "VH")
    ffHscaletos = functools.partial(getscaletos, "ggH", "bbH", "ttH")
    VBFscaletos = functools.partial(getscaletos, "VBF")
    notVBFscaletos = functools.partial(getscaletos, "VH", "ggH", "bbH", "ttH")
    VHscaletos = functools.partial(getscaletos, "VH")
    notVHscaletos = functools.partial(getscaletos, "VBF", "ggH", "bbH", "ttH")

    SMffHhistogramnormalization = makehistogramnormalization(
      name="SMffHnormalization",
      trees=ffHtrees,
      weightformulas=SMffHweights(),
      scaletos=ffHscaletos(),
    )
    SMVVHhistogramnormalization = makehistogramnormalization(
      name="SMVVHnormalization",
      trees=SMVVHtrees(),
      weightformulas=SMVVHweights(),
      scaletos=VVHscaletos(),
    )

    for hypothesis in hypothesislines:
      isL1Zg = Hypothesis(hypothesis.hypothesis) == "L1Zg"

      ffHhistogramnormalization = makehistogramnormalization(
        name=str(hypothesis.hypothesis)+"ffHnormalization",
        trees=ffHtrees,
        weightformulas=hypothesis.ffHweights(isL1Zg=isL1Zg),
        scaletos=ffHscaletos(isL1Zg=isL1Zg),
        normalizationweightformulas=SMffHweights(isL1Zg=isL1Zg),
      )
      VVHhistogramnormalization = makehistogramnormalization(
        name=str(hypothesis.hypothesis)+"VVHnormalization",
        trees=hypothesis.VVHtrees(isL1Zg=isL1Zg),
        weightformulas=hypothesis.VVHweights(isL1Zg=isL1Zg),
        normalizationweightformulas=SMVVHweights(isL1Zg=isL1Zg),
        normalizationtrees=SMVVHtrees(isL1Zg=isL1Zg),
        scaletos=VVHscaletos(isL1Zg=isL1Zg),
      )

      if hypothesis.muV is hypothesis.muf is None:
        ffHrescaling = makehistogramrescaling(numerator=SMffHhistogramnormalization, denominator=ffHhistogramnormalization)
        VVHrescaling = makehistogramrescaling(numerator=SMVVHhistogramnormalization, denominator=VVHhistogramnormalization)
      else:
        assert hypothesis.muV is not None is not hypothesis.muf
        ffHrescaling = VVHrescaling = None

      VVHkwargs = dict(
        name=name+"_VVH_"+str(hypothesis.hypothesis),
        trees=hypothesis.VVHtrees(),
        xformula=xformula,
        weightformulas=hypothesis.VVHweights(),
        normalizationweightformulas=SMVVHweights(isL1Zg=isL1Zg),
        normalizationtrees=SMVVHtrees(isL1Zg=isL1Zg),
        scaletos=VVHscaletos(isL1Zg=isL1Zg),
        cutformula=cutformula,
        binning=binning,
        linecolor=hypothesis.VVHlinecolor,
        linestyle=hypothesis.VVHlinestyle,
        linewidth=2,
        fillcolor=0,
        fillstyle=0,
        legendname="VBF+VH       "+hypothesis.legendname,
        legendlpf="f",
        addonbottom=[ZZhistogram],
        mirror=isDCP is not None and hypothesis.ispure,
        rescalings=[VVHrescaling for _ in hypothesis.VVHtrees()],
        inmainlegend=hypothesis.hypothesis in ("0+", "0-", "a2", "L1", "L1Zg", "BestFit19009"),
      )

      ffHkwargs = dict(
        name=name+"_ffH_"+str(hypothesis.hypothesis),
        trees=ffHtrees,
        xformula=xformula,
        weightformulas=hypothesis.ffHweights(),
        normalizationweightformulas=SMffHweights(isL1Zg=isL1Zg),
        scaletos=ffHscaletos(isL1Zg=isL1Zg),
        cutformula=cutformula,
        binning=binning,
        linecolor=hypothesis.ffHlinecolor,
        linestyle=hypothesis.ffHlinestyle,
        linewidth=2,
        fillcolor=0,
        fillstyle=0,
#        legendname="Total "+hypothesis.legendname,
        legendname="Total",
        legendlpf="l",
        mirror=isDCP == "prod" or (isDCP == "dec" and hypothesis.ispure),
        rescalings=[ffHrescaling for _ in ffHtrees],
        inmainlegend=hypothesis.hypothesis in ("0+", "0-", "a2", "L1", "L1Zg", "BestFit19009"),
      )

      if iscategorydiscriminant is None:
        pass
      elif iscategorydiscriminant == "VBF":
        VVHkwargs["trees"] = hypothesis.VBFtrees()
        ffHkwargs["trees"] = hypothesis.notVBFtrees()
        VVHkwargs["weightformulas"] = hypothesis.VBFweights()
        ffHkwargs["weightformulas"] = hypothesis.notVBFweights()
        VVHkwargs["scaletos"] = VBFscaletos(isL1Zg=isL1Zg)
        ffHkwargs["scaletos"] = notVBFscaletos(isL1Zg=isL1Zg)
        VVHkwargs["normalizationtrees"] = SMVBFtrees(isL1Zg=isL1Zg)
        ffHkwargs["normalizationtrees"] = SMnotVBFtrees(isL1Zg=isL1Zg)
        VVHkwargs["normalizationweightformulas"] = SMVBFweights(isL1Zg=isL1Zg)
        ffHkwargs["normalizationweightformulas"] = SMnotVBFweights(isL1Zg=isL1Zg)
        VVHkwargs["legendname"] = VVHkwargs["legendname"].replace("VBF+VH", "VBF")
        VVHkwargs["rescalings"] = [VVHrescaling if trees[0].isVVH else ffHrescaling for trees in hypothesis.VBFtrees()]
        ffHkwargs["rescalings"] = [VVHrescaling if trees[0].isVVH else ffHrescaling for trees in hypothesis.notVBFtrees()]
      elif iscategorydiscriminant == "VH":
        VVHkwargs["trees"] = hypothesis.VHtrees()
        ffHkwargs["trees"] = hypothesis.notVHtrees()
        VVHkwargs["weightformulas"] = hypothesis.VHweights()
        ffHkwargs["weightformulas"] = hypothesis.notVHweights()
        VVHkwargs["scaletos"] = VHscaletos(isL1Zg=isL1Zg)
        ffHkwargs["scaletos"] = notVHscaletos(isL1Zg=isL1Zg)
        VVHkwargs["normalizationtrees"] = SMVHtrees(isL1Zg=isL1Zg)
        ffHkwargs["normalizationtrees"] = SMnotVHtrees(isL1Zg=isL1Zg)
        VVHkwargs["normalizationweightformulas"] = SMVHweights(isL1Zg=isL1Zg)
        ffHkwargs["normalizationweightformulas"] = SMnotVHweights(isL1Zg=isL1Zg)
        VVHkwargs["legendname"] = VVHkwargs["legendname"].replace("VBF+VH", "VH")
        VVHkwargs["rescalings"] = [VVHrescaling if trees[0].isVVH else ffHrescaling for trees in hypothesis.VHtrees()]
        ffHkwargs["rescalings"] = [VVHrescaling if trees[0].isVVH else ffHrescaling for trees in hypothesis.notVHtrees()]
      else:
        assert False

      VVH = makehistogram(**VVHkwargs)
      ffHkwargs["addonbottom"] = [VVH]
      ffH = makehistogram(**ffHkwargs)

      histograms += [ffH, VVH]

    histograms += [ZZhistogram, ZXhistogram]

    graphs.append(datahistogram)

    self.__CMStext = CMStext

  def makeplot(self):
    c = self.__plotcopier.TCanvas("c_"+self.__name, "",  8, 30, 800, 800)
    style.applycanvasstyle(c)
    c.SetBottomMargin(0.14)

    for h in self.histograms+self.graphs:
      print self, h
      h.makefinalhistogram()

    hstack = ROOT.THStack(self.__name, "")
    l = ROOT.TLegend(*self.__legendargs)
    style.applylegendstyle(l)
    l.SetNColumns(self.__legendcolumns)

    for h in self.histograms:
      hstack.Add(h.histogram)
      if not h.inmainlegend:
        h.addtolegend(l)

    if self.__iscategorydiscriminant:
      dummycolor = 1
      solidstyle = 1
      dummywidth = 2
      dashedstyle = {h.histogram.GetLineStyle() for h in self.histograms if h.histogram.GetLineStyle() != 1}.pop()

      dummysolid = ROOT.TH1F("dummy_"+self.__name+"solid", "", 1, 0, 1)
      dummysolid.SetLineColor(dummycolor)
      dummysolid.SetLineWidth(dummywidth)
      dummysolid.SetLineStyle(solidstyle)
      l.AddEntry(dummysolid, "Total", "l")

      dummydashed = ROOT.TH1F("dummy_"+self.__name+"dashed", "", 1, 0, 1)
      dummydashed.SetLineColor(dummycolor)
      dummydashed.SetLineWidth(dummywidth)
      dummydashed.SetLineStyle(dashedstyle)
      l.AddEntry(dummydashed, self.__iscategorydiscriminant, "l")

    hstack.Draw("hist nostack")
    style.applyaxesstyle(hstack)
    hstack.GetXaxis().SetTitle(self.__xtitle)
    hstack.GetYaxis().SetTitle(self.__ytitle)
    hstack.SetMaximum(self.__ymax)

    for g in self.graphs:
      g.graph.Draw("PE")
      if not g.inmainlegend:
        g.addtolegend(l)

    l.Draw()

    lumi = sum(production.dataluminosity for production in config.productionsforcombine)

    style.CMS(self.__CMStext, lumi)

    """
    HVV4ltext = ROOT.TPaveText(.2, .87, .4, .92, "brNDC")
    HVV4ltext.SetBorderSize(0)
    HVV4ltext.SetTextAlign(12)
    HVV4ltext.SetTextSize(0.045)
    HVV4ltext.SetFillStyle(0)
    HVV4ltext.SetTextFont(42)
    HVV4ltext.AddText(0.01,0.01,"HVV")
    HVV4ltext.Draw()

    if self.__categorylabel is not None:
      categorytext = ROOT.TPaveText(.4, .87, .65, .92, "brNDC")
      categorytext.SetBorderSize(0)
      categorytext.SetTextAlign(12)
      categorytext.SetTextSize(0.045)
      categorytext.SetFillStyle(0)
      categorytext.SetTextFont(42)
      categorytext.AddText(0.01,0.01,self.__categorylabel)
      categorytext.Draw()

    if self.__Dbkglabel is not None:
      Dbkgtext = ROOT.TPaveText(.65, .87, .9, .92, "brNDC")
      Dbkgtext.SetBorderSize(0)
      Dbkgtext.SetTextAlign(12)
      Dbkgtext.SetTextSize(0.045)
      Dbkgtext.SetFillStyle(0)
      Dbkgtext.SetTextFont(42)
      Dbkgtext.AddText(0.01,0.01,self.__Dbkglabel)
      Dbkgtext.Draw()
    """

    textparts = ["HVV", "H#rightarrow4l"]
    if self.__categorylabel is not None:
      textparts.append(self.__categorylabel)
    if self.__Dbkglabel is not None:
      textparts.append(self.__Dbkglabel)

    labeltext = ROOT.TPaveText(.2, .87, .4, .92, "brNDC")
    labeltext.SetBorderSize(0)
    labeltext.SetTextAlign(12)
    labeltext.SetTextSize(0.045)
    labeltext.SetFillStyle(0)
    labeltext.SetTextFont(42)
    labeltext.AddText(0.01,0.01,", ".join(textparts))
    labeltext.Draw()

    if self.__iscategorydiscriminant:
      line = ROOT.TLine(0.5, 0, 0.5, self.__ymax * .4)
      line.SetLineWidth(4)
      line.SetLineColor(ROOT.kViolet-1)
      line.SetLineStyle(9)
      line.Draw()

    mkdir_p(self.__saveasdir)
    for ext in "png pdf root C".split():
      c.SaveAs(os.path.join(self.__saveasdir, self.__name+"."+ext))

    self.makemainlegend()

  def makemainlegend(self):
    c = self.__plotcopier.TCanvas("c_"+self.__name+"_legend", "",  8, 30, 800, 800)
    style.applycanvasstyle(c)
    l = ROOT.TLegend(*self.mainlegendargs)
    style.applylegendstyle(l)
    l.SetNColumns(self.__legendcolumns)

    for h in self.histograms:
      if h.inmainlegend:
        h.addtolegend(l)
    
    for g in self.graphs:
      g.makefinalhistogram()
      if g.inmainlegend:
        g.addtolegend(l)

    l.Draw()

    mkdir_p(self.__saveasdir)
    for ext in "png pdf root C".split():
      c.SaveAs(os.path.join(self.__saveasdir, "mainlegend."+ext))

categoryname = "category_0P_or_0M_or_a2_or_L1_or_L1Zg"

dijetcut = masscut + " && D_0minus_VBFdecay > -998"
untaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("Untagged").idnumbers) + ")"
VBFtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VBFtagged").idnumbers) + ")"
HadVHtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VHHadrtagged").idnumbers) + ")"
boostedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("Boosted").idnumbers) + ")"
VBF1jcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VBF1jtagged").idnumbers) + ")"
LepVHcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VHLepttagged").idnumbers) + ")"

dijetenrichcut = dijetcut + " && D_bkg > 0.7"
untaggedenrichcut = untaggedcut + " && D_bkg > 0.7"
VBFtaggedenrichcut = VBFtaggedcut + " && D_bkg_VBFdecay > 0.2"
HadVHtaggedenrichcut = HadVHtaggedcut + " && D_bkg_HadVHdecay > 0.2"
boostedenrichcut = boostedcut + " && D_bkg > 0.7"
VBF1jenrichcut = VBF1jcut + " && D_bkg > 0.7"
LepVHenrichcut = LepVHcut + " && D_bkg > 0.7"

purehypothesislines = [
  HypothesisLine("0+",   2,               1, 2,               2, "SM"),
  HypothesisLine("0-",   4,               1, 4,               2, "f_{a3}=1"),
  HypothesisLine("a2",   ROOT.kGreen+3,   1, ROOT.kGreen+3,   2, "f_{a2}=1"),
  HypothesisLine("L1",   ROOT.kMagenta+3, 1, ROOT.kMagenta+3, 2, "f_{#Lambda1}=1"),
  HypothesisLine("L1Zg", ROOT.kOrange+2,  1, ROOT.kOrange+2,  2, "f_{#Lambda1}^{Z#gamma}=1"),
]

a3mixdecay = HypothesisLine("fa30.5",  ROOT.kAzure-4, 1, ROOT.kAzure-4, 2, "f_{a3}=0.5")
a3mixVBF = HypothesisLine("fa3VBF0.5", ROOT.kAzure-4, 1, ROOT.kAzure-4, 2, "f_{a3}^{VBF}=0.5")
a3mixVH = HypothesisLine("fa3VH0.5",   ROOT.kAzure-4, 1, ROOT.kAzure-4, 2, "f_{a3}^{VH}=0.5")

a2mixdecay = HypothesisLine("fa2-0.5", ROOT.kGreen-3, 1, ROOT.kGreen-3, 2, "f_{a2}=#minus0.5")
a2mixVBF = HypothesisLine("fa2VBF0.5", ROOT.kGreen-3, 1, ROOT.kGreen-3, 2, "f_{a2}^{VBF}=0.5")
a2mixVH = HypothesisLine("fa2VH0.5",   ROOT.kGreen-3, 1, ROOT.kGreen-3, 2, "f_{a2}^{VH}=0.5")

L1mixdecay = HypothesisLine("fL10.5", ROOT.kMagenta-4, 1, ROOT.kMagenta-4, 2, "f_{#Lambda1}=0.5")

BestFit19009 = HypothesisLine("BestFit19009", 1, 1, 1, 2, "Best fit", muV=1.05094e-10, muf=5.70223)

def makeplots(filter):
  with PlotCopier() as pc:
    plotkwargses = [
      dict(
        name="D_0minus_decay",
        xtitle="D_{0-}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=130,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_0hplus_decay",
        xtitle="D_{0h+}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixdecay],
        xformula="D_0hplus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .5, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=150,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1_decay",
        xtitle="D_{#Lambda1}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [L1mixdecay],
        xformula="D_L1_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .55, .8, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=150,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1Zg_decay",
        xtitle="D_{#Lambda1}^{Z#gamma,dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .4, .55, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=140,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_CP_decay",
        xtitle="D_{CP}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixdecay],
        xformula="D_CP_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=120,
        plotcopier=pc,
        isDCP="dec",
        CMStext="",
      ),
      dict(
        name="D_int_decay",
        xtitle="D_{int}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixdecay],
        xformula="D_int_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([-1, .8, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=180,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_0minus_VBFdecay",
        xtitle="D_{0-}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=14,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_0hplus_VBFdecay",
        xtitle="D_{0h+}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=14,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1_VBFdecay",
        xtitle="D_{#Lambda1}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=14,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1Zg_VBFdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .8, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=14,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_CP_VBF",
        xtitle="D_{CP}^{VBF}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixVBF],
        xformula="D_CP_VBF",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=17,
        plotcopier=pc,
        isDCP="prod",
        CMStext="",
      ),
      dict(
        name="D_int_VBF",
        xtitle="D_{int}^{VBF}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixVBF],
        xformula="D_int_VBF",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([-1., 0, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=15,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_0minus_HadVHdecay",
        xtitle="D_{0-}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=8,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_0hplus_HadVHdecay",
        xtitle="D_{0h+}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=7,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1_HadVHdecay",
        xtitle="D_{#Lambda1}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=10,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_L1Zg_HadVHdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=10,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_CP_HadVH",
        xtitle="D_{CP}^{VH}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixVH],
        xformula="D_CP_HadVH",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=9,
        plotcopier=pc,
        isDCP="prod",
        CMStext="",
      ),
      dict(
        name="D_int_HadVH",
        xtitle="D_{int}^{VH}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixVH],
        xformula="D_int_HadVH",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([-1, -.6, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel="D_{bkg} > 0.2",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=13,
        plotcopier=pc,
        CMStext="",
      ),

      dict(
        name="D_bkg",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=untaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Untagged",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=270,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_bkg_VBFdecay",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_VBFdecay",
        cutformula=VBFtaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-2jet",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=16,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_bkg_HadVHdecay",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_HadVHdecay",
        cutformula=HadVHtaggedcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-hadronic",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=12,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_2jet_VBF",
        xtitle="max#left(#vec{D}_{2jet}^{VBF}#right)",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="max(D_2jet_0plus, max(D_2jet_0minus, max(D_2jet_a2, max(D_2jet_L1, D_2jet_L1Zg))))",
        cutformula=dijetenrichcut,
        binning=np.arange(0, 1+1/8., 1/8.),
        legendargs=(.2, .79, .6, .84),
        categorylabel="dijet events",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=13,
        plotcopier=pc,
        CMStext="",
        iscategorydiscriminant="VBF",
      ),
      dict(
        name="D_2jet_VH",
        xtitle="max#left(#vec{D}_{2jet}^{WH}, #vec{D}_{2jet}^{ZH}#right)",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="max(max(D_HadZH_0plus, max(D_HadZH_0minus, max(D_HadZH_a2, max(D_HadZH_L1, D_HadZH_L1Zg)))), max(D_HadWH_0plus, max(D_HadWH_0minus, max(D_HadWH_a2, max(D_HadWH_L1, D_HadWH_L1Zg)))))",
        cutformula=dijetenrichcut,
        binning=np.arange(0, 1+1/8., 1/8.),
        legendargs=(.2, .79, .6, .84),
        categorylabel="dijet events",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=20,
        plotcopier=pc,
        CMStext="",
        iscategorydiscriminant="VH",
      ),

      dict(
        name="ZZPt_boosted",
        xtitle="p_{T}^{4l}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="min(ZZPt, 650)",
        cutformula=boostedenrichcut,
        binning = np.array([120, 200., 300, 400, 500, 600, 700]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Boosted",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=13,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="ZZPt_VBF1j",
        xtitle="p_{T}^{4l}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="min(ZZPt, 170)",
        cutformula=VBF1jenrichcut,
        binning = np.array([0., 60, 120., 180]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-1jet",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=12,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="ZZPt_VHLep",
        xtitle="p_{T}^{4l}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="min(ZZPt, 350)",
        cutformula=LepVHenrichcut,
        binning = np.array([0, 100, 200., 300, 400]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-leptonic",
        Dbkglabel="D_{bkg} > 0.7",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=4,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_bkg_boosted",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=boostedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="Boosted",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=17,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_bkg_VBF1j",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=VBF1jcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VBF-1jet",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=19,
        plotcopier=pc,
        CMStext="",
      ),
      dict(
        name="D_bkg_VHLep",
        xtitle="D_{bkg}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=LepVHcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .79, .8, .84),
        categorylabel="VH-leptonic",
        Dbkglabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=11,
        plotcopier=pc,
        CMStext="",
      ),
    ]

    for kwargs in plotkwargses[:]:
      withbestfitkwargs = kwargs.copy()
      withbestfitkwargs["hypothesislines"] = kwargs["hypothesislines"] + [BestFit19009]
      withbestfitkwargs["saveasdir"] = os.path.join(kwargs["saveasdir"], "withbestfit")
      withbestfitkwargs["name"] += "_withbestfit"
      plotkwargses.append(withbestfitkwargs)

    for kwargs in plotkwargses[:]:
      preliminarykwargs = kwargs.copy()
      preliminarykwargs["CMStext"] = "Preliminary"
      preliminarykwargs["saveasdir"] = os.path.join(kwargs["saveasdir"], "preliminary")
      preliminarykwargs["name"] += "_preliminary"
      #plotkwargses.append(preliminarykwargs)

      workinprogresskwargs = kwargs.copy()
      workinprogresskwargs["CMStext"] = "Work in progress"
      workinprogresskwargs["saveasdir"] = os.path.join(kwargs["saveasdir"], "workinprogress")
      workinprogresskwargs["name"] += "_workinprogress"
      #plotkwargses.append(workinprogresskwargs)

    for kwargs in plotkwargses: print kwargs["name"]

    plots = [
      Plot(**kwargs) for kwargs in plotkwargses if filter(kwargs)
    ]

    for plot in plots:
      plot.makeplot()

if __name__ == "__main__":
  makeplots(**args.__dict__)
