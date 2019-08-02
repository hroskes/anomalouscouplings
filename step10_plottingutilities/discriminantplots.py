#!/usr/bin/env python

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
      if i % 10000 == 0 or i == nentries: print i, "/", nentries
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
  def __init__(self, name, trees, xformula, weightformula, cutformula, binning, scaleto, mirror=False, normalizationtrees=None, normalizationweightformula=None):
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
      self.__normalization = makenormalization(name=name+"_normalization", trees=normalizationtrees, weightformula=normalizationweightformula, scaleto=scaleto)

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
    print self.__histogram.GetName(), self.__histogram.Integral(),
    if self.__normalization is not None:
      self.__histogram.Scale(self.__normalization.scaleby)
    print self.__histogram.Integral()

  def SetBinContentError(self, x, value):
    assert not self.__finalized
    self.__histogram.SetBinContent(x, value.n)
    self.__histogram.SetBinError(x, value.s)

class Normalization(HistogramComponent):
  def __init__(self, name, trees, weightformula, scaleto):
    super(Normalization, self).__init__(
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
def makenormalization(**kwargs): return Normalization(**kwargs)

class Histogram(object):
  def __init__(self, name, trees, xformula, weightformulas, scaletos, cutformula, binning, linecolor, linestyle, linewidth, fillcolor, fillstyle, legendname, legendlpf, addonbottom, mirror, normalizationtrees=None, normalizationweightformulas=None):
    if normalizationweightformulas is None: normalizationweightformulas = weightformulas
    if normalizationtrees is None: normalizationtrees = trees
    if scaletos is None:
      scaletos = []
    else:
      assert len(scaletos) == len(weightformulas) == len(trees) == len(normalizationtrees) == len(normalizationweightformulas), (name, len(scaletos), len(weightformulas), len(trees), len(normalizationtrees), len(normalizationweightformulas))
    self.__components = [
      HistogramComponent(
        name=name+"_"+str(i), trees=componenttrees, xformula=xformula, weightformula=weightformula,
        cutformula=cutformula, binning=binning, mirror=mirror, scaleto=scaleto,
        normalizationtrees=componentnormalizationtrees, normalizationweightformula=normalizationweightformula,
      ) for i, (componenttrees, weightformula, scaleto, componentnormalizationtrees, normalizationweightformula) in enumerate(itertools.izip_longest(trees, weightformulas, scaletos, normalizationtrees, normalizationweightformulas))
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
    return legend.AddEntry(self.histogram, self.__legendname, self.__legendlpf)

@cache_keys(
  name=lambda name: None,
  trees=lambda trees: tuple(tuple(_) for _ in trees),
  addonbottom=tuple,
  weightformulas=tuple,
  binning=tuple,
  normalizationtrees=lambda normalizationtrees: tuple(tuple(_) for _ in normalizationtrees),
  normalizationweightformulas=tuple,
  scaletos=lambda scaletos: tuple(scaletos) if scaletos is not None else scaletos,
)
def makehistogram(**kwargs): return Histogram(**kwargs)

@cache
def gettree(*args, **kwargs):
  return Tree(*args, **kwargs)

@cache_file("trees_tmp.pkl")
def getfilenames(*productionmodesandhypotheses):
  print "Finding trees for", productionmodesandhypotheses
  return [
    [
      sample.withdiscriminantsfile()
      for sample in sorted(st, key=lambda x: x.withdiscriminantsfile())
    ]
    for production in config.productionsforcombine
    for otherargs in productionmodesandhypotheses
    for st in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg", *otherargs).reweightfrom()
  ]

def gettrees(*productionmodesandhypotheses):
  return [
    [
      gettree(
        filename,
        "candTree",
      ) for filename in lst
    ] for lst in getfilenames(*productionmodesandhypotheses)
  ]


def getweights(*productionmodesandhypotheses, **kwargs):
  scalebys = kwargs.pop("scalebys", [1 for _ in productionmodesandhypotheses])
  uselumi = kwargs.pop("uselumi", False)
  assert not kwargs
  return [
    "({}) * ({}) * ({})".format(weight, scaleby, production.dataluminosity if uselumi else 1)
    for production in config.productionsforcombine
    for otherargs, scaleby in itertools.izip_longest(productionmodesandhypotheses, scalebys)
    for weight in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg", *otherargs).weightname()
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
        "fa3fa2fL1fL1Zg",
        (
          productionmode
          if productionmode != "VH"
          else weight.split("p_Gen_")[1].split("_")[0]
        )
      )
      for ch in channels
      for ca in categories if ca != "Boosted"
    )
    for production in config.productionsforcombine
    for productionmode in productionmodes
    for weight in Template(production, "2e2mu", "Untagged", "fa3fa2fL1fL1Zg", productionmode, "0+", "Hff0+" if productionmode == "ttH" else None).weightname()
  ]

class HypothesisLine(object):
  def __init__(self, hypothesis, ffHlinecolor, ffHlinestyle, VVHlinecolor, VVHlinestyle, legendname):
    self.hypothesis = hypothesis
    self.ffHlinecolor = ffHlinecolor
    self.ffHlinestyle = ffHlinestyle
    self.VVHlinecolor = VVHlinecolor
    self.VVHlinestyle = VVHlinestyle
    self.legendname = legendname

  @property
  def ispure(self):
    return Hypothesis(self.hypothesis).ispure

  def getweights(self, *otherargs):
    scalebys = []
    for args in otherargs:
      productionmode = args[0]
      productionmodes = {
        "ggH": ("ggH", "bbH", "ttH"),
        "bbH": ("ggH", "bbH", "ttH"),
        "ttH": ("ggH", "bbH", "ttH"),
        "VBF": ("VBF", "VH"),
        "VH": ("VBF", "VH"),
        "ZH": ("VBF", "VH"),
        "ZH": ("VBF", "VH"),
      }[productionmode]
      numerator = sum(
        Template(production, "fa3fa2fL1fL1Zg", ca, ch, "0+", pm).gettemplate().Integral() * production.dataluminosity
        for ca in categories if ca != "Boosted"
        for ch in channels
        for pm in productionmodes
        for production in config.productionsforcombine
      )
      if Hypothesis(self.hypothesis).ispure:
        denominator = sum(
          Template(production, "fa3fa2fL1fL1Zg", ca, ch, self.hypothesis, pm).gettemplate().Integral() * production.dataluminosity
          for ca in categories if ca != "Boosted"
          for ch in channels
          for pm in productionmodes
          for production in config.productionsforcombine
        )
      else:
        s = ReweightingSample("ggH", self.hypothesis)
        if   s.g2 == s.g1prime2 == s.ghzgs1prime2 == 0:
          intletter = "i"
          ACname = "g4"
          AChypothesis = "0-"
        elif s.g4 == s.g1prime2 == s.ghzgs1prime2 == 0:
          intletter = "j"
          ACname = "g2"
          AChypothesis = "0h+"
        elif s.g2 == s.g4       == s.ghzgs1prime2 == 0:
          intletter = "k"
          ACname = "g1prime2"
          AChypothesis = "L1"

        gi = getattr(s, ACname)
        assert gi == getattr(ReweightingSample("VBF", self.hypothesis), ACname)
        assert gi != 0
        g1 = s.g1
        assert g1 == ReweightingSample("VBF", self.hypothesis).g1
        assert g1 != 0

        if productionmode in ("ggH", "bbH", "ttH"):
          inttype = "g11g{}1".format(intletter)
          denominator = sum(
            (
                 Template(production, "fa3fa2fL1fL1Zg", ca, ch, "0+",         pm).gettemplate().Integral() * g1*g1
            + IntTemplate(production, "fa3fa2fL1fL1Zg", ca, ch, inttype,      pm).gettemplate().Integral() * g1*gi
            +    Template(production, "fa3fa2fL1fL1Zg", ca, ch, AChypothesis, pm).gettemplate().Integral() * gi*gi
            ) * production.dataluminosity
            for ca in categories if ca != "Boosted"
            for ch in channels
            for pm in productionmodes
            for production in config.productionsforcombine
          )

        elif productionmode in ("VBF", "VH", "ZH", "WH"):
          def inttype(i): return "g1{}g{}{}".format(4-i, intletter, i)
          denominator = sum(
            (
                 Template(production, "fa3fa2fL1fL1Zg", ca, ch, "0+",         pm).gettemplate().Integral() * g1*g1*g1*g1
            + IntTemplate(production, "fa3fa2fL1fL1Zg", ca, ch, inttype(1),   pm).gettemplate().Integral() * g1*g1*g1*gi
            + IntTemplate(production, "fa3fa2fL1fL1Zg", ca, ch, inttype(2),   pm).gettemplate().Integral() * g1*g1*gi*gi
            + IntTemplate(production, "fa3fa2fL1fL1Zg", ca, ch, inttype(3),   pm).gettemplate().Integral() * g1*gi*gi*gi
            +    Template(production, "fa3fa2fL1fL1Zg", ca, ch, AChypothesis, pm).gettemplate().Integral() * gi*gi*gi*gi
            ) * production.dataluminosity
            for ca in categories if ca != "Boosted"
            for ch in channels
            for pm in productionmodes
            for production in config.productionsforcombine
          )
      scalebys.append(numerator/denominator)

    otherargs = [tuple(args) + (self.hypothesis,) for args in otherargs]
    assert len(otherargs) == len(scalebys)
    return getweights(*otherargs, scalebys=scalebys)

  def gettrees(self, *otherargs):
    hypothesis = self.hypothesis
    otherargs = (tuple(args) + (hypothesis,) for args in otherargs)
    return gettrees(*otherargs)

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
  def VVHtrees(self, isL1Zg=False):
    return self.gettrees(
      ("VBF",),
      ("VH" if not isL1Zg else "ZH",),
    )
  def VBFweights(self, isL1Zg=False):
    return self.getweights(
      ("VBF",),
    )
  def VBFtrees(self, isL1Zg=False):
    return self.gettrees(
      ("VBF",),
    )
  def notVBFweights(self, isL1Zg=False):
    return self.getweights(
      ("VH" if not isL1Zg else "ZH",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )
  def notVBFtrees(self, isL1Zg=False):
    return self.gettrees(
      ("VH" if not isL1Zg else "ZH",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )
  def VHweights(self, isL1Zg=False):
    return self.getweights(
      ("VH" if not isL1Zg else "ZH",),
    )
  def VHtrees(self, isL1Zg=False):
    return self.gettrees(
      ("VH" if not isL1Zg else "ZH",),
    )
  def notVHweights(self, isL1Zg=False):
    return self.getweights(
      ("VBF",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )
  def notVHtrees(self, isL1Zg=False):
    return self.gettrees(
      ("VBF",),
      ("ggH",),
      ("bbH",),
      ("ttH", "Hff0+"),
    )

class Plot(object):
  def __init__(self, **kwargs):
    preliminarykwargs = kwargs.copy()

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
    self.__saveasdir = saveasdir
    self.__plotcopier = plotcopier
    self.__iscategorydiscriminant = iscategorydiscriminant

    assert isDCP in (None, "prod", "dec")
    assert iscategorydiscriminant in (None, "VBF", "VH")

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

    for hypothesis in hypothesislines:
      isL1Zg = Hypothesis(hypothesis.hypothesis) == "L1Zg"

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
        legendname="VBF+VH "+hypothesis.legendname,
        legendlpf="f",
        addonbottom=[ZZhistogram],
        mirror=isDCP is not None and hypothesis.ispure,
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
        legendname="Total "+hypothesis.legendname,
        legendlpf="l",
        mirror=isDCP == "prod" or (isDCP == "dec" and hypothesis.ispure)
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
      else:
        assert False

      VVH = makehistogram(**VVHkwargs)
      ffHkwargs["addonbottom"] = [VVH]
      ffH = makehistogram(**ffHkwargs)

      histograms += [ffH, VVH]

    histograms += [ZZhistogram, ZXhistogram]

    self.__CMStext = CMStext

    if self.__CMStext != "Preliminary":
      preliminarykwargs["CMStext"] = "Preliminary"
      preliminarykwargs["saveasdir"] = os.path.join(saveasdir, "preliminary")
      preliminarykwargs["name"] += "_preliminary"
      self.__preliminary = type(self)(**preliminarykwargs)



  def makeplot(self):
    c = self.__plotcopier.TCanvas("c_"+self.__name, "",  8, 30, 800, 800)
    style.applycanvasstyle(c)

    for h in self.histograms:
      h.makefinalhistogram()

    hstack = ROOT.THStack(self.__name, "")
    l = ROOT.TLegend(*self.__legendargs)
    style.applylegendstyle(l)
    l.SetNColumns(self.__legendcolumns)

    for h in self.histograms:
      hstack.Add(h.histogram)
      h.addtolegend(l)

    hstack.Draw("hist nostack")
    style.applyaxesstyle(hstack)
    hstack.GetXaxis().SetTitle(self.__xtitle)
    hstack.GetYaxis().SetTitle(self.__ytitle)
    hstack.SetMaximum(self.__ymax)
    l.Draw()

    lumi = sum(production.dataluminosity for production in config.productionsforcombine)

    style.CMS(self.__CMStext, lumi)

    if self.__categorylabel is not None:
      categorytext = ROOT.TPaveText(.71, .5, .91, .6, "brNDC")
      categorytext.SetBorderSize(0)
      categorytext.SetTextAlign(12)
      categorytext.SetTextSize(0.045)
      categorytext.SetFillStyle(0)
      categorytext.SetTextFont(42)
      categorytext.AddText(0.01,0.01,self.__categorylabel)
      categorytext.Draw()

    if self.__iscategorydiscriminant:
      line = ROOT.TLine(0.5, 0, 0.5, self.__ymax * .4)
      line.SetLineWidth(4)
      line.SetLineColor(ROOT.kViolet-1)
      line.SetLineStyle(9)
      line.Draw()

    mkdir_p(self.__saveasdir)
    for ext in "png pdf root C".split():
      c.SaveAs(os.path.join(self.__saveasdir, self.__name+"."+ext))

    if self.__CMStext != "Preliminary":
      self.__preliminary.makeplot()

categoryname = "category_0P_or_0M_or_a2_or_L1_or_L1Zg"

dijetcut = masscut + " && D_0minus_VBFdecay > -998"
untaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("Untagged").idnumbers) + ")"
VBFtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VBFtagged").idnumbers) + ")"
HadVHtaggedcut = masscut + " && (" + " || ".join("{} == {}".format(categoryname, c) for c in Category("VHHadrtagged").idnumbers) + ")"

dijetenrichcut = dijetcut + " && D_bkg > 0.7"
untaggedenrichcut = untaggedcut + " && D_bkg > 0.7"
VBFtaggedenrichcut = VBFtaggedcut + " && D_bkg_VBFdecay > 0.2"
HadVHtaggedenrichcut = HadVHtaggedcut + " && D_bkg_HadVHdecay > 0.2"

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


def makeplots():
  with PlotCopier() as pc:
    plots = [
      Plot(
        name="D_0minus_decay",
        xtitle="D_{0-}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_0hplus_decay",
        xtitle="D_{0h+}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixdecay],
        xformula="D_0hplus_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .5, .7, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_L1_decay",
        xtitle="D_{#Lambda1}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [L1mixdecay],
        xformula="D_L1_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .55, .8, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_L1Zg_decay",
        xtitle="D_{#Lambda1}^{Z#gamma,dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([0, .4, .55, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_CP_decay",
        xtitle="D_{CP}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixdecay],
        xformula="D_CP_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        isDCP="dec",
        CMStext="Supplementary",
      ),
      Plot(
        name="D_int_decay",
        xtitle="D_{int}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixdecay],
        xformula="D_int_decay",
        cutformula=untaggedenrichcut,
        binning=np.array([-1, .8, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_0minus_VBFdecay",
        xtitle="D_{0-}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_0hplus_VBFdecay",
        xtitle="D_{0h+}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_L1_VBFdecay",
        xtitle="D_{#Lambda1}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_L1Zg_VBFdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_VBFdecay",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([0, .1, .8, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_CP_VBF",
        xtitle="D_{CP}^{VBF}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixVBF],
        xformula="D_CP_VBF",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        isDCP="prod",
        CMStext="",
      ),
      Plot(
        name="D_int_VBF",
        xtitle="D_{int}^{VBF}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixVBF],
        xformula="D_int_VBF",
        cutformula=VBFtaggedenrichcut,
        binning=np.array([-1., 0, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_0minus_HadVHdecay",
        xtitle="D_{0-}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0minus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_0hplus_HadVHdecay",
        xtitle="D_{0h+}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_0hplus_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_L1_HadVHdecay",
        xtitle="D_{#Lambda1}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, 1./3, 2./3, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_L1Zg_HadVHdecay",
        xtitle="D_{#Lambda1}^{Z#gamma,VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_L1Zg_HadVHdecay",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([0, .1, .9, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),
      Plot(
        name="D_CP_HadVH",
        xtitle="D_{CP}^{VH}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a3mixVH],
        xformula="D_CP_HadVH",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([-1, 0., 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        isDCP="prod",
        CMStext="Supplementary",
      ),
      Plot(
        name="D_int_HadVH",
        xtitle="D_{int}^{VH}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines + [a2mixVH],
        xformula="D_int_HadVH",
        cutformula=HadVHtaggedenrichcut,
        binning=np.array([-1, -.6, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="Supplementary",
      ),

      Plot(
        name="D_bkg",
        xtitle="D_{bkg}^{dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg",
        cutformula=untaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="Untagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=500,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_bkg_VBFdecay",
        xtitle="D_{bkg}^{VBF+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_VBFdecay",
        cutformula=VBFtaggedcut,
        binning=np.array([0, .2, .7, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VBF-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_bkg_HadVHdecay",
        xtitle="D_{bkg}^{VH+dec}",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="D_bkg_HadVHdecay",
        cutformula=HadVHtaggedcut,
        binning=np.array([0, .2, .8, 1]),
        legendargs=(.2, .5, .9, .9),
        categorylabel="VH-tagged",
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=50,
        plotcopier=pc,
        CMStext="",
      ),
      Plot(
        name="D_2jet_VBF",
        xtitle="max#left(D_{2jet}^{VBF}, #vec{D}_{2jet}^{VBF, BSM}#right)",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="max(D_2jet_0plus, max(D_2jet_0minus, max(D_2jet_a2, max(D_2jet_L1, D_2jet_L1Zg))))",
        cutformula=dijetenrichcut,
        binning=np.arange(0, 1+1/8., 1/8.),
        legendargs=(.2, .5, .9, .9),
        categorylabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=30,
        plotcopier=pc,
        CMStext="",
        iscategorydiscriminant="VBF",
      ),
      Plot(
        name="D_2jet_VH",
        xtitle="max#left(D_{2jet}^{WH}, #vec{D}_{2jet}^{WH, BSM}, D_{2jet}^{ZH}, #vec{D}_{2jet}^{ZH, BSM}#right)",
        ytitle="Events / bin",
        hypothesislines=purehypothesislines,
        xformula="max(max(D_HadZH_0plus, max(D_HadZH_0minus, max(D_HadZH_a2, max(D_HadZH_L1, D_HadZH_L1Zg)))), max(D_HadWH_0plus, max(D_HadWH_0minus, max(D_HadWH_a2, max(D_HadWH_L1, D_HadWH_L1Zg)))))",
        cutformula=dijetenrichcut,
        binning=np.arange(0, 1+1/8., 1/8.),
        legendargs=(.2, .5, .9, .9),
        categorylabel=None,
        legendcolumns=2,
        saveasdir=os.path.join(config.plotsbasedir, "templateprojections", "niceplots"),
        ymax=30,
        plotcopier=pc,
        CMStext="",
        iscategorydiscriminant="VH",
      ),
    ]

    for plot in plots:
      plot.makeplot()

if __name__ == "__main__":
  makeplots()