from collections import namedtuple
import re
import random

import ROOT

from helperstuff import config
from helperstuff import constants
from helperstuff.discriminants import discriminant
from helperstuff.enums import Analysis, Category, Channel
import helperstuff.rootoverloads.histogramfloor
from helperstuff.samples import ReweightingSample, Sample, SampleBase
from helperstuff.utilities import TFile

mandatory = object()

class Options(dict):
  def __getattr__(self, attr):
    return self[attr]
  def __setattr__(self, attr, value):
    self[attr] = value

def plotfromtree(**kwargs):
  o = Options({
    "reweightfrom":   mandatory,
    "reweightto":     None,
    "production":     None,

    "disc":           mandatory,
    "bins":           None,
    "min":            None,
    "max":            None,

    "disc2":          None,
    "bins2":          None,
    "min2":           None,
    "max2":           None,

    "transformation": None,  #should be a string with {disc} in it

    "enrich":         False,
    "masscut":        True,
    "normalizeto1":   False,

    "channel":        None,
    "category":       None,
    "categorization": None,  #for which categorization to use, or set the next argument
    "analysis":       None,  #for which categorization to use, or set the previous argument
    "cut":            None,

    "color":          1,
    "linestyle":      1,
    "hname":          None,
    "xaxislabel":     None,
  })
  for kw, kwarg in kwargs.iteritems():
    if kw in o:
      o[kw] = kwarg
    else:
      raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

  if not isinstance(o.reweightfrom, Sample):
    if o.production is not None: o.reweightfrom = Sample(o.reweightfrom, o.production)
    elif len(config.productionsforcombine) == 1: o.reweightfrom = Sample(o.reweightfrom, config.productionsforcombine[0])
    else: raise ValueError("More than one productionforcombine in config, so you have to provide it")

  if o.reweightfrom.production != o.production is not None:
    raise ValueError("production and reweightfrom.production are inconsistent\n{}, {}".format(o.production, o.reweightfrom.production))

  for kw, kwarg in o.iteritems():
    if kwarg is mandatory:
      raise TypeError("kwarg {} is mandatory!".format(kw))

  discname, title, discbins, discmin, discmax, _, discformula = discriminant(o.disc)
  if o.transformation is not None:
    title = "f(" + title + ")"
  if o.xaxislabel is not None:
    title = o.xaxislabel
  if o.bins is None:
    o.bins = discbins
  if o.min is None:
    o.min = discmin
  if o.max is None:
    o.max = discmax

  if o.disc2 is not None:
    disc2name, disc2title, disc2bins, disc2min, disc2max, _, disc2formula = discriminant(o.disc2)
    if o.bins2 is None:
      o.bins2 = disc2bins
    if o.min2 is None:
      o.min2 = disc2min
    if o.max2 is None:
      o.max2 = disc2max

  if o.reweightto is None:
    o.reweightto = o.reweightfrom
  if isinstance(o.reweightto, SampleBase):
    o.reweightto = o.reweightto.MC_weight

  for name, value in constants.__dict__.iteritems():
    try:
      value = "{:g}".format(value)
    except ValueError:
      continue

    if o.transformation is not None:
      o.transformation = re.sub(r"\b"+name+r"\b", value, o.transformation)
    if o.reweightto is not None:
      o.reweightto = re.sub(r"\b"+name+r"\b", value, o.reweightto)

  t = ROOT.TChain("candTree", "candTree")

  t.Add(o.reweightfrom.withdiscriminantsfile())
  if o.hname is None:
    o.hname = "h{}".format(random.getrandbits(100))
  weightname = o.reweightto if o.reweightto is not None else o.reweightfrom.weightname()

  weightfactors = [
                   weightname,
                   "{}>-98".format(discformula),
                  ]
  if o.category is not None:
    o.category = Category(o.category)
    if o.analysis is o.categorization is None: raise TypeError("analysis or categorization is mandatory if category is provided!")
    if o.analysis is not None is not o.categorization: raise TypeError("Can't provide both analysis and categorization!")
    if o.analysis is not None: categoryname = Analysis(o.analysis).categoryname
    else:                      categoryname = o.categorization
    weightfactors.append(" || ".join("(category_{}=={})".format(categoryname, _) for _ in Category(o.category).idnumbers))
  if o.enrich:
    weightfactors.append("D_bkg>0.5")
  if o.channel is not None:
    weightfactors.append("Z1Flav*Z2Flav=={}".format(Channel(o.channel).ZZFlav))
  if o.masscut:
    weightfactors.append("ZZMass > {}".format(config.m4lmin))
    weightfactors.append("ZZMass < {}".format(config.m4lmax))
  if o.disc2 is not None:
    weightfactors.append("{}>-98".format(disc2formula))
  if o.cut is not None:
    weightfactors.append(o.cut)

  wt = "*".join("("+_+")" for _ in weightfactors)

  formula = "min(max({}, {}), {})".format(discformula, o.min, o.max - (o.max-o.min)/100)
  if o.transformation is not None:
    formula = o.transformation.format(formula)

  todraw = ""
  if o.disc2 is not None:
    todraw += "min(max({}, {}), {}):".format(disc2formula, o.min2, o.max2 - (o.max2-o.min2)/100)
  todraw += "{}".format(formula)
  todraw += ">>{}".format(o.hname)
  if o.bins != 0:
    todraw += "({},{},{}".format(o.bins, o.min, o.max)
    if o.disc2 is not None:
      todraw += ",{},{},{}".format(o.bins2, o.min2, o.max2)
    todraw += ")"

  t.Draw(todraw, wt, "hist")
  try:
    h = getattr(ROOT, o.hname)
  except:
    print
    print "using file:"
    print o.reweightfrom.withdiscriminantsfile()
    print
    raise
  h.SetTitle(title)
  h.GetXaxis().SetTitle(title)
  if o.disc2 is not None:
    h.GetYaxis().SetTitle(disc2title)
  h.Floor()
  if o.normalizeto1:
    try:
      h.Scale(1/h.Integral())
    except ZeroDivisionError:
      pass
  h.SetLineColor(o.color)
  h.SetLineStyle(o.linestyle)
  h.SetMarkerColor(o.color)
  h.SetMarkerStyle(1)

  return h

class Line(namedtuple("Line", "sample title color reweightfrom")):
  """useful namedtuple, no special interaction with anything else here"""
  def __new__(cls, sample, title, color, reweightfrom=None, bkpreweightfrom=None, production=None):
    if reweightfrom is None: reweightfrom = sample
    if not isinstance(reweightfrom, Sample):
      if production is None:
        assert len(config.productionsforcombine) == 1
        production = config.productionsforcombine[0]
      try:
        reweightfrom = Sample(reweightfrom, production)
        with TFile(reweightfrom.withdiscriminantsfile()) as f:
          if f.candTree.GetEntries() == 0: raise ValueError("No events in the tree: "+reweightfrom.withdiscriminantsfile())
      except ValueError:
        if bkpreweightfrom is None: raise
        reweightfrom = Sample(bkpreweightfrom, production)
      with TFile(reweightfrom.withdiscriminantsfile()) as f:
        if f.candTree.GetEntries() == 0: raise ValueError("No events in the tree: "+reweightfrom.withdiscriminantsfile())
    return super(Line, cls).__new__(cls, sample, title, color, reweightfrom)
