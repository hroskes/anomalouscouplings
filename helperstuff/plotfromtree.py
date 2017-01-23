import re
import random

import ROOT

from helperstuff import config
from helperstuff import constants
from helperstuff.discriminants import discriminant
from helperstuff.enums import Category, Channel
import helperstuff.rootoverloads.histogramfloor
from helperstuff.samples import Sample, SampleBase

mandatory = object()

class Options(dict):
  def __getattr__(self, attr):
    return self[attr]
  def __setattr__(self, attr, value):
    self[attr] = value

def plotfromtree(**kwargs):
  o = Options({
    "productionmode": mandatory,
    "hypothesis":     mandatory,
    "weight":         None,

    "disc":           mandatory,
    "bins":           None,
    "min":            None,
    "max":            None,

    "disc2":          None,
    "bins2":           None,
    "min2":            None,
    "max2":            None,

    "transformation": None,  #should be a string with {disc} in it

    "enrich":         False,
    "masscut":        True,
    "normalizeto1":   False,

    "channel":        None,
    "category":       None,

    "color":          1,
    "hname":          None,
    "xaxislabel":     None,
  })
  for kw, kwarg in kwargs.iteritems():
    if kw in o:
      o[kw] = kwarg
    else:
      raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))
  for kw, kwarg in o.iteritems():
    if kwarg is mandatory:
      raise TypeError("kwarg {} is mandatory!".format(kw))

  discname, title, discbins, discmin, discmax = discriminant(o.disc)
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
    disc2name, disc2title, disc2bins, disc2min, disc2max = discriminant(o.disc2)
    if o.bins2 is None:
      o.bins2 = disc2bins
    if o.min2 is None:
      o.min2 = disc2min
    if o.max2 is None:
      o.max2 = disc2max

  if isinstance(o.weight, SampleBase):
    o.weight = o.weight.MC_weight

  for name, value in constants.__dict__.iteritems():
    try:
      value = "{:g}".format(value)
    except ValueError:
      continue

    if o.transformation is not None:
      o.transformation = re.sub(r"\b"+name+r"\b", value, o.transformation)
    if o.weight is not None:
      o.weight = re.sub(r"\b"+name+r"\b", value, o.weight)

  t = ROOT.TChain("candTree", "candTree")
  assert len(config.productionsforcombine) == 1
  sample = Sample(o.productionmode, o.hypothesis, config.productionsforcombine[0])

  t.Add(sample.withdiscriminantsfile())
  if o.hname is None:
    o.hname = "h{}".format(random.getrandbits(100))
  weightname = o.weight if o.weight is not None else sample.weightname()

  weightfactors = [
                   weightname,
                   "{}>-98".format(discname),
                  ]
  if o.category is not None:
    weightfactors.append(" || ".join("(category=={})".format(_) for _ in Category(o.category).idnumbers))
  if o.enrich:
    weightfactors.append("D_bkg>0.5")
  if o.channel is not None:
    weightfactors.append("Z1Flav*Z2Flav=={}".format(Channel(o.channel).ZZFlav))
  if o.masscut:
    weightfactors.append("ZZMass > {}".format(config.m4lmin))
    weightfactors.append("ZZMass < {}".format(config.m4lmax))
  if o.disc2 is not None:
    weightfactors.append("{}>-98".format(disc2name))

  wt = "*".join("("+_+")" for _ in weightfactors)

  if o.transformation is None:
    formula = discname
  else:
    formula = o.transformation.format(disc=discname)

  todraw = ""
  if o.disc2 is not None:
    todraw += "{}:".format(disc2name)
  todraw += "{}".format(formula)
  todraw += ">>{}".format(o.hname)
  if o.bins != 0:
    todraw += "({},{},{}".format(o.bins, o.min, o.max)
    if o.disc2 is not None:
      todraw += ",{},{},{}".format(o.bins2, o.min2, o.max2)
    todraw += ")"

  t.Draw(todraw, wt, "hist")
  h = getattr(ROOT, o.hname)
  h.GetXaxis().SetTitle(title)
  if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
    h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()+1) + h.GetBinContent(h.GetNbinsX()))
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
  if o.disc2 is not None:
    h.GetYaxis().SetTitle(disc2title)
  h.Floor()
  if o.normalizeto1:
    try:
      h.Scale(1/h.Integral())
    except ZeroDivisionError:
      pass
  h.SetLineColor(o.color)
  h.SetMarkerStyle(1)

  return h
