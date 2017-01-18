from helperstuff import config
from helperstuff import constants
from helperstuff.discriminants import discriminant
from helperstuff.enums import Category, Channel
import helperstuff.rootoverloads.histogramfloor
from helperstuff.samples import Sample, SampleBase
import re
import ROOT

class Options(dict):
  def __getattr__(self, attr):
    return self[attr]
  def __setattr__(self, attr, value):
    self[attr] = value

def plotfromtree(**kwargs):
  o = Options({
    "productionmode": None,
    "hypothesis":     None,
    "weight":         None,

    "disc":           None,
    "bins":           None,
    "min":            None,
    "max":            None,

    "enrich":         False,
    "masscut":        True,
    "normalizeto1":   False,

    "channel":        None,
    "category":       None,

    "color":          1,
    "hname":          None,
  })
  for kw, kwarg in kwargs.iteritems():
    if kw in o:
      o[kw] = kwarg
    else:
      raise ValueError("Unknown kwarg {}={}".format(kw, kwarg))

  discname, title, discbins, discmin, discmax = discriminant(o.disc)
  if o.bins is None:
    o.bins = discbins
  if o.min is None:
    o.min = discmin
  if o.max is None:
    o.max = discmax

  if isinstance(o.weight, SampleBase):
    o.weight = o.weight.MC_weight

  for name, value in constants.__dict__.iteritems():
    try:
      value = "{:g}".format(value)
    except ValueError:
      continue

    o.disc = re.sub(r"\b"+name+r"\b", value, o.disc)
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
                   "{}>-998".format(discname),
                  ]
  if o.category is not None:
    weightfactors.append(" || ".join("(category=={})".format(_) for _ in Category(o.category).idnumbers))
  if o.enrich:
    weightfactors.append("D_bkg_0plus>0.5")
  if o.channel is not None:
    weightfactors.append("Z1Flav*Z2Flav=={}".format(Channel(o.channel).ZZFlav))
  if o.masscut:
    weightfactors.append("ZZMass > {}".format(config.m4lmin))
    weightfactors.append("ZZMass < {}".format(config.m4lmax))

  wt = "*".join("("+_+")" for _ in weightfactors)

  t.Draw("{}>>{}({},{},{})".format(discname, o.hname, o.bins, o.min, o.max), wt, "hist")
  h = getattr(ROOT, o.hname)
  h.GetXaxis().SetTitle(title)
  if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
    h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()+1) + h.GetBinContent(h.GetNbinsX()))
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
  h.Floor()
  if o.normalizeto1:
    try:
      h.Scale(1/h.Integral())
    except ZeroDivisionError:
      pass
  h.SetLineColor(o.color)
  h.SetMarkerStyle(1)

  return h
