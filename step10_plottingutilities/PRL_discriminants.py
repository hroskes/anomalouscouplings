#!/usr/bin/env python

from functools import wraps
import os

import ROOT

from helperstuff import config
from helperstuff import stylefunctions as style
from helperstuff.enums import EnumItem, Category, MultiEnum, MyEnum
from helperstuff.utilities import cache, tfiles

from projections import Projections

class DiscriminantAxis(MyEnum):
  enumname = "discriminantaxis"
  enumitems = (
    EnumItem("X"),
    EnumItem("Y"),
    EnumItem("Z"),
  )

class ExtendedCategory(MyEnum):
  enumname = "extendedcategory"
  enumitems = Category.enumitems + (
    EnumItem("allevents"),
    EnumItem("Untagged_with2015"),
  )

def SetLegendPosition(legend, x1, y1, x2, y2):
  legend.SetX1NDC(x1)
  legend.SetX2NDC(x2)
  legend.SetY1NDC(y1)
  legend.SetY2NDC(y2)

class Discriminant(MultiEnum):
  enums = (Projections, DiscriminantAxis, ExtendedCategory)

  def __init__(self, *args, **kwargs):
    super(Discriminant, self).__init__(*args)
    hstack = self.hstack
    graph = self.graph
    legend = self.legend

    fulllegend = False
    for kw, kwarg in kwargs.iteritems():
      if kw == "maximum":
        hstack.SetMaximum(kwarg)
      elif kw == "legendposition":
        SetLegendPosition(self.legend, *kwarg)
      elif kw == "fulllegend":
        fulllegend = kwarg
      else:
        raise ValueError("Unknown kwarg {}={}!".format(kw, kwarg))

    if not fulllegend:
      lst = legend.GetListOfPrimitives()
      for entry in lst:
        if entry.GetLineStyle() == 2:
          if entry.GetLineColor() == 810: #total BSM
            label = entry.GetLabel()
            label = label.replace("total ", "#lower[0.5]{") + "}"
            entry.SetLabel(label)
          elif entry.GetLineColor() == 4:
            entry.SetLabel("")
          else:
            raise ValueError("color = {}???".format(entry.GetLineColor()))
        else:
          lst.Remove(entry)

    print legend.GetTextSize()
    legend.SetTextSize(.08)

  def check(self, *args):
    if self.extendedcategory == "allevents" and self.discriminantaxis != "Z":
      raise ValueError("Can only make a D_bkg plot with allevents!\n{}".format(args))
    super(Discriminant, self).check(*args)

  @property
  def discriminant(self):
    return self.projections.discriminants(self.category)["XYZ".index(str(self.discriminantaxis))]

  @property
  def category(self):
    if self.extendedcategory == "allevents":
      return Category("Untagged")
    else:
      if self.extendedcategory == "Untagged_with2015":
        return Category("Untagged")
      else:
        return Category(str(self.extendedcategory))

  @property
  def with2015(self):
    return self.extendedcategory in ("allevents", "Untagged_with2015")

  @property
  @cache
  def canvas(self):
    if self.extendedcategory == "allevents":
      folder = self.projections.saveasdir_Dbkgsum()
    else:
      folder = self.projections.saveasdir_niceplots(self.category, self.with2015)

    rootfile = tfiles[os.path.join(folder, self.discriminant.name+".root")]
    return rootfile.c1

  @cache
  def GetListOfPrimitives(self):
    return list(self.canvas.GetListOfPrimitives())

  @property
  @cache
  def hstack(self):
    hstack = self.GetListOfPrimitives()[1]
    if not isinstance(hstack, ROOT.THStack):
      raise ValueError("hstack has the wrong type:\n{}".format(hstack))
    return hstack.Clone()

  @property
  @cache
  def graph(self):
    graph = self.GetListOfPrimitives()[-1]
    if not isinstance(graph, ROOT.TGraph):
      raise ValueError("graph has the wrong type:\n{}".format(graph))
    graph.SetMarkerSize(7)
    return graph.Clone()

  @property
  @cache
  def legend(self):
    result = {_ for _ in self.GetListOfPrimitives() if isinstance(_, ROOT.TLegend)}
    assert len(result) == 1, result
    result = result.pop().Clone()
    return result

@cache
def TPad(*args, **kwargs):
  return ROOT.TPad(*args, **kwargs)

def PRL_discriminants():
  rows = [
    [
      Discriminant("fa3", "fullrange", "Z", "allevents", config.productionforcombine, fulllegend=True, legendposition=(.2, .35, .8, .95)),
      Discriminant("fa2", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=22, legendposition=(.5,.7,.9,.85)),
      Discriminant("fL1", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=39, legendposition=(.53,.7,.93,.85)),
      Discriminant("fL1Zg", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=35, legendposition=(.6,.7,1,.85)),
    ],
    [
      Discriminant("fa3", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=14, legendposition=(.5,.82,.9,.97)),
      Discriminant("fa3", "enrich", "X", "VBFtagged", config.productionforcombine, maximum=4.8, legendposition=(.6,.82,1,.97)),
      Discriminant("fa3", "enrich", "X", "VHHadrtagged", config.productionforcombine, maximum=3.7, legendposition=(.6,.82,1,.97)),
      Discriminant("fa3", "enrich", "Y", "Untagged_with2015", config.productionforcombine, maximum=19, legendposition=(.56,.82,.9,.97)),
    ],
  ]

  nrows = len(rows)
  ncolumns = {len(row) for row in rows}
  assert len(ncolumns) == 1; ncolumns = ncolumns.pop()

  c = ROOT.TCanvas("c", "c", 1600*ncolumns, 1600*nrows)
  style.applycanvasstyle(c)
  c.SetLeftMargin(.01)
  c.SetRightMargin(.1)
  c.SetTopMargin(.13)
  c.SetBottomMargin(0)
  c.Divide(ncolumns, nrows, 0, 0)

  for iy, row in enumerate(rows):
    for ix, plot in enumerate(row, start=1):
      pad = c.cd(iy*len(row) + ix)
      style.applycanvasstyle(pad)
      pad.SetLeftMargin(.08)
      pad.SetRightMargin(.02)
      pad.SetTopMargin(0)
      pad.SetBottomMargin(.18)
      hstack, graph, legend = plot.hstack, plot.graph, plot.legend
      hstack.Draw("nostack")
      hstack.GetYaxis().SetTitle("")
      hstack.GetXaxis().SetLabelSize(.08)
      hstack.GetYaxis().SetLabelSize(.08)
      hstack.GetXaxis().SetTitleSize(.08)
      hstack.GetYaxis().SetTitleSize(.08)
      graph.Draw("P")
      legend.Draw()

  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "PRL.{}".format(ext)))

if __name__ == "__main__":
    PRL_discriminants()
