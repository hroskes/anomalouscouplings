#!/usr/bin/env python

from functools import wraps
import os

print "importing root"
import ROOT

print "config"
from helperstuff import config
print "style"
from helperstuff import stylefunctions as style
print "enums"
from helperstuff.enums import EnumItem, Category, MultiEnum, MyEnum
print "utilities"
from helperstuff.utilities import cache, tfiles

print "PRL_loglinear"
from PRL_loglinear import yaxislabel
print "projections"
from projections import Projections
print "done"

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
    print self
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
def makearrow():
  color = ROOT.kOrange+7
  width = 2
  line = ROOT.TLine(0.5, 0, 0.5, 20)
  arrow = ROOT.TArrow(0.5, 20, 0.8, 20, .008, "|>")
  line.SetLineColor(color)
  line.SetLineWidth(width)
  arrow.SetLineColor(color)
  arrow.SetLineWidth(width)
  return line, arrow

def drawarrow():
  for _ in makearrow(): _.Draw()

def PRL_discriminants():
  rows = [
    [
      Discriminant("fa3", "fullrange", "Z", "allevents", config.productionforcombine, fulllegend=True, legendposition=(.15, .38, .75, .98)),
      Discriminant("fa2", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=22, legendposition=(.5,.7,.9,.85)),
      Discriminant("fL1", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=39, legendposition=(.53,.7,.93,.85)),
      Discriminant("fL1Zg", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=35, legendposition=(.6,.7,1,.85)),
    ],
    [
      Discriminant("fa3", "enrich", "X", "Untagged_with2015", config.productionforcombine, maximum=14, legendposition=(.5,.85,.9,1)),
      Discriminant("fa3", "enrich", "X", "VBFtagged", config.productionforcombine, maximum=4.8, legendposition=(.6,.85,1,1)),
      Discriminant("fa3", "enrich", "X", "VHHadrtagged", config.productionforcombine, maximum=3.7, legendposition=(.6,.85,1,1)),
      Discriminant("fa3", "enrich", "Y", "Untagged_with2015", config.productionforcombine, maximum=19, legendposition=(.56,.85,.9,1)),
    ],
  ]

  nrows = len(rows)
  ncolumns = {len(row) for row in rows}
  assert len(ncolumns) == 1; ncolumns = ncolumns.pop()

  c = ROOT.TCanvas("c", "c", 1600*ncolumns, 1600*nrows)
  style.applycanvasstyle(c)

  leftmargin = .024
  rightmargin = .02
  topmargin = .07
  bottommargin = 0
  c.SetLeftMargin(leftmargin)
  c.SetRightMargin(rightmargin)
  c.SetTopMargin(topmargin)
  c.SetBottomMargin(bottommargin)

  cachelist = []

  for iy, row in enumerate(rows):
    ymin = bottommargin + (1-topmargin-bottommargin) * (nrows-iy-1) / nrows
    ymax = ymin + (1-topmargin-bottommargin) / nrows
    for ix, plot in enumerate(row, start=1):
      xmin = leftmargin + (1-rightmargin-leftmargin) * (ix-1) / ncolumns
      xmax = xmin + (1-rightmargin-leftmargin) / ncolumns
      print xmin, xmax, ymin, ymax
      idx = iy*ncolumns + ix
      c.cd()
      pad = ROOT.TPad("pad{}".format(idx), "", xmin, ymin, xmax, ymax)
      pad.SetCanvas(c)
      pad.Draw()
      pad.cd()
      cachelist.append(pad)
      style.applycanvasstyle(pad)
      letter = "abcdefgh"[idx-1]
      pad.SetLeftMargin(.08)
      pad.SetRightMargin(.02)
      pad.SetTopMargin(0)
      pad.SetBottomMargin(.23)
      hstack, graph, legend = plot.hstack, plot.graph, plot.legend
      hstack.Draw("nostack")
      hstack.GetYaxis().SetTitle("")
      hstack.GetXaxis().SetLabelSize(.08)
      hstack.GetYaxis().SetLabelSize(.08)
      hstack.GetXaxis().SetTitleSize(.10)
      hstack.GetYaxis().SetTitleSize(.10)
      graph.Draw("P")
      legend.Draw()
      if ix == 1 and iy == 0: drawarrow()
      if iy == 1:
        style.subfig(letter, textsize=.08, x1=.16, x2=.2, y1=.92, y2=.98)
      else:
        style.subfig(letter, textsize=.08, x1=.86, x2=.9, y1=.92, y2=.98)


  c.cd()
  style.CMS("", lumi=None, lumitext="{:.1f} fb^{{-1}} (13 TeV)"
                                    .format(config.productionforcombine.dataluminosity+config.lumi2015),
                x1=0.04, x2=1.015, #???
                drawCMS=False, extratextsize=.045)
  style.CMS("", x1=0.022, x2=1.015, CMStextsize=.055)
  yaxislabel("Events / bin", textsize=.045).Draw()

  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "PRL.{}".format(ext)))

if __name__ == "__main__":
    PRL_discriminants()
