#!/usr/bin/env python
import argparse


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("analysis")
    p.add_argument("foldername")
    p.add_argument("--plotname")
    p.add_argument("--printformat")
    p.add_argument("--subdirectory")
    p.add_argument("--poi", default="CMS_zz4l_fai1")
    p.add_argument("--airatio", action="store_true")
    p.add_argument("--usemg", action="store_true")
    args = p.parse_args()

import math
import os
import sys
from collections import namedtuple, OrderedDict

import ROOT

from helperstuff import config
from helperstuff.enums import Analysis, EnumItem, MyEnum

class Point(object):
  def __init__(self, x, y):
    self.x = x
    self.y = y
  def __isub__(self, other):
    self.y -= other.y

class PrintFormat(MyEnum):
    enumitems = (
                 EnumItem("latex"),
                 EnumItem("powerpoint", "ppt"),
                )

    def printformat(self, minimum, pluscl, minuscl, str95):
        try:
          digits = max(
            min(-math.ceil(math.log(pluscl, 10)), -math.ceil(math.log(-minuscl, 10)))+2,
            -math.ceil(math.log(pluscl/2., 10))+1,
            -math.ceil(math.log(-minuscl/2., 10))+1,
            0,
          )
        except ValueError as e:
          print "digits value error", e
          digits = 2
        formatdict = {
            "minimum": minimum,
            "pluscl": pluscl,
            "minuscl": minuscl,
            "str95": str95,
        }
        percentdict = {
            "digits": digits
        }
        if float(("{minimum:.%(digits)df}" % percentdict).format(**formatdict)) == 0:
            percentdict["plus"] = ""
            formatdict["minimum"] = 0  #to avoid -0.0000
        else:
            percentdict["plus"] = ""

        if self == "latex":
            fmt = "${minimum:%(plus)s.%(digits)df}^{{{pluscl:+.%(digits)df}}}_{{{minuscl:+.%(digits)df}}}$ ${str95}$" % percentdict
        if self == "ppt":
            fmt = "${minimum:%(plus)s.%(digits)df}^({pluscl:+.%(digits)df})_({minuscl:+.%(digits)df}) {str95}" % percentdict
        try:
            return fmt.format(**formatdict)
        except:
            print fmt
            raise

def str95(*ranges):
  result = []
  for lo, hi in ranges:
    if lo < 0 < hi:
      digits = max(
        min(-math.ceil(math.log(hi, 10)), -math.ceil(math.log(-lo, 10)))+1,
        -math.ceil(math.log(hi, 10))+1,
        -math.ceil(math.log(-lo, 10))+1,
        0,
      )+1
    else:
      digits = 2
    result.append(("[{:.%(digits)df}, {:.%(digits)df}]" % {"digits": digits}).format(lo, hi))
  
  return " \cup ".join(result)

def findwhereyequals(y, p1, p2):
    m = (p2.y-p1.y) / (p2.x - p1.x)
    b =  p1.y - m*p1.x
    return (y-b)/m

def getgraphs(legend, mg, usemg):
  if usemg:
    for i, g in enumerate(mg.GetListOfGraphs()):
      yield g, str(i)
  else:
    for entry in legend.GetListOfPrimitives():
      yield entry.GetObject(), entry.GetLabel()

def getlimits(filename, poi, domirror=False, airatio=False, analysis=None, usemg=False):
    f = ROOT.TFile(filename)
    c = f.c1
    legend = c.GetListOfPrimitives()[2]
    mg = c.GetListOfPrimitives()[1]

    poimin, poimax = {
        "CMS_zz4l_fai1": (-float("inf"), float("inf")) if airatio else (-1, 1),
        "RV": (0, float("inf")),
        "RF": (0, float("inf")),
    }[poi]

    thresholds = [1, 3.84]
    NLL = {}

    finalresult = OrderedDict()

    if airatio:
        assert poi == "CMS_zz4l_fai1", poi
        assert analysis
        def transform(fa3):
            from helperstuff.combinehelpers import sigmaioversigma1
            if abs(fa3) == 1: return math.copysign(float("inf"), fa3)
            return math.copysign(math.sqrt(abs(fa3) / (1-abs(fa3)) / sigmaioversigma1(analysis, "ggH")), fa3)
        if analysis in ("fa3", "fa2"):
            pass
        elif analysis in ("fL1", "fL1Zg"):
            fa3transform = transform
            def transform(fL1):
                result = fa3transform(fL1)
                #result is g1prime2 / a1
                L1JHU = 10000 #GeV
                #set g1prime2 = L1JHU^2 / L1^2
                #then result is L1JHU^2 / (a1 L1^2)
                if result == 0: return float("inf")
                return math.copysign(math.sqrt(abs(L1JHU**2 / result)), result)
                #which is L1 sqrt(a1)
    else:
        def transform(fa3): return fa3

    for g, label in getgraphs(legend, mg, usemg):
        minimum = Point(float("nan"), float("infinity"))
        points = [Point(transform(x), y) for x, y, n in zip(g.GetX(), g.GetY(), range(g.GetN()))]
        if domirror and "Expected" in label:
            def findy(x):
                y = {yy for xx, yy in points if x==xx or abs(xx - x) == 0}   #x==xx to cover float("inf")
                assert len(y) == 1, (x, y)
                return y.pop()
            points = [Point(x, (y+findy(-x))/2 if any(xx==-x for xx, yy in points) else y) for x, y in points]
        points.sort()
        isabove = {threshold: points[0].y > threshold for threshold in thresholds}
        results = {threshold: [[poimin]] if not isabove[threshold] else [] for threshold in thresholds}
        for point in points:
            if point.y < minimum.y: minimum = point

        for point in sorted(points, key=lambda x: x is minimum):
          point -= minimum
        assert minimum.y == 0

        for point in points:
            for threshold in thresholds:
                if isabove[threshold] and point.y < threshold:
                    results[threshold].append([findwhereyequals(threshold, point, lastpoint)])
                    isabove[threshold]=False
                if not isabove[threshold] and point.y > threshold:
                    results[threshold][-1].append(findwhereyequals(threshold, point, lastpoint))
                    isabove[threshold]=True

            NLL[point.x] = point.y

            lastpoint = point

        for threshold in thresholds:
            if not isabove[threshold]:
                results[threshold][-1].append(poimax)

        finalresult[label] = minimum.x, results[1.0], results[3.84]

    return finalresult

def printlimits(analysis, foldername, **kwargs):
    if str(analysis) == "table":
      print r"\begin{tabular}{ccc}"
      print r"\hline"
      print r" Parameter                                   &  {Observed} & {Expected}    \\"
      print r"\hline"
      for analysis in Analysis.items(lambda x: x in ("fa3", "fa2", "fL1", "fL1Zg")):
        print analysis.latexfai,
        printlimits(analysis, foldername, fortable=True, **kwargs)
        print r"\\"
      print r"\hline"
      print r"\end{tabular}"
      return

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    plotname = "limit"
    printformat = PrintFormat("latex")
    subdirectory = ""
    poi = "CMS_zz4l_fai1"
    fortable = False
    airatio = False
    for kw, kwarg in kwargs.iteritems():
        if kw == "plotname":
            plotname = kwarg
            for ext in "png eps root pdf".split():
                plotname = plotname.replace("."+ext, "")
        elif kw == "printformat":
            printformat = PrintFormat(kwarg)
        elif kw == "subdirectory":
            subdirectory = kwarg
        elif kw == "poi":
            poi = kwarg
        elif kw == "fortable":
            fortable = kwarg
        elif kw == "airatio":
            airatio = kwarg
        elif kw == "usemg":
            usemg = kwarg
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    filename = os.path.join(config.plotsbasedir, "limits", subdirectory, foldername, plotname+".root")
    allresults = getlimits(filename, poi=poi, domirror=(analysis=="fa3"), airatio=airatio, analysis=analysis, usemg=usemg)
    if fortable:
        assert len(allresults.keys()) == 2 and allresults.keys()[0].startswith("Observed") and allresults.keys()[1].startswith("Expected"), allresults.keys()

    if not fortable and airatio:
        if analysis in ("fa3", "fa2"):
            print "Printing a_i / a_1"
            print
        elif analysis in ("fL1", "fL1Zg"):
            print "Printing Lambda_1^VV sqrt(a_1)"
            print

    for label, (minimum, c68, c95) in allresults.iteritems():
        repmap = {}

        if not fortable:
            print label
            print

        repmap["minimum"] = minimum
        if len(c68) == 1:
            repmap["pluscl"] = c68[0][1] - minimum
            repmap["minuscl"] = c68[0][0] - minimum
        else:
            if len(c68) == 2 and -c68[0][0] == c68[1][1] == 1:
                if minimum < c68[0][1]:
                    repmap["pluscl"] = c68[0][1] - minimum
                    repmap["minuscl"] = c68[1][0] - minimum - 2
                elif minimum > c68[1][0]:
                    repmap["pluscl"] = c68[0][1] - minimum + 2
                    repmap["minuscl"] = c68[1][0] - minimum
                else:
                    assert False
            else:
                print "Need something more complicated for 68% CL!"
                print "ranges: ", str95(*c68)
                for range_ in c68:
                    if range_[0] < minimum < range_[1]:
                        break
                repmap["pluscl"] = range_[1] - minimum
                repmap["minuscl"] = range_[0] - minimum

        repmap["str95"] = str95(*c95)

        if fortable:
            print "&", printformat.printformat(**repmap),
        else:
            print printformat.printformat(**repmap)
            print
            print

    return minimum, allresults

if __name__ == "__main__":
    for k, v in args.__dict__.copy().iteritems():
        if v is None: del args.__dict__[k]
    printlimits(**args.__dict__)
