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
    args = p.parse_args()

import math
import os
import sys
from collections import namedtuple, OrderedDict

import ROOT

from helperstuff import config
from helperstuff.enums import Analysis, EnumItem, MyEnum

Point = namedtuple("Point", "x y")

class PrintFormat(MyEnum):
    enumitems = (
                 EnumItem("latex"),
                 EnumItem("powerpoint", "ppt"),
                )

    def printformat(self, minimum, pluscl, minuscl, str95):
        digits = max(
          min(-math.ceil(math.log(pluscl, 10)), -math.ceil(math.log(-minuscl, 10)))+2,
          -math.ceil(math.log(pluscl/2, 10))+1,
          -math.ceil(math.log(-minuscl/2, 10))+1,
        )
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
        min(-math.ceil(math.log(hi, 10)), -math.ceil(math.log(-lo, 10)))+2,
        -math.ceil(math.log(hi/2, 10))+1,
        -math.ceil(math.log(-lo/2, 10))+1,
      )
    else:
      digits = 2
    result.append(("[{:.%(digits)df}, {:.%(digits)df}]" % {"digits": digits}).format(lo, hi))
  
  return " \cup ".join(result)

def findwhereyequals(y, p1, p2):
    m = (p2.y-p1.y) / (p2.x - p1.x)
    b =  p1.y - m*p1.x
    return (y-b)/m

def getlimits(filename, poi, domirror=False):
    f = ROOT.TFile(filename)
    c = f.c1
    legend = c.GetListOfPrimitives()[2]

    poimin, poimax = {
        "CMS_zz4l_fai1": (-1, 1),
        "RV": (0, float("inf")),
        "RF": (0, float("inf")),
    }[poi]

    thresholds = [1, 3.84]
    NLL = {}

    finalresult = OrderedDict()

    for entry in legend.GetListOfPrimitives():
        minimum = Point(float("nan"), float("infinity"))
        g, label = entry.GetObject(), entry.GetLabel()
        points = [Point(x, y) for x, y, n in zip(g.GetX(), g.GetY(), range(g.GetN()))]
        if domirror and "Expected" in label:
            def findy(x):
                y = {yy for xx, yy in points if abs(xx - x) == 0}
                assert len(y) == 1, (x, y)
                return y.pop()
            points = [Point(x, (y+findy(-x))/2 if any(xx==-x for xx, yy in points) else y) for x, y in points]
        isabove = {threshold: points[0].y > threshold for threshold in thresholds}
        results = {threshold: [[poimin]] if not isabove[threshold] else [] for threshold in thresholds}
        for point in points:
            for threshold in thresholds:
                if isabove[threshold] and point.y < threshold:
                    results[threshold].append([findwhereyequals(threshold, point, lastpoint)])
                    isabove[threshold]=False
                if not isabove[threshold] and point.y > threshold:
                    results[threshold][-1].append(findwhereyequals(threshold, point, lastpoint))
                    isabove[threshold]=True

            if point.y < minimum.y: minimum = point

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
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    filename = os.path.join(config.plotsbasedir, "limits", subdirectory, foldername, plotname+".root")
    allresults = getlimits(filename, poi=poi, domirror=(analysis=="fa3"))
    if fortable:
        assert len(allresults.keys()) == 2 and allresults.keys()[0].startswith("Observed") and allresults.keys()[1].startswith("Expected"), allresults.keys()

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
