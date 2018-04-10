#!/usr/bin/env python
from collections import namedtuple, OrderedDict
from helperstuff import config
from helperstuff.enums import Analysis, EnumItem, MyEnum
import os
import ROOT
import sys

Point = namedtuple("Point", "x y")

class PrintFormat(MyEnum):
    enumitems = (
                 EnumItem("latex"),
                 EnumItem("powerpoint", "ppt"),
                )

    @property
    def printformat(self):
        if self == "latex":
            return "${min:.2g}^{{{pluscl:+.2g}}}_{{{minuscl:+.2g}}}$ ${95%}$"
        if self == "ppt":
            return "{min:.2g}^({pluscl:+.2g})_({minuscl:+.2g}) {95%}"
        assert False

def findwhereyequals(y, p1, p2):
    m = (p2.y-p1.y) / (p2.x - p1.x)
    b =  p1.y - m*p1.x
    return (y-b)/m

def getlimits(filename, domirror=False):
    f = ROOT.TFile(filename)
    c = f.c1
    legend = c.GetListOfPrimitives()[2]

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
        results = {threshold: [[-1]] if not isabove[threshold] else [] for threshold in thresholds}
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
                results[threshold][-1].append(1)

        finalresult[label] = minimum.x, results[1.0], results[3.84]

    return finalresult

def printlimits(analysis, foldername, **kwargs):
    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    plotname = "limit"
    printformat = PrintFormat("latex")
    subdirectory = ""
    for kw, kwarg in kwargs.iteritems():
        if kw == "plotname":
            plotname = kwarg
            for ext in "png eps root pdf".split():
                plotname = plotname.replace("."+ext, "")
        elif kw == "printformat":
            printformat = PrintFormat(kwarg)
        elif kw == "subdirectory":
            subdirector = kwarg
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    filename = os.path.join(config.plotsbasedir, "limits", subdirectory, foldername, plotname+".root")
    allresults = getlimits(filename, domirror=(analysis=="fa3"))

    for label, (minimum, c68, c95) in allresults.iteritems():
        repmap = {}

        print label
        print

        repmap["min"] = minimum
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
                print "ranges: ", " \cup ".join("[{:.2g},{:.2g}]".format(range_[0],range_[1]) for range_ in c68)
                for range_ in c68:
                    if range_[0] < minimum < range_[1]:
                        break
                repmap["pluscl"] = range_[1] - minimum
                repmap["minuscl"] = range_[0] - minimum

        repmap["95%"] = " \cup ".join("[{:.2g},{:.2g}]".format(range_[0],range_[1]) for range_ in c95)

        print printformat.printformat.format(**repmap)

        print
        print

        return minimum, results

if __name__ == "__main__":
    analysis = Analysis(sys.argv[1])
    foldername = sys.argv[2]
    kwargs = {}
    for arg in sys.argv[3:]:
        kw, kwarg = arg.split("=")
        if kw in kwargs:
            raise TypeError("Duplicate kwarg {}!".format(kw))
        kwargs[kw] = kwarg
    printlimits(analysis, foldername, **kwargs)
