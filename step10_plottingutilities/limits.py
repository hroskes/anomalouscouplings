from collections import namedtuple
from helperstuff import config
from helperstuff.enums import Analysis
import os
import ROOT
import sys

Point = namedtuple("Point", "x y")

def findwhereyequals(y, p1, p2):
    m = (p2.y-p1.y) / (p2.x - p1.x)
    b =  p1.y - m*p1.x
    return (y-b)/m

def printlimits(analysis, foldername, **kwargs):
    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    plotname = "limit"
    for kw, kwarg in kwargs.iteritems():
        if kw == "plotname":
            plotname = kwarg
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))
    filename = os.path.join(config.plotsbasedir, "limits", foldername, plotname+".root")

    f = ROOT.TFile(filename)
    c = f.c1
    legend = c.GetListOfPrimitives()[2]

    thresholds = [1, 3.84]

    for entry in legend.GetListOfPrimitives():
        g, label = entry.GetObject(), entry.GetLabel()
        points = [Point(x, y) for x, y, n in zip(g.GetX(), g.GetY(), range(g.GetN()))]
        isabove = {threshold: points[0].y > threshold for threshold in thresholds}
        results = {threshold: ["[-1,"] if not isabove[threshold] else [] for threshold in thresholds}
        for point in points:
            for threshold in thresholds:
                if isabove[threshold] and point.y < threshold:
                    results[threshold].append("[{:.2f},".format(findwhereyequals(threshold, point, lastpoint)))
                    isabove[threshold]=False
                if not isabove[threshold] and point.y > threshold:
                    results[threshold][-1] += "{:.2f}]".format(findwhereyequals(threshold, point, lastpoint))
                    isabove[threshold]=True
            lastpoint = point

        if not isabove[threshold]:
            results[threshold][-1] += "1]"

        for threshold in thresholds:
            print label, threshold, " U ".join(results[threshold])

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
