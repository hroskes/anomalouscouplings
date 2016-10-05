from collections import namedtuple
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
            return "${min:.2g}^{{{pluscl:+.2g}}}_{{{minuscl:+.2g}}}$ $ {95%}$"
        if self == "ppt":
            return "{min:.2g}^({pluscl:+.2g})_({minuscl:+.2g}) {95%}"
        assert False

def findwhereyequals(y, p1, p2):
    m = (p2.y-p1.y) / (p2.x - p1.x)
    b =  p1.y - m*p1.x
    return (y-b)/m

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

    f = ROOT.TFile(filename)
    c = f.c1
    legend = c.GetListOfPrimitives()[2]

    thresholds = [1, 3.84]
    NLL = {}

    for entry in legend.GetListOfPrimitives():
        minimum = Point(float("nan"), float("infinity"))
        g, label = entry.GetObject(), entry.GetLabel()
        points = [Point(x, y) for x, y, n in zip(g.GetX(), g.GetY(), range(g.GetN()))]
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

        repmap = {}

        print label
        print

        repmap["min"] = minimum.x
        if len(results[1.0]) == 1:
            repmap["pluscl"] = results[1.0][0][1] - minimum.x
            repmap["minuscl"] = results[1.0][0][0] - minimum.x
        else:
            if len(results[1.0]) == 2 and -results[1.0][0][0] == results[1.0][1][1] == 1:
                if minimum.x < results[1.0][0][1]:
                    repmap["pluscl"] = results[1.0][0][1] - minimum.x
                    repmap["minuscl"] = results[1.0][1][0] - minimum.x - 2
                elif minimum.x > results[1.0][1][0]:
                    repmap["pluscl"] = results[1.0][0][1] - minimum.x + 2
                    repmap["minuscl"] = results[1.0][1][0] - minimum.x
                else:
                    assert False
            else:
                print "Need something more complicated for 68% CL!"
                print "ranges: ", " \cup ".join("[{:.2g},{:.2g}]".format(range_[0],range_[1]) for range_ in results[1])
                for range_ in results[1.0]:
                    if range_[0] < minimum.x < range_[1]:
                        break
                repmap["pluscl"] = range_[1] - minimum.x
                repmap["minuscl"] = range_[0] - minimum.x

        repmap["95%"] = " \cup ".join("[{:.2g},{:.2g}]".format(range_[0],range_[1]) for range_ in results[3.84])

        print printformat.printformat.format(**repmap)

        if 0 in NLL and 1 in NLL:
            if NLL[1] >= NLL[0]:
                prob = ROOT.TMath.Prob(NLL[1]-NLL[0], 1)/2
            else:
                prob = 1 - ROOT.TMath.Prob(NLL[0]-NLL[1], 1)/2

            print "Probability for pure BSM vs. pure SM: {:.2g}{}".format(prob*100, "%" if printformat == "ppt" else r"\%")
        print
        print

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
