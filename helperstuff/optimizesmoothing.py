from collections import namedtuple
from discriminants import discriminant
from filemanager import tfiles
import itertools
from math import sqrt
import numpy
import ROOT
from templates import Template

Bin = namedtuple("Bin", "low hi contents error")
Derivative = namedtuple("Derivative", "center difference differror")

#http://stackoverflow.com/a/5434936/5228524
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def sign(x):
    return cmp(x, 0)

def cache(function):
    cachename = "__cache_{}".format(function.__name__)
    def newfunction(self):
        try:
            return getattr(self, cachename)
        except AttributeError:
            setattr(self, cachename, function(self))
            return newfunction(self)
    newfunction.__name__ = function.__name__
    return newfunction

def generatortolist(function):
    def newfunction(*args, **kwargs):
        return list(function(*args, **kwargs))
    newfunction.__name__ = function.__name__
    return newfunction

class Range(namedtuple("Range", "low hi lowval hival lowerr hierr")):
    def __init__(self, *args, **kwargs):
        super(Range, self).__init__(*args, **kwargs)
        if self.low > self.hi:
            raise ValueError("Can't have a range that goes from {} to {}!\nMaybe try switching the order?".format(self.low, self.hi))

    @property
    @cache
    def significance(self):
        """
        This function is supposed to return a number
        The bigger this number is, the greater the probability
          that this range is a real thing rather than a statistical fluctuation
        """
        return (self.hi - self.low) * abs(self.hival - self.lowval) / sqrt(self.lowerr**2 + self.hierr**2)

    def __contains__(self, other):
        return self.low <= other.low < other.hi <= self.hi

    def overlap(self, other):
        """
        http://stackoverflow.com/a/2953979/5228524
        """
        return max(0, min(self.hi, other.hi) - max(self.low, other.low))

    def samedirection(self, other):
        return sign(self.hival - self.lowval) == sign(other.hival - other.lowval)

    @property
    def length(self):
        return self.hi - self.low

class HistInfo(object):
    def __init__(self, h):
        self.h = h

    @property
    @cache
    def bins(self):
        return [Bin(self.GetBinLowEdge(i), self.GetBinLowEdge(i+1), self.GetBinContent(i), self.GetBinError(i))
                       for i in range(1, self.GetNbinsX()+1)]

    @property
    @cache
    def derivative(self):
        result = []
        for bin1, bin2 in pairwise(self.bins):
            assert bin1.hi == bin2.low
            result.append(Derivative(bin1.hi, bin2.contents - bin1.contents, sqrt(bin2.error**2+bin1.error**2)))
        return result

    @property
    @cache
    def increasingordecreasingranges(self):
        ranges = []
        direction = 0
        rangebegin = 1
        nextrangebegin = rangeend = nextrangebeginval = rangeendval = None
        for i, (middle, difference, differror) in enumerate(self.derivative, start=2):
            if direction == 0:
                direction = sign(difference)
                if rangebegin is None:
                    rangebegin = i
            if sign(difference) != direction:
                if rangeend is None:
                    rangeend = i-1
                if sign(difference) == 0:
                    nextrangebegin = i-1
                if sign(difference) != 0:
                    if nextrangebegin is None:
                        nextrangebegin = i-1
                    direction = sign(difference)
                    ranges.append(self.getrange(rangebegin, rangeend))
                    rangebegin, rangeend, nextrangebegin = nextrangebegin, None, None
            else:
                nextrangebegin = rangeend = None
        rangeend = self.GetNbinsX()
        ranges.append(self.getrange(rangebegin, rangeend))
        return ranges

    def getrange(self, begin, end):
        return Range(
                     self.GetBinLowEdge(begin) + self.binwidth/2, self.GetBinLowEdge(end) + self.binwidth/2,
                     self.GetBinContent(begin),                   self.GetBinContent(end),
                     self.GetBinError(begin),                     self.GetBinError(end),
                    )

    @property
    @cache
    def binwidth(self):
        assert len(set(self.h.GetBinWidth(i) for i in range(1, self.GetNbinsX()+1))) == 1
        return self.h.GetBinWidth(1)

    def __getattr__(self, attr):
        try:
            return getattr(self.h, attr)
        except AttributeError as e:
            if "'{}'".format(attr) in str(e):
                raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, attr))
            else:
                raise

class ControlPlot(object):
    def __init__(self, disc, *args):
        self.template = Template(*args)
        f = tfiles[self.template.templatesfile.templatesfile()]
        controlPlots = f.controlPlots
        axis = self.template.discriminants.index(discriminant(disc))

        self.h = {}

        canvas = controlPlots.Get("control_{}_projAxis{}_afterReweight".format(self.template.templatename(final=False), axis))
        self.h["raw"], self.h["reweight"] = (HistInfo(h) for h in canvas.GetListOfPrimitives())

        canvas = controlPlots.Get("control_{}_projAxis{}_afterSmooth".format(self.template.templatename(final=False), axis))
        _, self.h["smooth"] = (HistInfo(h) for h in canvas.GetListOfPrimitives())

    @property
    @cache
    def binwidth(self):
        _ = set(v.binwidth for k, v in self.h.iteritems())
        assert len(_) == 1
        return _.pop()

    @generatortolist
    def rangesthataresmoothedaway(self, step):
        result = []

        for rawrange in self.ranges("raw"):
            overlaps = [smoothrange for smoothrange in self.ranges(step) if rawrange.overlap(smoothrange)]
            if max(smoothrange.overlap(rawrange) for smoothrange in overlaps) <= rawrange.length/2:  #if it's broken in half or more, with no majority, it's basically gone
                yield rawrange
                continue
            biggestoverlap = max(overlaps, key = lambda x: x.overlap(rawrange))
            if not biggestoverlap.samedirection(rawrange):
                yield rawrange
                continue

    @property
    @generatortolist
    def rangesthataresmoothedawaybutreweightedback(self):
        smoothedaway = self.rangesthataresmoothedaway("smooth")
        stillgoneafterreweight = self.rangesthataresmoothedaway("reweight")
        return [range_ for range_ in smoothedaway if range_ not in stillgoneafterreweight]

    def ranges(self, step):
        return self.h[step].increasingordecreasingranges


def printranges(disc, *args):
    controlplot = ControlPlot(disc, *args)

    print controlplot.rangesthataresmoothedaway("reweight")

    fmt = "    {:8.3g} {:8.3g} {:8.3g} {:8}"

    print "raw ranges:"
    for _ in controlplot.ranges("raw"):
        print fmt.format(_.low, _.hi, _.significance, _ in controlplot.rangesthataresmoothedawaybutreweightedback)
    for step in "smooth", "reweight":
        print step, "ranges:"
        for _ in controlplot.ranges(step):
            print fmt.format(_.low, _.hi, _.significance, "")

if __name__ == "__main__":
    t = Template("WH", "fa3", "D_int_prod", "2e2mu", "VHHadrtagged", "160928", "fa3prod0.5")
    for d in t.discriminants:
        print d
        printranges(d, t)
        print
        print
        print
