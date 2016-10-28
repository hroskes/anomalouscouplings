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

class HistInfo(object):
    def __init__(self, h):
        self.h = h

    @property
    @cache
    def bins(self):
        return [Bin(self.h.GetBinLowEdge(i), self.h.GetBinLowEdge(i+1), self.h.GetBinContent(i), self.h.GetBinError(i))
                       for i in range(1, self.h.GetNbinsX()+1)]

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
        rangebegin = self.bins[0].low
        nextrangebegin = rangeend = None
        for middle, difference, differror in self.derivative:
            if direction == 0:
                direction = sign(difference)
                if rangebegin is None:
                    rangebegin = middle
            if sign(difference) != direction:
                if rangeend is None:
                    rangeend = middle
                if sign(difference) == 0:
                    nextrangebegin = middle
                if sign(difference) != 0:
                    if nextrangebegin is None:
                        nextrangebegin = middle
                    direction = sign(difference)
                    ranges.append((rangebegin, rangeend))
                    rangebegin, rangeend, nextrangebegin = nextrangebegin, None, None
            else:
                nextrangebegin = rangeend = None
        rangeend = self.bins[-1].hi
        ranges.append((rangebegin, rangeend))
        return ranges

def printranges(disc, *args):
    template = Template(*args)
    f = tfiles[template.templatesfile.templatesfile()]
    controlPlots = f.controlPlots
    axis = template.discriminants.index(discriminant(disc))

    canvas = controlPlots.Get("control_{}_projAxis{}_afterNormalization".format(template.templatename(final=False), axis))
    hraw, hproj = (HistInfo(h) for h in canvas.GetListOfPrimitives())

    print "raw ranges:"
    for low, hi in hraw.increasingordecreasingranges:
        print "   ", low, hi
    print "proj ranges:"
    for low, hi in hproj.increasingordecreasingranges:
        print "   ", low, hi

if __name__ == "__main__":
    t = Template("WH", "fa3", "D_int_prod", "2e2mu", "VHHadrtagged", "160928", "0+")
    for d in t.discriminants:
        print d
        printstuff(d, t)
        print
        print
        print
