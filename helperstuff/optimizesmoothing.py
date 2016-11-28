#!/usr/bin/env python
from abc import ABCMeta, abstractmethod, abstractproperty
from collections import namedtuple
from discriminants import discriminant
from filemanager import cache, tfiles
import itertools
from math import sqrt
import numbers
import numpy
import os
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

def generatortolist(function):
    def newfunction(*args, **kwargs):
        return list(function(*args, **kwargs))
    newfunction.__name__ = function.__name__
    return newfunction

class Interval(namedtuple("Interval", "low hi tolerance")):
    def __init__(self, *args, **kwargs):
        super(Interval, self).__init__(*args, **kwargs)
        if self.low > self.hi:
            raise ValueError("Can't have an interval that goes from {} to {}!\nMaybe try switching the order?".format(self.low, self.hi))

    def __contains__(self, other):
        if isinstance(other, numbers.Number):
            return self.low-self.tolerance <= other <= self.hi+self.tolerance
        return self.low-self.tolerance <= other.low < other.hi <= self.hi+self.tolerance

    def overlapinterval(self, other):
        return Interval(max(self.low, other.low), min(self.hi, other.hi), max(self.tolerance, other.tolerance))

    def overlap(self, other):
        """
        http://stackoverflow.com/a/2953979/5228524
        """
        result = max(0, min(self.hi, other.hi) - max(self.low, other.low))
        if result <= self.tolerance:
            return 0  #so that __bool__ returns the right thing (False)
        return result

    @property
    def length(self):
        return self.hi - self.low

class Range(namedtuple("Range", "low hi tolerance lowval hival lowerr hierr"), Interval):
    @property
    @cache
    def deltax(self):
        return self.hi - self.low
    @property
    @cache
    def deltay(self):
        return abs(self.hival - self.lowval)
    @property
    @cache
    def deltayerror(self):
        return sqrt(self.lowerr**2 + self.hierr**2)

    @property
    @cache
    def significance(self):
        """
        This function is supposed to return a number
        The bigger this number is, the greater the probability
          that this range is a real thing rather than a statistical fluctuation
        """
        return self.deltax * self.deltay / self.deltayerror

    @property
    @cache
    def significanceerror(self):
        """
        Error on significance would have to involve error on error on deltay
        This involves the kurtosis https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
        I think this could in principle be calculated from here https://arxiv.org/pdf/1309.1287.pdf
        But there's no such thing as Sumw4() or something like that
        I'll let 15% stand in, for now
        """
        return .15 * self.significance

    def samedirection(self, other):
        return sign(self.hival - self.lowval) == sign(other.hival - other.lowval)

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
                     self.binwidth * ReweightBinning.tolerancefactor,
                     self.GetBinContent(begin),                   self.GetBinContent(end),
                     self.GetBinError(begin),                     self.GetBinError(end),
                    )

    @property
    @cache
    def binwidth(self):
        assert len(set(self.h.GetBinWidth(i) for i in range(1, self.GetNbinsX()+1))) == 1
        return self.h.GetBinWidth(1)

    @property
    @cache
    def xmin(self):
        return self.h.GetXaxis().GetXmin()

    @property
    @cache
    def xmax(self):
        return self.h.GetXaxis().GetXmax()

    def __getattr__(self, attr):
        try:
            return getattr(self.h, attr)
        except AttributeError as e:
            if "'{}'".format(attr) in str(e):
                raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, attr))
            else:
                raise

class ControlPlotBase(object):
    __metaclass__ = ABCMeta
    def __init__(self, disc, *args, **kwargs):
        self.template = Template(*args)
        self.disc = discriminant(disc)
        f = tfiles[self.template.templatefile(**kwargs)]
        controlPlots = f.controlPlots
        axis = self.template.discriminants.index(self.disc)

        self.h = {}

        for step in "Fill", "Floor", "Smooth", "Reweight", "Normalization":

            canvas = controlPlots.Get("control_{}_projAxis{}_after{}".format(self.template.templatename(final=False), axis, step))
            if canvas:
                _, self.h[step.lower()] = (HistInfo(h) for h in canvas.GetListOfPrimitives())

        self.h["raw"] = _

        self.checkthatthisworks()

    @abstractmethod
    def checkthatthisworks(self):
        """
        Called by __init__.  This should check the parameters of self.template, and maybe other things
        raise WrongControlPlotError if this particular subclass doesn't apply to this template
        """

    @property
    @cache
    def binwidth(self):
        _ = set(v.binwidth for k, v in self.h.iteritems())
        assert len(_) == 1
        return _.pop()

    @property
    @cache
    def xmin(self):
        _ = set(v.xmin for k, v in self.h.iteritems())
        assert len(_) == 1
        return _.pop()

    @property
    @cache
    def xmax(self):
        _ = set(v.xmin for k, v in self.h.iteritems())
        assert len(_) == 1
        return _.pop()

    @property
    @cache
    def nbins(self):
        result = float(self.xmax - self.xmin) / self.binwidth
        if abs(result - int(result)) < 1e-9: return int(result)
        if 1-1e-9 < result - int(result) < 1: return int(result)+1
        assert False

    def isontheedge(self, range):
        #print "=================================="
        #print self.xmin, range.low - self.binwidth/2, self.xmax, range.hi+self.binwidth/2
        #print "=================================="
        return (
                   abs(range.low - self.binwidth/2 - self.xmin) < range.tolerance
                or abs(range.hi  + self.binwidth/2 - self.xmax) < range.tolerance
               )

    @cache
    def overlaps(self, therange, step):
        return sorted(
                      [_ for _ in self.ranges(step) if therange.overlap(_) > therange.tolerance],
                      key = lambda x: therange.overlap(x)
                     )

    def whathappenstoit(self, therange, fromstep, instep):
        """
        returns:
          1 if therange from fromstep is absorbed in a bigger range in instep
         -1 if therange from fromstep is smoothed away in instep
          0 if neither
        """
        if therange not in self.ranges(fromstep):
            raise ValueError("therange {}-{} does not come from fromstep {}!".format(therange.low, therange.hi, fromstep))

        overlaps = self.overlaps(therange, instep)
        if abs(therange.length - self.binwidth*2) < therange.tolerance and len(overlaps) == 2:
            #special case --> 3 bins wide, split in half
            #could just mean it was moved over a bit
            samedirection = [_ for _ in overlaps if _.samedirection(therange)]
            assert len(samedirection) == 1
            samedirection = samedirection[0]
            if samedirection.length+samedirection.tolerance > 5*self.binwidth:
                return 1

        if not overlaps:
            return -1
        biggestoverlap = overlaps[-1]
        if biggestoverlap.overlap(therange) <= therange.length/2 + therange.tolerance:
            #if it's broken in half or more, with no majority, it's basically gone
            #except for the special case above
            return -1
        if not biggestoverlap.samedirection(therange):
            return -1

        if len(overlaps) >= 2:
            if therange.length - biggestoverlap.overlap(therange) + therange.tolerance > therange.length / 5:
                return 0

        overlapoverlaps = [_ for _ in self.overlaps(biggestoverlap, fromstep) if _ != therange]
        for otherrange in overlapoverlaps:
            if otherrange in biggestoverlap: #if biggestoverlap contains at least one other range
                return 1

        return 0

    @cache
    def isrelevant(self, therange):
        maxsignificance = max(_.significance for step in self.h.keys() for _ in self.ranges(step))
        if therange.significance < maxsignificance / 100:
            return False
        return True

    @cache
    @generatortolist
    def rangesinonebutnotintheother(self, inthisone, notinthisone):
        for rawrange in self.ranges(inthisone):
            if self.whathappenstoit(rawrange, inthisone, notinthisone) == -1:
                yield rawrange

    @cache
    def rangesthataresmoothedaway(self, step):
        return self.rangesinonebutnotintheother("raw", step)

    @property
    @cache
    def rangesthataresmoothedawaybutreweightedback(self):
        smoothedaway = self.rangesthataresmoothedaway("smooth")
        stillgoneafterreweight = self.rangesthataresmoothedaway("reweight")
        therawranges = [range_ for range_ in smoothedaway if range_ not in stillgoneafterreweight]
        return therawranges

    @property
    @cache
    def rangesintroducedbyreweighting(self):
        """
        Not entirely sure how these get introduced
        but I've seen some of them
        """
        smoothedaway = self.rangesthataresmoothedaway("smooth")
        stillgoneafterreweight = self.rangesthataresmoothedaway("reweight")
        therawranges = [range_ for range_ in stillgoneafterreweight if range_ not in smoothedaway]
        reweightedranges = [
                            max(
                                (reweightrange for reweightrange in self.ranges("reweight")),
                                key = lambda x: x.overlap(rawrange)
                               ) for rawrange in therawranges
                           ]
        return reweightedranges

    @cache
    @generatortolist
    def rangesthatareabsorbed(self, inthisone, absorbedinthisone):
        """
        Ranges that are absorbed in a bigger range going in the same direction
        """
        for rawrange in self.ranges(inthisone):
            if self.whathappenstoit(rawrange, inthisone, absorbedinthisone) == 1:
                yield rawrange

    @property
    @cache
    def rangesabsorbedbysmoothing(self):
        return self.rangesthatareabsorbed("raw", "smooth")


    def ranges(self, step):
        return self.h[step].increasingordecreasingranges

    @abstractproperty
    @cache
    @generatortolist
    def rangesthatshouldnotbereweighted(self):
        try:
            maxsmoothedsignificance = max(_.significance for _ in self.rangesthataresmoothedaway("reweight"))
        except ValueError:#no ranges are smoothed away
            try:
                maxsmoothedsignificance = max(_.significance for _ in self.rangesthataresmoothedaway("smooth"))
            except ValueError:
                #in this case smoothing didn't do too much, so it's safe to say that reweighting won't ruin it
                return

        for therange in self.rangesthataresmoothedawaybutreweightedback:
            if not self.isrelevant(therange): continue
            significance = therange.significance
            if self.isontheedge(therange):
                significance *= 2
            if significance < maxsmoothedsignificance:
                yield therange

        for therange in self.rangesintroducedbyreweighting:
            if not self.isrelevant(therange): continue
            yield therange

    def GetEffectiveEntries(self, step, *args, **kwargs):
        return self.h[step].GetEffectiveEntries(*args, **kwargs)

    @abstractproperty
    @cache
    @generatortolist
    def rangesthatshouldhavebeensmoothed(self):
        def maybeprint(*stuff):
            for a in stuff:
                #print a,
                pass

        try:
            maxsmoothedsignificance_pluserror = max(
                                                    range.significance + range.significanceerror
                                                      for range in self.rangesthataresmoothedaway("smooth")
                                                   )
        except ValueError:#nothing was smoothed away!  ...I have no idea!
            significances = [_.significance for _ in self.ranges("raw")]
            average = numpy.average(significances)
            stdev = numpy.std(significances)
            for therange in self.ranges("raw"):
                if self.isrelevant(therange) and therange.significance < average-3*stdev:
                    yield therange
            return

        for i, therange in enumerate(self.ranges("raw")):
            if not self.isrelevant(therange): maybeprint(i, therange.low, therange.hi, 0); continue
            if therange in self.rangesthataresmoothedaway("smooth"): maybeprint(i, therange.low, therange.hi, 1); continue
            if therange in self.rangesabsorbedbysmoothing: maybeprint(i, therange.low, therange.hi, 2); continue
            if not self.isontheedge(therange) and therange.significance < maxsmoothedsignificance_pluserror: maybeprint(i, therange.low, therange.hi, 4); yield therange; continue
            if self.isontheedge(therange) and 2*therange.significance < maxsmoothedsignificance_pluserror: maybeprint(i, therange.low, therange.hi, 5); yield therange; continue
            maybeprint(i, therange.low, therange.hi, 3); continue

    @property
    def sufficientsmoothing(self):
        return not self.rangesthatshouldhavebeensmoothed

class WrongControlPlotException(Exception): pass

class ControlPlotSimple(ControlPlotBase):
    def checkthatthisworks(self):
        pass
    @property
    def rangesthatshouldnotbereweighted(self):
        return super(ControlPlotSimple, self).rangesthatshouldnotbereweighted
    @property
    def rangesthatshouldhavebeensmoothed(self):
        return super(ControlPlotSimple, self).rangesthatshouldhavebeensmoothed

class ControlPlotOneBigRange(ControlPlotBase):
    """
    For plots that are just increasing, or just decreasing
    (with possibly a bin or two at each end doing something else)
    """
    @property
    @cache
    @generatortolist
    def rangesthatshouldnotbereweighted(self):
        result = super(ControlPlotOneBigRange, self).rangesthatshouldnotbereweighted
        if result: pass
            #for _ in result: yield _
            #return
        smoothranges = self.ranges("smooth")
        reweightranges = self.ranges("reweight")
        biggestsmoothrange = max((_ for _ in smoothranges), key=lambda x: x.significance)
        if biggestsmoothrange.length <= self.binwidth * (self.nbins/2 + .0001):
            raise WrongControlPlotException
        for reweightrange in reweightranges:
            if (
                reweightrange.overlap(biggestsmoothrange)
                and not reweightrange.samedirection(biggestsmoothrange)
                and not self.isontheedge(reweightrange)
               ):
                print reweightrange.overlapinterval(biggestsmoothrange)
                yield reweightrange.overlapinterval(biggestsmoothrange)

    @property
    def rangesthatshouldhavebeensmoothed(self):
        return super(ControlPlotOneBigRange, self).rangesthatshouldhavebeensmoothed

    def checkthatthisworks(self):
        if "reweight" in self.h:
            self.rangesthatshouldnotbereweighted  #might raise WrongControlPlotException
        t, d = self.template, self.disc
        if d == discriminant("D_bkg_0plus"): return
        if d == discriminant("D_0minus_VBFdecay"):
            if t.productionmode == "VBF" and t.hypothesis not in ("0+", "0-"):
                raise WrongControlPlotException
            else:
                return
        if d == discriminant("D_g2_VBFdecay"):
            if t.productionmode == "VBF" and t.hypothesis not in ("0+", "a2"): #L1 also doesn't work
                raise WrongControlPlotException
            else:
                return
        if d == discriminant("D_g1prime2_VBFdecay"):
            if (
                t.productionmode == "VBF" and t.hypothesis not in ("0+", "L1")
                or t.productionmode in ("ZH", "WH") and t.hypothesis != "0+"
               ):
                raise WrongControlPlotException
            else:
                return
        raise WrongControlPlotException

def ControlPlot(*args, **kwargs):
    for cls in ControlPlotOneBigRange, ControlPlotSimple:
        #ControlPlotSimple should be the last, since it works for everything
        try:
            return cls(*args, **kwargs)
        except WrongControlPlotException:
            continue


def printranges(disc, *args, **kwargs):
    controlplot = ControlPlot(disc, *args, **kwargs)
    print controlplot

    print "smoothed sufficiently:", controlplot.sufficientsmoothing

    fmt = "    {:8.3g} {:8.3g} {:8.3g} {:5.2g} {:5.2g} {:5.2g}"

    print "raw ranges:"
    for _ in controlplot.ranges("raw"):
        print fmt.format(
                         _.low, _.hi, _.significance,
                         _ in controlplot.rangesthataresmoothedaway("smooth"),
                         _ in controlplot.rangesabsorbedbysmoothing,
                         _ in controlplot.rangesthatshouldhavebeensmoothed,
                        )

    for step in "smooth", "reweight":
        if controlplot.rangesthataresmoothedaway(step):
            maxsmoothedsignificance = max(range.significance for range in controlplot.rangesthataresmoothedaway(step))
            maxsmoothedsignificance_pluserror = max(range.significance+range.significanceerror for range in controlplot.rangesthataresmoothedaway(step))
            print step, "maxsmoothedsignificance =", maxsmoothedsignificance, "+/-", maxsmoothedsignificance_pluserror-maxsmoothedsignificance
        else:
            print step, "nothing smoothed away"
        print step, "ranges:"
        for _ in controlplot.ranges(step):
            print fmt.format(
                             _.low, _.hi, _.significance,
                             min(sum(_.overlap(_2) for _2 in controlplot.rangesthataresmoothedawaybutreweightedback) / _.length, 1),
                             _ in controlplot.rangesintroducedbyreweighting,
                             min(sum(_.overlap(_2) for _2 in controlplot.rangesthatshouldnotbereweighted) / _.length, 1),
                            )

#    print
#    print "controlplot.rangesthataresmoothedawaybutreweightedback"
#    for _ in controlplot.rangesthataresmoothedawaybutreweightedback:
#        print "    {:8.3g} {:8.3g}".format(_.low, _.hi)

class ReweightBinning(object):
    tolerancefactor = 1.0 / 100000

    def __init__(self, *args):
        if len(args) == 4:
            #combine intervals
            self.combineintervals, self.nbins, self.xmin, self.xmax = args
            self.binwidth = float(self.xmax - self.xmin) / self.nbins
            self.tolerance = self.binwidth*self.tolerancefactor
            self.sortcombineintervals()
            self.reweightbinningfromcombineintervals()
        elif len(args) == 2:
            #reweight binning
            self.__reweightbinning, self.nbins = args
            self.xmin, self.xmax = self.__reweightbinning[0], self.__reweightbinning[-1]
            self.binwidth = float(self.xmax - self.xmin) / self.nbins
            self.tolerance = self.binwidth*self.tolerancefactor
            self.combineintervalsfromreweightbinning()
        else:
            raise TypeError("Wrong number of args for {}.__init__".format(type(self).__name__))

    def sortcombineintervals(self):
        self.combineintervals.sort(key=lambda x: x.low)
        while True:
            for (i1, interval1), (i2, interval2) in pairwise(enumerate(self.combineintervals[:])):
                if interval1.hi > interval2.low:
                    if interval1.hi >= interval2.hi:  #interval2 is unnecessary
                        self.combineintervals.remove(interval2)
                    else:                             #merge them
                        assert interval1.tolerance == interval2.tolerance
                        self.combineintervals.remove(interval2)
                        self.combineintervals[i1] = Interval(interval1.low, interval2.hi, interval1.tolerance)
                    break
            else:
                break

    def combineintervalsfromreweightbinning(self):
        for bin in self.__reweightbinning:
            if abs((bin - self.xmin) % self.binwidth) > self.tolerance and abs(abs((bin - self.xmin) % self.binwidth) - self.binwidth) > self.tolerance:
                raise ValueError("({!r} - xmin) % binwidth = {!r} != 0".format(bin, (bin - self.xmin) % self.binwidth))

        self.combineintervals = []
        combinefrom = combineto = None
        for i in range(self.nbins+1):
            currentvalue = self.xmin + i*self.binwidth
            if any(abs(bin - currentvalue) < self.tolerance for bin in self.__reweightbinning):   #there is a bin boundary with this value
                if combinefrom is not None:
                    combineto = currentvalue
                    self.combineintervals.append(Interval(combinefrom, combineto, self.tolerance))
                    combinefrom = combineto = None
            else:                                               #no bin boundary with this value
                if combinefrom is None:
                    combinefrom = currentvalue - self.binwidth

        assert combinefrom is None

    def reweightbinningfromcombineintervals(self):
        self.sortcombineintervals()

        for interval in self.combineintervals:
            for _ in interval.low, interval.hi:
                if abs((_ - self.xmin) % self.binwidth) > self.tolerance and abs((_ - self.xmin) % self.binwidth - self.binwidth) > self.tolerance:
                    raise ValueError("({!r} - xmin) % binwidth = {!r} != 0".format(_, (_ - self.xmin) % self.binwidth))
            if interval.tolerance != self.tolerance:
                raise ValueError("Inconsistent tolerances {!r} {!r}".format(self.tolerance, interval.tolerance))

        self.__reweightbinning = []

        combineintervals = iter(self.combineintervals+[None])
        nextcombineinterval = next(combineintervals)
        currentcombineinterval = None
        for i in range(self.nbins+1):
            if currentcombineinterval is not None:
                if self.xmin + (i+1)*self.binwidth not in currentcombineinterval:  #then we are at the upper edge
                    currentcombineinterval = None

            if currentcombineinterval is None:    #not else!!
                self.__reweightbinning.append(self.xmin + i*self.binwidth)
                if nextcombineinterval is not None and self.xmin + i*self.binwidth in nextcombineinterval:
                    currentcombineinterval = nextcombineinterval
                    nextcombineinterval = next(combineintervals)

    def addcombineinterval(self, low, hi):
        self.combineintervals.append(Interval(low, hi, self.tolerance))
        self.reweightbinningfromcombineintervals()

    @property
    def reweightbinning(self):
        if self.combineintervals:
            return self.__reweightbinning
        else:
            return None

class TemplateIterate(Template):
    @property
    @cache
    def controlplots(self):
        return [ControlPlot(disc, self) for disc in self.discriminants]

    def getnextiteration(self):
        message = ""

        if not os.path.exists(self.templatefile()):
            return [None, None, None], "0th iteration-->run with no smoothing"

        if self.smoothingparameters[0] is None:
            #first iteration --> try smoothing with reweight on all axes, with no rebinning
            #                    200 bins if possible, as in the TemplateBuilder examples
            #                    but maximum of 4 bins in the template
            neffectiveentries = self.controlplots[0].GetEffectiveEntries("raw")  #effective entries should be the same for any axis
            entriesperbin = min(neffectiveentries/4, 200)
            nbins = neffectiveentries / entriesperbin
            message += "First iteration-->{} entries per bin ({} bins), ".format(entriesperbin, nbins)
            message += "reweight on all axes with JB's automatic binning procedure."
            return [entriesperbin, [0, 1, 2], None], message

        baddiscriminants = ", ".join(
                                     disc.name for disc, controlplot in zip(self.discriminants, self.controlplots)
                                       if not controlplot.sufficientsmoothing
                                    )
        if baddiscriminants:
            cansmoothmore = True
            entriesperbin = self.smoothingparameters[0]
            neffectiveentries = self.controlplots[0].GetEffectiveEntries("raw")  #effective entries should be the same for any axis
            nbins = neffectiveentries / entriesperbin
            if nbins >= 20:
                entriesperbin *= 2
                nbins = neffectiveentries / entriesperbin
            elif nbins >= 4:
                entriesperbin *= 1.5
                nbins = neffectiveentries / entriesperbin
            else:
                newnbins = max(int(nbins+.5) - 1, 2)
                entriesperbin = int(neffectiveentries / newnbins)
                newnbins = neffectiveentries / entriesperbin #to be correct in the message
                if newnbins == nbins:
                    #um... not much we can do...
                    message += "Would like to smooth more (due to {}), but can't, already at {} bins.  ".format(baddiscriminants, nbins)
                    cansmoothmore = False
                else:
                    nbins = newnbins

            if cansmoothmore:
                message += "Insufficient smoothing ({}), going to {} entries per bin ({} bins)".format(baddiscriminants, entriesperbin, nbins)
                return [entriesperbin, [0, 1, 2], None], message

        #find increasing/decreasing ranges that were not properly smoothed away
        if self.smoothingparameters[2] is None:
            self.smoothingparameters[2] = [None, None, None]
        reweightbinnings = [
                            ReweightBinning(axis, disc.bins)
                              if axis is not None
                              else
                            ReweightBinning([], disc.bins, disc.min, disc.max)
                              for axis, disc in zip(self.smoothingparameters[2], self.discriminants)
                           ]
        for controlplot, binning in zip(self.controlplots, reweightbinnings):
            for range in controlplot.rangesthatshouldnotbereweighted:
                binning.addcombineinterval(range.low - binning.binwidth/2, range.hi + binning.binwidth/2)

        newreweightbinning = [binning.reweightbinning for binning in reweightbinnings]
        baddiscriminants = ", ".join(
                                     disc.name for disc, oldbinning, newbinning in zip(self.discriminants, self.smoothingparameters[2], newreweightbinning)
                                       if oldbinning != newbinning
                                    )
        if baddiscriminants:
            message += "Reweighting with automatic rebinning did not work for {}.  Removing the bad ranges.".format(baddiscriminants)
            return [self.smoothingparameters[0], self.smoothingparameters[1], newreweightbinning], message

        return None

if __name__ == "__main__":
    t = Template(*"4mu  vbf fL1 160928 Untagged VBF L1".split())
    print t
    print
    for d in t.discriminants:
        print d
        printranges(d, t, iteration=1)
        print
        print
        print
