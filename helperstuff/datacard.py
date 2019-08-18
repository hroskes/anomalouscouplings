import collections
import functools
import inspect
import multiprocessing
import numbers
import os
import pprint
import re
import sys
import time
import traceback

from collections import Counter
from itertools import chain, izip, product
from math import pi

import ROOT

from rootoverloads import histogramfloor

import config
import utilities

from combinehelpers import discriminants, getdatatree, gettemplate, getnobserved, getrate, Luminosity, mixturesign, sigmaioversigma1, zerotemplate
from enums import Analysis, categories, Category, Channel, channels, EnumItem, Hypothesis, MultiEnum, MyEnum, Production, ProductionMode, ShapeSystematic, SystematicDirection, WorkspaceShapeSystematic
from samples import ReweightingSample
from templates import TemplatesFile
from utilities import cache, callclassinitfunctions, cd, deprecate, generatortolist, mkdir_p, multienumcache, OneAtATime, Tee, WriteOnceDict
from yields import YieldSystematic, YieldSystematicValue

names = set()

def makename(name):
    if not isinstance(name, basestring): raise ValueError("I think you are confused, '{!r}' is not a string!".format(name))
    if name in names: raise ValueError("Name '{}' is already taken!".format(name))
    names.add(name)
    return name

class Section(object):
    def __init__(self, *labels):
        self.labels = []
        for label in labels:
             if isinstance(label, SystematicFromEnums_BaseClass):
                 self.labels += label.allnames
             else:
                 self.labels.append(label)

    def __get__(self, obj, objtype):
        return "\n".join(self.getlines(obj, objtype))
    def getlines(self, obj, objtype):
        for label in self.labels:
            if not isinstance(label, basestring) or label.startswith("#"):
                yield label
            else:
                value = getattr(obj, label)
                if value is None: continue
                if isinstance(value, (basestring, numbers.Number)):
                    yield "{} {}".format(label, value)
                elif isinstance(value, (tuple, list)):
                    yield label, value

class SystematicsSection(Section):
    def __init__(self, *labels):
        self.systematicnames = collections.defaultdict(list)
        super(SystematicsSection, self).__init__(*labels)
    def getlines(self, obj, objtype):
        for line in super(SystematicsSection, self).getlines(obj, objtype):
            if isinstance(line, (tuple, list)):
                if line[1][0] == "group" and len(line) == len(line[1]) == 2:
                    #syntax for group: property should return "group", lambda systematicname: return should_this_systematic_be_in_the_group(systematicname)
                    line = line[0] + " group = " + " ".join(name for name in self.systematicnames[obj] if line[1][1](name))
                else:
                    raise ValueError("Unknown result for "+line[0]+":\n"+str(line))
            if len(line.split()) > 2 and all(systematicvalue == "-" for systematicvalue in line.split()[2:]):
                continue
            if len(line.split()) > 2:
                if re.match("ln[NU]|gmM|trG|shape.*|unif|dfD2?|(p|(flat|rate)P)aram|extArg|discrete|group", line.split()[1]):
                    self.systematicnames[obj].append(line.split()[0])
                else:
                    raise ValueError("Unknown pdf type for line:\n"+line)
            yield line

class SystematicFromEnums_BaseClass(object):
    pass

def MakeSystematicFromEnums(*theenums, **kwargs):
    theenums = list(theenums)
    class SystematicFromEnums(SystematicFromEnums_BaseClass):
        def __init__(self, function):
            self.name = self.origname = function.__name__
            for i, enum in enumerate(theenums[:]):
                if isinstance(enum, int):
                    class Index(MyEnum):
                        enumitems = list(EnumItem(str(i)) for i in xrange(1, enum+1))
                        enumname = "index"
                    enum = theenums[i] = Index
                self.name = self.name.replace(enum.enumname, "{"+enum.enumname+"}")

            for kw, kwarg in kwargs.iteritems():
                if kw == "name":
                    self.name = kwarg
                else:
                    raise ValueError("Unknown kwarg {}={}".format(kw, kwarg))

            class SystematicFromEnumsWithValues(MultiEnum):
                enums = theenums

                def __get__(self_systematic, self_datacard, objtype):
                    if self_datacard is None: return self_systematic #when you get it directly from the Datacard class
                    return function(
                                    self_datacard,
                                    **{enum.enumname: getattr(self_systematic, enum.enumname)
                                             for enum in theenums}
                                   )

            SystematicFromEnumsWithValues.__name__ = self.name

            self.systematicfromenumswithvalues = SystematicFromEnumsWithValues

        @property
        @generatortolist
        def allnames(self):
            cartesianproduct = product(*(enum.items() for enum in theenums))
            for enumvalues in cartesianproduct:
                formatdict = {enum.enumname: enumvalue for enum, enumvalue in zip(theenums, enumvalues)}
                yield self.name.format(**formatdict)


        def applyallfunctions(self_systematic, datacardcls):
            cartesianproduct = product(*(enum.items() for enum in theenums))
            for enumvalues in cartesianproduct:
                formatdict = {enum.enumname: enumvalue for enum, enumvalue in zip(theenums, enumvalues)}
                name = self_systematic.name.format(**formatdict)
                newsystematic = self_systematic.systematicfromenumswithvalues(*enumvalues)
                if hasattr(datacardcls, name):
                    oldsystematic = getattr(datacardcls, name)
                    def doublesystematic(self_datacard, counter=Counter(), oldsystematic=oldsystematic, newsystematic=newsystematic):
                        counter[self_datacard] += 1
                        if counter[self_datacard] == 1: return oldsystematic.__get__(self_datacard, type(self_datacard))
                        if counter[self_datacard] == 2: return newsystematic.__get__(self_datacard, type(self_datacard))
                    doublesystematic.__name__ = name  #does this do anything when I'm going to property it anyway? not sure
                    doublesystematic = property(doublesystematic)
                    setattr(datacardcls, name, doublesystematic)
                else:
                    setattr(datacardcls, name, newsystematic)

            delattr(datacardcls, self_systematic.origname)

    return SystematicFromEnums

class _Datacard(MultiEnum):
    enums = (Analysis, Category, Channel, Luminosity)
    enumname = "datacard"
    def __init__(self, *args):
        super(_Datacard, self).__init__(*args)
        self.pdfs = None
    @property
    def year(self):
        return self.production.year
    @property
    def txtfile(self):
        return "hzz4l_{}S_{}_{}.lumi{:.2f}.txt".format(self.channel, self.category, self.year, float(self.luminosity))
    @property
    def rootfile(self):
        return self.rootfile_base
    @property
    def rootfile_base(self):
        return "hzz4l_{}S_{}_{}.input.root".format(self.channel, self.category, self.year)

    @property
    def productionmodes(self):
        result = [ProductionMode(p) for p in ("ggH", "qqH", "VH", "ttH", "bbH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets")]
        if self.analysis.isdecayonly:
            result.remove("VBF")
            result.remove("VH")
            result.remove("ttH")
            result.remove("bbH")
        if self.production.LHE or self.production.GEN:
            result.remove("ggZZ")
            result.remove("VBFbkg")
        if not config.usedata:
            result.remove("ZX")
        for _ in result[:]:
            if _ == "VBFbkg": deprecate(result.remove(_), 2019, 8, 20)
        return result

    @property
    def imax(self):
        return 1
    @property
    def jmax(self):
        return len(self.getprocesses(Counter()))-1
    @property
    def kmax(self):
        return "*"

    section1 = Section("imax", "jmax", "kmax")

    @property
    def shapes(self):
        return "* * $CHANNEL.input.root $PROCESS $PROCESS_$SYSTEMATIC"

    section2 = Section("shapes")

    @property
    def bin(self, counter=Counter()):
        counter[self] += 1

        bin = "hzz4l_{}S_{}_{}".format(self.channel, self.category, self.year)

        if counter[self] == 1:
            return bin
        elif counter[self] == 2:
            return " ".join([str(bin)]*len(self.getprocesses(Counter())))
        assert False

    @property
    def observation(self):
        return getnobserved(self.channel, self.production, self.category, self.analysis)

    section3 = Section("bin", "observation")

    @property
    @cache
    @generatortolist
    def histograms(self):
        for histogramname in self.allhistograms:
            if self.histogramintegrals[histogramname] == 0: continue
            yield histogramname

    @property
    @cache
    @generatortolist
    def allhistograms(self):
        for p in self.productionmodes:
            if p.isbkg:
                yield p.combinename
            elif p.issignal:
                templategroup = str(p).lower()
                templatesfile = TemplatesFile(templategroup, self.analysis, self.production, self.channel, self.category)
                for t in templatesfile.templates():
                    if not t.hypothesis.ispure: continue
                    name = p.combinename
                    if t.productionmode in ("ggH", "ttH"):
                        if t.productionmode == "ttH" or t.category in ("VBFtagged", "VHHadrtagged"):
                            name += "_" + t.hffhypothesis.combinename
                        else:
                            assert t.hffhypothesis == "Hff0+", t
                    name += "_"+t.hypothesis.combinename
                    yield name
                for t in templatesfile.inttemplates():
                    for sign in "positive", "negative":
                        templatenamepart = (
                          t.templatename().replace("template", "").replace("Mirror", "").replace("Int", "").replace("Hff0MinusHVV", "")
                                          .replace("a1", "g11").replace("a3", "g41").replace("a2", "g21").replace("L1Zg", "ghzgs1prime21").replace("L1", "g1prime21")
                        )
                        if t.templatename().startswith("templateInt"): templatenamepart = "g11"+self.analysis.couplingname+"1"

                        name = p.combinename
                        if t.productionmode in ("ggH", "ttH"):
                            if t.productionmode == "ttH" or t.category in ("VBFtagged", "VHHadrtagged"):
                                name += "_" + t.hffhypothesis.combinename
                            else:
                                assert t.hffhypothesis == "Hff0+", t
                        name += "_"+templatenamepart+"_"+sign
                        yield name
            else:
                assert False

    @property
    def process(self): return " ".join(self.getprocesses())

    def getprocesses(self, counter=Counter()):
        counter[self] += 1

        result = {1: [], 2: []}

        nsignal = nbkg = 0
        for h in self.histograms:
            if "bkg_" in h:
                assert ProductionMode(h).isbkg
                result[1].append(h)
                result[2].append(1 + nbkg)
                nbkg += 1
            else:
                assert ProductionMode(h.split("_")[0]).issignal
                result[1].append(h)
                result[2].append(-nsignal)
                nsignal += 1
        return [str(_) for _ in result[counter[self]]]

    @property
    def rate(self):
        return " ".join(str(self.histogramintegrals[h]) for h in self.histograms)

    section4 = Section("## mass window [{},{}]".format(config.m4lmin, config.m4lmax),
                       "bin", "process", "process", "rate")

    @MakeSystematicFromEnums(YieldSystematic)
    def yieldsystematic(self, yieldsystematic):
        productionmodes = [ProductionMode(h if "bkg_" in h else h.split("_")[0]) for h in self.histograms]

        if self.production.LHE or self.production.GEN:
            if yieldsystematic == "QCDscale_ren_VV":
                result = " ".join(["lnN"] + ["1.1" if h=="bkg_qqzz" else "-" for h in productionmodes])
                return result
            else: return None

        lst = ["lnN"]

        for p in productionmodes:
            ysv = str(YieldSystematicValue(yieldsystematic, self.channel, self.category, self.analysis, self.production, p))
            try:
                wss = WorkspaceShapeSystematic(str(yieldsystematic))
            except ValueError:
                pass
            else:
                if wss in p.workspaceshapesystematics(self.category) and ysv != "-":
                    if (yieldsystematic == "JEC" and self.category in ("VBFtagged", "VHHadrtagged")):
                        ysv = "1"
                        lst[0] = "shape1?"
                    else:
                        raise ValueError("{} has both a yield and shape systematic in {} {}".format(wss, p, self.category))
            lst.append(ysv)

        return " ".join(lst)

    @MakeSystematicFromEnums(WorkspaceShapeSystematic, Channel)
    def workspaceshapesystematic_channel(self, workspaceshapesystematic, channel):
      if self.year not in workspaceshapesystematic.years: return None
      if self.production.LHE or self.production.GEN: return None
      if workspaceshapesystematic.isperchannel and channel == self.channel:
        return " ".join(
          ["shape1"] + [
            "1" if workspaceshapesystematic in
              ProductionMode(h if "bkg_" in h else h.split("_")[0]).workspaceshapesystematics(self.category)
            else "-"
            for h in self.histograms
          ]
        )

    @MakeSystematicFromEnums(WorkspaceShapeSystematic)
    def workspaceshapesystematic(self, workspaceshapesystematic):
      if self.year not in workspaceshapesystematic.years: return None
      if workspaceshapesystematic.isperchannel: return None
      try:
          YieldSystematic(str(workspaceshapesystematic))
      except ValueError:
          pass
      else:
          return None  #in that case shape1? is taken care of in yieldsystematic
      return " ".join(
        ["shape1"] + [
          "1" if workspaceshapesystematic in
            ProductionMode(h if "bkg_" in h else h.split("_")[0]).workspaceshapesystematics(self.category)
          else "-"
          for h in self.histograms
        ]
      )

    @MakeSystematicFromEnums(Category, Channel, config.staticmaxbins)
    def binbybin_category_channel_background_index(self, category, channel, index):
        if category != self.category or channel != self.channel:
          return None
        if "binbybin_{}_{}_background_{}".format(category, channel, index) not in self.binbybinuncertainties:
          return None
        result = ["shape1"]
        return " ".join(["shape1"] + ["1" if "bkg_" in h else "-" for h in self.histograms])

    @property
    def everything_but_binbybin(self):
        return "group", lambda systematicname: "binbybin" not in systematicname

    @property
    def binbybin(self):
        return "group", lambda systematicname: "binbybin" in systematicname

    systematicssection = section5 = SystematicsSection(yieldsystematic, workspaceshapesystematic_channel, workspaceshapesystematic, binbybin_category_channel_background_index, "binbybin", "everything_but_binbybin")

    divider = "\n------------\n"

    __initedsystematics = False

    @classmethod
    def initsystematicsfromenums(cls):
        if cls.__initedsystematics: return
        cls.__initedsystematics = True
        for name, systematic in inspect.getmembers(cls, predicate=lambda x: isinstance(x, SystematicFromEnums_BaseClass)):
            systematic.applyallfunctions(cls)

    def writedatacard(self):
        self.initsystematicsfromenums()
        sections = self.section1, self.section2, self.section3, self.section4, self.section5
        if not os.path.exists(self.rootfile):
            raise IOError("workspace file {} should exist first!".format(self.rootfile))
        with open(self.txtfile, "w") as f:
            f.write(self.divider.join(sections)+"\n")

    @property
    @cache
    def histogramintegrals(self):
        """
        Gets modified during makehistograms
        """
        return {}

    def makehistograms(self):
        f = ROOT.TFile(self.rootfile_base, "RECREATE")
        cache = WriteOnceDict()
        self.binbybinuncertainties = []
        print self
        domirror = False  #will be set to true
        for h in chain(self.allhistograms, ["data"]):
            if "bkg_" in h or h == "data":
                p = h
                hypothesis = None
                hffhypothesis = None
                sign = None
            elif "positive" in h or "negative" in h:
                if "ff" in h:
                    p, hffhypothesis, inttype, sign = h.split("_")
                else:
                    p, inttype, sign = h.split("_")
                for letter, name in zip("ijkl", self.analysis.couplingnames):
                    inttype = inttype.replace(name, "g"+letter)
                hypothesis = inttype
            elif "ff" in h:
                p, hffhypothesis, hypothesis = h.split("_")
                sign = None
            else:
                p, hypothesis = h.split("_")
                hffhypothesis = sign = None

            p = ProductionMode(p)

            for shapesystematic, direction in chain(product(p.workspaceshapesystematics(self.category), ("Up", "Down")), [(None, None)]):
                if shapesystematic is None is direction:
                    systematic = None
                else:
                    if self.production.GEN: continue
                    if p == "data": continue
                    systematic = ShapeSystematic(shapesystematic.nickname+direction)

                name = h
                if systematic: name += "_"+shapesystematic.combinename(self.channel)+direction
                if p == "data": name = "data_obs"

                if p == "data":
                    scaleby = 1 if config.unblindscans else 0
                else:
                    numerator = getrate(p, self.channel, self.category, self.analysis, self.production, 1)
                    denominator = gettemplate(p, self.analysis, self.production, self.category, "0+" if hypothesis else None, "Hff0+" if hffhypothesis else None, self.channel, systematic).Integral()
                    if systematic:
                        try:
                            yieldsystematic = YieldSystematic(str(shapesystematic))
                        except ValueError:
                            pass
                        else:
                            ysv = YieldSystematicValue(p, yieldsystematic, self.analysis, self.production, self.category, self.channel)
                            numerator *= {"Up": ysv.upvalue, "Down": ysv.downvalue}[direction]
                    scaleby = numerator / denominator * float(self.luminosity)

                assert scaleby >= 0, scaleby

                originaltemplate = gettemplate(p, self.analysis, self.production, self.category, hypothesis, hffhypothesis, self.channel, systematic)
                t3D = originaltemplate.Clone(name+"_3D")
                if sign is not None:
                    if sign == "positive": pass
                    elif sign == "negative": t3D.Scale(-1)
                    else: assert False
                    t3D.Floor(0)
                if "Mirror" in originaltemplate.GetName(): domirror = True

                t3D.Scale(scaleby)

                nbinsx = t3D.GetNbinsX()
                nbinsy = t3D.GetNbinsY()
                nbinsz = t3D.GetNbinsZ()
                nbinsxyz = nbinsx*nbinsy*nbinsz

                t = getattr(ROOT, type(t3D).__name__.replace("3", "1"))(name, name, nbinsxyz, 0, nbinsxyz)

                binindices = {}

                xyz = 1
                for x, y, z in product(xrange(1, nbinsx+1), xrange(1, nbinsy+1), xrange(1, nbinsz+1)):
                    binindices[t3D.GetXaxis().GetBinCenter(x), t3D.GetYaxis().GetBinCenter(y), t3D.GetZaxis().GetBinCenter(z)] = xyz
                    binindices[xyz] = t3D.GetXaxis().GetBinCenter(x), t3D.GetYaxis().GetBinCenter(y), t3D.GetZaxis().GetBinCenter(z)
                    t.SetBinContent(xyz, t3D.GetBinContent(x, y, z))
                    t.SetBinError(xyz, t3D.GetBinError(x, y, z))
                    xyz += 1

                t3D.SetDirectory(f)
                assert t3D.GetName() not in cache, t3D.GetName()
                cache[t3D.GetName()] = t3D
                t.SetDirectory(f)
                assert t.GetName() not in cache, t.GetName()
                cache[t.GetName()] = t

                assert t.Integral() == t3D.Integral(), (t.Integral(), t3D.Integral())
                if systematic is not None: self.histogramintegrals[name] = t.Integral()

                if systematic is None and config.usebinbybin and p != "data" and p.isbkg:
                    for i in xrange(1, t.GetNbinsX()+1):
                        if i > config.staticmaxbins:
                            raise ValueError("config.staticmaxbins is not big enough.  If you want to use bin by bin uncertainties, increase it.")
                        xcenter, ycenter, zcenter = binindices[i]
                        systname = "binbybin_{self.category}_{self.channel}_background_{}".format(i, self=self)
                        newname = name + "_" + systname
                        if domirror:
                            if ycenter<0: continue
                            otheri = binindices[xcenter, -ycenter, zcenter]

                        up = t.Clone(newname+"Up")
                        dn = t.Clone(newname+"Down")

                        up.SetBinContent(i, t.GetBinContent(i) + t.GetBinError(i))
                        dn.SetBinContent(i, max(t.GetBinContent(i) - t.GetBinError(i), t.GetBinContent(i)/2))

                        if domirror:
                            up.SetBinContent(otheri, t.GetBinContent(otheri) + t.GetBinError(otheri))
                            dn.SetBinContent(otheri, max(t.GetBinContent(otheri) - t.GetBinError(otheri), t.GetBinContent(otheri)/2))

                        up.SetDirectory(f)
                        cache[up.GetName()] = up
                        dn.SetDirectory(f)
                        cache[dn.GetName()] = dn

                        self.binbybinuncertainties.append(systname)

        for name in sorted(cache):
            if name == "data_obs" or name.endswith("Up") or name.endswith("Down") or name.endswith("_3D"): continue
            if "_negative" in name: continue  #handled by positive

            if "_positive" in name:
                namepositive = name
                namenegative = namepositive.replace("_positive", "_negative")

                systematicnamespositive = sorted([k for k in cache if k.startswith(namepositive) and (k.endswith("Up") or k.endswith("Down"))])
                systematicnamesnegative = sorted([k for k in cache if k.startswith(namenegative) and (k.endswith("Up") or k.endswith("Down"))])

                assert all(neg == pos.replace("_positive", "_negative") for pos, neg in izip(systematicnamespositive, systematicnamesnegative)), (systematicnamespositive, systematicnamesnegative)

                histogramspositive = [cache[namepositive]] + [cache[_] for _ in systematicnamespositive]
                histogramsnegative = [cache[namenegative]] + [cache[_] for _ in systematicnamesnegative]

                for x in xrange(1, histogramspositive[0].GetNbinsX()+1):
                    bincontentpositive = [_.GetBinContent(x) for _ in histogramspositive]
                    bincontentnegative = [_.GetBinContent(x) for _ in histogramsnegative]

                    if any(bincontentpositive) and not all(bincontentpositive) or any(bincontentnegative) and not all(bincontentnegative):
                        print namepositive, x
                        for hpos, hneg, pos, neg in izip(histogramspositive, histogramsnegative, bincontentpositive, bincontentnegative):
                            if not pos or not neg:
                                hpos.Fill(hpos.GetBinCenter(x), 1e-10)
                                hneg.Fill(hneg.GetBinCenter(x), 1e-10)

                self.histogramintegrals[namepositive] = cache[namepositive].Integral()
                self.histogramintegrals[namenegative] = cache[namenegative].Integral()

            else:
                systematicnames = sorted([k for k in cache if k.startswith(name) and (k.endswith("Up") or k.endswith("Down"))])
                histograms = [cache[name]] + [cache[_] for _ in systematicnames]

                for x in xrange(1, histograms[0].GetNbinsX()+1):
                    bincontent = [_.GetBinContent(x) for _ in histograms]
                    if any(bincontent) and not all(bincontent):
                        print namepositive, x
                        for h, val in izip(histograms, bincontent):
                            if not val:
                                h.Fill(x, 1e-10)

                self.histogramintegrals[name] = cache[name].Integral()

        f.Write()
        f.Close()

    def makeCardsWorkspaces(self, outdir="."):
        mkdir_p(outdir)
        with cd(outdir):
            Datacard(self.channel, self.category, self.analysis, self.luminosity).makehistograms()
            self.writedatacard()

Datacard = multienumcache(_Datacard)

def makeDCsandWSs(productions, channels, categories, *otherargs, **kwargs):
    with OneAtATime("makeDCsandWSs.tmp", 30):
        fullarglist = [
          (production, channel, category) + otherargs
          for production, channel, category
          in product(productions, channels, categories)
        ]
        arglist = []
        for dcargs in fullarglist:
          dc = Datacard(*dcargs)
          if all(os.path.exists(thing) for thing in (dc.rootfile_base, dc.rootfile, dc.txtfile)): continue
          if dc.production.LHE and dc.channel != "2e2mu": continue
          if dc.analysis.isdecayonly and dc.category != "Untagged": continue
          if not dc.analysis.useboosted and dc.category == "Boosted": continue
          if not dc.analysis.usemorecategories and dc.category in ("VHLeptTagged", "VBF1jtagged"): continue
          arglist.append(dcargs)
        pool = multiprocessing.Pool(processes=8 if utilities.LSB_JOBID() else 1)
        mapresult = pool.map(functools.partial(makeDCandWS, **kwargs), arglist)
        pool.close()

def makeDCandWS(dcargs, **kwargs):
    try:
        dc = Datacard(*dcargs)
        dc.makeCardsWorkspaces(**kwargs)
        for thing in dc.rootfile_base, dc.rootfile, dc.txtfile:
            if not os.path.exists(thing):
                raise ValueError("{} was not created.  Something is wrong.".format(thing))
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
