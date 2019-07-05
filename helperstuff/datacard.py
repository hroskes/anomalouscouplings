import collections
import inspect
import numbers
import os
import re
import time

from collections import Counter
from itertools import chain, product
from math import pi

import ROOT

from rootoverloads import histogramfloor

import config
import utilities

from combinehelpers import discriminants, getdatatree, gettemplate, getnobserved, getrate, Luminosity, mixturesign, sigmaioversigma1, zerotemplate
from enums import Analysis, categories, Category, Channel, channels, EnumItem, Hypothesis, MultiEnum, MyEnum, Production, ProductionMode, ShapeSystematic, SystematicDirection, WorkspaceShapeSystematic
from samples import ReweightingSample
from templates import TemplatesFile
from utilities import cache, callclassinitfunctions, cd, deprecate, generatortolist, mkdir_p, multienumcache, OneAtATime, Tee
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
        if self.analysis.usehistogramsforcombine: return self.rootfile_base
        return "hzz4l_{}S_{}_{}.lumi{:.2f}.input.root".format(self.channel, self.category, self.year, float(self.luminosity))
    @property
    def rootfile_base(self):
        return "hzz4l_{}S_{}_{}.input.root".format(self.channel, self.category, self.year)

    @property
    def productionmodes(self):
        result = [ProductionMode(p) for p in ("ggH", "qqH", "WH", "ZH", "ttH", "bbH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets")]
        if self.analysis.isdecayonly:
            result.remove("VBF")
            result.remove("WH")
            result.remove("ZH")
            result.remove("ttH")
            result.remove("bbH")
        if self.production.LHE or self.production.GEN:
            result.remove("ggZZ")
            result.remove("VBFbkg")
        if not config.usedata:
            result.remove("ZX")
        for _ in result[:]:
            if _ == "WH": deprecate(result.remove(_), 2019, 7, 5, hour=10)
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
        if self.analysis.usehistogramsforcombine:
            return "* * $CHANNEL.input.root $PROCESS $PROCESS_$SYSTEMATIC"
        else:
            return "* * $CHANNEL.input.root w:$PROCESS w:$PROCESS_$SYSTEMATIC"

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
            if self.analysis.mergeallVVHandallffH:
                if "ttH" in histogramname or "bbH" in histogramname: continue
                histogramname = histogramname.replace("ggH", "ffH")
                if "ZH" in histogramname or "WH" in histogramname: continue
                histogramname = histogramname.replace("qqH", "VVH")
            if self.histogramintegrals[histogramname] == 0: continue
            yield histogramname

    @property
    @cache
    @generatortolist
    def allhistograms(self):
        if not self.analysis.usehistogramsforcombine:
            raise ValueError("Should not be calling this function for {}".format(self.analysis))
        for p in self.productionmodes:
            if p.isbkg:
                yield p.combinename
            elif p.issignal:
                templategroup = str(p).lower()
                templatesfile = TemplatesFile(templategroup, self.analysis, self.production, self.channel, self.category)
                for t in templatesfile.templates():
                    if not t.hypothesis.ispure: continue
                    yield p.combinename+"_"+t.hypothesis.combinename
                for t in templatesfile.inttemplates():
                    for sign in "positive", "negative":
                        templatenamepart = (
                          t.templatename().replace("template", "").replace("AdapSmooth", "").replace("Mirror", "").replace("Int", "")
                                          .replace("a1", "g11").replace("a3", "g41").replace("a2", "g21").replace("L1Zg", "ghzgs1prime21").replace("L1", "g1prime21")
                        )
                        if t.templatename().startswith("templateIntAdapSmooth"): templatenamepart = "g11"+self.analysis.couplingname+"1"

                        yield (p.combinename+"_"+templatenamepart+"_"+sign)
            else:
                assert False

    @property
    def process(self): return " ".join(self.getprocesses())

    def getprocesses(self, counter=Counter()):
        counter[self] += 1

        if self.analysis.usehistogramsforcombine:
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
        else:
            if counter[self] == 1:
                return [_.combinename for _ in self.productionmodes]

            if counter[self] == 2:
                nsignal = sum(_.issignal for _ in self.productionmodes)
                nbkg = sum(_.isbkg for _ in self.productionmodes)
                assert nsignal+nbkg == len(self.productionmodes)

                return [str(_) for _ in range(-nsignal+1, nbkg+1)]

        assert False

    @property
    def rate(self):
        if self.analysis.usehistogramsforcombine:
            return " ".join(str(self.histogramintegrals[h]) for h in self.histograms)
        else:
            return " ".join(str(getrate(p, self.channel, self.category, self.analysis, self.luminosity)) for p in self.productionmodes)

    section4 = Section("## mass window [{},{}]".format(config.m4lmin, config.m4lmax),
                       "bin", "process", "process", "rate")

    @MakeSystematicFromEnums(YieldSystematic)
    def yieldsystematic(self, yieldsystematic):
        if self.analysis.usehistogramsforcombine:
            productionmodes = [h if "bkg_" in h else h.split("_")[0] for h in self.histograms]
        else:
            productionmodes = self.productionmodes

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
                    if (yieldsystematic in ("QCDscale_ggH2in", "CMS_scale_j_13TeV_2016", "CMS_scale_j_13TeV_2017") and p == "ggH" and self.category in ("VBFtagged", "VHHadrtagged")):
                        ysv = "1"
                        lst[0] = "shape1?"
                    else:
                        raise ValueError("{} has both a yield and shape systematic in {} {}".format(wss, p, self.category))
            lst.append(ysv)

        return " ".join(lst)

    @MakeSystematicFromEnums(WorkspaceShapeSystematic, Channel)
    def workspaceshapesystematicchannel(self, workspaceshapesystematic, channel):
      if self.year not in workspaceshapesystematic.years: return None
      if self.production.LHE or self.production.GEN: return None
      if workspaceshapesystematic.isperchannel and channel == self.channel:
        if self.analysis.usehistogramsforcombine:
          return " ".join(
                          ["shape1"] +
                          ["1" if workspaceshapesystematic in
                                  ProductionMode(h if "bkg_" in h else h.split("_")[0]).workspaceshapesystematics(self.category)
                               else "-"
                              for h in self.histograms]
                         )
        return " ".join(
                        ["shape1"] +
                        ["1" if workspaceshapesystematic in p.workspaceshapesystematics(self.category)
                             else "-"
                            for p in self.productionmodes]
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
      if self.analysis.usehistogramsforcombine:
        return " ".join(
                        ["shape1"] +
                        ["1" if workspaceshapesystematic in
                                ProductionMode(h if "bkg_" in h else h.split("_")[0]).workspaceshapesystematics(self.category)
                             else "-"
                            for h in self.histograms]
                       )
      return " ".join(
                      ["shape1"] +
                      ["1" if workspaceshapesystematic in p.workspaceshapesystematics(self.category) else "-"
                          for p in self.productionmodes]
                     )

    @MakeSystematicFromEnums(Channel)
    def CMS_zz4l_smd_zjets_bkg_channel(self, channel):
        if config.usenewZXsystematics: return None
        if self.production.LHE or self.production.GEN: return None
        category = "Untagged"
        if channel == self.channel and category == self.category:
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category != "Untagged":
                return "param 0 1 [-3,3]"

    @MakeSystematicFromEnums(Category)
    def CMS_zz4l_smd_zjets_bkg_category(self, category):
        if config.usenewZXsystematics: return None
        if category == "Untagged": return None
        if not config.mergeZXVBFVHsystematics: return None
        if self.production.LHE or self.production.GEN: return None
        if category == self.category:
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category != "Untagged":
                return "param 0 1 [-3,3]"

    @MakeSystematicFromEnums(Category, Channel)
    def CMS_zz4l_smd_zjets_bkg_category_channel(self, category, channel):
        if config.usenewZXsystematics: return None
        if category == "Untagged": return None
        if config.mergeZXVBFVHsystematics: return None
        if self.production.LHE or self.production.GEN: return None
        if category == self.category and channel == self.channel:
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category != "Untagged":
                return "param 0 1 [-3,3]"

    @MakeSystematicFromEnums(Channel)
    def CMS_fake_channel(self, channel):
        if not config.usenewZXsystematics: return None
        if self.production.LHE or self.production.GEN: return None
        if channel == self.channel:
          return "param 0 1 [-3,3]"

    @MakeSystematicFromEnums(ProductionMode, Category, Channel, config.staticmaxbins)
    def binbybin_category_channel_productionmode_index(self, productionmode, category, channel, index):
        if productionmode == "data" or not productionmode.isbkg: return None
        if not self.analysis.usehistogramsforcombine:
          return None
        if category != self.category or channel != self.channel:
          return None
        if "binbybin_{}_{}_{}_{}".format(category, channel, productionmode, index) not in self.binbybinuncertainties:
          return None
        result = ["shape1"]
        return " ".join(["shape1"] + ["1" if productionmode == (h if "bkg_" in h else h.split("_")[0]) else "-" for h in self.histograms])

    @property
    def kbkg_gg(self):
        if self.production.LHE or self.production.GEN: return None
        return "param 1 0.1 [0,2]"

    @property
    def everything_but_binbybin(self):
        return "group", lambda systematicname: "binbybin" not in systematicname

    @property
    def binbybin(self):
        return "group", lambda systematicname: "binbybin" in systematicname

    systematicssection = section5 = SystematicsSection(yieldsystematic, workspaceshapesystematicchannel, workspaceshapesystematic, CMS_zz4l_smd_zjets_bkg_channel, CMS_zz4l_smd_zjets_bkg_category, CMS_zz4l_smd_zjets_bkg_category_channel, binbybin_category_channel_productionmode_index, CMS_fake_channel, "kbkg_gg", "binbybin", "everything_but_binbybin")

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

    @classmethod
    @cache
    def fai(cls):
        name = makename("CMS_zz4l_fai1")
        fai = ROOT.RooRealVar(name, name, -1., 1.)
        fai.setBins(1000)
        return fai
    @classmethod
    @cache
    def faj(cls):
        name = makename("CMS_zz4l_fai2")
        faj = ROOT.RooRealVar(name, name, -1., 1.)
        faj.setBins(1000)
        return faj
    @classmethod
    @cache
    def fak(cls):
        name = makename("CMS_zz4l_fai3")
        fak = ROOT.RooRealVar(name, name, -1., 1.)
        fak.setBins(1000)
        return fak
    @classmethod
    @cache
    def fal(cls):
        name = makename("CMS_zz4l_fai4")
        fal = ROOT.RooRealVar(name, name, -1., 1.)
        fal.setBins(1000)
        return fal
    @classmethod
    @cache
    def phiai(cls):
        name = makename("CMS_zz4l_phiai1")
        phiai = ROOT.RooRealVar(name, name, -pi, pi)
        phiai.setBins(1000)
        return phiai
    @classmethod
    @cache
    def phiaj(cls):
        name = makename("CMS_zz4l_phiai2")
        phiaj = ROOT.RooRealVar(name, name, -pi, pi)
        phiaj.setBins(1000)
        return phiaj
    @classmethod
    @cache
    def a1(cls, analysis):
        if analysis.dimensions == 1:
            name = makename("a1")
            return ROOT.RooFormulaVar(name, name, "sqrt(1-abs(@0))", ROOT.RooArgList(cls.fai()))
        elif analysis.dimensions == 4:
            name = makename("a1")
            return ROOT.RooFormulaVar(name, name, "sqrt(1-abs(@0)-abs(@1)-abs(@2)-abs(@3)-abs(@4))", ROOT.RooArgList(
                cls.fai(), cls.faj(), cls.fak(), cls.fal()
            ))
    @classmethod
    @cache
    def aletter(cls, analysis, letter):
        index = "ijkl".index(letter)
        name = makename("a"+letter)
        return ROOT.RooFormulaVar(
            name, name,
            "{mixturesign} * (@0>0 ? 1 : -1) * sqrt(abs(@0)/{sigmaioversigma1})".format(
                sigmaioversigma1=sigmaioversigma1(analysis.fais[index], "ggH"),
                mixturesign=mixturesign(analysis.fais[index])
            ),
            ROOT.RooArgList(getattr(cls, "fa"+letter)())
        )

    @classmethod
    def ai(cls, analysis): return cls.aletter(analysis, "i")
    @classmethod
    def aj(cls, analysis): return cls.aletter(analysis, "j")
    @classmethod
    def ak(cls, analysis): return cls.aletter(analysis, "k")
    @classmethod
    def al(cls, analysis): return cls.aletter(analysis, "l")

    def makepdfs(self):
        if self.analysis.usehistogramsforcombine:
            raise ValueError("Should not be calling this function for {}".format(self.analysis))
        if self.pdfs is not None: return

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##

        fai = self.fai()
        if self.analysis.dimensions == 2:
            faj = self.faj()
            phiai = self.phiai()
            phiaj = self.phiaj()
        else:
            faj = phiai = phiaj = None

        if self.analysis.isdecayonly:
            a1 = ai = None
        else:
            a1 = self.a1(self.analysis)
            ai = self.ai(self.analysis)

        discs = discriminants(self.analysis, self.category, self.production)
        D1Name, D2Name, D3Name = (d.name for d in discs)
        dBinsX, dBinsY, dBinsZ = (d.bins for d in discs)
        dLowX, dLowY, dLowZ = (d.min for d in discs)
        dHighX, dHighY, dHighZ = (d.max for d in discs)

        D1 = ROOT.RooRealVar(D1Name, D1Name, dLowX, dHighX)
        D2 = ROOT.RooRealVar(D2Name, D2Name, dLowY, dHighY)
        D3 = ROOT.RooRealVar(D3Name, D3Name, dLowZ, dHighZ)
        D1.setBins(dBinsX)
        D2.setBins(dBinsY)
        D3.setBins(dBinsZ)

        self.pdfs = []

        for p in self.productionmodes:
            pdfkwargs = {"fai": fai, "a1": a1, "ai": ai, "D1": D1, "D2": D2, "D3": D3, "faj": faj, "phiai": phiai, "phiaj": phiaj}
            self.pdfs.append(Pdf(self, p, **pdfkwargs))
            for systematic in p.workspaceshapesystematics(self.category):
                if self.production.year not in systematic.years: continue
                if self.production.LHE or self.production.GEN: continue
                self.pdfs.append(Pdf(self, p, systematic, "Up", **pdfkwargs))
                self.pdfs.append(Pdf(self, p, systematic, "Down", **pdfkwargs))

        ## --------------------------- DATASET --------------------------- ##

        data_obs_tree = getdatatree(self.channel, self.production, self.category, self.analysis)
        datasetName = makename("data_obs_{}_{}_{}".format(self.channel, self.category, self.production))


        self.data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(D1,D2,D3))


        ## --------------------------- WORKSPACE -------------------------- ##

    def writeworkspace(self):
        if os.path.exists(self.rootfile_base): return
        w = ROOT.RooWorkspace("w","w")

        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf_1D_fast_a1ai.Class(),True)
        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf_2D_fast.Class(),True)
        w.importClassCode(ROOT.VVHZZ4L_RooSpinZeroPdf_1D_fast_a1ai.Class(),True)

        getattr(w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        for pdf in self.pdfs:
            getattr(w, 'import')(pdf.pdf, ROOT.RooFit.RecycleConflictNodes())
            if (pdf.productionmode.issignal or pdf.productionmode in ("ggZZ", "ZX")) and pdf.shapesystematic == "":
                getattr(w, 'import')(pdf.norm, ROOT.RooFit.RecycleConflictNodes())

        w.writeToFile(self.rootfile_base)

    def linkworkspace(self):
        if not os.path.exists(self.rootfile):
            os.symlink(self.rootfile_base, self.rootfile)

    @property
    @cache
    def histogramintegrals(self):
        """
        Gets modified during makehistograms
        """
        return {}

    def makehistograms(self):
        if not self.analysis.usehistogramsforcombine:
            raise ValueError("Should not be calling this function for {}".format(self.analysis))

        f = ROOT.TFile(self.rootfile_base, "RECREATE")
        cache = {}
        self.binbybinuncertainties = []
        print self
        domirror = False  #will be set to true
        for h in chain(self.allhistograms, ["data"]):
            if "bkg_" in h or h == "data":
                p = h
                hypothesis = None
                sign = None
            elif "positive" in h or "negative" in h:
                p, inttype, sign = h.split("_")
                for letter, name in zip("ijkl", self.analysis.couplingnames):
                    inttype = inttype.replace(name, "g"+letter)
                hypothesis = inttype
            else:
                p, hypothesis = h.split("_")
                sign = None

            p = ProductionMode(p)

            for systematic, direction in chain(product(p.workspaceshapesystematics(self.category), ("Up", "Down")), [(None, None)]):
                if systematic is not None is not direction:
                    if self.production.GEN: continue
                    if p == "data": continue
                    systematic = ShapeSystematic(str(systematic)+direction)

                name = h
                if systematic: name += "_"+str(systematic)
                if p == "data": name = "data_obs"

                if p == "data":
                    scaleby = 1 if config.unblindscans else 0
                else:
                    scaleby = (
                                getrate(p, self.channel, self.category, self.analysis, self.production, 1)
                              /
                                gettemplate(p, self.analysis, self.production, self.category, "0+" if hypothesis else None, self.channel).Integral()
                              ) * float(self.luminosity)
                assert scaleby >= 0, scaleby

                originaltemplate = gettemplate(p, self.analysis, self.production, self.category, hypothesis, self.channel, systematic)
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
                cache[t3D.GetName()] = t3D
                t.SetDirectory(f)
                cache[t.GetName()] = t

                assert t.Integral() == t3D.Integral(), (t.Integral(), t3D.Integral())
                self.histogramintegrals[name] = t.Integral()

                if systematic is None and config.usebinbybin and p != "data" and p.isbkg:
                    for i in xrange(1, t.GetNbinsX()+1):
                        if i > config.staticmaxbins:
                            raise ValueError("config.staticmaxbins is not big enough.  If you want to use bin by bin uncertainties, increase it.")
                        xcenter, ycenter, zcenter = binindices[i]
                        systname = "binbybin_{self.category}_{self.channel}_{}_{}".format(p, i, self=self)
                        newname = name + "_" + systname
                        if domirror:
                            if ycenter<0: continue
                            otheri = binindices[xcenter, -ycenter, zcenter]

                        up = t.Clone(newname+"Up")
                        dn = t.Clone(newname+"Down")

                        up.SetBinContent(i, t.GetBinContent(i) + t.GetBinError(i))
                        dn.SetBinContent(i, max(t.GetBinContent(i) - t.GetBinError(i), 1e-9))

                        if domirror:
                            up.SetBinContent(otheri, t.GetBinContent(otheri) + t.GetBinError(otheri))
                            dn.SetBinContent(otheri, max(t.GetBinContent(otheri) - t.GetBinError(otheri), 1e-9))

                        up.SetDirectory(f)
                        cache[up.GetName()] = up
                        dn.SetDirectory(f)
                        cache[dn.GetName()] = dn

                        self.binbybinuncertainties.append(systname)

        for k in sorted(cache):
            v = cache[k]
            newname = None
            if self.analysis.mergeallVVHandallffH:
                if "qqH" in k or "ZH" in k or "WH" in k:
                    newname = k.replace("qqH", "VVH").replace("ZH", "VVH").replace("WH", "VVH")
                    assert newname.count("VVH") == 1 and newname.startswith("VVH"), newname
                if "ggH" in k or "ttH" in k or "bbH" in k:
                    newname = k.replace("ggH", "ffH").replace("ttH", "ffH").replace("bbH", "ffH")
                    assert newname.count("ffH") == 1 and newname.startswith("ffH"), newname

            if newname is not None:
                if newname in cache:
                    cache[newname].Add(v)
                else:
                    cache[newname] = v.Clone(newname)
                    cache[newname].SetDirectory(f)

                if "_3D" not in newname:
                    self.histogramintegrals[newname] = cache[newname].Integral()

        nominalnames = set(t.GetName() for t in cache.itervalues() if "binbybin" not in t.GetName())
        assert len(nominalnames) + 2*len(self.binbybinuncertainties) == len(cache)

        expectedhistnames = set(self.allhistograms) | {"data_obs"}
        expectedhistnames |= {_+"_3D" for _ in expectedhistnames}
        if self.analysis.mergeallVVHandallffH:
            expectedhistnames |= {_.replace("qqH", "VVH").replace("ggH", "ffH") for _ in expectedhistnames}
        assert nominalnames == expectedhistnames, (nominalnames, expectedhistnames, nominalnames ^ expectedhistnames)

        f.Write()
        f.Close()

    def makeCardsWorkspaces(self, outdir="."):
        mkdir_p(outdir)
        with cd(outdir):
            if self.analysis.usehistogramsforcombine:
                Datacard(self.channel, self.category, self.analysis, self.luminosity).makehistograms()
            else:    
                for channel, category in product(channels, categories):
                    if self.production.LHE and channel != "2e2mu": continue
                    if self.analysis.isdecayonly and category != "Untagged": continue
                    Datacard(channel, category, self.analysis, self.luminosity).makepdfs()
                self.writeworkspace()
                self.linkworkspace()
            self.writedatacard()

class PdfBase(MultiEnum):
    enums = (_Datacard, ProductionMode, WorkspaceShapeSystematic, SystematicDirection)

    def __init__(self, *args):
        super(PdfBase, self).__init__(*args)
        if self.workspaceshapesystematic is None:
          self.shapesystematic = ShapeSystematic("")
        else:
          self.shapesystematic = ShapeSystematic(self.workspaceshapesystematic.nickname + str(self.systematicdirection))

    def check(self, *args):
        dontcheck = []
        if self.workspaceshapesystematic is None:
            dontcheck.append(WorkspaceShapeSystematic)
            dontcheck.append(SystematicDirection)
        super(PdfBase, self).check(self, *args, dontcheck=dontcheck)

class _Pdf(PdfBase):
    def __init__(self, *args, **kwargs):
        for thing in "fai", "a1", "ai", "D1", "D2", "D3", "faj", "phiai", "phiaj":
            assert not hasattr(self, thing)
            setattr(self, thing, kwargs[thing])
            del kwargs[thing]
        super(_Pdf, self).__init__(*args, **kwargs)
        self.__pdf = None
        self.__norm = None

    @property
    def pdf(self):
        return self.getpdf()
    @property
    def norm(self):
        if self.shapesystematic != "" or self.productionmode.isbkg and self.productionmode not in ("ggZZ", "ZX"):
            raise ValueError("Can't get norm for systematic or bkg pdf!\n{!r}".format(self))
        self.makepdf()
        return self.__norm
    @property
    def muscaled(self):
        return self.getmuscaled()

    @property
    def appendname(self):
        nameparts = [self.productionmode, self.production.year, self.channel, self.category, self.shapesystematic]
        nameparts = [str(_) for _ in nameparts]
        return "_".join(_ for _ in nameparts if _)

    def templatename(self, i=None):
        if i is None:
            return makename("T_{}".format(self.appendname))
        else:
            return makename("T_{}_{}".format(self.appendname, i))
    def integralname(self, i):
        return makename("integral_T_{}_{}".format(self.appendname, i))
    def datahistname(self, i=None):
        if i is None:
            return makename("T_datahist_{}".format(self.appendname))
        else:
            return makename("T_datahist_{}_{}".format(self.appendname, i))
    def histfuncname(self, i):
        return makename("T_histfunc_{}_{}".format(self.appendname, i))
    def ZXsinglepdfname(self, i):
        return makename("ZX_{}_{}".format(self.appendname, i))
    @property
    def pdftmpname(self):
        return makename("pdf_tmp_{}".format(self.appendname))
    @property
    def pdfname(self):
        result = self.productionmode.combinename
        if self.workspaceshapesystematic is not None:
            result += "_" + str(self.workspaceshapesystematic)
            if self.workspaceshapesystematic.isperchannel:
                result += str(self.channel)
            result += str(self.systematicdirection)
        return result

    @property
    def puresnormname(self):
        return makename("puresNorm_{}".format(self.appendname))
    @property
    def realintsnormname(self):
        return makename("realIntsNorm_{}".format(self.appendname))
    @property
    def individualnormname(self):
        if self.shapesystematic != "":
            raise ValueError("Can't get norm for systematic pdf!\n{!r}".format(self))
        return makename("norm_{}".format(self.appendname))
    @property
    def normtmpname(self):
        return makename("norm_tmp_{}".format(self.appendname))
    @property
    def normname(self):
        if self.shapesystematic != "":
            raise ValueError("Can't get norm for systematic pdf!\n{!r}".format(self))
        return "{}_norm".format(self.productionmode.combinename)  #no makename!

    @classmethod
    def murationame(cls, productionmodes, year):
        if productionmodes == ("ggH", "ttH", "bbH"):
            return makename("mufratio_{}".format(year))
        if productionmodes == ("VBF", "ZH", "WH"):
            return makename("muVratio_{}".format(year))
    @classmethod
    def muscaledname(cls, productionmodes, year):
        if productionmodes == ("ggH", "ttH", "bbH"):
            return makename("muf_scaled_{}".format(year))
        if productionmodes == ("VBF", "ZH", "WH"):
            return makename("muV_scaled_{}".format(year))

    @classmethod
    @cache
    def alphaMorphBkg_perchannel(cls, channel):
        channel = Channel(channel)
        name = makename(("CMS_fake_{}" if config.usenewZXsystematics else "CMS_zz4l_smd_zjets_bkg_{}").format(channel))
        result = ROOT.RooRealVar(name,name,0,-20,20)
        result.setConstant(False)
        return result

    @classmethod
    @cache
    def alphaMorphBkg_percategory(cls, category):
        assert not config.usenewZXsystematics
        category = Category(category)
        name = makename("CMS_zz4l_smd_zjets_bkg_{}".format(category))
        result = ROOT.RooRealVar(name,name,0,-20,20)
        result.setConstant(False)
        return result

    @classmethod
    @cache
    def alphaMorphBkg_percategorychannel(cls, category, channel):
        assert not config.usenewZXsystematics
        category = Category(category)
        channel = Channel(channel)
        name = makename("CMS_zz4l_smd_zjets_bkg_{}_{}".format(category, channel))
        result = ROOT.RooRealVar(name,name,0,-20,20)
        result.setConstant(False)
        return result

    @classmethod
    def alphaMorphBkg(cls, channel, category):
        class _(MultiEnum): enums = (Channel, Category)
        _ = _(channel, category)
        channel, category = _.channel, _.category
        if category == "Untagged" or config.usenewZXsystematics:
            return cls.alphaMorphBkg_perchannel(channel)
        else:
            if config.mergeZXVBFVHsystematics:
                return cls.alphaMorphBkg_percategory(category)
            else:
                return cls.alphaMorphBkg_percategorychannel(category, channel)


    @classmethod
    @cache
    def RV(cls):
        result = ROOT.RooRealVar(makename("RV"), "RV", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def RF(cls):
        result = ROOT.RooRealVar(makename("RF"), "RF", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def R(cls):
        result = ROOT.RooRealVar(makename("R"), "R", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def RV_13TeV(cls):
        result = ROOT.RooRealVar(makename("RV_13TeV"), "RV_13TeV", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def RF_13TeV(cls):
        result = ROOT.RooRealVar(makename("RF_13TeV"), "RF_13TeV", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def R_13TeV(cls):
        result = ROOT.RooRealVar(makename("R_13TeV"), "R_13TeV", 1, 0, 400)
        result.setConstant()
        return result
    @classmethod
    @cache
    def muV(cls):
        return ROOT.RooFormulaVar(makename("muV"), "muV", "@0*@1*@2*@3", ROOT.RooArgList(cls.RV(), cls.RV_13TeV(), cls.R(), cls.R_13TeV()))
    @classmethod
    @cache
    def muf(cls):
        return ROOT.RooFormulaVar(makename("muf"), "muf", "@0*@1*@2*@3", ROOT.RooArgList(cls.RF(), cls.RF_13TeV(), cls.R(), cls.R_13TeV()))

    @cache
    def makepdf_decayonly(self):
        self.T = [
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                 ]
        for i, t in enumerate(self.T, start=1):
            t.SetName(self.templatename(i))

        self.T_datahist = [ROOT.FastHisto3D_f(t) for i, t in enumerate(self.T, start=1)]
        self.T_histfunc = [ROOT.FastHisto3DFunc_f(self.histfuncname(i), "", ROOT.RooArgList(self.D1,self.D2,self.D3), datahist) for i, datahist in enumerate(self.T_datahist, start=1)]

        if self.shapesystematic == "":
            self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
            for i, integral in enumerate(self.T_integral, start=1):
                print "{} T{}".format(self.productionmode, i), integral.getVal()

            self.r_fai_pures_norm = ROOT.RooFormulaVar(self.puresnormname, "", "( (1-abs(@0))*@1+abs(@0)*@2 )",ROOT.RooArgList(self.fai, self.T_integral[0], self.T_integral[1]))
            self.r_fai_realints_norm = ROOT.RooFormulaVar(self.realintsnormname, "", "(sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1)",ROOT.RooArgList(self.fai, self.T_integral[2]))
            self.individualnorm = ROOT.RooFormulaVar(self.individualnormname, "", "(abs(@2))>1 ? 0. : TMath::Max((@0+@1),0)", ROOT.RooArgList(self.r_fai_pures_norm, self.r_fai_realints_norm, self.fai))
            self.norm_SM = self.T_integral[0]
            self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "@0/@1 * @2", ROOT.RooArgList(self.individualnorm, self.norm_SM, self.muf()))

    @cache
    def makepdf_decayonly_2D(self):
        #https://github.com/hroskes/HiggsAnalysis-CombinedLimit/blob/e828bd35a5b27dda1ae31a00154053e3d2b6fe89/src/HZZ4L_RooSpinZeroPdf_2D_fast.cc#L89-L97
        self.T = [
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[2], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g11gj1", self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "gi1gj1", self.channel, self.shapesystematic),
                  zerotemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                  zerotemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                  zerotemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                 ]
        for i, t in enumerate(self.T, start=1):
            t.SetName(self.templatename(i))

        self.T_datahist = [ROOT.FastHisto3D_f(t) for i, t in enumerate(self.T, start=1)]
        self.T_histfunc = [ROOT.FastHisto3DFunc_f(self.histfuncname(i), "", ROOT.RooArgList(self.D1,self.D2,self.D3), datahist) for i, datahist in enumerate(self.T_datahist, start=1)]

        if self.shapesystematic == "":
            self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
            for i, integral in enumerate(self.T_integral, start=1):
                print "{} T{}".format(self.productionmode, i), integral.getVal()

            self.individualnorm = ROOT.RooFormulaVar(self.individualnormname, "", "(abs(@0)+abs(@1))>1 ? 0. : TMath::Max(0, "
                                                                                    "(1-abs(@0)-abs(@1)) * @2 + "
                                                                                    "abs(@0) * @3 + "
                                                                                    "abs(@1) * @4 + "
                                                                                    "(@0>0?1:-1)   * sqrt((1-abs(@0)-abs(@1))*abs(@0)) * @5 + "
                                                                                    "(@1>0?1:-1)   * sqrt((1-abs(@0)-abs(@1))*abs(@1)) * @6 + "
                                                                                    "(@0*@1>0?1:-1)* sqrt(abs(@0*@1)) * @7"
                                                                                  ")",
                                                     utilities.RooArgList(self.fai, self.faj, *self.T_integral))
            self.norm_SM = self.T_integral[0]
            self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "@0/@1 * @2", ROOT.RooArgList(self.individualnorm, self.norm_SM, self.muf()))

    @cache
    def makepdf_proddec(self):
        self.T = [
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g13gi1", self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g12gi2", self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi3", self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, self.shapesystematic),
                 ]
        for i, t in enumerate(self.T, start=1):
            t.SetName(self.templatename(i))

        self.T_datahist = [ROOT.FastHisto3D_f(t) for i, t in enumerate(self.T, start=1)]
        self.T_histfunc = [ROOT.FastHisto3DFunc_f(self.histfuncname(i), "", ROOT.RooArgList(self.D1,self.D2,self.D3), datahist) for i, datahist in enumerate(self.T_datahist, start=1)]

        if self.shapesystematic == "":
            self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
            for i, integral in enumerate(self.T_integral, start=1):
                print "{} T{}".format(self.productionmode, i), integral.getVal()

            formula = " + ".join("@0**{}*@1**{}*@{}".format(4-i, i, i+2) for i in range(5))
            self.individualnorm = ROOT.RooFormulaVar(self.individualnormname, "", formula, ROOT.RooArgList(self.a1, self.ai, *self.T_integral))
            self.norm_SM = self.T_integral[0]
            self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "@0/@1 * @2", ROOT.RooArgList(self.individualnorm, self.norm_SM, self.muV()))

    def getmuscaled(self):
        """
        The individualnorm, set by makepdf, models how the number of events changes
        as a function of fai.  Now we want to divide by the sum of the individualnorms
        so that constant muV corresponds to the same number of observed VBF+ZH+WH events
        for any fai, and constant muf corresponds to the same number of observed ggH+ttH+bbH
        events for any fai.
        """

        if self.productionmode in ("VBF", "ZH", "WH"):
            productionmodes = [ProductionMode(_) for _ in ("VBF", "ZH", "WH")]
            mu = self.muV()
        if self.productionmode in ("ggH", "ttH", "bbH"):
            productionmodes = [ProductionMode(_) for _ in ("ggH", "ttH", "bbH")]
            mu = self.muf()
        return self.makemuscaled(tuple(productionmodes), self.luminosity, self.analysis, mu)

    @classmethod
    @cache
    def makemuratio(cls, productionmodes, luminosity, analysis):
        pdfswithsamemu = [Pdf(Datacard(analysis, luminosity, ca, ch), p) for p in productionmodes
                                                                         for ca in categories
                                                                         for ch in channels]
        for pdf in pdfswithsamemu:
            pdf.makepdf()  #does nothing if pdf is already made

        numerator = " + ".join("@{}".format(i) for i, pdf in enumerate(pdfswithsamemu, start=len(pdfswithsamemu)))
        denominator = " + ".join("@{}".format(i) for i, pdf in enumerate(pdfswithsamemu))
        formula = "({}) / ({})".format(numerator, denominator)
        arglist = utilities.RooArgList(*[_.individualnorm for _ in pdfswithsamemu]
                                       +[_.norm_SM for _ in pdfswithsamemu])
        return ROOT.RooFormulaVar(cls.murationame(productionmodes, luminosity.production.year), "", formula, arglist)

    @classmethod
    @cache
    def makemuscaled(cls, productionmodes, luminosity, analysis, mu):
        muratio = cls.makemuratio(productionmodes, luminosity, analysis)
        return ROOT.RooFormulaVar(cls.muscaledname(productionmodes, luminosity.production.year), "", "@0/@1", ROOT.RooArgList(mu, muratio))

    @classmethod
    @cache
    def one(cls):
        return ROOT.RooConstVar(makename("one"), "one", 1)

    @classmethod
    @cache
    def ZXnormstuff(cls, channel, category, year, nominal, up, dn):
        return (
            cls.one(),
            ROOT.RooConstVar(makename("zxup_{}_{}_{}".format(channel, category, year)), "", up / nominal),
            ROOT.RooConstVar(makename("zxdn_{}_{}_{}".format(channel, category, year)), "", dn / nominal),
        )

    @cache
    def makepdf_ZX(self):
        if hasattr(self, "T"): return
        if self.shapesystematic != "":
            raise ValueError("Do not give shape systematics to Z+X pdf!  They are handled internally.\n{!r}".format(self))

        shapesystematics = ShapeSystematic.items(lambda x: x in ("", "ZXUp", "ZXDn"))

        self.T = {shapesystematic: gettemplate(self.productionmode, self.analysis, self.production, self.category, self.channel, shapesystematic) for shapesystematic in shapesystematics}
        for shapesystematic, t in self.T.iteritems():
            t.SetName(self.templatename(shapesystematic))
        self.T_datahist = {shapesystematic: ROOT.RooDataHist(self.datahistname(shapesystematic), "", ROOT.RooArgList(self.D1,self.D2,self.D3), t) for shapesystematic, t in self.T.iteritems()}
        self.ZXpdfs = {shapesystematic: ROOT.RooHistPdf(self.ZXsinglepdfname(shapesystematic), "", ROOT.RooArgSet(self.D1,self.D2,self.D3), datahist) for shapesystematic, datahist in self.T_datahist.iteritems()}

        self.funcList_zjets = ROOT.RooArgList()
        self.morphVarListBkg = ROOT.RooArgList()

        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("")])
        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXUp")])
        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXDn")])
        self.morphVarListBkg.add(self.alphaMorphBkg(self.category, self.channel))

        self.__norm = ROOT.AsymQuad(self.normname, self.normname, 
          ROOT.RooArgList(*self.ZXnormstuff(self.channel, self.category, self.production.year, self.T[ShapeSystematic("")].Integral(), self.T[ShapeSystematic("ZXUp")].Integral(), self.T[ShapeSystematic("ZXDn")].Integral())),
          ROOT.RooArgList(self.alphaMorphBkg(self.category, self.channel)),
          1, 2
        )

    @classmethod
    @cache
    def kbkg_gg(cls):
        k = ROOT.RooRealVar(makename("kbkg_gg"), "kbkg_gg", 1, 0, 2)
        k.setBins(200)
        return k

    @cache
    def makepdf_bkg(self):
        if hasattr(self, "T"): return
        self.T = gettemplate(self.productionmode, self.analysis, self.production, self.category, self.channel, self.shapesystematic)
        self.T.SetName(self.templatename())
        self.T_datahist = ROOT.RooDataHist(self.datahistname(), "", ROOT.RooArgList(self.D1,self.D2,self.D3), self.T)
        self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "abs(@0)", ROOT.RooArgList(self.kbkg_gg()))

    def makepdf(self):
        if self.productionmode in ("ggH", "ttH", "bbH"):
            if self.analysis.dimensions == 2:
                self.makepdf_decayonly_2D()
            else:
                self.makepdf_decayonly()
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.analysis.dimensions == 2:
                self.makepdf_proddec_2D()
            else:
                self.makepdf_proddec()
        elif self.productionmode == "ZX":
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category in ("VBFtagged" ,"VHHadrtagged"):
                self.makepdf_ZX()
            else:
                self.makepdf_bkg()
        elif self.productionmode in ("ggZZ", "qqZZ", "VBFbkg"):
            self.makepdf_bkg()
        else:
            raise ValueError("Unknown productionmode {}".format(self.productionmode))

    @cache
    def getpdf_decayonly(self):
        return ROOT.HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(self.pdfname, self.pdfname, self.a1, self.ai, ROOT.RooArgList(self.D1, self.D2, self.D3), ROOT.RooArgList(*self.T_histfunc))
    @cache
    def getpdf_decayonly_2D(self):
        return ROOT.HZZ4L_RooSpinZeroPdf_2D_fast(self.pdfname, self.pdfname, self.fai, self.faj, self.phiai, self.phiaj, ROOT.RooArgList(self.D1, self.D2, self.D3), ROOT.RooArgList(*self.T_histfunc))
    @cache
    def getpdf_proddec(self):
        return ROOT.VVHZZ4L_RooSpinZeroPdf_1D_fast_a1ai(self.pdfname, self.pdfname, self.a1, self.ai, ROOT.RooArgList(self.D1, self.D2, self.D3), ROOT.RooArgList(*self.T_histfunc))
    @cache
    def getpdf_ZX(self):
        return ROOT.FastVerticalInterpHistPdf3D(self.pdfname,self.pdfname,self.D1,self.D2,self.D3,False,self.funcList_zjets,self.morphVarListBkg,1.0,1)
    @cache
    def getpdf_bkg(self):
        return ROOT.RooHistPdf(self.pdfname, self.pdfname, ROOT.RooArgSet(self.D1,self.D2,self.D3), self.T_datahist)

    def getpdf(self):
        self.makepdf()
        if self.productionmode in ("ggH", "ttH", "bbH"):
            if self.analysis.dimensions == 2:
                return self.getpdf_decayonly_2D()
            else:
                return self.getpdf_decayonly()
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.analysis.dimensions == 2:
                return self.getpdf_proddec_2D()
            else:
                return self.getpdf_proddec()
        elif self.productionmode == "ZX":
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category in ("VBFtagged" ,"VHHadrtagged"):
                return self.getpdf_ZX()
            else:
                return self.getpdf_bkg()
        elif self.productionmode in ("ggZZ", "qqZZ", "VBFbkg"):
            return self.getpdf_bkg()
        else:
            raise ValueError("Unknown productionmode {}".format(self.productionmode))

Datacard = multienumcache(_Datacard)
Pdf = multienumcache(_Pdf, haskwargs=True, multienumforkey=PdfBase)

def makeDCsandWSs(productions, channels, categories, *otherargs, **kwargs):
    with OneAtATime("makeDCsandWSs.tmp", 30):
        dcs = [
               Datacard(production, channel, category, *otherargs)
                       for production, channel, category
                       in product(productions, channels, categories)
              ]
        if all(os.path.exists(thing) for dc in dcs for thing in (dc.rootfile_base, dc.rootfile, dc.txtfile)):
            return
        for dc in dcs:
            if dc.production.LHE and dc.channel != "2e2mu": continue
            if dc.analysis.isdecayonly and dc.category != "Untagged": continue
            dc.makeCardsWorkspaces(**kwargs)
            for thing in dc.rootfile_base, dc.rootfile, dc.txtfile:
                if not os.path.exists(thing):
                    raise ValueError("{} was not created.  Something is wrong.".format(thing))
