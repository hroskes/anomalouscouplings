from collections import Counter
import inspect
import itertools
import os
import time

import ROOT

import combineinclude
from combinehelpers import discriminants, getdatatree, gettemplate, getnobserved, getrate, Luminosity, mixturesign, sigmaioversigma1
import config
from enums import Analysis, categories, Category, Channel, channels, MultiEnum, Production, ProductionMode, ShapeSystematic, SystematicDirection, WorkspaceShapeSystematic
import utilities
from utilities import cache, callclassinitfunctions, cd, generatortolist, mkdir_p, multienumcache, OneAtATime, Tee
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
        result = []
        for label in self.labels:
            if label.startswith("#"):
                result.append(label)
            else:
                value = getattr(obj, label)
                if value is None: continue
                result.append("{} {}".format(label, value))
        return result

class SystematicsSection(Section):
    def getlines(self, obj, objtype):
        result = super(SystematicsSection, self).getlines(obj, objtype)
        for line in result[:]:
            if len(line.split()) > 2 and all(systematicvalue == "-" for systematicvalue in line.split()[2:]):
                result.remove(line)
        return result

class SystematicFromEnums_BaseClass(object):
    pass

def MakeSystematicFromEnums(*theenums, **kwargs):
    class SystematicFromEnums(SystematicFromEnums_BaseClass):
        def __init__(self, function):
            self.name = self.origname = function.__name__
            for enum in theenums:
                self.name = self.name.replace(enum.enumname, "{"+enum.enumname+"}")

            for kw, kwarg in kwargs.iteritems():
                if kw == "name":
                    self.name = kwarg
                else:
                    raise ValueError("Unknown kwarg {}={}".format(kw, kwarg))

            class SystematicFromEnumsWithValues(MultiEnum):
                enums = theenums

                def __get__(self_systematic, self_datacard, objtype):
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
            cartesianproduct = itertools.product(*(enum.items() for enum in theenums))
            for enumvalues in cartesianproduct:
                formatdict = {enum.enumname: enumvalue for enum, enumvalue in zip(theenums, enumvalues)}
                yield self.name.format(**formatdict)


        def applyallfunctions(self, datacardcls):
            cartesianproduct = itertools.product(*(enum.items() for enum in theenums))
            for enumvalues in cartesianproduct:
                formatdict = {enum.enumname: enumvalue for enum, enumvalue in zip(theenums, enumvalues)}
                name = self.name.format(**formatdict)
                setattr(datacardcls, name, self.systematicfromenumswithvalues(*enumvalues))

            delattr(datacardcls, self.origname)

    return SystematicFromEnums

@callclassinitfunctions("initsystematicsfromenums")
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
        return "hzz4l_{}S_{}_{}.lumi{}.txt".format(self.channel, self.category, self.year, float(self.luminosity))
    @property
    def rootfile(self):
        return "hzz4l_{}S_{}_{}.lumi{}.input.root".format(self.channel, self.category, self.year, float(self.luminosity))
    @property
    def rootfile_base(self):
        return "hzz4l_{}S_{}_{}.input.root".format(self.channel, self.category, self.year)
    @property
    def logfile(self):
        return "log_hzz4l_{}S_{}_{}.input".format(self.channel, self.category, self.year)

    @property
    def productionmodes(self):
        return [ProductionMode(p) for p in ("ggH", "qqH", "WH", "ZH", "ttH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets")]

    @property
    def imax(self):
        return 1
    @property
    def jmax(self):
        return len(self.productionmodes)-1
    @property
    def kmax(self):
        return "*"

    section1 = Section("imax", "jmax", "kmax")

    @property
    def shapes(self):
        return "* * {} w:$PROCESS w:$PROCESS_$SYSTEMATIC".format(self.rootfile)

    section2 = Section("shapes")

    @property
    def bin(self, counter=Counter()):
        counter[self] += 1

        bin = "a{}".format(len(channels) * channels.index(self.channel) + categories.index(self.category))

        if counter[self] == 1:
            return bin
        elif counter[self] == 2:
            return " ".join([str(bin)]*len(self.productionmodes))
        assert False

    @property
    def observation(self):
        return getnobserved(self.channel, self.production, self.category, self.analysis)

    section3 = Section("bin", "observation")

    @property
    def process(self, counter=Counter()):
        counter[self] += 1

        if counter[self] == 1:
            return " ".join(_.combinename for _ in self.productionmodes)

        if counter[self] == 2:
            nsignal = sum(_.issignal for _ in self.productionmodes)
            nbkg = sum(_.isbkg for _ in self.productionmodes)
            assert nsignal+nbkg == len(self.productionmodes)

            return " ".join(str(_) for _ in range(-nsignal+1, nbkg+1))

        assert False

    @property
    def rate(self):
        return " ".join(str(getrate(p, self.channel, self.category, self.analysis, self.luminosity)) for p in self.productionmodes)

    section4 = Section("## mass window [{},{}]".format(config.m4lmin, config.m4lmax),
                       "bin", "process", "process", "rate")

    @MakeSystematicFromEnums(YieldSystematic)
    def yieldsystematic(self, yieldsystematic):
        return " ".join(
                        ["lnN"] +
                        [str(YieldSystematicValue(yieldsystematic, self.channel, self.category, self.analysis, p))
                            for p in self.productionmodes]
                       )

    @MakeSystematicFromEnums(WorkspaceShapeSystematic, Channel)
    def workspaceshapesystematicchannel(self, workspaceshapesystematic, channel):
      if workspaceshapesystematic.isperchannel and channel == self.channel:
        return " ".join(
                        ["shape1"] +
                        ["1" if workspaceshapesystematic in p.workspaceshapesystematics(self.category) else "-"
                            for p in self.productionmodes]
                       )

    @MakeSystematicFromEnums(WorkspaceShapeSystematic)
    def workspaceshapesystematic(self, workspaceshapesystematic):
      if workspaceshapesystematic.isperchannel: return None
      return " ".join(
                      ["shape1"] +
                      ["1" if workspaceshapesystematic in p.workspaceshapesystematics(self.category) else "-"
                          for p in self.productionmodes]
                     )

    @MakeSystematicFromEnums(Channel, Category)
    def CMS_zz4l_smd_zjets_bkg_channel_category(self, channel, category):
        if channel == self.channel and category == self.category:
            if config.applyZXshapesystematicsUntagged and category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and category != "Untagged":
                return "param 0 1 [-3,3]"

    @property
    def muVratio(self): return "extArg {}:w:RecycleConflictNodes".format(self.rootfile)
    @property
    def mufratio(self): return "extArg {}:w:RecycleConflictNodes".format(self.rootfile)

    section5 = SystematicsSection(yieldsystematic, workspaceshapesystematicchannel, workspaceshapesystematic, CMS_zz4l_smd_zjets_bkg_channel_category, "muVratio", "mufratio")

    divider = "\n------------\n"

    @classmethod
    def initsystematicsfromenums(cls):
        for name, systematic in inspect.getmembers(cls, predicate=lambda x: isinstance(x, SystematicFromEnums_BaseClass)):
            systematic.applyallfunctions(cls)

    def writedatacard(self):
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
    def mixturesign_constvar(cls, analysis):
        name = makename("mixturesign_{}".format(analysis))
        return ROOT.RooConstVar(name, name, mixturesign(analysis))
    @classmethod
    @cache
    def sigmaioversigma1_constvar(cls, analysis):
        name = makename("sigmaioversigma1_{}".format(analysis))
        return ROOT.RooConstVar(name, name, sigmaioversigma1(analysis, "ggH"))
    @classmethod
    @cache
    def a1(cls):
        name = makename("a1")
        return ROOT.RooFormulaVar(name, name, "sqrt(1-abs(@0))", ROOT.RooArgList(cls.fai()))
    @classmethod
    @cache
    def ai(cls, analysis):
        name = makename("ai")
        return ROOT.RooFormulaVar(name, name, "@2 * (@0>0 ? 1 : -1) * sqrt(abs(@0)/@1)", ROOT.RooArgList(cls.fai(), cls.sigmaioversigma1_constvar(analysis), cls.mixturesign_constvar(analysis)))

    def makepdfs(self):
        if self.pdfs is not None: return

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##

        fai = self.fai()
        mixturesign_constvar = self.mixturesign_constvar(self.analysis)
        sigmaioversigma1_constvar = self.sigmaioversigma1_constvar(self.analysis)
        a1 = self.a1()
        ai = self.ai(self.analysis)

        #add category name in case the same discriminant is used in multiple categories
        discs = discriminants(self.analysis, self.category)
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
            pdfkwargs = {"fai": fai, "a1": a1, "ai": ai, "D1": D1, "D2": D2, "D3": D3}
            self.pdfs.append(Pdf(self, p, **pdfkwargs))
            for systematic in p.workspaceshapesystematics(self.category):
                self.pdfs.append(Pdf(self, p, systematic, "Up", **pdfkwargs))
                self.pdfs.append(Pdf(self, p, systematic, "Down", **pdfkwargs))

        ## --------------------------- DATASET --------------------------- ##

        data_obs_tree = getdatatree(self.channel, self.production, self.category, self.analysis)
        datasetName = makename("data_obs_{}_{}".format(self.channel, self.category))


        self.data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(D1,D2,D3))


        ## --------------------------- WORKSPACE -------------------------- ##

    def writeworkspace(self):
        if os.path.exists(self.rootfile_base): return
        w = ROOT.RooWorkspace("w","w")

        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf.Class(),True)
        w.importClassCode(ROOT.VBFHZZ4L_RooSpinZeroPdf.Class(),True)

        getattr(w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        for pdf in self.pdfs:
            getattr(w, 'import')(pdf.pdf, ROOT.RooFit.RecycleConflictNodes())
            if pdf.productionmode.issignal and pdf.shapesystematic == "":
                getattr(w, 'import')(pdf.norm, ROOT.RooFit.RecycleConflictNodes())
                getattr(w, 'import')(pdf.muratio, ROOT.RooFit.RecycleConflictNodes())

        w.writeToFile(self.rootfile_base)

    def linkworkspace(self):
        if not os.path.exists(self.rootfile):
            os.symlink(self.rootfile_base, self.rootfile)

    def makeCardsWorkspaces(self, outdir="."):
        mkdir_p(outdir)
        with cd(outdir), Tee(self.logfile, 'w'):
            for channel, category in itertools.product(channels, categories):
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
          self.shapesystematic = ShapeSystematic(str(self.workspaceshapesystematic) + str(self.systematicdirection))

    def check(self, *args):
        dontcheck = []
        if self.workspaceshapesystematic is None:
            dontcheck.append(WorkspaceShapeSystematic)
            dontcheck.append(SystematicDirection)
        super(PdfBase, self).check(self, *args, dontcheck=dontcheck)

class _Pdf(PdfBase):
    def __init__(self, *args, **kwargs):
        for thing in "fai", "a1", "ai", "D1", "D2", "D3":
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
        if self.shapesystematic != "" or self.productionmode.isbkg:
            raise ValueError("Can't get norm for systematic or bkg pdf!\n{!r}".format(self))
        self.makepdf()
        return self.__norm
    @property
    def muratio(self):
        return self.getmuratio()

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
    def murationame(cls, productionmodes):
        if productionmodes == ("ggH", "ttH"):
            return makename("mufratio")
        if productionmodes == ("VBF", "ZH", "WH"):
            return makename("muVratio")

    @cache
    def makepdf_decayonly(self):
        self.T = [
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, self.shapesystematic),
                  gettemplate(self.productionmode, self.analysis, self.production, self.category, "g11gi1", self.channel, self.shapesystematic),
                 ]
        for i, t in enumerate(self.T, start=1):
            t.SetName(self.templatename(i))

        self.T_datahist = [ROOT.RooDataHist(self.datahistname(i), "", ROOT.RooArgList(self.D1,self.D2,self.D3), t) for i, t in enumerate(self.T, start=1)]
        self.T_histfunc = [ROOT.RooHistFunc(self.histfuncname(i), "", ROOT.RooArgSet(self.D1,self.D2,self.D3), datahist) for i, datahist in enumerate(self.T_datahist, start=1)]

        if self.shapesystematic == "":
            self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
            for i, integral in enumerate(self.T_integral, start=1):
                print "{} T{}".format(self.productionmode, i), integral.getVal()

            self.r_fai_pures_norm = ROOT.RooFormulaVar(self.puresnormname, "", "( (1-abs(@0))*@1+abs(@0)*@2 )",ROOT.RooArgList(self.fai, self.T_integral[0], self.T_integral[1]))
            self.r_fai_realints_norm = ROOT.RooFormulaVar(self.realintsnormname, "", "(sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1)",ROOT.RooArgList(self.fai, self.T_integral[2]))
            self.individualnorm = ROOT.RooFormulaVar(self.individualnormname, "", "(abs(@2))>1 ? 0. : TMath::Max((@0+@1),0)", ROOT.RooArgList(self.r_fai_pures_norm, self.r_fai_realints_norm, self.fai))
            self.norm_SM = self.T_integral[0]
            self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "@0/@1", ROOT.RooArgList(self.individualnorm, self.norm_SM))

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

        self.T_datahist = [ROOT.RooDataHist(self.datahistname(i), "", ROOT.RooArgList(self.D1,self.D2,self.D3), t) for i, t in enumerate(self.T, start=1)]
        self.T_histfunc = [ROOT.RooHistFunc(self.histfuncname(i), "", ROOT.RooArgSet(self.D1,self.D2,self.D3), datahist) for i, datahist in enumerate(self.T_datahist, start=1)]

        if self.shapesystematic == "":
            self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
            for i, integral in enumerate(self.T_integral, start=1):
                print "{} T{}".format(self.productionmode, i), integral.getVal()

            formula = " + ".join("@0**{}*@1**{}*@{}".format(4-i, i, i+2) for i in range(5))
            self.individualnorm = ROOT.RooFormulaVar(self.individualnormname, "", formula, ROOT.RooArgList(self.a1, self.ai, *self.T_integral))
            self.norm_SM = self.T_integral[0]
            self.__norm = ROOT.RooFormulaVar(self.normname, self.normname, "@0/@1", ROOT.RooArgList(self.individualnorm, self.norm_SM))

    def getmuratio(self):
        """
        The individualnorm, set by makepdf, models how the number of events changes
        as a function of fai.  Now we want to divide by the sum of the individualnorms
        so that constant muV corresponds to the same number of observed VBF+ZH+WH events
        for any fai, and constant muf corresponds to the same number of observed ggH+ttH
        events for any fai.
        """

        if self.productionmode in ("VBF", "ZH", "WH"):
            productionmodes = [ProductionMode(_) for _ in ("VBF", "ZH", "WH")]
        if self.productionmode in ("ggH", "ttH"):
            productionmodes = [ProductionMode(_) for _ in ("ggH", "ttH")]
        return self.makemuratio(tuple(productionmodes), self.luminosity, self.analysis)

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
        return ROOT.RooFormulaVar(cls.murationame(productionmodes), "", formula, arglist)

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
        self.morphBkgVarName =  makename("CMS_zz4l_smd_zjets_bkg_{}_{}".format(self.channel, self.category))
        self.alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        self.morphVarListBkg = ROOT.RooArgList()

        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("")])
        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXUp")])
        self.funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXDn")])
        self.alphaMorphBkg.setConstant(False)
        self.morphVarListBkg.add(self.alphaMorphBkg)

    @cache
    def makepdf_bkg(self):
        if hasattr(self, "T"): return
        self.T = gettemplate(self.productionmode, self.analysis, self.production, self.category, self.channel, self.shapesystematic)
        self.T.SetName(self.templatename())
        self.T_datahist = ROOT.RooDataHist(self.datahistname(), "", ROOT.RooArgList(self.D1,self.D2,self.D3), self.T)

    def makepdf(self):
        if self.productionmode in ("ggH", "ttH"):
            self.makepdf_decayonly()
        elif self.productionmode in ("VBF", "ZH", "WH"):
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
        return ROOT.HZZ4L_RooSpinZeroPdf(self.pdfname, self.pdfname, self.D1, self.D2, self.D3, self.fai, ROOT.RooArgList(*self.T_histfunc))
    @cache
    def getpdf_proddec(self):
        return ROOT.VBFHZZ4L_RooSpinZeroPdf(self.pdfname, self.pdfname, self.D1, self.D2, self.D3, self.a1, self.ai, ROOT.RooArgList(*self.T_histfunc))
    @cache
    def getpdf_ZX(self):
        return ROOT.FastVerticalInterpHistPdf3D(self.pdfname,self.pdfname,self.D1,self.D2,self.D3,False,self.funcList_zjets,self.morphVarListBkg,1.0,1)
    @cache
    def getpdf_bkg(self):
        return ROOT.RooHistPdf(self.pdfname, self.pdfname, ROOT.RooArgSet(self.D1,self.D2,self.D3), self.T_datahist)

    def getpdf(self):
        self.makepdf()
        if self.productionmode in ("ggH", "ttH"):
            return self.getpdf_decayonly()
        elif self.productionmode in ("VBF", "ZH", "WH"):
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
        for production, channel, category in itertools.product(productions, channels, categories):
            dc = Datacard(production, channel, category, *otherargs)
            dc.makeCardsWorkspaces(**kwargs)
            for thing in dc.rootfile_base, dc.rootfile, dc.txtfile:
                if not os.path.exists(thing):
                    raise ValueError("{} was not created.  Something is wrong.".format(thing))
