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
from utilities import callclassinitfunctions, cd, generatortolist, KeepWhileOpenFile, mkdir_p, Tee
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
            if all(systematicvalue == "-" for systematicvalue in line.split()[2:]):
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
class Datacard(MultiEnum):
    enums = (Analysis, Category, Channel, Luminosity)
    enumname = "datacard"
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

    section5 = SystematicsSection(yieldsystematic, workspaceshapesystematicchannel, workspaceshapesystematic, CMS_zz4l_smd_zjets_bkg_channel_category)

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

    def writeworkspace(self):
        if os.path.exists(self.rootfile_base): return

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##

        bins = 1000

        fai_name = "CMS_zz4l_fai1"

        fai = ROOT.RooRealVar(fai_name,fai_name,-1.,1.)
        mixturesign_constvar = ROOT.RooConstVar("mixturesign", "mixturesign", mixturesign(self.analysis))
        sigmaioversigma1_constvar = ROOT.RooConstVar("sigmaioversigma1", "sigmaioversigma1", sigmaioversigma1(self.analysis, "ggH"))
        a1 = ROOT.RooFormulaVar("a1", "a1", "sqrt(1-abs(@0))", ROOT.RooArgList(fai))
        ai = ROOT.RooFormulaVar("ai", "ai", "@2 * (@0>0 ? 1 : -1) * sqrt(abs(@0)/@1)", ROOT.RooArgList(fai, sigmaioversigma1_constvar, mixturesign_constvar))
        fai.setBins(bins)

        #add category name in case the same discriminant is used in multiple categories
        discs = discriminants(self.analysis, self.category)
        D1Name, D2Name, D3Name = ("{}_{}".format(d.name, self.category) for d in discs)
        dBinsX, dBinsY, dBinsZ = (d.bins for d in discs)
        dLowX, dLowY, dLowZ = (d.min for d in discs)
        dHighX, dHighY, dHighZ = (d.max for d in discs)

        D1 = ROOT.RooRealVar(D1Name, D1Name, dLowX, dHighX)
        D2 = ROOT.RooRealVar(D2Name, D2Name, dLowY, dHighY)
        D3 = ROOT.RooRealVar(D3Name, D3Name, dLowZ, dHighZ)
        D1.setBins(dBinsX)
        D2.setBins(dBinsY)
        D3.setBins(dBinsZ)

        pdfs = []

        for p in self.productionmodes:
            pdfkwargs = {"fai": fai, "a1": a1, "ai": ai, "D1": D1, "D2": D2, "D3": D3}
            pdfs.append(Pdf(self, p, **pdfkwargs))
            for systematic in p.workspaceshapesystematics(self.category):
                pdfs.append(Pdf(self, p, systematic, "Up", **pdfkwargs))
                pdfs.append(Pdf(self, p, systematic, "Down", **pdfkwargs))

        ## --------------------------- DATASET --------------------------- ##

        data_obs_tree = getdatatree(self.channel, self.production, self.category, self.analysis)
        data_obs = ROOT.RooDataSet()
        datasetName = makename("data_obs_{}_{}".format(self.channel, self.category))


        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(D1,D2,D3))


        ## --------------------------- WORKSPACE -------------------------- ##

        w = ROOT.RooWorkspace("w","w")

        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf.Class(),True)
        w.importClassCode(ROOT.VBFHZZ4L_RooSpinZeroPdf.Class(),True)

        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        for pdf in pdfs:
            getattr(w, 'import')(pdf.pdf, ROOT.RooFit.RecycleConflictNodes())
            if pdf.productionmode.issignal and pdf.shapesystematic == "":
                getattr(w, 'import')(pdf.norm, ROOT.RooFit.RecycleConflictNodes())

        w.writeToFile(self.rootfile_base)

    def makeworkspace(self):
        self.writeworkspace()
        if not os.path.exists(self.rootfile):
            os.symlink(self.rootfile_base, self.rootfile)

    def makeCardsWorkspaces(self, outdir="."):
        mkdir_p(outdir)
        with cd(outdir), Tee(self.logfile, 'w'):
            self.makeworkspace()
            self.writedatacard()

class Pdf(MultiEnum):
    enums = (Datacard, ProductionMode, WorkspaceShapeSystematic, SystematicDirection)
    def __init__(self, *args, **kwargs):
        for thing in "fai", "a1", "ai", "D1", "D2", "D3":
            setattr(self, thing, kwargs[thing])
            del kwargs[thing]
        super(Pdf, self).__init__(*args, **kwargs)
        self.__pdf = None
        self.__norm = None
        if self.workspaceshapesystematic is None:
          self.shapesystematic = ShapeSystematic("")
        else:
          self.shapesystematic = ShapeSystematic(str(self.workspaceshapesystematic) + str(self.systematicdirection))

    def check(self, *args):
        dontcheck = []
        if self.workspaceshapesystematic is None:
            dontcheck.append(WorkspaceShapeSystematic)
            dontcheck.append(SystematicDirection)
        super(Pdf, self).check(self, *args, dontcheck=dontcheck)

    @property
    def pdf(self):
        if self.__pdf is None: self.makepdf()
        return self.__pdf
    @property
    def norm(self):
        if self.shapesystematic != "":
            raise ValueError("Can't get norm for systematic pdf!\n{!r}".format(self))
        if self.__norm is None: self.makepdf()
        return self.__norm

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
    def normname(self):
        if self.shapesystematic != "":
            raise ValueError("Can't get norm for systematic pdf!\n{!r}".format(self))
        return "{}_norm".format(self.productionmode.combinename)  #no makename!

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
            self.__pdf = ROOT.HZZ4L_RooSpinZeroPdf(self.pdfname, self.pdfname, self.D1, self.D2, self.D3, self.fai, ROOT.RooArgList(*self.T_histfunc))

            if self.shapesystematic == "":
                self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
                for i, integral in enumerate(self.T_integral, start=1):
                    print "{} T{}".format(self.productionmode, i), integral.getVal()

                self.r_fai_pures_norm = ROOT.RooFormulaVar(self.puresnormname, "", "( (1-abs(@0))*@1+abs(@0)*@2 )/@1",ROOT.RooArgList(self.fai, self.T_integral[0], self.T_integral[1]))
                self.r_fai_realints_norm = ROOT.RooFormulaVar(self.realintsnormname, "", "(sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1)/@2",ROOT.RooArgList(self.fai, self.T_integral[2], self.T_integral[0]))
                self.__norm = ROOT.RooFormulaVar(self.normname, "", "(abs(@2))>1 ? 0. : TMath::Max((@0+@1),0)", ROOT.RooArgList(self.r_fai_pures_norm, self.r_fai_realints_norm, self.fai))

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
            self.__pdf = ROOT.VBFHZZ4L_RooSpinZeroPdf(self.pdfname, self.pdfname, self.D1, self.D2, self.D3, self.a1, self.ai, ROOT.RooArgList(*self.T_histfunc))

            if self.shapesystematic == "":
                self.T_integral = [ROOT.RooConstVar(self.integralname(i), "", t.Integral()) for i, t in enumerate(self.T, start=1)]
                for i, integral in enumerate(self.T_integral, start=1):
                    print "{} T{}".format(self.productionmode, i), integral.getVal()

                formula = " + ".join("@0**{}*@1**{}*@{}".format(4-i, i, i+2) for i in range(5))
                formula = "("+formula+") / @2"
                self.__norm = ROOT.RooFormulaVar(self.normname, formula, ROOT.RooArgList(self.a1, self.ai, *self.T_integral))

    def makepdf_ZX(self):
        if self.shapesystematic != "":
            raise ValueError("Do not give shape systematics to Z+X pdf!  They are handled internally.\n{!r}".format(self))

        shapesystematics = ShapeSystematic.items(lambda x: x in ("", "ZXUp", "ZXDn"))

        self.T = {shapesystematic: gettemplate(self.productionmode, self.analysis, self.production, self.category, self.channel, shapesystematic) for shapesystematic in shapesystematics}
        for shapesystematic, t in self.T.iteritems():
            t.SetName(self.templatename(shapesystematic))
        self.T_datahist = {shapesystematic: ROOT.RooDataHist(self.datahistname(shapesystematic), "", ROOT.RooArgList(self.D1,self.D2,self.D3), t) for shapesystematic, t in self.T.iteritems()}
        self.ZXpdfs = {shapesystematic: ROOT.RooHistPdf(self.ZXsinglepdfname(shapesystematic), "", ROOT.RooArgSet(self.D1,self.D2,self.D3), datahist) for shapesystematic, datahist in self.T_datahist.iteritems()}

        funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  makename("CMS_zz4l_smd_zjets_bkg_{}_{}".format(self.channel, self.category))
        self.alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()

        funcList_zjets.add(self.ZXpdfs[ShapeSystematic("")])
        funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXUp")])
        funcList_zjets.add(self.ZXpdfs[ShapeSystematic("ZXDn")])
        self.alphaMorphBkg.setConstant(False)
        morphVarListBkg.add(self.alphaMorphBkg)

        self.__pdf = ROOT.FastVerticalInterpHistPdf3D(self.pdfname,self.pdfname,self.D1,self.D2,self.D3,False,funcList_zjets,morphVarListBkg,1.0,1)

    def makepdf_bkg(self):
        self.T = gettemplate(self.productionmode, self.analysis, self.production, self.category, self.channel, self.shapesystematic)
        self.T.SetName(self.templatename())
        self.T_datahist = ROOT.RooDataHist(self.datahistname(), "", ROOT.RooArgList(self.D1,self.D2,self.D3), self.T)
        self.__pdf = ROOT.RooHistPdf(self.pdfname, self.pdfname, ROOT.RooArgSet(self.D1,self.D2,self.D3), self.T_datahist)

    def makepdf(self):
        if self.productionmode in ("ggH", "ttH"):
            self.makepdf_decayonly()
        elif self.productionmode in ("VBF", "ZH", "WH"):
            self.makepdf_proddec()
        elif self.productionmode == "ZX":
            if config.applyZXshapesystematicsUntagged and self.category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and self.category != "Untagged":
                self.makepdf_ZX()
            else:
                self.makepdf_bkg()
        elif self.productionmode in ("ggZZ", "qqZZ", "VBFbkg"):
            self.makepdf_bkg()
        else:
            raise ValueError("Unknown productionmode {}".format(self.productionmode))

def makeDCsandWSs(productions, channels, categories, *otherargs, **kwargs):
    done = {(production, channel, category): False for production, channel, category in itertools.product(productions, channels, categories)}
    while not all(done.values()):
        anychanged = False
        for production, channel, category in itertools.product(productions, channels, categories):
            if done[production,channel,category]: continue
            dc = Datacard(production, channel, category, *otherargs)
            with KeepWhileOpenFile(dc.txtfile+".tmp") as f:
                if not f:
                    continue
                if os.path.exists(dc.rootfile_base) and os.path.exists(dc.rootfile) and os.path.exists(dc.txtfile):
                    done[production,channel,category] = True
                    continue
                dc.makeCardsWorkspaces(**kwargs)
                for thing in dc.rootfile_base, dc.rootfile, dc.txtfile:
                    if not os.path.exists(thing):
                        raise ValueError("{} was not created.  Something is wrong.".format(thing))
                anychanged = done[production,channel,category] = True
        if not anychanged and not all(done.values()):
            print "Some datacards are being created by other processes.  Waiting 30 seconds..."
            time.sleep(30)
