from collections import Counter
import os

import combineinclude
from combinehelpers import discriminants, getdatatree, gettemplate, getnobserved, getrate, Luminosity, mixturesign, sigmaioversigma1
import config
from enums import Analysis, categories, Category, Channel, channels, MultiEnum, Production, ProductionMode
from utilities import cd, mkdir_p

import ROOT

class Section(object):
    def __init__(self, *labels):
        self.labels = labels
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

class Datacard(MultiEnum):
    enums = (Analysis, Category, Channel, Production, Luminosity)
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
    def productionmodes(self):
        return [ProductionMode(p) for p in ("ggH", "qqH", "WH", "ZH", "bkg_qqzz", "bkg_ggzz", "bkg_vbf", "bkg_zjets")]

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

    #just one dummy systematic for now
    @property
    def lumi_13TeV_common(self):
        return " ".join(["lnN"] + ["1.023" for p in self.productionmodes])

    @property
    def CMS_zz4l_smd_zjets_bkg_2e2mu_Untagged(self):
        channel, category = "2e2mu", "Untagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_2e2mu_VBFtagged(self):
        channel, category = "2e2mu", "VBFtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_2e2mu_VHHadrtagged(self):
        channel, category = "2e2mu", "VHHadrtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4e_Untagged(self):
        channel, category = "4e", "Untagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4e_VBFtagged(self):
        channel, category = "4e", "VBFtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4e_VHHadrtagged(self):
        channel, category = "4e", "VHHadrtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4mu_Untagged(self):
        channel, category = "4mu", "Untagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4mu_VBFtagged(self):
        channel, category = "4mu", "VBFtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    @property
    def CMS_zz4l_smd_zjets_bkg_4mu_VHHadrtagged(self):
        channel, category = "4mu", "VHHadrtagged"
        if self.channel == channel and self.category == category:
            return "param 0 1 [-3,3]"
        return None

    section5 = SystematicsSection(*["lumi_13TeV_common"] + ["CMS_zz4l_smd_zjets_bkg_{}_{}".format(channel, category) for channel in channels for category in categories])

    divider = "\n------------\n"

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

        x_name = "CMS_zz4l_fai1"

        x = ROOT.RooRealVar(x_name,x_name,-1.,1.)
        mixturesign_constvar = ROOT.RooConstVar("mixturesign", "mixturesign", mixturesign(self.analysis))
        sigmaioversigma1_constvar = ROOT.RooConstVar("sigmaioversigma1", "sigmaioversigma1", sigmaioversigma1(self.analysis, "ggH"))
        a1 = ROOT.RooFormulaVar("a1", "a1", "sqrt(1-abs(@0))", ROOT.RooArgList(x))
        ai = ROOT.RooFormulaVar("ai", "ai", "@2 * (@0>0 ? 1 : -1) * sqrt(abs(@0)/@1)", ROOT.RooArgList(x, sigmaioversigma1_constvar, mixturesign_constvar))
        x.setBins(bins)

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

        T = {}
        T_ScaleResUp = {}
        T_ScaleResDown = {}
        T_integral = {}
        T_hist = {}
        T_ScaleResUp_hist = {}
        T_ScaleResDown_hist = {}
        T_ScaleResDown_hist = {}
        T_histfunc = {}
        T_ScaleResUp_histfunc = {}
        T_ScaleResDown_histfunc = {}
        pdf = {}
        pdf_syst1Up = {}
        pdf_syst1Down = {}
        norm = {}


        #for ggH, the order is SM, BSM, int
        T["ggH"] = [
                    gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel),
                    gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel),
                    gettemplate("ggH", self.analysis, self.production, self.category, "g11gi1", self.channel),
                   ]
        for i, t in enumerate(T["ggH"], start=1):
            t.SetName("ggH_T_ZZ_{}_{}_3D_{}".format(self.production.year,self.channel, i))

        T_ScaleResUp["ggH"] = [
                               gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, "ScaleResUp"),
                               gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, "ScaleResUp"),
                               gettemplate("ggH", self.analysis, self.production, self.category, "g11gi1", self.channel, "ScaleResUp"),
                              ]
        for i, t in enumerate(T_ScaleResUp["ggH"], start=1):
            t.SetName("ggH_T_ZZ_{}_{}_3D_{}_ScaleResUp".format(self.production.year,self.channel, i))

        T_ScaleResDown["ggH"] = [
                                 gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, "ScaleResDown"),
                                 gettemplate("ggH", self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, "ScaleResDown"),
                                 gettemplate("ggH", self.analysis, self.production, self.category, "g11gi1", self.channel, "ScaleResDown"),
                                ]
        for i, t in enumerate(T_ScaleResDown["ggH"], start=1):
            t.SetName("ggH_T_ZZ_{}_{}_3D_{}_ScaleResDown".format(self.production.year,self.channel, i))

        T_integralName = ["ggH_normt{}_{}_{}_{}".format(i, self.channel, self.category, self.production.year) for i in range(1, 4)]
        T_integral["ggH"] = [ROOT.RooConstVar(integralName, integralName, t.Integral()) for t, integralName in zip(T["ggH"], T_integralName)]
        for i, integral in enumerate(T_integral["ggH"], start=1):
            print "ggH T{}".format(i), integral.getVal()

        r_fai_pures_norm_Name = "ggH_PuresNorm_{}_{}_{}".format(self.channel, self.category, self.production.year)
        r_fai_realints_norm_Name = "ggH_RealIntsNorm_{}_{}_{}".format(self.channel, self.category, self.production.year)
        r_fai_pures_norm = ROOT.RooFormulaVar(r_fai_pures_norm_Name,r_fai_pures_norm_Name,"( (1-abs(@0))*@1+abs(@0)*@2 )/@1",ROOT.RooArgList(x,T_integral["ggH"][0],T_integral["ggH"][1]))
        r_fai_realints_norm = ROOT.RooFormulaVar(r_fai_realints_norm_Name,r_fai_realints_norm_Name,"( sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1 )/@2",ROOT.RooArgList(x,T_integral["ggH"][2],T_integral["ggH"][0]))
        norm["ggH"] = ROOT.RooFormulaVar("ggH_norm","ggH_norm","(abs(@2))>1 ? 0. : TMath::Max((@0+@1),0)",ROOT.RooArgList(r_fai_pures_norm,r_fai_realints_norm,x))

        T_hist["ggH"] = [ROOT.RooDataHist("ggH_T_{}_hist".format(i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T["ggH"], start=1)]
        T_ScaleResUp_hist["ggH"] = [ROOT.RooDataHist("ggH_T_{}_ScaleResUp_hist".format(i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T_ScaleResUp["ggH"], start=1)]
        T_ScaleResDown_hist["ggH"] = [ROOT.RooDataHist("ggH_T_{}_ScaleResDown_hist".format(i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T_ScaleResDown["ggH"], start=1)]

        T_histfunc["ggH"] = [ROOT.RooHistFunc("ggH_T_{}_histfunc".format(i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_hist["ggH"], start=1)]
        T_ScaleResUp_histfunc["ggH"] = [ROOT.RooHistFunc("ggH_T_{}_histfunc_ScaleResUp".format(i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_ScaleResUp_hist["ggH"], start=1)]
        T_ScaleResDown_histfunc["ggH"] = [ROOT.RooHistFunc("ggH_T_{}_histfunc_ScaleResDown".format(i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_ScaleResDown_hist["ggH"], start=1)]

        pdf["ggH"] = ROOT.HZZ4L_RooSpinZeroPdf("ggH", "ggH", D1, D2, D3, x, ROOT.RooArgList(*T_histfunc["ggH"]))

        pdfName_syst1Up = "ggH_ScaleRes{}Up".format(self.channel)
        pdfName_syst1Down = "ggH_ScaleRes{}Down".format(self.channel)
        pdf_syst1Up["ggH"] = ROOT.HZZ4L_RooSpinZeroPdf(pdfName_syst1Up, pdfName_syst1Up, D1, D2, D3, x, ROOT.RooArgList(*T_ScaleResUp_histfunc["ggH"]))
        pdf_syst1Down["ggH"] = ROOT.HZZ4L_RooSpinZeroPdf(pdfName_syst1Down, pdfName_syst1Down, D1, D2, D3, x, ROOT.RooArgList(*T_ScaleResDown_histfunc["ggH"]))

        for p in "qqH", "ZH", "WH":
            #for prod+dec the order is a1^4 a1^3ai a1^2ai^2 a1ai^3 ai^4
            T[p] = [
                    gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel),
                    gettemplate(p, self.analysis, self.production, self.category, "g13gi1", self.channel),
                    gettemplate(p, self.analysis, self.production, self.category, "g12gi2", self.channel),
                    gettemplate(p, self.analysis, self.production, self.category, "g11gi3", self.channel),
                    gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel),
                   ]
            for i, t in enumerate(T[p], start=1):
                t.SetName("{}_T_ZZ_{}_{}_{}_3D_{}".format(p, self.production.year, self.channel, self.category, i))

            T_ScaleResUp[p] = [
                               gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, "ScaleResUp"),
                               gettemplate(p, self.analysis, self.production, self.category, "g13gi1", self.channel, "ScaleResUp"),
                               gettemplate(p, self.analysis, self.production, self.category, "g12gi2", self.channel, "ScaleResUp"),
                               gettemplate(p, self.analysis, self.production, self.category, "g11gi3", self.channel, "ScaleResUp"),
                               gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, "ScaleResUp"),
                              ]
            for i, t in enumerate(T_ScaleResUp[p], start=1):
                t.SetName("{}_T_ZZ_{}_{}_{}_3D_{}_ScaleResUp".format(p, self.production.year, self.channel, self.category, i))

            T_ScaleResDown[p] = [
                                 gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[0], self.channel, "ScaleResDown"),
                                 gettemplate(p, self.analysis, self.production, self.category, "g13gi1", self.channel, "ScaleResDown"),
                                 gettemplate(p, self.analysis, self.production, self.category, "g12gi2", self.channel, "ScaleResDown"),
                                 gettemplate(p, self.analysis, self.production, self.category, "g11gi3", self.channel, "ScaleResDown"),
                                 gettemplate(p, self.analysis, self.production, self.category, self.analysis.purehypotheses[1], self.channel, "ScaleResDown"),
                                ]
            for i, t in enumerate(T_ScaleResDown[p], start=1):
                t.SetName("{}_T_ZZ_{}_{}_{}_3D_{}_ScaleResDown".format(p, self.production.year, self.channel, self.category, i))

            T_integralName = ["{}_normt{}_{}_{}_{}".format(p, i, self.channel, self.category, self.production.year) for i in range(1, 6)]
            T_integral[p] = [ROOT.RooConstVar(T_integralName[i], T_integralName[i], t.Integral()) for i, t in enumerate(T[p])]
            for i, integral in enumerate(T_integral[p], start=1):
                print "{} T{}".format(p, i), integral.getVal()

            normname = "{}_norm".format(p)
            formula = " + ".join("@0**{}*@1**{}*@{}".format(4-i, i, i+2) for i in range(5))
            formula = "("+formula+") / @2"
            norm[p] = ROOT.RooFormulaVar(normname, formula, ROOT.RooArgList(a1, ai, *T_integral[p]))

            T_hist[p] = [ROOT.RooDataHist("{}_T_{}_hist".format(p, i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T[p], start=1)]
            T_ScaleResUp_hist[p] = [ROOT.RooDataHist("{}_T_{}_ScaleResUp_hist".format(p, i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T_ScaleResUp[p], start=1)]
            T_ScaleResDown_hist[p] = [ROOT.RooDataHist("{}_T_{}_ScaleResDown_hist".format(p, i), "", ROOT.RooArgList(D1,D2,D3), t) for i, t in enumerate(T_ScaleResDown[p], start=1)]

            T_histfunc[p] = [ROOT.RooHistFunc("{}_T_{}_histfunc".format(p, i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_hist[p], start=1)]
            T_ScaleResUp_histfunc[p] = [ROOT.RooHistFunc("{}_T_{}_histfunc_ScaleResUp".format(p, i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_ScaleResUp_hist[p], start=1)]
            T_ScaleResDown_histfunc[p] = [ROOT.RooHistFunc("{}_T_{}_histfunc_ScaleResDown".format(p, i), "", ROOT.RooArgSet(D1,D2,D3), datahist) for i, datahist in enumerate(T_ScaleResDown_hist[p], start=1)]

            pdf[p] = ROOT.VBFHZZ4L_RooSpinZeroPdf(p, p, D1, D2, D3, a1, ai, ROOT.RooArgList(*T_histfunc[p]))

            pdfName_syst1Up = "{}_ScaleRes{}{}Up".format(p, self.category, self.channel)
            pdfName_syst1Down = "{}_ScaleRes{}{}Down".format(p, self.category, self.channel)
            pdf_syst1Up[p] = ROOT.VBFHZZ4L_RooSpinZeroPdf(pdfName_syst1Up, pdfName_syst1Up, D1, D2, D3, a1, ai, ROOT.RooArgList(*T_ScaleResUp_histfunc[p]))
            pdf_syst1Down[p] = ROOT.VBFHZZ4L_RooSpinZeroPdf(pdfName_syst1Down, pdfName_syst1Down, D1, D2, D3, a1, ai, ROOT.RooArgList(*T_ScaleResDown_histfunc[p]))

        ## ------------------ END 2D SIGNAL SHAPES FOR PROPERTIES ------------------------ ##


        ## ------------------ 2D BACKGROUND SHAPES FOR PROPERTIES ------------------- ##

        qqZZTemplate = gettemplate(self.analysis, self.production, self.category, "qqZZ", self.channel)

        TemplateName = "qqZZTempDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        qqZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),qqZZTemplate)
        PdfName = "bkg_qqzz"
        qqZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),qqZZTempDataHist)

        ggZZTemplate = gettemplate(self.analysis, self.production, self.category, "ggZZ", self.channel)

        TemplateName = "ggZZTempDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ggZZTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ggZZTemplate)
        PdfName = "bkg_ggzz"
        ggZZTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ggZZTempDataHist)

        VBFbkgTemplate = gettemplate(self.analysis, self.production, self.category, "VBFbkg", self.channel)

        TemplateName = "VBFbkgTempDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        VBFbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),VBFbkgTemplate)
        PdfName = "bkg_vbf"
        VBFbkgTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),VBFbkgTempDataHist)

        ZjetsTemplate = gettemplate(self.analysis, self.production, self.category, "ZX", self.channel)
        TemplateName = "ZjetsTempDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplate)
        PdfName = "Zjets_TemplatePdf_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTemplatePdf = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHist)

        ZjetsTemplateDown = gettemplate(self.analysis, self.production, self.category, "ZX", self.channel, "ZXDown")
        TemplateName = "ZjetsTempDownDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTempDataHistDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplateDown)
        PdfName = "Zjets_TemplateDownPdf_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTemplatePdfDown = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHistDown)

        ZjetsTemplateUp = gettemplate(self.analysis, self.production, self.category, "ZX", self.channel, "ZXUp")
        TemplateName = "ZjetsTempUpDataHist_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTempDataHistUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(D1,D2,D3),ZjetsTemplateUp)
        PdfName = "Zjets_TemplateUpPdf_{}_{}_{}".format(self.channel,self.category,self.production.year)
        ZjetsTemplatePdfUp = ROOT.RooHistPdf(PdfName,PdfName,ROOT.RooArgSet(D1,D2,D3),ZjetsTempDataHistUp)

        funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{}_{}".format(self.channel, self.category)
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()

        funcList_zjets.add(ZjetsTemplatePdf)
        funcList_zjets.add(ZjetsTemplatePdfUp)
        funcList_zjets.add(ZjetsTemplatePdfDown)
        alphaMorphBkg.setConstant(False)
        morphVarListBkg.add(alphaMorphBkg)

        MorphName = "bkg_zjets"
        ZjetsTemplateMorphPdf = ROOT.FastVerticalInterpHistPdf3D(MorphName,MorphName,D1,D2,D3,False,funcList_zjets,morphVarListBkg,1.0,1)


        ## ---------------- END 2D BACKGROUND SHAPES FOR PROPERTIES ----------------- ##

        ## --------------------------- DATASET --------------------------- ##

        data_obs_tree = getdatatree(self.channel, self.production, self.category, self.analysis)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{}".format(self.channel)


        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(D1,D2,D3))


        ## --------------------------- WORKSPACE -------------------------- ##

        w = ROOT.RooWorkspace("w","w")

        w.importClassCode(ROOT.HZZ4L_RooSpinZeroPdf.Class(),True)
        w.importClassCode(ROOT.VBFHZZ4L_RooSpinZeroPdf.Class(),True)

        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
        #getattr(w,'import')(r_fai_norm) ### Should this be renamed?


        for p in ("ggH", "qqH", "ZH", "WH"):
            getattr(w,'import')(pdf[p], ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(norm[p], ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(pdf_syst1Up[p], ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(pdf_syst1Down[p], ROOT.RooFit.RecycleConflictNodes())

        getattr(w,'import')(qqZZTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(ggZZTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(VBFbkgTemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(ZjetsTemplateMorphPdf, ROOT.RooFit.RecycleConflictNodes())


        w.writeToFile(self.rootfile_base)

    def makeworkspace(self):
        self.writeworkspace()
        if not os.path.exists(self.rootfile):
            os.symlink(self.rootfile_base, self.rootfile)

    def makeCardsWorkspaces(self, outdir="."):
        mkdir_p(outdir)
        with cd(outdir):
            self.makeworkspace()
            self.writedatacard()

def makeDCandWS(*args, **kwargs):
    Datacard(*args).makeCardsWorkspaces(**kwargs)
