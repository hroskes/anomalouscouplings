from abc import ABCMeta, abstractproperty
import config
import constants
from enums import BlindStatus, Flavor, decayonlyhypotheses, prodonlyhypotheses, proddechypotheses, hffhypotheses, Hypothesis, MultiEnum, MultiEnumABCMeta, ProductionMode, Production
from math import sqrt
import os
import ROOT


class SampleBase(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def g1(self):
        pass
    @abstractproperty
    def g2(self):
        pass
    @abstractproperty
    def g4(self):
        pass
    @abstractproperty
    def g1prime2(self):
        pass

    @property
    def JHUxsec(self):
        if self.productionmode == "ggH":
            return (
                      constants.JHUXSggH2L2la1*self.g1**2
                    + constants.JHUXSggH2L2la2*self.g2**2
                    + constants.JHUXSggH2L2la3*self.g4**2
                    + constants.JHUXSggH2L2lL1*self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )

        if self.productionmode == "VBF":
            return (
                      constants.JHUXSVBFa1 * self.g1**2
                    + constants.JHUXSVBFa2 * self.g2**2
                    + constants.JHUXSVBFa3 * self.g4**2
                    + constants.JHUXSVBFL1 * self.g1prime2**2
                    + constants.JHUXSVBFa1a2 * self.g1*self.g2 / constants.g2VBF
                    + constants.JHUXSVBFa1a3 * self.g1*self.g4 / constants.g4VBF
                    + constants.JHUXSVBFa1L1 * self.g1*self.g1prime2 / constants.g1prime2VBF_gen
                   ) * (
                      constants.JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )

        if self.productionmode == "ZH":
            return (
                      constants.JHUXSZHa1 * self.g1**2
                    + constants.JHUXSZHa2 * self.g2**2
                    + constants.JHUXSZHa3 * self.g4**2
                    + constants.JHUXSZHL1 * self.g1prime2**2
                    + constants.JHUXSZHa1a2 * self.g1*self.g2 / constants.g2ZH
                    + constants.JHUXSZHa1a3 * self.g1*self.g4 / constants.g4ZH
                    + constants.JHUXSZHa1L1 * self.g1*self.g1prime2 / constants.g1prime2ZH_gen
                   ) * (
                      constants.JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )

        if self.productionmode == "WH":
            return (
                      constants.JHUXSWHa1 * self.g1**2
                    + constants.JHUXSWHa2 * self.g2**2
                    + constants.JHUXSWHa3 * self.g4**2
                    + constants.JHUXSWHL1 * self.g1prime2**2
                    + constants.JHUXSWHa1a2 * self.g1*self.g2 / constants.g2WH
                    + constants.JHUXSWHa1a3 * self.g1*self.g4 / constants.g4WH
                    + constants.JHUXSWHa1L1 * self.g1*self.g1prime2 / constants.g1prime2WH_gen
                   ) * (
                      constants.JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )
        if self.productionmode == "HJJ":
            return (
                      constants.JHUXSHJJa2 * self.ghg2**2
                    + constants.JHUXSHJJa3 * self.ghg4**2
                    + constants.JHUXSHJJa2a3 * self.ghg2*self.ghg4 / constants.ghg4HJJ
                   ) * (
                      constants.JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )
        if self.productionmode == "ttH":
            return (
                      constants.JHUXSttHkappa * self.kappa**2
                    + constants.JHUXSttHkappatilde * self.kappa_tilde**2
                    + constants.JHUXSttHkappakappatilde * self.kappa*self.kappa_tilde / constants.kappa_tilde_ttH
                   ) * (
                      constants.JHUXSggH2L2la1 * self.g1**2
                    + constants.JHUXSggH2L2la2 * self.g2**2
                    + constants.JHUXSggH2L2la3 * self.g4**2
                    + constants.JHUXSggH2L2lL1 * self.g1prime2**2
                    + constants.JHUXSggH2L2la1a2 * self.g1*self.g2 / constants.g2decay
                    + constants.JHUXSggH2L2la1a3 * self.g1*self.g4 / constants.g4decay
                    + constants.JHUXSggH2L2la1L1 * self.g1*self.g1prime2 / constants.g1prime2decay_gen
                   )
        assert False

    @property
    def xsec(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "SM": return constants.SMXSggH2L2l
            return constants.SMXSggH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "VBF":
            if self.hypothesis == "SM": return constants.SMXSVBF2L2l
            return constants.SMXSVBF2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "ZH":
            if self.hypothesis == "SM": return constants.SMXSZH2L2l
            return constants.SMXSZH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "WH":
            if self.hypothesis == "SM": return constants.SMXSWH2L2l
            return constants.SMXSWH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "HJJ":
            if self.hypothesis == "SM": return constants.SMXSHJJ2L2l
            return constants.SMXSHJJ2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "ttH":
            if self.hypothesis == "SM": return constants.SMXSttH2L2l
            return constants.SMXSttH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        assert False


class ArbitraryCouplingsSample(SampleBase):
    def __init__(self, productionmode, g1, g2, g4, g1prime2, ghg2=None, ghg4=None, kappa=None, kappa_tilde=None):
        self.productionmode = ProductionMode(productionmode)
        self.__g1, self.__g2, self.__g4, self.__g1prime2 = g1, g2, g4, g1prime2
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH"):
            raise ValueError("Bad productionmode {}".format(self.productionmode))
        if sum(bool(g) for g in (g2, g4, g1prime2)) > 1:
            raise ValueError("Can only set at most one of g2, g4, g1prime2")

        if self.productionmode == "HJJ":
            self.__ghg2, self.__ghg4 = ghg2, ghg4
            if ghg2 is None or ghg4 is None:
                raise ValueError("Have to set ghg2 and ghg4 for HJJ")
        else:
            if ghg2 is not None or ghg4 is not None:
                raise ValueError("Can't set ghg2 or ghg4 for {}".format(self.productionmode))

        if self.productionmode == "ttH":
            self.__kappa, self.__kappa_tilde = kappa, kappa_tilde
            if kappa is None or kappa_tilde is None:
                raise ValueError("Have to set kappa and kappa_tilde for ttH")
        else:
            if kappa is not None or kappa_tilde is not None:
                raise ValueError("Can't set kappa or kappa_tilde for {}".format(self.productionmode))

    @property
    def g1(self):
        return self.__g1
    @property
    def g2(self):
        return self.__g2
    @property
    def g4(self):
        return self.__g4
    @property
    def g1prime2(self):
        return self.__g1prime2
    @property
    def ghg2(self):
        return self.__ghg2
    @property
    def ghg4(self):
        return self.__ghg4
    @property
    def kappa(self):
        return self.__kappa
    @property
    def kappa_tilde(self):
        return self.__kappa_tilde

class ReweightingSample(MultiEnum, SampleBase):
    __metaclass__ = MultiEnumABCMeta
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, Flavor]

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if (
                   self.productionmode == "ggH" and self.hypothesis not in decayonlyhypotheses
                or self.productionmode in ("HJJ", "ttH") and self.hypothesis not in hffhypotheses
               ):
                raise ValueError("{} hypothesis can't be {}\n{}".format(self.productionmode, self.hypothesis, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ggZZ", "VBF bkg"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is None:
                raise ValueError("No flavor provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.productionmode == "VBF bkg" and self.flavor.hastaus:
                raise ValueError("No {} samples with taus\n{}".format(self.productionmode, args))
        elif self.productionmode in ("qqZZ", "ZX", "data"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    def reweightingsamples(self):
        if self.productionmode == "ggH":
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in decayonlyhypotheses]
        elif self.productionmode in ("VBF", "ZH", "WH"):
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in proddechypotheses]
        elif self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX", "ttH", "HJJ"):
            return [self]
        elif self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    def isbkg(self):
        if self.productionmode in ("ggH", "data", "VBF", "ZH", "WH", "ttH", "HJJ"):
            return False
        elif self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX"):
            return True
        raise self.ValueError("isbkg")

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ"):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "VBF", "ZH", "WH", "ttH", "HJJ"):
            return False
        elif self.productionmode == "data":
            return True
        raise self.ValueError("isdata")
        
    def weightname(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return "MC_weight_ggH_g1"
            elif self.hypothesis == "a2":
                return "MC_weight_ggH_g2"
            elif self.hypothesis == "0-":
                return "MC_weight_ggH_g4"
            elif self.hypothesis == "L1":
                return "MC_weight_ggH_g1prime2"
            elif self.hypothesis == "fa20.5":
                return "MC_weight_ggH_g1g2"
            elif self.hypothesis == "fa30.5":
                return "MC_weight_ggH_g1g4"
            elif self.hypothesis == "fL10.5":
                return "MC_weight_ggH_g1g1prime2"
        elif self.productionmode in ("VBF", "ZH", "WH"):
            if self.hypothesis == "0+":
                return "MC_weight_{}_g1".format(self.productionmode)
            elif self.hypothesis == "a2":
                return "MC_weight_{}_g2".format(self.productionmode)
            elif self.hypothesis == "0-":
                return "MC_weight_{}_g4".format(self.productionmode)
            elif self.hypothesis == "L1":
                return "MC_weight_{}_g1prime2".format(self.productionmode)
            elif self.hypothesis == "fa2dec0.5":
                return "MC_weight_{}_g1g2_dec".format(self.productionmode)
            elif self.hypothesis == "fa3dec0.5":
                return "MC_weight_{}_g1g4_dec".format(self.productionmode)
            elif self.hypothesis == "fL1dec0.5":
                return "MC_weight_{}_g1g1prime2_dec".format(self.productionmode)
            elif self.hypothesis == "fa2prod0.5":
                return "MC_weight_{}_g1g2_prod".format(self.productionmode)
            elif self.hypothesis == "fa3prod0.5":
                return "MC_weight_{}_g1g4_prod".format(self.productionmode)
            elif self.hypothesis == "fL1prod0.5":
                return "MC_weight_{}_g1g1prime2_prod".format(self.productionmode)
            elif self.hypothesis == "fa2proddec-0.5":
                return "MC_weight_{}_g1g2_proddec_pi".format(self.productionmode)
            elif self.hypothesis == "fa3proddec-0.5":
                return "MC_weight_{}_g1g4_proddec_pi".format(self.productionmode)
            elif self.hypothesis == "fL1proddec-0.5":
                return "MC_weight_{}_g1g1prime2_proddec_pi".format(self.productionmode)
        elif self.productionmode == "HJJ":
            if self.hypothesis == "0+":
                return "MC_weight_HJJ_g2"
            elif self.hypothesis == "0-":
                return "MC_weight_HJJ_g4"
            elif self.hypothesis == "fCP0.5":
                return "MC_weight_HJJ_g2g4"
        elif self.productionmode == "ttH":
            if self.hypothesis == "0+":
                return "MC_weight_ttH_kappa"
            elif self.hypothesis == "0-":
                return "MC_weight_ttH_kappatilde"
            elif self.hypothesis == "fCP0.5":
                return "MC_weight_ttH_kappakappatilde"
        elif self.productionmode == "ggZZ":
            return "MC_weight_ggZZ"
        elif self.productionmode == "qqZZ":
            return "MC_weight_qqZZ"
        elif self.productionmode == "VBF bkg":
            return "MC_weight_VBFbkg"
        elif self.productionmode == "ZX":
            return "MC_weight_ZX"
        raise self.ValueError("weightname")

    def TDirectoryname(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBFbkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ") or self.productionmode == "ZX" and not config.usedata:
            return "ZZTree"
        if self.productionmode == "ZX":
            return "CRZLLTree"
        raise self.ValueError("TDirectoryname")

    def color(self):
        if self.productionmode == "ggH":
            if self.hypothesis == "0+":
                return 1
            elif self.hypothesis == "0-":
                return 2
            elif self.hypothesis == "fa30.5":
                return ROOT.kRed+3
            elif self.hypothesis == "a2":
                return ROOT.kCyan
            elif self.hypothesis == "fa20.5":
                return 4
            elif self.hypothesis == "L1":
                return 3
            elif self.hypothesis == "fL10.5":
                return ROOT.kGreen+3
        raise self.ValueError("color")

    def onlyweights(self):
        """True if this sample is not needed for making templates,
           and only the weight and ZZMass should be recorded in the tree"""
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "data", "HJJ", "ttH"):
            return False
        if self.productionmode in ():
            return True
        raise self.ValueError("onlyweights")

    def weightingredients(self):
        """only needs to be defined if self.onlyweights() is True.
           gives a list of variable names that are used to calculate
           the weight"""
        raise self.ValueError("weightingredients")

    @property
    def g1(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 1
        if self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            if self.hypothesis in ["0+"] + ["{}{}0.5".format(a, b) for a in ("fa2", "fa3", "fL1") for b in ("dec", "prod", "proddec-")]:
                return 1
            if self.hypothesis in ("a2", "0-", "L1"):
                return 0
        raise self.ValueError("g1")

    @property
    def g2(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 0
        if self.hypothesis in ["0+", "0-", "L1"] + ["{}{}0.5".format(a, b) for a in ("fa3", "fL1") for b in ("dec", "prod", "proddec-")]:
            return 0
        if self.hypothesis == "a2":
            if self.productionmode == "ggH":
                return constants.g2decay
            else:
                return 1
        if self.hypothesis == "fa2dec0.5":
            return constants.g2decay

        if self.productionmode == "VBF":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2VBF
            if self.hypothesis == "fa2proddec-0.5":
                return -sqrt(constants.g2VBF*constants.g2decay)

        if self.productionmode == "ZH":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2ZH
            if self.hypothesis == "fa2proddec-0.5":
                return -sqrt(constants.g2ZH*constants.g2decay)

        if self.productionmode == "WH":
            if self.hypothesis == "fa2prod0.5":
                return constants.g2WH
            if self.hypothesis == "fa2proddec-0.5":
                return -sqrt(constants.g2WH*constants.g2decay)

        raise self.ValueError("g2")

    @property
    def g4(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 0
        if self.hypothesis in ["0+", "a2", "L1"] + ["{}{}0.5".format(a, b) for a in ("fa2", "fL1") for b in ("dec", "prod", "proddec-")]:
            return 0
        if self.hypothesis == "0-":
            if self.productionmode == "ggH":
                return constants.g4decay
            else:
                return 1
        if self.hypothesis == "fa3dec0.5":
            return constants.g4decay

        if self.productionmode == "VBF":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4VBF
            if self.hypothesis == "fa3proddec-0.5":
                return -sqrt(constants.g4VBF*constants.g4decay)

        if self.productionmode == "ZH":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4ZH
            if self.hypothesis == "fa3proddec-0.5":
                return -sqrt(constants.g4ZH*constants.g4decay)

        if self.productionmode == "WH":
            if self.hypothesis == "fa3prod0.5":
                return constants.g4WH
            if self.hypothesis == "fa3proddec-0.5":
                return -sqrt(constants.g4WH*constants.g4decay)

        raise self.ValueError("g4")

    @property
    def g1prime2(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 0
        if self.hypothesis in ["0+", "a2", "0-"] + ["{}{}0.5".format(a, b) for a in ("fa2", "fa3") for b in ("dec", "prod", "proddec-")]:
            return 0
        if self.hypothesis == "L1":
            if self.productionmode == "ggH":
                return constants.g1prime2decay_gen
            else:
                return 1
        if self.hypothesis == "fL1dec0.5":
            return constants.g1prime2decay_gen

        if self.productionmode == "VBF":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2VBF_gen
            if self.hypothesis == "fL1proddec-0.5":
                return -sqrt(constants.g1prime2VBF_gen*constants.g1prime2decay_gen)

        if self.productionmode == "ZH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2ZH_gen
            if self.hypothesis == "fL1proddec-0.5":
                return -sqrt(constants.g1prime2ZH_gen*constants.g1prime2decay_gen)

        if self.productionmode == "WH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2WH_gen
            if self.hypothesis == "fL1proddec-0.5":
                return -sqrt(constants.g1prime2WH_gen*constants.g1prime2decay_gen)

        raise self.ValueError("g1prime2")

    @property
    def ghg2(self):
        if self.productionmode == "HJJ":
            if self.hypothesis in ("0+", "fCP0.5"):
                return 1
            if self.hypothesis == "0-":
                return 0
        raise self.ValueError("ghg2")

    @property
    def ghg4(self):
        if self.productionmode == "HJJ":
            if self.hypothesis == "0+":
                return 0
            if self.hypothesis == "0-":
                return 1
            if self.hypothesis == "fCP0.5":
                return constants.ghg4HJJ

    @property
    def kappa(self):
        if self.productionmode == "ttH":
            if self.hypothesis in ("0+", "fCP0.5"):
                return 1
            if self.hypothesis == "0-":
                return 0
        raise self.ValueError("kappa")

    @property
    def kappa_tilde(self):
        if self.productionmode == "ttH":
            if self.hypothesis == "0+":
                return 0
            if self.hypothesis == "0-":
                return 1
            if self.hypothesis == "fCP0.5":
                return constants.kappa_tilde_ttH

class Sample(ReweightingSample):
    enums = [ReweightingSample, Production, BlindStatus]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))
        if self.blindstatus is None and self.productionmode == "data":
            raise ValueError("No blindstatus provided for data sample!\n{}".format(args))
        if self.blindstatus is not None and self.productionmode != "data":
            raise ValueError("blindstatus provided for MC sample!\n{}".format(args))
        if self.productionmode in ("VBF", "ZH", "WH") and self.hypothesis not in prodonlyhypotheses:
            raise ValueError("No {} sample produced with hypothesis {}!\n{}".format(self.productionmode, self.hypothesis, args))
        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
        if self.productionmode == "ggH":
            return self.production.CJLSTdir_anomalous()
        if self.productionmode == "VBF":
            return self.production.CJLSTdir_anomalous_VBF()
        if self.productionmode in ("ZH", "WH"):
            return self.production.CJLSTdir_anomalous_VH()
        if self.productionmode in ("HJJ", "ttH"):
            return self.production.CJLSTdir_HJJttH()
        if self.productionmode in ("data", "ZX"):
            return self.production.CJLSTdir_data()
        return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.productionmode == "ggH" and self.production <= "160624":
            if self.hypothesis == "0+": return "0PM"
            if self.hypothesis == "a2": return "0PH"
            if self.hypothesis == "0-": return "0M"
            if self.hypothesis == "L1": return "0L1"
            if self.hypothesis == "fa20.5": return "0PHf05ph0"
            if self.hypothesis == "fa30.5": return "0Mf05ph0"
            if self.hypothesis == "fL10.5": return "0L1f05ph0"
        if (self.productionmode in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH") and self.production >= "160714"):
            s = str(self.productionmode)
            if self.productionmode == "VBF": s = "VBFH"
            if self.hypothesis == "0+": return "{}0PM_M125".format(s)
            if self.hypothesis == "a2": return "{}0PH_M125".format(s)
            if self.hypothesis == "0-": return "{}0M_M125".format(s)
            if self.hypothesis == "L1": return "{}0L1_M125".format(s)
            if self.hypothesis in ("fa20.5", "fa2prod0.5"): return "{}0PHf05ph0_M125".format(s)
            if self.hypothesis in ("fa30.5", "fa3prod0.5"): return "{}0Mf05ph0_M125".format(s)
            if self.hypothesis in ("fL10.5", "fL1prod0.5"): return "{}0L1f05ph0_M125".format(s)
        if self.productionmode == "ggZZ":
            if self.production == "160714" and self.flavor in ("4e", "4mu", "4tau", "2mu2tau") or self.production >= "160720":
                return "ggTo{}_Contin_MCFM701".format(self.flavor)
            elif self.production <= "160714":
                return "ggZZ{}".format(self.flavor)
        if self.productionmode == "VBF bkg":
            return "VBFTo{}JJ_Contin_phantom128".format(self.flavor)
        if self.productionmode == "qqZZ":
            return "ZZTo4l"
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")

    def withdiscriminantsfile(self):
        result = os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        return result
        raise self.ValueError("withdiscriminantsfile")

    @property
    def useMELAv2(self):
        return self.production.useMELAv2

    @property
    def unblind(self):
        if self.productionmode == "data":
            return self.blindstatus == "unblind"
        raise self.ValueError("unblind")
