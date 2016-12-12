#!/usr/bin/env python
from abc import ABCMeta, abstractproperty
import config
import constants
from enums import AlternateGenerator, Analysis, BlindStatus, Flavor, decayonlyhypotheses, prodonlyhypotheses, proddechypotheses, purehypotheses, hffhypotheses, Hypothesis, MultiEnum, MultiEnumABCMeta, ProductionMode, Production
from math import copysign, sqrt
import numpy
import os
import ROOT
from utilities import cache


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
    @abstractproperty
    def ghg2(self):
        pass
    @abstractproperty
    def ghg4(self):
        pass
    @abstractproperty
    def kappa(self):
        pass
    @abstractproperty
    def kappa_tilde(self):
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
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSggH2L2l
            return constants.SMXSggH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "VBF":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSVBF2L2l
            return constants.SMXSVBF2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "ZH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSZH2L2l
            return constants.SMXSZH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "WH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWH2L2l
            return constants.SMXSWH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "WplusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWpH2L2l
        if self.productionmode == "WminusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWmH2L2l
        if self.productionmode == "HJJ":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSHJJ2L2l
            return constants.SMXSHJJ2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        if self.productionmode == "ttH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSttH2L2l
            return constants.SMXSttH2L2l * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec
        assert False

    def __eq__(self, other):
        return all((
                   self.productionmode == other.productionmode,
                   self.g1 == other.g1,
                   self.g2 == other.g2,
                   self.g4 == other.g4,
                   self.g1prime2 == other.g1prime2,
               )) and (
                   self.productionmode != "HJJ" or all((
                       self.ghg2 == other.ghg2,
                       self.ghg4 == other.ghg4,
                   ))
               ) and (
                   self.productionmode != "ttH" or all((
                       self.kappa == other.kappa,
                       self.kappa_tilde == other.kappa_tilde,
                   ))
               )
    def __ne__(self, other):
        return not self == other

    def linearcombinationofreweightingsamples(self, hypotheses):
        g1 = self.g1
        if   self.g2 == self.g1prime2 == 0: analysis = Analysis("fa3"); gi = self.g4
        elif self.g4 == self.g1prime2 == 0: analysis = Analysis("fa2"); gi = self.g2
        elif self.g2 == self.g4 == 0:       analysis = Analysis("fL1"); gi = self.g1prime2
        else: assert False

        if hypotheses == "templates":
            from templates import TemplatesFile
            templatesfile = TemplatesFile(analysis, str(self.productionmode).lower(), "2e2mu", "Untagged", config.productionsforcombine[0])
            hypotheses = [_.hypothesis for _ in templatesfile.templates()]
        if hypotheses == "directlyreweighted":
            if self.productionmode == "ggH":
                if analysis == "fa3": hypotheses = ["0+", "0-", "fa30.5"]
                if analysis == "fa2": hypotheses = ["0+", "a2", "fa20.5"]
                if analysis == "fL1": hypotheses = ["0+", "L1", "fL10.5"]
            elif self.productionmode in ("VBF", "WH", "ZH"):
                if analysis == "fa3": hypotheses = ["0+", "0-", "fa3dec0.5", "fa3prod0.5", "fa3proddec-0.5"]
                if analysis == "fa2": hypotheses = ["0+", "a2", "fa2dec0.5", "fa2prod0.5", "fa2proddec-0.5"]
                if analysis == "fL1": hypotheses = ["0+", "L1", "fL1dec0.5", "fL1prod0.5", "fL1proddec0.5"]

        basis = SampleBasis(hypotheses, self.productionmode, analysis)
        vectorofreweightingsamples = [ReweightingSample(self.productionmode, hypothesis) for hypothesis in hypotheses]
        invertedmatrix = basis.invertedmatrix
        maxpower = invertedmatrix.shape[0]-1
        vectorTofg1xgiy = numpy.matrix([[g1**(maxpower-i)*gi**i for i in range(maxpower+1)]])
        vectorTmatrix = vectorTofg1xgiy*invertedmatrix
        vectorTmatrixaslist = numpy.array(vectorTmatrix).flatten().tolist()
        return {
                reweightingsample: factor
                 for reweightingsample, factor in zip(vectorofreweightingsamples, vectorTmatrixaslist)
                 if factor
               }

    @property
    def MC_weight(self):
        return " + ".join(
                          "{}*{}".format(reweightingsample.weightname(), factor*reweightingsample.xsec/self.xsec)
                            for reweightingsample, factor in self.linearcombinationofreweightingsamples("directlyreweighted").iteritems()
                         )

    def get_MC_weight_function(self_sample):
        terms = tuple(
                      (reweightingsample.weightname(), factor*reweightingsample.xsec/self_sample.xsec)
                               for reweightingsample, factor in self_sample.linearcombinationofreweightingsamples("directlyreweighted").iteritems()
                     )
        def MC_weight_function(self_tree):
            return sum(factor*getattr(self_tree, componentweightname) for componentweightname, factor in terms)
        MC_weight_function.__name__ = self_sample.weightname()
        return MC_weight_function

    def fai(self, productionmode, analysis, withdecay=False):
        from combinehelpers import mixturesign
        analysis = Analysis(analysis)
        hypothesis = analysis.purehypotheses[1]
        mixhypothesis = analysis.mixdecayhypothesis
        if productionmode == "VH":
            productionmodes = [ProductionMode("ZH"), ProductionMode("WH")]
        else:
            productionmodes = [ProductionMode(productionmode)]

        if productionmode == "ggH":
            withdecay = True

        power = 4 if productionmode != "ggH" and withdecay else 2
        numerator = None
        denominator = 0
        for h in purehypotheses:
            kwargs = {_.couplingname: 1 if _ == h else 0 for _ in purehypotheses}
            term = 0
            for productionmode in productionmodes:
                term += getattr(self, h.couplingname)**power * ArbitraryCouplingsSample(productionmode, **kwargs).xsec
            if not withdecay:
                term /= ArbitraryCouplingsSample("ggH", **kwargs).xsec
            denominator += term
            if h == hypothesis:
                numerator = term
        return copysign(
                        numerator / denominator,
                        #mixturesign: multiply by -1 for fL1
                        mixturesign(analysis)*getattr(self, hypothesis.couplingname)
                       )

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

    def __repr__(self):
        couplings = ("g1", "g2", "g4", "g1prime2")
        if self.productionmode == "HJJ": couplings += ("ghg2", "ghg4")
        if self.productionmode == "ttH": couplings += ("kappa", "kappa_tilde")
        return "{}({}, {})".format(
                                   type(self).__name__,
                                   repr(self.productionmode.item.name),
                                   ", ".join(
                                             "{}={}".format(coupling, getattr(self, coupling))
                                             for coupling in couplings
                                            ),
                                  )

def samplewithfai(productionmode, analysis, fai, withdecay=False, productionmodeforfai=None):
    from combinehelpers import mixturesign, sigmaioversigma1
    analysis = Analysis(analysis)
    productionmode = ProductionMode(productionmode)
    if productionmodeforfai is None:
        productionmodeforfai = productionmode
    productionmodeforfai = ProductionMode(productionmodeforfai)
    if productionmodeforfai == "ggH":
        withdecay = True
    kwargs = {coupling: 0 for coupling in ("g1", "g2", "g4", "g1prime2")}
    power = (.25 if productionmodeforfai != "ggH" and withdecay else .5)
    kwargs["g1"] = (1-abs(fai))**power
    xsecratio = sigmaioversigma1(analysis, productionmodeforfai)
    if not withdecay:
        xsecratio /= sigmaioversigma1(analysis, "ggH")
    mixhypothesis = analysis.mixdecayhypothesis
    kwargs[analysis.couplingname] = copysign(
                                             (abs(fai) / xsecratio)**power,
                                             #mixturesign: multiply by -1 for fL1
                                             mixturesign(analysis)*fai
                                            )
    return ArbitraryCouplingsSample(productionmode, **kwargs)

class ReweightingSample(MultiEnum, SampleBase):
    __metaclass__ = MultiEnumABCMeta
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, Flavor]

    def __eq__(self, other):
        try:
            return MultiEnum.__eq__(self, other)
        except BaseException as e1:
            try:
                return SampleBase.__eq__(self, other)
            except BaseException as e2:
                raise ValueError(
                                 "Can't compare {} == {}!  Trying to compare as MultiEnums:\n"
                                 "  {}\n"
                                 "Trying to compare couplings:\n"
                                 "  {}"
                                 .format(self, other, e1, e2)
                                )

    def check(self, *args):
        if self.productionmode is None:
            raise ValueError("No option provided for productionmode\n{}".format(args))
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "HJJ", "ttH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hypothesis not in self.productionmode.validhypotheses:
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
        if self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX", "ttH", "HJJ") or self.alternategenerator == "POWHEG":
            return [self]
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in self.productionmode.validhypotheses]
        elif self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    def directreweightingsamples(self):
        if self.productionmode in ("ggZZ", "qqZZ", "VBF bkg", "ZX", "ttH", "HJJ") or self.alternategenerator == "POWHEG":
            return [self]
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in self.productionmode.directlyreweightedhypotheses]
        elif self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    def isbkg(self):
        return self.productionmode.isbkg

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH"):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH"):
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
        elif self.productionmode in ("VBF", "ZH", "WH", "WplusH", "WminusH"):
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
            elif self.hypothesis == "fL1proddec0.5":
                return "MC_weight_{}_g1g1prime2_proddec".format(self.productionmode)

            elif self.hypothesis == "fa2dec-0.5":
                return "MC_weight_{}_g1g2_dec_pi".format(self.productionmode)
            elif self.hypothesis == "fa2prod-0.5":
                return "MC_weight_{}_g1g2_prod_pi".format(self.productionmode)
            elif self.hypothesis == "fa2proddec0.5":
                return "MC_weight_{}_g1g2_proddec".format(self.productionmode)

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
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBFbkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH") or self.productionmode == "ZX" and not config.usedata:
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "data", "HJJ", "ttH", "WplusH", "WminusH"):
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH"):
            if self.hypothesis in (
                                   ["0+"]
                                   + ["{}{}0.5".format(a, b) for a in ("fa2", "fa3") for b in ("dec", "prod", "proddec-")]
                                   + ["{}{}0.5".format(a, b) for a in ("fL1",) for b in ("dec", "prod", "proddec")]
                                   + ["{}{}0.5".format(a, b) for a in ("fa2",) for b in ("dec-", "prod-", "proddec")]
                                  ):
                return 1
            if self.hypothesis in ("a2", "0-", "L1"):
                return 0
        raise self.ValueError("g1")

    @property
    def g2(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 0
        if self.hypothesis in (
                               ["0+", "0-", "L1"]
                             + ["{}{}0.5".format(a, b) for a in ("fa3",) for b in ("dec", "prod", "proddec-")]
                             + ["{}{}0.5".format(a, b) for a in ("fL1",) for b in ("dec", "prod", "proddec")]
                              ):
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

        if self.hypothesis == "fa2prod-0.5": return -ReweightingSample(self.productionmode, "fa2prod0.5").g2
        if self.hypothesis == "fa2dec-0.5": return -ReweightingSample(self.productionmode, "fa2dec0.5").g2
        if self.hypothesis == "fa2proddec0.5": return -ReweightingSample(self.productionmode, "fa2proddec-0.5").g2

        raise self.ValueError("g2")

    @property
    def g4(self):
        if self.productionmode in ("ttH", "HJJ"):
            return 0
        if self.hypothesis in (
                               ["0+", "a2", "L1"]
                             + ["{}{}0.5".format(a, b) for a in ("fa2",) for b in ("dec", "prod", "proddec-")]
                             + ["{}{}0.5".format(a, b) for a in ("fL1",) for b in ("dec", "prod", "proddec")]
                             + ["{}{}0.5".format(a, b) for a in ("fa2",) for b in ("dec-", "prod-", "proddec")]
                              ):
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
        if self.hypothesis in (
                               ["0+", "a2", "0-"]
                             + ["{}{}0.5".format(a, b) for a in ("fa2", "fa3") for b in ("dec", "prod", "proddec-")]
                             + ["{}{}0.5".format(a, b) for a in ("fa2",) for b in ("dec-", "prod-", "proddec")]
                              ):
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
            if self.hypothesis == "fL1proddec0.5":
                return -sqrt(constants.g1prime2VBF_gen*constants.g1prime2decay_gen)

        if self.productionmode == "ZH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2ZH_gen
            if self.hypothesis == "fL1proddec0.5":
                return -sqrt(constants.g1prime2ZH_gen*constants.g1prime2decay_gen)

        if self.productionmode == "WH":
            if self.hypothesis == "fL1prod0.5":
                return constants.g1prime2WH_gen
            if self.hypothesis == "fL1proddec0.5":
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
    enums = [ReweightingSample, Production, BlindStatus, AlternateGenerator]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))

        if self.blindstatus is None and self.productionmode == "data":
            raise ValueError("No blindstatus provided for data sample!\n{}".format(args))
        if self.blindstatus is not None and self.productionmode != "data":
            raise ValueError("blindstatus provided for MC sample!\n{}".format(args))

        if self.hypothesis is not None and self.hypothesis not in self.productionmode.generatedhypotheses:
            raise ValueError("No {} sample produced with hypothesis {}!\n{}".format(self.productionmode, self.hypothesis, args))

        if self.productionmode in ("WplusH", "WminusH") and self.alternategenerator != "POWHEG":
            raise ValueError("Separate {} sample is produced with POWHEG.  Maybe you meant to specify POWHEG, or WH?\n{}".format(self.productionmode, args))

        if (
            self.alternategenerator == "POWHEG"
            and (
                 self.hypothesis != "0+"
                 or self.productionmode not in ("VBF", "ZH", "WplusH", "WminusH", "ttH")
                )
           ):
            raise ValueError("No {} {} sample produced with {}\n{}".format(self.productionmode, self.hypothesis, self.alternategenerator, args))

        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
        if self.productionmode == "ggH":
            return self.production.CJLSTdir_anomalous()
        if self.productionmode == "VBF":
            return self.production.CJLSTdir_anomalous_VBF()
        if self.productionmode in ("ZH", "WH"):
            return self.production.CJLSTdir_anomalous_VH()
        if self.productionmode in ("data", "ZX"):
            return self.production.CJLSTdir_data()
        return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.alternategenerator == "POWHEG":
            if self.productionmode in ("VBF", "ZH", "WplusH", "WminusH", "ttH"):
                s = str(self.productionmode)
                if self.productionmode == "VBF": s = "VBFH"
                if self.hypothesis == "0+": return "{}125".format(s)
            raise self.ValueError("CJLSTdirname")
        if self.productionmode == "ggH" and self.production <= "160624":
            if self.hypothesis == "0+": return "0PM"
            if self.hypothesis == "a2": return "0PH"
            if self.hypothesis == "0-": return "0M"
            if self.hypothesis == "L1": return "0L1"
            if self.hypothesis == "fa20.5": return "0PHf05ph0"
            if self.hypothesis == "fa30.5": return "0Mf05ph0"
            if self.hypothesis == "fL10.5": return "0L1f05ph0"
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH") and self.production >= "160714":
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

class SampleBasis(MultiEnum):
    enums = [ProductionMode, Analysis]
    def __init__(self, hypotheses, *args):
        self.hypotheses = [Hypothesis(_) for _ in hypotheses]
        super(SampleBasis, self).__init__(*args)

    def check(self, *args):
        args = (self.hypotheses,)+args
        if self.productionmode == "ggH":
            dimension = 3
        elif self.productionmode in ("VBF", "WH", "ZH"):
            dimension = 5
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

        if len(self.hypotheses) != dimension:
            raise ValueError("Wrong number of hypotheses ({}, should be {})\n{}".format(len(self.hypotheses), dimension, args))
        if len(set(self.hypotheses)) != len(self.hypotheses):
            raise ValueError("Duplicate hypothesis\n{}".format(args))

    @property
    @cache
    def matrix(self):
        dimension = len(self.hypotheses)
        maxpower = dimension-1
        samples = [ReweightingSample(self.productionmode, _) for _ in self.hypotheses]
        return numpy.matrix(
                            [
                             [
                              sample.g1**(maxpower-i) * getattr(sample, self.analysis.couplingname)**i
                                  for i in range(dimension)
                             ]
                                 for sample in samples
                            ]
                           )
    @property
    @cache
    def invertedmatrix(self):
        return self.matrix.I


if __name__ == "__main__":
    for productionmode in "ggH", "VBF", "ZH":
        for analysis in "fa3", "fa2", "fL1":
            for fai in 0.5, -0.2:
                assert abs(samplewithfai(productionmode, analysis, fai).fai(productionmode, analysis)/fai - 1) < 1e-15
    f = ROOT.TFile(Sample("VBF", "0+", "160928").withdiscriminantsfile())
    t = f.candTree
    s = ReweightingSample("VBF", "fa2prod-0.5")
    length = min(t.GetEntries(), 100)
    for i, entry in enumerate(t, start=1):
        if s.get_MC_weight_function()(t) < 0:
            t.Show()
            assert False
        print i, "/", length
        if i == length:
            break
