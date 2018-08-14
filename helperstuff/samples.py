#!/usr/bin/env python
from abc import ABCMeta, abstractproperty
from collections import Counter
from itertools import combinations, permutations, product as cartesianproduct
from math import copysign, sqrt
import numpy
import os

import ROOT

import config
import constants
from enums import AlternateGenerator, analyses, Analysis, Extension, Flavor, flavors, purehypotheses, HffHypothesis, hffhypotheses, Hypothesis, MultiEnum, MultiEnumABCMeta, Production, ProductionMode, productions, PythiaSystematic, pythiasystematics
from utilities import cache, deprecate, mkdtemp, product, TFile, tlvfromptetaphim
from weightshelper import WeightsHelper


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
    def ghzgs1prime2(self):
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
    @abstractproperty
    def pdf(self):
        pass

    def ValueError(self, functionname):
        return ValueError("invalid sample {} in function {}".format(self, functionname))

    @property
    def nominalJHUxsec(self):
        result = constants.JHUXSHZZ2e2mua1.nominal_value
        if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g1**2
        elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g2**2 / constants.g2HZZ**2
        elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
            result *= self.g4**2 / constants.g4HZZ**2
        elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
            result *= self.g1prime2**2 / constants.g1prime2HZZ**2
        elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
            result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2HZZ**2
        else:
            raise ValueError("Nominal xsec is only defined for samples with no interference")

        if self.productionmode in ("ggH", "bbH", "tqH"):
            pass
        elif self.productionmode == "VBF":
            result *= constants.NNPDF30_lo_as_0130.JHUXSVBFa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2VBF**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4VBF**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2VBF**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2VBF**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "ZH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSZHa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2ZH**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4ZH**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2ZH**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2ZH**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "WH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSWHa1.nominal_value
            if self.g2 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g1**2
            elif self.g1 == self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g2**2 / constants.g2WH**2
            elif self.g1 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0:
                result *= self.g4**2 / constants.g4WH**2
            elif self.g1 == self.g2 == self.g4 == self.ghzgs1prime2 == 0:
                result *= self.g1prime2**2 / constants.g1prime2WH**2
            elif self.g1 == self.g2 == self.g4 == self.g1prime2 == 0:
                result *= self.ghzgs1prime2**2 / constants.ghzgs1prime2WH**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "HJJ":
            result *= constants.NNPDF30_lo_as_0130.JHUXSHJJa2.nominal_value
            if self.ghg4 == 0:
                result *= self.ghg2**2
            elif self.ghg2 == 0:
                result *= self.ghg4**2 / constants.ghg4HJJ**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        elif self.productionmode == "ttH":
            result *= constants.NNPDF30_lo_as_0130.JHUXSttHkappa.nominal_value
            if self.kappa_tilde == 0:
                result *= self.kappa**2
            elif self.kappa == 0:
                result *= self.kappa_tilde**2 / constants.kappa_tilde_ttH**2
            else:
                raise ValueError("Nominal xsec is only defined for samples with no interference")
        else:
            raise self.ValueError("nominalJHUxsec")

        if not (self.productionmode == "WH" and self.g1 == self.g2 == self.g4 == self.g1prime2 == 0) and result == 0 != self.JHUxsec:
            raise ValueError("Something is wrong")

        return result

    @property
    def JHUxsec(self):
        try:
            return self.JHUxsec_pdf(self.pdf)
        except ValueError as e:
            if "with no pdf" in str(e):
                raise ValueError("{} has no pdf, so can't get JHUxsec".format(self))
            raise

    def JHUxsec_pdf(self, pdf):
        from constants import JHUXSHZZ2e2mua1, JHUXSHZZ2e2mua2, JHUXSHZZ2e2mua3, JHUXSHZZ2e2muL1, JHUXSHZZ2e2muL1Zg,\
                              JHUXSHZZ2e2mua1a2, JHUXSHZZ2e2mua1a3, JHUXSHZZ2e2mua1L1, JHUXSHZZ2e2mua1L1Zg,         \
                              JHUXSHZZ2e2mua2a3, JHUXSHZZ2e2mua2L1, JHUXSHZZ2e2mua2L1Zg, JHUXSHZZ2e2mua3L1,         \
                              JHUXSHZZ2e2mua3L1Zg, JHUXSHZZ2e2muL1L1Zg

        if pdf is not None:
            exec 'from helperstuff.constants.'+pdf+' import ' + """                                            \
                              JHUXSVBFa1, JHUXSVBFa2, JHUXSVBFa3, JHUXSVBFL1, JHUXSVBFL1Zg,                    \
                              JHUXSVBFa1a2, JHUXSVBFa1a3, JHUXSVBFa1L1, JHUXSVBFa1L1Zg,                        \
                              JHUXSVBFa2a3, JHUXSVBFa2L1, JHUXSVBFa2L1Zg, JHUXSVBFa3L1,                        \
                              JHUXSVBFa3L1Zg, JHUXSVBFL1L1Zg,                                                  \
                              JHUXSZHa1, JHUXSZHa2, JHUXSZHa3, JHUXSZHL1, JHUXSZHL1Zg,                         \
                              JHUXSZHa1a2, JHUXSZHa1a3, JHUXSZHa1L1, JHUXSZHa1L1Zg,                            \
                              JHUXSZHa2a3, JHUXSZHa2L1, JHUXSZHa2L1Zg, JHUXSZHa3L1,                            \
                              JHUXSZHa3L1Zg, JHUXSZHL1L1Zg,                                                    \
                              JHUXSWHa1, JHUXSWHa2, JHUXSWHa3, JHUXSWHL1, JHUXSWHL1Zg,                         \
                              JHUXSWHa1a2, JHUXSWHa1a3, JHUXSWHa1L1, JHUXSWHa1L1Zg,                            \
                              JHUXSWHa2a3, JHUXSWHa2L1, JHUXSWHa2L1Zg, JHUXSWHa3L1,                            \
                              JHUXSWHa3L1Zg, JHUXSWHL1L1Zg,                                                    \
                              JHUXSHJJa2, JHUXSHJJa3, JHUXSHJJa2a3,                                            \
                              JHUXSttHkappa, JHUXSttHkappatilde, JHUXSttHkappakappatilde
            """
        else:
            if self.productionmode not in ("ggH", "bbH"):
                raise ValueError("Can't get xsec for {} with no pdf".format(self.productionmode))

        if self.productionmode in ("ggH", "bbH"):
            return (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "VBF":
            return (
                      JHUXSVBFa1*self.g1**2
                    + JHUXSVBFa2*self.g2**2
                    + JHUXSVBFa3*self.g4**2
                    + JHUXSVBFL1*self.g1prime2**2
                    + JHUXSVBFL1Zg*self.ghzgs1prime2**2

                    + JHUXSVBFa1a2 * self.g1*self.g2
                    + JHUXSVBFa1a3 * self.g1*self.g4
                    + JHUXSVBFa1L1 * self.g1*self.g1prime2
                    + JHUXSVBFa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSVBFa2a3 * self.g4*self.g2
                    + JHUXSVBFa2L1 * self.g2*self.g1prime2
                    + JHUXSVBFa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSVBFa3L1 * self.g4*self.g1prime2
                    + JHUXSVBFa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSVBFL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "ZH":
            return (
                      JHUXSZHa1*self.g1**2
                    + JHUXSZHa2*self.g2**2
                    + JHUXSZHa3*self.g4**2
                    + JHUXSZHL1*self.g1prime2**2
                    + JHUXSZHL1Zg*self.ghzgs1prime2**2

                    + JHUXSZHa1a2 * self.g1*self.g2
                    + JHUXSZHa1a3 * self.g1*self.g4
                    + JHUXSZHa1L1 * self.g1*self.g1prime2
                    + JHUXSZHa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSZHa2a3 * self.g4*self.g2
                    + JHUXSZHa2L1 * self.g2*self.g1prime2
                    + JHUXSZHa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSZHa3L1 * self.g4*self.g1prime2
                    + JHUXSZHa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSZHL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "WH":
            return (
                      JHUXSWHa1*self.g1**2
                    + JHUXSWHa2*self.g2**2
                    + JHUXSWHa3*self.g4**2
                    + JHUXSWHL1*self.g1prime2**2
                    + JHUXSWHL1Zg*self.ghzgs1prime2**2

                    + JHUXSWHa1a2 * self.g1*self.g2
                    + JHUXSWHa1a3 * self.g1*self.g4
                    + JHUXSWHa1L1 * self.g1*self.g1prime2
                    + JHUXSWHa1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSWHa2a3 * self.g4*self.g2
                    + JHUXSWHa2L1 * self.g2*self.g1prime2
                    + JHUXSWHa2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSWHa3L1 * self.g4*self.g1prime2
                    + JHUXSWHa3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSWHL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )
        if self.productionmode == "HJJ":
            return (
                      JHUXSHJJa2 * self.ghg2**2
                    + JHUXSHJJa3 * self.ghg4**2
                    + JHUXSHJJa2a3 * self.ghg2*self.ghg4
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )
        if self.productionmode == "ttH":
            return (
                      JHUXSttHkappa * self.kappa**2
                    + JHUXSttHkappatilde * self.kappa_tilde**2
                    + JHUXSttHkappakappatilde * self.kappa*self.kappa_tilde
                   ) * (
                      JHUXSHZZ2e2mua1*self.g1**2
                    + JHUXSHZZ2e2mua2*self.g2**2
                    + JHUXSHZZ2e2mua3*self.g4**2
                    + JHUXSHZZ2e2muL1*self.g1prime2**2
                    + JHUXSHZZ2e2muL1Zg*self.ghzgs1prime2**2

                    + JHUXSHZZ2e2mua1a2 * self.g1*self.g2
                    + JHUXSHZZ2e2mua1a3 * self.g1*self.g4
                    + JHUXSHZZ2e2mua1L1 * self.g1*self.g1prime2
                    + JHUXSHZZ2e2mua1L1Zg * self.g1*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua2a3 * self.g4*self.g2
                    + JHUXSHZZ2e2mua2L1 * self.g2*self.g1prime2
                    + JHUXSHZZ2e2mua2L1Zg * self.g2*self.ghzgs1prime2

                    + JHUXSHZZ2e2mua3L1 * self.g4*self.g1prime2
                    + JHUXSHZZ2e2mua3L1Zg * self.g4*self.ghzgs1prime2

                    + JHUXSHZZ2e2muL1L1Zg * self.g1prime2*self.ghzgs1prime2
                   )

        if self.productionmode == "tqH":
            assert self.kappa_tilde == self.g4 == self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0, self
            return 1

        assert False

    @property
    def xsec(self):
        if self.productionmode == "ggH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSggH2e2mu
            return constants.SMXSggH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "VBF":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSVBF2e2mu
            return constants.SMXSVBF2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "ZH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSZH2e2mu
            return constants.SMXSZH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "WH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWH2e2mu
            return constants.SMXSWH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "WplusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWpH2e2mu
        if self.productionmode == "WminusH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSWmH2e2mu
        if self.productionmode == "HJJ":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return 1
            return constants.SMXSHJJ2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "ttH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM" and self.hffhypothesis == "Hff0+": return constants.SMXSttH2e2mu
            return constants.SMXSttH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "bbH":
            if hasattr(self, "hypothesis") and self.hypothesis == "SM": return constants.SMXSbbH2e2mu
            return constants.SMXSbbH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+").JHUxsec_pdf(self.pdf)
        if self.productionmode == "tqH":
            if hasattr(self, hypothesis) and self.hypothesis == "SM": return constants.SMXStqH2e2mu
            return constants.SMXStqH2e2mu * self.JHUxsec / ReweightingSample(self.productionmode, "0+", "Hff0+").JHUxsec_pdf(self.pdf)
        raise self.ValueError("xsec")

    @property
    def SMxsec(self):
        hffhypothesis = "Hff0+" if self.hffhypothesis else None
        return ReweightingSample(self.productionmode, "SM", hffhypothesis).xsec

    def __eq__(self, other):
        return all((
                   self.productionmode == other.productionmode,
                   self.g1 == other.g1,
                   self.g2 == other.g2,
                   self.g4 == other.g4,
                   self.g1prime2 == other.g1prime2,
                   self.ghzgs1prime2 == other.ghzgs1prime2,
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
               ) and (
                   self.pdf is None or other.pdf is None or self.pdf == other.pdf
               )
    def __ne__(self, other):
        return not self == other

    def linearcombinationofreweightingsamples(self, hypotheses):
        g1 = self.g1
        if   self.g2 == self.g1prime2 == self.ghzgs1prime2 == 0: analysis = Analysis("fa3");   gi = self.g4
        elif self.g4 == self.g1prime2 == self.ghzgs1prime2 == 0: analysis = Analysis("fa2");   gi = self.g2
        elif self.g2 == self.g4       == self.ghzgs1prime2 == 0: analysis = Analysis("fL1");   gi = self.g1prime2
        elif self.g2 == self.g4       == self.g1prime2     == 0: analysis = Analysis("fL1Zg"); gi = self.ghzgs1prime2
        else: assert False

        if hypotheses == "templates":
            from templates import TemplatesFile
            templatesfile = TemplatesFile(analysis, str(self.productionmode).lower(), "2e2mu", "Untagged", config.productionsforcombine[0])
            hypotheses = [_.hypothesis for _ in templatesfile.templates()]

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
    def MC_weight_terms(self):
        factors = []

        weightshelper = WeightsHelper(self)

        for prodordec in "prod", "dec":
          if weightshelper.useproddec(prodordec):
            counter = Counter()
            for couplingname, couplingvalue, weight in weightshelper.couplingsandweights(prodordec, mix=False):
              coupling = getattr(self, couplingname) / float(couplingvalue)
              counter[weight] += coupling**2
            for ((coupling1name, coupling1value, weight1),
                 (coupling2name, coupling2value, weight2),
                 weightint) in weightshelper.couplingsandweights(prodordec, mix=True):
              coupling1 = getattr(self, coupling1name ) / float(coupling1value)
              coupling2 = getattr(self, coupling2name ) / float(coupling2value)
              counter[weight1]   -= coupling1*coupling2
              counter[weight2]   -= coupling1*coupling2
              counter[weightint] += coupling1*coupling2

            factors.append([
                            (weightname, couplingsq)
                                  for weightname, couplingsq in counter.iteritems()
                                   if couplingsq
                                         and weightname is not None
                           ])

        assert factors
        return factors

    @property
    def MC_weight(self):
        result = "*".join(
                          "("+
                          "+".join(
                                   "({}*{})".format(weightname, couplingsq)
                                         for weightname, couplingsq in factor
                                  )
                          +")" for factor in self.MC_weight_terms
                         )
        return result

    def get_MC_weight_function(self_sample, functionname=None, reweightingonly=False):
        if functionname is None:
            functionname = self_sample.weightname()

        strsample = str(self_sample)

        if not config.LHE:

            if self_sample.productionmode == "ggZZ":
                def MC_weight_function(self_tree):
                    KFactor = self_tree.KFactor_QCD_ggZZ_Nominal
                    return self_tree.overallEventWeight * self_tree.xsec * KFactor / self_tree.nevents

            elif self_sample.productionmode == "qqZZ":
                def MC_weight_function(self_tree):
                    if hasattr(self_tree, "KFactor_EW_qqZZ"): #qqZZ->itself
                        KFactor = self_tree.KFactor_EW_qqZZ * self_tree.KFactor_QCD_qqZZ_M
                        return self_tree.overallEventWeight * self_tree.xsec * KFactor #/ self_tree.nevents
                                                                                       #Do NOT divide by nevents.  This will happen
                                                                                       #in creating the templates when we divide by
                                                                                       #sum(effectiveentries[==nevents])
                    else:                                     #VBFbkg->qqZZ
                        result = self_tree.p_Gen_JJQCD_BKG_MCFM
                        if not reweightingonly and result != 0:
                            cutoff = self_tree.cutoffs[strsample]
                            if result > cutoff:
                                result = cutoff**2 / result
                            result *= (
                                         self_tree.overallEventWeight
                                       * self_tree.multiplyweight[strsample]
                                      )
                        return result


            elif self_sample.productionmode == "ZX":
                import ZX
                def MC_weight_function(self_tree):
                    LepPt, LepEta, LepLepId = self_tree.tree.LepPt, self_tree.tree.LepEta, self_tree.tree.LepLepId
                    return ZX.fakeRate13TeV(LepPt[2],LepEta[2],LepLepId[2]) * ZX.fakeRate13TeV(LepPt[3],LepEta[3],LepLepId[3])

            elif (self_sample.productionmode in ("VBF bkg", "tqH")
                       or hasattr(self_sample, "alternategenerator") and self_sample.alternategenerator in ("POWHEG", "MINLO", "NNLOPS")
                       or hasattr(self_sample, "pythiasystematic") and self_sample.pythiasystematic is not None):
                def MC_weight_function(self_tree):
                    if reweightingonly: return 1
                    return self_tree.overallEventWeight * self_tree.xsec / self_tree.nevents

            elif self_sample.issignal:

                factors = self_sample.MC_weight_terms

                photoncut_decay = False
                photoncut_ZH = False and self_sample.productionmode == "ZH"
                photoncut_VBF = False and self_sample.productionmode == "VBF"

                def MC_weight_function(self_tree):
                    if photoncut_decay:
                        leptons = [(id, tlvfromptetaphim(pt, eta, phi, m)) for id, pt, eta, phi, m in zip(self_tree.LHEDaughterId, self_tree.LHEDaughterPt, self_tree.LHEDaughterEta, self_tree.LHEDaughterPhi, self_tree.LHEDaughterMass)]
                        for (id1, p1), (id2, p2) in cartesianproduct(leptons, leptons):
                            if id1 == -id2 and (p1+p2).M() < 4: return 0
                    if photoncut_ZH:
                        Zdecay = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(self_tree.LHEAssociatedParticlePt, self_tree.LHEAssociatedParticleEta, self_tree.LHEAssociatedParticlePhi, self_tree.LHEAssociatedParticleMass)]
                        assert len(Zdecay) == 2
                        if (Zdecay[0]+Zdecay[1]).M() < 4: return 0
                    if photoncut_VBF:
                        jets = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(self_tree.LHEAssociatedParticlePt, self_tree.LHEAssociatedParticleEta, self_tree.LHEAssociatedParticlePhi, self_tree.LHEAssociatedParticleMass)]
                        assert len(jets) == 2
                        if any(jet.Pt() < 15 for jet in jets): return 0
                    result = product(
                                     sum(
                                         getattr(self_tree, weightname) * couplingsq
                                               for weightname, couplingsq in factor
                                        ) for factor in factors
                                    )

                    if not reweightingonly and result != 0:
                        cutoff = self_tree.cutoffs[strsample]
                        if result > cutoff:
                            result = cutoff**2 / result
                        result *= (
                                     self_tree.overallEventWeight
                                   * self_tree.multiplyweight[strsample]
                                  )
                    return result

            else:
                raise ValueError("No MC weight function defined for {}".format(self_sample))

        else: #config.LHE
            photoncut_decay = False
            if self_sample.productionmode == "ggH":
                def MC_weight_function(self_tree):
                    if photoncut_decay:
                        leptons = [(id, tlv) for id, tlv in self_tree.event.gendaughters]
                        for (id1, p1), (id2, p2) in cartesianproduct(leptons, leptons):
                            if id1 == -id2 and (p1+p2).M() < 4: return 0
                    result = self_tree.event.weight
                    if not reweightingonly and result != 0:
                        result *= self_sample.SMxsec / self_tree.sumofweights
                    return result
            elif self_sample.productionmode == "qqZZ":
                xsec = 1.256 * 1000
                def MC_weight_function(self_tree):
                    result = self_tree.event.weight
                    if not reweightingonly and result != 0:
                        result *= xsec / self_tree.sumofweights
                    return result
            else:
                raise ValueError("No MC weight function defined for {}".format(self_sample))

        MC_weight_function.__name__ = functionname
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
            if not getattr(self, h.couplingname) and h != hypothesis: continue
            kwargs = {_.couplingname: 0 for _ in purehypotheses}
            kwargs["pdf"] = "NNPDF30_lo_as_0130"
            kwargs[h.couplingname] = 1
            term = 0
            for productionmode in productionmodes:
                term += getattr(self, h.couplingname)**power * ArbitraryCouplingsSample(productionmode, **kwargs).nominalJHUxsec \
                         * (constants.nominal_normalize_WH_to_ZH if productionmode == "WH" else 1)
            if not withdecay:
                term /= ArbitraryCouplingsSample("ggH", **kwargs).nominalJHUxsec
            denominator += term
            if h == hypothesis:
                numerator = term
        return copysign(
                        numerator / denominator,
                        #mixturesign: multiply by -1 for fL1
                        mixturesign(analysis)*getattr(self, hypothesis.couplingname)*self.g1
                       )

    @property
    def a1eLeR(self):
        assert self.g2 == self.g4 == 0
        ghz1 = self.g1 + 2 * self.g1prime2 * (constants.M_Z**2 - 1j*constants.M_Z*constants.Ga_Z)/constants.L1**2
        ghz1 = ghz1.real  #<-- !!!
        ghzzp1 = constants.M_Z**2/constants.L1**2
        ezp_L = constants.aL * self.g1prime2 + constants.e * self.ghzgs1prime2
        ezp_R = constants.aR * self.g1prime2 + constants.e * self.ghzgs1prime2
        return ghz1, ghzzp1*ezp_L, ghzzp1*ezp_R

    @property
    def feL(self):
        a1, eL, eR = self.a1eLeR
        return (1 if a1*eL>0 else -1) * eL**2 * constants.JHUXSHZZ2e2mueL / (eL**2 * constants.JHUXSHZZ2e2mueL + eR**2 * constants.JHUXSHZZ2e2mueR + abs(a1)**2 * constants.JHUXSHZZ2e2mua1)

    @property
    def feR(self):
        a1, eL, eR = self.a1eLeR
        return (1 if a1*eR>0 else -1) * eR**2 * constants.JHUXSHZZ2e2mueR / (eL**2 * constants.JHUXSHZZ2e2mueL + eR**2 * constants.JHUXSHZZ2e2mueR + abs(a1)**2 * constants.JHUXSHZZ2e2mua1)

class ArbitraryCouplingsSample(SampleBase):
    def __init__(self, productionmode, g1, g2, g4, g1prime2, ghzgs1prime2, ghg2=None, ghg4=None, kappa=None, kappa_tilde=None, pdf=None):
        self.productionmode = ProductionMode(productionmode)
        self.__g1, self.__g2, self.__g4, self.__g1prime2, self.__ghzgs1prime2 = g1, g2, g4, g1prime2, ghzgs1prime2
        if self.productionmode not in ("ggH", "VBF", "ZH", "WH", "HJJ", "ttH", "bbH"):
            raise ValueError("Bad productionmode {}".format(self.productionmode))
        if sum(bool(g) for g in (g2, g4, g1prime2 or ghzgs1prime2)) > 1:
            raise ValueError("Can only set at most one of g2, g4, or g1prime2 and/or ghzgs1prime2")

        if self.productionmode == "HJJ":
            self.__ghg2, self.__ghg4 = ghg2, ghg4
            if ghg2 is None or ghg4 is None:
                raise ValueError("Have to set ghg2 and ghg4 for HJJ")
        elif self.productionmode == "ggH":
            if ghg2 is None:
                ghg2 = 1
            if ghg4 is None:
                ghg4 = 0
            if ghg2 != 1:
                raise ValueError("ghg2 has to = 1 for {}".format(self.productionmode))
            if ghg4 != 0:
                raise ValueError("ghg4 has to = 0 for {}".format(self.productionmode))
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

        self.__pdf = pdf
        if pdf is not None:
            try:
                exec "from constants import "+pdf in globals(), locals()
            except ImportError:
                raise ValueError("Unknown pdf "+pdf)

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
    def ghzgs1prime2(self):
        return self.__ghzgs1prime2
    @property
    def ghg2(self):
        if self.__ghg2 is None: raise ValueError("{} doesn't have ghg2!".format(self))
        return self.__ghg2
    @property
    def ghg4(self):
        if self.__ghg4 is None: raise ValueError("{} doesn't have ghg4!".format(self))
        return self.__ghg4
    @property
    def kappa(self):
        return self.__kappa
    @property
    def kappa_tilde(self):
        return self.__kappa_tilde
    @property
    def pdf(self):
        return self.__pdf

    def __repr__(self):
        kwargs = ()
        if self.productionmode == "HJJ": kwargs += ("ghg2", "ghg4")
        if self.productionmode == "ttH": kwargs += ("kappa", "kappa_tilde")
        if self.pdf is not None: kwargs += ("pdf",)
        kwargs += ("g1", "g2", "g4", "g1prime2", "ghzgs1prime2")
        return "{}({}, {})".format(
                                   type(self).__name__,
                                   repr(self.productionmode.item.name),
                                   ", ".join(
                                             "{}={}".format(kwarg, getattr(self, kwarg))
                                             for kwarg in kwargs
                                            ),
                                  )

def samplewithfai(productionmode, analysis, fai, withdecay=False, productionmodeforfai=None, pdf=None):
    from combinehelpers import mixturesign, sigmaioversigma1
    analysis = Analysis(analysis)
    productionmode = ProductionMode(productionmode)

    if isinstance(withdecay, (basestring, ProductionMode)) and productionmodeforfai is None:
        raise TypeError(
                        "withdecay is a {thetype}!  Maybe you meant\n"
                        "  samplewithfai({productionmode!r}, {analysis!r}, {fai!r}, productionmodeforfai={withdecay!r})\n"
                        "or\n"
                        "  samplewithfai({productionmode!r}, {analysis!r}, {fai!r}, bool({withdecay!r}))"
                        .format(thetype=type(withdecay).__name__, productionmode=productionmode, analysis=analysis, fai=fai, withdecay=withdecay)
                       )

    if productionmodeforfai is None:
        productionmodeforfai = productionmode

    if productionmodeforfai == "ggH":
        withdecay = True

    if str(productionmodeforfai) == "VH":
        pass
    else:
        productionmodeforfai = ProductionMode(productionmodeforfai)

    kwargs = {coupling: 0 for coupling in ("g1", "g2", "g4", "g1prime2", "ghzgs1prime2")}
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
    kwargs["pdf"] = pdf
    return ArbitraryCouplingsSample(productionmode, **kwargs)

def samplewitha1eLeR(a1, eL, eR):
    kwargs = {
      "g1": a1 - 2/(constants.aL - constants.aR) * (eL - eR),
      "g1prime2": constants.L1**2 / (constants.M_Z**2 * (constants.aL - constants.aR)) * (eL-eR),
      "ghzgs1prime2": -constants.L1**2 / (constants.M_Z**2 * constants.e * (constants.aL - constants.aR)) * (constants.aR*eL-constants.aL*eR),
      "g2": 0,
      "g4": 0,
    }
    result = ArbitraryCouplingsSample("ggH", **kwargs)
    return result

def samplewithfeLfeR(feL, feR):
    fa1 = 1-abs(feL)-abs(feR)
    a1 = sqrt(fa1)
    eL = (1 if feL>0 else -1) * sqrt(abs(feL) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2mueL)
    eR = (1 if feR>0 else -1) * sqrt(abs(feR) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2mueR)
    return samplewitha1eLeR(a1, eL, eR)

def samplewithfL1fL1Zg(fL1, fL1Zg):
    fa1 = 1-abs(fL1)-abs(fL1Zg)
    a1 = sqrt(fa1)
    L1 = -(1 if fL1>0 else -1) * sqrt(abs(fL1) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2muL1)
    L1Zg = -(1 if fL1Zg>0 else -1) * sqrt(abs(fL1Zg) * constants.JHUXSHZZ2e2mua1 / constants.JHUXSHZZ2e2muL1Zg)
    return ArbitraryCouplingsSample("ggH", g1=a1, g1prime2=L1, ghzgs1prime2=L1Zg, g2=0, g4=0)

class ReweightingSample(MultiEnum, SampleBase):
    __metaclass__ = MultiEnumABCMeta
    enumname = "reweightingsample"
    enums = [ProductionMode, Hypothesis, HffHypothesis]

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
        elif self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "HJJ", "ttH", "bbH", "tqH"):
            if self.hypothesis is None:
                raise ValueError("No hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hypothesis not in self.productionmode.validhypotheses:
                raise ValueError("{} hypothesis can't be {}\n{}".format(self.productionmode, self.hypothesis, args))
            if self.productionmode in ("HJJ", "ttH", "tqH"):
                if self.hffhypothesis is None:
                    raise ValueError("Hff hypothesis not provided for {} productionmode\n{}".format(self.productionmode, args))
                if self.productionmode == "tqH" and self.hffhypothesis != "Hff0+":
                    raise ValueError("Bad hffhypothesis {} for {} productionmode\n{}".format(self.hffhypothesis, self.productionmode, args))
            else:
                if self.hffhypothesis is not None:
                    raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ggZZ", "VBF bkg"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode == "qqZZ":
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        elif self.productionmode in ("ZX", "data"):
            if self.hypothesis is not None:
                raise ValueError("Hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.hffhypothesis is not None:
                raise ValueError("Hff hypothesis provided for {} productionmode\n{}".format(self.productionmode, args))
        else:
            raise ValueError("Bad productionmode {}\n{}".format(self.productionmode, args))

    def reweightingsamples(self):
        if self.productionmode == "qqZZ":      #self == ReweightingSample("qqZZ") or ReweightingSamplePlus("qqZZ", "ext")
            return [ReweightingSample("qqZZ")]
        if self.productionmode in ("ggZZ", "qqZZ", "ZX"):
            return [self]
        if self.alternategenerator in ("POWHEG", "MINLO", "NNLOPS") or self.pythiasystematic is not None:
            return [ReweightingSamplePlus(self.reweightingsample, self.alternategenerator, self.pythiasystematic)] #but not self.extension
        if self.productionmode == "VBF bkg":
            return [self, ReweightingSample("qqZZ")]
        if self.productionmode == "ttH":  #ttH reweighting doesn't work unfortunately
            return [
                    ReweightingSample(self.productionmode, hypothesis, hffhypothesis)
                            for hypothesis in self.productionmode.validhypotheses
                            for hffhypothesis in hffhypotheses
                            if hffhypothesis == self.hffhypothesis
                   ]
        if self.productionmode in ("ttH", "HJJ"):
            return [
                    ReweightingSample(self.productionmode, hypothesis, hffhypothesis)
                            for hypothesis in self.productionmode.validhypotheses
                            for hffhypothesis in hffhypotheses
                   ]
        if self.productionmode == "ggH" and config.LHE:
            return [self]
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "bbH"):
            return [ReweightingSample(self.productionmode, hypothesis) for hypothesis in self.productionmode.validhypotheses
                     if not deprecate(self.productionmode == "bbH" and hypothesis in ("fa30.5fa20.5", "fa20.5fL10.5", "fa30.5fL10.5", "fa30.5fL1Zg0.5", "fL10.5fL1Zg0.5", "fa20.5fL1Zg0.5"), 2018, 9, 1)]
        if self.productionmode == "tqH":
            return [self]
        if self.productionmode == "data":
            return []
        raise self.ValueError("reweightingsamples")

    @property
    def isbkg(self):
        return self.productionmode.isbkg
    @property
    def issignal(self):
        return not self.isbkg

    def isZX(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH", "bbH", "tqH"):
            return False
        elif self.productionmode == "ZX":
            return True
        raise self.ValueError("isZX")

    def isdata(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH", "bbH", "tqH"):
            return False
        elif self.productionmode == "data":
            return True
        raise self.ValueError("isdata")

    def weightname(self):
        if self.productionmode in ("ggH", "bbH"):
            if self.hypothesis == "0+":
                return "MC_weight_{}_g1".format(self.productionmode)
            elif self.hypothesis == "a2":
                return "MC_weight_{}_g2".format(self.productionmode)
            elif self.hypothesis == "0-":
                return "MC_weight_{}_g4".format(self.productionmode)
            elif self.hypothesis == "L1":
                return "MC_weight_{}_g1prime2".format(self.productionmode)
            elif self.hypothesis == "L1Zg":
                return "MC_weight_{}_ghzgs1prime2".format(self.productionmode)
            elif self.hypothesis == "fa20.5":
                return "MC_weight_{}_g1g2".format(self.productionmode)
            elif self.hypothesis == "fa2dec-0.9":
                return "MC_weight_{}_g1g2_minuspoint9".format(self.productionmode)
            elif self.hypothesis == "fa30.5":
                return "MC_weight_{}_g1g4".format(self.productionmode)
            elif self.hypothesis == "fL10.5":
                return "MC_weight_{}_g1g1prime2".format(self.productionmode)
            elif self.hypothesis == "fL1Zg0.5":
                return "MC_weight_{}_g1ghzgs1prime2".format(self.productionmode)
            elif self.hypothesis == "fL10.5fL1Zg0.5":
                return "MC_weight_{}_g1prime2ghzgs1prime2".format(self.productionmode)
            for a in "fa3", "fa2", "fL1", "fL1Zg":
                if self.hypothesis == "{}-0.5".format(a):
                    return ReweightingSample(self.productionmode, "{}0.5".format(a)).weightname()+"_pi"
            for a, b in combinations(("fa3", "fa2", "fL1", "fL1Zg"), 2):
                if self.hypothesis == a+"dec0.5"+b+"dec0.5":
                    return "MC_weight_{}_{}{}".format(self.productionmode, Analysis(a).couplingname, Analysis(b).couplingname)

        elif self.productionmode in ("VBF", "ZH", "WH", "WplusH", "WminusH"):
            if self.hypothesis == "0+":
                return "MC_weight_{}_g1".format(self.productionmode)
            elif self.hypothesis == "a2":
                return "MC_weight_{}_g2".format(self.productionmode)
            elif self.hypothesis == "0-":
                return "MC_weight_{}_g4".format(self.productionmode)
            elif self.hypothesis == "L1":
                return "MC_weight_{}_g1prime2".format(self.productionmode)
            elif self.hypothesis == "L1Zg":
                return "MC_weight_{}_ghzgs1prime2".format(self.productionmode)

            elif self.hypothesis == "fa2dec0.5":
                return "MC_weight_{}_g1g2_dec".format(self.productionmode)
            elif self.hypothesis == "fa3dec0.5":
                return "MC_weight_{}_g1g4_dec".format(self.productionmode)
            elif self.hypothesis == "fL1dec0.5":
                return "MC_weight_{}_g1g1prime2_dec".format(self.productionmode)
            elif self.hypothesis == "fL1Zgdec0.5":
                return "MC_weight_{}_g1ghzgs1prime2_dec".format(self.productionmode)

            elif self.hypothesis == "fa2prod0.5":
                return "MC_weight_{}_g1g2_prod".format(self.productionmode)
            elif self.hypothesis == "fa3prod0.5":
                return "MC_weight_{}_g1g4_prod".format(self.productionmode)
            elif self.hypothesis == "fL1prod0.5":
                return "MC_weight_{}_g1g1prime2_prod".format(self.productionmode)
            elif self.hypothesis == "fL1Zgprod0.5":
                return "MC_weight_{}_g1ghzgs1prime2_prod".format(self.productionmode)

            elif self.hypothesis == "fa2proddec0.5":
                return "MC_weight_{}_g1g2_proddec".format(self.productionmode)
            elif self.hypothesis == "fa3proddec0.5":
                return "MC_weight_{}_g1g4_proddec".format(self.productionmode)
            elif self.hypothesis == "fL1proddec0.5":
                return "MC_weight_{}_g1g1prime2_proddec".format(self.productionmode)
            elif self.hypothesis == "fL1Zgproddec0.5":
                return "MC_weight_{}_g1ghzgs1prime2_proddec".format(self.productionmode)

            elif self.hypothesis == "fa2dec-0.9":
                return "MC_weight_{}_g1g2_dec_minuspoint9".format(self.productionmode)

            for a in "fa3", "fa2", "fL1", "fL1Zg":
                for b in "prod", "dec", "proddec":
                    if self.hypothesis == "{}{}-0.5".format(a, b):
                        return ReweightingSample(self.productionmode, "{}{}0.5".format(a, b)).weightname()+"_pi"

            for proddec in "prod", "dec", "proddec":
                sign = "-" if proddec == "proddec" else ""
                for a, b in combinations(("fa3", "fa2", "fL1", "fL1Zg"), 2):
                    hypsuffix = wtsuffix = ""
                    if sign == "-": wtsuffix += "_pi"
                    if self.hypothesis == a+proddec+"0.5"+b+proddec+sign+"0.5"+hypsuffix:
                        return "MC_weight_{}_{}{}_{}".format(
                          self.productionmode,
                          Analysis(a).couplingname,
                          Analysis(b).couplingname,
                          proddec,
                        ) + wtsuffix
                    if self.hypothesis == a+proddec+"0.33"+b+proddec+sign+"0.33"+hypsuffix:
                        return "MC_weight_{}_g1{}{}_{}".format(
                          self.productionmode,
                          Analysis(a).couplingname,
                          Analysis(b).couplingname,
                          proddec,
                        ) + wtsuffix
                for a, b, c in combinations(("fa3", "fa2", "fL1", "fL1Zg"), 3):
                    hypsuffix = wtsuffix = ""
                    if sign == "-": wtsuffix += "_pi"
                    if self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+sign+"0.33"+hypsuffix:
                        return "MC_weight_{}_{}{}{}_{}".format(
                          self.productionmode,
                          Analysis(a).couplingname,
                          Analysis(b).couplingname,
                          Analysis(c).couplingname,
                          proddec,
                        ) + wtsuffix
                    if proddec != "proddec": continue
                    wtsuffix = hypsuffix
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+hypsuffix:
                        return "MC_weight_{}_g1{}{}{}_{}".format(
                          self.productionmode,
                          Analysis(a).couplingname,
                          Analysis(b).couplingname,
                          Analysis(c).couplingname,
                          proddec,
                        ) + wtsuffix
                for a, b, c, d in combinations(("fa3", "fa2", "fL1", "fL1Zg"), 4):
                    suffix = ""
                    if proddec != "proddec": continue
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25"+suffix:
                        return "MC_weight_{}_{}{}{}{}_{}".format(
                          self.productionmode,
                          Analysis(a).couplingname,
                          Analysis(b).couplingname,
                          Analysis(c).couplingname,
                          Analysis(d).couplingname,
                          proddec,
                        ) + suffix

        elif self.productionmode in ("HJJ", "ttH"):
            if self.productionmode == "HJJ":
                if self.hffhypothesis == "Hff0+":
                    result = "MC_weight_HJJ_ghg2"
                elif self.hffhypothesis == "Hff0-":
                    result = "MC_weight_HJJ_ghg4"
                elif self.hffhypothesis == "fCP0.5":
                    result = "MC_weight_HJJ_ghg2ghg4"
            elif self.productionmode == "ttH":
                if self.hffhypothesis == "Hff0+":
                    result = "MC_weight_ttH_kappa"
                elif self.hffhypothesis == "Hff0-":
                    result = "MC_weight_ttH_kappatilde"
                elif self.hffhypothesis == "fCP0.5":
                    result = "MC_weight_ttH_kappakappatilde"
            result += ReweightingSample("ggH", self.hypothesis).weightname().replace("MC_weight_ggH", "")
            return result

        elif self.productionmode == "tqH":
            if self.hypothesis == "0+" and self.hffhypothesis == "Hff0+":
                return "MC_weight_tqH_kappa_g1"

        elif self.productionmode == "ggZZ":
            return "MC_weight_ggZZ"
        elif self.productionmode == "qqZZ":
            return "MC_weight_qqZZ"
        elif self.productionmode == "VBF bkg":
            return "MC_weight_VBFbkg"
        elif self.productionmode == "ZX":
            return "MC_weight_ZX"
        raise self.ValueError("weightname")

    @property
    def MC_weight(self):
        return self.weightname()

    def TDirectoryname(self):
        if self.productionmode in ("ggH", "ggZZ", "qqZZ", "VBFbkg", "data", "VBF", "ZH", "WH", "ttH", "HJJ", "WplusH", "WminusH", "bbH", "tqH") or self.productionmode == "ZX" and not config.usedata:
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ", "VBF bkg", "ZX", "data", "HJJ", "ttH", "WplusH", "WminusH", "bbH", "tqH"):
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
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "HJJ", "bbH", "tqH"):
            if self.hypothesis == "0+": return 1
            if self.hypothesis in ("0-", "a2", "L1", "L1Zg"): return 0
            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                minus = "-" if proddec == "proddec" else ""
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        return 1
                    if self.hypothesis == a+proddec+"-0.5":
                        return 1
                for a, b in combinations(couplings, 2):
                    if self.hypothesis == a+proddec+"0.5"+b+proddec+minus+"0.5":
                        return 0
                    if self.hypothesis == a+proddec+"0.33"+b+proddec+minus+"0.33":
                        return 1
                for a, b, c in combinations(couplings, 3):
                    if self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33":
                        return 0
                    if proddec != "proddec": continue
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25":
                        return 1
                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        return 0

            if self.hypothesis == "fa2dec-0.9": return 1

        raise self.ValueError("g1")

    @property
    def g2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "HJJ"):
            if self.hypothesis == "a2": return constants.g2HZZ if self.productionmode in ("ggH", "ttH", "bbH", "HJJ") else 1
            if self.hypothesis in ("0+", "0-", "L1", "L1Zg"): return 0

            g2 = {"dec": constants.g2HZZ}
            if self.productionmode == "VBF": g2["prod"] = constants.g2VBF
            elif self.productionmode == "ZH":  g2["prod"] = constants.g2ZH
            elif self.productionmode == "WH":  g2["prod"] = constants.g2WH
            else: g2["prod"] = None
            g2["proddec"] = None if g2["prod"] is None else sqrt(g2["prod"]*g2["dec"])

            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        if a == "fa2":
                            return g2[proddec]
                        return 0
                    if self.hypothesis == a+proddec+"-0.5":
                        if a == "fa2":
                            return -g2[proddec]
                        return 0

                minus = "-" if proddec == "proddec" else ""
                for a, b in combinations(couplings, 2):
                    if self.hypothesis in (a+proddec+"0.5"+b+proddec+minus+"0.5",
                                           a+proddec+"0.33"+b+proddec+minus+"0.33"):
                        if a == "fa2":
                            return g2[proddec]
                        if b == "fa2":
                            return float(minus+"1") * g2[proddec]
                        return 0
                for a, b, c in combinations(couplings, 3):
                    if (self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33"
                          or proddec=="proddec"
                          and self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"
                    ):
                        if a == "fa2" or b == "fa2" or "0.25" in str(self.hypothesis) and c == "fa2":
                            return g2[proddec]
                        if c == "fa2":
                            return float(minus+"1") * g2[proddec]
                        return 0

                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        if a == "fa2" or b == "fa2" or c == "fa2" or d == "fa2":
                            return g2[proddec]
                        return 0

            if self.hypothesis == "fa2dec-0.9": return -constants.g2HZZ * 3

        raise self.ValueError("g2")

    @property
    def g4(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "HJJ"):
            if self.hypothesis == "0-": return constants.g4HZZ if self.productionmode in ("ggH", "ttH", "bbH", "HJJ") else 1
            if self.hypothesis in ("0+", "a2", "L1", "L1Zg"): return 0

            g4 = {"dec": constants.g4HZZ}
            if self.productionmode == "VBF": g4["prod"] = constants.g4VBF
            elif self.productionmode == "ZH":  g4["prod"] = constants.g4ZH
            elif self.productionmode == "WH":  g4["prod"] = constants.g4WH
            else: g4["prod"] = None
            g4["proddec"] = None if g4["prod"] is None else sqrt(g4["prod"]*g4["dec"])

            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        if a == "fa3":
                            return g4[proddec]
                        return 0
                    if self.hypothesis == a+proddec+"-0.5":
                        if a == "fa3":
                            return -g4[proddec]
                        return 0

                minus = "-" if proddec == "proddec" else ""
                for a, b in combinations(couplings, 2):
                    if self.hypothesis in (a+proddec+"0.5"+b+proddec+minus+"0.5",
                                           a+proddec+"0.33"+b+proddec+minus+"0.33"):
                        if a == "fa3":
                            return g4[proddec]
                        if b == "fa3":
                            return float(minus+"1") * g4[proddec]
                        return 0
                for a, b, c in combinations(couplings, 3):
                    if (self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33"
                          or proddec=="proddec"
                          and self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"
                    ):
                        if a == "fa3" or b == "fa3" or "0.25" in str(self.hypothesis) and c == "fa3":
                            return g4[proddec]
                        if c == "fa3":
                            return float(minus+"1") * g4[proddec]
                        return 0

                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        if a == "fa3" or b == "fa3" or c == "fa3" or d == "fa3":
                            return g4[proddec]
                        return 0

            if self.hypothesis == "fa2dec-0.9": return 0

        raise self.ValueError("g4")

    @property
    def g1prime2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "HJJ"):
            if self.hypothesis == "L1": return constants.g1prime2HZZ if self.productionmode in ("ggH", "ttH", "bbH", "HJJ") else 1
            if self.hypothesis in ("0+", "a2", "0-", "L1Zg"): return 0

            g1prime2 = {"dec": constants.g1prime2HZZ}
            if self.productionmode == "VBF": g1prime2["prod"] = constants.g1prime2VBF
            elif self.productionmode == "ZH":  g1prime2["prod"] = constants.g1prime2ZH
            elif self.productionmode == "WH":  g1prime2["prod"] = constants.g1prime2WH
            else: g1prime2["prod"] = None
            g1prime2["proddec"] = None if g1prime2["prod"] is None else -sqrt(g1prime2["prod"]*g1prime2["dec"])

            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        if a == "fL1":
                            return g1prime2[proddec]
                        return 0
                    if self.hypothesis == a+proddec+"-0.5":
                        if a == "fL1":
                            return -g1prime2[proddec]
                        return 0

                minus = "-" if proddec == "proddec" else ""
                for a, b in combinations(couplings, 2):
                    if self.hypothesis in (a+proddec+"0.5"+b+proddec+minus+"0.5",
                                           a+proddec+"0.33"+b+proddec+minus+"0.33"):
                        if a == "fL1":
                            return g1prime2[proddec]
                        if b == "fL1":
                            return float(minus+"1") * g1prime2[proddec]
                        return 0
                for a, b, c in combinations(couplings, 3):
                    if (self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33"
                          or proddec=="proddec"
                          and self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"
                    ):
                        if a == "fL1" or b == "fL1" or "0.25" in str(self.hypothesis) and c == "fL1":
                            return g1prime2[proddec]
                        if c == "fL1":
                            return float(minus+"1") * g1prime2[proddec]
                        return 0

                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        if a == "fL1" or b == "fL1" or c == "fL1" or d == "fL1":
                            return g1prime2[proddec]
                        return 0

            if self.hypothesis == "fa2dec-0.9": return 0

        raise self.ValueError("g1prime2")

    @property
    def ghzgs1prime2(self):
        if self.productionmode in ("ggH", "VBF", "ZH", "WH", "WplusH", "WminusH", "ttH", "bbH", "HJJ"):
            if self.hypothesis == "L1Zg": return constants.ghzgs1prime2HZZ if self.productionmode in ("ggH", "ttH", "bbH", "HJJ") else 1
            if self.hypothesis in ("0+", "a2", "0-", "L1"): return 0

            ghzgs1prime2 = {"dec": constants.ghzgs1prime2HZZ}
            if self.productionmode == "VBF": ghzgs1prime2["prod"] = constants.ghzgs1prime2VBF
            elif self.productionmode == "ZH":  ghzgs1prime2["prod"] = constants.ghzgs1prime2ZH
            elif self.productionmode == "WH":  ghzgs1prime2["prod"] = constants.ghzgs1prime2WH
            else: ghzgs1prime2["prod"] = None
            ghzgs1prime2["proddec"] = None if ghzgs1prime2["prod"] is None else -sqrt(ghzgs1prime2["prod"]*ghzgs1prime2["dec"])

            couplings = "fa3", "fa2", "fL1", "fL1Zg"
            for proddec in "prod", "dec", "proddec":
                for a in couplings:
                    if self.hypothesis == a+proddec+"0.5":
                        if a == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        return 0
                    if self.hypothesis == a+proddec+"-0.5":
                        if a == "fL1Zg":
                            return -ghzgs1prime2[proddec]
                        return 0

                minus = "-" if proddec == "proddec" else ""
                for a, b in combinations(couplings, 2):
                    if self.hypothesis in (a+proddec+"0.5"+b+proddec+minus+"0.5",
                                           a+proddec+"0.33"+b+proddec+minus+"0.33"):
                        if a == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        if b == "fL1Zg":
                            return float(minus+"1") * ghzgs1prime2[proddec]
                        return 0
                for a, b, c in combinations(couplings, 3):
                    if (self.hypothesis == a+proddec+"0.33"+b+proddec+"0.33"+c+proddec+minus+"0.33"
                          or proddec=="proddec"
                          and self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"
                    ):
                        if a == "fL1Zg" or b == "fL1Zg" or "0.25" in str(self.hypothesis) and c == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        if c == "fL1Zg":
                            return float(minus+"1") * ghzgs1prime2[proddec]
                        return 0

                if proddec != "proddec": continue
                for a, b, c, d in combinations(couplings, 4):
                    if self.hypothesis == a+proddec+"0.25"+b+proddec+"0.25"+c+proddec+"0.25"+d+proddec+"0.25":
                        if a == "fL1Zg" or b == "fL1Zg" or c == "fL1Zg" or d == "fL1Zg":
                            return ghzgs1prime2[proddec]
                        return 0

            if self.hypothesis == "fa2dec-0.9": return 0

        raise self.ValueError("ghzgs1prime2")

    @property
    def ghg2(self):
        if self.productionmode == "ggH":
            return 1
        if self.productionmode == "HJJ":
            if self.hffhypothesis in ("Hff0+", "fCP0.5"):
                return 1
            if self.hffhypothesis == "Hff0-":
                return 0
        raise self.ValueError("ghg2")

    @property
    def ghg4(self):
        if self.productionmode == "ggH":
            return 0
        if self.productionmode == "HJJ":
            if self.hffhypothesis == "Hff0+":
                return 0
            if self.hffhypothesis == "Hff0-":
                return 1
            if self.hffhypothesis == "fCP0.5":
                return constants.ghg4HJJ
        raise self.ValueError("ghg4")

    @property
    def kappa(self):
        if self.productionmode == "ttH":
            if self.hffhypothesis in ("Hff0+", "fCP0.5"):
                return 1
            if self.hffhypothesis == "Hff0-":
                return 0
        if self.productionmode == "tqH":
            assert self.hffhypothesis == "Hff0+"
            return 1
        raise self.ValueError("kappa")

    @property
    def kappa_tilde(self):
        if self.productionmode == "ttH":
            if self.hffhypothesis == "Hff0+":
                return 0
            if self.hffhypothesis == "Hff0-":
                return 1
            if self.hffhypothesis == "fCP0.5":
                return constants.kappa_tilde_ttH
        if self.productionmode == "tqH":
            assert self.hffhypothesis == "Hff0+"
            return 0
    @property
    def pdf(self): return None

class ReweightingSampleWithPdf(ReweightingSample):
    enums = [ReweightingSample, Production]
    @property
    def pdf(self): return self.production.pdf

class ReweightingSampleWithFlavor(ReweightingSample):
    enums = [ReweightingSample, Flavor]
    def check(self, *args):
        if self.productionmode in ("ggZZ", "VBF bkg"):
            if self.flavor is None:
                raise ValueError("No flavor provided for {} productionmode\n{}".format(self.productionmode, args))
            if self.productionmode == "VBF bkg" and self.flavor.hastaus:
                raise ValueError("No {} samples with taus\n{}".format(self.productionmode, args))
        else:
            if self.flavor is not None:
                raise ValueError("Flavor provided for {} productionmode\n{}".format(self.productionmode, args))

class ReweightingSamplePlus(ReweightingSample):
    enums = [ReweightingSample, AlternateGenerator, PythiaSystematic, Extension]
    enumname = "reweightingsampleplus"

    def check(self, *args):
        if (
            self.alternategenerator == "POWHEG"
            and (
                 self.hypothesis != "0+"
                 or self.productionmode not in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH")
                )
           ) or (
            self.alternategenerator in ("MINLO", "NNLOPS")
            and (
                 self.hypothesis != "0+"
                 or self.productionmode != "ggH"
                )
           ):
            raise ValueError("No {} {} sample produced with {}\n{}".format(self.productionmode, self.hypothesis, self.alternategenerator, args))

        if (
            self.alternategenerator == "POWHEG"
            and (
                 self.hffhypothesis != "Hff0+"
                 and self.hffhypothesis is not None
                )
           ):
            raise ValueError("No {} {} sample produced with {}\n{}".format(self.productionmode, self.hffhypothesis, self.alternategenerator, args))

        if (
            self.extension == "ext"
            and not (
                     self.productionmode == "qqZZ"
                     or self.productionmode in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH") and self.alternategenerator == "POWHEG"
                    )
           ):
            raise ValueError("No extra {} sample produced\n{}".format(self.productionmode, args))

        if self.pythiasystematic is not None:
            if (
                self.hypothesis != "0+"
                or self.productionmode not in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH")
               ):
                raise ValueError("No {} {} sample produced with pythia {}\n{}".format(self.productionmode, self.hypothesis, self.pythiasystematic, args))
            if self.alternategenerator not in ("POWHEG", "MINLO"):
                raise ValueError("{} sample with pythia {} is produced with POWHEG{}\n{}".format(self.productionmode, self.pythiasystematic, " or MINLO" if self.productionmode == "ggH" else "", args))

        super(ReweightingSamplePlus, self).check(*args)

    def weightname(self):
        result = super(ReweightingSamplePlus, self).weightname()
        if self.alternategenerator:
            result += "_" + str(self.alternategenerator)
        if self.pythiasystematic:
            result += "_" + str(self.pythiasystematic)
        return result

class Sample(ReweightingSamplePlus):
    enums = [ReweightingSamplePlus, Flavor, Production]

    def check(self, *args):
        if self.production is None:
            raise ValueError("No option provided for production\n{}".format(args))

        if self.hypothesis is not None and self.hypothesis not in self.productionmode.generatedhypotheses:
            raise ValueError("No {} sample produced with hypothesis {}!\n{}".format(self.productionmode, self.hypothesis, args))

        if self.productionmode in ("WplusH", "WminusH") and self.alternategenerator != "POWHEG":
            raise ValueError("Separate {} sample is produced with POWHEG.  Maybe you meant to specify POWHEG, or WH?\n{}".format(self.productionmode, args))

        self.reweightingsamplewithflavor = ReweightingSampleWithFlavor(self.reweightingsample, self.flavor)
        self.reweightingsamplewithpdf = ReweightingSampleWithPdf(self.reweightingsample, self.production)

        if self.pythiasystematic is not None and self.alternategenerator == "NNLOPS":
            raise ValueError("No NNLOPS samples with systematics!\n{}".format(args))

        super(Sample, self).check(*args)

    def CJLSTmaindir(self):
      if self.alternategenerator is None:
        if self.productionmode == "ggH":
          return self.production.CJLSTdir_anomalous()
        if self.productionmode == "VBF":
          return self.production.CJLSTdir_anomalous_VBF()
        if self.productionmode in ("ZH", "WH"):
          return self.production.CJLSTdir_anomalous_VH()
        if self.productionmode in ("data", "ZX"):
          return self.production.CJLSTdir_data()
      if self.productionmode == "ggH" and self.alternategenerator == "MINLO" and self.pythiasystematic is None:
        return self.production.CJLSTdir_MINLO()
      return self.production.CJLSTdir()

    def CJLSTdirname(self):
        if self.alternategenerator is not None:
            if self.alternategenerator == "POWHEG":
                if self.productionmode in ("ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH"):
                    s = str(self.productionmode)
                    if self.productionmode == "VBF": s = "VBFH"
                    if self.hypothesis == "0+": result = "{}125".format(s)
                if self.pythiasystematic is not None:
                    result += self.pythiasystematic.appendname
                if self.extension == "ext": result += "ext"
                return result
        if self.alternategenerator in ("MINLO", "NNLOPS"):
            if self.productionmode == "ggH":
                result = "ggH125"
                if self.pythiasystematic is not None:
                    result += self.pythiasystematic.appendname
                if self.alternategenerator == "MINLO":
                    result += "_minloHJJ"
                if self.alternategenerator == "NNLOPS":
                    result += "_NNLOPS"
                return result
            raise self.ValueError("CJLSTdirname")
        if self.productionmode in ("ggH", "VBF", "ZH", "WH"):
            s = str(self.productionmode)
            if self.productionmode == "VBF": s = "VBFH"
            if self.hypothesis == "0+": result = "{}0PM_M125".format(s)
            if self.hypothesis == "a2": result = "{}0PH_M125".format(s)
            if self.hypothesis == "0-": result = "{}0M_M125".format(s)
            if self.hypothesis == "L1": result = "{}0L1_M125".format(s)
            if self.hypothesis in ("fa20.5", "fa2prod0.5"): result = "{}0PHf05ph0_M125".format(s)
            if self.hypothesis in ("fa30.5", "fa3prod0.5"): result = "{}0Mf05ph0_M125".format(s)
            if self.hypothesis in ("fL10.5", "fL1prod0.5"): result = "{}0L1f05ph0_M125".format(s)
            return result
        if self.productionmode in ("HJJ", "ttH"):
            s = str(self.productionmode)
            if self.hffhypothesis == "Hff0+": return "{}0PM_M125".format(s)
            if self.hffhypothesis == "Hff0-": return "{}0M_M125".format(s)
            if self.hffhypothesis == "fCP0.5": return "{}0Mf05ph0_M125".format(s)
        if self.productionmode == "bbH": return "bbH125"
        if self.productionmode == "tqH": return "tqH125"
        if self.productionmode == "ggZZ":
            return "ggTo{}_Contin_MCFM701".format(self.flavor)
        if self.productionmode == "VBF bkg":
            return "VBFTo{}JJ_Contin_phantom128".format(self.flavor)
        if self.productionmode == "qqZZ":
            if self.extension == "ext":
                return "ZZTo4lext"
            else:
                return "ZZTo4l"
        if self.productionmode == "ZX" or self.productionmode == "data":
            return "AllData"
        raise self.ValueError("CJLSTdirname")

    def CJLSTfile(self):
        if config.LHE: raise ValueError("Can't get the CJLST file when in LHE mode!")
        return os.path.join(self.CJLSTmaindir(), self.CJLSTdirname(), "ZZ4lAnalysis.root")

    @property
    def LHEfile(self):
        if not config.LHE: raise ValueError("Can't get the lhe file when not in LHE mode!")
        if self.production == "LHE_170509":
            if self.productionmode == "ggH":
                folder = "/work-zfs/lhc/heshy/LHEanomalouscouplings/processlhe/fL1fL1Zg"
                if self.hypothesis == "0+": filename = "ggHa1.lhe"
                if self.hypothesis == "L1": filename = "ggHL1.lhe"
                if self.hypothesis == "fL10.5": filename = "ggHa1L1.lhe"
                if self.hypothesis == "L1Zg": filename = "ggHL1Zg.lhe"
                if self.hypothesis == "fL1Zg0.5": filename = "ggHa1L1Zg.lhe"
                if self.hypothesis == "fL10.5fL1Zg0.5": filename = "ggHL1L1Zg.lhe"
                return os.path.join(folder, filename)
            if self.productionmode == "qqZZ":
                return "/work-zfs/lhc/ianderso/hep/LHEFiles/qqZZ/MG/8T/ZZJetsTo4L_TuneZ2_8TeV-madgraph-tauola.lhe"
            if self.productionmode == "data":
                tmpdir = mkdtemp()
                filename = os.path.join(tmpdir, "empty.lhe")
                with open(filename, 'w') as f:
                    pass
                return filename
        raise self.ValueError("LHEfile")

    def withdiscriminantsfile(self):
        if self.copyfromothersample: return self.copyfromothersample.withdiscriminantsfile()
        result = os.path.join(config.repositorydir, "step3_withdiscriminants", "{}.root".format(self).replace(" ", ""))
        return result
        raise self.ValueError("withdiscriminantsfile")

    @staticmethod
    @cache
    def effectiveentries(reweightfrom, reweightto):
        from utilities import tfiles
        if (reweightto.productionmode in ("ggZZ", "VBF bkg", "ZX", "data")
                or reweightfrom.alternategenerator is not None
                or config.LHE):
            assert reweightfrom.reweightingsample == reweightto
            return 1
        f = tfiles[reweightfrom.withdiscriminantsfile()]
        t = f.effectiveentries
        assert t.GetEntries() == 1, "{}???".format(t.GetEntries())
        t.GetEntry(0)
        try:
            return getattr(t, reweightto.weightname())
        except:
            t.Show()
            raise

    @property
    def xsec(self):
        try:
            return super(Sample, self).xsec
        except ValueError:
            t = ROOT.TChain("candTree")
            t.Add(self.withdiscriminantsfile())
            t.GetEntry(1)
            result = t.xsec
            if t.xsec <= 0:
                raise ValueError("Can't get xsec for {!r}".format(self))
            return result

    @property
    def copyfromothersample(self):
        if not config.LHE and self == Sample("170203", "ggH", "MINLO", "0+"):
            return Sample("170222", "ggH", "MINLO", "0+")
        if self.production == "170712" and self.productionmode in ("ggH", "HJJ", "VBF", "WH", "ZH", "ttH") and self.alternategenerator is None:
            args = [self.reweightingsampleplus, self.flavor, "170222"]
            return Sample(*args)
        if self.production in ("180224_10bins", "180224_newdiscriminants") or "_Ulascan" in str(self):
            production = str(self.production).split("_")[0]
            args = [self.flavor, self.reweightingsampleplus, production]
            return Sample(*args)
        return None

    @property
    def pdf(self): return self.reweightingsamplewithpdf.pdf

    @cache
    def alternateweightxsec(self, alternateweight):
        if alternateweight in ("1", "PythiaScaleUp", "PythiaScaleDown"):
            return 1
        with TFile(self.withdiscriminantsfile()) as f:
             t = f.alternateweightxsecstree
             t.GetEntry(0)
             return getattr(t, alternateweight.weightname)

class SampleBasis(MultiEnum):
    enums = [ProductionMode, Analysis]
    def __init__(self, hypotheses, *args):
        self.hypotheses = [Hypothesis(_) for _ in hypotheses]
        super(SampleBasis, self).__init__(*args)

    def __hash__(self):
        return hash((super(SampleBasis, self).__hash__(), tuple(self.hypotheses)))
    def __eq__(self, other):
        return super(SampleBasis, self).__eq__(other) and self.hypotheses == other.hypotheses
    #__ne__ is already defined as not == in MultiEnum

    def check(self, *args):
        args = (self.hypotheses,)+args
        if self.productionmode in ("ggH", "ttH", "bbH"):
            if self.analysis.dimensions == 4:
                dimension = 15
            elif self.analysis.dimensions == 2:
                dimension = 6
            else:
                dimension = 3
        elif self.productionmode in ("VBF", "WH", "ZH"):
            if self.analysis.dimensions == 4:
                dimension = 70
            elif self.analysis.dimensions == 2:
                assert False  #have to figure this out
            else:
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
        hffhypothesis = None
        if self.productionmode == "ttH": hffhypothesis = "Hff0+"
        dimension = len(self.hypotheses)
        samples = [ReweightingSample(self.productionmode, _, hffhypothesis) for _ in self.hypotheses]
        if self.analysis.dimensions == 4:
            if dimension == 15:
                maxpower = 2
            elif dimension == 70:
                maxpower = 4
            else:
                assert False, dimension
            assert len(self.analysis.couplingnames) == 4
            return numpy.matrix(
                                [
                                 [
                                  sample.g1**(maxpower-i-j-k-l)
                                  * getattr(sample, self.analysis.couplingnames[0])**i
                                  * getattr(sample, self.analysis.couplingnames[1])**j
                                  * getattr(sample, self.analysis.couplingnames[2])**k
                                  * getattr(sample, self.analysis.couplingnames[3])**l
                                     for l in range(maxpower+1)
                                     for k in range(maxpower+1-l)
                                     for j in range(maxpower+1-k-l)
                                     for i in range(maxpower+1-j-k-l)
                                 ]
                                    for sample in samples
                                ]
                               )

        elif self.analysis.dimensions == 2:
            assert dimension == 6  #not VBF or VH
            maxpower = 2
            return numpy.matrix(
                                [
                                 [
                                  sample.g1**(maxpower-i-j)
                                  * getattr(sample, self.analysis.couplingnames[0])**i
                                  * getattr(sample, self.analysis.couplingnames[1])**j
                                     for j in range(maxpower+1)
                                     for i in range(maxpower+1-j)
                                 ]
                                    for sample in samples
                                ]
                               )
        else:
            maxpower = dimension-1
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

def allsamples():
    __xcheck()
    for production in productions:
        if "Ulascan" in str(production): continue
        for productionmode in "ggH", "VBF", "ZH", "WH", "bbH":
            for hypothesis in ProductionMode(productionmode).generatedhypotheses:
                yield Sample(productionmode, hypothesis, production)
        for hypothesis in hffhypotheses:
            yield Sample("HJJ", hypothesis, "0+", production)
            yield Sample("ttH", hypothesis, "0+", production)
        yield Sample("tqH", "0+", "Hff0+", production)
        for flavor in flavors:
            yield Sample("ggZZ", flavor, production)
            if not flavor.hastaus:
                yield Sample("VBF bkg", flavor, production)
        yield Sample("ggH", "0+", "NNLOPS", production)
        yield Sample("ggH", "0+", "MINLO", production)
        for productionmode in "ggH", "VBF", "ZH", "WplusH", "WminusH":
            yield Sample(productionmode, "0+", "POWHEG", production)
            if production.year == 2017: yield Sample(productionmode, "0+", "POWHEG", "ext", production)
        yield Sample("ttH", "Hff0+", "0+", "POWHEG", production)
        if production.year == 2017: yield Sample("ttH", "Hff0+", "0+", "POWHEG", "ext", production)
        for systematic in pythiasystematics:
            if systematic in ("ScaleUp", "ScaleDown") and production.year >= 2017: continue
            for productionmode in "ggH", "VBF", "ZH", "WplusH", "WminusH":
                yield Sample(productionmode, "0+", "POWHEG", production, systematic)
            yield Sample("ttH", "Hff0+", "0+", "POWHEG", production, systematic)
            yield Sample("ggH", "0+", "MINLO", production, systematic)
        yield Sample("qqZZ", production)
        yield Sample("qqZZ", "ext", production)
        yield Sample("ZX", production)
        yield Sample("data", production)

if config.LHE:
    def allsamples():
        __xcheck()
        for production in productions:
            for hypothesis in "0+", "L1Zg", "fL1Zg0.5", "fL10.5fL1Zg0.5":
                yield Sample("ggH", hypothesis, production)
            yield Sample("qqZZ", production)
            yield Sample("data", production)

def __xcheck():
    global __xcheck
    def __xcheck(): pass

    class WriteOnceDict(dict):
        def __init__(self, messagefmt="{key} has already been set"):
            self.messagefmt = messagefmt
        def __setitem__(self, key, value):
            if key in self:
                if value == self[key]: return
                raise KeyError(self.messagefmt.format(key=key, newvalue=value, oldvalue=self[key]))
            super(WriteOnceDict, self).__setitem__(key, value)

    if config.LHE:
        def CJLSTfile(sample):
             return sample.LHEfile
    else:
        def CJLSTfile(sample):
             return sample.CJLSTfile()

    CJ2wd = WriteOnceDict("Multiple withdiscriminantsfiles:\n  {oldvalue}\n  {newvalue}\nfrom the same CJLST file:\n  {key}")
    wd2CJ = WriteOnceDict("Multiple CJLST files:\n  {oldvalue}\n  {newvalue}\nthat give the same withdiscriminantsfile:\n  {key}")

    for sample in allsamples():
        if sample.productionmode == "data": continue
        try:
            CJ2wd[CJLSTfile(sample)] = sample.withdiscriminantsfile()
            wd2CJ[sample.withdiscriminantsfile()] = CJLSTfile(sample)
        except KeyError as e:
            print e
            raise

if __name__ == "__main__":
    for s in (
      ReweightingSample("ggH", "fa2dec0.5fa3dec0.5"),
      ReweightingSample("VBF", "fa2prod0.5fa3prod0.5"),
      ReweightingSample("ZH",  "fa2proddec-0.33fa3proddec0.33"),
      ReweightingSample("VBF", "fa2proddec0.25fa3proddec0.25fL1proddec0.25"),
    ):
        print s, s.g1, s.g2, s.g4, s.g1prime2, s.ghzgs1prime2
