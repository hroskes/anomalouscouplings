from abc import ABCMeta, abstractmethod, abstractproperty
import inspect
import re

class MakeSystematics(object):
    __metaclass__ = ABCMeta
    def __init__(self, function):
        if isinstance(function, MakeSystematics):
            self.otheralternates = function.getalternates()
            self.nominal = function.getnominal()
        else:
            self.otheralternates = []
            self.nominal = function

    @property
    def name(self): return self.nominal.__name__
    @abstractproperty
    def upname(self): pass
    @abstractproperty
    def dnname(self): pass

    def getnominal(self): return self.nominal

    def getalternates(self):
        #python is magic!
        sourcelines = inspect.getsourcelines(self.nominal)[0]
        for i, line in enumerate(sourcelines):
            if "@" not in line:
                break
        del sourcelines[0:i]
        lookfor = "def "+self.name+"("
        replacewith = "def {UpDn}("
        if lookfor not in sourcelines[0]:
            raise ValueError("function '{}' is defined with weird syntax:\n\n{}\n\nWant to have '{}' in the first line".format(self.name, "".join(sourcelines), lookfor))
        sourcelines[0] = sourcelines[0].replace(lookfor, replacewith)
        code = "".join(sourcelines) #they already have '\n' in them

        code = self.doreplacements(code)

        code = re.sub("^ *(def )", r"\1", code) #so that we don't get IndentationError in exec

        exec code.format(UpDn="Up")
        exec code.format(UpDn="Dn")

        Up.__name__ = self.upname
        Dn.__name__ = self.dnname

        return self.otheralternates + [Up, Dn]

    @abstractmethod
    def doreplacements(self, code):
        return code

class MakeJetSystematicsBase(MakeSystematics):
    @abstractproperty
    def JEX(self): pass
    @property
    def JEXlower(self): return self.JEX.lower()
    @property
    def JEXlowerforMET(self): return self.JEXlower.replace("jec", "jes")
    @property
    def upname(self): return self.name+"_"+self.JEX+"Up"
    @property
    def dnname(self): return self.name+"_"+self.JEX+"Dn"

    def doreplacements(self, code):
        for thing in (
            r"(self[.][\w]*)_JECNominal\b",
            r"(self[.]M2(?:(?:g1)?(?:g(?:|hzgs|hgsgs)(?:2|4|1prime2))?(?:Zg|gg)?_(?:VBF|HadZH|HadWH))|qqZZJJ)\b",
            r"(self[.]notdijet)\b",
            r"(self[.](?:binning_4couplings(?:_photons)?|D_bkg_kin)_(?:HadVH|VBF|)(?:decay)?)\b",
        ):
            code = re.sub(thing, r"\1_"+self.JEX+"{UpDn}", code)
        for thing in (
            r"(self[.](Hjj|Jet1)Pt)\b",
            r"(self[.]DiJet(Mass|DEta))\b",
            r"(self[.]nCleanedJetsPt30(BTagged_bTagSF)?)\b",
            r"(self[.]jet(QGLikelihood|Phi))\b",
        ):
            code = re.sub(thing, r"\1_"+self.JEXlower+"{UpDn}", code)
        for thing in (
            r"(self[.]PFMET_corrected)\b",
        ):
            code = re.sub(thing, r"\1_"+self.JEXlowerforMET+"{UpDn}", code)

        result = re.findall("(self[.][\w]*)[^\w{]", code)
        for variable in result:
            if re.match(
              "self[.]("
                "M2(?:g1)?(?:g(?:|hzgs|hgsgs)(?:2|4|1prime2))?(?:Zg|gg)?_decay"
                "|ZZ(?:Eta|Pt|Mass)|D_4couplings_general(?:_raw|)|foldbins_4couplings_(?:VBF|HadVH)decay"
                "|p_m4l_(?:SIG|BKG)|cconstantforDbkg(?:kin)?|flavor"
                "|g(?:(?:|hzgs|hgsgs)(?:2|4|1prime2))(?:Zg|gg)?(?:HZZ|VBF|VH|ZH|WH|HZg|Hgg)_m4l"
                "|nExtra(?:Lep|Z)"
              ")$",
              variable
            ): continue
            raise ValueError("Unknown self.variable for '{}' in function '{}':\n\n{}\n{}".format(type(self).__name__, self.name, code, variable))

        return code

class MakeJECSystematics(MakeJetSystematicsBase):
  @property
  def JEX(self): return "JEC"
class MakeJESSystematics(MakeJetSystematicsBase):
  @property
  def JEX(self): return "JES"
class MakeJERSystematics(MakeJetSystematicsBase):
  @property
  def JEX(self): return "JER"

def MakeJetSystematics(function):
  return MakeJECSystematics(MakeJESSystematics(MakeJERSystematics(function)))

class MakeBtagSystematics(MakeSystematics):
    @property
    def upname(self): return self.name+"_bTagSFUp"
    @property
    def dnname(self): return self.name+"_bTagSFDn"

    def doreplacements(self, code):
        for thing in (
            r"(self[.][\w]*)_bTagSF\b",
        ):
            code = re.sub(thing, r"\1_bTagSF{UpDn}", code)

        result = re.findall("(self[.][\w]*)[^\w{]", code)
        for variable in result:
            if re.match(
              "self[.]("
              "nCleanedJetsPt30|nExtra(Lep|Z)|jet(QGLikelihood|Phi)"
              "|p(Aux)?_(J+(QCD|VBF)|Had[ZW]H)_(SIG_gh[vzwg][12]_1_JHUGen|mavjj(_true)?)_JECNominal"
              "|PFMET_corrected|ZZ(Mass|Pt)|DiJetMass|HjjPt"
              ")$",
              variable
            ): continue
            raise ValueError("Unknown self.variable for '{}' in function '{}':\n\n{}\n{}".format(type(self).__name__, self.name, code, variable))

        return code


#The following imports are not needed for any of the code written explicitly here,
#but when the functions are defined with exec, it looks for the global variables
#in this module, rather than where the original function was defined.

from math import sqrt

import CJLSTscripts
import config
import constants
import STXS
