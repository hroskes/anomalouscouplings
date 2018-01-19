from abc import ABCMeta, abstractmethod, abstractproperty
import inspect
import re

class MakeSystematics(object):
    __metaclass__ = ABCMeta
    def __init__(self, function):
        while isinstance(function, MakeSystematics):
            function = function.function
        self.function = function

    @property
    def name(self): return self.function.__name__
    @abstractproperty
    def upname(self): pass
    @abstractproperty
    def dnname(self): pass

    def getnominal(self): return self.function

    def getupdn(self):
        #python is magic!
        sourcelines = inspect.getsourcelines(self.function)[0]
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

        return Up, Dn

    @abstractmethod
    def doreplacements(self, code):
        return code

class MakeJECSystematics(MakeSystematics):
    @property
    def upname(self): return self.name+"_JECUp"
    @property
    def dnname(self): return self.name+"_JECDn"

    def doreplacements(self, code):
        for thing in (
            r"(self[.]M2(?:g1)?(?:g2|g4|g1prime2|ghzgs1prime2)?_(?:VBF|HadZH|HadWH))",
            r"(self[.]notdijet)",
            r"(self[.]binning_4couplings_(HadVH|VBF|)decay)",
        ):
            code = re.sub(thing, r"\1_JEC{UpDn}", code)
        for thing in (
            r"(self[.](Hjj|Jet1)Pt)",
            r"(self[.]DiJet(Mass|DEta))",
            r"(self[.]nCleanedJetsPt30)",
        ):
            code = re.sub(thing, r"\1_jec{UpDn}", code)

        result = re.findall("(self[.][\w]*)[^\w{]", code)
        for variable in result:
            if re.match("self[.](M2(?:g1)?(?:g2|g4|g1prime2|ghzgs1prime2)?_decay|ZZ(Eta|Pt)|D_4couplings_general)", variable): continue
            raise ValueError("Unknown self.variable in function '{}':\n\n{}\n{}".format(self.name, code, variable))

        return code
