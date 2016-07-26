import config
from enums import Channel, channels, EnumItem, MultiEnum, MyEnum, Production
import os
import yaml

class SystematicsType(MyEnum):
    enumitems = (
                 EnumItem("Moriond"),
                 EnumItem("ICHEP"),
                )
    enumname = "systematicstype"

class ReadSystematics(MultiEnum):
    enums = (Production, Channel, SystematicsType)
    def __init__(self, *args):
        super(ReadSystematics, self).__init__(*args)
        self.systematics = self.getsystematics()

    def check(self, *args, **kwargs):
        assert self.systematicstype is None
        if self.production.year == 2015:
            self.systematicstype = "Moriond"
        elif self.production.year == 2016:
            self.systematicstype = "ICHEP"
        else:
            assert False
        super(ReadSystematics, self).check(*args, **kwargs)

    def getsystematics(self):
        if self.systematicstype == "Moriond":
            with open(self.channel.moriondcardfile()) as f:
                self.contents = f.read()

            sections = self.contents.split("------------------------------------------------------------")
            systematicssection = sections[-1]
            return {line.split()[0]: line for line in systematicssection.split("\n") if line.strip()}
        elif self.systematicstype == "ICHEP":
            filenames = [os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_ICHEP2016", "LegoCards", "configs", "inputs", f)
                           for f in [
                                     "systematics_13TeV_{}.yaml".format(self.channel),
                                     "systematics_expt_13TeV.yaml",
                                     "systematics_leptonScaleAndResol_13TeV.yaml",
                                     "systematics_theory_13TeV.yaml",
                                    ]
                        ]
            self.systematics = {}
            for filename in filenames:
                with open(filename) as f:
                    y = yaml.load(f)
                assert not set(y).intersection(set(self.systematics))
                self.systematics.update(y)
            return self.systematics

    def __getitem__(self, item):
        return self.systematics[item]

    @classmethod
    def singlesystematictype(cls, production, *args, **kwargs):
        self = cls(production, "2e2mu")
        if self.systematicstype == "Moriond":
            return MoriondSystematic(*args, **kwargs)
        if self.systematicstype == "ICHEP":
            return ICHEPSystematic(*args, **kwargs)
        assert False

class Systematic(object):
    def __init__(self, name, useprocesses=None):
        self.name = name
        if useprocesses is None:
            useprocesses = ["ggH", "qqZZ", "ggZZ", "zjets"]
        self.useprocesses = useprocesses

class MoriondSystematic(Systematic):
    def __init__(self, name, useprocesses=None):
        super(MoriondSystematic, self).__init__(name, useprocesses)
        self.useindices = [self.getindex(processname) for processname in self.useprocesses]

    @staticmethod
    def getindex(processname):
        return "WH ZH ggH qqH ttH ggZZ qqZZ zjets".split().index(processname)

    def getline(self, name, readsystematics):
        line = readsystematics[self.name]
        split = line.split()
        if split[0] != self.name:
            raise ValueError("Wrong systematic passed to {}({}).getline():\n{}".format(type(self).__name__, self.name, line))
        if split[1] != "lnN":
            raise ValueError("Non-normalization systematic passed to {}.getline():\n{}".format(type(self).__name__, self.name, line))
        uselist = [name, split[1]] + [split[i+2] for i in self.useindices]
        return " ".join(uselist)

class ICHEPSystematic(Systematic):
    def getline(self, name, readsystematics):
        if self.name == "QCDscale_ggVV":   #included in ggH
            uselist = [name, "lnN"] + ["-" for p in self.useprocesses]
            return " ".join(uselist)
        try:
            s = readsystematics[self.name]["UnTagged"]
        except KeyError:
            s = readsystematics[self.name]["Any"]
        if s["type"] != "lnN":
            raise ValueError("Non-normalization systematic passed to {}.getline():\n{}".format(type(self).__name__, self.name, line))
        uselist = [name, s["type"]] + [str(s[p]) if p in s else "-" for p in self.useprocesses]
        return " ".join(uselist)

class Systematic(object):
    def __init__(self, line, production):
        self.name = self.type = None
        try:
            self.name = line.split()[0]
            self.type = line.split()[1]
        except IndexError:
            pass
        self.line = line
        self.production = Production(production)

    def newname(self):
        if self.type != "lnN": return None
        if self.name == "lumi_8TeV": return "lumi_13TeV"
        if self.name == "pdf_gg": return "pdf_Higgs_gg"
        if self.name == "pdf_qqbar": return "pdf_qq"
        if self.name == "pdf_hzz4l_accept": return None
        if self.name == "QCDscale_ggH": return "QCDscale_ggH"
        if self.name == "QCDscale_qqH": return "QCDscale_qqH"  #these evaluate to all dashes
        if self.name == "QCDscale_VH": return "QCDscale_VH"    #need to figure out how to add them
        if self.name == "QCDscale_ttH": return "QCDscale_ttH"  #to QCDscale_ggH properly
        if self.name == "QCDscale_ggVV": return "QCDscale_ggVV"
        if self.name == "QCDscale_VV": return "QCDscale_VV"
        if self.name == "BRhiggs_hzz4l": return "BRhiggs_hzz4l"
        if self.name == "CMS_eff_m": return "CMS_eff_m"
        if self.name == "CMS_eff_e": return "CMS_eff_e"
        for channel in channels:
            if self.name == "CMS_hzz{}_Zjets".format(channel): return "CMS_zz{}_zjets".format(channel)
        raise ValueError("Unknown systematic name {}!".format(self.name))

    def readsystematic(self):
        newname = self.newname()
        if newname is None: return None
        return ReadSystematics.singlesystematictype(self.production, newname)

def replacesystematics(channel, production):
    channel = Channel(channel)
    production = Production(production)
    cardfilename = "hzz4l_{}S_{}.txt".format(channel, production.year)
    with open(cardfilename) as f:
        contents = f.read()
    sections = contents.split("------------\n")
    systematicssection = sections[-1]
    systematicslines = systematicssection.split("\n")
    systematics = [Systematic(line, production) for line in systematicslines]

    readsystematics = ReadSystematics(channel, production)

    for i, systematic in enumerate(systematics[:]):
        readsys = systematic.readsystematic()
        if readsys is None: continue
        systematicslines[i] = readsys.getline(systematic.name, readsystematics)

    sections[-1] = systematicssection = "\n".join(systematicslines)
    contents = "------------\n".join(sections)
    with open(cardfilename, "w") as f:
        f.write(contents)
