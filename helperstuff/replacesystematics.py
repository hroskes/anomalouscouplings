import config
from enums import Channel, channels, EnumItem, MultiEnum, MyEnum, Production
from math import sqrt
import os
import yaml

class SystematicsType(MyEnum):
    enumitems = (
                 EnumItem("Moriond"),
                 EnumItem("ICHEP"),
                )
    enumname = "systematicstype"

class SingleNormalizationSystematic(object):
    def __init__(self, value):
        self.value = float(value)
    def __mul__(self, other):
        return type(self)(1 + (self.value-1)*other)
    def __div__(self, other):
        return self * (1.0/other)
    def __float__(self):
        return self.value
    def __str__(self):
        return str(float(self))

class UpDownSystematic(object):
    def __init__(self, valueorup, down=None):
        if down is None:
            try:
                self.up, self.down = (float(_) for _ in valueorup.split("/"))
            except ValueError:
                raise ValueError("{} is not a valid up/down systematic.  Has to be two floats separated by /.".format(value))
        else:
            self.up, self.down = float(valueorup), float(down)

    def __mul__(self, other):
        up = 1 + (self.up-1)*other
        down = 1 + (self.down-1)*other
        return type(self)(up, down)

    def __div__(self, other):
        return self * (1.0/other)

    def __str__(self):
        return "{}/{}".format(self.up, self.down)

def SystematicValue(value):
    try:
        return SingleNormalizationSystematic(value)
    except ValueError as e:
        try:
            return UpDownSystematic(value)
        except ValueError as e2:
            raise ValueError("Can't construct systematic from {}:\n\nTrying as a single value:\n{}\n\nTrying as up/down:\n{}".format(value, e.message, e2.message))

def ReadSystematics(*args, **kwargs):
    class _(MultiEnum):
        enums = (Production, Channel)
    year = _(*args).production.year
    return {2015: MoriondSystematics, 2016: ICHEPSystematics}[year](*args, **kwargs)

class ReadSystematics_BaseClass(MultiEnum):
    enums = (Production, Channel)
    def __init__(self, *args, **kwargs):
        super(ReadSystematics_BaseClass, self).__init__(*args)
        self.scenario = 1
        self.luminosity = self.production.dataluminosity
        for kw, kwarg in kwargs.iteritems():
            if kw == "scenario":
                self.scenario = kwarg
            elif kw == "luminosity":
                self.luminosity = kwarg
            else:
                raise TypeError("Bad kwarg {}={}".format(kw, kwarg))
        self.systematics = self.getsystematics()
    def check(self, *args):
        if self.production.year != self.year:
            raise ValueError("{} can only take productions from {}!".format(type(self).__name__, self.year))
        super(ReadSystematics_BaseClass, self).check(*args)
    def __getitem__(self, item):
        return self.systematics[item]

    def getline(self, oldline):
        if not oldline.strip(): return oldline
        useprocesses = ["ggH", "qqZZ", "ggZZ", "zjets"]

        oldname = oldline.split()[0]
        systype = oldline.split()[1]

        if systype != "lnN": return oldline

        #special cases here
        if oldname == "pdf_hzz4l_accept":
            uselist = [oldname, "lnN"] + ["1.02" if p == "ggH" else "-" for p in useprocesses]
            return " ".join(uselist)
        #if oldname == "lumi_8TeV":
            #https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
            #correlation = " ".join(["lumi_13TeV_common", "lnN"] + ["1.023" for p in useprocesses])
            #if self.year == 2015:
            #    mysys = " ".join(["lumi_13TeV_2015", "lnN"] + ["1.015" for p in useprocesses])
            #if self.year == 2016:
            #    mysys = " ".join(["lumi_13TeV_2016", "lnN"] + ["1.058" for p in useprocesses])
            #return "\n".join((correlation, mysys))
        #if "lumi_13TeV" in oldname:
            #return oldline

        thedict = {
                   "pdf_gg": "pdf_Higgs_gg",
                   "pdf_qqbar": "pdf_qq",
                   "QCDscale_ggH": "QCDscale_ggH",
                   "QCDscale_qqH": "QCDscale_qqH",  #these evaluate to all dashes
                   "QCDscale_VH": "QCDscale_VH",    #need to figure out how to add them
                   "QCDscale_ttH": "QCDscale_ttH",  #to QCDscale_ggH properly
                   "QCDscale_ggVV": "QCDscale_ggVV",
                   "QCDscale_ggVV_bonly": "QCDscale_ggVV_bonly",
                   "QCDscale_VV": "QCDscale_VV",
                   "BRhiggs_hzz4l": "BRhiggs_hzz4l",
                   "CMS_eff_m": "CMS_eff_m",
                   "CMS_eff_e": "CMS_eff_e",
                   "lumi_8TeV": "lumi_13TeV",
                  }
        thedict.update({"CMS_hzz{}_Zjets".format(channel): "CMS_zz{}_zjets".format(channel) for channel in channels})
        if oldname in thedict.values():
            name = oldname
        else:
            try:
                name = thedict[oldname]
            except KeyError:
                raise ValueError("Unknown systematic name {}!".format(oldname))

        return self.readline(name, systype, useprocesses)

class MoriondSystematics(ReadSystematics_BaseClass):
    year = 2015
    def __init__(self, *args, **kwargs):
        super(MoriondSystematics, self).__init__(*args, **kwargs)
        if self.scenario != 1:
            raise ValueError("Bad scenario {}".format(self.scenario))
        #for theory systematics
        #self.ichepsystematics = ICHEPSystematics([_ for _ in config.productionsforcombine if _.year == 2016][0], self.channel)
        self.ichepsystematics = ICHEPSystematics("160729", self.channel)

    def getsystematics(self):
        with open(self.channel.moriondcardfile()) as f:
            self.contents = f.read()

        sections = self.contents.split("------------------------------------------------------------")
        systematicssection = sections[-1]
        return {line.split()[0]: line for line in systematicssection.split("\n") if line.strip()}

    @staticmethod
    def getindex(processname):
        return "WH ZH ggH qqH ttH ggZZ qqZZ zjets".split().index(processname)

    def readline(self, name, systype, useprocesses):
        if "QCDscale" in name or "pdf_Higgs" in name:
            return self.ichepsystematics.readline(name, systype, useprocesses)

        useindices = [self.getindex(processname) for processname in useprocesses]

        line = self[name]
        split = line.split()
        if split[0] != name:
            raise ValueError("Wrong systematic passed to {}({}).getline():\n{}".format(type(self).__name__, name, line))
        if split[1] != "lnN":
            raise ValueError("Non-normalization systematic passed to {}.getline():\n{}".format(type(self).__name__, name, line))
        uselist = [name, split[1]] + [split[i+2] for i in useindices]
        return " ".join(uselist)

class ICHEPSystematics(ReadSystematics_BaseClass):
    year = 2016
    def __init__(self, *args, **kwargs):
        super(ICHEPSystematics, self).__init__(*args, **kwargs)
        if self.scenario not in (1, 2):
            raise ValueError("Bad scenario {}".format(self.scenario))

    def getsystematics(self):
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
            if "leptonScaleAndResol" in filename: continue

            if self.scenario == 2:
                for sysname, sys in y.iteritems():
                    for catname, category in sys.iteritems():
                        for productionmode, value in category.iteritems():
                            if productionmode == "type": continue
                            value = SystematicValue(value)
                            if "theory" in filename:
                                value /= 2
                            elif sysname.lower() in ("cms_eff_e", "cms_eff_m"):
                                if "2e2mu.yaml" in filename:
                                    value = SystematicValue(1.02)
                                else:
                                    value = SystematicValue(1.04)
                            elif sysname.endswith("_zjets") or sysname in ("JES", "bTagSF"):
                                value /= sqrt(self.luminosity / self.production.dataluminosity)
                            elif sysname == "lumi_13TeV":
                                pass
                            else:
                                raise ValueError("Don't know what do do with {}".format(sysname))
                            category[productionmode] = value

            self.systematics.update(y)
        return self.systematics

    def readline(self, name, systype, useprocesses):
        if name == "QCDscale_ggVV":   #included in ggH
            return None
        try:
            s = self[name]["UnTagged"]
        except KeyError:
            s = self[name]["Any"]
        if s["type"] != "lnN":
            raise ValueError("Non-normalization systematic passed to {}.getline():\n{}".format(type(self).__name__, name, line))
        uselist = [name, s["type"]] + [str(s[p]) if p in s else "-" for p in useprocesses]
        return " ".join(uselist)

def replacesystematics(channel, production, scenario=1, luminosity=None):
    channel = Channel(channel)
    production = Production(production)
    cardfilename = "hzz4l_{}S_{}.txt".format(channel, production.year)
    with open(cardfilename) as f:
        contents = f.read()
    sections = contents.split("------------\n")
    systematicssection = sections[-1]
    systematicslines = systematicssection.split("\n")

    readsystematics = ReadSystematics(channel, production, scenario=scenario, luminosity=luminosity)

    for i, systematic in enumerate(systematicslines[:]):
        systematicslines[i] = readsystematics.getline(systematic)

    sections[-1] = systematicssection = "\n".join(line for line in systematicslines if line)
    contents = "------------\n".join(sections)
    with open(cardfilename, "w") as f:
        f.write(contents)
