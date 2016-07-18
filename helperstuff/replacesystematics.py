from enums import Channel, channels

class MoriondSystematics(object):
    def __init__(self, channel):
        self.channel = Channel(channel)
        with open(self.channel.moriondcardfile()) as f:
            self.contents = f.read()

        sections = self.contents.split("------------------------------------------------------------")
        systematicssection = sections[-1]
        self.systematics = {line.split()[0]: line for line in systematicssection.split("\n") if line.strip()}

    def __getitem__(self, item):
        return self.systematics[item]

class MoriondSystematic(object):
    def __init__(self, name, useprocesses=None):
        self.name = name
        if useprocesses is None:
            useprocesses = ["ggH", "qqZZ", "ggZZ", "zjets"]
        self.useprocesses = useprocesses
        self.useindices = [self.getindex(processname) for processname in useprocesses]

    @staticmethod
    def getindex(processname):
        return "WH ZH ggH qqH ttH ggZZ qqZZ zjets".split().index(processname)

    def readline(self, name, moriondsystematics):
        line = moriondsystematics[self.name]
        split = line.split()
        if split[0] != self.name:
            raise ValueError("Wrong systematic passed to {}({}).readline():\n{}".format(type(self).__name__, self.name, line))
        if split[1] != "lnN":
            raise ValueError("Non-normalization systematic passed to {}.readline():\n{}".format(type(self).__name__, self.name, line))
        uselist = [name, split[1]] + [split[i+2] for i in self.useindices]
        return " ".join(uselist)

class Systematic(object):
    def __init__(self, line):
        self.name = self.type = None
        try:
            self.name = line.split()[0]
            self.type = line.split()[1]
        except IndexError:
            pass
        self.line = line

    def moriondsystematic(self):
        if self.type != "lnN": return None
        if self.name == "lumi_8TeV": return MoriondSystematic("lumi_13TeV")
        if self.name == "pdf_gg": return MoriondSystematic("pdf_Higgs_gg")
        if self.name == "pdf_qqbar": return MoriondSystematic("pdf_qq")
        if self.name == "pdf_hzz4l_accept": return None
        if self.name == "QCDscale_ggH": return MoriondSystematic("QCDscale_ggH")
        if self.name == "QCDscale_qqH": return MoriondSystematic("QCDscale_qqH")  #these evaluate to all dashes
        if self.name == "QCDscale_VH": return MoriondSystematic("QCDscale_VH")    #need to figure out how to add them
        if self.name == "QCDscale_ttH": return MoriondSystematic("QCDscale_ttH")  #to QCDscale_ggH properly
        if self.name == "QCDscale_ggVV": return MoriondSystematic("QCDscale_ggVV")
        if self.name == "QCDscale_VV": return MoriondSystematic("QCDscale_VV")
        if self.name == "BRhiggs_hzz4l": return MoriondSystematic("BRhiggs_hzz4l")
        if self.name == "CMS_eff_m": return MoriondSystematic("CMS_eff_m")
        if self.name == "CMS_eff_e": return MoriondSystematic("CMS_eff_e")
        for channel in channels:
            if self.name == "CMS_hzz{}_Zjets".format(channel): return MoriondSystematic("CMS_zz{}_zjets".format(channel))
        raise ValueError("Unknown systematic name {}!".format(self.name))

def replacesystematics(channel):
    channel = channel
    cardfilename = "hzz4l_{}S_8TeV.txt".format(channel)
    with open(cardfilename) as f:
        contents = f.read()
    sections = contents.split("------------\n")
    systematicssection = sections[-1]
    systematicslines = systematicssection.split("\n")
    systematics = [Systematic(line) for line in systematicslines]

    moriondsystematics = MoriondSystematics(channel)

    for i, systematic in enumerate(systematics[:]):
        moriond = systematic.moriondsystematic()
        if moriond is None: continue
        systematicslines[i] = moriond.readline(systematic.name, moriondsystematics)

    sections[-1] = systematicssection = "\n".join(systematicslines)
    contents = "------------\n".join(sections)
    with open(cardfilename, "w") as f:
        f.write(contents)
