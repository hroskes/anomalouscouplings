from collections import Counter
import os

from combinehelpers import Luminosity
from enums import Analysis, Category, Channel, MultiEnum, Production

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
                result.append("{} {}".format(label, getattr(obj, label)))
        return result

class SystematicsSection(Section):
    def getlines(self, obj, objtype):
        result = super(SystematicsSection, self).__get__(obj, objtype)
        for line in result[:]:
            if all(systematicvalue == "-" for systematicvalue in line.split()[1:]):
                result.remove(line)
        return result

class Datacard(MultiEnum):
    enums = (Analysis, Category, Channel, Production, Luminosity)
    @property
    def year(self):
        return self.production.year
    @property
    def txtfile(self):
        return "hzz4l_{}S_{}_{}.lumi{:f}.txt".format(self.channel, self.category, self.year, self.luminosity)
    @property
    def rootfile(self):
        return "hzz4l_{}S_{}_{}.lumi{:f}.input.root".format(self.channel, self.category, self.year, self.luminosity)

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

        bin = len(channels) * channels.index(self.channel) + categories.index(self.category)

        if counter[self] == 1:
            return bin
        elif counter[self] == 2:
            return " ".join([str(bin)]*len(productionmodes))
        assert False

    @property
    def observation(self):
        return getnobserved(self.channel, self.production, self.category, self.analysis)

    section3 = Section("bin", "observation")

    masswindowcomment = "## mass window [{},{}]".format(config.m4lmin, config.m4lmax)

    @property
    def process(self, counter=Counter()):
        counter[self] += 1
 
        if counter[self] == 1:
            return " ".join(_.name for _ in self.productionmodes)

        if counter[self] == 2:
            nsignal = len(_ for _ in self.productionmodes if _.issignal())
            nbkg = len(_ for _ in self.productionmodes if _.isbkg())
            assert nsignal+nbkg == len(self.productionmodes)

            return " ".join(str(_) for _ in range(-nsignal+1, nbkg+1))

        assert False

    @property
    def rate(self):
        return " ".join(str(getrate(p, self.channel, self.category, self.analysis, self.luminosity)) for p in self.productionmodes)

    section4 = Section("masswindowcomment", "bin", "process", "process", "rate")

    #just one dummy systematic for now
    @property
    def lumi_13TeV_common(self):
        return " ".join(["lnN"] + ["1.023" for p in self.productionmodes])

    section5 = SystematicsSection("lumi_13TeV_common")

    divider = "\n------------\n"

    def writedatacard(self):
        sections = self.section1, self.section2, self.section3, self.section4, self.section5
        if not os.path.exists(self.rootfile):
            raise IOError("workspace file {} should exist first!".format(self.rootfile))
        with open(self.txtfile, "w") as f:
            f.write(self.divider.join(sections)+"\n")

def writedatacard(*args):
    return Datacard(*args).writedatacard()
