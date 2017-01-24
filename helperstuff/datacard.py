from collections import Counter

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
#            elif "," in label:
#                name, args = label.split(",", 1)
#                args = args.split(",")
#                result += "{} {}".format(name, getattr(obj, name)(*args))
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
    def jmax(self):
        return len(self.productionmodes)-1
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
        raise NotImplementedError

    section3 = Section("bin", "observation")

    masswindowcomment = "## mass window [105.0,140.0]"

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
        raise NotImplementedError

    section4 = Section("masswindowcomment", "bin", "process", "process", "rate")

    section5 = SystematicsSection()

    divider = "\n------------\n"

    def writedatacard(self):
        sections = self.section1, self.section2, self.section3, self.section4, self.section5
        with open(self.txtfile, "w") as f:
            f.write(self.divider.join(sections)+"\n")
