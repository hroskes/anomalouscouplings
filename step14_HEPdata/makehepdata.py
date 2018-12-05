#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("cadi", choices=["HIG-17-011"])
  g = p.add_mutually_exclusive_group()
  g.add_argument("--force", "-f", action="store_true")
  g.add_argument("--check", "-c", action="store_true")
  args = p.parse_args()

import abc, contextlib, os, shutil, subprocess, tempfile, urllib2

import yaml

class HEPData(object):
  def run(self, force=False, check=False):
    if check:
      self.check(self.cadi)
      return
    if os.path.exists(self.cadi):
      if not force:
        raise RuntimeError(self.cadi+" exists, try removing it or running with -f")
      shutil.rmtree(self.cadi)
    os.mkdir(self.cadi)
    with cd(self.cadi):
      for table in self.tables:
        table.writeyaml()
        table.getplots()
      with open("submission.yaml", "w") as f:
        yaml.dump_all(self.submission, f, default_flow_style=False)
    self.check(self.cadi)

  @classmethod
  def check(cls, directory):
    directory = os.path.abspath(directory)
    with cd(tempfile.mkdtemp()):
      with open("check.py", "w") as newf, contextlib.closing(urllib2.urlopen("https://github.com/HEPData/hepdata-submission/raw/master/scripts/check.py")) as oldf:
        shutil.copyfileobj(oldf, newf)
      output = subprocess.check_output(["python", "check.py", directory])
      print output
      if "error" in output: raise RuntimeError("check script failed")

  @property
  def submission(self):
    submission = [{"comment": "Constraints on anomalous Higgs boson couplings using production and decay information in the four-lepton final state"}]
    for table in self.tables:
      submission.append(table.submissionyaml)
    return submission

  @abc.abstractproperty
  def tables(self): pass
  @abc.abstractproperty
  def cadi(self): pass

class PaperInfo(object):
  def __init__(self, cadi):
    self.__cadi = cadi
  @property
  def cadi(self):
    return self.__cadi
  @property
  def reactions(self):
    return {"HIG-17-011": ["P P --> HIGGS --> Z0 Z0 --> LEPTON+ LEPTON- LEPTON+ LEPTON-"]}[self.cadi]
  @property
  def sqrts(self):
    return {"HIG-17-011": [13000]}[self.cadi]
  @property
  def rootfile(self):
    return {"HIG-17-011": "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/limits/{self.fai}_limit_lumi35.8671.root"}[self.cadi].format(self=self)
  @property
  def lumidict(self):
    return {"HIG-17-011": {7: 5.1, 8: 19.7, 13: 38.6}}[self.cadi]

class Table(PaperInfo):
  __metaclass__ = abc.ABCMeta
  @property
  def figurename(self):
    result = "CMS-" + self.cadi + "_Figure_{:03d}".format(self.figurenumber)
    if self.figureletter is not None: result += "-"+self.figureletter
    result += ".png"
    return result
  @property
  def thumbfigurename(self):
    return "thumb_"+self.figurename
  @property
  def figureurl(self):
    return "http://cms-results.web.cern.ch/cms-results/public-results/publications/"+self.cadi+"/"+self.figurename
  @property
  def location(self):
    return "Fig. {:d}{} on page {:d} of the paper".format(self.figurenumber, self.figureletter if self.figureletter else None, self.figurepage)
  @abc.abstractproperty
  def figurepage(self): pass
  @property
  def submissionyaml(self):
    return {
      "additional_resources": [
         {"description": "Image file", "location": self.figurename},
         {"description": "Thumbnail image file", "location": self.thumbfigurename},
      ],
      "data_file": self.yamlfile,
      "description": self.description,
      "keywords": [
        {"name": "reactions", "values": self.reactions},
        {"name": "observables", "values": self.observables},
        {"name": "cmenergies", "values": self.sqrts},
      ],
      "location": self.location,
      "name": self.name,
    }

  @abc.abstractproperty
  def yaml(self):
    pass
  @abc.abstractproperty
  def name(self):
    pass
  @abc.abstractproperty
  def figurenumber(self):
    pass
  @abc.abstractproperty
  def figureletter(self):
    pass
  def writeyaml(self):
    with open(self.yamlfile, "w") as f:
      yaml.dump(self.yaml, f, default_flow_style=False)
  @property
  def yamlfile(self):
    return self.name.replace(" ", "")
  def getplots(self):
    try:
      with open(self.figurename, "w") as newf, contextlib.closing(urllib2.urlopen(self.figureurl)) as oldf:
        shutil.copyfileobj(oldf, newf)
    except urllib2.HTTPError:
      print self.figureurl
      raise
    subprocess.check_call(["convert", "-thumbnail", "200", self.figurename, self.thumbfigurename])
    assert os.path.exists(self.thumbfigurename)
  @abc.abstractproperty
  def description(self): pass

class FaiLimitTable(Table):
  def __init__(self, cadi, fai):
    self.fai = fai
    super(FaiLimitTable, self).__init__(cadi)
  @property
  def faifancy(self):
    return {"fa3": "f_{a3}", "fa2": "f_{a2}", "fL1": "f_{\Lambda1}", "fL1Zg": "f_{\Lambda1}^{Z\gamma}"}[self.fai]
  @property
  def phiaifancy(self):
    return self.faifancy.replace("f", "\phi", 1)
  @property
  def aifancy(self):
    return {"fa3": "a_{3}", "fa2": "a_{2}", "fL1": "\Lambda_{1}", "fL1Zg": "\Lambda_{1}^{Z\gamma}"}[self.fai]
  @property
  def xheader(self):
    return {"name": "${}\cos({})$".format(self.faifancy, self.phiaifancy)}
  @property
  def observables(self):
    return ["SIG/SIG"]
  @property
  def yaml(self):
    import ROOT
    x = {"header": self.xheader}
    ys = []
    result = {"name": self.yamlfile, "independent_variables": [x], "dependent_variables": ys}
    xvalues = None
    f = ROOT.TFile(self.rootfile)
    canvas = f.GetListOfKeys()
    assert len(canvas) == 1
    canvas = canvas[0]
    canvas = canvas.ReadObj()
    pad = canvas.GetListOfPrimitives()[0]
    legend = pad.GetListOfPrimitives()[-1]
    for entry in self.legendentries(legend):
      if xvalues is None:
        xvalues = entry.xvalues
        x["values"] = [{"value": _} for _ in xvalues]
      else:
        assert xvalues == entry.xvalues, (set(xvalues) ^ set(entry.xvalues))
    ys.append(entry.ydict)
    return result

  def legendentries(self, legend):
    for entry in legend.GetListOfPrimitives():
      yield FaiLegendEntry(entry, self.lumidict, self.aifancy, self.removepoints)

  @property
  def name(self):
    return self.fai+".yaml"

  @property
  def removepoints(self):
    if self.cadi == "HIG-17-011":
      if self.fai == "fa3": return {-6.821699116699165e-06, 7.081189892232942e-07, 1.2522907582024345e-06, -4.0629518480272964e-05, -3.221497172489762e-05, 0.008799999952316284, 0.006000000052154064, 0.009999999776482582, 0.011599999852478504, 0.012400000356137753, 0.015599999576807022, 0.012000000104308128, 0.00800000037997961, 0.00559999980032444, 0.007199999876320362, 0.010400000028312206, 0.012799999676644802, 0.006800000090152025, 0.013199999928474426, 0.009600000455975533, 0.01080000028014183, 0.006399999838322401, 0.01360000018030405, 0.01119999960064888, 0.01600000075995922, 6.609242291233386e-07, 0.005200000014156103, 1.3594632264357642e-06, 0.00839999970048666, 0.007600000128149986, 0.009200000204145908}
      if self.fai == "fa2": return {0.007500944659113884, 0.012208748608827591, -6.613622982598599e-09, 3.1140121592443393e-09, 0.0242480281740427, 0.012350913137197495, 1.719479225670284e-09, -5.649837819809989e-10}
      if self.fai == "fL1": return {0.0001509769936092198, 9.811160452954937e-06, 0.00015089577937033027, 7.448184291547477e-10, 6.218647467903793e-05, 6.269038567552343e-05, 6.28244652034482e-06, 9.799659892451018e-06}
      if self.fai == "fL1Zg": return {0.00014273269334807992, -0.00011004077532561496, -1.194288233818952e-05, -9.83606241788948e-06, 0.00010419703176012263, 1.7338144971290603e-05, -0.00010066330287372693, 8.42108711367473e-05}
    return set()

  @property
  def figurenumber(self):
    if self.cadi == "HIG-17-011": return 3
    assert False
  @property
  def figureletter(self):
    return next(letter for letter, fai in zip("abcd", ("fa3", "fa2", "fL1", "fL1Zg")) if fai == self.fai)

  @property
  def description(self):
    return "Observed and expected likelihood scans ${self.fancyfai}\cos{self.fancyphiai}$.  See Section 2 of the paper for more details."

  @property
  def figurepage(self):
    if self.cadi == "HIG-17-011": return 6
    assert False

class FaiLegendEntry(object):
  def __init__(self, entry, lumidict, aifancy, removepoints):
    self.entry = entry
    self.lumidict = lumidict
    self.aifancy = aifancy
    self.removepoints = removepoints
  @property
  def tgraph(self):
    return self.entry.GetObject()
  @property
  def header(self):
    return {"name": "-2\Delta\ln L "+self.entry.GetLabel()}
  @property
  def xyvalues(self):
    return [(x, y) for i, x, y in zip(range(self.tgraph.GetN()), self.tgraph.GetX(), self.tgraph.GetY()) if x not in self.removepoints]
  @property
  def xvalues(self):
    return [x for x, y in self.xyvalues]
  @property
  def yvalues(self):
    return [y for x, y in self.xyvalues]
  @property
  def ydict(self):
    return {
      "header": self.header,
      "qualifiers": [
        {"name": "LUMINOSITY_{:d}TeV".format(sqrts), "units": 'fb$^{-1}$', "value": lumi}
          for sqrts, lumi in self.lumidict.iteritems() if not ("13 TeV" in self.header and sqrts != 13)
      ] + [
        {"name": "HVV couplings", "value": "$%s$ real, other $f_{ai}=0$" % self.aifancy}
      ],
      "values": [{"value": y} for y in self.yvalues]
    }

class HIG17011(HEPData):
  @property
  def cadi(self):
    return "HIG-17-011"
  @property
  def tables(self):
    for fai in "fa3", "fa2", "fL1", "fL1Zg":
      yield FaiLimitTable(self.cadi, fai)

@contextlib.contextmanager
def cd(newdir):
    """http://stackoverflow.com/a/24176022/5228524"""
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

if __name__ == "__main__":
  globals()[args.cadi.replace("-", "")]().run(force=args.force, check=args.check)
