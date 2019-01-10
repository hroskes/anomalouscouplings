#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("cadi", choices=["HIG-17-011", "HIG-18-002", "HIG-17-034"])
  g = p.add_mutually_exclusive_group()
  g.add_argument("--force", "-f", action="store_true")
  g.add_argument("--check", "-c", action="store_true")
  g.add_argument("--tar", "-t", action="store_true")
  args = p.parse_args()

import abc, contextlib, os, shutil, subprocess, tempfile, urllib2

import yaml

class HEPData(object):
  def run(self, force=False, check=False, tar=False):
    if check or tar:
      self.check(self.cadi)
      if tar:
        subprocess.check_call(["tar", "cvzf", self.cadi+".tgz", "-C", self.cadi] + os.listdir(self.cadi))
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
    return {
      "HIG-17-011": ["P P --> HIGGS --> Z0 Z0 --> LEPTON+ LEPTON- LEPTON+ LEPTON-"],
      "HIG-18-002": ["P P --> Z0 Z0 --> LEPTON+ LEPTON- LEPTON+ LEPTON-"],
      "HIG-17-034": ["P P --> HIGGS --> TAU+ TAU-"],
    }[self.cadi]
  @property
  def sqrts(self):
    return {
      "HIG-17-011": [13000],
      "HIG-18-002": [13000],
      "HIG-17-034": [13000],
    }[self.cadi]
  @property
  def rootfile(self):
    return {
      "HIG-17-011": "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/limits/{self.fai}_limit_lumi35.8671.root",
      "HIG-18-002": "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_HIG18002/limits/{self.fai}_limit_lumi80.15.root",
      "HIG-17-034": "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_HIG18002/limits/HTT/{self.fai}_limit_lumi80.15.root",
    }[self.cadi].format(self=self)
  @property
  def lumidict(self):
    return {
      "HIG-17-011": {7: 5.1, 8: 19.7, 13: 38.6},
      "HIG-18-002": {7: 5.1, 8: 19.7, 13: 80.2},
      "HIG-17-034": {7: 5.1, 8: 19.7, 13: 80.2},
    }[self.cadi]

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
  def publishedinjournal(self):
    if self.cadi == "HIG-17-034": return False
    if self.cadi == "HIG-18-002": return False
    return True
  @property
  def publishedonarxiv(self):
    if self.cadi == "HIG-17-034": return False
    return True
  @property
  def figureurl(self):
    if not self.publishedonarxiv:
      result = "file:/afs/cern.ch/work/h/hroskes/AN/papers/{self.cadi}/trunk/Figure_{self.figurenumber:03d}".format(self=self)
      if self.figureletter is not None: result += "-"+self.figureletter
      result += ".pdf"
      return result
    return "http://cms-results.web.cern.ch/cms-results/public-results/publications/{self.cadi}/{self.figurename}".format(self=self)
  @property
  def location(self):
    if self.publishedinjournal:
      return "Fig. {:d}{} on page {:d} of the paper".format(self.figurenumber, self.figureletter if self.figureletter else None, self.figurepage)
    else:
      return "Fig. {:d}{} of the paper".format(self.figurenumber, self.figureletter if self.figureletter else None)
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
    figureurl = self.figureurl
    if figureurl.endswith(".pdf"):
      try:
        with tempfile.NamedTemporaryFile(suffix=".pdf", bufsize=0) as newf, contextlib.closing(urllib2.urlopen(figureurl)) as oldf:
          shutil.copyfileobj(oldf, newf)
          subprocess.check_call(["convert", "-density", "150", "-trim", newf.name, "-quality", "100", self.figurename])
      except urllib2.HTTPError:
        print self.figureurl
        raise
    elif figureurl.endswith(".png"):
      try:
        with open(self.figurename, "w") as newf, contextlib.closing(urllib2.urlopen(figureurl)) as oldf:
          shutil.copyfileobj(oldf, newf)
      except urllib2.HTTPError:
        print self.figureurl
        raise
    else:
      assert False, self.figureurl
    subprocess.check_call(["convert", "-thumbnail", "200", self.figurename, self.thumbfigurename])
    assert os.path.exists(self.thumbfigurename)
  @abc.abstractproperty
  def description(self): pass

class FaiLimitTable(Table):
  def __init__(self, cadi, fai):
    self.fai = fai
    super(FaiLimitTable, self).__init__(cadi)
  @property
  def fancyfai(self):
    return {"fa3": "f_{a3}", "fa2": "f_{a2}", "fL1": "f_{\Lambda1}", "fL1Zg": "f_{\Lambda1}^{Z\gamma}"}[self.fai]
  @property
  def fancyphiai(self):
    return self.fancyfai.replace("f", "\phi", 1)
  @property
  def fancyai(self):
    return {"fa3": "a_{3}", "fa2": "a_{2}", "fL1": "\Lambda_{1}", "fL1Zg": "\Lambda_{1}^{Z\gamma}"}[self.fai]
  @property
  def xheader(self):
    return {"name": "${}\cos({})$".format(self.fancyfai, self.fancyphiai)}
  @property
  def observables(self):
    return ["SIG/SIG"]

  @abc.abstractmethod
  def getlegendfromrootfile(self, tfile):
    pass

  @property
  def yaml(self):
    import ROOT
    x = {"header": self.xheader}
    ys = []
    result = {"name": self.yamlfile, "independent_variables": [x], "dependent_variables": ys}
    xvalues = None
    f = ROOT.TFile(self.rootfile)
    legend = self.getlegendfromrootfile(f)
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
      yield FaiLegendEntry(entry, self.lumidict, self.fancyphiai, self.removepoints)

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
    if self.cadi == "HIG-18-002":
      if self.fai == "fa3": return {-3.469446951953614e-17, -5.048824459663592e-05, -7.212763011921197e-05, -4.85241471324116e-05, -5.592270326815196e-07, 0.1599999964237213, -0.1599999964237213}
      if self.fai == "fa2": return {0.0003515805583447218, -1.5517587215185813e-12, -2.256646258746997e-10, 0.019600000232458115, 0.00041023208177648485, 1.6157702702912502e-05, 6.066570676921401e-07, -2.392181785992875e-09, 2.2096777740898688e-09, 0.5600000023841858, 0.01600000075995922}
      if self.fai == "fL1": return {3.816104435827583e-05, 9.086481567166516e-15, -0.3400000035762787, -0.6800000071525574, -0.3199999928474426, 4.4846903620054945e-05, 6.056595225345518e-10, 6.253050059967791e-07, 1.2750886526191607e-05, -5.451978868364904e-10, 2.5288371396925413e-10}
      if self.fai == "fL1Zg": return {0.7200000286102295, 0.7799999713897705, 5.376144329716226e-08, -0.3799999952316284, -0.0002512808714527637, 3.570632856053635e-08, -0.0002523555886000395, -0.8799999952316284, 1.2302773484407226e-06, 2.031392796197906e-05, -0.00026075352798216045, 0.800000011920929, -0.699999988079071, -0.0002693945134524256}
    if self.cadi == "HIG-17-034":
      if self.fai == "fa3": return {-3.248867841421088e-11, 1.760848959975192e-07, 0.0012000000569969416, 0.6800000071525574, 0.017999999225139618, -0.014000000432133675, 0.0020000000949949026, 0.3799999952316284, 1.6212986508890026e-07, -1.9253953986719807e-10, -7.212763011921197e-05, -5.048824459663592e-05, -3.469446951953614e-17, 5.018843012294383e-07, 4.964784920957754e-07, 0.012400000356137753, -5.114797474448096e-13, -2.177543423353967e-10}
      if self.fai == "fa2": return {7.216371159302071e-05, 8.422934479312971e-05, -0.018400000408291817, -0.30000001192092896, 1.234568003383174e-13, -5.851388817923464e-12, -0.014800000004470348, 0.005200000014156103, -0.4399999976158142, 0.007199999876320362, 0.00041023208177648485, 0.019600000232458115, 0.0003515805583447218, -2.256646258746997e-10, -1.5517587215185813e-12, 0.7799999713897705, 0.006800000090152025, 2.4923905584728345e-05}
      if self.fai == "fL1": return {5.741680979554076e-07, 1.7215561456396244e-06, -0.07999999821186066, -0.01080000028014183, 0.10000000149011612, 7.998490758609478e-11, 0.4399999976158142, 0.4000000059604645, 0.6200000047683716, 0.6000000238418579, -1.2600177845545346e-11, 0.3799999952316284, 0.5600000023841858, 0.3400000035762787, -0.8600000143051147, 3.816104435827583e-05, -0.6800000071525574, 4.4846903620054945e-05, 9.086481567166516e-15, -0.3400000035762787, -0.3199999928474426, 6.056595225345518e-10, -0.3799999952316284, 0.9200000166893005, -0.5799999833106995, -0.00839999970048666, -0.1599999964237213, -0.6600000262260437, 0.11999999731779099, 4.1646595150268695e-07, 0.009200000204145908}
      if self.fai == "fL1Zg": return {-4.05151867610698e-09, -1.4713416931044776e-05, -3.2008043490350246e-05, -0.3799999952316284, 3.570632856053635e-08, 5.376144329716226e-08, 0.7200000286102295, 0.7799999713897705, -0.0002512808714527637, -0.0002523555886000395, -0.009600000455975533, -0.003599999938160181, 7.022480303930934e-07, -0.014800000004470348, 0.002400000113993883, 0.699999988079071, -0.03999999910593033, -4.3368086899420177e-16}
    return set()

  @property
  def figurenumber(self):
    return {
      "HIG-17-011": 3,
      "HIG-18-002": 4,
      "HIG-17-034": 10,
    }[self.cadi]
    assert False
  @property
  def figureletter(self):
    return next(letter for letter, fai in zip("abcd", ("fa3", "fa2", "fL1", "fL1Zg")) if fai == self.fai)

  @property
  def description(self):
    return "Observed and expected likelihood scans of ${self.fancyfai}\cos{self.fancyphiai}$.  See Section {self.paperphenosection} of the paper for more details.".format(self=self)

  @property
  def figurepage(self):
    if self.cadi == "HIG-17-011": return 6
    assert False

  @property
  def paperphenosection(self):
    return {
      "HIG-17-011": 2,
      "HIG-17-034": 4,
    }[self.cadi]

class FaiLimitTable_TwoPartPlot(FaiLimitTable):
  def getlegendfromrootfile(self, tfile):
    import ROOT
    canvas = tfile.GetListOfKeys()
    assert len(canvas) == 1
    canvas = canvas[0]
    canvas = canvas.ReadObj()
    pad = canvas.GetListOfPrimitives()[0]
    legend = pad.GetListOfPrimitives()[-1]
    assert isinstance(legend, ROOT.TLegend), legend
    return legend

class FaiLimitTable_SixPartPlot(FaiLimitTable):
  def getlegendfromrootfile(self, tfile):
    import ROOT
    canvas = tfile.GetListOfKeys()
    assert len(canvas) == 1
    canvas = canvas[0]
    canvas = canvas.ReadObj()
    legend = canvas.GetListOfPrimitives()[7]
    assert isinstance(legend, ROOT.TLegend), legend
    return legend

  @property
  def figureurl(self):
    if self.cadi == "HIG-17-034":
      if os.path.exists("/afs/cern.ch/work/h/hroskes/AN/papers/{self.cadi}/trunk/Figure_010-a.pdf".format(self=self)) or self.publishedonarxiv:
        raise RuntimeError("Remove this function")
      return "file:/afs/cern.ch/work/h/hroskes/AN/papers/{self.cadi}/trunk/Figures/{self.fai}_limit_lumi80.pdf".format(self=self)
    return super(FaiLimitTable_SixPartPlot, self).figureurl

class FaiLimitTable_HIG18002(FaiLimitTable_SixPartPlot):
  @property
  def name(self):
    return self.fai+"_onshell.yaml"

  @property
  def description(self):
    return "Observed and expected likelihood scans of ${self.fancyfai}\cos{self.fancyphiai}$ using on-shell events only.  See Section 2 of the paper for more details.".format(self=self)

class FaiLegendEntry(object):
  def __init__(self, entry, lumidict, fancyphiai, removepoints):
    self.entry = entry
    self.lumidict = lumidict
    self.fancyphiai = fancyphiai
    self.removepoints = removepoints
  @property
  def tgraph(self):
    return self.entry.GetObject()
  @property
  def header(self):
    return {"name": "$-2\Delta\ln L$ "+self.entry.GetLabel()}
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
          for sqrts, lumi in self.lumidict.iteritems() if not ("13 TeV" in self.header["name"] and sqrts != 13)
      ] + [
        {"name": "HVV couplings", "value": "${self.fancyphiai}=0$ or $\pi$, other $f_{{ai}}=0$".format(self=self)}
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
      yield FaiLimitTable_TwoPartPlot(self.cadi, fai)

class HIG18002(HEPData):
  @property
  def cadi(self):
    return "HIG-18-002"
  @property
  def tables(self):
    for fai in "fa3", "fa2", "fL1", "fL1Zg":
      yield FaiLimitTable_HIG18002(self.cadi, fai)

class HIG17034(HEPData):
  @property
  def cadi(self):
    return "HIG-17-034"
  @property
  def tables(self):
    for fai in "fa3", "fa2", "fL1", "fL1Zg":
      yield FaiLimitTable_SixPartPlot(self.cadi, fai)

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
  globals()[args.cadi.replace("-", "")]().run(force=args.force, check=args.check, tar=args.tar)
