#!/usr/bin/env python

from collections import OrderedDict
from itertools import izip
import os
import shutil
import subprocess

from helperstuff.enums import analyses
from helperstuff.utilities import cd, mkdir_p
from helperstuff import config

from limits import getlimits

class Model(object):
  def __init__(self, name, exp68, exp68_up, exp95, exp95_2, exp95_up_2, exp95_up, bestfit, obs68, obs95, obs95_2, obs68_up, obs95_up, obs95_up_2):
    self.name = name
    self.exp68 = exp68
    self.exp68_up = exp68_up
    self.exp95 = exp95
    self.exp95_2 = exp95_2
    self.exp95_up_2 = exp95_up_2
    self.exp95_up = exp95_up
    self.bestfit = bestfit
    self.obs68 = obs68
    self.obs95 = obs95
    self.obs95_2 = obs95_2
    self.obs68_up = obs68_up
    self.obs95_up = obs95_up
    self.obs95_up_2 = obs95_up_2

  @property
  def stuff(self):
    return (
      self.name,
      self.exp68,
      self.exp68_up,
      self.exp95,
      self.exp95_2,
      self.exp95_up_2,
      self.exp95_up,
      self.bestfit,
      self.obs68,
      self.obs95,
      self.obs95_2,
      self.obs68_up,
      self.obs95_up,
      self.obs95_up_2,
    )

  def __iter__(self):
    return iter(self.stuff)

def setupdats(plotid):
  with cd(os.path.join(config.repositorydir, "step10_plottingutilities", "NIS_summary", "Input_"+plotid)):
    files = Model(None, None, None, None, None, None, None, None, None, None, None, None, None, None)
    models = OrderedDict()
    with \
         open("model_names.dat") as files.name, \
         open("model_1sig.dat") as files.exp68, \
         open("model_1sig_up.dat") as files.exp68_up, \
         open("model_95CL.dat") as files.exp95, \
         open("model_95CL_2.dat") as files.exp95_2, \
         open("model_95CL_up_2.dat") as files.exp95_up_2, \
         open("model_95CL_up.dat") as files.exp95_up, \
         open("obs_BF.dat") as files.bestfit, \
         open("obs_model_1sig.dat") as files.obs68, \
         open("obs_model_95CL.dat") as files.obs95, \
         open("obs_model_95CL_2.dat") as files.obs95_2, \
         open("obs_model_1sig_up.dat") as files.obs68_up, \
         open("obs_model_95CL_up.dat") as files.obs95_up, \
         open("obs_model_95CL_up_2.dat") as files.obs95_up_2:
      for data in izip(*files):
        data = [data[0].strip()] + [float(_) if float(_) != 0 else int(float(_)) for _ in data[1:]]
        assert data[0] not in models, data[0]
        models[data[0]] = Model(*data)

    def thekey(kv):
      k, v = kv
      if "Z#gamma" in k: k = "yy"+k
      if "#gamma#gamma" in k: k = "zz"+k
      k = k.replace("3", "1")
      k = k.replace("#", "")
      k = k.lower()
      return k
    models = OrderedDict((k, v) for (k, v) in sorted(models.items(), key=thekey))

    for analysis in analyses:
      filename = os.path.join(config.plotsbasedir, "limits", "{}_September7combination".format(analysis), "limit_lumi77.45_7813_101,-1.0,1.0_101,-0.02,0.02.root")
      allresults = getlimits(filename, poi="CMS_zz4l_fai1", domirror=analysis=="fa3")

      m = models[analysis.couplingtitle]
      expname = "Expected, {} = 0 or #pi".format(analysis.phi)
      obsname = "Observed, {} = 0 or #pi".format(analysis.phi)
      (expmin, exp68, exp95), (obsmin, obs68, obs95) = allresults[expname], allresults[obsname]

      assert abs(expmin) < .001, expmin

      if len(exp68) == 1:
          m.exp68, m.exp68_up = -exp68[0][0], exp68[0][1]
      else:
          assert False, exp68

      if len(exp95) == 1:
         m.exp95, m.exp95_up = exp95[0][0], exp95[0][1]  #different sign convention than 68
         m.exp95_2 = m.exp95_up_2 = 0
      elif len(exp95) == 2:
         m.exp95, m.exp95_2, m.exp95_up_2, m.exp95_up = exp95[0][0], exp95[0][1], exp95[1][0], exp95[1][1]
      else:
          assert False, exp95

      m.bestfit = obsmin

      if len(obs68) == 1:
          m.obs68, m.obs68_up = obsmin-obs68[0][0], obs68[0][1]-obsmin
      else:
          assert False, obs68

      if len(obs95) == 1:
         m.obs95, m.obs95_up = obs95[0][0], obs95[0][1]  #different sign convention than 68
         m.obs95_2 = m.obs95_up_2 = 0
      elif len(obs95) == 2:
         m.obs95, m.obs95_2, m.obs95_up_2, m.obs95_up = obs95[0][0], obs95[0][1], obs95[1][0], obs95[1][1]
      else:
          assert False, obs95

    with \
         open("model_names.dat", "w") as files.name, \
         open("model_1sig.dat", "w") as files.exp68, \
         open("model_1sig_up.dat", "w") as files.exp68_up, \
         open("model_95CL.dat", "w") as files.exp95, \
         open("model_95CL_2.dat", "w") as files.exp95_2, \
         open("model_95CL_up_2.dat", "w") as files.exp95_up_2, \
         open("model_95CL_up.dat", "w") as files.exp95_up, \
         open("obs_BF.dat", "w") as files.bestfit, \
         open("obs_model_1sig.dat", "w") as files.obs68, \
         open("obs_model_95CL.dat", "w") as files.obs95, \
         open("obs_model_95CL_2.dat", "w") as files.obs95_2, \
         open("obs_model_1sig_up.dat", "w") as files.obs68_up, \
         open("obs_model_95CL_up.dat", "w") as files.obs95_up, \
         open("obs_model_95CL_up_2.dat", "w") as files.obs95_up_2:
      for model in models.values():
        for data, f in zip(model, files):
          f.write(str(data).strip()+"\n")

def makeplot(plotid):
  with cd(os.path.join(config.repositorydir, "step10_plottingutilities", "NIS_summary")):
    subprocess.check_call(["make"])
    subprocess.check_call(["./NIS_summary_4", plotid])
    for ext in "png eps root pdf C".split():
      if plotid == "HIG18002onshell": shutil.copy("Summary_HIG18002onshell."+ext, os.path.join(config.plotsbasedir, "limits", "summary."+ext))
      if plotid == "HIG18002onshellPAS": shutil.copy("Summary_HIG18002onshellPAS."+ext, os.path.join(config.plotsbasedir, "limits", "forPAS", "summary."+ext))
      if plotid == "HIG18002onshelloffshell": shutil.copy("Summary_HIG18002onshelloffshell."+ext, os.path.join(config.plotsbasedir, "limits", "summary_offshell."+ext))
      if plotid == "HIG18002onshelloffshellPAS": shutil.copy("Summary_HIG18002onshelloffshellPAS."+ext, os.path.join(config.plotsbasedir, "limits", "forPAS", "summary_offshell."+ext))

if __name__ == "__main__":
#  setupdats("HIG18002onshell")
#  makeplot("HIG18002onshell")
#  setupdats("HIG18002onshellPAS")
#  makeplot("HIG18002onshellPAS")
#  setupdats("HIG18002onshelloffshell")
  makeplot("HIG18002onshelloffshell")
  makeplot("HIG18002onshelloffshellPAS")
