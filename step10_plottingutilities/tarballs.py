#!/usr/bin/env python

import glob
import os
import shutil
import subprocess

from helperstuff import config
from helperstuff.enums import analyses, categories

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

if os.path.exists(os.path.join(config.plotsbasedir, "templateprojections", "projections")):
  tmpdir = utilities.mkdtemp()
  for analysis in analyses:
    for category in categories:
      if category == "Untagged": productionmode = "ggH"
      if category == "VBFtagged": productionmode = "VBF"
      if category == "VHHadrtagged": productionmode = "VH"
      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "projections", "fullrange", "{}_{}".format(analysis, production), str(category), "2e2mu", productionmode, "*.pdf"))
      assert theglob, (analysis, category)
      for filename in theglob:
         shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(analysis, category, os.path.basename(filename))))

  subprocess.check_call(["tar", "-cvzf", "templateprojections.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

if os.path.exists(os.path.join(config.plotsbasedir, "xchecks", "compare_POWHEG_JHUGen")):
  tmpdir = utilities.mkdtemp()
  for analysis in analyses:
    for productionmode in "VBF", "ZH", "WH", "HJJ", "ttH":
      theglob = glob.glob(os.path.join(config.plotsbasedir, "xchecks", "compare_POWHEG_JHUGen", str(analysis), productionmode, "*.pdf"))
      assert theglob, (analysis, productionmode)
      for filename in theglob:
         shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(analysis, productionmode, os.path.basename(filename))))

  subprocess.check_call(["tar", "-cvzf", "POWHEGvsJHUGen.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

if os.path.exists(os.path.join(config.plotsbasedir, "categorization")):
  tmpdir = utilities.mkdtemp()
  for a in analyses:
    for h in a.purehypotheses:
      h = str(h).replace("0+", "0plus").replace("0-", "0minus").replace("_photoncut", "")
      for p in "2jet", "HadWH", "HadZH":
        if a == "fL1Zg" and p == "HadWH" and h == "L1Zg": continue
        filename = os.path.join(config.plotsbasedir, "categorization", str(a), "D_{}_{}.pdf".format(p, h))
        shutil.copy(filename, os.path.join(tmpdir, "{}_{}".format(a, os.path.basename(filename))))

  subprocess.check_call(["tar", "-cvzf", "categorydiscriminants.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

if os.path.exists(os.path.join(config.plotsbasedir, "templateprojections", "niceplots")):
  tmpdir = utilities.mkdtemp()
  for a in analyses:
    for c in categories:
      if c == "Untagged": c = "Untagged_with2015"
      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", str(a), str(c), "*.pdf"))
      assert theglob, (a, c)
      for filename in theglob:
        if "D_bkg" in filename:
          filename = filename.replace("enrich", "fullrange")
        shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, c, os.path.basename(filename))))
      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "fullrange", str(a), "D_bkg_with2015.pdf"))
      assert theglob, (a, c)
      for filename in theglob:
        shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, c, os.path.basename(filename))))

  subprocess.check_call(["tar", "-cvzf", "discriminantdistributions.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

if glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", "*", "*", "animation")):
  tmpdir = utilities.mkdtemp()
  for a in analyses:
    for c in categories:
      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", str(a), str(c), "animation", "*.gif"))
      assert theglob, (a, c)
      for filename in theglob:
        if "D_bkg" in filename:
          filename = filename.replace("enrich", "fullrange")
        shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, c, os.path.basename(filename))))

  subprocess.check_call(["tar", "-cvzf", "discriminantdistributions_animations.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))
