#!/usr/bin/env python

import argparse
import glob
import os
import shutil
import subprocess

def tar_cvzf(tarballname, tmpdir, plotcopier=None):
  subprocess.check_call(["tar", "-cvzf", tarballname, "-C", tmpdir] + os.listdir(tmpdir))
  if plotcopier: plotcopier.copy(tarballname)

def templateprojections(plotcopier=None):
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

    tar_cvzf("templateprojections.tar.gz", tmpdir, plotcopier)

def POWHEGvsJHUGen(plotcopier=None):
  if os.path.exists(os.path.join(config.plotsbasedir, "xchecks", "compare_POWHEG_JHUGen")):
    tmpdir = utilities.mkdtemp()
    for analysis in analyses:
      for productionmode in "VBF", "ZH", "WH", "HJJ", "ttH":
        theglob = glob.glob(os.path.join(config.plotsbasedir, "xchecks", "compare_POWHEG_JHUGen", str(analysis), productionmode, "*.pdf"))
        assert theglob, (analysis, productionmode)
        for filename in theglob:
           shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(analysis, productionmode, os.path.basename(filename))))

    tar_cvzf("POWHEGvsJHUGen.tar.gz", tmpdir, plotcopier)

def categorydiscriminants(plotcopier=None):
  if os.path.exists(os.path.join(config.plotsbasedir, "categorization")):
    tmpdir = utilities.mkdtemp()
    for a in analyses:
      for h in a.purehypotheses:
        h = str(h).replace("0+", "0plus").replace("0-", "0minus")
        for p in "2jet", "HadWH", "HadZH":
          if a == "fL1Zg" and p == "HadWH" and h == "L1Zg": continue
          filename = os.path.join(config.plotsbasedir, "categorization", str(a), "D_{}_{}.pdf".format(p, h))
          shutil.copy(filename, os.path.join(tmpdir, "{}_{}".format(a, os.path.basename(filename))))

    tar_cvzf("categorydiscriminants.tar.gz", tmpdir, plotcopier)

def discriminantdistributions(plotcopier=None):
  if os.path.exists(os.path.join(config.plotsbasedir, "templateprojections", "niceplots")):
    tmpdir = utilities.mkdtemp()
    for a in analyses:
      for c in categories:
        cfolder = str(c)
#        if c == "Untagged": cfolder = "Untagged_with2015"
        theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", str(a), cfolder, "*.pdf"))
        assert theglob, (a, c)
        for filename in theglob:
          if "D_bkg" in filename:
            filename = filename.replace("enrich", "fullrange")
          shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, c, os.path.basename(filename).replace("_new", ""))))

#      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "fullrange", str(a), "D_bkg_with2015.pdf"))
#      assert theglob, a
#      for filename in theglob:
#        shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, "all", os.path.basename(filename.replace("_with2015", "")))))

    tar_cvzf("discriminantdistributions.tar.gz", tmpdir, plotcopier)

def animations(plotcopier=None):
  if glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", "*", "*", "animation")):
    tmpdir = utilities.mkdtemp()
    for a in analyses:
      for c in categories:
        cfolder = str(c)
        if c == "Untagged": cfolder = "Untagged_with2015"
        theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "enrich", str(a), cfolder, "animation", "*.gif"))
        assert theglob, (a, c)
        for filename in theglob:
          if "D_bkg" in filename:
            filename = filename.replace("enrich", "fullrange")
          shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, c, os.path.basename(filename))))

      theglob = glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "niceplots", "fullrange", str(a), "D_bkg_with2015.gif"))
      assert theglob, a
      for filename in theglob:
        shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(a, "all", os.path.basename(filename.replace("_with2015", "")))))

      theglob = glob.glob(os.path.join(config.plotsbasedir, "limits", "{}_limit_lumi35.8671.gif".format(a)))
      assert theglob, a
      for filename in theglob:
        shutil.copy(filename, os.path.join(tmpdir, "{}_limit.gif".format(a)))

    tar_cvzf("discriminantdistributions_animations.tar.gz", tmpdir, plotcopier)

if __name__ == "__main__":
  functions = {_.__name__: _ for _ in (templateprojections, POWHEGvsJHUGen, categorydiscriminants, discriminantdistributions, animations)}
  parser = argparse.ArgumentParser()
  parser.add_argument("tarball", choices=functions.values(), type=lambda _: functions[_], action="append")
  args = parser.parse_args()

from helperstuff import config, utilities
from helperstuff.enums import analyses, categories
from helperstuff.utilities import PlotCopier

if __name__ == "__main__":
  with PlotCopier(copyfromfolder=".", copytofolder="/afs/cern.ch/work/h/hroskes/AN/papers/HIG-18-002/trunk/Figures/stacked") as pc:
    for tarball in args.tarball:
      tarball(pc)
