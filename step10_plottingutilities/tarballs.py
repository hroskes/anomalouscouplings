#!/usr/bin/env python

import glob
import os
import shutil
import subprocess
import tempfile

from helperstuff import config
from helperstuff.enums import analyses, categories

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

tmpdir = tempfile.mkdtemp()
for analysis in analyses:
    for category in categories:
        if category == "Untagged": productionmode = "ggH"
        if category == "VBFtagged": productionmode = "VBF"
        if category == "VHHadrtagged": productionmode = "VH"
        for filename in glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "fullrange", "{}_{}".format(analysis, production), str(category), "2e2mu", productionmode, "*.pdf")):
             shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(analysis, category, os.path.basename(filename))))

subprocess.check_call(["tar", "-cvzf", "templateprojections.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

tmpdir = tempfile.mkdtemp()
for analysis in analyses:
    for productionmode in "VBF", "ZH", "WH", "HJJ", "ttH":
        for filename in glob.glob(os.path.join(config.plotsbasedir, "templateprojections", "compare_POWHEG_JHUGen", str(analysis), productionmode, "*.pdf")):
             shutil.copy(filename, os.path.join(tmpdir, "{}_{}_{}".format(analysis, productionmode, os.path.basename(filename))))

subprocess.check_call(["tar", "-cvzf", "POWHEGvsJHUGen.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))

tmpdir = tempfile.mkdtemp()
for a in analyses:
    for h in a.purehypotheses:
        h = str(h).replace("0+", "0plus").replace("0-", "0minus").replace("_photoncut", "")
        for p in "2jet", "HadWH", "HadZH":
            if a == "fL1Zg" and p == "HadWH" and h == "L1Zg": continue
            filename = os.path.join(config.plotsbasedir, "categorization", str(a), "D_{}_{}.pdf".format(p, h))
            shutil.copy(filename, os.path.join(tmpdir, "{}_{}".format(a, os.path.basename(filename))))

subprocess.check_call(["tar", "-cvzf", "categorydiscriminants.tar.gz", "-C", tmpdir] + os.listdir(tmpdir))
