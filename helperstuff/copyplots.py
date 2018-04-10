#!/usr/bin/env python

import config
import os
import subprocess
import sys

from utilities import cd

def copyplots(folder):
  if config.host == "lxplus":
    return
  elif config.host == "MARCC":
    with cd(config.plotsbasedir):
      folder = os.path.join(folder, "")
      if not os.path.exists(folder):
        raise ValueError("{} does not exist in {}".format(folder, config.plotsbasedir))
      lxplusconfig = config.getconfiguration("lxplus.cern.ch", config.lxplususername)
      copyto = "{}@lxplus.cern.ch:{}".format(config.lxplususername, lxplusconfig["plotsbasedir"])
      #https://stackoverflow.com/a/22908437/5228524
      copyfrom = folder
      subprocess.check_call(["rsync", "-azvP", "--relative", copyfrom, copyto])
  else:
    assert False, config.host

if __name__ == "__main__":
  copyplots(sys.argv[1])