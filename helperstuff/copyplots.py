#!/usr/bin/env python

import config
import os
import pipes
import subprocess
import sys

from utilities import cd, LSB_JOBID

def copyplots(folder):
  if LSB_JOBID():
    return
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
      command = ["rsync", "-azvP", "--relative", copyfrom, copyto]
      try:
        subprocess.check_call(command)
      except:
        print
        print "Failed to copy plots.  To do it yourself, try:"
        print
        print "("
        print "set -e"
        print "cd "+os.getcwd()
        print " ".join(pipes.quote(_) for _ in command)
        print ")"
        print
  else:
    assert False, config.host

if __name__ == "__main__":
  copyplots(sys.argv[1])
