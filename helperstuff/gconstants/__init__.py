#!/usr/bin/env python
import json
import os
import shutil
import urllib

from ..utilities import cd, LSB_JOBID, OneAtATime

def __init():
  from gconstantclass import GConstant, Process, gconstants

  def wget(url, dest=None):
    if dest is None: dest = os.path.basename(url)
    tmpfilename, message = urllib.urlretrieve(url)
    shutil.move(tmpfilename, dest)

  thingstodownload = []
  for gconstant in gconstants():
    thingstodownload.append([gconstant.url, gconstant.filename])
  thingstodownload.sort()

  with cd(os.path.dirname(__file__)), OneAtATime("download_info.txt.tmp", 5, task="downloading gconstants"):
    try:
      with open("download_info.txt") as f:
        if [GConstant.commit, thingstodownload] != json.load(f):
          raise IOError
      for url, dest in thingstodownload:
        if not os.path.exists(dest):
          raise IOError
    except (IOError, ValueError):
      for url, dest in thingstodownload:
        try:
          os.remove(dest)
        except OSError:
          pass
        wget(url, dest)
        if not os.path.exists(dest):
          raise OSError(dest+" was not downloaded successfully")
      with open("download_info.txt", "w") as f:
        json.dump([GConstant.commit, thingstodownload], f)

__init()

def gconstant(process, coupling, m4l):
  from gconstantclass import GConstant
  return GConstant(process, coupling).getvalue(m4l)
