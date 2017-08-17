import os
import subprocess

import ROOT

import config
from utilities import cache, tfiles

class CConstantMeta(type):
  __tocopy = []

  def __new__(cls, clsname, bases, dct):
    try:
      filenames = dct["filenames"]
    except KeyError:
      for base in bases:
        try:
          filenames = base.filenames
          break
        except AttributeError:
          pass
      else:
        raise TypeError("CConstant class {} does not define filenames!".format(clsname))

    cls.__tocopy += filenames.values()

    if config.host == "lxplus":
      pass
    elif config.host == "MARCC":
      filenames = dct["filenames"] = {k: os.path.join(cls.newcconstantsdir(), os.path.basename(v)) for k, v in filenames.iteritems()}

    return super(CConstantMeta, cls).__new__(cls, clsname, bases, dct)

  @classmethod
  def newcconstantsdir(cls):
    if config.host == "lxplus": pass
    elif config.host == "MARCC": return os.path.join(os.path.dirname(__file__), "cconstantsplines")

  def copyfiles(self):
    if config.host == "lxplus":
      pass
    elif config.host == "MARCC":
      newdir = self.newcconstantsdir()
      tocopy = set(self.__tocopy)
      if not all(os.path.exists(os.path.join(newdir, os.path.basename(filename))) for filename in self.__tocopy):
        subprocess.check_call(["rsync", "-azvP"] + ["{}@lxplus.cern.ch:{}".format(config.lxplususername, filename) for filename in tocopy] + [newdir])
      assert all(os.path.exists(os.path.join(newdir, os.path.basename(filename))) for filename in self.__tocopy)
    else:
      assert False

class CConstant(object):
  __metaclass__ = CConstantMeta
  __slots__ = []
  filenames = {}
  def __init__(self):
    type(self).copyfiles()

  def __hash__(self):
    return hash(type(self))


class DBkgKin(CConstant):
  filenames = {
    11*11*11*11: "/afs/cern.ch/work/u/usarica/public/forMELAv206AndAbove/SmoothKDConstant_m4l_Dbkgkin_4e13TeV.root",
    11*11*13*13: "/afs/cern.ch/work/u/usarica/public/forMELAv206AndAbove/SmoothKDConstant_m4l_Dbkgkin_2e2mu13TeV.root",
    13*13*13*13: "/afs/cern.ch/work/u/usarica/public/forMELAv206AndAbove/SmoothKDConstant_m4l_Dbkgkin_4mu13TeV.root",
  }
  for k, v in filenames.items():
    for factor in 2, 4, -1, -2, -4:
      filenames[k*factor] = filenames[k]
  del k, v, factor

  @cache
  def getspline(self, ZZFlav):
    filename = self.filenames[ZZFlav]
    f = tfiles[filename]
    return f.sp_gr_varTrue_Constant_Smooth

  def __call__(self, ZZFlav, ZZMass):
    return self.getspline(ZZFlav).Eval(ZZMass)

getDbkgkinConstant = DBkgKin()

def getDbkgConstant(*args, **kwargs): return getDbkgkinConstant(*args, **kwargs)

def getDVBF2jetsConstant(*args, **kwargs): return None
def getDVBF1jetConstant(*args, **kwargs): return None
def getDZHhConstant(*args, **kwargs): return None
def getDWHhConstant(*args, **kwargs): return None

def getDVBF2jetsConstant_shiftWP(*args, **kwargs): return None
def getDVBF1jetConstant_shiftWP(*args, **kwargs): return None
def getDZHhConstant_shiftWP(*args, **kwargs): return None
def getDWHhConstant_shiftWP(*args, **kwargs): return None

def getDVBF2jetsWP(*args, **kwargs): return None
def getDVBF1jetWP(*args, **kwargs): return None
def getDZHhWP(*args, **kwargs): return None
def getDWHhWP(*args, **kwargs): return None
