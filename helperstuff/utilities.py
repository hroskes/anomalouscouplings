import collections
import contextlib
from itertools import tee, izip
import operator
import json
import os
import time

import ROOT

class KeyDefaultDict(collections.defaultdict):
    """
    http://stackoverflow.com/a/2912455
    """
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

class TFilesDict(KeyDefaultDict):
    def __init__(self):
        return super(TFilesDict, self).__init__(ROOT.TFile.Open)
    def __delitem__(self, key):
        self[key].Close()
        return super(TFilesDict, self).__delitem__(key)

tfiles = TFilesDict()

class MultiplyCounter(collections.Counter):
    def __add__(self, other):
        return type(self)(super(MultiplyCounter, self).__add__(other))
    def __sub__(self, other):
        return type(self)(super(MultiplyCounter, self).__sub__(other))
    def __mul__(self, other):
        return type(self)({k: v*other for k, v in self.iteritems()})
    def __rmul__(self, other):
        return type(self)({k: other*v for k, v in self.iteritems()})
    def __div__(self, other):
        return type(self)({k: v/other for k, v in self.iteritems()})

    def __imul__(self, other):
        for key in self:
            self[key] *= other
        return self
    def __idiv__(self, other):
        for key in self:
            self[key] /= other
        return self

def cache(function):
    cachename = "__cache_{}".format(function.__name__)
    def newfunction(self, *args):
        try:
            return getattr(self, cachename)[args]
        except AttributeError:
            setattr(self, cachename, {})
            return newfunction(self, *args)
        except KeyError:
            getattr(self, cachename)[args] = function(self, *args)
            return newfunction(self, *args)
    newfunction.__name__ = function.__name__
    return newfunction

@contextlib.contextmanager
def cd(newdir):
    """http://stackoverflow.com/a/24176022/5228524"""
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

class KeepWhileOpenFile(object):
    def __init__(self, name):
        self.filename = name
        self.pwd = os.getcwd()
        self.fd = self.f = None
    def __enter__(self):
        with cd(self.pwd):
            try:
                self.fd = os.open(self.filename, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except OSError, e:
                if e.errno == 17:
                    return None
                else:
                    raise
            else:
                self.f = os.fdopen(self.fd, 'w')
                return self.f
    def __exit__(self, *args):
        if self:
            self.f.close()
            with cd(self.pwd):
                os.remove(self.filename)
            self.fd = self.f = None
    def __nonzero__(self):
        return bool(self.f)

def jsonloads(jsonstring):
    try:
        return json.loads(jsonstring)
    except:
        print jsonstring
        raise

def getnesteddictvalue(thedict, *keys, **kwargs):
    hasdefault = False
    for kw, kwarg in kwargs.iteritems():
       if kw == "default":
           hasdefault = True
           default = kwarg
       else:
           raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    if len(keys) == 0:
        return thedict

    if hasdefault and keys[0] not in thedict:
        if len(keys) == 1:
            thedict[keys[0]] = default
        else:
            thedict[keys[0]] = {}

    return getnesteddictvalue(thedict[keys[0]], *keys[1:], **kwargs)

def pairwise(iterable):
    """
    https://docs.python.org/2/library/itertools.html#recipes
    """
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def callclassinitfunctions(*names):
    def callfunctions(cls):
        for name in names:
            getattr(cls, name)()
        return cls
    return callfunctions

def product(iterable):
    return reduce(operator.mul, iterable, 1)

def LoadMacro(filename):
    done = False
    while not done:
        with KeepWhileOpenFile(filename.rstrip("+")+".tmp") as kwof:
            if not kwof:
                print "Another process is already loading {}, waiting 5 seconds...".format(filename)
                time.sleep(5)
                continue
            error = ROOT.gROOT.LoadMacro(filename)
            if error:
                raise IOError("Couldn't load "+filename+"!")
            done = True
