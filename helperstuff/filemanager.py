import collections
import contextlib
import json
import os
import ROOT

class keydefaultdict(collections.defaultdict):
    """
    http://stackoverflow.com/a/2912455
    """
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret
tfiles = keydefaultdict(ROOT.TFile.Open)

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
