import abc
import collections
import contextlib
import errno
from itertools import tee, izip
import operator
import json
import os
import shutil
import sys
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
    cache = {}
    def newfunction(*args, **kwargs):
        try:
            return cache[args, tuple(sorted(kwargs.iteritems()))]
        except KeyError:
            cache[args, tuple(sorted(kwargs.iteritems()))] = function(*args, **kwargs)
            return newfunction(*args, **kwargs)
    newfunction.__name__ = function.__name__
    return newfunction

def cache_instancemethod(function):
    """
    for when self doesn't support __hash__
    """
    cachename = "__cache_{}".format(function.__name__)
    def newfunction(self, *args, **kwargs):
        try:
            return getattr(self, cachename)[args, tuple(sorted(kwargs.iteritems()))]
        except AttributeError:
            setattr(self, cachename, {})
            return newfunction(self, *args, **kwargs)
        except KeyError:
            getattr(self, cachename)[args, tuple(sorted(kwargs.iteritems()))] = function(self, *args, **kwargs)
            return newfunction(self, *args, **kwargs)
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
    def __init__(self, name, message=None):
        self.filename = name
        self.message = message
        self.pwd = os.getcwd()
        self.fd = self.f = None
        self.bool = False

    def __enter__(self):
        with cd(self.pwd):
            try:
                self.fd = os.open(self.filename, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except OSError as e:
                return None
            else:
                self.f = os.fdopen(self.fd, 'w')
                try:
                    if self.message is not None:
                        self.f.write(self.message+"\n")
                except:
                    self.__exit__()
                    raise
                finally:
                    self.f.close()
                self.bool = True
                return True

    def __exit__(self, *args):
        try:
            with cd(self.pwd):
                os.remove(self.filename)
        except OSError:
            pass #ignore it
        self.fd = self.f = None
        self.bool = False

    def __nonzero__(self):
        return self.bool

class Tee(object):
    """http://stackoverflow.com/a/616686/5228524"""
    def __init__(self, *openargs, **openkwargs):
        self.openargs = openargs
        self.openkwargs = openkwargs
    def __enter__(self):
        self.file = open(*self.openargs, **self.openkwargs)
        self.stdout = sys.stdout
        sys.stdout = self
    def __exit__(self, *args):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

class OneAtATime(KeepWhileOpenFile):
    def __init__(self, name, delay, message=None, task="doing this"):
        super(OneAtATime, self).__init__(name)
        self.delay = delay
        if message is None:
            message = "Another process is already {task}!  Waiting {delay} seconds."
        message = message.format(delay=delay, task=task)
        self.message = message

    def __enter__(self):
        while True:
            result = super(OneAtATime, self).__enter__()
            if result:
                return result
            print self.message
            time.sleep(self.delay)

def jsonloads(jsonstring):
    try:
        return json.loads(jsonstring)
    except:
        print jsonstring
        raise

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

def tlvfromptetaphim(pt, eta, phi, m):
    result = ROOT.TLorentzVector()
    result.SetPtEtaPhiM(pt, eta, phi, m)
    return result

def sign(x):
    return cmp(x, 0)

def generatortolist(function):
    return generatortolist_condition(lambda x: True)(function)

def generatortolist_condition(condition):
    def generatortolist(function):
        def newfunction(*args, **kwargs):
            return [_ for _ in function(*args, **kwargs) if condition(_)]
        newfunction.__name__ = function.__name__
        return newfunction
    return generatortolist


def rreplace(s, old, new, occurrence):
    """http://stackoverflow.com/a/2556252/5228524"""
    li = s.rsplit(old, occurrence)
    return new.join(li)

def mkdir_p(path):
    """http://stackoverflow.com/a/600612/5228524"""
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def is_almost_integer(flt):
    if isinstance(flt, (int, long)) or flt.is_integer(): return True
    if float("{:.8g}".format(flt)).is_integer(): return True
    return False

class JsonDict(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def keys(self): pass

    @property
    def default(self):
        return JsonDict.__nodefault

    @abc.abstractmethod
    def dictfile(self):
        """should be a member, not a method"""





    __nodefault = object()
    __dictscache = collections.defaultdict(lambda: None)

    def setvalue(self, value):
        self.setnesteddictvalue(self.getdict(), *self.keys, value=value)
        assert self.value == value

    def getvalue(self):
        try:
            return self.getnesteddictvalue(self.getdict(), *self.keys, default=self.default)
        except:
            print "Error while getting value of\n{!r}".format(self)
            raise

    @property
    def value(self):
        return self.getvalue()

    @value.setter
    def value(self, value):
        self.setvalue(value)

    @classmethod
    def getdict(cls, trycache=True):
      import globals
      if cls.__dictscache[cls] is None or not trycache:
        try:
          with open(cls.dictfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(cls.dictfile))
          except OSError:
            pass
          with open(cls.dictfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        cls.__dictscache[cls] = json.loads(jsonstring)
      return cls.__dictscache[cls]

    @classmethod
    def writedict(cls):
      dct = cls.getdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(cls.dictfile, "w") as f:
        f.write(jsonstring)

    @classmethod
    def getnesteddictvalue(cls, thedict, *keys, **kwargs):
        hasdefault = False
        for kw, kwarg in kwargs.iteritems():
           if kw == "default":
               if kwarg is not JsonDict.__nodefault:
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

        return cls.getnesteddictvalue(thedict[keys[0]], *keys[1:], **kwargs)

    @classmethod
    def setnesteddictvalue(cls, thedict, *keys, **kwargs):
        for kw, kwarg in kwargs.iteritems():
            if kw == "value":
                value = kwarg
            else:
                raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

        try:
            value
        except NameError:
            raise TypeError("Didn't provide value kwarg!")

        if len(keys) == 1:
            thedict[keys[0]] = value
            return

        if keys[0] not in thedict:
            thedict[keys[0]] = {}

        return cls.setnesteddictvalue(thedict[keys[0]], *keys[1:], **kwargs)

def LSB_JOBID():
    return os.environ.get("LSB_JOBID", None)

class LSF_creating(object):
    def __init__(self, *files, **kwargs):
        self.files = files
        for filename in files:
            if not filename.startswith("/"):
                raise ValueError("{} should be an absolute path!".format(filename))
        self.jsonfile = None
        for kw, kwarg in kwargs.iteritems():
            if kw == "jsonfile":
                self.jsonfile = kwarg
                if not self.jsonfile.startswith("/"): raise ValueError("jsonfile={} should be an absolute path!".format(self.jsonfile))
            else:
                raise TypeError("Unknown kwarg {}={}!".format(kw, kwarg))

    def __enter__(self):
        if not LSB_JOBID(): return self
        if self.jsonfile is not None:
            shutil.copy(self.jsonfile, "./")
            self.jsonfile = os.path.basename(self.jsonfile)
            if len(self.files) != 1: raise ValueError("only know how to handle 1 file")

            with open(os.path.basename(self.jsonfile)) as f:
                content = f.read()

            if content.count(self.files[0]) != 1:
                raise ValueError("{} is not in {}".format(self.files[0], self.jsonfile))
            content = content.replace(self.files[0], os.path.basename(self.files[0]))

            with open(os.path.basename(self.jsonfile), "w") as f:
                f.write(content)

        return self

    def __exit__(self, *errorinfo):
        if not LSB_JOBID(): return

        notcreated = []

        for filename in self.files:
            if os.path.exists(os.path.basename(filename)):
                shutil.copy(os.path.basename(filename), filename)
            else:
                notcreated.append(os.path.basename(filename))

        if notcreated:
            raise RuntimeError("\n".join("{} was not created!".format(os.path.basename(filename)) for filename in filenames))
