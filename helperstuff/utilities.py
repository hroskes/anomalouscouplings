import abc
import collections
import contextlib
import errno
from functools import wraps
import inspect
from itertools import tee, izip
import logging
import operator
import json
import os
import pipes
import shutil
import subprocess
import sys
import tempfile
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
    def clear(self):
        for key in self.keys(): del self[key]
        return super(TFilesDict, self).clear()

tfiles = TFilesDict()

class MultiplyCounter(collections.Counter):
    def __init__(self, *args, **kwargs):
        self.__frozen = False
        super(MultiplyCounter, self).__init__(*args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        super(MultiplyCounter, self).__setitem__(*args, **kwargs)

    def __add__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        return type(self)(super(MultiplyCounter, self).__add__(other))
    def __sub__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        return type(self)(super(MultiplyCounter, self).__sub__(other))
    def __mul__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        return type(self)({k: v*other for k, v in self.iteritems()})
    def __rmul__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        return type(self)({k: other*v for k, v in self.iteritems()})
    def __div__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        return type(self)({k: v/other for k, v in self.iteritems()})

    def __imul__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        for key in self:
            self[key] *= other
        return self
    def __idiv__(self, other):
        if self.__frozen: raise TypeError("MultiplyCounter is already frozen!")
        for key in self:
            self[key] /= other
        return self

    def freeze(self):
        self.__frozen = True

def cache(function):
    cache = {}
    @wraps(function)
    def newfunction(*args, **kwargs):
        try:
            return cache[args, tuple(sorted(kwargs.iteritems()))]
        except TypeError:
            print args, tuple(sorted(kwargs.iteritems()))
            raise
        except KeyError:
            cache[args, tuple(sorted(kwargs.iteritems()))] = function(*args, **kwargs)
            return newfunction(*args, **kwargs)
    return newfunction

def cache_instancemethod(function):
    """
    This one can't take arguments.
    But the cache clears when self is deleted (as opposed to the cache keeping self alive).
    Probably could be modified to take arguments without too much trouble.
    """
    @wraps(function)
    def newfunction(self):
        if not hasattr(self, "__cache_instancemethod_{}".format(function.__name__)):
            setattr(self, "__cache_instancemethod_{}".format(function.__name__), function(self))
        return getattr(self, "__cache_instancemethod_{}".format(function.__name__))

def multienumcache(function, haskwargs=False, multienumforkey=None):
    from enums import MultiEnum
    if multienumforkey is None:
        multienumforkey = function
    assert issubclass(function, MultiEnum)
    assert issubclass(multienumforkey, MultiEnum)
    cache = {}
    def newfunction(*args, **kwargs):
        if kwargs and not haskwargs:
            raise TypeError("{} has no kwargs!".format(function.__name__))
        key = multienumforkey(*args)
        try:
            oldkwargs, result = cache[key]
            if kwargs and kwargs != oldkwargs:
                raise ValueError("{}({}, **kwargs) called with 2 different kwargs:\n{}\n{}".format(function.__name__, ", ".join(repr(_) for _ in args), oldkwargs, kwargs))
            return result
        except KeyError:
            if haskwargs and not kwargs:
                raise ValueError("Have to give kwargs the first time you call {}({}, **kwargs)".format(function.__name__, ", ".join(repr(_) for _ in args)))
            cache[key] = kwargs, function(*args, **kwargs)
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
        logging.debug("creating KeepWhileOpenFile {}".format(name))
        self.filename = name
        self.__message = message
        self.pwd = os.getcwd()
        self.fd = self.f = None
        self.bool = False

    def __enter__(self):
        logging.debug("entering KeepWhileOpenFile {}".format(self.filename))
        with cd(self.pwd):
            logging.debug("does it exist? {}".format(os.path.exists(self.filename)))
            try:
                logging.debug("trying to open")
                self.fd = os.open(self.filename, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except OSError:
                logging.debug("failed: it already exists")
                return None
            else:
                logging.debug("succeeded: it didn't exist")
                logging.debug("does it now? {}".format(os.path.exists(self.filename)))
                if not os.path.exists(self.filename):
                    logging.warning("{} doesn't exist!??".format(self.filename))
                self.f = os.fdopen(self.fd, 'w')
                try:
                    if self.__message is not None:
                        logging.debug("writing message")
                        self.f.write(self.__message+"\n")
                        logging.debug("wrote message")
                except IOError:
                    logging.debug("failed to write message")
                    pass
                try:
                    logging.debug("trying to close")
                    self.f.close()
                    logging.debug("closed")
                except IOError:
                    logging.debug("failed to close")
                    pass
                self.bool = True
                return True

    def __exit__(self, *args):
        logging.debug("exiting")
        if self:
            try:
                with cd(self.pwd):
                    logging.debug("trying to remove")
                    os.remove(self.filename)
                    logging.debug("removed")
            except OSError:
                logging.debug("failed")
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
        self.__message = message

    def __enter__(self):
        while True:
            result = super(OneAtATime, self).__enter__()
            if result:
                return result
            print self.__message
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
        if isinstance(value, list): value = tuple(value)
        self.setnesteddictvalue(self.getdict(), *self.keys, value=value)
        assert self.value == value, (self.value, value)

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
    import config
    if config.host == "lxplus":
        return os.environ.get("LSB_JOBID", None)
    if config.host == "MARCC":
        return os.environ.get("SLURM_JOBID", None)
    assert False, config.host

class LSF_creating(object):
    def __init__(self, *files, **kwargs):
        self.files = files
        for filename in files:
            if not filename.startswith("/"):
                raise ValueError("{} should be an absolute path!".format(filename))

        self.jsonfile = None
        self.ignorefailure = False
        for kw, kwarg in kwargs.iteritems():
            if kw == "jsonfile":
                self.jsonfile = kwarg
                if not self.jsonfile.startswith("/"): raise ValueError("jsonfile={} should be an absolute path!".format(self.jsonfile))
            elif kw == "ignorefailure":
                self.ignorefailure = kwarg
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

    def basename(self, filename):
        if filename not in self.files: raise ValueError("Unknown filename {}!".format(filename))
        if LSB_JOBID():
            return os.path.basename(filename)
        else:
            return filename

    def __exit__(self, *errorinfo):
        if not LSB_JOBID(): return

        notcreated = []

        for filename in self.files:
            if os.path.exists(os.path.basename(filename)):
                shutil.move(os.path.basename(filename), filename)
            else:
                notcreated.append(os.path.basename(filename))

        if notcreated and not self.ignorefailure:
            raise RuntimeError("\n".join("{} was not created!".format(os.path.basename(filename)) for filename in filenames))

def RooArgList(*args, **kwargs):
    name = None
    for kw, kwarg in kwargs.iteritems():
        if kw == "name":
            name = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}!".format(kw, kwarg))
    args = list(args)
    if name is None and isinstance(args[-1], basestring):
        name = args[-1]
        args = args[:-1]

    if len(args) < 4:
        if name is not None:
            args.append(name)
        return ROOT.RooArgList(*args)

    nameargs = [name] if name is not None else []
    result = ROOT.RooArgList(*nameargs)
    for arg in args:
        result.add(arg)
    return result

def inscreen():
    return bool(os.environ.get("STY"))

class DummyContextManager(object):
    def __enter__(self): return self
    def __exit__(*stuff): pass

def mkdtemp(**kwargs):
    import config
    if "dir" not in kwargs:
        if LSB_JOBID() is not None:
            if config.host == "lxplus":
                kwargs["dir"] = os.environ["LSB_JOB_TMPDIR"]
            elif config.host == "MARCC":
                pass
            else:
                assert False, config.host
    return tempfile.mkdtemp(**kwargs)

def getmembernames(*args, **kwargs):
    return [_[0] for _ in inspect.getmembers(*args, **kwargs)]

lastcmsswbase = None

"""
this doesn't work
def cmsenv(folder="."):
    global lastcmsswbase
    #pythonpath needs special handling
    oldpythonpath = os.environ["PYTHONPATH"].split(":")
    indexinsyspath = sys.path.index(oldpythonpath[0])

    with cd(folder):
        scram = subprocess.check_output(["scram", "ru", "-sh"])
        for line in scram.splitlines():
            potentialerror = ValueError("Unknown scram b output:\n{}".format(line))
            if line.split()[0] == "unset":
                if line[-1] != ';': raise potentialerror
                line = line[:-1]
                for variable in line.split()[1:]:
                    del os.environ[variable]
            elif line.split()[0] == "export":
                afterexport = line.split(None, 1)[1]
                variable, value = afterexport.split("=", 1)
                if value[0] != '"' or value[-2:] != '";': raise potentialerror
                value = value[1:-2]
                if "\\" in value or '"' in value or "'" in line: raise potentialerror
                os.environ[variable] = value
            else:
                raise potentialerror

    newpythonpath = os.environ["PYTHONPATH"].split(":")
    for _ in oldpythonpath: sys.path.remove(_)
    sys.path[indexinsyspath:indexinsyspath] = newpythonpath

    if lastcmsswbase is not None and os.environ["CMSSW_BASE"] != lastcmsswbase:
        raise ValueError("Can't cmsenv in both {} and {}!".format(lastcmsswbase, os.environ["CMSSW_BASE"]))
    lastcmsswbase = os.environ["CMSSW_BASE"]
"""

def requirecmsenv(folder):
    needcmsswbase = subprocess.check_output("cd {} && eval $(scram ru -sh) >& /dev/null && echo $CMSSW_BASE".format(pipes.quote(folder)), shell=True).strip()
    cmsswbase = os.environ["CMSSW_BASE"]
    if cmsswbase != needcmsswbase:
        raise ValueError("Need to cmsenv in {}!".format(needcmsswbase))

def deletemelastuff():
    if os.path.exists("Pdfdata"):
        shutil.rmtree("Pdfdata")
    for thing in "br.sm1", "br.sm2", "ffwarn.dat", "input.DAT", "process.DAT":
        if os.path.exists(thing):
            os.remove(thing)

class cdtemp_slurm(object):
    def __enter__(self):
        import config
        self.cd = None
        if config.host == "lxplus":
            return
        elif config.host == "MARCC":
            if LSB_JOBID() is not None:
                self.cd = cd(mkdtemp())
                return self.cd.__enter__()
            else:
                return
        else:
            assert False, config.host

    def __exit__(self, *args, **kwargs):
        if self.cd is not None:
            return self.cd.__exit__(*args, **kwargs)

def recursivesubclasses(cls):
    result = [cls]
    for subcls in cls.__subclasses__():
        result += recursivesubclasses(subcls)
    return result
