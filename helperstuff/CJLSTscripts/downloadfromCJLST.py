"""
This is here because I really do not want to include CJLST and all its dependencies as submodules
in CMSSW, and then do a huge scram b.
"""

from abc import ABCMeta, abstractmethod, abstractproperty
import os
import shutil
import urllib


class CJLSTScript_base(object):
    __metaclass__ = ABCMeta
    class WrongExtensionError(Exception): pass

    @abstractproperty
    def validextensions(self):
        """give file extensions supported by this class"""
        return []

    def __init__(self, locationinZZAnalysis):
        self.locationinZZAnalysis = locationinZZAnalysis
        self.filename = os.path.basename(locationinZZAnalysis)
        extension = os.path.splitext(self.filename)[1]
        if extension not in self.validextensions:
            raise self.WrongExtensionError("{} does not support {} which has extension {}".format(type(self).__name__, locationinZZAnalysis, extension))

    def download(self, SHA1, force=False):
        if self.existsandishappy() and not force:
            return

        url = os.path.join("https://github.com/CJLST/ZZAnalysis/raw/", SHA1, self.locationinZZAnalysis)
        tmpfilename, message = urllib.urlretrieve(url)

        self.fixandmove(tmpfilename)

        if not self.existsandishappy():
            raise OSError("file wasn't downloaded! "+url)

    @abstractmethod
    def fixandmove(self, tmpfilename):
        shutil.move(tmpfilename, self.filename)

    def existsandishappy(self):
        try:
            with open(self.filename) as f:
                for line in f:
                    if "<!DOCTYPE html>" in line:  #doesn't exist on github, wget gives failure but urlretrieve doesn't
                        return False
                    break
        except IOError:
            return False
        return True

    def __str__(self):
        return self.locationinZZAnalysis

class CJLSTScript_Cpp(CJLSTScript_base):
    @property
    def validextensions(self):
        return [".cc", ".C", ".h"]
    def fixandmove(self, tmpfilename):
        with open(tmpfilename) as tmpf, open(self.filename, "w") as f:
            for line in tmpf:
                if line.startswith("#include "):
                    include = line.split("#include ")[1].strip()[1:-1]
                    if "ZZAnalysis" in include or ".." in include:
                        if include.split("/")[-2] in ("interface", "src"):
                            newinclude = include.split("/")[-1]
                            newinclude = newinclude.replace(".cc", ".h")
                        try:
                            line = line.replace(include, newinclude).replace("<", '"').replace(">", '"')
                            del newinclude
                        except UnboundLocalError:
                            raise OSError("Can't handle #including " + include)
                line = line.replace('extern "C" ', "")
                f.write(line)

class CJLSTScript_other(CJLSTScript_base):
    @property
    def validextensions(self):
        return [".root"]
    def fixandmove(self, tmpfilename):
        super(CJLSTScript_other, self).fixandmove(tmpfilename)

def CJLSTScript(*args, **kwargs):
    result = []
    exceptions = []
    for subclass in CJLSTScript_base.__subclasses__():
        try:
            result.append(subclass(*args, **kwargs))
        except CJLSTScript_base.WrongExtensionError as e:
            exceptions.append(e)
    if not result:
        raise CJLSTScript_base.WrongExtensionError("\n".join(str(e) for e in exceptions))
    if len(result) > 1:
        raise TypeError("Multiple classes ({}) allow the same extension!!".format(", ".join(str(type(_)) for _ in result)))
    return result[0]

class Downloader(object):
    downloadinfofilename = "download_info.txt"

    def __init__(self, SHA1, *thingstodownload):
        self.thingstodownload = []
        for _ in thingstodownload:
            self.add(_)
        self.SHA1 = SHA1

    def add(self, thingtodownload):
        self.thingstodownload.append(CJLSTScript(thingtodownload))

    def downloadinfocorrect(self):
        try:
            with open(self.downloadinfofilename) as f:
                contents = f.read()
        except IOError:
            return False
        if contents.strip() == self.downloadinfo.strip():
            return True

        #download info is wrong: remove old files
        filenames = contents.split("\n")[1:]
        for filename in filenames:
            if os.path.exists(os.path.basename(filename)):
                os.remove(os.path.basename(filename))
        return False

    def writedownloadinfo(self):
        with open(self.downloadinfofilename, "w") as f:
            f.write(self.downloadinfo)

    def download(self):
        for _ in self.thingstodownload:
            _.download(self.SHA1, force=not self.downloadinfocorrect())
        self.writedownloadinfo()

    @property
    def downloadinfo(self):
        return "\n".join([self.SHA1] + sorted(str(thing) for thing in self.thingstodownload))
