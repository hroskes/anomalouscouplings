"""
This is here because I really do not want to include CJLST and all its dependencies as submodules
in CMSSW, and then do a huge scram b.
"""

import os
import urllib

SHA1 = "5178a02a67b631e4bc657e03924ea0e09dc08369"

def existsandishappy(filename):
    if not os.path.exists(filename):
        return False
    with open(filename) as f:
        for line in f:
            if "<!DOCTYPE html>" in line:  #doesn't exist on github, wget gives failure but urlretrieve doesn't
                return False
            break
    return True

def download(locationinZZAnalysis):
    filename = os.path.basename(locationinZZAnalysis)
    if existsandishappy(filename):
        return

    url = os.path.join("https://github.com/CJLST/ZZAnalysis/raw/", SHA1, locationinZZAnalysis)
    tmpfilename, message = urllib.urlretrieve(url)
    if not existsandishappy(tmpfilename):
        raise OSError("file wasn't downloaded! "+url)

    with open(tmpfilename) as tmpf, open(filename, "w") as f:
        for line in tmpf:
            if line.startswith("#include "):
                include = line.split("#include ")[1].strip()[1:-1]
                if "ZZAnalysis" in include:
                    if "interface" in include:
                        newinclude = include.split("/")[-1]
                    try:
                        line = line.replace(include, newinclude).replace("<", '"').replace(">", '"')
                        del newinclude
                    except UnboundLocalError:
                        raise OSError("Can't handle #including " + include)
            line = line.replace('extern "C" ', "")
            f.write(line)

if __name__ == "__main__":
    download("AnalysisStep/src/Category.cc")
    download("AnalysisStep/interface/Category.h")
