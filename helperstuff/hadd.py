#!/usr/bin/env python

import os
import subprocess
import sys

def hadd(finalrootfile, *inputrootfiles):
    if os.path.exists(finalrootfile):
        raise OSError("{} exists!".format(finalrootfile))
    for _ in inputrootfiles:
        if not os.path.exists(_):
            raise OSError("{} does not exist!".format(_))
    try:
        subprocess.check_call(["hadd", finalrootfile] + list(inputrootfiles))
    except:
        try:
            raise
        finally:
            try:
                os.remove(finalrootfile)
            except:
                pass

if __name__ == "__main__":
    hadd(*sys.argv[1:])
