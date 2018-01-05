from datetime import datetime, timedelta
import os
from pipes import quote
import re
import stat
import subprocess
import tempfile
import textwrap

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap

import config

if config.host == "lxplus":
    def submitjob(jobtext, jobname=None, jobtime=None, queue=None, interactive=False, waitids=[], waitsuccessids=[], outputfile=None, errorfile=None, docd=False, morerepmap=None, email=True):
        if not email: raise RuntimeError("Not sure how to turn off email for lxplus")
        if outputfile is not None:
            outputfile = outputfile.format(jobid="%J")
        if errorfile is not None:
            errorfile = errorfile.format(jobid="%J")
        if interactive: queue = "cmsinter"
        if queue is None:
            if jobtime is None:
                raise ValueError("Have to provide either queue or jobtime (or set interactive=True)!")
            #http://stackoverflow.com/a/4628148/5228524
            match = re.match(r'((?P<days>\d+?)-)?((?P<hours>\d+?):)?((?P<minutes>\d+?):)?((?P<seconds>\d+?))?', jobtime)
            if not match: raise ValueError("Bad time string {}".format(jobtime))
            kwargs = {key: int(value) for key, value in match.groupdict().iteritems() if value is not None}
            delta = timedelta(**kwargs)
            if not delta: raise ValueError("Job takes no time! {}".format(jobtime))
            if delta <= timedelta(minutes=8): queue = "8nm"
            elif delta <= timedelta(hours=1): queue = "1nh"
            elif delta <= timedelta(hours=8): queue = "8nh"
            elif delta <= timedelta(days=1): queue = "1nd"
            elif delta <= timedelta(days=2): queue = "2nd"
            elif delta <= timedelta(weeks=1): queue = "1nw"
            elif delta <= timedelta(weeks=2): queue = "2nw"
            else: raise ValueError("time {} is more than 2 weeks!".format(jobtime))
        run = """
                 echo .oO[jobtext]Oo. | bsub .oO[options]Oo.
              """
        optionsdict = {
                       "-q": queue,
                       "-J": jobname,
                       "-o": outputfile,
                       "-e": errorfile,
                       "-w": " && ".join(["ended({:d})".format(id) for id in waitids]+["done({:d})".format(id) for id in waitsuccessids])
                      }
        repmap = {
                  "jobtext": quote("cd .oO[CMSSW_BASE]Oo. && eval $(scram ru -sh) && .oO[docd]Oo. && "+jobtext),
                  "CMSSW_BASE": os.environ["CMSSW_BASE"],
                  "pwd": os.getcwd(),
                  "options": " ".join("{} {}".format(key, quote(value)) for key, value in optionsdict.iteritems() if value),
                  "docd": "cd .oO[pwd]Oo." if docd else "cd -"
                 }
        if morerepmap:
            assert not (set(repmap) & set(morerepmap))
            repmap.update(morerepmap)

        if interactive:
            repmap["options"] += " -I"
            subprocess.check_call(replaceByMap(run, repmap), shell=True)
        else:
            bsubout = subprocess.check_output(replaceByMap(run, repmap), shell=True)
            jobid = int(bsubout.split("<")[1].split(">")[0])
            print bsubout,
            return jobid

elif config.host == "MARCC":
    def submitjob(jobtext, jobname=None, jobtime=None, queue=None, interactive=False, waitids=[], outputfile=None, errorfile=None, docd=False, morerepmap=None, email=False, memory="3000M"):
        if queue is None:
            queue = "shared"
        if outputfile is not None:
            outputfile = outputfile.format(jobid="%j")
        if errorfile is not None:
            errorfile = errorfile.format(jobid="%j")
        jobtemplate = textwrap.dedent("""
            #!/bin/bash

            . /work-zfs/lhc/cms/cmsset_default.sh &&
            cd .oO[CMSSW_BASE]Oo.                 &&
            eval $(scram ru -sh)                  &&
            .oO[docd]Oo.                          &&
            echo "SLURM job running in: " `pwd`   &&

            #commands
            .oO[jobtext]Oo.
        """).strip()

        repmap = {
                  "jobtext": jobtext,
                  "CMSSW_BASE": os.environ["CMSSW_BASE"],
                  "pwd": os.getcwd(),
                  "jobname": jobname,
                  "jobtime": jobtime,
                  "queue": queue,
                  "waitids": ":".join("{:d}".format(id) for id in waitids),
                  "docd": "cd .oO[pwd]Oo." if docd else "cd $(mktemp -d)",
                 }
        if morerepmap:
            assert not (set(repmap) & set(morerepmap))
            repmap.update(morerepmap)

        options = {
                   "--job-name": ".oO[jobname]Oo.",
                   "--time": ".oO[jobtime]Oo.",
                   "--nodes": "1",
                   "--ntasks-per-node": "1",
                   "--partition": ".oO[queue]Oo.",
                   "--mem": memory,
                   "--output": outputfile,
                   "--error": errorfile,
                  }
        if waitids:
            options["--dependency"] = "afterany:.oO[waitids]Oo."
        if email:
            options["--mail-user"] = config.email
            options["--mail-type"] = "end"

        options = {replaceByMap(k, repmap): replaceByMap(v, repmap) for k, v in options.iteritems() if v is not None}
        options = ["{}={}".format(k, quote(v)) for k, v in options.iteritems()]

        if interactive:
            subprocess.check_call(
                                  ["srun"] + options +
                                  ["bash", "-c", replaceByMap(".oO[jobtext]Oo.", repmap)]
                                 )
        else:
            with tempfile.NamedTemporaryFile(bufsize=0) as f:
                f.write(replaceByMap(jobtemplate, repmap))
                os.chmod(f.name, os.stat(f.name).st_mode | stat.S_IEXEC)

                sbatchout = subprocess.check_output(["sbatch"] + options + [f.name])
                jobid = int(sbatchout.split()[-1])
                print sbatchout,
                return jobid

else:
    assert False
