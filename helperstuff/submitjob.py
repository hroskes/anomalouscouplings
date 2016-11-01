from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap
import config
import os
from pipes import quote
import stat
import subprocess
import tempfile
import textwrap

if config.host == "lxplus":
    def submitjob(jobtext, jobname=None, jobtime=None, queue="1nd", interactive=False, waitids=[]):
        run = """
                 echo .oO[jobtext]Oo. | bsub .oO[options]Oo.
              """
        if interactive and queue != "cmsinter":
            raise ValueError("Interactive jobs have to be submitted to cmsinter")
        optionsdict = {
                       "-q": queue,
                       "-J": jobname,
                      }
        if waitids:
            optionsdict["-w"] = " && ".join("ended({:d})".format(id) for id in waitids)
        repmap = {
                  "jobtext": quote("cd .oO[CMSSW_BASE]Oo. && eval $(scram ru -sh) && cd .oO[pwd]Oo. && "+jobtext),
                  "CMSSW_BASE": os.environ["CMSSW_BASE"],
                  "pwd": os.getcwd(),
                  "options": " ".join("{} {}".format(key, quote(value)) for key, value in optionsdict.iteritems() if value is not None),
                 }
        if interactive:
            repmap["options"] += " -I"
            subprocess.check_call(replaceByMap(run, repmap), shell=True)
        else:
            bsubout = subprocess.check_output(replaceByMap(run, repmap), shell=True)
            jobid = int(bsubout.split("<")[1].split(">")[0])
            print bsubout,
            return jobid

elif config.host == "MARCC":
    def submitjob(jobtext, jobname=None, jobtime=None, queue="shared", interactive=False, waitids=[]):
        jobtemplate = textwrap.dedent("""
            #!/bin/bash

            . /work-zfs/lhc/cms/cmsset_default.sh &&
            cd .oO[CMSSW_BASE]Oo.                 &&
            eval $(scram ru -sh)                  &&
            cd .oO[pwd]Oo.                        &&
            echo "SLURM job running in: " `pwd`   &&

            #commands
            .oO[jobtext]Oo.
        """).strip()

        repmap = {
                  "jobtext": "cd .oO[CMSSW_BASE]Oo. && eval $(scram ru -sh) && cd .oO[pwd]Oo. && "+jobtext,
                  "CMSSW_BASE": os.environ["CMSSW_BASE"],
                  "pwd": os.getcwd(),
                  "jobname": jobname,
                  "jobtime": jobtime,
                  "queue": queue,
                  "waitids": ":".join("{:d}".format(id) for id in waitids)
                 }

        options = {
                   "--job-name": ".oO[jobname]Oo.",
                   "--time": ".oO[jobtime]Oo.",
                   "--nodes": "1",
                   "--ntasks-per-node": "1",
                   "--partition": ".oO[queue]Oo.",
                   "--mem": "3000",
                  }
        if waitids:
            options["--dependency"] = "afterany:.oO[waitids]Oo."

        options = {replaceByMap(k, repmap): replaceByMap(v, repmap) for k, v in options.iteritems()}
        options = [quote("{}={}".format(k, v)) for k, v in options.iteritems()]

        if interactive:
            subprocess.check_call(
                                  ["srun"] + options +
                                  ["bash", "-c", quote(replaceByMap(".oO[jobtext]Oo.", repmap))]
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
