from Alignment.OfflineValidation.TkAlAllInOneTool.helperfunctions import replaceByMap
import config
import os
from pipes import quote
import subprocess
import tempfile

if config.host == "lxplus":
    def submitjob(jobtext, jobname=None, jobtime=None, queue="1nd", interactive=False):
        run = """
                 echo .oO[jobtext]Oo. | bsub .oO[options]Oo.
              """
        if interactive and queue != "cmsinter":
            raise ValueError("Interactive jobs have to be submitted to cmsinter")
        optionsdict = {
                       "-q": queue,
                       "-J": jobname,
                      }
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
    def submitjob(jobtext, jobname=None, jobtime=None, queue="shared", interactive=False):
        jobtemplate = textwrap.dedent("""
            #!/bin/bash
            #SBATCH --job-name=.oO[jobname]Oo.
            #SBATCH --time=.oO[jobtime]Oo.
            #SBATCH --nodes=1
            #SBATCH --ntasks-per-node=1
            #SBATCH --partition=.oO[queue]Oo.
            #SBATCH --mem=3000

            . /work-zfs/lhc/cms/cmsset_default.sh &&
            cd .oO[CMSSW_BASE]Oo.                 &&
            eval $(scram ru -sh)                  &&
            cd .oO[pwd]Oo.                        &&
            echo "SLURM job running in: " `pwd`   &&

            #commands
            .oO[jobtext]Oo.
        """).strip()

        repmap = {
                  "jobtext": quote("cd .oO[CMSSW_BASE]Oo. && eval $(scram ru -sh) && cd .oO[pwd]Oo. && "+jobtext),
                  "CMSSW_BASE": os.environ["CMSSW_BASE"],
                  "pwd": os.getcwd(),
                  "jobname": jobname,
                  "jobtime": jobtime,
                  "queue": queue,
                 }

        with tempfile.NamedTemporaryFile(bufsize=0) as f:
            f.write(replaceByMap(jobtemplate, repmap))
            os.chmod(f.name, os.stat(f.name).st_mode | stat.S_IEXEC)

            if interactive:
                subprocess.check_call(["srun", f.name])
            else:
                sbatchout = subprocess.check_output(["sbatch", f.name])
                jobid = int(sbatchout.split()[-1])
                print sbatchout,
                return jobid

else:
    assert False
