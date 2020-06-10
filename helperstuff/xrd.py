import os
from pipes import quote
import subprocess
from utilities import cache, LSB_JOBID

xrdlstemplate = """
echo '
connect {host}
cd {folder}
ls
' | xrd
"""

@cache
def xrdls(host, folder):
    folder = os.path.join(folder, "")
    output = subprocess.check_output(xrdlstemplate.format(host=quote(host), folder=quote(folder)), shell=True)
    if "no such file or directory" in output:
        raise IOError("root://{}/{} doesn't exist!".format(host, folder))
    result = []
    for line in output.split("\n"):
        line = line.split(folder+"> ")[-1]
        if folder not in line:
            continue
        line = line.split(folder)[-1]
        result.append(line)
    return result
    

def exists(filename):
    if filename.startswith("/eos/cms/") and LSB_JOBID():
        filename = "root://eoscms/"+filename
    if filename.startswith("root://"):
        filename = filename.rstrip("/")
        filename = filename.replace("root://", "")
        host, location = filename.split("//")
        location = "/"+location
        dirname, basename = os.path.split(location)
        dirdirname, basedirname = os.path.split(dirname)
        try:
            return basename in xrdls(host, dirname)
        except IOError:
            return False
    else:
        try:
            with open(filename): return True
        except IOError:
            return False

if __name__ == '__main__':
    print "\n".join(xrdls["lxcms03", "/data3/Higgs/160624/"])
