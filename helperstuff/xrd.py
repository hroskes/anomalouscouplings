import os
from pipes import quote
import subprocess
from utilities import keydefaultdict

xrdlstemplate = """
echo '
connect {host}
cd {folder}
ls
' | xrd
"""

def xrdls(hostfolder):
    host, folder = hostfolder
    folder = os.path.join(folder, "")
    output = subprocess.check_output(xrdlstemplate.format(host=quote(host), folder=quote(folder)), shell=True)
    if "no such file or directory" in output:
        raise OSError("root://{}/{} doesn't exist!".format(host, folder))
    result = []
    for line in output.split("\n"):
        line = line.split(folder+"> ")[-1]
        if folder not in line:
            continue
        line = line.split(folder)[-1]
        result.append(line)
    return result
xrdls = keydefaultdict(xrdls)
    

def exists(filename):
    if filename.startswith("root://"):
        filename = filename.rstrip("/")
        filename = filename.replace("root://", "")
        host, location = filename.split("//")
        location = "/"+location
        dirname, basename = os.path.split(location)
        dirdirname, basedirname = os.path.split(dirname)
        if basedirname not in xrdls[host, dirdirname]: return False
        return basename in xrdls[host, dirname]
    else:
        return os.path.exists(filename)

if __name__ == '__main__':
    print "\n".join(xrdls["lxcms03", "/data3/Higgs/160624/"])
