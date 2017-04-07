import getpass
import os
import re
import socket


if (".cern.ch" in socket.gethostname() or "lxplus" in socket.gethostname()) and getpass.getuser() == "hroskes":
    host = "lxplus"
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings_production/"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/"
    svndir = "/afs/cern.ch/work/h/hroskes/AN/notes"

elif ("login-node" in socket.gethostname() or "compute" in socket.gethostname() or "bigmem" in socket.gethostname()) and getpass.getuser() == "jroskes1@jhu.edu":
    host = "MARCC"
    repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/"

repositorydir = os.path.realpath(repositorydir)

try:
    repositorydir
    plotsbasedir
    host
except NameError:
    raise ValueError("Who/where are you?\n{}\n{}".format(socket.gethostname(), getpass.getuser()))

usedata = True
showblinddistributions = True
unblinddistributions = True
unblindscans = True
useQGTagging = False
useVHMETTagged = True

usefastpdf = True
usefastpdfdouble = False

applym4lshapesystematicsUntagged = True
applym4lshapesystematicsVBFVHtagged = True
applym4lshapesystematicsggH = True
applym4lshapesystematicsggHUntagged = True
applym4lshapesystematicsdiagonal = True #VBF in VBFtagged, VH in VHHadrtagged
combinem4lshapesystematics = True

applyZXshapesystematicsVBFVHtagged = True
applyZXshapesystematicsUntagged = True
mergeZXVBFVHsystematics = True

applyJECshapesystematics = False
applyMINLOsystematics = True

assert unblindscans <= unblinddistributions <= showblinddistributions <= usedata

expectedscanluminosity = 35.8671

lumi2015 = 2.7

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg() < 0.5

productionsforcombine = ["170222"]
if len(productionsforcombine) == 1:
    productionforcombine = productionsforcombine[0]

defaultnbins = 40
