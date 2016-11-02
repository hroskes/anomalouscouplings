import getpass
import socket


if (".cern.ch" in socket.gethostname() or "lxplus" in socket.gethostname()) and getpass.getuser() in ["hroskes","rbarr"] :
    host = "lxplus"
    if getpass.getuser() == "hroskes":
        repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings_production"
        plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/"
    if getpass.getuser() == "rbarr":
        repositorydir = "/afs/cern.ch/work/r/rbarr/anomalouscouplings"
        plotsbasedir = "/afs/cern.ch/user/r/rbarr/www/anomalouscouplings_production/"

elif "login-node" in socket.gethostname() and getpass.getuser() == "jroskes1@jhu.edu":
    host = "MARCC"
    repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/"

try:
    repositorydir
    plotsbasedir
    host
except NameError:
    raise ValueError("Who/where are you?\n{}\n{}".format(socket.gethostname(), getpass.getuser()))

usedata = False
unblinddistributions = False
unblindscans = False
applyshapesystematics = False
useQGTagging = False

assert unblindscans >= unblinddistributions >= usedata

expectedscanluminosity = 15.7

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg_0plus() < 0.5

productionsforcombine = ["160928"]

defaultnbins = 40
