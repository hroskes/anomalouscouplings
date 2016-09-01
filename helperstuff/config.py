import getpass
import socket

if (".cern.ch" in socket.gethostname() or "lxplus" in socket.gethostname()) and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings_production"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/"
elif "login-node" in socket.gethostname() and getpass.getuser() == "jroskes1@jhu.edu":
    repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/"
else:
    raise ValueError("Who/where are you?")

unblinddistributions = False
unblindscans = False
applyshapesystematics = False
useQGTagging = False

expectedscanluminosity = 10

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg_0plus() < 0.5

productionsforcombine = ["160901"]
