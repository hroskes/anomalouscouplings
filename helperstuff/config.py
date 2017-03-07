import getpass
import socket

if (".cern.ch" in socket.gethostname() or "lxplus" in socket.gethostname()) and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"
elif "login-node" in socket.gethostname() and getpass.getuser() == "jroskes1@jhu.edu":
    repositorydir = "/work-zfs/lhc/heshy/ICHEPanomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/ICHEPanomalouscouplings/plots/"
else:
    raise ValueError("Who/where are you?")

unblinddistributions = True
unblindscans = True
expectedscanluminosity = 30
m4lmin, m4lmax = 105, 140
blindcut = lambda self: self.D_bkg_0plus() < 0.5
productionsforcombine = ("160729",)
