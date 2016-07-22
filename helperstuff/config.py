import getpass
import socket

if "lxplus" in socket.gethostname() and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"
elif "login-node" in socket.gethostname() and getpass.getuser() == "jroskes1@jhu.edu":
    repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/"
else:
    raise ValueError("Who/where are you?")

unblinddistributions = True
unblindscans = True
expectedscanluminosity = 10
m4lmin, m4lmax = 105, 140
blindcut = lambda self: self.D_bkg_0plus() < 0.5
productionforcombine = "160720"
