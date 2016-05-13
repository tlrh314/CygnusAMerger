import glob
import os
import numpy
from matplotlib import pyplot

from amuse.lab import *

from parser import Gadget2BinaryF77UnformattedType2Parser
from parser import Toycluster2RuntimeOutputParser
from initial import Cluster
from analysis import plot_individual_cluster_mass
from analysis import plot_individual_cluster_density

IClog = Toycluster2RuntimeOutputParser("runToycluster.log")
# print IClog
pyplot.figure(figsize=(12,9))
cluster_at_time = Cluster("./", "./", "IC_single_0")
plot_individual_cluster_density(cluster_at_time)
pyplot.show()

import sys; sys.exit(0)

# assume order of writing is earlier -> later times
snaps = sorted(glob.glob("../runs/20160423T2219/snaps/snapshot_*"),  key=os.path.getmtime)
for snap in snaps:
    print "Checking snapshot:", snap
    # data = Gadget2BinaryF77UnformattedType2Parser(snap)
    # print data

    # path/to/snapshot_snapnr
    fname = snap.split('/')[-1]
    snapnr = fname.split('_')[-1]
    fname = "density_"+snapnr

    pyplot.figure(figsize=(12,9))
    cluster_at_time = Cluster("../runs/20160423T2219/ICs/", "", snap)
    plot_individual_cluster_mass(cluster_at_time)
    # pyplot.show()
    pyplot.savefig("../runs/20160423T2219/analysis/"+fname, dpi=500)


    print "Done checking snapshot:", snap
