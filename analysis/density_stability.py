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

run = "/Volumes/Taurus/Toycluster/"
IClog = Toycluster2RuntimeOutputParser(run+"runToycluster.log")
print IClog

cluster_at_time = Cluster(run, run, "IC_double_0")
# cluster_at_time = Cluster("./", "./", "IC_single_0")
plot_individual_cluster_density(cluster_at_time)
pyplot.show()

import sys; sys.exit(0)

# assume order of writing is earlier -> later times
snaps = sorted(glob.glob(run+"/snaps/snapshot_*"),  key=os.path.getmtime)
for i, snap in enumerate(snaps):
    print "Checking snapshot:", snap

    # path/to/snapshot_snapnr
    fname = snap.split('/')[-1]
    snapnr = fname.split('_')[-1]
    fname = "density_"+snapnr

    pyplot.figure(figsize=(24,18))
    cluster_at_time = Cluster(run+"/ICs/", "", snap)
    plot_individual_cluster_density(cluster_at_time)
    # pyplot.show()
    pyplot.savefig(run+"/analysis/"+fname, dpi=500)

    print "Done checking snapshot:", snap
