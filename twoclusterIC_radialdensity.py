import glob
import os
import numpy
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 22})
# pyplot.rcParams.update({"text.usetex": True})

from amuse.units import units
from amuse.units.quantities import VectorQuantity
import amuse.plot as amuse_plot

import globals
from ioparser import Gadget2BinaryF77UnformattedType2Parser
from ioparser import Toycluster2RuntimeOutputParser
from cluster import NumericalCluster
from cluster import AnalyticalCluster


if __name__ == "__main__":
    # 20160526T1354 should have correct Mtot, Xm and rc_0, rc_1...
    numerical = NumericalCluster(
        icdir="../runs/20160526T1354/ICs/",
        snapdir="../runs/20160526T1354/snaps/",
        logfile="runToycluster.log",
        icfile="IC_single_0")

    # numerical.perform_sanity_checks()

    print numerical.raw_data.boxSize
    print numerical.toyclusterlog.kinematics['D_CoM_0']
    print numerical.toyclusterlog.kinematics['D_CoM_1']

    pyplot.figure(figsize=(12, 12))
    pyplot.title("Gas")
    amuse_plot.scatter(numerical.gas.x, numerical.gas.y, edgecolor="face", s=1)
    amuse_plot.xlabel("x")
    amuse_plot.ylabel("y")

    pyplot.figure(figsize=(12, 12))
    pyplot.title("Gas")
    amuse_plot.scatter(numerical.gas.z, numerical.gas.y, edgecolor="face", s=1)
    amuse_plot.xlabel("z")
    amuse_plot.ylabel("y")

    pyplot.figure(figsize=(12, 12))
    pyplot.title("Gas")
    amuse_plot.scatter(numerical.gas.x, numerical.gas.z, edgecolor="face", s=1)
    amuse_plot.xlabel("x")
    amuse_plot.ylabel("z")

    pyplot.figure(figsize=(12, 12))
    pyplot.title("DM")
    amuse_plot.scatter(numerical.dm.x, numerical.dm.y, edgecolor="face", s=1)
    amuse_plot.xlabel("x")
    amuse_plot.ylabel("y")

    pyplot.figure(figsize=(12, 12))
    pyplot.title("DM")
    amuse_plot.scatter(numerical.dm.z, numerical.dm.y, edgecolor="face", s=1)
    amuse_plot.xlabel("z")
    amuse_plot.ylabel("y")

    pyplot.figure(figsize=(12, 12))
    pyplot.title("DM")
    amuse_plot.scatter(numerical.dm.x, numerical.dm.z, edgecolor="face", s=1)
    amuse_plot.xlabel("x")
    amuse_plot.ylabel("z")

    pyplot.show()

    #pyplot.gca().set_yscale("log")
