import os
import glob
import argparse

import numpy
import matplotlib
from matplotlib import pyplot
matplotlib.use("Agg", warn=False)

from deco import concurrent, synchronized
import amuse.plot as amuse_plot
from cluster import NumericalCluster
from ioparser import parse_gadget_parms


@concurrent(processes=8)
def plot_hist_gasvx(rundir, snapnr, TimeBetSnapshot, ic=False):
    simulation = NumericalCluster(
        icdir=rundir+"ICs/",
        snapdir=rundir+"ICs/" if ic else rundir+"snaps/",
        logfile="runToycluster.log",
        icfile="IC_single_0" if ic else "snapshot_"+snapnr,)
    simulation.gas, simulation.dm = simulation.place_ic_data_in_datamodel(simulation.raw_data)

    pyplot.figure(figsize=(12,9))
    amuse_plot.hist(simulation.gas.vx, normed=True, stacked=True,
                    bins=int(numpy.sqrt(simulation.raw_data.Ngas)))
    pyplot.xlabel("gas velocity x-direction")
    pyplot.ylabel("normalized count")
    pyplot.xlim(-2000, 1200)
    pyplot.suptitle("T = {0:04.2f} Gyr".format(TimeBetSnapshot*int(snapnr)),
                color="black", size=18, y=0.95)
    pyplot.savefig(rundir+"out/gasvx_"+snapnr)

@synchronized
def make_all_hists(timestamp):
    rundir="../runs/{0}/".format(timestamp)
    if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
        os.mkdir(rundir+"out")

    gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
    TimeBetSnapshot = gadgetparms["TimeBetSnapshot"]

    snaps = sorted(glob.glob(rundir+"snaps/snapshot_*"),
                   key=os.path.getmtime)
    for snap in snaps:
        # e.g. "../runs/20160727T1112/snaps/snapshot_051"
        snapnr = snap.split('/')[-1].split('_')[-1]  # gets the 051
        plot_hist_gasvx(rundir, snapnr, TimeBetSnapshot)

def new_argument_parser():
    description="Plot histogram of gas vx for all simulation snapshots"
    parser = argparse.ArgumentParser(description=description)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="simulationID", nargs=1,
        help="string of the Simulation ID")

    return parser


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    simulationID = arguments.simulationID[0]
    make_all_hists(simulationID)
