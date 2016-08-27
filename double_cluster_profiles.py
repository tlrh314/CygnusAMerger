import os
import numpy
import scipy
import argparse

import matplotlib
from matplotlib import pyplot
matplotlib.use("Agg", warn=False)
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('text', usetex=True)

from amuse.units import units
from amuse.units import constants
import amuse.plot as amuse_plot

from macro import *
# from cosmology import CosmologyCalculator
from ioparser import Toycluster2RuntimeOutputParser
from ioparser import Gadget2BinaryF77UnformattedType2Parser
from cluster import ObservedCluster
from cluster import AnalyticalCluster
from cluster import SampledBox
import convert

from density_stability import plot_individual_cluster_density
from density_stability import plot_individual_cluster_mass
from density_stability import plot_individual_cluster_temperature

def new_argument_parser():
    parser = argparse.ArgumentParser(description="""Plot rho(r), M(<r),
        T(r) for given SimulationID with two clusters.""")
    parser.add_argument("-s", "--show", dest="save",
        action="store_false", default=True,
        help="show instead of save plots; default is save")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="simulationID", nargs=1,
        help="string of the Simulation ID")

    return parser


if __name__ == "__main__":
    """ Parse Toycluster output for two clusters living in a box.
    Obtain the radial density and radial mass profile for the DM and gas
    for both clusters and check that it matches the Chandra observation.
    """

    arguments = new_argument_parser().parse_args()
    myRun = arguments.simulationID[0]
    print myRun

    print 80*'-'
    print "Parsing Toycluster output: double cluster (Xm != 0)"
    print 80*'-'

    snapnr = "000-ICOnly"

    outdir="../runs/{0}/out/".format(myRun)
    if not (os.path.isdir(outdir) or os.path.exists(outdir)):
        os.mkdir(outdir)

    world = SampledBox(myRun)

    world.halo0_numerical.get_gas_mass_via_density(DESNNGB=50)
    world.halo0_numerical.get_dm_mass_via_number_density()
    world.halo0_numerical.set_dm_density()

    world.halo1_numerical.get_gas_mass_via_density(DESNNGB=50)
    world.halo1_numerical.get_dm_mass_via_number_density()
    world.halo1_numerical.set_dm_density()

    mass0_filename = outdir+"{0}-mass-halo0-{1}.png".format(myRun, snapnr)
    density0_filename = outdir+"{0}-density-hal0-{1}.png".format(myRun, snapnr)
    temperature0_filename = outdir+"{0}-temperature-halo0-{1}.png".format(myRun, snapnr)
    mass1_filename = outdir+"{0}-mass-halo1-{1}.png".format(myRun, snapnr)
    density1_filename = outdir+"{0}-density-halo1-{1}.png".format(myRun, snapnr)
    temperature1_filename = outdir+"{0}-temperature-halo1-{1}.png".format(myRun, snapnr)

    observed_cygA = ObservedCluster("cygA")
    plot_individual_cluster_density(world.halo0_numerical,
        world.halo0_analytical, observed_cygA)

    if arguments.save:
        pyplot.savefig(density0_filename)
        pyplot.close()

    observed_cygB = ObservedCluster("cygB")
    plot_individual_cluster_density(world.halo1_numerical,
        world.halo1_analytical, observed_cygB)
    if arguments.save:
        pyplot.savefig(density1_filename)
        pyplot.close()

    plot_individual_cluster_mass(world.halo0_numerical,
                                 world.halo0_analytical)
    if arguments.save:
        pyplot.savefig(mass0_filename)
        pyplot.close()

    plot_individual_cluster_mass(world.halo1_numerical,
                                 world.halo1_analytical)
    if arguments.save:
        pyplot.savefig(mass1_filename)
        pyplot.close()

    plot_individual_cluster_temperature(world.halo0_numerical,
                                        world.halo0_analytical)
    if arguments.save:
        pyplot.savefig(temperature0_filename)
        pyplot.close()

    plot_individual_cluster_temperature(world.halo1_numerical,
                                        world.halo1_analytical)
    if arguments.save:
        pyplot.savefig(temperature1_filename)
        pyplot.close()
    else:
        pyplot.show()

    print 80*'-'
