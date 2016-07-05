import os
import numpy
import scipy

import matplotlib
matplotlib.use("Qt4Agg")
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('text', usetex=True)
from matplotlib import pyplot

from amuse.units import units
from amuse.units import constants
import amuse.plot as amuse_plot

from macro import *
# from cosmology import CosmologyCalculator
from ioparser import Toycluster2RuntimeOutputParser
from ioparser import Gadget2BinaryF77UnformattedType2Parser
# from cluster import ObservedCluster
from cluster import AnalyticalCluster
from cluster import SampledBox
import convert

from density_stability import plot_individual_cluster_density
from density_stability import plot_individual_cluster_mass
from density_stability import plot_individual_cluster_temperature


if __name__ == "__main__":
    """ Parse Toycluster output for two clusters living in a box.
    Obtain the radial density and radial mass profile for the DM and gas
    for both clusters and check that it matches the Chandra observation.
    """

    print 80*'-'
    print "Parsing Toycluster output: double cluster (Xm != 0)"
    print 80*'-'

    myRun="20160704T2243"
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

    plot_individual_cluster_density(world.halo0_numerical,
                                    world.halo0_analytical)
    pyplot.savefig(density0_filename)

    plot_individual_cluster_mass(world.halo0_numerical,
                                 world.halo0_analytical)
    pyplot.savefig(density1_filename)

    plot_individual_cluster_density(world.halo1_numerical,
                                    world.halo1_analytical)
    pyplot.savefig(mass0_filename)

    plot_individual_cluster_mass(world.halo1_numerical,
                                 world.halo1_analytical)
    pyplot.savefig(mass1_filename)

    plot_individual_cluster_temperature(world.halo1_numerical,
                                        world.halo1_analytical)
    pyplot.savefig(temperature0_filename)

    plot_individual_cluster_temperature(world.halo1_numerical,
                                        world.halo1_analytical)
    pyplot.savefig(temperature1_filename)

    print 80*'-'
