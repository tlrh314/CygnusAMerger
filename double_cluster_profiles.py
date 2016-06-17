import os
import numpy
import scipy

import matplotlib
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


if __name__ == "__main__":
    """ Parse Toycluster output for two clusters living in a box.
    Obtain the radial density and radial mass profile for the DM and gas
    for both clusters and check that it matches the Chandra observation.
    """

    print 80*'-'
    print "Parsing Toycluster output: double cluster (Xm != 0)"
    print 80*'-'

    timestamp="20160617T1544"
    # timestamp="20160617T1936"

    world = SampledBox(timestamp)

    world.halo0_numerical.get_gas_mass_via_density(DESNNGB=295)
    world.halo0_numerical.get_dm_mass_via_number_density()
    world.halo0_numerical.set_dm_density()

    world.halo1_numerical.get_gas_mass_via_density(DESNNGB=295)
    world.halo1_numerical.get_dm_mass_via_number_density()
    world.halo1_numerical.set_dm_density()

    plot_individual_cluster_density(world.halo0_numerical,
                                    world.halo0_analytical)
    pyplot.show()

    plot_individual_cluster_mass(world.halo0_numerical,
                                 world.halo0_analytical)
    pyplot.show()

    plot_individual_cluster_density(world.halo1_numerical,
                                    world.halo1_analytical)
    pyplot.show()

    plot_individual_cluster_mass(world.halo1_numerical,
                                 world.halo1_analytical)
    pyplot.show()

    print 80*'-'
