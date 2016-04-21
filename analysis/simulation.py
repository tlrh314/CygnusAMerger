"""
File: simulation.py
Author: Timo L. R. Halbesma <timohalbesma@gmail.com>
Date created: Tue Apr 19, 2016 04:22 pm
Last modified: Tue Apr 19, 2016 04:59 pm


"""

import numpy
# Bigger fontsize is better *_*
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})
# pyplot.rcParams.update({"text.usetex": True})

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import VectorQuantity
import amuse.plot as amuse_plot
from amuse.community.gadget2.interface import Gadget2

from initial import Cluster


def make_plot(particles):
    pyplot.gca().set_aspect("equal")
    # pyplot.xlim(-150, 150)
    # pyplot.ylim(-150, 150)
    amuse_plot.plot(particles.x.as_quantity_in(units.kpc), particles.y.as_quantity_in(units.kpc), 'r.')

    #pyplot.savefig(os.path.join("plots", "plot_galaxy_merger_{0}_{1:=04}.png".format("disk" if j else "total", i)))
    #pyplot.close()
    pyplot.show()


def integrate(cluster):
    converter = nbody_system.nbody_to_si(1.0e10 | units.MSun, 1 | units.kpc)
    print converter
    dynamics = Gadget2(converter, number_of_workers=4, redirection='none')
    dynamics.parameters.epsilon_squared = 0.0000001 | nbody_system.length**2
    print dynamics.parameters
    dm = dynamics.dm_particles.add_particles(cluster.dm)

    make_plot(dynamics.particles)

    for i in range(1, 101):
        break
        dynamics.evolve_model(i * (1 | units.yr))
        print dynamics.model_time.as_quantity_in(units.Myr)
        make_plots(dynamics.particles)

    dynamics.stop()


if __name__ == "__main__":
    cluster = Cluster("../input/IC_single_0_2e5")
    print "Parsed Cluster IC"
    integrate(cluster)
