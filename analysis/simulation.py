"""
File: simulation.py
Author: Timo L. R. Halbesma <timohalbesma@gmail.com>
Date created: Tue Apr 19, 2016 04:22 pm
Last modified: Fri Apr 22, 2016 10:07 am


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


def make_plot(particles, outdir, i=0):
    pyplot.gca().set_aspect("equal")
    # pyplot.xlim(-150, 150)
    # pyplot.ylim(-150, 150)
    amuse_plot.plot(particles.x.as_quantity_in(units.kpc), particles.y.as_quantity_in(units.kpc), 'r.')

    pyplot.savefig(os.path.join(outdir, "integrated_dm{0:=04}.png".format(i)))
    pyplot.close()
    #pyplot.show()


def integrate(cluster, outdir):
    # TODO: check what converter actually does.
    # If I set the converter to the same units as Gadget2/Toycluster2, are the internal units the same?
    # NB Toycluster is cgs, the converter is SI
    converter = nbody_system.nbody_to_si(1.0e10 | units.MSun, 1 | units.kpc)

    dynamics = Gadget2(converter, number_of_workers=4, redirection='none')
    dynamics.parameters.epsilon_squared = 0.0000001 | nbody_system.length**2
    print dynamics.parameters
    dm = dynamics.dm_particles.add_particles(cluster.dm)

    make_plot(dynamics.particles, outdir)

    for i in range(1, 101):
        dynamics.evolve_model(i * (1 | units.Myr))
        print dynamics.model_time.as_quantity_in(units.Myr)
        make_plot(dynamics.particles, outdir, i)

    dynamics.stop()



if __name__ == "__main__":
    datadir="../../workdir/ToyclusterICs/20160420T1852/"
    logfile="runToycluster.log"
    icfile="IC_single_0"

    import os
    if not os.path.exists("out"):
        os.mkdir("out")
    outdir=datadir+"out"

    cluster = Cluster(datadir)
    integrate(cluster, outdir)
