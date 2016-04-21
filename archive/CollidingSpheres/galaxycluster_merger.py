"""
File: galaxycluster_merger.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Version: 0.01 (Initial)
Date created: Wed Mar 09, 2016 01:01 PM
Last modified:

Description: CygA Cluster Merger Simulation

"""

import os
import glob
import pickle
from datetime import datetime

import numpy

from PIL import Image
# https://raw.githubusercontent.com/rec/echomesh/master/code/python/external/images2gif.py
# import images2gif

from matplotlib import pyplot
pyplot.rcParams.update({'font.size': 22})
# from mpl_toolkits.mplot3d import Axes3D

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units.quantities import VectorQuantity
from amuse.ic.kingmodel import new_king_model
from amuse.ic.plummer import new_plummer_model
from amuse.community.gadget2.interface import Gadget2
from amuse.plot import plot
from amuse.plot import scatter
from amuse.plot import xlabel
from amuse.plot import ylabel
from amuse.plot import xlim
from amuse.plot import ylim
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.support.console import set_printing_strategy


class SubCluster(object):
    """ TODO: write docstring """

    def __init__(self, name, Mtot=(1e15 | units.MSun), Rvir=(500 | units.kpc),
                 Ndm=1e2, Ngas=1e3):
        self.name = name
        self.converter = nbody_system.nbody_to_si(Mtot, Rvir)
        self.Mtot = Mtot
        self.Rvir = Rvir
        self.Ndm = int(Ndm)
        self.Ngas = int(Ngas)

        # Set up numerical smoothing fractions.
        self.dm_smoothing_fraction = 0.001
        self.gas_smoothing_fraction = 0.05
        self.dm_epsilon = self.dm_smoothing_fraction * self.Rvir
        self.gas_epsilon = self.gas_smoothing_fraction * self.Rvir

        # Set up gas and dark matter fractions and gas/dm mass.
        self.gas_fraction = 0.1
        self.dm_fraction = 1.0 - self.gas_fraction
        self.dm_mass = self.dm_fraction * self.Mtot
        self.gas_mass = self.gas_fraction * self.Mtot

        # Set up King Sphere for the gas
        # Second parameter W0: Dimension-less depth of the King potential
        W0 = 3
        self.gas = new_king_model(self.Ngas, W0, convert_nbody=self.converter)
        self.gas.h_smooth = self.gas_epsilon
        self.gas.mass = (1.0/self.Ngas) * self.gas_mass
        self.gas.move_to_center()

        # Set up NFW profile for the Dark Matter (TODO: now plummer sphere!)
        self.dm = new_plummer_model(self.Ndm, convert_nbody=self.converter)
        self.dm.radius = self.dm_epsilon
        self.dm.mass = (1.0/self.Ndm) * self.dm_mass
        self.dm.move_to_center()

        print "Succesfuly created a Sub Cluster with the following properties"
        print self

        print self.gas

    def __str__(self):
        tmp = "Sub Cluster {0}\n".format(self.name)
        tmp += "\tMtot:\t\t{0}\n\tRvir\t\t{1}\n".format(self.Mtot, self.Rvir)
        tmp += "\tNgas:\t\t{0}\n\tNdm\t\t{1}\n".format(self.Ngas, self.Ndm)
        tmp += "\tGas fraction:\t{0}\n\tDM fraction:\t{1}".format(self.gas_fraction, self.dm_fraction)

        return tmp

    def dm_rvir_gas_sph_subplot(self):
        fig = pyplot.figure(figsize=(20, 10))
        ax_dm = fig.add_subplot(121, aspect='equal')
        ax_gas = fig.add_subplot(122, aspect='equal',
            sharex=ax_dm, sharey=ax_dm)

        # plot dark matter
        center_of_mass = self.dm.center_of_mass()
        virial_radius = self.dm.virial_radius().as_quantity_in(units.kpc)
        innersphere = self.dm.select(lambda r: (center_of_mass-r).length()<virial_radius,["position"])
        outersphere = self.dm.select(lambda r: (center_of_mass-r).length()>= virial_radius,["position"])
        pyplot.gcf().sca(ax_dm)
        x = outersphere.x.as_quantity_in(units.kpc)
        y = outersphere.y.as_quantity_in(units.kpc)
        scatter(x, y, c='red', edgecolor='red', label=r'$r \geq r_{\rm vir}$')
        x = innersphere.x.as_quantity_in(units.kpc)
        y = innersphere.y.as_quantity_in(units.kpc)
        scatter(x, y, c='green', edgecolor='green', label=r'$r < r_{\rm vir}$')
        xlabel(r'$x$')
        ylabel(r'$y$')
        pyplot.legend()

        # plot gas as sph plot

        # # Adjusted code from amuse.plot.sph_particles_plot
        # pyplot.gcf().sca(ax_gas)
        # min_size = 100
        # max_size = 10000
        # alpha = 0.1
        # x = self.gas.x.as_quantity_in(units.kpc)
        # y = self.gas.y.as_quantity_in(units.kpc)
        # z = self.gas.z.as_quantity_in(units.kpc)
        # z, x, y, us, h_smooths = z.sorted_with(x, y, self.gas.u, self.gas.h_smooth)
        # u_min, u_max = min(us), max(us)

        # log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
        # clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

        # red = 1.0 - clipped_log_u**4
        # blue = clipped_log_u**4
        # green = numpy.minimum(red, blue)

        # colors = numpy.transpose(numpy.array([red, green, blue]))
        # n_pixels = pyplot.gcf().get_dpi() * pyplot.gcf().get_size_inches()

        # ax_gas.set_axis_bgcolor('#101010')
        # ax_gas.set_aspect("equal", adjustable="datalim")
        # phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
        # sizes = numpy.minimum(numpy.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

        # scatter(x, y, color=colors, s=sizes, edgecolors="none", alpha=alpha)
        # xlabel(r'$x$')
        # ylabel(r'$y$')

        xlim(-2.*virial_radius, 2*virial_radius)
        ylim(-2.*virial_radius, 2*virial_radius)
        pyplot.tight_layout()
        pyplot.show()


class World(object):
    """ TODO: write docstring """
    def __init__(self, massratio=1./3):
        # Set up directories to store data in
        self.timestamp = datetime.today().strftime('%Y%m%dT%H%M')
        if not os.path.exists('out/{0}'.format(self.timestamp)):
            os.mkdir('out/{0}'.format(self.timestamp))
        if not os.path.exists('out/{0}/plots'.format(self.timestamp)):
            os.mkdir('out/{0}/plots'.format(self.timestamp))
        if not os.path.exists('out/{0}/data'.format(self.timestamp)):
            os.mkdir('out/{0}/data'.format(self.timestamp))

        # Set up sub clusters
        self.subClusterA = SubCluster(name="Sub Cluster A")
        self.subClusterB = SubCluster(name="Sub Cluster B",
            Mtot=massratio*(1e15 | units.MSun), Rvir=(200 | units.kpc))
        self.converter = self.subClusterA.converter

        self.timesteps = VectorQuantity.arange(0 | units.Myr, 1 | units.Gyr, 50 | units.Myr)

        # Set up world and gravity/hydro solvers
        self.place_clusters_in_world()

        # Write simulation parameters to text file
        filename = "out/{0}/data/merger.dat".format(self.timestamp)
        print "Dumping ClusterMerger instance to", filename, "\n"
        pickle.dump(self, open(filename, 'wb'))

        # Set up simulation codes
        self.setup_codes()

        print "Created subclusters.\n", str(self)

        hydro_plot(
            [-2.0, 2.0, -2.0, 2.0] | units.Mpc,
            self.code,
            (100, 100),
            base_output_file_name + "_hydro_image.png"
        )

    def __str__(self):
        tmp = "Initial Conditions"
        tmp += "\n" + str(self.subClusterA)
        tmp += "dm cm:\t\t\t" + str(self.subClusterA.dm.center_of_mass())
        tmp += "\ngas cm:\t\t\t" + str(self.subClusterA.gas.center_of_mass())
        tmp += "\ngas dm velocity:\t" + str(self.subClusterA.dm.center_of_mass_velocity())
        tmp += "\ngas cm velocity:\t" + str(self.subClusterA.gas.center_of_mass_velocity())

        tmp += "\n\n" + str(self.subClusterB)
        tmp += "dm cm:\t\t\t" + str(self.subClusterB.dm.center_of_mass())
        tmp += "\ngas cm:\t\t\t" + str(self.subClusterB.gas.center_of_mass())
        tmp += "\ndm cm velocity:\t\t" + str(self.subClusterB.dm.center_of_mass_velocity())
        tmp += "\ngas cm velocity\t\t" + str(self.subClusterB.gas.center_of_mass_velocity())

        return tmp

    def place_clusters_in_world(self):
        self.subClusterA.dm.rotate(0.0, numpy.pi/4, 0.0)
        self.subClusterA.dm.position += [0.0, 0.0, 0.0] | units.kpc
        self.subClusterA.dm.velocity -= [0.0, 0.0, 0.0] | units.km/units.s
        self.subClusterA.gas.rotate(0.0, numpy.pi/4, 0.0)
        self.subClusterA.gas.position += [0.0, 0.0, 0.0] | units.kpc
        self.subClusterA.gas.velocity -= [0.0, 0.0, 0.0] | units.km/units.s

        self.subClusterB.dm.rotate(numpy.pi/6, 0.0, 0.0)
        self.subClusterB.dm.position += [5.0, 0.0, 0.0] | units.Mpc
        self.subClusterB.dm.velocity -= [2000.0, 0.0, 0.0] | units.km/units.s
        self.subClusterB.gas.rotate(numpy.pi/6, 0.0, 0.0)
        self.subClusterB.gas.position += [5.0, 0.0, 0.0] | units.Mpc
        self.subClusterB.gas.velocity -= [2000.0, 0.0, 0.0] | units.km/units.s

    def setup_codes(self):
        # 1 code is created and started
        self.code = Gadget2(self.converter, redirection='none', number_of_workers=4)
        # 2 parameters are set
        self.code.parameters.epsilon_squared = 0.0000001 | nbody_system.length**2
        self.dmA = self.code.dm_particles.add_particles(self.subClusterA.dm)
        self.gasA = self.code.gas_particles.add_particles(self.subClusterA.gas)
        self.dmB = self.code.dm_particles.add_particles(self.subClusterB.dm)
        self.gasB = self.code.gas_particles.add_particles(self.subClusterB.gas)
        self.code.commit_particles()

        # self.Etot_init = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        # print "Ekin:", self.code.kinetic_energy
        # print "Epot:", self.code.potential_energy
        # print "Eth:", self.code.thermal_energy

    def evolve(self):
        tot = len(self.timesteps) - 1
        end_time = self.timesteps[-1]

        print "Starting Simulation :-)"
        print "Generating plots on the fly :-)"

        for i, time in enumerate(self.timesteps):
            print_progressbar(i, tot)
            self.code.evolve_model(time)
            self.dm_gas_sph_subplot(int(time.value_in(units.Myr)))

            write_set_to_file(self.gasA,
                "out/{0}/data/gasA_{1}.amuse"
                    .format(self.timestamp, int(time.value_in(units.Myr))),
                "amuse")
            write_set_to_file(self.dmA,
                "out/{0}/data/dmA_{1}.amuse"
                    .format(self.timestamp, int(time.value_in(units.Myr))),
                "amuse")
            write_set_to_file(self.gasB,
                "out/{0}/data/gasB_{1}.amuse"
                    .format(self.timestamp, int(time.value_in(units.Myr))),
                "amuse")
            write_set_to_file(self.dmB,
                "out/{0}/data/dmB_{1}.amuse"
                    .format(self.timestamp, int(time.value_in(units.Myr))),
                "amuse")
        print

        self.stop()

def hydro_plot(view, hydro_code, image_size, figname):
    """
    view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    hydro_code: hydrodynamics code in which the gas to be plotted is defined
    image_size: size of the output image in pixels (x, y)
    """
    if not HAS_MATPLOTLIB:
        return
    shape = (image_size[0], image_size[1], 1)
    size = image_size[0] * image_size[1]
    axis_lengths = [0.0, 0.0, 0.0] | units.m
    axis_lengths[0] = view[1] - view[0]
    axis_lengths[1] = view[3] - view[2]
    grid = Grid.create(shape, axis_lengths)
    grid.x += view[0]
    grid.y += view[2]
    speed = grid.z.reshape(size) * (0 | 1/units.s)
    rho, rhovx, rhovy, rhovz, rhoe = hydro_code.get_hydro_state_at_point(grid.x.reshape(size),
        grid.y.reshape(size), grid.z.reshape(size), speed, speed, speed)

    min_v =  800.0 | units.km / units.s
    max_v = 3000.0 | units.km / units.s
    min_rho = 3.0e-4 | units.g / units.cm**3
    max_rho = 0.1 | units.g / units.cm**3
    min_E = 1.0e11 | units.J / units.kg
    max_E = 1.0e13 | units.J / units.kg

    v_sqr = (rhovx**2 + rhovy**2 + rhovz**2) / rho**2
    E = rhoe / rho
    log_v = numpy.log((v_sqr / min_v**2)) / numpy.log((max_v**2 / min_v**2))
    log_rho = numpy.log((rho / min_rho)) / numpy.log((max_rho / min_rho))
    log_E = numpy.log((E / min_E)) / numpy.log((max_E / min_E))

    red   = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_rho)).reshape(shape)
    green = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_v)).reshape(shape)
    blue  = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_E)).reshape(shape)
    alpha = numpy.minimum(numpy.ones_like(log_v), numpy.maximum(numpy.zeros_like(log_v),
        numpy.log((rho / (10*min_rho))))).reshape(shape)

    rgba = numpy.concatenate((red, green, blue, alpha), axis = 2)

    pyplot.figure(figsize = (image_size[0]/100.0, image_size[1]/100.0), dpi = 100)
    im = pyplot.figimage(rgba, origin='lower')

    # pyplot.savefig(figname, transparent=True, dpi = 100)
    # print "\nHydroplot was saved to: ", figname
    pyplot.show()
    pyplot.close()

if __name__ == "__main__":
    set_printing_strategy("custom", preferred_units=[units.MSun, units.Mpc,
        units.Myr, units.kms, units.erg])
    print "hello, world"
    # subcluster = SubCluster("Sub Cluster A")
    merger = World()

