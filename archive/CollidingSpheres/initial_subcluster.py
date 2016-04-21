"""
File: initial_subcluster.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Thu Feb 25, 2016 10:03 am
Last modified: Mon Mar 07, 2016 07:09 pm

Description: Initial conidtions for two subclusters in Cygnus A cluster.

"""

import numpy

from matplotlib import pyplot
pyplot.rcParams.update({'font.size': 22})
from mpl_toolkits.mplot3d import Axes3D

from amuse.units import units
from amuse.units import nbody_system
# TODO: change Plummer sphere to.. something else?
from amuse.ic.plummer import new_plummer_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.plot import plot
from amuse.plot import scatter
from amuse.plot import xlabel
from amuse.plot import ylabel
from amuse.plot import xlim
from amuse.plot import ylim

class SubCluster(object):
    """ Set up a virialized cluster of virial radius Rvir, and total mass Mtot,
        where gas_fraction of the mass is in the cluster gas and 1-gas_fraction
        of the mass is dark matter. Both the gas and dark matter are initially
        at rest and located at the center. For now we assume Plummer spheres. """
    def __init__(self, name, Mtot=(1e15 | units.MSun), Rvir=(500 | units.kpc)):
        self.name = name
        self.converter = nbody_system.nbody_to_si(Mtot, Rvir)
        self.number_of_dm_particles = 1e2
        self.number_of_gas_particles = 1e3

        # Set up numerical smoothing fractions
        self.dm_smoothing_fraction = 0.001
        self.gas_smoothing_fraction = 0.05
        self.dm_epsilon = self.dm_smoothing_fraction * Rvir
        self.gas_epsilon = self.gas_smoothing_fraction * Rvir

        # Setup up gas and dark matter fractions and gas/dm mass
        self.gas_fraction = 0.1
        self.dm_fraction = 1.0 - self.gas_fraction
        self.dm_mass = self.dm_fraction * Mtot
        self.gas_mass = self.gas_fraction * Mtot

        # Set up dark matter particles
        # TODO: probably not use Plummer sphere
        self.dm = new_plummer_model(self.number_of_dm_particles, convert_nbody=self.converter)
        self.dm.radius = self.dm_epsilon
        self.dm.mass = (1.0/self.number_of_dm_particles) * self.dm_mass
        self.dm.move_to_center()

        # Set up virialized gas sphere of N gas particles, Plummer sphere
        # TODO: probably not use Plummer sphere
        self.gas = new_plummer_gas_model(self.number_of_gas_particles, convert_nbody=self.converter)
        self.gas.h_smooth = self.gas_epsilon
        self.gas.move_to_center()
        self.gas.mass = (1.0/self.number_of_gas_particles) * self.gas_mass

        # print "Created ", str(self)

    def __str__(self):
        tmp = "{0}\n".format(self.name)
        tmp += "Mtot:\t\t\t" + str((self.dm.mass.sum() + self.gas.mass.sum()))
        tmp += " (gas fraction: {0}, dm fraction: {1})\n".format(self.gas_fraction, self.dm_fraction)
        tmp += "Mgas:\t\t\t" + str((self.gas.mass.sum()))
        tmp += "\nMdm:\t\t\t" + str((self.dm.mass.sum()))
        tmp += "\nRvir dm:\t\t" + str(self.dm.virial_radius())
        tmp += "\nRvir gas:\t\t" + str(self.gas.virial_radius())
        tmp += "\n"

        return tmp

    def dm_rvir_gas_sph_3dsubplot(self):
        fig = pyplot.figure(figsize=(20, 10))
        ax_dm = fig.add_subplot(121, aspect='equal', projection='3d')
        ax_gas = fig.add_subplot(122, aspect='equal', projection='3d',
            sharex=ax_dm, sharey=ax_dm)

        # plot dark matter
        center_of_mass = self.dm.center_of_mass()
        virial_radius = self.dm.virial_radius().as_quantity_in(units.kpc)
        innersphere = self.dm.select(lambda r: (center_of_mass-r).length()<virial_radius,["position"])
        outersphere = self.dm.select(lambda r: (center_of_mass-r).length()>= virial_radius,["position"])
        pyplot.gcf().sca(ax_dm)
        x = outersphere.x.as_quantity_in(units.kpc)
        y = outersphere.y.as_quantity_in(units.kpc)
        z = outersphere.z.as_quantity_in(units.kpc)
        plot(x, y, z, 'o', c='red', label=r'$r \geq r_{\rm vir}$')
        x = innersphere.x.as_quantity_in(units.kpc)
        y = innersphere.y.as_quantity_in(units.kpc)
        z = innersphere.z.as_quantity_in(units.kpc)
        plot(x, y, z, 'o', c='green', label=r'$r < r_{\rm vir}$')
        xlabel(r'$x$')
        ylabel(r'$y$')
        ax_dm.set_zlabel(r'$z$ [{0}]'.format(virial_radius.unit))
        pyplot.legend()

        # plot gas as sph plot
        # Adjusted code from amuse.plot.sph_particles_plot
        pyplot.gcf().sca(ax_gas)
        min_size = 100
        max_size = 10000
        alpha = 0.1
        x = self.gas.x.as_quantity_in(units.kpc)
        y = self.gas.y.as_quantity_in(units.kpc)
        z = self.gas.z.as_quantity_in(units.kpc)
        z, x, y, us, h_smooths = z.sorted_with(x, y, self.gas.u, self.gas.h_smooth)
        u_min, u_max = min(us), max(us)

        log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
        clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

        red = 1.0 - clipped_log_u**4
        blue = clipped_log_u**4
        green = numpy.minimum(red, blue)

        colors = numpy.transpose(numpy.array([red, green, blue]))
        n_pixels = pyplot.gcf().get_dpi() * pyplot.gcf().get_size_inches()

        ax_gas.set_axis_bgcolor('#101010')
        ax_gas.set_aspect("equal", adjustable = "datalim")
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
        sizes = numpy.minimum(numpy.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

        ax_gas.scatter(x.number, y.number, z.number, color=colors, s=sizes, edgecolors="none", alpha=alpha)
        xlabel(r'$x$')
        ylabel(r'$y$')
        ax_gas.set_zlabel(r'$z$ [{0}]'.format(virial_radius.unit))

        xlim(-2.*virial_radius, 2*virial_radius)
        ylim(-2.*virial_radius, 2*virial_radius)
        ax_dm.set_zlim(-2.*virial_radius.number, 2*virial_radius.number)
        ax_gas.set_zlim(-2.*virial_radius.number, 2*virial_radius.number)
        pyplot.tight_layout()
        pyplot.show()

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

        # Adjusted code from amuse.plot.sph_particles_plot
        pyplot.gcf().sca(ax_gas)
        min_size = 100
        max_size = 10000
        alpha = 0.1
        x = self.gas.x.as_quantity_in(units.kpc)
        y = self.gas.y.as_quantity_in(units.kpc)
        z = self.gas.z.as_quantity_in(units.kpc)
        z, x, y, us, h_smooths = z.sorted_with(x, y, self.gas.u, self.gas.h_smooth)
        u_min, u_max = min(us), max(us)

        log_u = numpy.log((us / u_min)) / numpy.log((u_max / u_min))
        clipped_log_u = numpy.minimum(numpy.ones_like(log_u), numpy.maximum(numpy.zeros_like(log_u), log_u))

        red = 1.0 - clipped_log_u**4
        blue = clipped_log_u**4
        green = numpy.minimum(red, blue)

        colors = numpy.transpose(numpy.array([red, green, blue]))
        n_pixels = pyplot.gcf().get_dpi() * pyplot.gcf().get_size_inches()

        ax_gas.set_axis_bgcolor('#101010')
        ax_gas.set_aspect("equal", adjustable="datalim")
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
        sizes = numpy.minimum(numpy.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

        scatter(x, y, color=colors, s=sizes, edgecolors="none", alpha=alpha)
        xlabel(r'$x$')
        ylabel(r'$y$')

        xlim(-2.*virial_radius, 2*virial_radius)
        ylim(-2.*virial_radius, 2*virial_radius)
        pyplot.tight_layout()
        pyplot.show()
