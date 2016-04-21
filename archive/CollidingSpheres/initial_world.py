"""
File: initial_world.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Fri Mar 04, 2016 03:45 pm
Last modified: Wed Mar 09, 2016 12:44 pm

Description: Initial conidtions for Cygnus A cluster, and routines for merger.

"""

import os
import glob
from datetime import datetime

import pickle

import numpy

from matplotlib import pyplot
from matplotlib import animation
import matplotlib.gridspec as gridspec
pyplot.rcParams.update({'font.size': 22})

from PIL import Image
import images2gif # https://raw.githubusercontent.com/rec/echomesh/master/code/python/external/images2gif.py

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import VectorQuantity
from amuse.community.gadget2.interface import Gadget2
from amuse.plot import plot
from amuse.plot import scatter
from amuse.plot import xlabel
from amuse.plot import ylabel
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file

from initial_subcluster import SubCluster
from helper_functions import print_progressbar
from helper_functions import smart_length_units_for_vector_quantity


class ClusterMerger(object):
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

    def dm_gas_sph_subplot(self, time=0):
        fig = pyplot.figure(figsize=(20, 12))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4])
        # gs.update(left=0.05, right=0.48, wspace=0.05)

        ax_text = pyplot.subplot(gs[0, :])
        ax_text.axis('off')
        time_text = ax_text.text(0.02, 1.0, '', transform=ax_text.transAxes, fontsize=42)
        time_text.set_text('Time: {0:.1f} Myr'.format(self.code.model_time.value_in(units.Myr)))
        energy_text = ax_text.text(0.02, -0.2, '', transform=ax_text.transAxes, fontsize=42)

        Ekin = self.code.kinetic_energy.value_in(units.erg)
        Epot = self.code.potential_energy.value_in(units.erg)
        Eth = self.code.thermal_energy.value_in(units.erg)

        energy_text.set_text('Ekin: {0:.3e} erg\nEpot: {1:.3e} erg\nEth: {2:.3e} erg'
            .format(Ekin, Epot, Eth))


        # lim = max(abs((self.dmA.center_of_mass() - self.dmB.center_of_mass()).value_in(units.Mpc)))
        lim = 2.5
        ax_dm = pyplot.subplot(gs[1, 0], xlim=(-4*lim, 4*lim), ylim=(-4*lim, 4*lim))
        ax_gas = pyplot.subplot(gs[1, 1], aspect='equal',
            sharex=ax_dm, sharey=ax_dm, xlim=(-4*lim, 4*lim), ylim=(-4*lim, 4*lim))
        # ax_dm = fig.add_subplot(121, aspect='equal')
        # ax_gas = fig.add_subplot(122, aspect='equal',
        #    sharex=ax_dm, sharey=ax_dm)

        # plot dark matter
        pyplot.gcf().sca(ax_dm)
        x = self.dmA.x.as_quantity_in(units.Mpc)
        y = self.dmA.y.as_quantity_in(units.Mpc)
        scatter(x, y, c='red', edgecolor='red', label=str(self.subClusterA.name))
        x = self.dmB.x.as_quantity_in(units.Mpc)
        y = self.dmB.y.as_quantity_in(units.Mpc)
        scatter(x, y, c='green', edgecolor='green', label=str(self.subClusterB.name))
        xlabel(r'$x$')
        ylabel(r'$y$')
        pyplot.legend()

        # plot gas as sph plot
        def plot_sph(gas):
            # Adjusted code from amuse.plot.sph_particles_plot
            pyplot.gcf().sca(ax_gas)
            min_size = 100
            max_size = 10000
            alpha = 0.1
            x = gas.x
            y = gas.y
            z = gas.z
            z, x, y, us, h_smooths = z.sorted_with(x, y, gas.u, gas.h_smooth)
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
            length_unit = smart_length_units_for_vector_quantity(x)
            phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
            sizes = numpy.minimum(numpy.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

            scatter(x.as_quantity_in(length_unit), y.as_quantity_in(length_unit),
                    color=colors, s=sizes, edgecolors="none", alpha=alpha)

        plot_sph(self.gasA)
        plot_sph(self.gasB)
        xlabel(r'$x$')
        ylabel(r'$y$')

        pyplot.tight_layout()
        pyplot.savefig('out/{0}/plots/dm_gas_sph_subplot_{1}.png'
            .format(self.timestamp, time), dpi=50)
        # pyplot.show()
        pyplot.close()

    def create_gif(self):
        print "Generating gif :-)"

        # Convert matplotlib generated plots to a gif using imges2gif (requires Pillow!)
        images = []
        sorted_plots = sorted(glob.glob('out/{0}/plots/*.png'.format(self.timestamp)),
                              key=os.path.getmtime)
        for plot in sorted_plots:
            images.append(Image.open(plot))

        #with NamedTemporaryFile(suffix='.gif') as f:
            #print f.name
        self.gifname = 'out/{0}/plots/clustermerger.gif'.format(self.timestamp)
        images2gif.writeGif(self.gifname, images)

    def gif_to_html(self):
        IMG_TAG = """<img src="data:image/gif;base64,{0}"></img>"""

        # convert gif to html
        data = open(self.gifname, 'rb').read()
        data = data.encode('base64')

        return IMG_TAG.format(data)

    def plot_velocities(self):
        print "Plotting velocity as function of time"

        gasA_vel_list = [] | (units.km/units.s)
        dmA_vel_list = [] | (units.km/units.s)
        gasB_vel_list = [] | (units.km/units.s)
        dmB_vel_list = [] | (units.km/units.s)
        time_list = [] | units.Gyr

        for i, time in enumerate(self.timesteps):
            gasA = read_set_from_file("out/{0}/data/gasA_{1}.amuse"
                .format(self.timestamp, int(time.value_in(units.Myr))), "amuse")
            dmA = read_set_from_file("out/{0}/data/dmA_{1}.amuse"
                .format(self.timestamp, int(time.value_in(units.Myr))), "amuse")
            gasB = read_set_from_file("out/{0}/data/gasB_{1}.amuse"
                .format(self.timestamp, int(time.value_in(units.Myr))), "amuse")
            dmB = read_set_from_file("out/{0}/data/dmB_{1}.amuse"
                .format(self.timestamp, int(time.value_in(units.Myr))), "amuse")

            gasA_vel_list.append(gasA.center_of_mass_velocity().x)
            dmA_vel_list.append(dmA.center_of_mass_velocity().x)
            gasB_vel_list.append(gasB.center_of_mass_velocity().x)
            dmB_vel_list.append(dmB.center_of_mass_velocity().x)
            time_list.append(time)
        print

        fig = pyplot.figure(figsize=(12, 10), dpi=50)
        plot(time_list, gasA_vel_list, label="gasA", c="r", lw=10, ls="solid")
        plot(time_list, dmA_vel_list, label="dmA", c="r", lw=10, ls="dashed")
        plot(time_list, gasB_vel_list, label="gasB", c="g",lw=10, ls="solid")
        plot(time_list, dmB_vel_list, label="dmB", c="g", lw=10, ls="dashed")
        xlabel("Time")
        ylabel("Velocity")
        pyplot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        pyplot.show()

        print gasA_vel_list
        print dmA_vel_list
        print gasB_vel_list
        print dmB_vel_list
        print time_list

    def plot_energies(self):
        print "Plotting energies as function of time"

        thermal_list = [] | (units.erg)
        kinetic_list = [] | (units.erg)
        potential_list = [] | (units.erg)
        total_list = [] | (units.erg)
        time_list = [] | units.Gyr

        paths = sorted(glob.glob('out/{0}/data/gasA_*'.format(self.timestamp)), key=os.path.getmtime)
        tot = len(paths) - 1
        for i, path in enumerate(paths):
            time = path[28:-6]  # fix due to round error when storing the data.
            print_progressbar(i, tot)
            gasA = read_set_from_file("out/{0}/data/gasA_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            dmA = read_set_from_file("out/{0}/data/dmA_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            gasB = read_set_from_file("out/{0}/data/gasB_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            dmB = read_set_from_file("out/{0}/data/dmB_{1}.amuse"
                .format(self.timestamp, time), "amuse")

            thermal = gasA.thermal_energy() + gasB.thermal_energy()
            kinetic = gasA.kinetic_energy() + dmA.kinetic_energy() + gasB.kinetic_energy() + dmB.kinetic_energy()
            potential = gasA.potential_energy() + dmA.potential_energy() + gasB.potential_energy() + dmB.potential_energy()

            thermal_list.append(thermal)
            kinetic_list.append(kinetic)
            potential_list.append(potential)
            total_list.append(thermal+kinetic+potential)
            time_list.append((int(time) | units.Myr))
        print

        fig = pyplot.figure(figsize=(12, 10), dpi=50)
        plot(time_list, thermal_list, label="thermal", c="m", lw=10, ls="solid")
        plot(time_list, kinetic_list, label="kinetic", c="g", lw=10, ls="solid")
        plot(time_list, potential_list, label="potential", c="y", lw=10, ls="solid")
        plot(time_list, total_list, label="total", c="r", lw=10, ls="solid")
        xlabel("Time")
        ylabel("Energy")
        pyplot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        pyplot.show()

        print thermal_list
        print kinetic_list
        print potential_list
        print total_list
        print time_list

    def plot_densities(self):
        print "Plotting x-density as function of time"

        if not os.path.exists('out/{0}/plots/density'.format(self.timestamp)):
            os.mkdir('out/{0}/plots/density'.format(self.timestamp))

        density_gasA_list = [] | (units.MSun/units.Mpc**3)
        density_dmA_list = [] | (units.MSun/units.Mpc**3)
        density_gasB_list = [] | (units.MSun/units.Mpc**3)
        density_dmB_list = [] | (units.MSun/units.Mpc**3)
        time_list = [] | units.Gyr

        paths = sorted(glob.glob('out/{0}/data/gasA_*'.format(self.timestamp)), key=os.path.getmtime)
        tot = len(paths) - 1
        for i, path in enumerate(paths):
            time = path[28:-6]  # fix due to round error when storing the data.
            print_progressbar(i, tot)
            gasA = read_set_from_file("out/{0}/data/gasA_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            dmA = read_set_from_file("out/{0}/data/dmA_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            gasB = read_set_from_file("out/{0}/data/gasB_{1}.amuse"
                .format(self.timestamp, time), "amuse")
            dmB = read_set_from_file("out/{0}/data/dmB_{1}.amuse"
                .format(self.timestamp, time), "amuse")

            density_gasA = gasA.mass/(4./3 * numpy.pi * gasA.h_smooth**3)
            density_dmA = dmA.mass/(4./3 * numpy.pi * dmA.radius**3)
            density_gasB = gasB.mass/(4./3 * numpy.pi * gasB.h_smooth**3)
            density_dmB = dmA.mass/(4./3 * numpy.pi * dmA.radius**3)

            fig = pyplot.figure(figsize=(12, 10), dpi=50)
            plot(gasA.x, density_gasA, label="gasA", c="r", ls="solid")
            plot(dmA.x, density_dmA, label="dmA", c="r", ls="dashed")
            plot(gasB.x, density_gasB, label="gasB", c="g", ls="solid")
            plot(dmB.x, density_dmB, label="dmB", c="g", ls="dashed")
            xlabel("x")
            ylabel("Density")
            pyplot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            pyplot.show()
            pyplot.savefig('out/{0}/plots/density/density_{1}'.format(self.timestamp, time), dpi=50)
            pyplot.close()

            if i > 10:
                break
        print

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

        pyplot.savefig(figname, transparent=True, dpi = 100)
        print "\nHydroplot was saved to: ", figname
        pyplot.close()

    def stop(self):
        print "Stopping the code. End of pipeline :-)"
        self.code.stop()


if __name__ == "__main__":
    from amuse.support.console import set_printing_strategy
    set_printing_strategy("custom", preferred_units=[units.MSun, units.kpc, units.Myr, units.kms],
                          precision=4, prefix = "", separator = " [", suffix = "]")

    merger = ClusterMerger()
    merger.evolve()
