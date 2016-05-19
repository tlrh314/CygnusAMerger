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


def plot_individual_cluster_density(numerical, analytical):
    """ Plot the particles' density radial profile and compare to model """
    pyplot.figure(figsize=(16, 12))
    pyplot.style.use(["dark_background"])
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]
    fit_colour = "white"


    # Gas RHO and RHOm from Toycluster living in AMUSE datamodel
    amuse_plot.scatter(numerical.gas.r,
        numerical.gas.rho,
        c=data_colour[0], edgecolor="face", s=40, label=r"Gas")

    # amuse_plot.scatter(numerical.gas.r,
    #    numerical.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")

    # DM density obtained from mass profile (through number density)
    amuse_plot.scatter(numerical.dm_radii, numerical.rho_dm_below_r,
       c=data_colour[2], edgecolor="face", s=40, label=r"DM")

    # Analytical solutions.

    # Plot analytical cut-off beta model (Donnert 2014) for the gas density
    amuse_plot.plot(analytical.radius, analytical.gas_density(), c=fit_colour, lw=4, ls="dotted",
        label=analytical.modelname+"\n"+\
            r"$\rho_0$ = {0:.3e} g/cm$^3$; $rc = ${1:.2f} kpc"\
            .format(analytical.rho0.number, analytical.rc.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(analytical.radius, analytical.dm_density(), c=fit_colour, lw=4, ls="solid",
        label=r"Analytical, Hernquist-model" "\n" r"$M_{{\rm dm}}= ${0:.2e} MSun; $a = $ {1} kpc".format(
        analytical.M_dm.number, analytical.a.number))

    amuse_plot.xlabel(r"$r$")
    amuse_plot.ylabel(r"$\rho$")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.gca().set_ylim(ymin=1e-30, ymax=9e-24)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    pyplot.axvline(x=numerical.R200.value_in(units.kpc), lw=1, c=fit_colour)
    pyplot.text(numerical.R200.value_in(units.kpc), 5e-24, r"$r_{{cut}} =$ {0}".format(numerical.rcut))
    pyplot.axvline(x=numerical.rc.value_in(units.kpc), lw=1, c=fit_colour)
    pyplot.text(numerical.rc.value_in(units.kpc), 5e-24, r"$rc =$ {0}".format(numerical.rc))
    pyplot.axvline(x=numerical.a.value_in(units.kpc), lw=1, c=fit_colour)
    pyplot.text(numerical.a.value_in(units.kpc), 1e-24, r"$a =$ {0}".format(numerical.a))

    pyplot.legend(loc=3, prop={'size': 12})


def plot_individual_cluster_mass(numerical, analytical):
    pyplot.figure(figsize=(16, 12))
    pyplot.style.use(["dark_background"])
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]
    fit_colour = "white"


    amuse_plot.scatter(numerical.gas_radii, numerical.M_gas_below_r,
        c=data_colour[0], edgecolor="face", s=40, label="Gas")

    amuse_plot.scatter(numerical.dm_radii, numerical.M_dm_below_r,
        c=data_colour[2], edgecolor="face", s=40, label="DM")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.dm_cummulative_mass(), ls="solid", c=fit_colour, lw=4, label="DM")

    # Plot analytical cut-off beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.gas_mass(),
        ls="dotted", c=fit_colour, lw=4, label=analytical.modelname)

    pyplot.xlabel(r"$r$ [kpc]")
    pyplot.ylabel(r"$M (<r)$ [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=4)


if __name__ == "__main__":
    snaps = sorted(glob.glob("../runs/yes_wvt_relax/snaps/snapshot_*"),  key=os.path.getmtime)

    for i, snap in enumerate(snaps):
        fname = snap.split('/')[-1]
        snapnr = fname.split('_')[-1]

        numerical_cluster = NumericalCluster(
            icdir="../runs/yes_wvt_relax/ICs/",
            snapdir="../runs/yes_wvt_relax/snaps/",
            logfile="runToycluster.log",
            icfile="snapshot_"+snapnr)

        # Caution: parms[0] is number density! Caution: use unitsless numbers!
        ne0 = (numerical_cluster.rho0gas.value_in(units.g/units.cm**3)) / \
            (globals.mu*globals.m_p)
        parms = (ne0, numerical_cluster.rc.value_in(units.kpc),
            numerical_cluster.R200.value_in(units.kpc))
        # TODO: not use this, but if doing so then use halo0_sampling...
        dm_parms = numerical_cluster.Mass_in_DM, numerical_cluster.a

        analytical_cluster = AnalyticalCluster(parms, dm_parms)

        fast_dm_density_computation = False
        if fast_dm_density_computation:
            # magenta, dark blue, orange, green, light blue (?)
            data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                           (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                           (41./255, 239./255, 239./255)]
            print numerical_cluster.dm
            radii = numpy.arange(1, 10000, 1)
            particles = numpy.zeros(len(radii))
            dr = radii[1] - radii[0]
            print dr
            for i, r in enumerate(radii):
                particles[i] = ((numpy.where(numerical_cluster.dm.r.value_in(units.kpc) < r)[0]).size)

                print i, r, particles[i]

            dparticle = numpy.zeros(len(particles))
            for i in range(1, len(particles)):
                dparticle[i-1] = particles[i] - particles[i-1]

            print dparticle

            pyplot.figure(figsize=(12,9))
            volume2 = 4./3 * numpy.pi * radii**3
            volume = 4 * numpy.pi * radii**2 * dr
            density = dparticle/(numerical_cluster.raw_data.Ndm*volume*numerical_cluster.M_dm.number)
            density2 = particles/(numerical_cluster.raw_data.Ndm*volume2*numerical_cluster.M_dm.number)
            pyplot.scatter(radii, 0.17*density, c=data_colour[0], edgecolor="face", label="numerical")
            pyplot.scatter(radii, 0.17*density2, c=data_colour[2], edgecolor="face", label="numerical cubed")
            amuse_plot.plot(analytical_cluster.radius, analytical_cluster.dm_density(), label="analytical")
            # amuse_plot.plot(analytical_cluster.radius, analytical_cluster.gas_density(), label="gas")
            pyplot.gca().set_xscale("log")
            pyplot.gca().set_yscale("log")
            pyplot.legend(loc=3)
            pyplot.show()

            import sys; sys.exit(0)

        numerical_cluster.get_gas_mass_via_density()
        numerical_cluster.get_dm_mass_via_number_density()
        numerical_cluster.set_dm_density()

        plot_individual_cluster_density(numerical_cluster, analytical_cluster)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*0.25))
        pyplot.savefig("out/yes_wvt_relax-density-"+snapnr+".png")
        pyplot.close()
        plot_individual_cluster_mass(numerical_cluster, analytical_cluster)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*0.25))
        pyplot.savefig("out/yes_wvt_relax-mass-"+snapnr+".png")
        pyplot.close()

        print "Done checking snapshot:", snap
