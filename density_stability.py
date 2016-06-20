import glob
import os
import numpy
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 22})
# pyplot.rcParams.update({"text.usetex": True})

from amuse.units import units
from amuse.units.quantities import VectorQuantity
import amuse.plot as amuse_plot

from ioparser import Gadget2BinaryF77UnformattedType2Parser
from ioparser import Toycluster2RuntimeOutputParser
from ioparser import parse_gadget_parms
from cluster import NumericalCluster
from cluster import AnalyticalCluster
import convert


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
        c=data_colour[0], edgecolor="face", s=20, label=r"Gas")

    # amuse_plot.scatter(numerical.gas.r,
    #    numerical.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")

    # DM density obtained from mass profile (through number density)
    amuse_plot.scatter(numerical.dm_radii, numerical.rho_dm_below_r,
       c=data_colour[2], edgecolor="face", s=20, label=r"DM")

    # Analytical solutions.

    # Plot analytical cut-off beta model (Donnert 2014) for the gas density
    amuse_plot.plot(analytical.radius, analytical.gas_density(), c=fit_colour, lw=2, ls="dotted",
        label=analytical.modelname+"\n"+\
            r"$\rho_0$ = {0:.3e} g/cm$^3$; $rc = ${1:.2f} kpc"\
            .format(analytical.rho0.number, analytical.rc.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(analytical.radius, analytical.dm_density(), c=fit_colour, lw=2, ls="solid",
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
        c=data_colour[0], edgecolor="face", s=20, label="Gas")

    amuse_plot.scatter(numerical.dm_radii, numerical.M_dm_below_r,
        c=data_colour[2], edgecolor="face", s=20, label="DM")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.dm_cummulative_mass(), ls="solid", c=fit_colour, lw=2, label="DM")

    # Plot analytical cut-off beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.gas_mass(),
        ls="dotted", c=fit_colour, lw=2, label=analytical.modelname)

    pyplot.xlabel(r"$r$ [kpc]")
    pyplot.ylabel(r"$M (<r)$ [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=4)


if __name__ == "__main__":
    IC_only = False

    myRun = "yes_wvt_relax"
    myRun = "no_wvt_relax"
    myRun = "20160620T1152"
    if IC_only:
        # TODO: there is something different in the Toycluster rho
        # than in the Gadget output :(
        icdir = "../runs/{0}/ICs/".format(myRun)
        snapdir = "../runs/{0}/ICs/".format(myRun)
        snaps = ["../runs/{0}/ICs/IC_single_0".format(myRun)]
    else:
        icdir = "../runs/{0}/ICs/".format(myRun)
        snapdir = "../runs/{0}/snaps/".format(myRun)
        snaps = sorted(glob.glob("../runs/{0}/snaps/snapshot_*".format(myRun)),  key=os.path.getmtime)

    for i, snap in enumerate(snaps):
        fname = snap.split('/')[-1]
        snapnr = fname.split('_')[-1]

        numerical = NumericalCluster(
            icdir=icdir,
            snapdir=snapdir,
            logfile="runToycluster.log",
            icfile="IC_single_0" if IC_only else "snapshot_"+snapnr)

        """ TODO: can use cluster.py setup_analytical_cluster
                  to obtain AnalyticalCluster from NumericalCluster"""
        # Caution: parms[0] is number density! Caution: use unitsless numbers!
        ne0 = convert.rho_to_ne(numerical.rho0gas.value_in(units.g/units.cm**3),
            numerical.z)
        parms = (ne0, numerical.rc.value_in(units.kpc),
            numerical.R200.value_in(units.kpc))
        # TODO: not use this, but if doing so then use halo0_sampling...
        dm_parms = (numerical.Mass_in_DM, numerical.a)

        analytical = AnalyticalCluster(parms, dm_parms, z=numerical.z)

        fast_dm_density_computation = False
        if fast_dm_density_computation:
            # magenta, dark blue, orange, green, light blue (?)
            data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                           (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                           (41./255, 239./255, 239./255)]
            print numerical.dm
            radii = numpy.arange(1, 10000, 100)
            particles = numpy.zeros(len(radii))
            dr = radii[1] - radii[0]
            print dr
            for i, r in enumerate(radii):
                particles[i] = ((numpy.where(numerical.dm.r.value_in(units.kpc) < r)[0]).size)

                print i, r, particles[i]

            dparticle = numpy.zeros(len(particles))
            for i in range(1, len(particles)):
                dparticle[i-1] = particles[i] - particles[i-1]

            print dparticle

            pyplot.figure(figsize=(12,9))
            volume = 4 * numpy.pi * radii**2 * dr
            density = dparticle/(numerical.raw_data.Ndm*volume*numerical.M_dm.number)
            density = density
            volume2 = 4./3 * numpy.pi * radii**3
            density2 = particles/(numerical.raw_data.Ndm*volume2*numerical.M_dm.number)
            # TODO: why 0.17? Is this density gas + dm??
            pyplot.scatter(radii, 0.17*density, c=data_colour[0], edgecolor="face", label="numerical")
            pyplot.scatter(radii, 0.17*density2, c=data_colour[2], edgecolor="face", label="numerical cubed")
            amuse_plot.plot(analytical.radius, analytical.dm_density(), label="analytical")
            # amuse_plot.plot(analytical.radius, analytical.gas_density(), label="gas")
            pyplot.gca().set_xscale("log")
            pyplot.gca().set_yscale("log")
            pyplot.xlim(5, 1e5)
            pyplot.ylim(1e-30, 1e-21)
            pyplot.legend(loc=3)
            pyplot.show()

            import sys; sys.exit(0)

        numerical.get_gas_mass_via_density()
        numerical.get_dm_mass_via_number_density()
        numerical.set_dm_density()

        plot_individual_cluster_density(numerical, analytical)

        if IC_only:
            TimeBetSnapshot = 0
        else:
            # For time counter
            gadgetparms = parse_gadget_parms(snapdir+"gadget2.par")
            TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

        pyplot.title("Time = {0:1.2f} Gyr".format(i*TimeBetSnapshot))
        pyplot.savefig("out/{0}-density-".format(myRun)+snapnr+".png")
        pyplot.close()
        plot_individual_cluster_mass(numerical, analytical)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*TimeBetSnapshot))
        pyplot.savefig("out/{0}-mass-".format(myRun)+snapnr+".png")
        pyplot.close()

        # break

        print "Done checking snapshot:", snap
