import glob
import os
import numpy
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 28})
pyplot.rcParams.update({"xtick.major.size": 8})
pyplot.rcParams.update({"xtick.minor.size": 4})
pyplot.rcParams.update({"ytick.major.size": 8})
pyplot.rcParams.update({"ytick.minor.size": 4})
pyplot.rcParams.update({"xtick.major.width": 4})
pyplot.rcParams.update({"xtick.minor.width": 2})
pyplot.rcParams.update({"ytick.major.width": 4})
pyplot.rcParams.update({"ytick.minor.width": 2})
pyplot.rcParams.update({"xtick.major.pad": 8})
pyplot.rcParams.update({"xtick.minor.pad": 8})
pyplot.rcParams.update({"ytick.major.pad": 8})
pyplot.rcParams.update({"ytick.minor.pad": 8})
pyplot.rcParams.update({"legend.loc": "best"})
pyplot.rcParams.update({"figure.autolayout": True})

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import VectorQuantity
import amuse.plot as amuse_plot

from ioparser import Gadget2BinaryF77UnformattedType2Parser
from ioparser import Toycluster2RuntimeOutputParser
from ioparser import parse_gadget_parms
from cluster import NumericalCluster
from cluster import AnalyticalCluster
from cluster import ObservedCluster
import convert


def plot_individual_cluster_density(numerical, analytical, observed=None):
    """ Plot the particles' density radial profile and compare to model """
    pyplot.figure(figsize=(16, 12))
    pyplot.style.use(["dark_background"])
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]
    fit_colour = "white"

    if observed:
        pyplot.errorbar(observed.radius+observed.binsize/2,
                        observed.density, xerr=observed.binsize/2,
                        yerr=observed.density_std, marker="o",
                        ms=7, elinewidth=5, ls="", c=data_colour[3],
                        label=observed.name)


    # Gas RHO and RHOm from Toycluster living in AMUSE datamodel
    # cut = numpy.where(numerical.gas.r < numerical.R200)
    amuse_plot.scatter(numerical.gas.r, numerical.gas.rho,
        c=data_colour[0], edgecolor="face", s=1, label=r"Gas, sampled")

    # amuse_plot.scatter(numerical.gas.r,
    #    numerical.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")

    # DM density obtained from mass profile (through number density)
    amuse_plot.scatter(numerical.dm_radii, numerical.rho_dm_below_r,
       c=data_colour[2], edgecolor="face", s=10, label=r"DM, sampled")

    # Analytical solutions.

    # Plot analytical cut-off beta model (Donnert 2014) for the gas density
    amuse_plot.plot(analytical.radius, analytical.gas_density(), c=fit_colour, lw=1, ls="dotted",
        label="Gas: "+analytical.modelname)#+"\n"+\
            #r"$\rho_0$ = {0:.3e} g/cm$^3$; $rc = ${1:.2f} kpc"\
            #.format(analytical.rho0.number, analytical.rc.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(analytical.radius, analytical.dm_density(), c=fit_colour, lw=1, ls="solid",
        label=r"DM: Hernquist")# "\n" r"$M_{{\rm dm}}= ${0:.2e} MSun; $a = $ {1} kpc".format(
        #analytical.M_dm.number, analytical.a.number))

    amuse_plot.ylabel(r"Density [g/cm$**$3]")
    amuse_plot.xlabel(r"Radius [kpc]")

    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.gca().set_ylim(ymin=1e-32, ymax=9e-24)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    pyplot.axvline(x=numerical.R200.value_in(units.kpc), lw=1, c=fit_colour)
    pyplot.text(numerical.R200.value_in(units.kpc), 3e-24, r"r200 = {0:.2f} kpc".format(numerical.rcut.value_in(units.kpc)))
    pyplot.axvline(x=numerical.rc.value_in(units.kpc), lw=1, c=fit_colour)
    pyplot.text(numerical.rc.value_in(units.kpc), 1e-24, r"rc = {0:.2f} kpc".format(numerical.rc.value_in(units.kpc)))
    # pyplot.axvline(x=numerical.a.value_in(units.kpc), lw=1, c=fit_colour)
    # pyplot.text(numerical.a.value_in(units.kpc), 1e-24, r"$a =$ {0}".format(numerical.a))

    # pyplot.savefig("out/actuallyran.png")
    pyplot.legend(loc=3)#, prop={"size": 22})


def plot_individual_cluster_mass(numerical, analytical):
    pyplot.figure(figsize=(16, 12))
    pyplot.style.use(["dark_background"])
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]
    fit_colour = "white"

    amuse_plot.scatter(numerical.gas_radii, numerical.M_gas_below_r,
        c=data_colour[0], edgecolor="face", s=1, label="Gas, sampled")

    amuse_plot.scatter(numerical.dm_radii, numerical.M_dm_below_r,
        c=data_colour[2], edgecolor="face", s=1, label="DM, sampled")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.dm_cummulative_mass(), ls="solid", c=fit_colour, lw=2, label="DM: Hernquist")

    # Plot analytical cut-off beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.gas_mass(),
        ls="dotted", c=fit_colour, lw=2, label="Gas: "+analytical.modelname)

    pyplot.xlabel(r"Radius [kpc]")
    pyplot.ylabel(r"Cummulative Mass [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=4)


def plot_individual_cluster_temperature(numerical, analytical):
    pyplot.figure(figsize=(16, 12))
    pyplot.style.use(["dark_background"])
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]
    fit_colour = "white"

    numerical.set_gas_temperature()
    pyplot.scatter(numerical.gas.r.value_in(units.kpc),
        numerical.gas.T.value_in(units.K),
        c=data_colour[0], edgecolor="face", s=1, label="Gas, sampled")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical temperature (Donnert 2014)
    # NB the temperature is set by both gas and dark matter
    # because it follows from hydrostatic equation, which contains M(<r),
    # where M(<r) is the total gravitating mass delivered by gas+dm!
    amuse_plot.plot(analytical.radius, analytical.temperature(),
        ls="dotted", c=fit_colour, lw=2, label="Analytical")
    pyplot.axhline(analytical.characteristic_temperature().value_in(units.K))
    pyplot.axvline(analytical.rc.value_in(units.kpc))

    pyplot.xlabel(r"Radius [kpc]")
    pyplot.ylabel(r"Temperature [K]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=3)


if __name__ == "__main__":
    IC_only = True
    save = True
    replot = True

    myRun = "yes_wvt_relax"
    myRun = "no_wvt_relax"
    myRun = "20160712T2017"
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
        print snaps

    outdir="../runs/{0}/out/".format(myRun)
    if not (os.path.isdir(outdir) or os.path.exists(outdir)):
        os.mkdir(outdir)

    for i, snap in enumerate(snaps):
        fname = snap.split('/')[-1]
        snapnr = fname.split('_')[-1]

        if IC_only:
            TimeBetSnapshot = 0
            snapnr += "00-ICOnly"
        else:
            # For time counter
            gadgetparms = parse_gadget_parms(snapdir+"gadget2.par")
            TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

        mass_filename = outdir+"{0}-mass-{1}.png".format(myRun, snapnr)
        density_filename = outdir+"{0}-density-{1}.png".format(myRun, snapnr)
        temperature_filename = outdir+"{0}-temperature-{1}.png".format(myRun, snapnr)

        if (os.path.isfile(mass_filename) and os.path.exists(mass_filename)
                and os.path.isfile(density_filename) \
                and os.path.exists(density_filename)) \
                and not replot:
            print "Plots exist. Skipping", snap
            continue

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

        # 295 for WC6, 50 for Cubic Spline kernel
        numerical.get_gas_mass_via_density(DESNNGB=50 if IC_only else 50)
        numerical.get_dm_mass_via_number_density(log_binning=True)
        numerical.set_dm_density()

        plot_individual_cluster_density(numerical, analytical)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*TimeBetSnapshot), y=1.03)
        if save:
            pyplot.savefig(density_filename)
            pyplot.close()

        plot_individual_cluster_mass(numerical, analytical)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*TimeBetSnapshot), y=1.03)
        if save:
            pyplot.savefig(mass_filename)
            pyplot.close()

        plot_individual_cluster_temperature(numerical, analytical)
        pyplot.title("Time = {0:1.2f} Gyr".format(i*TimeBetSnapshot), y=1.03)
        if save:
            pyplot.savefig(temperature_filename)
            pyplot.close()
        else:
            pyplot.show()

        # break

        print "Done checking snapshot:", snap
