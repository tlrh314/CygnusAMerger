import glob
import os
import numpy
import argparse

from deco import concurrent, synchronized

import matplotlib
from matplotlib import pyplot
from plotsettings import PlotSettings
style = PlotSettings()
pyplot.rcParams.update({"xtick.major.width": 4})
pyplot.rcParams.update({"ytick.major.width": 4})

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


class Run(object):
    def __init__(self, arguments, discard_firstbins=True):
        self.timestamp = arguments.simulationID[0]
        self.observed = ObservedCluster(arguments.clustername) if arguments.clustername else None
        if arguments.clustername and arguments.clustername == "cygA" and discard_firstbins:
            print "WARNING: Discarding first three CygA bins.\n"
            self.observed.radius = self.observed.radius[3:]
            self.observed.binsize = self.observed.binsize[3:]
            self.observed.density = self.observed.density[3:]
            self.observed.density_std = self.observed.density_std[3:]
            self.observed.number_density = self.observed.number_density[3:]
            self.observed.number_density_std = self.observed.number_density_std[3:]
        if arguments.IC_only:
            # TODO: there is something different in the Toycluster rho
            # than in the Gadget output :(
            self.IC_only = True
            self.icdir = "../runs/{0}/ICs/".format(self.timestamp)
            self.snapdir = "../runs/{0}/ICs/".format(self.timestamp)
            self.snaps = ["../runs/{0}/ICs/IC_single_0".format(self.timestamp)]
        else:
            self.IC_only = False
            self.icdir = "../runs/{0}/ICs/".format(self.timestamp)
            self.snapdir = "../runs/{0}/snaps/".format(self.timestamp)
            self.snaps = sorted(glob.glob("../runs/{0}/snaps/snapshot_*"
                .format(self.timestamp)),  key=os.path.getmtime)
            print self.snaps, "\n"

        self.outdir="../runs/{0}/out/".format(self.timestamp)
        if not (os.path.isdir(self.outdir) or os.path.exists(self.outdir)):
            os.mkdir(self.outdir)

        if self.IC_only:
            self.TimeBetSnapshot = 0
        else:
            # For time counter
            gadgetparms = parse_gadget_parms(self.snapdir+"gadget2.par")
            self.TimeBetSnapshot = gadgetparms['TimeBetSnapshot']


def plot_individual_cluster_density(numerical, analytical, observed=None,
                                    poster_style=False):
    """ Plot the particles' density radial profile and compare to model """
    pyplot.figure(figsize=(12, 9))
    if poster_style:
        pyplot.style.use(["dark_background"])
        pyplot.rcParams.update({"font.weight": "bold"})
        # magenta, dark blue, orange, green, light blue (?)
        fit_colour = "white"
    else:
        fit_colour = "k"
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]

    if observed:
        obs_colour = "k" if observed.name == "cygA" else "k"
        pyplot.errorbar(observed.radius+observed.binsize/2, observed.density,
                        xerr=observed.binsize/2, yerr=observed.density_std,
                        marker="o", ms=7 if poster_style else 3, alpha=1,
                        elinewidth=5 if poster_style else 1, ls="",
                        c=obs_colour, label=observed.name)

    # Gas RHO and RHOm from Toycluster living in AMUSE datamodel
    # cut = numpy.where(numerical.gas.r < numerical.R200)
    amuse_plot.scatter(numerical.gas.r, numerical.gas.rho,
        c=data_colour[3], edgecolor="face", s=10 if poster_style else 4,
        label=r"Gas, sampled")

    # amuse_plot.scatter(numerical.gas.r,
    #    numerical.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")

    # DM density obtained from mass profile (through number density)
    amuse_plot.scatter(numerical.dm_radii, numerical.rho_dm_below_r,
       c=data_colour[1], edgecolor="face", s=10 if poster_style else 4,
       label=r"DM, sampled")

    # Analytical solutions.

    # Plot analytical cut-off beta model (Donnert 2014) for the gas density
    amuse_plot.plot(analytical.radius, analytical.gas_density(), c=fit_colour,
            lw=5 if poster_style else 1, ls="dashed",
        label="Gas: "+analytical.modelname)#+"\n"+\
            #r"$\rho_0$ = {0:.3e} g/cm$^3$; $rc = ${1:.2f} kpc"\
            #.format(analytical.rho0.number, analytical.rc.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(analytical.radius, analytical.dm_density(), c=fit_colour,
            lw=5 if poster_style else 1, ls="solid",
        label=r"DM: Hernquist")# "\n" r"$M_{{\rm dm}}= ${0:.2e} MSun; $a = $ {1} kpc".format(
        #analytical.M_dm.number, analytical.a.number))

    amuse_plot.ylabel(r"Density [g/cm$^3$]")
    amuse_plot.xlabel(r"Radius [kpc]")

    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.gca().set_ylim(ymin=1e-32, ymax=9e-24)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    pyplot.axvline(x=numerical.R200.value_in(units.kpc),
                   lw=1, c=fit_colour)
    pyplot.text(numerical.R200.value_in(units.kpc)+100, 4e-24, r"$r_{200}$",
                ha="left", fontsize=22)
    pyplot.fill_between(numpy.arange(numerical.R200.value_in(units.kpc), 1e4, 0.01),
            1e-32, 9e-24, facecolor="grey", edgecolor="grey", alpha=0.1)

    ymin = analytical.gas_density(numerical.rc.as_quantity_in(units.kpc))
    pyplot.vlines(x=numerical.rc.value_in(units.kpc),
                  ymin=ymin.value_in(units.g/units.cm**3), ymax=9e-24,
                  lw=1, linestyles="dashed", color=fit_colour)
    pyplot.text(numerical.rc.value_in(units.kpc)+25
                if (observed and observed.name == "cygB") else
                numerical.rc.value_in(units.kpc), 4e-24, r"$r_c$",
                ha="left", fontsize=22)

    ymin = analytical.dm_density(numerical.a.as_quantity_in(units.kpc))
    pyplot.vlines(x=numerical.a.value_in(units.kpc),
                  ymin=ymin.value_in(units.g/units.cm**3), ymax=9e-24,
                  lw=1, linestyles="solid", color=fit_colour)
    pyplot.text(numerical.a.value_in(units.kpc)-25, 4e-24, r"$a$",
                ha="right", fontsize=22)

    # Maximum value in histogram of numerical.gas.h
    ymin = analytical.gas_density(107.5 | units.kpc)
    pyplot.vlines(x=107.5, ymin=ymin.value_in(units.g/units.cm**3), ymax=9e-24,
                  lw=4, linestyles="dashed", color=data_colour[3])
    pyplot.text(107.5-10, 4e-24, r"$2h_{\rm sml}$", color=data_colour[3],
                ha="right", fontsize=22)

    # pyplot.savefig("out/actuallyran.png")
    # pyplot.legend(loc=3, fontsize=28)
    pyplot.tight_layout()
    # pyplot.show()
    # pyplot.savefig("out/no_wvt_relax.pdf", dpi=300)

def plot_individual_cluster_mass(numerical, analytical, poster_style=False):
    pyplot.figure(figsize=(16, 12))
    if poster_style:
        pyplot.style.use(["dark_background"])
        pyplot.rcParams.update({"font.weight": "bold"})
        fit_colour = "white"
    else:
        fit_colour = "k"
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]

    amuse_plot.scatter(numerical.gas_radii, numerical.M_gas_below_r,
        c=data_colour[0], edgecolor="face", s=4 if poster_style else 1,
        label="Gas, sampled")

    amuse_plot.scatter(numerical.dm_radii, numerical.M_dm_below_r,
        c=data_colour[2], edgecolor="face", s=4 if poster_style else 1,
        label="DM, sampled")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.dm_cummulative_mass(),
        ls="solid", c=fit_colour, lw=5 if poster_style else 1, label="DM: Hernquist")

    # Plot analytical cut-off beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.gas_mass(), ls="dotted",
            c=fit_colour, lw=5 if poster_style else 1, label="Gas: "+analytical.modelname)

    pyplot.xlabel(r"Radius [kpc]")
    pyplot.ylabel(r"Cummulative Mass [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=4)

def plot_individual_cluster_temperature(numerical, analytical, poster_style=False):
    pyplot.figure(figsize=(16, 12))
    if poster_style:
        pyplot.style.use(["dark_background"])
        pyplot.rcParams.update({"font.weight": "bold"})
        fit_colour = "white"
    else:
        fit_colour = "k"
    # magenta, dark blue, orange, green, light blue (?)
    data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                   (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                   (41./255, 239./255, 239./255)]

    numerical.set_gas_temperature()
    pyplot.scatter(numerical.gas.r.value_in(units.kpc),
        numerical.gas.T.value_in(units.K), c=data_colour[0], edgecolor="face",
        s=5 if poster_style else 1, label="Gas, sampled")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical temperature (Donnert 2014)
    # NB the temperature is set by both gas and dark matter
    # because it follows from hydrostatic equation, which contains M(<r),
    # where M(<r) is the total gravitating mass delivered by gas+dm!
    amuse_plot.plot(analytical.radius, analytical.temperature(), ls="dotted",
        c=fit_colour, lw=5 if poster_style else 1, label="Analytical")
    pyplot.axhline(analytical.characteristic_temperature().value_in(units.K))
    pyplot.axvline(analytical.rc.value_in(units.kpc))

    pyplot.xlabel(r"Radius [kpc]")
    pyplot.ylabel(r"Temperature [K]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.gca().set_xlim(xmin=1, xmax=1e4)
    pyplot.legend(loc=3)


#@concurrent(processes=4)
def make_all_plots(run, snap, i, replot=False):
    print "Generating plots for:", snap

    fname = snap.split('/')[-1]
    snapnr = fname.split('_')[-1]
    if run.IC_only:
        snapnr += "00-ICOnly"

    mass_filename = run.outdir+"{0}-mass-{1}.png".format(run.timestamp, snapnr)
    density_filename = run.outdir+"{0}-density-{1}.png".format(run.timestamp, snapnr)
    temperature_filename = run.outdir+"{0}-temperature-{1}.png".format(run.timestamp, snapnr)

    if (os.path.isfile(mass_filename) and os.path.exists(mass_filename)
            and os.path.isfile(density_filename) \
            and os.path.exists(density_filename)) \
            and not replot:
        print "Plots exist. Skipping", snap
        return

    numerical = NumericalCluster(
        icdir=run.icdir,
        snapdir=run.snapdir,
        logfile="runToycluster.log",
        icfile="IC_single_0" if run.IC_only else "snapshot_"+snapnr,
        verbose=False)

    """ TODO: can use cluster.py setup_analytical_cluster
              to obtain AnalyticalCluster from NumericalCluster"""
    # Caution: parms[0] is number density! Caution: use unitsless numbers!
    ne0 = convert.rho_to_ne(numerical.rho0gas.value_in(units.g/units.cm**3),
        numerical.z)
    parms = (ne0, numerical.rc.value_in(units.kpc),
        numerical.R200.value_in(units.kpc), float(numerical.beta))
    # TODO: not use this, but if doing so then use halo0_sampling...
    dm_parms = (numerical.Mass_in_DM, numerical.a)

    analytical = AnalyticalCluster(parms, dm_parms, z=numerical.z, free_beta=True)

    # 295 for WC6, 50 for Cubic Spline kernel
    numerical.get_gas_mass_via_density(
        DESNNGB=50 if arguments.IC_only else 50)
    numerical.get_dm_mass_via_number_density(log_binning=True)
    numerical.set_dm_density()

    plot_individual_cluster_density(numerical, analytical, run.observed)
    # pyplot.title("Time = {0:1.2f} Gyr".format(i*run.TimeBetSnapshot), y=1.03)
    if arguments.save:
        pyplot.savefig(density_filename)
        pyplot.close()

    pyplot.show()
    # plot_individual_cluster_mass(numerical, analytical)
    # pyplot.title("Time = {0:1.2f} Gyr".format(i*run.TimeBetSnapshot), y=1.03)
    # if arguments.save:
    #     pyplot.savefig(mass_filename)
    #     pyplot.close()

    # plot_individual_cluster_temperature(numerical, analytical)
    # pyplot.title("Time = {0:1.2f} Gyr".format(i*run.TimeBetSnapshot), y=1.03)
    # if arguments.save:
    #     pyplot.savefig(temperature_filename)
    #     pyplot.close()
    # else:
    #     pyplot.show()

    # break
    print "Done with snapshot:  ", snap


#@synchronized
def generate_plots(arguments):
    run = Run(arguments)

    for i, snap in enumerate(run.snaps):
        if i != 5:
            continue
        make_all_plots(run, snap, i, replot=arguments.replot)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        description="Plot rho(r), M(<r), T(r) for given SimulationID.")
    parser.add_argument("-o", "--onlyic", dest="IC_only",
        action="store_true", default=False,
        help="only use ICfile, else use Gadget snapshots")
    parser.add_argument("-s", "--show", dest="save",
        action="store_false", default=True,
        help="show instead of save plots; default is save")
    parser.add_argument("-r", "--replot", dest="replot",
        action="store_true", default=False,
        help="if plot already exists, plot it anyway")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="simulationID", nargs=1,
        help="string of the Simulation ID")
    parser.add_argument("-c", "--cluster", dest="clustername", default=None,
        choices=["cygA", "cygB"], help="also plot observed cluster")

    return parser

if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    generate_plots(arguments)
