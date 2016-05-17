# import glob
# import os
import numpy
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})
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

    # Gas RHO and RHOm from Toycluster living in AMUSE datamodel
    amuse_plot.scatter(numerical.gas.r,
        numerical.gas.rho,
        c="g", edgecolor="face", s=1, label=r"Gas")

    # amuse_plot.scatter(numerical.gas.r,
    #    numerical.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")

    # DM density obtained from mass profile (through number density)
    amuse_plot.scatter(numerical.dm_radii, numerical.rho_dm_below_r,
       c="b", edgecolor="none", label=r"DM")

    # Analytical solutions.

    # Plot analytical cut-off beta model (Donnert 2014) for the gas density
    amuse_plot.plot(analytical.radius, analytical.gas_density(), c="k", ls="dashed",
        label=analytical.modelname+"\n"+\
            r"$\rho_0$ = {0:.3e} g/cm$^3$; $rc = ${1:.2f} kpc"\
            .format(analytical.rho0.number, analytical.rc.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(analytical.radius, analytical.dm_density(), c="k", ls="solid",
        label=r"Analytical, Hernquist-model" "\n" r"$M_{{\rm dm}}= ${0:.2e} MSun; $a = $ {1} kpc".format(
        analytical.M_dm.number, analytical.a.number))

    pyplot.legend(loc=3, prop={'size': 12})
    amuse_plot.xlabel(r"$r$")
    amuse_plot.ylabel(r"$\rho$")
    pyplot.gca().set_xlim(xmin=10, xmax=1e4)
    pyplot.gca().set_ylim(ymin=1e-30, ymax=9e-24)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    pyplot.axvline(x=numerical.R200.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(numerical.R200.value_in(units.kpc), 5e-24, r"$r_{{cut}} =$ {0}".format(numerical.rcut))
    pyplot.axvline(x=numerical.rc.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(numerical.rc.value_in(units.kpc), 5e-24, r"$rc =$ {0}".format(numerical.rc))
    pyplot.axvline(x=numerical.a.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(numerical.a.value_in(units.kpc), 1e-24, r"$a =$ {0}".format(numerical.a))


def plot_individual_cluster_mass(numerical, analytical):
    pyplot.figure(figsize=(12, 9))

    amuse_plot.scatter(numerical.gas_radii, numerical.M_gas_below_r,
        c="g", edgecolor="face", label="Gas")

    amuse_plot.scatter(numerical.dm_radii, numerical.M_dm_below_r,
        c="b", edgecolor="face", label="DM")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # Analytical solutions. Sample radii and plug into analytical expression.

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.dm_cummulative_mass(), ls="solid", c="k", label="DM")

    # Plot analytical cut-off beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(analytical.radius, analytical.gas_mass(),
        ls="dotted", c="k", label=analytical.modelname)

    pyplot.xlabel(r"$r$ [kpc]")
    pyplot.ylabel(r"$M (<r)$ [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.legend(loc=4)


if __name__ == "__main__":
    icdir = "../runs/no_wvt_relax/ICs/"
    snapdir = "../runs/no_wvt_relax/ICs/"
    numerical_cluster = NumericalCluster(
        icdir="../runs/no_wvt_relax/ICs/",
        snapdir="../runs/no_wvt_relax/ICs/",
        logfile="runToycluster.log",
        icfile="IC_single_0")

    # Caution: parms[0] is number density! Caution: use unitsless numbers!
    ne0 = (numerical_cluster.rho0gas.value_in(units.g/units.cm**3)) / \
        (globals.mu*globals.m_p)
    parms = (ne0, numerical_cluster.rc.value_in(units.kpc),
        numerical_cluster.R200.value_in(units.kpc))
    # TODO: not use this, but if doing so then use halo0_sampling...
    dm_parms = numerical_cluster.Mass_in_DM, numerical_cluster.a

    analytical_cluster = AnalyticalCluster(parms, dm_parms)

    numerical_cluster.get_gas_mass_via_density()
    numerical_cluster.get_dm_mass_via_number_density()
    numerical_cluster.set_dm_density()

    plot_individual_cluster_density(numerical_cluster, analytical_cluster)
    plot_individual_cluster_mass(numerical_cluster, analytical_cluster)
    pyplot.show()

