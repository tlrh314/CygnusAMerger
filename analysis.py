"""
File: analysis.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Mon Apr 18, 2016 02:28 pm

Check if initial conditions file was parsed correctly, and check
the stability of the pas in the dark matter halo.
"""

import numpy
# Bigger fontsize is better *_*
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})
# pyplot.rcParams.update({"text.usetex": True})

from amuse.units import units
from amuse.units.quantities import VectorQuantity
import amuse.plot as amuse_plot

from initial import Cluster


def plot_individual_cluster_density(cluster):
    """ Plot the particles' density radial profile and compare to model """
    pyplot.figure(figsize=(24, 18))

    # AMUSE datamodel particles. Gas has RHO and RHOm; dm rho from model.
    amuse_plot.scatter(cluster.gas.r,
        cluster.gas.rho,
        c="g", edgecolor="face", s=1, label=r"Generated IC: gas $\rho$")
    # amuse_plot.scatter(cluster.gas.r,
    #    cluster.gas.rhom,
    #    c="r", edgecolor="face", s=1, label=r"Generated IC: gas $\rho_{\rm model}$")
    amuse_plot.scatter(cluster.dm.r,
        cluster.dm.rho,
        c="b", edgecolor="none", label=r"Generated IC: DM $\rho$")

    # Analytical solutions. Sample radii and plug into analytical expression.
    r = VectorQuantity.arange(units.kpc(1), units.kpc(10000), units.parsec(100))

    # Plot analytical beta model (Donnert 2014) for the gas density
    amuse_plot.plot(r, cluster.gas_density_double_beta(r), c="k", ls="dashed",
        label=r"Analytical, $\beta$-model:" "\n" r"$\rho_0$ = {0} g/cm$^3$; $rc = ${1} kpc".format(
        cluster.rho0gas.number, cluster.rc.number))
    # Plot analytical double beta model (Donnert et al. 2016, in prep) for gas
    amuse_plot.plot(r, cluster.gas_density_beta(r), c="k", ls="dotted",
        label=r"Analytical, double $\beta$-model:" "\n" r"$\rho_0$ = {0} g/cm$^3$; $rc =$ {1} kpc; $r_{{\rm cut}}$ = {2} kpc".format(
        cluster.rho0gas.number, cluster.rc.number, cluster.rcut.number))

    # Plot analytical Hernquist model for the DM density
    amuse_plot.plot(r, cluster.dm_density(r), c="k", ls="solid",
        label=r"Analytical, Hernquist-model" "\n" r"$M_{{\rm dm}}= ${0:.2e} MSun; $a = $ {1} kpc".format(
        cluster.M_dm.number, cluster.a.number))

    pyplot.legend(loc=3)
    amuse_plot.xlabel(r"$r$")
    amuse_plot.ylabel(r"$\rho$")
    pyplot.gca().set_xlim(xmin=10, xmax=1e4)
    pyplot.gca().set_ylim(ymin=1e-30, ymax=9e-24)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    pyplot.axvline(x=cluster.R200.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(cluster.R200.value_in(units.kpc), 5e-24, r"$r_{{cut}} =$ {0}".format(cluster.rcut))
    pyplot.axvline(x=cluster.rc.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(cluster.rc.value_in(units.kpc), 5e-24, r"$rc =$ {0}".format(cluster.rc))
    pyplot.axvline(x=cluster.a.value_in(units.kpc), lw=1, c="grey")
    pyplot.text(cluster.a.value_in(units.kpc), 1e-24, r"$a =$ {0}".format(cluster.a))


def get_mass_via_number_density(cluster):
    # progressbar
    import sys
    pbwidth = 42

    # TODO: find different method to calculate M(<r)
    print "Counting particles for which radii < r to obtain M(<r)"
    radii = VectorQuantity.arange(units.kpc(1), units.kpc(10000), units.parsec(1000))
    M_gas_below_r = numpy.zeros(len(radii))
    M_dm_below_r = numpy.zeros(len(radii))
    N = len(radii)
    for i, r in enumerate(radii):
        M_gas_below_r[i] = ((numpy.where(cluster.gas.r < r)[0]).size)
        M_dm_below_r[i] = ((numpy.where(cluster.dm.r < r)[0]).size)
        if i%1000 == 0:
            # update the bar
            progress = float(i+1000)/N
            block = int(round(pbwidth*progress))
            text = "\rProgress: [{0}] {1:.1f}%".format( "#"*block + "-"*(pbwidth-block), progress*100)
            sys.stdout.write(text)
            sys.stdout.flush()

    sys.stdout.write("\n")
    print "Done counting particles :-)... TODO: improve this method?!"

    M_gas_below_r *= cluster.M_gas/cluster.raw_data.Ngas
    M_dm_below_r *= cluster.M_dm/cluster.raw_data.Ndm

    amuse_plot.scatter(radii, M_gas_below_r, c="g", label="Gas")
    amuse_plot.scatter(radii, M_dm_below_r, c="b", label="DM")


def get_mass_via_density(cluster):
    gas_i = numpy.argsort(cluster.gas.r.value_in(units.kpc))
    gas_r = cluster.gas.r[gas_i].value_in(units.kpc)
    gas_rho = (cluster.gas.rho[gas_i]).value_in(units.MSun/units.kpc**3)
    gas_mass = (4./3*numpy.pi*gas_r**3*gas_rho)
    gas_mass_cummulative = numpy.cumsum(gas_mass)

    dm_i = numpy.argsort(cluster.dm.r.value_in(units.kpc))
    dm_r = cluster.dm.r[dm_i].value_in(units.kpc)
    dm_rho = (cluster.dm.rho[dm_i]).value_in(units.MSun/units.kpc**3)
    dm_mass = (4./3*numpy.pi*dm_r**3*dm_rho)
    dm_mass_cummulative = numpy.cumsum(dm_mass)

    pyplot.scatter(gas_r, gas_mass_cummulative)
    pyplot.scatter(dm_r, dm_mass_cummulative)
    # pyplot.gca().set_xscale("log")
    # pyplot.gca().set_yscale("log")

    pyplot.scatter(gas_r, gas_mass)
    pyplot.scatter(dm_r, dm_mass)

    #amuse_plot.scatter(cluster.gas.r, M_gas_below_r, c="g", label="Gas")
    #amuse_plot.scatter(cluster.dm.r, M_dm_below_r, c="b", label="DM")


def plot_individual_cluster_mass(cluster):
    pyplot.figure(figsize=(24,18))

    # TODO: use different method :-)...
    get_mass_via_number_density(cluster)
    # get_mass_via_number_density_parallel(cluster)
    # get_mass_via_density(cluster)

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    # M_gas = (cluster.gas.mass.value_in(units.MSun))
    # M_dm = (cluster.dm.mass.value_in(units.MSun))
    # amuse_plot.scatter(cluster.gas.r, units.MSun(M_gas), c="g", label="Gas")
    # amuse_plot.scatter(cluster.dm.r, units.MSun(M_dm), c="b", label="DM")

    # Analytical solutions. Sample radii and plug into analytical expression.
    r = VectorQuantity.arange(units.kpc(1), units.kpc(10000), units.parsec(100))

    # Plot analytical Hernquist model for the DM mass M(<r)
    amuse_plot.plot(r, cluster.dm_cummulative_mass(r), ls="solid", c="k", label="DM")
    # Plot analytical beta model (Donnert 2014) for the gas mass M(<r)
    amuse_plot.plot(r, cluster.gas_cummulative_mass_beta(r),
        ls="dotted", c="k", label=r"$\beta$-model")
    # Plot analytical double beta model (Donnert et al. 2016, in prep) gas M(<r)
    amuse_plot.plot(r, cluster.gas_cummulative_mass_double_beta(r),
        ls="dashed", c="k", label=r"double $\beta$-model")


    pyplot.xlabel(r"$r$ [kpc]")
    pyplot.ylabel(r"$M (<r)$ [MSun]")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.legend(loc=4)


def cutoff_beta_model(r, rho0, rc, rcut):
    beta = 2./3
    return rho0 * (1 + (r/rc)**2)**(-3*beta/2) / (1 + (r/rcut)**4);


if __name__ == "__main__":
    datadir="../../workdir/ToyclusterICs/20160420T1852/"
    logfile="runToycluster.log"
    icfile="IC_single_0"

    cluster = Cluster(datadir)
    cluster.perform_sanity_checks()
    plot_individual_cluster_density(cluster)
    plot_individual_cluster_mass(cluster)
    pyplot.show()
