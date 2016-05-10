import numpy
import pandas
import csv
import scipy
from scipy import stats

import matplotlib
matplotlib.rcParams.update({'font.size': 22})
from matplotlib import pyplot

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import VectorQuantity
from amuse.support.console import set_printing_strategy
# set_printing_strategy("custom", preferred_units = [units.RSun, units.MSun, units.Myr, units.MSun/units.yr, units.RSun/units.yr])
# set_printing_strategy("cgs")
import amuse.plot as amuse_plot

from macro import *
from cosmology import CosmologyCalculator

# Bug in the version I have currently installed. Nice...
# Latest amuse.plot.errorbar
def errorbar(*args, **kwargs):
    for label in ['yerr', 'xerr']:
        if label in kwargs:
            args += (kwargs.pop(label),)
        else:
            args += (None,)

    yerr, xerr = args[2:4]

    args1 = amuse_plot.UnitlessArgs.strip(*args[:2])
    if xerr is not None:
        xerr = amuse_plot.UnitlessArgs.value_in_x_unit(xerr)
    if yerr is not None:
        yerr = amuse_plot.UnitlessArgs.value_in_y_unit(yerr)
    args = args1 + [yerr, xerr]
    result = pyplot.errorbar(*args, **kwargs)
    pyplot.xlabel(amuse_plot.UnitlessArgs.x_label())
    pyplot.ylabel(amuse_plot.UnitlessArgs.y_label())
    return result


def Gas_density_profile(parm, r):
    """ Double beta profile at rc and rcut """
    rho0 = parm[0]
    rc = parm[1]
    rcut = parm[2]

    return rho0 / (1 + p2(r/rc)) / (1 + p3(r/rcut) * (r/rcut))

# Define the statistical model, in this case we shall use a chi-squared distribution, assuming normality in the errors
def stat(parm, x, y, dy):
    ymodel = Gas_density_profile(parm, x)
    chisq = numpy.sum((y - ymodel)**2 / dy**2)
    return(chisq)

# Read Martijn Data
raw_cygA = pandas.read_csv("data/cygA_1T_pressureprofile.dat", delimiter="|")
cygA_bin_number = raw_cygA["Bin number "].as_matrix()
cygA_bin_volume = raw_cygA[" Bin volume (cm^3) "].as_matrix()
cygA_density = raw_cygA["   Density (cm^-3) "].as_matrix()
cygA_density_std = raw_cygA["     Sigma density "].as_matrix()
cygA_pressure = raw_cygA[" Pressure (erg cm^-3) "].as_matrix()
cygA_pressure_std = raw_cygA["    Sigma Pressure "].as_matrix()
cygA_compton_y = raw_cygA[" Compton Y parameter"].as_matrix()

raw_cygA_sn100 = pandas.read_csv("data/cygA_sn100_sbprofile.dat", delimiter="|")
cygA_inner_radius = raw_cygA_sn100["Inner radius (arcsec) "].as_matrix()
cygA_outer_radius = raw_cygA_sn100[" Outer radius (arcsec) "].as_matrix()
cygA_source_sb = raw_cygA_sn100[" Source SB (counts/cm^2/arcsec^2/s) "].as_matrix()
cygA_source_sb_std = raw_cygA_sn100["   Sigma source SB "].as_matrix()
cygA_background_sb = raw_cygA_sn100[" Background SB(counts/cm^2/arcsec^2/s) "].as_matrix()
cygA_background_sb_std = raw_cygA_sn100[" Sigma Background SB "].as_matrix()
cygA_bin_number_sn100 = raw_cygA_sn100[" Bin number"].as_matrix()

raw_cygB = pandas.read_csv("data/cygB_1T_pressureprofile.dat", delimiter="|")
cygB_bin_number = raw_cygB["Bin number "].as_matrix()
cygB_bin_volume = raw_cygB[" Bin volume (cm^3) "].as_matrix()
cygB_density = raw_cygB["   Density (cm^-3) "].as_matrix()
cygB_density_std = raw_cygB["     Sigma density "].as_matrix()
cygB_pressure = raw_cygB[" Pressure (erg cm^-3) "].as_matrix()
cygB_pressure_std = raw_cygB["    Sigma Pressure "].as_matrix()
cygB_compton_y = raw_cygB[" Compton Y parameter"].as_matrix()

raw_cygB_sn100 = pandas.read_csv("data/cygB_sn100_sbprofile.dat", delimiter="|")
cygB_inner_radius = raw_cygB_sn100["Inner radius (arcsec) "].as_matrix()
cygB_outer_radius = raw_cygB_sn100[" Outer radius (arcsec) "].as_matrix()
cygB_source_sb = raw_cygB_sn100[" SB (cnts/cm^2/arcsec^2/s) "].as_matrix()
cygB_source_sb_std = raw_cygB_sn100["   Sigma source SB "].as_matrix()
cygB_background_sb = raw_cygB_sn100[" Bkg SB(cnts/cm^2/arcsec^2/s) "].as_matrix()
cygB_background_sb_std = raw_cygB_sn100["      Sigma Bkg SB "].as_matrix()
cygB_bin_number_sn100 = raw_cygB_sn100[" Bin number"].as_matrix()

# print cygA_bin_number == cygA_bin_number_sn100
# print cygB_bin_number == cygB_bin_number_sn100

m_p = constants.proton_mass.as_quantity_in(units.g)
mu = 0.17
cc = CosmologyCalculator()
arcsec2kpc = cc.kpc_DA | units.kpc

cygA_radius = cygA_outer_radius * arcsec2kpc
cygA_mass_density = (cygA_density | units.cm**-3) * mu * m_p
cygA_mass_density_std = (cygA_density_std | units.cm**-3) * mu * m_p
cygB_radius = (cygB_outer_radius * arcsec2kpc).as_quantity_in(units.kpc)
cygB_mass_density = (cygB_density | units.cm**-3) * mu * m_p
cygB_mass_density_std = (cygB_density_std | units.cm**-3) * mu * m_p

# For analytical purposes
radius = VectorQuantity.arange(units.kpc(1), max(cygA_radius), units.parsec(100))

# plot Martijn data
fit_density = True
if fit_density:
    # set_printing_strategy("custom", preferred_units = [units.kpc, units.cm**-3*units.g])
    f, (axA, axB) = pyplot.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 8))

    # CygA
    pyplot.sca(axA)
    pyplot.errorbar(cygA_radius.value_in(units.kpc), cygA_density,
        yerr=cygA_density_std, label="CygA")
    # amuse_plot.errorbar(cygA_radius.as_quantity_in(units.kpc), cygA_mass_density,
    #     yerr=cygA_mass_density_std, label="CygA")
    # amuse_plot.plot(radius, Gas_density_profile(((5e-26 | units.g/units.cm**3), units.kpc(30), units.kpc(1000)), radius),
    #     label="eyebal 'fit' CygA")

    # Fit
    parm = [9e-28, 200, 1000]
    # parm = [1., 1., 1.]
    bounds = [(None, None), (None, None), (None, None)]
    result = scipy.optimize.minimize(stat, parm,
        args=(cygA_radius.value_in(units.kpc), cygA_density,
              cygA_density_std),
        method='L-BFGS-B', bounds=bounds)
    print result
    rho0_fit = result["x"][0] | units.g/units.cm**3
    rc_fit = result["x"][1] | units.kpc
    rcut_fit = result["x"][2] | units.kpc
    amuse_plot.plot(radius, Gas_density_profile((rho0_fit, rc_fit, rcut_fit), radius), label="CygB Fit")

    # CygB
    pyplot.sca(axB)
    pyplot.errorbar(cygB_radius.value_in(units.kpc), cygB_density,
        yerr=cygB_density_std, label="CygB")
    # pyplot.errorbar(cygB_radius.as_quantity_in(units.kpc), cygB_mass_density,
    #     yerr=cygB_mass_density_std, label="CygB")
    # amuse_plot.plot(radius, Gas_density_profile(((9e-28 | units.g/units.cm**3), units.kpc(200), units.kpc(1000)), radius),
    #    label="eyebal 'fit' CygB")

    # Fit
    parm = [2e-4, 200, 1000]
    # parm = [1., 1., 1.]
    bounds = [(None, None), (None, None), (None, None)]
    result = scipy.optimize.minimize(stat, parm,
        args=(cygB_radius.value_in(units.kpc), cygB_density,
              cygB_density_std),
        method='L-BFGS-B', bounds=bounds)
    print result
    rho0_fit = result["x"][0] | units.g/units.cm**3
    rc_fit = result["x"][1] | units.kpc
    rcut_fit = result["x"][2] | units.kpc
    amuse_plot.plot(radius, Gas_density_profile((rho0_fit, rc_fit, rcut_fit), radius), label="CygB Fit")


    # ml_func = result["fun"]
    # dof = len(dens) - len(ml_vals)
    # print "[MLEs], chisq/dof:", ml_vals, ml_func/dof
    # ch = scipy.stats.chi2(dof)
    # pval = 1.0 - ch.cdf(ml_func)
    # print "The p-value for the {0} distribution is: {1:.5f}".format(dists[i], pval)

    for ax in [axA, axB]:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$n_e$")
        ax.legend(loc=1)
    pyplot.show()
    # set_printing_strategy("cgs")

plot_pressure = False
if plot_pressure:
    pyplot.figure(figsize=(12, 9))
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    amuse_plot.errorbar(cygA_bin_number, (cygA_pressure | units.erg/units.cm**3),
        yerr=(cygA_pressure_std | units.erg/units.cm**3), label="CygA")
    amuse_plot.errorbar(cygB_bin_number, (cygB_pressure | units.erg/units.cm**3),
        yerr=(cygB_pressure_std | units.erg/units.cm**3), label="CygB")
    pyplot.xlabel("Bin number")  # TODO: change to radius. How?
    pyplot.ylabel(r'$P$')
    pyplot.legend()
    pyplot.show()

plot_bin_volume = False
if plot_bin_volume:
    pyplot.figure(figsize=(12, 9))
    pyplot.title("Bin Volume")
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.plot(cygA_bin_number, cygA_bin_volume, label="CygA")
    pyplot.plot(cygB_bin_number, cygB_bin_volume, label="CygB")
    pyplot.legend()
    pyplot.savefig("cygAB_bin_volume.png")
    pyplot.show()

plot_compton_y = False
if plot_compton_y:
    pyplot.figure(figsize=(12, 9))
    pyplot.title("Compton Y")
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.plot(cygA_bin_number, cygA_compton_y, label="CygA")
    pyplot.plot(cygB_bin_number, cygB_compton_y, label="CygB")
    pyplot.legend()
    pyplot.show()

# pyplot.figure(figsize=(12, 9))
# pyplot.gca().set_xscale("log")
# pyplot.gca().set_yscale("log")
# amuse_plot.errorbar(cygA_bin_number, cygA_density*cygA_bin_volume,
#     yerr=cygA_density_std, label="CygA")
# amuse_plot.errorbar(cygB_bin_number, cygB_density*cygB_bin_volume,
#     yerr=cygB_density_std, label="CygB")
# pyplot.xlabel("Bin number")  # TODO: change to radius. How?
# pyplot.ylabel(r'$n_e\cdot V$')
# pyplot.legend()
# pyplot.show()
