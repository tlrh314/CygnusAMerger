import os
import numpy
import pandas
import csv
import scipy
from scipy import stats
from collections import OrderedDict

import matplotlib
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('text', usetex=True)
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
from ioparser import parse_toycluster_parms
from ioparser import Toycluster2RuntimeOutputParser
from cluster import ObservedCluster
from cluster import AnalyticalCluster
from cluster import NumericalCluster
import convert


def Gas_density_profile(parm, r):
    """ Beta model if len(parm)==2 (Donnert 2014);
        Cut-off beta profile at rc and rcut if len(parm)==3 (Donnert et al. 2016, in prep)
        Double cut-off beta profile if len(parm)==5. See Toycluster/setup.c """
    DOUBLE_BETA_COOL_CORES = False
    rho0 = parm[0]
    rc = parm[1]

    rho = rho0 / (1 + p2(r/rc))   # not cut-off (Donnert 2014)
    if len(parm) >= 3:
        rcut = parm[2]
        rho /= (1 + p3(r/rcut) * (r/rcut))

    if len(parm) == 5:  # Double Beta model (for Coolcore clusters)
        Rho0_Fac = parm[3]
        Rc_Fac = parm[4]

        rho0_cc = rho0 * Rho0_Fac
        rc_cc = rc / Rc_Fac

        rho += rho0_cc / (1 + p2(r/rc_cc)) / (1 + p3(r/rcut) * (r/rcut))

    return rho


def Gas_density_profile_wrapper(r, parm0, parm1, parm2=None, parm3=None, parm4=None):
    """ Same same, but scipy.optimize.minimize expects different form of
        function than scipy.optimize.curve_fit does. """
    if parm3 and parm4:  # coolcore, double beta
        return Gas_density_profile((parm0, parm1, parm2, parm3, parm4), r)
    elif parm2:
        return Gas_density_profile((parm0, parm1, parm2), r)
    else:
        return Gas_density_profile((parm0, parm1), r)


def Mass_profile(r, rho0, rc, rcut):
    """ return M(<= R) of a double beta profile with beta=2/3  """
    r2 = p2(r)
    rc2 = p2(rc)
    rcut2 = p2(rcut)
    sqrt2 = numpy.sqrt(2)

    A = (rc2 - rcut2)*(numpy.log(rcut2 - sqrt2*rcut*r + r2) \
        - numpy.log(rcut2 + sqrt2*rcut*r + r2))
    Bplus = 2*(rc2 + rcut2)*numpy.arctan(1 + sqrt2*r/rcut)
    Bmin = 2*(rc2 + rcut2)*numpy.arctan(1 - sqrt2*r/rcut)

    # NB the paper is slightly different: equation 2 does not contain 4*pi*rho
    M_gas_below_r = (4 * numpy.pi * rho0) *\
        rc2*p3(rcut)/(8*(p2(rcut2)+p2(rc2))) *\
        (sqrt2*(A - Bmin + Bplus) - 8*rc*rcut*numpy.arctan(r/rc))
    return M_gas_below_r


def concentration_parameter(cluster, gas_mass):
    A = 5.74
    B = -0.097
    C = -0.47
    Mpivot = 2e12 / cluster.cc.h

    UnitMass_in_g = 1.989e43  # 1e10 Msun
    Msol2cgs = 1.98892e33

    bf = 0.17
    total_mass = gas_mass * (1+1/bf)  #Mgas+Mdm, Mdm = Mgas/bf

    # total mass has to be in 1e10 Msun

    mass = total_mass*UnitMass_in_g/Msol2cgs

    c_NFW = A * pow(mass/Mpivot, B) * pow(1 + cluster.cc.z, C)

    return c_NFW


# Define the statistical model, in this case we shall use a chi-squared distribution, assuming normality in the errors
def stat(parm, x, y, dy):
    ymodel = Gas_density_profile(parm, x)
    chisq = numpy.sum((y - ymodel)**2 / dy**2)
    return(chisq)


def fit_betamodel_to_chandra(observed, parm=[1., 1., 1.], double=False):
    """ Fit betamodel to cluster observation, or
        double betamodel when len(parm)==5. """

    # Fit to data
    if observed.name == "cygA":
        if len(parm) == 3:
            bounds = [(None, None), (None, None), (800, 1400)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1400), (None, None), (None, None)]
        else:
            bounds = [(None, None), (100, 400)]
    elif observed.name == "cygB":
        if len(parm) == 3:
            bounds = [(None, None), (None, None), (600, 800)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1400), (None, None), (None, None)]
        else:
            bounds = [(None, None), (None, None)]
    result = scipy.optimize.minimize(stat, parm,
            args=(observed.radius, observed.number_density, observed.number_density_std),
            method='L-BFGS-B', bounds=bounds)

    # Obtain and print MLEs
    obtain_mles(observed, result)

    return result


def obtain_mles(observed, result):
    ml_vals = result["x"]
    ml_func = result["fun"]

    moddof = len(ml_vals)  # Model degrees of freedom; nr of fit parameters
    dof = len(observed.number_density) - moddof  # degrees of freedom

    ch = scipy.stats.chi2(dof)
    pval = 1.0 - ch.cdf(ml_func)

    print observed.name
    print "Results for the '{0}' model:".format("Cut-off beta")
    print "  Using scipy.optimize.minimize to minimize chi^2 yields:"
    print "    n_e,0       = {0:.3f}".format(ml_vals[0])
    print "    r_c         = {0:.3f}".format(ml_vals[1])
    if len(ml_vals) == 3:
        print "    r_cut       = {0:.3f}".format(ml_vals[2])
    print "    chisq/dof   = {0:.3f}".format(ml_func/dof)
    print "    p-value     = {0:.5f}".format(pval)

    ml_vals, ml_covar = scipy.optimize.curve_fit(
            Gas_density_profile_wrapper,
            observed.radius, observed.number_density, p0=ml_vals,
            sigma=observed.number_density_std,)
    # ml_funcval = stat(ml_vals, edges, dens, err, model)

    if not result["success"]:
        print "  scipy.optimize.curve_fit broke down!\n    Reason: '{0}'".format(result["message"])
        print "  No confidence intervals have been calculated."

    err = numpy.sqrt(numpy.diag(ml_covar))
    print "  Using scipy.optimize.curve_fit to obtain confidence intervals yields:"
    print "    n_e,0       = {0:.3f} +/- {1:.3f}".format(ml_vals[0], err[0])
    print "    r_c         = {0:.3f} +/- {1:.3f}".format(ml_vals[1], err[1])
    if len(ml_vals) == 3:
        print "    r_cut       = {0:.3f} +/- {1:.3f}".format(ml_vals[2], err[2])
    print


def plot_fit_results(observed, analytical, numerical=None, mass_density=False, save=False):
    poster_style = True
    if poster_style:
        pyplot.style.use(["dark_background"])
        data_colour = (255./255, 64./255, 255./255)
        fit_colour = "white"
    else:
        data_colour = "r"
        fit_colour = "k"

    if mass_density:
        observed_density = observed.density
        observed_density_std = observed.density_std
        analytical_density = analytical.gas_density().value_in(units.g/units.cm**3)
    else:  # using number density as density
        observed_density = observed.number_density
        observed_density_std = observed.number_density_std
        analytical_density = analytical.gas_number_density().value_in(1/units.cm**3)

    fig, (ax, ax_r) = pyplot.subplots(2, 2, sharex=True, figsize=(16, 12))
    gs1 = matplotlib.gridspec.GridSpec(3, 3)
    gs1.update(hspace=0)
    ax = pyplot.subplot(gs1[:-1,:])
    ax_r = pyplot.subplot(gs1[-1,:])  # residuals

    # Plot data
    pyplot.sca(ax)
    pyplot.title(observed.name)
    pyplot.errorbar(observed.radius+observed.binsize/2, observed_density, xerr=observed.binsize/2,
            yerr=observed_density_std, marker='o', ms=6, ls='', c=data_colour,
            label="800 ks Chandra\n(De Vries, 2016)",)

    label = analytical.modelname+"\n\t"
    label +=  r"$n_{{e,0}} \,$ = {0:.2e}".format(analytical.ne0.number) if not mass_density else r"$\rho_{{0}} \,$ = {0:.2e}".format(analytical.rho0.number)
    label += "\n\t"+ r"$r_c \,$ = {0:.2f}".format(analytical.rc.number)
    if analytical.rcut is not None:
        label += "\n\t"+r"$r_{{\rm cut}}$ = {0:.2f}".format(analytical.rcut.number)
    if analytical.rho0_cc is not None:
        label += "\n\t"+r"$n_{{e,0,cc}}$ = {0:.2f}".format(analytical.ne0_cc.number) if not mass_density else "\n\t"+r"$\rho_{{0,cc}}$ = {0:.2f}".format(analytical.rho0_cc.number)
        label += "\n\t"+r"$r_{{c,cc}}$ = {0:.2f}".format(analytical.rc_cc.number)
    amuse_plot.plot(analytical.radius, analytical_density,
            c=fit_colour, lw=1, label=label)  #, drawstyle="steps-mid")

    if numerical:
        pyplot.scatter(numerical.gas.r.value_in(units.kpc),
            numerical.gas.rho.value_in(units.g / (units.cm**3)),
            c="g", edgecolor="face", s=4, label=r"Gas, numerical")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$ [kpc]")
    ax.set_ylabel(r"$n_e$ [cm$^{-3}$]")
    ax.legend(loc=3, prop={'size':20})

    # Plot Residuals
    pyplot.sca(ax_r)
    residual_density = (observed_density - analytical_density)/observed_density
    pyplot.errorbar(observed.radius+observed.binsize/2, residual_density,
            yerr=observed_density_std, c=fit_colour, drawstyle="steps-mid")
    ax_r.axhline(y=0, lw=2, ls="dashed", c="white")

    ax_r.set_xscale("log")
    # ax_r.set_yscale("log")

    ax.set_xlim(min(observed.radius)-0.3, max(observed.radius)+2000)
    ax_r.set_xlim(min(observed.radius)-0.3, max(observed.radius)+2000)

    if observed.name == "cygA" and not mass_density:
        ax_r.set_ylim(-0.3, 0.3)
    if observed.name == "cygB" and not mass_density:
        ax_r.set_ylim(-0.3, 0.3)

    # TODO: set residuals limit for mass density and avoid overlap...
    if observed.name == "cygA" and mass_density:
        ax_r.set_ylim(-0.3, 0.3)
    if observed.name == "cygB" and mass_density:
        ax_r.set_ylim(-0.3, 0.3)

    # Fix for overlapping y-axis markers
    from matplotlib.ticker import MaxNLocator
    ax.tick_params(labelbottom='off')
    nbins = len(ax_r.get_yticklabels())
    ax_r.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

    if save:
        pyplot.savefig("out/density_betamodel_fit_{0}.png".format(observed.name), dpi=600)
    # pyplot.show()


# TODO: place this in AnalyticalCluster?
def get_cluster_mass_analytical(observed, result):
    poster_style = True
    if poster_style:
        pyplot.style.use(["dark_background"])
        data_colour = (255./255, 64./255, 255./255)
        fit_colour = "white"
    else:
        data_colour = "g"
        fit_colour = "b"

    ml_vals = result["x"]

    ne0_fit = ml_vals[0]  # | 1/units.cm**3
    rc_fit = ml_vals[1]  # | units.kpc
    if len(ml_vals) == 3:
        rcut_fit = ml_vals[2]  # | units.kpc
    if len(ml_vals) == 5:
        Rho0_Fac = ml_vals[3]
        Rc_Fac = ml_vals[4]

    kpc2cm = units.kpc(1).value_in(units.cm)
    rho0_analytical = convert.ne_to_rho(ne0_fit) # g/cm**3
    rc_analytical = rc_fit * kpc2cm  # cm
    rhocrit200 = 200*observed.cc.rho_crit()

    analytical_radius = numpy.arange(1, 1e4, 0.001) * kpc2cm
    analytical_density = Gas_density_profile((rho0_analytical, rc_analytical), analytical_radius)
    analytical_mass = cummulative_mass_profile(analytical_radius, rc_analytical,  rho0_analytical)
    rho_average = (analytical_mass / (4./3*numpy.pi*analytical_radius**3))
    r200_analytical = analytical_radius[(numpy.abs(rho_average-rhocrit200)).argmin()]

    print observed.name
    print "  r200 = {0:.4e} kpc".format((r200_analytical | units.cm).value_in(units.kpc))
    # M200_analytical = analytical_mass[(numpy.abs(analytical_radius-r200_analytical)).argmin()]
    # print "M200 =", (M200_analytical | units.g).value_in(units.MSun)
    M200_analytical = cummulative_mass_profile(r200_analytical, rc_analytical, rho0_analytical)
    print "  M200 = {0:.4e} MSun".format((M200_analytical | units.g).value_in(units.MSun))

    fig, (ax1, ax2) = pyplot.subplots(1, 2, sharex=True, figsize=(18, 9))
    pyplot.suptitle("Obtaining M200 from analytical fit to observed density profile")
    # Step 1: calculate average density
    # Step 2: find where average density equals 200 times critical density at given z (0.0562)
    # Step 3: find r200
    # Step 4: find M200 given the above obtained r200
    pyplot.sca(ax1)
    amuse_plot.loglog((analytical_radius | units.cm).as_quantity_in(units.kpc),
        (analytical_density | units.g/units.cm**3), c=fit_colour, label="analytical")
    amuse_plot.loglog((analytical_radius | units.cm).as_quantity_in(units.kpc),
        (rho_average | units.g/units.cm**3), c=data_colour, label="average")
    pyplot.axhline(rhocrit200, ls="dashed", c=data_colour)
    pyplot.axvline((r200_analytical | units.cm).value_in(units.kpc), ls="dashed", c=data_colour)
    amuse_plot.ylabel(r"$\rho_{\rm gas}(r)$")
    amuse_plot.xlabel(r"$r$")
    pyplot.legend(loc=2)

    pyplot.sca(ax2)
    amuse_plot.loglog((analytical_radius | units.cm).as_quantity_in(units.kpc),
        (analytical_mass | units.g).as_quantity_in(units.MSun), c=fit_colour)
    pyplot.axhline((M200_analytical | units.g).value_in(units.MSun), ls="dashed", c=data_colour)
    pyplot.axvline((r200_analytical | units.cm).value_in(units.kpc), ls="dashed", c=data_colour)
    amuse_plot.ylabel(r"$M(<r)$")
    amuse_plot.xlabel(r"$r$")
    pyplot.savefig("out/m200_analytical_{0}.png".format(observed.name), dpi=600)


def cummulative_mass_profile(r, r_c, rho_0):
    return 4*numpy.pi*r_c**3*rho_0 * (r/r_c - numpy.arctan(r/r_c))


# TODO: Place this in ObservedCluster and get it working
def get_mass_profile(observed, result):
    ml_vals = result["x"]

    ne0_fit = ml_vals[0]  # | 1/units.cm**3
    rc_fit = ml_vals[1]  # | units.kpc
    if len(ml_vals) == 3:
        rcut_fit = ml_vals[2]  # | units.kpc
    if len(ml_vals) == 5:
        Rho0_Fac = ml_vals[3]
        Rc_Fac = ml_vals[4]

    kpc2cm = units.kpc(1).value_in(units.cm)
    rho0_analytical = convert.ne_to_rho(ne0_fit)  # g/cm**3
    rc_analytical = rc_fit * kpc2cm  # cm

    analytical_radius = numpy.arange(min(observed.radius), max(observed.radius), 0.01) * kpc2cm
    analytical_mass = cummulative_mass_profile(analytical_radius, rc_analytical,  rho0_analytical)

    mass = observed.density * observed.bin_volume
    volume = 4./3*numpy.pi*(observed.radius*kpc2cm)**3  # cm cubed
    mass2 = observed.density * volume
    # numpy.cumsum() actually works :-)....
    # mass_summed = numpy.zeros(len(mass))
    # mass_summed[0] = mass[0]
    # mass2_summed = numpy.zeros(len(mass))
    # mass2_summed[0] = mass2[0]
    # for i in range(1, len(mass)):
    #     mass_summed[i] = mass_summed[i-1] + mass[i]
    #     mass2_summed[i] = mass2_summed[i-1] + mass2[i]

    pyplot.figure(figsize=(12, 9))
    amuse_plot.plot((analytical_radius | units.cm).as_quantity_in(units.kpc),
        (analytical_mass | units.g).as_quantity_in(units.MSun),
        label="Analytical")
    # amuse_plot.plot((observed.radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
    #     (mass_summed | units.g).as_quantity_in(units.MSun),
    #     label="Martijn selfsummed")
    amuse_plot.plot((observed.radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
        (mass.cumsum() | units.g).as_quantity_in(units.MSun),
        label="Martijn")
    # amuse_plot.plot((observed.outer_radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
    #     (mass2_summed | units.g).as_quantity_in(units.MSun), label="Timo selfsummed")
    amuse_plot.plot((observed.outer_radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
        (mass2.cumsum() | units.g).as_quantity_in(units.MSun), label="Timo")
    pyplot.legend(loc=2)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    # pyplot.show()


def write_parms_to_textfile(filename, parms):
    """ Write Toycluster parms (dict) to plain text file """

    with open(filename, "w") as f:
        str=""
        for k, v in parms.iteritems():
            str += "{0} {1}\n".format(k, v)

        # print "Toycluster parameters"
        # print str

        f.write(str)

def plot_toycluster_rho0(observed, analytical):
    parms = OrderedDict({
        ('Output_file', './IC_single_0'),
        ('Ntotal', 2000000),
        ('Mtotal', 1000000),
        ('Mass_Ratio', 0),
        ('ImpactParam', 0),
        ('ZeroEOrbitFrac', 0.1),
        ('Cuspy', 1),
        ('Redshift', 0.0562),
        ('Bfld_Norm', 0),
        ('Bfld_Eta', 0),
        ('Bfld_Scale', 0),
        ('bf', 0.17),
        ('h_100', 0.7),
        ('UnitLength_in_cm', "3.085678e+21"),
        ('UnitMass_in_g', "1.989e+43"),
        ('UnitVelocity_in_cm_per_s', 100000),
        # ('c_nfw_0', 3.0),
        # ('v_com_0', 0.0),
        # ('rc_0', 28.06),
        # ('c_nfw_1', 3.0),
        # ('v_com_1', -1100),
        # ('rc_1', 276.26)
    })

    chandra_rho0 = analytical.rho0.value_in(units.g/units.cm**3)
    toycluster_rho0 = []
    mass = range(100, 1000000)
    N = len(mass)
    for i, m in enumerate(mass):
        parms['Mtotal'] = m
        write_parms_to_textfile("ToyclusterTrial.par", parms)
        os.system("./NoParticles.sh")
        tc = Toycluster2RuntimeOutputParser("NoParticles.txt")
        toycluster_rho0.append((tc.halosetup[0]['rho0gas_cgs']).value_in(units.g/units.cm**3))
        print_progressbar(i, N)

        if toycluster_rho0[i]/chandra_rho0 < 0.1:
            break

    pyplot.figure()
    pyplot.plot(mass, toycluster_rho0)
    pyplot.axhline(chandra_rho0)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.show()

    print "goal (chandra):", chandra_rho0
    print "Chandra    :", chandra_rho0
    print "Toycluster :", toycluster_rho0
    print "Ratio      :", toycluster_rho0/chandra_rho0


def get_cluster_mass_from_toycluster(observed, analytical, verbose=True):
    """ We have rho0 and r_c from Chandra data. We want Toycuster to give us
        the same rho0. For which mass is this? We brute force which mass
        to plug into Toycluster for a given r_c to match rho0.

        Note that in order to run fast Toycluster should exit after
        the setup only. (i.e. hack Toycluster executable!!). """
    parms = OrderedDict({
        ('Output_file', './IC_single_0'),
        ('Ntotal', 2000000),
        ('Mtotal', 100000),
        ('Mass_Ratio', 0),
        ('ImpactParam', 0),
        ('ZeroEOrbitFrac', 0.1),
        ('Cuspy', 1),
        ('Redshift', 0.0562),
        ('Bfld_Norm', 0),
        ('Bfld_Eta', 0),
        ('Bfld_Scale', 0),
        ('bf', 0.17),
        ('h_100', 0.7),
        ('UnitLength_in_cm', "3.085678e+21"),
        ('UnitMass_in_g', "1.989e+43"),
        ('UnitVelocity_in_cm_per_s', 100000),
        ('c_nfw_0', 3.0),
        ('v_com_0', 0.0),
        ('rc_0', 0.0),
        ('c_nfw_1', 3.0),
        ('v_com_1', 0),
        ('rc_1', 0.0)
    })

    parms['rc_0'] = analytical.rc.value_in(units.kpc)

    if observed.name == "cygA":
        parms['Cuspy'] = 1
        parms['Output_file'] = "CygA_trial"
    elif observed.name == "cygB":
        parms['Cuspy'] = 0
        parms['Output_file'] = "CygB_trial"
        # parms['c_nfw_0'] = 4

    # Bisection the mass. Here: initial guess 1e15. NB: code units = 1e10 MSun
    lower = 100         # 1e12
    upper = 1000000  # 1e16
    #upper = 100000000000  # 1e21

    # start loop
    chandra_rho0 = analytical.rho0.value_in(units.g/units.cm**3)
    ratio = 0
    epsilon = 1
    print "\nObtaining mass trough bisection method for", observed.name
    while upper-lower > epsilon:
        current = (lower+upper)/2.
        parms['Mtotal'] = current
        write_parms_to_textfile("ToyclusterTrial.par", parms)
        os.system("./NoParticles.sh")
        tc = Toycluster2RuntimeOutputParser("NoParticles.txt")
        toycluster_rho0 = (tc.halosetup[0]['rho0gas_cgs']).value_in(units.g/units.cm**3)
        ratio = toycluster_rho0/chandra_rho0

        print "Lower      :", lower, "10^10 Msun"
        print "Current    :", parms['Mtotal'], "10^10 Msun"
        print "Upper      :", upper, "10^10 Msun"
        print "Chandra    :", chandra_rho0, "g/cm**3"
        print "Toycluster :", toycluster_rho0, "g/cm**3"
        print "Ratio      :", ratio
        print "R200       :", tc.halosetup[0]['R200']
        print "Rc         :", tc.halosetup[0]['rc']
        print "c_nfw      :", tc.halosetup[0]['c_nfw']
        print

        if toycluster_rho0 > chandra_rho0:
            upper = current
        if toycluster_rho0 < chandra_rho0:
            lower = current

        # raw_input("Press enter to continue")

    return tc
    os.system("./YesParticles.sh")
    if observed.name == "cygA":
        os.system("mv YesParticles.txt RunToyClusterCygA.log")
    elif observed.name == "cygB":
        os.system("mv YesParticles.txt RunToyClusterCygB.log")
    return tc


if __name__ == "__main__":
    print "Reading Chandra Observed Density Profiles"
    print 80*"-"
    cygA_observed = ObservedCluster("cygA")
    cygB_observed = ObservedCluster("cygB")
    print cygA_observed
    print cygB_observed
    print 80*"-"

    fit = True
    if fit:
        print "Obtaining parameters from observed density profiles"
        print 80*"-"
        discard_firstbins = False
        if discard_firstbins:
            cygA_fit = fit_betamodel_to_chandra(cygA_observed, parm=[0.135, 27, 1.])
            cygA_observed.radius = cygA_observed.radius[4:]
            cygA_observed.binsize = cygA_observed.binsize[4:]
            cygA_observed.density= cygA_observed.density[4:]
            cygA_observed.density_std = cygA_observed.density_std[4:]
        cygA_fit = fit_betamodel_to_chandra(cygA_observed, parm=[0.135, 27, 1.])
        cygB_fit = fit_betamodel_to_chandra(cygB_observed, parm=[1., 1., 1.])
        cygB_fit = fit_betamodel_to_chandra(cygB_observed, parm=[0.002, 200, 700])
        cygB_fit = fit_betamodel_to_chandra(cygB_observed, parm=[0.002, 1.])

        # Giving radius --> discrete (observed) radius; not 'continuous'
        cygA_analytical = AnalyticalCluster(cygA_fit["x"], None, cygA_observed.radius)
        cygB_analytical = AnalyticalCluster(cygB_fit["x"], None, cygB_observed.radius)
        print "CygA n_e0  = {0:1.4e} 1/cm**3".format(
            cygA_analytical.ne0.value_in(1/units.cm**3))
        print "CygA rho_0 = {0:1.4e} g/cm**3".format(
            cygA_analytical.rho0.value_in(units.g/units.cm**3))
        print "CygB n_e0  = {0:1.4e} 1/cm**3".format(
            cygB_analytical.ne0.value_in(1/units.cm**3))
        print "CygB rho_0 = {0:1.4e} g/cm**3".format(
            cygB_analytical.rho0.value_in(units.g/units.cm**3))
        print 80*"-"

    get_total_mass = False
    if get_total_mass:
        print "Obtaining total mass and mass ratio"
        print 80*"-"
        print cygA_analytical.rc.value_in(units.kpc)
        print cygA_analytical.rho0.value_in(units.g/units.cm**3)

        rhocrit200 = 200*cygA_observed.cc.rho_crit()

        bf = 0.17

        # TODO

        print 80*"-"


    fit_toycluster = False
    if fit_toycluster:
        print "Fitting Toycluster calculated rho_0, rc to Chandra data"
        print 80*"-"
        tcA = get_cluster_mass_from_toycluster(cygA_observed, cygA_analytical)
        tcB = get_cluster_mass_from_toycluster(cygB_observed, cygB_analytical)

        cygA_mass = tcA.halosetup[0]['Mass_in_R200']
        cygB_mass = tcB.halosetup[0]['Mass_in_R200']

        print "M_CygA          =", cygA_mass
        print "rc              =", tcA.halosetup[0]['rc']
        print "c_nfw           =", tcA.halosetup[0]['c_nfw']
        print
        print "M_CygB          =", cygB_mass
        print "rc              =", tcB.halosetup[0]['rc']
        print "c_nfw           =", tcB.halosetup[0]['c_nfw']

        print
        print "M_CygA/M_CygB   =", cygA_mass/cygB_mass
        print "M_CygB/M_CygA   =", cygB_mass/cygA_mass
        print "Mtotal          =", cygA_mass + cygB_mass

        print 80*"-"

    plot_fit = True
    if plot_fit:
        # Number density
        # plot_fit_results(cygA_observed, cygA_analytical,
        #                  numerical=None, mass_density=False)
        # plot_fit_results(cygB_observed, cygB_analytical,
        #                  numerical=None, mass_density=False)

        # Mass density, including numerical best-fit model
        # cygA_numerical = NumericalCluster(
        #     icdir="./",
        #     snapdir="./",
        #     logfile="RunToyClusterCygA.log",
        #     icfile="CygA_trial")
        # cygB_numerical = NumericalCluster(
        #     icdir="./",
        #     snapdir="./",
        #     logfile="RunToyClusterCygB.log",
        #     icfile="CygB_trial")

        # Works, but residuals need some attention...
        plot_fit_results(cygA_observed, cygA_analytical,
                         mass_density=True)
        pyplot.show()
        plot_fit_results(cygB_observed, cygB_analytical,
                         mass_density=True)
        pyplot.show()

    obtain_mass = False
    if obtain_mass:
        print "Obtaining observed mass profiles"
        print 80*"-"
        cygA_mass = get_mass_profile(cygA_observed, cygA_fit)
        cygB_mass = get_mass_profile(cygB_observed, cygA_fit)
        print 80*"-"

    find_r200_observed = False
    if find_r200_observed:
        print "Obtaining r_200 from Chandra observation"
        print 80*"-"
        rhocrit200A = 200*cygA_observed.cc.rho_crit()
        rhocrit200B = 200*cygB_observed.cc.rho_crit()

        mass_densityA = cygA_observed.density
        volumeA = 4*numpy.pi*cygA_observed.radius**2 * cygA_observed.binsize
        # cm**-3 --> kpc**-3 = 2.9379989455e+64 ?
        rho_averageA = (cygA_observed.bin_volume/2.9379989455e+64)*cygA_observed.density / volumeA
        pyplot.figure(figsize=(12,9))
        pyplot.plot(cygA_observed.radius, cygA_observed.density, label="obs")
        pyplot.plot(cygA_observed.radius, rho_averageA, label="avg")
        pyplot.axhline(rhocrit200A)
        pyplot.gca().set_xscale("log")
        pyplot.gca().set_yscale("log")
        pyplot.legend()
        pyplot.show()
        print rho_averageA
        r_200A = cygA_observed.radius[(numpy.abs(rho_averageA-rhocrit200A)).argmin()]

        mass_densityB = cygB_observed.density
        rho_averageB = (cygB_observed.density / (4./3*numpy.pi*cygB_observed.radius**3))
        r_200B = cygB_observed.radius[(numpy.abs(rho_averageB-rhocrit200B)).argmin()]

        print "CygA"
        print "  200*rhocrit =", rhocrit200A
        print "  Indices of rho_average where rho_average > 200*rhocrit"
        print "   ", numpy.where(rho_averageA > rhocrit200A)
        print "  r200        =", r_200A
        print
        print "CygB"
        print "  200*rhocrit =", rhocrit200B
        print "  Indices of rho_average where rho_average > 200*rhocrit"
        print "   ", numpy.where(rho_averageB > rhocrit200B)
        print "  r200        =", r_200B
        print 80*"-"
        import sys; sys.exit(0)

    mass_analytical = False
    if mass_analytical:
        print "Obtaining cluster mass from analytical density profile"
        print 80*"-"
        get_cluster_mass_analytical(cygA_observed, cygA_fit)
        get_cluster_mass_analytical(cygB_observed, cygB_fit)
        print 80*"-"

    pyplot.show()
