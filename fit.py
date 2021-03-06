import os
import numpy
import pandas
import csv
import scipy
from scipy import stats
from collections import OrderedDict

import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
matplotlib.rc("text", usetex=True)
pyplot.rcParams.update({"font.size": 42})
pyplot.rcParams.update({"xtick.major.size": 8})
pyplot.rcParams.update({"xtick.minor.size": 4})
pyplot.rcParams.update({"ytick.major.size": 8})
pyplot.rcParams.update({"ytick.minor.size": 4})
pyplot.rcParams.update({"xtick.major.width": 4})
pyplot.rcParams.update({"xtick.minor.width": 2})
pyplot.rcParams.update({"ytick.major.width": 4})
pyplot.rcParams.update({"ytick.minor.width": 2})
pyplot.rcParams.update({"xtick.major.pad": 16})
pyplot.rcParams.update({"xtick.minor.pad": 16})
pyplot.rcParams.update({"ytick.major.pad": 16})
pyplot.rcParams.update({"ytick.minor.pad": 16})
pyplot.rcParams.update({"legend.loc": "best"})
pyplot.rcParams.update({"figure.autolayout": True})


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


def Gas_density_profile(parm, r, free_beta=False, fix_cygA=False):
    """ Beta model if len(parm)==2 (Donnert 2014);
        Cut-off beta profile at rc and rcut if len(parm)==3 (Donnert et al. 2016, in prep)
        Double cut-off beta profile if len(parm)==5. See Toycluster/setup.c """
    DOUBLE_BETA_COOL_CORES = False
    rho0 = parm[0]
    rc = parm[1]
    beta = 2.0/3

    if fix_cygA:
        rho0 = 0.09162
        rc = 26.01621

    if len(parm) >= 3 and free_beta:
        beta = parm[2]

    rho = rho0 * (1 + p2(r/rc))**(-3.0/2*beta)   # not cut-off (Donnert 2014)
    if len(parm) >= 3 and not free_beta:
        rcut = parm[2]
        rho /= (1 + p3(r/rcut) * (r/rcut))

    if len(parm) == 5:  # Double Beta model (for Coolcore clusters)
        Rho0_Fac = parm[3]
        Rc_Fac = parm[4]

        rho0_cc = rho0 * Rho0_Fac
        rc_cc = rc / Rc_Fac

        rho += rho0_cc / (1 + p2(r/rc_cc)) / (1 + p3(r/rcut) * (r/rcut))

    return rho


def Gas_density_profile_wrapper(r, parm0, parm1, parm2=None, parm3=None,
        parm4=None, free_beta=False, fix_cygA=False):
    """ Same same, but scipy.optimize.minimize expects different form of
        function than scipy.optimize.curve_fit does. """
    if parm3 and parm4:  # coolcore, double beta
        return Gas_density_profile((parm0, parm1, parm2, parm3, parm4), r)
    elif parm2 and not free_beta:
        return Gas_density_profile((parm0, parm1, parm2), r)
    elif parm2 and free_beta:
        return Gas_density_profile((parm0, parm1, parm2), r, free_beta=True)
    else:
        return Gas_density_profile((parm0, parm1), r, fix_cygA=fix_cygA)


# Define the statistical model, in this case we shall use a chi-squared distribution, assuming normality in the errors
def stat(parm, x, y, dy, free_beta=False, fix_cygA=False):
    ymodel = Gas_density_profile(parm, x, free_beta=free_beta, fix_cygA=fix_cygA)
    chisq = numpy.sum((y - ymodel)**2 / dy**2)
    return(chisq)


def fit_betamodel_to_chandra(observed, parm=[1., 1., 1.],
                             free_beta=False, fix_cygA=False):
    """ Fit betamodel to cluster observation
        if len(parm) == 2: single beta-model with beta = 2/3
           len(parm) == 3: cut-off beta-model with beta = 2/3
           len(parm) == 5: cut-off double beta-model with beta = 2/3
        if free_beta: use single beta-model with beta as fit parameter
        if fix_cygA : "run" fit with beta, rc, and rho0 fixed for compatibility
                      with routines that import this function externally.

        returns scipy.curve_fit fit result
    """

    # Fit to data
    if observed.name == "cygA":
        if len(parm) == 2:
            if fix_cygA:
                bounds = [(0.091, 0.09162), (26.0160, 26.0162)]
            else:
                bounds = [(None, None), (None, None)]
        elif len(parm) == 3 and not free_beta:
            bounds = [(None, None), (None, None), (1700, 1900)]
        elif len(parm) == 3 and free_beta:
            bounds = [(None, None), (None, None), (0.0, 1.0)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1800), (None, None), (None, None)]
    elif observed.name == "cygB":
        if len(parm) == 2:
            bounds = [(None, None), (None, None)]
        if len(parm) == 3 and not free_beta:
            bounds = [(None, None), (None, None), (1200, 1350)]
        if len(parm) == 3 and free_beta:
            bounds = [(0.002, 0.003), (None, None), (0.0, 1.0)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1800), (None, None), (None, None)]

    result = scipy.optimize.minimize(stat, parm,
            args=(observed.radius, observed.number_density,
                  observed.number_density_std, free_beta, fix_cygA),
            method='L-BFGS-B', bounds=bounds)

    # Obtain and print MLEs
    ml_vals, ml_covar = obtain_mles(observed, result,
        free_beta=free_beta, fix_cygA=fix_cygA)

    return result, ml_vals, ml_covar


def obtain_mles(observed, result, free_beta=False, fix_cygA=False):
    ml_vals = result["x"]
    ml_func = result["fun"]

    moddof = len(ml_vals)  # Model degrees of freedom; nr of fit parameters
    dof = len(observed.number_density) - moddof  # degrees of freedom

    ch = scipy.stats.chi2(dof)
    pval = 1.0 - ch.cdf(ml_func)

    if len(ml_vals) == 2:
        model = 0
    elif len(ml_vals) == 3 and not free_beta:
        model = 1
    elif len(ml_vals) == 3 and free_beta:
        model = 3
    elif len(ml_vals) == 5:
        model = 2

    modelnames = {0: "beta-model",
                  1: "cut-off beta-model",
                  2: "cut-off double beta-model",
                  3: "free beta-model"}

    print observed.name
    print "Results for the '{0}' model:".format(modelnames[model])
    print "  Using scipy.optimize.minimize to minimize chi^2 yields:"
    print "    n_e,0       = {0:.5f}".format(ml_vals[0])
    print "    r_c         = {0:.5f}".format(ml_vals[1])
    if len(ml_vals) == 3 and not free_beta:
        print "    r_cut       = {0:.5f}".format(ml_vals[2])
    if len(ml_vals) == 3 and free_beta:
        print "    beta        = {0:.5f}".format(ml_vals[2])
    print "    chisq       = {0:.5f}".format(ml_func)
    print "    dof         = {0:.5f}".format(dof)
    print "    chisq/dof   = {0:.5f}".format(ml_func/dof)
    print "    p-value     = {0:.5f}".format(pval)

    if free_beta:
        ml_vals, ml_covar = scipy.optimize.curve_fit(
                lambda r, parm0, parm1, parm2: Gas_density_profile_wrapper(r,
                parm0, parm1, parm2, free_beta=free_beta, fix_cygA=fix_cygA),
                observed.radius, observed.number_density, p0=ml_vals,
                sigma=observed.number_density_std)
    elif fix_cygA:
        ml_vals, ml_covar = scipy.optimize.curve_fit(
                lambda r, parm0, parm1: Gas_density_profile_wrapper(r,
                parm0, parm1, free_beta=free_beta, fix_cygA=fix_cygA),
                observed.radius, observed.number_density, p0=ml_vals,
                sigma=observed.number_density_std)
    else:
        ml_vals, ml_covar = scipy.optimize.curve_fit(
                Gas_density_profile_wrapper,
                observed.radius, observed.number_density, p0=ml_vals,
                sigma=observed.number_density_std)
        # ml_funcval = stat(ml_vals, edges, dens, err, model)

    if not result["success"]:
        print "  scipy.optimize.curve_fit broke down!\n    Reason: '{0}'"\
            .format(result["message"])
        print "  No confidence intervals have been calculated."

    err = numpy.sqrt(numpy.diag(ml_covar))
    print "  Using scipy.optimize.curve_fit to obtain confidence intervals yields:"
    print "    n_e,0       = {0:.5f} +/- {1:.5f}".format(ml_vals[0], err[0])
    print "    r_c         = {0:.5f} +/- {1:.5f}".format(ml_vals[1], err[1])
    if len(ml_vals) == 3 and not free_beta:
        print "    r_cut       = {0:.5f} +/- {1:.5f}".format(ml_vals[2], err[2])
    if len(ml_vals) == 3 and free_beta:
        print "    beta        = {0:.5f} +/- {1:.5f}".format(ml_vals[2], err[2])
    print

    return ml_vals, ml_covar


def plot_fit_results(observed, analytical, numerical=None,
                     mass_density=False, save=False, poster_style=False,
                     fix_cygA=False):
    if poster_style:
        pyplot.style.use(["dark_background"])
        pyplot.rcParams.update({"font.weight": "bold"})
        # magenta, dark blue, orange, green, light blue (?)
        data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                       (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                       (41./255, 239./255, 239./255)]
        fit_colour = "white"
    else:
        data_colour = ["g", "r", "b"]
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
    # pyplot.title(observed.name)
    pyplot.errorbar(observed.radius+observed.binsize/2, observed_density,
            xerr=observed.binsize/2, yerr=observed_density_std,
            marker='o', ms=5 if poster_style else 3,
            elinewidth=3 if poster_style else 1,
            ls="", c=data_colour[0] if observed.name=="cygA" else data_colour[2],
            label="800" if observed.oldICs else "982"+\
                  " ks Chandra\n(Wise et al. in prep)",)

    label = "model: "+analytical.modelname+"\n\t"
    label +=  r"$n_{{e,0}} \,$ = {0:.2e} cm$^{{-3}}$".format(analytical.ne0.number) if not mass_density else r"$\rho_{{0}} \,$ = {0:.2e} g cm$^{{-3}}$".format(analytical.rho0.number)
    label += "\n\t"+ r"$r_c \,$ = {0:.2f} kpc".format(analytical.rc.number)
    if analytical.rcut is not None and not analytical.free_beta:
        label += "\n\t"+r"$r_{{\rm cut}}$ = {0:.2f} kpc".format(analytical.rcut.number)
    if analytical.free_beta:
        label += "\n\t"+r"$\beta$ = {0:.2f}".format(analytical.beta)
    if analytical.rho0_cc is not None:
        label += "\n\t"+r"$n_{{e,0,cc}}$ = {0:.2f} cm$^{{-3}}$".format(analytical.ne0_cc.number) if not mass_density else "\n\t"+r"$\rho_{{0,cc}}$ = {0:.2f} kpc".format(analytical.rho0_cc.number)
        label += "\n\t"+r"$r_{{c,cc}}$ = {0:.2f} kpc".format(analytical.rc_cc.number)

    label = r"\begin{tabular}{lll}"
    if analytical.free_beta:
        label += " model & = & free beta \\\\"
    else:
        label += " model & = & fixed beta \\\\"
    label += r" rho0 & = & {0:.2e} g$\cdot$cm$^{{-3}}$ \\".format(analytical.rho0.number)
    label += " rc & = & {0:.2f} kpc \\\\".format(analytical.rc.number)
    if analytical.free_beta:
        label += " beta & = & {0:.3f} kpc \\\\".format(analytical.beta)
    label += (" \end{tabular}")


    amuse_plot.plot(analytical.radius, analytical_density,
            c=fit_colour, lw=3 if poster_style else 1, label=label)

    if numerical:
        pyplot.scatter(numerical.gas.r.value_in(units.kpc),
            numerical.gas.rho.value_in(units.g / (units.cm**3)),
            c="g", edgecolor="face", s=4, label=r"Gas, numerical")

    ax.set_xscale("log")
    ax.set_yscale("log")
    if mass_density:
        ax.set_ylabel(r"Gas density [g/cm$**$3]", fontsize=38)
        ax.set_ylim(min(observed_density)/1.5, max(observed_density)*1.3)
    else:
        ax.set_ylabel(r"Gas density [1/cm$**$3]", fontsize=38)
    ax.legend(loc=3, prop={'size':30})

    # Plot Residuals
    pyplot.sca(ax_r)
    residual_density = (observed_density - analytical_density)/observed_density
    pyplot.errorbar(observed.radius+observed.binsize/2, 100*residual_density,
            yerr=100*observed_density_std/observed_density, c=fit_colour,
            lw=3 if poster_style else 1,
            elinewidth=1 if poster_style else 1, drawstyle="steps-mid")
    ax_r.axhline(y=0, lw=3 if poster_style else 1, ls="dashed", c=fit_colour)

    ax_r.set_xscale("log")
    # ax_r.set_yscale("log")

    ax.set_xlim(min(observed.radius)-0.3, max(observed.radius)+2000)
    ax_r.set_xlim(min(observed.radius)-0.3, max(observed.radius)+2000)

    # Show 50% deviations from the data
    if fix_cygA:
        ax_r.set_ylim(-20, 100)
    else:
        ax_r.set_ylim(-50, 50)

    # Fix for overlapping y-axis markers
    from matplotlib.ticker import MaxNLocator
    ax.tick_params(labelbottom='off')
    nbins = len(ax_r.get_yticklabels())
    ax_r.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

    ax_r.set_xlabel(r"Radius [kpc]", fontsize=38)
    ax_r.set_ylabel(r"Residuals [\%]", fontsize=38)

    ax.axvline(x=analytical.rc.value_in(units.kpc),
               lw=3 if poster_style else 1, ls="dashed", c=fit_colour)
    ax_r.axvline(x=analytical.rc.value_in(units.kpc),
                 lw=3 if poster_style else 1, ls="dashed", c=fit_colour)

    # Force axis labels to align
    ax.get_yaxis().set_label_coords(-0.15, 0.5)
    ax_r.get_yaxis().set_label_coords(-0.15, 0.5)

    pyplot.sca(ax)
    ylim = pyplot.gca().get_ylim()
    pyplot.text(1.05*analytical.rc.value_in(units.kpc), 0.95*ylim[1], r"$r_c$",
                ha="left", va="top", color="k", fontsize=38)

    if save and not analytical.free_beta:
        pyplot.savefig("out/density_betamodel_fit_{0}{1}{2}.png"\
            .format(observed.name, "_800ksec" if observed.oldICs else "_900ksec",
                "_dark" if poster_style else ""), dpi=300)
    if save and analytical.free_beta:
        pyplot.savefig("out/density_free_betamodel_fit_{0}{1}{2}.png"\
            .format(observed.name, "_800ksec" if observed.oldICs else "_900ksec",
                "_dark" if poster_style else ""), dpi=300)
    if not save:
        pyplot.show()


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
    # We obtain core radius and central density of CygA from fit where
    # beta is also a free parameter. But then we want Toycuster to reach
    # the same central density for which we need to continue the analysis
    # to infer the total mass and mass ratio given the model that only
    # fits the inner region, but has a systematic deviation from the data
    # because the observed slope is steeper but we are bound to beta 2/3
    # in Toycluster... for now. Hence, fix_cygA

    fix_cygA = False

    fit = True
    discard_firstbins = True
    discard_lastbins = False
    plot_fit = True
    poster_style = False
    mass_density = True
    save = True
    fit_toycluster = False



    print "Reading Chandra Observed Density Profiles"
    print 80*"-"
    cygA_observed = ObservedCluster("cygA", oldICs=False)
    cygB_observed = ObservedCluster("cygB", oldICs=False)
    print cygA_observed
    print cygB_observed
    print 80*"-"

    if fit:
        print "Obtaining parameters from observed density profiles"
        print 80*"-"
        # Due to pile up in inner region (?). Also inside kernel: cannot model
        # in a stable way
        if discard_firstbins:
            cygA_observed.radius = cygA_observed.radius[3:]
            cygA_observed.binsize = cygA_observed.binsize[3:]
            cygA_observed.density= cygA_observed.density[3:]
            cygA_observed.density_std = cygA_observed.density_std[3:]
            cygA_observed.number_density = cygA_observed.number_density[3:]
            cygA_observed.number_density_std = cygA_observed.number_density_std[3:]

        if discard_lastbins:
            cut = 40
            cygA_observed.radius = cygA_observed.radius[:-cut]
            cygA_observed.binsize = cygA_observed.binsize[:-cut]
            cygA_observed.density= cygA_observed.density[:-cut]
            cygA_observed.density_std = cygA_observed.density_std[:-cut]
            cygA_observed.number_density = cygA_observed.number_density[:-cut]
            cygA_observed.number_density_std = cygA_observed.number_density_std[:-cut]

        # Cut-off single beta (2/3)
        cygA_fit, cygA_ml_vals, cygA_ml_covar = \
            fit_betamodel_to_chandra(cygA_observed, parm=[0.135, 27, 1755])
        # Single beta (beta is fit parameter)
        cygA_fit_free, cygA_free_ml_vals, cygA_free_ml_covar = \
            fit_betamodel_to_chandra(cygA_observed,
                parm=[0.1, 10, 0.67], free_beta=True)
        # Single beta (2/3): use this for numerical setup
        cygA_fit, cygA_ml_vals, cygA_ml_covar = fit_betamodel_to_chandra(
            cygA_observed, parm=[0.1, 10], fix_cygA=fix_cygA)

        # Cut-off single beta (2/3)
        cygB_fit, cygB_ml_vals, cygB_ml_covar = \
            fit_betamodel_to_chandra(cygB_observed, parm=[0.001, 264.65, 1313.82])
        # Single beta (beta is fit parameter)
        cygB_fit_free, cygB_free_ml_vals, cygB_free_ml_covar = \
            fit_betamodel_to_chandra(cygB_observed,
                parm=[0.001, 100, 0.79], free_beta=True)
        # Single beta (2/3): use this for numerical setup
        cygB_fit, cygB_ml_vals, cygB_ml_covar = \
            fit_betamodel_to_chandra(cygB_observed, parm=[0.001, 100])

        # Giving radius --> discrete (observed) radius; not 'continuous'
        cygA_analytical = AnalyticalCluster(cygA_fit["x"], None, cygA_observed.radius)
        cygB_analytical = AnalyticalCluster(cygB_fit["x"], None, cygB_observed.radius)

        cygA_analytical_free = AnalyticalCluster(cygA_fit_free["x"], None, cygA_observed.radius, free_beta=True)
        cygB_analytical_free = AnalyticalCluster(cygB_fit_free["x"], None, cygB_observed.radius, free_beta=True)
        print "CygA n_e0  = {0:1.4e} 1/cm**3".format(
            cygA_analytical.ne0.value_in(1/units.cm**3))
        print "CygA rho_0 = {0:1.4e} g/cm**3".format(
            cygA_analytical.rho0.value_in(units.g/units.cm**3))
        print "CygB n_e0  = {0:1.4e} 1/cm**3".format(
            cygB_analytical.ne0.value_in(1/units.cm**3))
        print "CygB rho_0 = {0:1.4e} g/cm**3".format(
            cygB_analytical.rho0.value_in(units.g/units.cm**3))
        print 80*"-"

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

        plot_fit_results(cygA_observed, cygA_analytical,
            mass_density=mass_density, save=save, poster_style=poster_style)
        plot_fit_results(cygB_observed, cygB_analytical,
            mass_density=mass_density, save=save, poster_style=poster_style)

        plot_fit_results(cygA_observed, cygA_analytical_free,
            mass_density=mass_density, save=save, poster_style=poster_style)
        plot_fit_results(cygB_observed, cygB_analytical_free,
            mass_density=mass_density, save=save, poster_style=poster_style)

    #pyplot.show()

    #import sys; sys.exit(0)

    # None of the code below actually works (physically). This is a neat
    # toycluster bruteforce which we leave in cause it is kind of awesome tho.
    # See get_mass_plus_plot.py for a working (physical) implementation.
    if fit_toycluster:
        print "Fitting Toycluster calculated rho_0, rc to Chandra data"
        print 80*"-"
        tcA = get_cluster_mass_from_toycluster(cygA_observed, cygA_analytical)

        cygA_mass = tcA.halosetup[0]['Mass_in_R200']

        print "M_CygA          =", cygA_mass
        print "rc              =", tcA.halosetup[0]['rc']
        print "c_nfw           =", tcA.halosetup[0]['c_nfw']
