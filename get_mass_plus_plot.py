"""
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Mon May 16, 2016 05:23 pm

Script to obtain M200 for CygA, CygB, thus xM (mass ratio),
and cNFW for cygA and cygB using reverse of Toycluster setup.c.
Code based on Julius Donnert's 20160617 cyg.pro script.

"""

import numpy
from scipy import special
import pandas
import matplotlib
# Anaconda python gives annoying "setCanCycle: is deprecated" when using Tk
# matplotlib.use("TkAgg")
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

from cluster import ObservedCluster
from cluster import AnalyticalCluster
from fit import fit_betamodel_to_chandra
from macro import print_progressbar


# from amuse.units import constants
# from amuse.units import units
# g2msun = units.g(1).value_in(units.MSun)
# msun2g = units.MSun(1).value_in(units.g)
# kpc2cm = units.kpc(1).value_in(units.cm)
# cm2kpc = units.cm(1).value_in(units.kpc)

g2msun = 5.02785e-34
msun2g = 1.98892e+33
kpc2cm = 3.08568e+21
cm2kpc = 3.24078e-22

# Macro's
def p2(a):
    return ((a)*(a))
def p3(a):
    return ((a)*(a)*(a))


def ne_to_rho(ne):
    """ convert electron number density (1/cm**3) to mass density (g/cm**3) """
    # Hydrogen fraction
    xH = 0.76              # Gadget-2 uses 0.76 in allvars.h
    # Mean molecular weight
    umu = 4./(5.*xH+3.)    # Gadget-2 read_ic.c:129, assuming full ionisation (T>1e4)

    # Proton mass in g
    mp = 1.672621637e-24
    rho = umu*ne*mp
    return rho


def rho_to_ne(rho):
    """ convert mass density to electron number density """
    # Hydrogen fraction
    xH = 0.76              # Gadget-2 uses 0.76 in allvars.h
    # Mean molecular weight
    umu = 4./(5.*xH+3.)    # Gadget-2 read_ic.c:129, assuming full ionisation (T>1e4)

    # Proton mass in g
    mp = 1.672621637e-24
    ne = rho/(umu*mp)
    return ne


def concentration_parameter(M200):
    """
    Calculate concentration parameter (Duffy+ 2008).
    NB here assuming H0 = 70 km/s/Mpc.

    returns the concentration parameter cNFW
    """

    M200 *= g2msun  # because M200 should be in Msun here
    A = 5.74
    B = -0.097
    C = -0.47
    Mpivot = 2e12 / 0.7
    cNFW = A * numpy.power( M200/Mpivot, B)  # * numpy.power( 1+z, C)
    return cNFW


def hernquist_a(rs, cNFW):
    """ Calculate the Hernquist scalelength a given r_s and the concentration parameter cNFW """

    a = rs * numpy.sqrt(2*numpy.log(1+cNFW) - cNFW/(1+cNFW))
    return a


def dm_density_hernquist(r, Mdm, a):
    """ Hernquist model for DM density profile (Hernquist 1990).
    NB this is similar to the NFW profile but with a steeper cut-off
    at radii > 0.1 R200 (?)"""

    rho_dm_hernquist = Mdm/(2*numpy.pi) * a / (r*p3(r + a))
    return rho_dm_hernquist


def M_dm_below_r(r, Mdm, a):
    """ Hernquist (1990) dark matter mass profile """

    M_dm_below_r = Mdm * p2(r) / p2(r + a)
    return M_dm_below_r


def gas_density_beta(r, rho0, rc, beta=None):
    """ Beta-model for gas density profile (Cavaliere & Fusco-Femiano 1978)
    with a fixed value beta = 2/3 (Mastropietro & Burkert 2008)

    Donnert (2014): "Should any other value of beta be used, this step
                     would involve the confluent hypergeometric function,
                     which usually makes the analytical solution of the
                     hydrostatic equation impossible. Another elegant
                     solution would be beta=1, which however is unfavoured
                     from observations"
    """

    if not beta:
        beta = 2./3
    rho_gas_beta = rho0 * (1 + p2(r)/p2(rc))**(-3*beta/2.)
    return rho_gas_beta


def M_gas_below_r(r, rho0, rc, beta=None):
    """ M(<r) for Hydrostatic Equilibrium, Spherical Symmetry
    Beta-Model where beta is fixed to 2/3 """

    if not beta:
        M_gas_below_r = 4*numpy.pi*p3(rc)*rho0 * (r/rc - numpy.arctan(r/rc))
    else:
        # well, isn't she lovely... scipy built-in Gauss Hypergeometric function
        M_gas_below_r = special.hyp2f1(1.5, 1.5*beta, 2.5, -p2(r/rc))
        M_gas_below_r *= 4*numpy.pi*rho0*p3(r)/3
    return M_gas_below_r


def H(z=0, H0=70, WM=0.3, WV=0.7):
    """ Hubble constant as a function of redshift """

    return H0 * numpy.sqrt((WM*p3(1+z) + WV))


def rho_crit(z=0.0562):
    """ Critical density of the Universe as a function of redshift """

    G = 6.67428e-8  # cm**3 g**-1 s**-2
    Hz = H(z) * 3.24077928966e-20 # km/s/Mpc --> 1/s
    rho_crit = 3 * p2(Hz) / (8 * numpy.pi * G)
    return rho_crit


def plot_observed_cluster(observed, analytical_density, poster_style=False):
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

    data_colour = data_colour[0] if observed.name == "cygA" else data_colour[2]

    fig, (ax, ax_r) = pyplot.subplots(2, 2, sharex=True, figsize=(16, 12))
    gs1 = matplotlib.gridspec.GridSpec(3, 3)
    gs1.update(hspace=0)
    ax = pyplot.subplot(gs1[:-1,:])
    ax_r = pyplot.subplot(gs1[-1,:])  # residuals

    # Plot data
    pyplot.sca(ax)
    pyplot.errorbar(observed.radius+observed.binsize/2,
                    observed.density, xerr=observed.binsize/2,
                    yerr=observed.density_std, marker="o",
                    ms=5 if poster_style else 3, elinewidth=3 if poster_style else 1,
                    ls="", c=data_colour)
                    #label="800 ks Chandra\n(Wise+ 2016, in prep)")

    # Plot analytical gas profile
    pyplot.plot(observed.radius, analytical_density, lw=3 if poster_style else 1,
                c=fit_colour, ls="dashed")#, label="gas")

    # Plot Residuals
    pyplot.sca(ax_r)
    residual_density = (observed.density - analytical_density)/observed.density
    pyplot.errorbar(observed.radius+observed.binsize/2, 100*residual_density,
            yerr=100*observed.density_std/observed.density, c=fit_colour,
            lw=3 if poster_style else 1,
            elinewidth=1 if poster_style else 1, drawstyle="steps-mid")
    # Show dashed residuals zero line
    ax_r.axhline(y=0, lw=3 if poster_style else 1, ls="dashed", c=fit_colour)

    # Set axis labels
    ax.set_ylabel(r"Density [g cm$^{-3}$]")
    ax_r.set_xlabel(r"Radius [kpc]")
    ax_r.set_ylabel("Residuals [\%]")

    # Set logscale, but the residual y-axis is not logarithmic!
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax_r.set_xscale("log")

    # Set residual y-limits: show 50% deviations from the data
    ax_r.set_ylim(-50, 50)

    # Fix for overlapping y-axis markers
    from matplotlib.ticker import MaxNLocator
    ax.tick_params(labelbottom="off")
    nbins = len(ax_r.get_yticklabels())
    ax_r.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune="upper"))

    # Force axis labels to align
    ax.get_yaxis().set_label_coords(-0.1,0.5)
    ax_r.get_yaxis().set_label_coords(-0.1,0.5)

    return fig


def obtain_M200_bisection(rc, rho0, beta=None, verbose=False,
                          visualise=False, observed=None):
    """ We follow the Toycluster (Donnert 2014) setup.c method in reverse.
    If we assume a value for r200 and we assume the baryon fraction at r200
    is equal to 0.17 we are able to obtain a total mass M200 for which
    rho_average(r200) == 200 rho_crit."""

    if visualise:
        gas_rhom = gas_density_beta(observed.radius, rho0, rc*cm2kpc, beta)
        n=0

        if poster_style:
            pyplot.style.use(["dark_background"])
            # magenta, dark blue, orange, green, light blue (?)
            data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                           (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                           (41./255, 239./255, 239./255)]
            fit_colour = "white"
            accent_colour = (151./255, 24./255, 24./255)
        else:
            data_colour = ["g", "r", "b"]
            fit_colour = "k"
            accent_colour = "red"

    # Find r200 such that rho200 / rho_crit == 200
    lower = 10 * kpc2cm
    upper = 2000 * kpc2cm

    if visualise and observed.name == "cygB":
        lower = 100 * kpc2cm

    # bisection method
    epsilon = 0.0000001
    while upper/lower > 1+epsilon:
        # bisection
        r200 = (lower+upper)/2.

        bf = 0.17  # Mission critical assumption that bf = 0.17 at r200! Planelles+ 2013

        Mgas200 = M_gas_below_r(r200, rho0, rc, beta)
        Mdm200 = Mgas200 * (1/bf - 1)

        M200 = Mgas200 + Mdm200

        cNFW = concentration_parameter(M200)
        rs = r200 / cNFW
        a = hernquist_a(rs, cNFW)

        Mdm = Mdm200 / M_dm_below_r(r200, 1, a)

        """ Now rho_average(r200)/rhocrit should equal 200.
                If not? Try different r200"""
        rho200_over_rhocrit = ( M200 / (4./3 * numpy.pi * p3(r200))) / rho_crit()
        if verbose:
            print "Lower                  = {0:3.1f}".format(lower * cm2kpc)
            print "r200                   = {0:3.1f}".format(r200 * cm2kpc)
            print "Upper                  = {0:3.1f}".format(upper * cm2kpc)
            print "rho_avg(r200)/rho_crit = {0:.1f}".format(rho200_over_rhocrit)
            print "Ratio                  = {0:.1f}".format(rho200_over_rhocrit/200)
            print

        if visualise:
            fig = plot_observed_cluster(observed, gas_rhom)
            pyplot.figure(fig.number)
            ax, ax_r = fig.axes
            pyplot.sca(ax)

            dm_rhom = dm_density_hernquist(observed.radius*kpc2cm, Mdm, a)

            pyplot.plot(observed.radius, dm_rhom, c=fit_colour, lw=3 if poster_style else 1, ls="solid")
            pyplot.axhline(200*rho_crit(), c=fit_colour, lw=3 if poster_style else 1)
            rho_avg_200 = M200 / (4./3 * numpy.pi * p3(r200))
            pyplot.axhline(rho_avg_200, c=accent_colour, lw=3 if poster_style else 1)

            # Indicate bisection bounds
            pyplot.axvline(x=lower*cm2kpc, c=fit_colour, ls="dotted", lw=3 if poster_style else 1)
            pyplot.axvline(x=upper*cm2kpc, c=fit_colour, ls="dotted", lw=3 if poster_style else 1)
            pyplot.axvline(x=r200*cm2kpc, c=accent_colour, ls="solid", lw=3 if poster_style else 1)

            # Plot textbox with bisection info
            bisection_info = r"\begin{tabular}{lll}"
            bisection_info += " lower & : & {0:.1f} \\\\".format(lower*cm2kpc)
            bisection_info += " r200 & : & {0:.1f} \\\\".format(r200*cm2kpc)
            bisection_info += " upper & : & {0:.1f} \\\\".format(upper*cm2kpc)
            bisection_info += r"$\frac{{\rho(r200)}}{{\rho_{{\rm crit}}}}$ & : & {0:.1f}"\
                                .format(rho200_over_rhocrit)
            bisection_info += (" \end{tabular}")

            if observed.name == "cygA":
                textX = 3
                textY = 2e-29
            else:
                textX = 50
                textY = 2e-29

            pyplot.text(textX, textY, bisection_info, size=18,
                        ha="left", va="bottom", color=fit_colour,
                        bbox=dict(boxstyle="round",
                                  ec=accent_colour if poster_style else (1., 0.5, 0.5),
                                  fc=accent_colour if poster_style else (1., 0.8, 0.8),
                                  )
                       )

            # Set axis limits
            if observed.name == "cygA":
                ax_r.set_xlim(1, 3000)
                ax.set_xlim(1, 3000)
                ax.set_ylim(1e-29, 1e-22)
            else:
                ax_r.set_xlim(40, 5000)
                ax.set_xlim(40, 5000)
                ax.set_ylim(1e-29, 3e-25)

            pyplot.savefig("out/findmass_{1}_{0:03d}.png".format(n, observed.name))
            pyplot.close()
            n+=1
            # import sys; sys.exit(0)

        # bisection
        if rho200_over_rhocrit < 200:
            upper = r200
        if rho200_over_rhocrit > 200:
            lower = r200

    # r200, thus M200 found
    halo = dict()
    halo["r200"] = r200
    halo["rho200_over_rhocrit"] = rho200_over_rhocrit
    halo["rho0"] = rho0
    halo["ne0"] = rho_to_ne(rho0)
    halo["rc"] = rc
    halo["beta"] = beta
    halo["Mgas200"] = Mgas200
    halo["Mdm200"] = Mdm200
    halo["M200"] = M200
    halo["cNFW"] = cNFW
    halo["rs"] = rs
    halo["a"] = a
    halo["Mdm"] = Mdm

    return halo


def propagate_errors(halo, to_print=True):
    """ We have an error on rho0 and rc from the fit to the Chandra data

        We obtain M200 and the other parameters using
        rho0+rho0_sigma;rc+rc_sigma for +1 sigma parameters, and
        rho0-rho0_sigma;rc-rc_sigma for -1 sigma parameters.

        boolean to_print allows to switch between just the errors, or the
            errors plus(or min) the values.

        to_print False can be used to plot confidence intervals
        to_print True can be used to print just the errors
    """

    plus = obtain_M200_bisection(halo["rc"]+halo["rc_sigma"],
        halo["rho0"]+halo["rho0_sigma"],
        None if not halo["beta"] else halo["beta"]+halo["beta_sigma"],
        verbose=False, visualise=False)
    min = obtain_M200_bisection(halo["rc"]-halo["rc_sigma"],
        halo["rho0"]-halo["rho0_sigma"],
        None if not halo["beta"] else halo["beta"]-halo["beta_sigma"],
        verbose=False, visualise=False)

    if to_print:
        plus = pandas.Series(plus)
        halo = pandas.Series(halo)
        min = pandas.Series(min)

        plus1sigma = plus-halo
        min1sigma = min-halo

        return plus1sigma, min1sigma
    else:
        return plus, min


def print_inferred_values(halo):
    plus, min = propagate_errors(halo)

    bf_200 = halo["Mgas200"]/(halo["Mdm200"]+halo["Mgas200"])
    bf_200_plus = plus["Mgas200"]/(plus["Mdm200"]+plus["Mgas200"])
    bf_200_min = min["Mgas200"]/(min["Mdm200"]+min["Mgas200"])

    masses = ["Mgas200", "Mdm200", "M200", "Mdm"]
    radii = ["r200", "rc", "rs", "a"]

    print "{0:<20}   {1:<10} {2:<10} {3:<10}".format("quantity", "value", "+1sigma", "-1sigma")
    print "{0:<20} = {1:<10.3f} {2:<10.3f} {3:<10.3f}".format(
        "bf", bf_200, bf_200_plus-bf_200, bf_200_min-bf_200_min)
    for p in ["r200", "rho200_over_rhocrit", "rho0", "ne0", "rc", "beta",
            "Mgas200", "Mdm200", "M200", "cNFW", "rs", "a", "Mdm"]:
        if p in masses:
            print "{0:<20} = {1:<10.3e} {2:<10.3e} {3:<10.3e}".format(
                p, halo[p]*g2msun, plus[p]*g2msun, min[p]*g2msun)
        elif p in radii:
            print "{0:<20} = {1:<10.3f} {2:<10.3f} {3:<10.3f}".format(
                p, halo[p]*cm2kpc, plus[p]*cm2kpc, min[p]*cm2kpc)
        elif p == "rho0":
            print "{0:<20} = {1:<10.3e} {2:<10.3e} {3:<10.3e}".format(
                p, halo[p], plus[p], min[p])
        elif p == "ne0":
            print "{0:<20} = {1:<10.5f} {2:<10.5f} {3:<10.5f}".format(
                p, halo[p], plus[p], min[p])
        elif p == "beta" and not halo[p]:  # beta = 2/3
            print "{0:<20} = {1:<10.3f} {2:<10.3f} {3:<10.3f}".format(
                p, 2.0/3, 0, 0)
        else:
            print "{0:<20} = {1:<10.3f} {2:<10.3f} {3:<10.3f}".format(
                p, halo[p], plus[p], min[p])
    print


def make_plot(cygA, cygB, cygA_observed=None, cygB_observed=None, mode=""):
    """ Make plot of inferred profiles

    @param cyg*: dictionary with best-fit parameters (of gas and dm)
    @param mode: string to select plot type.
        NB If mode contains single, then only one cluster is plotted.
           To select which cluster: set one of cyg* dictionaries to None!

        Currently, options are:
        "ratio"       : create mass ratio plot as a function of radius
        "rhomassboth" : 2-panel subplot. left rho; right mass. A&B in same panel
        "massboth"    : 2-panel mass subplot. Left CygA, right CygB
        "rhoboth"     : 2-panel denistysubplot. Left CygA, right CygB
        "masssameplot": mass of both CygA and CygB in the same plot
        "rhosingle"   : plot Chandra gas and fit with residuals and DM bestfit
        "bfsingle"    : Baryon fraction as a function of radius (single cluster)
        "nfwsingle"   : Dark Matter density comparison of NFW and Hernquist

    """
    # Get continuous radius range. NB observed radius is discrete!
    r = numpy.arange(1, 1e4, 0.1)  # kpc!

    # Get gas and dark matter density and mass profiles
    # Density in g/cm**3, mass in MSun, radius in kpc. So units require care
    if cygA:
        cygA_gas_rhom = gas_density_beta(r, cygA["rho0"], cygA["rc"]*cm2kpc, cygA["beta"])
        cygA_gas_mass = M_gas_below_r(r, cygA["rho0"]*g2msun/p3(cm2kpc),
                cygA["rc"]*cm2kpc, cygA["beta"])
        cygA_dm_rhom = dm_density_hernquist(r*kpc2cm, cygA["Mdm"], cygA["a"])
        cygA_dm_mass = M_dm_below_r(r, cygA["Mdm"]*g2msun, cygA["a"]*cm2kpc)
    if cygB:
        cygB_gas_rhom = gas_density_beta(r, cygB["rho0"], cygB["rc"]*cm2kpc, cygB["beta"])
        cygB_gas_mass = M_gas_below_r(r, cygB["rho0"]*g2msun/p3(cm2kpc),
                cygB["rc"]*cm2kpc, cygB["beta"])

        cygB_dm_rhom = dm_density_hernquist(r*kpc2cm, cygB["Mdm"], cygB["a"])
        cygB_dm_mass = M_dm_below_r(r, cygB["Mdm"]*g2msun, cygB["a"]*cm2kpc)

    if mode == "ratio":
        print "Generating plot of the mass ratio"
        cygA_plus, cygA_min = propagate_errors(cygA, to_print=False)
        cygB_plus, cygB_min = propagate_errors(cygB, to_print=False)

        cygA_gas_plus_sigma = M_gas_below_r(r, cygA_plus["rho0"]*g2msun/p3(cm2kpc),
                cygA_plus["rc"]*cm2kpc, cygA_plus["beta"])
        cygA_dm_plus_sigma = M_dm_below_r(
            r, cygA_plus["Mdm"]*g2msun, cygA_plus["a"]*cm2kpc)
        cygB_gas_plus_sigma = M_gas_below_r(r, cygB_plus["rho0"]*g2msun/p3(cm2kpc),
                cygB_plus["rc"]*cm2kpc, cygB_plus["beta"])
        cygB_dm_plus_sigma = M_dm_below_r(
            r, cygB_plus["Mdm"]*g2msun, cygB_plus["a"]*cm2kpc)

        cygA_gas_min_sigma = M_gas_below_r(r, cygA_min["rho0"]*g2msun/p3(cm2kpc),
                cygA_min["rc"]*cm2kpc, cygA_min["beta"])
        cygA_dm_min_sigma = M_dm_below_r(
            r, cygA_min["Mdm"]*g2msun, cygA_min["a"]*cm2kpc)
        cygB_gas_min_sigma = M_gas_below_r(r, cygB_min["rho0"]*g2msun/p3(cm2kpc),
                cygB_min["rc"]*cm2kpc, cygB_min["beta"])
        cygB_dm_min_sigma = M_dm_below_r(
            r, cygB_min["Mdm"]*g2msun, cygB_min["a"]*cm2kpc)

        cygA_mass = cygA_gas_mass + cygA_dm_mass
        cygA_mass_plus_sigma = cygA_gas_plus_sigma + cygA_dm_plus_sigma
        cygA_mass_min_sigma = cygA_gas_min_sigma + cygA_dm_min_sigma
        cygB_mass = cygB_gas_mass + cygB_dm_mass
        cygB_mass_plus_sigma = cygB_gas_plus_sigma + cygB_dm_plus_sigma
        cygB_mass_min_sigma = cygB_gas_min_sigma + cygB_dm_min_sigma

        ratio = numpy.array(cygA_mass/cygB_mass)
        ratio_plus = numpy.array(cygA_mass_plus_sigma/cygB_mass_plus_sigma)
        ratio_min = numpy.array(cygA_mass_min_sigma/cygB_mass_min_sigma)

        # pyplot.figure(figsize=(12,9))
        # pyplot.loglog(r, cygA_gas_mass, label="A gas")
        # pyplot.loglog(r, cygA_dm_mass, label="A dm")
        # pyplot.loglog(r, cygA_mass, label="A total")

        # pyplot.loglog(r, cygA_mass, label="A val")
        # pyplot.loglog(r, cygA_mass_plus_sigma, label="A plus")
        # pyplot.loglog(r, cygA_mass_min_sigma, label="A min")
        # pyplot.loglog(r, cygB_mass, label="B val")
        # pyplot.loglog(r, cygB_mass_plus_sigma, label="B plus")
        # pyplot.loglog(r, cygB_mass_min_sigma, label="B min")

        # pyplot.plot(r, ratio, label="val")
        # pyplot.plot(r, ratio_plus, label="plus")
        # pyplot.plot(r, ratio_min, label="min")
        # pyplot.plot(loc=4)
        # pyplot.gca().set_xscale("log")

        # pyplot.show()

        pyplot.figure(figsize=(12, 9))
        pyplot.plot(r, ratio, c="b")
        #pyplot.plot(r, ratio_min, c="k")
        #pyplot.plot(r, ratio_plus, c="r")
        pyplot.fill_between(r, ratio_min, ratio_plus, facecolor="green", alpha=0.2)
        pyplot.gca().set_xscale("log")
        pyplot.axvline(cygA["r200"]*cm2kpc, c="r")
        pyplot.xlabel(r"$r$ [kpc]")
        pyplot.ylabel(r"Mass Ratio [Cyg$_{\rm A}$/Cyg$_{\rm B}$]")
        pyplot.xlim(1, 4000)
        pyplot.ylim(numpy.min(ratio_plus)-0.3, numpy.max(ratio_min)+0.3)
        pyplot.savefig("out/cygA_cygB_massRatio{0}{1}.png"\
            .format("_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"), dpi=300)

    if mode == "rhomassboth":
        print "Generating ugly two-panel plot. Density left, mass right."
        print "Both panels contain both gas and dark matter of both clusters"
        fig, (ax0, ax1) = pyplot.subplots(1, 2, sharex=True, figsize=(16, 8))

        pyplot.sca(ax0)
        pyplot.loglog(r, cygA_gas_rhom, label="cygA gas", c="r", ls="dotted")
        pyplot.loglog(r, cygA_dm_rhom, label="cygA dm", c="r")
        pyplot.loglog(r, cygB_gas_rhom, label="cygB gas", c="g", ls="dashed")
        pyplot.loglog(r, cygB_dm_rhom, label="cygB dm", c="g")

        pyplot.sca(ax1)
        pyplot.loglog(r, cygA_gas_mass, label="cygA gas", c="r", ls="dotted")
        pyplot.loglog(r, cygA_dm_mass, label="cygA dm", c="r")
        pyplot.loglog(r, cygB_gas_mass, label="cygB gas", c="g", ls="dashed")
        pyplot.loglog(r, cygB_dm_mass, label="cygB dm", c="g")

        for ax in [ax0, ax1]:
            pyplot.sca(ax)
            pyplot.xlabel(r"$r$ [kpc]")
        ax0.axvline(cygA["r200"]*cm2kpc, c="r")
        ax1.axvline(cygB["r200"]*cm2kpc, c="g")
        ax0.set_ylabel(r"$\rho$ [g/cm$^3$]")
        ax1.set_ylabel(r"$M(<r)$ [$M_{\odot}$]")

        ax0.legend(loc=3, fontsize=12)
        ax1.legend(loc=4, fontsize=12)
        pyplot.tight_layout()
        pyplot.savefig("out/cygA_cygB_massAndDensity{0}{1}.png"\
            .format("_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"), dpi=300)

    if mode == "massboth":
        print "Generating two-panel mass plot. CygA left, CygB right."
        fig, (ax0, ax1) = pyplot.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 8))

        pyplot.sca(ax0)
        pyplot.loglog(r, cygA_gas_mass+cygA_dm_mass, label="total", c="k")
        pyplot.loglog(r, cygA_dm_mass, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygA_gas_mass, label="gas", c="k", ls="dashed")
        ax0.axvline(cygA["r200"]*cm2kpc, c="k", ls="dotted")

        pyplot.sca(ax1)
        pyplot.loglog(r, cygB_gas_mass+cygB_dm_mass, label="total", c="k")
        pyplot.loglog(r, cygB_dm_mass, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygB_gas_mass, label="gas", c="k", ls="dashed")
        ax1.axvline(cygB["r200"]*cm2kpc, c="k", ls="dotted")

        for ax in [ax0, ax1]:
            pyplot.sca(ax)
            pyplot.xlabel(r"$r$ [kpc]")
            pyplot.legend(loc=4, fontsize=12)

        ax0.set_ylabel(r"$M(<r)$ [$M_{\odot}$]")
        pyplot.tight_layout()
        pyplot.savefig("out/cygA_cygB_mass{0}{1}.png"\
            .format("_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"), dpi=300)

    if mode == "masssameplot":
        print "Generating single plot with both CygA and CygB masses."
        fig = pyplot.figure(figsize=(12, 9))

        pyplot.loglog(r, cygA_gas_mass+cygA_dm_mass, label="CygA total", c="r")
        # pyplot.loglog(r, cygA_dm_mass, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygA_gas_mass, label="CygA gas", c="r", ls="dashed")
        pyplot.axvline(cygA["r200"]*cm2kpc, c="r", ls="dotted")

        pyplot.loglog(r, cygB_gas_mass+cygB_dm_mass, label="CygB total", c="k")
        # pyplot.loglog(r, cygB_dm_mass, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygB_gas_mass, label="CygB gas", c="k", ls="dashed")
        pyplot.axvline(cygB["r200"]*cm2kpc, c="k", ls="dashed")

        # Smith et al. quote a mass at 500 kpc, assuming H0 = 50.
        # So converting back to arcsec, and then assuming H0 = 70 to kpc
        pyplot.axvline(500/1.527*1.091, c="red", ls="solid")

        pyplot.xlim(5, 3000)
        pyplot.ylim(1e10, 5e14)
        pyplot.xlabel(r"$r$ [kpc]")
        pyplot.legend(loc=2, fontsize=12)

        pyplot.ylabel(r"$M(<r)$ [$M_{\odot}$]")
        pyplot.savefig("out/cygA_cygB_mass_sameplot{0}{1}.png"\
            .format("_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"), dpi=300)

    if mode == "rhoboth":
        print "Generating two-panel denisty plot. CygA left, CygB right."
        fig, (ax0, ax1) = pyplot.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 8))

        pyplot.sca(ax0)
        pyplot.loglog(r, cygA_dm_rhom, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygA_gas_rhom, label="gas", c="k", ls="dashed")
        ax0.axvline(cygA["r200"]*cm2kpc, c="k", ls="dotted")

        pyplot.sca(ax1)
        pyplot.loglog(r, cygB_dm_rhom, label="dm", c="k", ls="dotted")
        pyplot.loglog(r, cygB_gas_rhom, label="gas", c="k", ls="dashed")
        ax1.axvline(cygB["r200"]*cm2kpc, c="k", ls="dotted")

        for ax in [ax0, ax1]:
            pyplot.sca(ax)
            pyplot.xlabel(r"$r$ [kpc]")
            pyplot.legend(loc=3, fontsize=12)

        ax0.set_ylabel(r"$\rho$ [g/cm$^3$]")
        pyplot.tight_layout()
        pyplot.savefig("out/cygA_cygB_density{0}{1}.png"\
            .format("_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"), dpi=300)

    if "single" in mode:
        print "Plotting single cluster!"
        # Find which cluster we are plotting
        if cygA and not cygB:
            print "    sub cluster cygA selected"
            gas_rhom = cygA_gas_rhom
            gas_mass = cygA_gas_mass
            dm_rhom  = cygA_dm_rhom
            dm_mass  = cygA_dm_mass
            parms = cygA
            observed = cygA_observed
        elif cygB and not cygA:
            print "    sub cluster cygB selected"
            gas_rhom = cygB_gas_rhom
            gas_mass = cygB_gas_mass
            dm_rhom  = cygB_dm_rhom
            dm_mass  = cygB_dm_mass
            parms = cygB
            observed = cygB_observed
        else:
            print "Incorrect usage"
            return

    if mode == "rhosingle":
        print "Generating single cluster density plot"
        gas_rhom_discrete = gas_density_beta(observed.radius, parms["rho0"],
                parms["rc"]*cm2kpc, parms["beta"])
        fig = plot_observed_cluster(observed, gas_rhom_discrete)

        # Plot the analytical models: "continuous" radii and function values
        pyplot.figure(fig.number)
        ax, ax_r = fig.axes
        pyplot.sca(ax)

        pyplot.plot(r, dm_rhom, c="k", ls="solid")
        pyplot.axhline(200*rho_crit(), c="k")
        rho_avg_200 = parms["M200"] / (4./3 * numpy.pi * p3(parms["r200"]))
        pyplot.axhline(rho_avg_200, c="r")

        pyplot.axvline(x=parms["r200"]*cm2kpc, c="r", ls="solid")


        # Set axis limits
        if observed.name == "cygA":
            ax_r.set_xlim(1, 3000)
            ax.set_xlim(1, 3000)
            ax.set_ylim(1e-29, 1e-22)
        else:
            ax_r.set_xlim(40, 5000)
            ax.set_xlim(40, 5000)
            ax.set_ylim(1e-29, 3e-25)

        pyplot.savefig("out/density_profile_with_dm_{0}{1}{2}.png"\
            .format(observed.name, "_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"))

    if mode == "bfsingle":
        print "Generating plot of the baryon fraction as a function of radius"
        # Here we find r500
        lower = 10 * kpc2cm
        upper = 1000 * kpc2cm
        epsilon = 0.0001
        while upper/lower > 1+epsilon:
            r500 = (lower+upper)/2
            rho500_over_rhocrit = (parms["M200"] / (4./3 * numpy.pi * p3(r500))) / rho_crit()
            if rho500_over_rhocrit < 500:
                upper = r500
            if rho500_over_rhocrit > 500:
                lower = r500

        bf = gas_mass/(dm_mass+gas_mass)

        pyplot.figure(figsize=(12,9))

        pyplot.plot(r, bf, c="k")
        pyplot.gca().set_xscale("log")

        pyplot.axvline(parms["r200"]*cm2kpc, c="r", ls="dotted")
        pyplot.axvline(r500*cm2kpc, c="r", ls="dashed")
        pyplot.axhline(0.17, c="k", ls="dotted", label=r"$\overline{\rho(r_{200})}=0.17$")
        pyplot.axhline(0.14, c="k", ls="dashed", label=r"$\overline{\rho(r_{500})}=0.14$")
        pyplot.xlabel(r"$r$ [kpc]")
        pyplot.ylabel(r"$b_f$")

        pyplot.legend(loc=2, fontsize=12)
        pyplot.savefig("out/baryon_fraction_{0}{1}{2}.png"\
            .format(observed.name, "_freebeta" if free_beta else "",
                    "_800ksec" if oldICs else "_900ksec"))

    print


def is_solution_unique(rc, rho0, beta, observed):
    gas_rhom = gas_density_beta(observed.radius, rho0, rc*cm2kpc, beta)
    fig = plot_observed_cluster(observed, gas_rhom)
    ax, ax_r = fig.axes
    pyplot.sca(ax)
    # pyplot.ion()

    # Set axis limits
    pyplot.axhline(200*rho_crit(), c="k")
    if observed.name == "cygA":
        ax_r.set_xlim(1, 3000)
        ax.set_xlim(1, 3000)
        ax.set_ylim(1e-29, 1e-22)
    else:
        ax_r.set_xlim(40, 5000)
        ax.set_xlim(40, 5000)
        ax.set_ylim(1e-29, 3e-25)

    searchrange = numpy.arange(400, 1000, 1)
    N = len(searchrange)
    for n, r200 in enumerate(searchrange*kpc2cm):
        # bisection
        bf = 0.17  # Mission critical assumption that bf = 0.17 at r200! Planelles+ 2013

        Mgas200 = M_gas_below_r(r200, rho0, rc, beta)
        Mdm200 = Mgas200 * (1/bf - 1)

        M200 = Mgas200 + Mdm200

        cNFW = concentration_parameter(M200)
        rs = r200 / cNFW
        a = hernquist_a(rs, cNFW)

        Mdm = Mdm200 / M_dm_below_r(r200, 1, a)

        """ Now rho_average(r200)/rhocrit should equal 200.
                If not? Try different r200"""
        rho200_over_rhocrit = ( M200 / (4./3 * numpy.pi * p3(r200))) / rho_crit()

        pyplot.figure(fig.number)

        dm_rhom = dm_density_hernquist(observed.radius*kpc2cm, Mdm, a)

        dmrho_line = pyplot.plot(observed.radius, dm_rhom, c="k", ls="solid")
        rho_avg_200 = M200 / (4./3 * numpy.pi * p3(r200))
        avg_line = pyplot.axhline(y=rho_avg_200, c="r", ls="solid")

        # Indicate bisection bounds
        r200_line = pyplot.axvline(x=r200*cm2kpc, c="r", ls="solid")

        # Plot textbox with bisection info
        bisection_info = "r200 : {0:.1f}\n".format(r200*cm2kpc) \
            + r"$\frac{{\rho(r_{{200}})}}{{\rho_{{\rm crit}}}}$: {0:.1f}".format(rho200_over_rhocrit)

        if observed.name == "cygA":
            textX = 3
            textY = 2e-29
        else:
            textX = 50
            textY = 2e-29

        info = pyplot.text(textX, textY, bisection_info, size=18,
                           ha="left", va="bottom",
                           bbox=dict(boxstyle="round",
                                     ec=(1., 0.5, 0.5),
                                     fc=(1., 0.8, 0.8),
                                     )
                          )
        if 199.9 < rho200_over_rhocrit < 200.1:
            print "Solution found:", r200
        # pyplot.draw()
        # pyplot.pause(0.001)
    # convert -delay 100 -loop 0 out/uniqueness-check_cygA*.png out/uniqueness_cygA.gif
        pyplot.savefig("out/uniqueness-check_{1}_{0:03d}.png".format(n, observed.name))
        r200_line.remove()
        avg_line.remove()
        dmrho_line.pop().remove()
        info.remove()

        print_progressbar(n, N)


if __name__ == "__main__":
    # If visualise is True we create plots of the bisection method
    visualise=True
    oldICs=False
    discard_firstbins = True
    discard_lastbins = False
    free_beta = False
    poster_style = True

    print "Reading Chandra observed density profiles..."
    print 80*"-"
    if oldICs:
        print "\nWARNING: Using old density profiles (799.5 ksec)."
        raw_input("Are you sure you want to continue?\n")
    cygA_observed = ObservedCluster("cygA", oldICs=oldICs)
    cygB_observed = ObservedCluster("cygB", oldICs=oldICs)
    print ".... done reading Chandra observed density profiles."
    print 80*"-"

    print "Obtaining central density and core radius."
    print 80*"-"
    # Due to pile up in inner region (?). Also inside kernel: cannot model
    # in a stable way
    if discard_firstbins:
        print "WARNING: Discarding first three CygA bins."
        cygA_observed.radius = cygA_observed.radius[3:]
        cygA_observed.binsize = cygA_observed.binsize[3:]
        cygA_observed.density= cygA_observed.density[3:]
        cygA_observed.density_std = cygA_observed.density_std[3:]
        cygA_observed.number_density = cygA_observed.number_density[3:]
        cygA_observed.number_density_std = cygA_observed.number_density_std[3:]
    if discard_lastbins:
        print "WARNING: Discarding last CygA bins."
        cygA_observed.radius = cygA_observed.radius[:-40]
        cygA_observed.binsize = cygA_observed.binsize[:-40]
        cygA_observed.density= cygA_observed.density[:-40]
        cygA_observed.density_std = cygA_observed.density_std[:-40]
        cygA_observed.number_density = cygA_observed.number_density[:-40]
        cygA_observed.number_density_std = cygA_observed.number_density_std[:-40]

    if free_beta:
        print "INFO: Using free beta model.\n"
        cygA_fit, cygA_ml_vals, cygA_ml_covar = \
            fit_betamodel_to_chandra(cygA_observed,
                parm=[0.1, 10, 0.67], free_beta=free_beta)
        cygB_fit, cygB_ml_vals, cygB_ml_covar = \
            fit_betamodel_to_chandra(cygB_observed,
                parm=[0.001, 100, 0.79], free_beta=free_beta)
    else:
        print "INFO: Using beta=2/3 model.\n"
        cygA_fit, cygA_ml_vals, cygA_ml_covar = \
            fit_betamodel_to_chandra(cygA_observed, parm=[0.1, 10])
        cygB_fit, cygB_ml_vals, cygB_ml_covar = \
            fit_betamodel_to_chandra(cygB_observed, parm=[0.001, 100])
    print 80*"-"

    print "Obtaining dark matter density/mass profile, r200, cNFW, etc."
    print 80*"-"

    cygA_ne0 = cygA_ml_vals[0]
    cygA_rho0 = ne_to_rho(cygA_ne0)  # g/cm**3
    cygA_rc = cygA_ml_vals[1] * kpc2cm
    cygA_beta = None if not free_beta else cygA_ml_vals[2]

# convert -delay 100 -loop 0 out/findmass_cygA_*.png out/bisection_cygA.gif
    # is_solution_unique(cygA_rc, cygA_rho0, cygA_observed)

    cygA = obtain_M200_bisection(cygA_rc, cygA_rho0, cygA_beta, verbose=False,
                                 visualise=visualise, observed=cygA_observed)

    cygA_sigma = numpy.sqrt(numpy.diag(cygA_ml_covar))
    cygA_ne0_sigma = cygA_sigma[0]
    cygA_rho0_sigma = ne_to_rho(cygA_ne0_sigma)
    cygA_rc_sigma = cygA_sigma[1] * kpc2cm
    cygA_beta_sigma = None if not free_beta else cygA_sigma[2]

    cygA["rho0_sigma"] = cygA_rho0_sigma
    cygA["rc_sigma"] = cygA_rc_sigma
    cygA["beta_sigma"] = cygA_beta_sigma

    print "CygA"
    print_inferred_values(cygA)

    print "CygB"

    cygB_ne0 = cygB_ml_vals[0]
    cygB_rho0 = ne_to_rho(cygB_ne0)  # g/cm**3
    cygB_rc = cygB_ml_vals[1] * kpc2cm
    cygB_beta = None if not free_beta else cygB_ml_vals[2]

    #is_solution_unique(cygB_rc, cygB_rho0, cygB_observed)
    #import sys; sys.exit(0)

    cygB = obtain_M200_bisection(cygB_rc, cygB_rho0, cygB_beta, verbose=False,
                                 visualise=visualise, observed=cygB_observed)

    cygB_sigma = numpy.sqrt(numpy.diag(cygB_ml_covar))
    cygB_ne0_sigma = cygB_sigma[0]
    cygB_rho0_sigma = ne_to_rho(cygB_ne0_sigma)
    cygB_rc_sigma =  cygB_sigma[1] * kpc2cm
    cygB_beta_sigma = None if not free_beta else cygB_sigma[2]

    cygB["rho0_sigma"] = cygB_rho0_sigma
    cygB["rc_sigma"] = cygB_rc_sigma
    cygB["beta_sigma"] = cygB_beta_sigma

    print_inferred_values(cygB)

    print "cygA_M200            = {0:1.4e} MSun".format(cygA["M200"] * g2msun)
    print "cygB_M200            = {0:1.4e} MSun".format(cygB["M200"] * g2msun)
    print
    print "M_CygA/M_CygB        = {0:1.4f}".format(cygA["M200"]/cygB["M200"])
    print "M_CygB/M_CygA        = {0:1.4f}".format(cygB["M200"]/cygA["M200"])
    print "Mtotal               = {0:1.4e} MSun".format((cygA["M200"] + cygB["M200"]) * g2msun)
    print

    print 80*"-"
    print "Plotting the results..."
    print 80*"-"

    make_plot(cygA, None, cygA_observed, None, mode="nfwsingle")
    pyplot.show()
    raw_input("Press enter to continue...\n")

    # Plot density+mass profiles (gas + dm in same plot); density left, mass right
    make_plot(cygA, cygB, mode="massboth")
    make_plot(cygA, cygB, mode="rhoboth")
    make_plot(cygA, cygB, mode="masssameplot")

    make_plot(cygA, None, cygA_observed, None, mode="rhosingle")
    make_plot(None, cygB, None, cygB_observed, mode="rhosingle")

    # Plot mass ratio
    make_plot(cygA, cygB, mode="ratio")

    # Plot baryon fraction
    make_plot(cygA, None, cygA_observed, None, mode="bfsingle")
    make_plot(None, cygB, None, cygB_observed, mode="bfsingle")

    #pyplot.show()

    print 80*"-"
    print "End of pipeline."
    print 80*"-"
