"""
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Mon May 16, 2016 05:23 pm

Script to obtain M200 for CygA, CygB, thus xM (mass ratio),
and cNFW for cygA and cygB using reverse of Toycluster setup.c.
Code based on Julius Donnert's 20160617 cyg.pro script.

Script is designed to stand alone from the rest of my code base.
Only depends on Numpy, but could easily be rewritten to use plain Python

"""

import numpy
from matplotlib import pyplot
import matplotlib
pyplot.rcParams.update({'font.size': 22})

from cluster import ObservedCluster
from cluster import AnalyticalCluster


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


def gas_density_beta(r, rho0, rc):
    """ Beta-model for gas density profile (Cavaliere & Fusco-Femiano 1978)
    with a fixed value beta = 2/3 (Mastropietro & Burkert 2008)

    Donnert (2014): "Should any other value of beta be used, this step
                     would involve the confluent hypergeometric function,
                     which usually makes the analytical solution of the
                     hydrostatic equation impossible. Another elegant
                     solution would be beta=1, which however is unfavoured
                     from observations"
    """

    beta = 2./3
    rho_gas_beta = rho0 * (1 + p2(r)/p2(rc))**(-3*beta/2.)
    return rho_gas_beta


def M_gas_below_r(r, rho0, rc):
    """ M(<r) for Hydrostatic Equilibrium, Spherical Symmetry
    Beta-Model where beta is fixed to 2/3 """

    M_gas_below_r = 4*numpy.pi*p3(rc)*rho0 * (r/rc - numpy.arctan(r/rc))
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


def plot_observed_cluster(observed, analytical_density):
    fig, (ax, ax_r) = pyplot.subplots(2, 2, sharex=True, figsize=(16, 12))
    gs1 = matplotlib.gridspec.GridSpec(3, 3)
    gs1.update(hspace=0)
    ax = pyplot.subplot(gs1[:-1,:])
    ax_r = pyplot.subplot(gs1[-1,:])  # residuals

    # Plot data
    pyplot.sca(ax)
    pyplot.errorbar(observed.radius+observed.binsize/2,
                    observed.density, xerr=observed.binsize/2,
                    yerr=observed.density_std, marker="o", ms=6, ls="", c="g")
                    #label="800 ks Chandra\n(Wise+ 2016, in prep)")

    # Plot analytical gas profile
    pyplot.plot(observed.radius, analytical_density,
                c="k", ls="dashed")#, label="gas")

    # Plot Residuals
    pyplot.sca(ax_r)
    residual_density = (observed.density - analytical_density)/observed.density
    pyplot.errorbar(observed.radius+observed.binsize/2, residual_density,
            yerr=observed.density_std/observed.density, c="k", drawstyle="steps-mid")
    # Show dashed residuals zero line
    ax_r.axhline(y=0, lw=2, ls="dashed", c="k")

    # Set axis labels
    ax.set_xlabel(r"$\rho(r)$ [g/cm$^{-3}$]")
    ax_r.set_xlabel(r"$r$ [kpc]")
    ax_r.set_ylabel("Residuals")

    # Set logscale, but the residual y-axis is not logarithmic!
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax_r.set_xscale("log")

    # Set residual y-limits: show 50% deviations from the data
    ax_r.set_ylim(-0.5, 0.5)

    # Fix for overlapping y-axis markers
    from matplotlib.ticker import MaxNLocator
    ax.tick_params(labelbottom='off')
    nbins = len(ax_r.get_yticklabels())
    ax_r.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

    return fig


def obtain_M200_bisection(rc, rho0, verbose=False,
                          visualise=False, observed=None):
    """ We follow the Toycluster (Donnert 2014) setup.c method in reverse.
    If we assume a value for r200 and we assume the baryon fraction at r200
    is equal to 0.17 we are able to obtain a total mass M200 for which
    rho_average(r200) == 200 rho_crit."""

    if visualise:
        gas_rhom = gas_density_beta(observed.radius, rho0, rc*cm2kpc)
        n=0

    # Find r200 such that rho200 / rho_crit == 200
    lower = 10 * kpc2cm
    upper = 2000 * kpc2cm

    if visualise and observed.name == "cygB":
        lower = 100 * kpc2cm

    # bisection method
    epsilon = 0.0001
    while upper/lower > 1+epsilon:
        # bisection
        r200 = (lower+upper)/2.

        bf = 0.17  # Mission critical assumption that bf = 0.17 at r200! Planelles+ 2013

        Mgas200 = M_gas_below_r(r200, rho0, rc)
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

            pyplot.plot(observed.radius, dm_rhom, c="k", ls="solid")
            pyplot.axhline(200*rho_crit(), c="k")
            rho_avg_200 = M200 / (4./3 * numpy.pi * p3(r200))
            pyplot.axhline(rho_avg_200, c="r")

            # Indicate bisection bounds
            pyplot.axvline(x=lower*cm2kpc, c="k", ls="dotted")
            pyplot.axvline(x=upper*cm2kpc, c="k", ls="dotted")
            pyplot.axvline(x=r200*cm2kpc, c="r", ls="solid")

            # Plot textbox with bisection info
            bisection_info = "\nlower: {0:.1f}\n".format(lower*cm2kpc) \
                + "r200 : {0:.1f}\n".format(r200*cm2kpc) \
                + "upper: {0:.1f}\n".format(upper*cm2kpc) \
                + r"$\frac{{\rho(r_{{200}})}}{{\rho_{{\rm crit}}}}$: {0:.1f}".format(rho200_over_rhocrit)

            if observed.name == "cygA":
                textX = 3
                textY = 2e-29
            else:
                textX = 50
                textY = 2e-29

            pyplot.text(textX, textY, bisection_info, size=18,
                        ha="left", va="bottom",
                        bbox=dict(boxstyle="round",
                                  ec=(1., 0.5, 0.5),
                                  fc=(1., 0.8, 0.8),
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


            # pyplot.show()
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
    halo['r200'] = r200
    halo['rho200_over_rhocrit'] = rho200_over_rhocrit
    halo['rho0'] = rho0
    halo['rc'] = rc
    halo['Mgas200'] = Mgas200
    halo['Mdm200'] = Mdm200
    halo['M200'] = M200
    halo['cNFW'] = cNFW
    halo['rs'] = rs
    halo['a'] = a
    halo['Mdm'] = Mdm

    return halo


def print_inferred_values(halo):
    bf_200 = halo["Mgas200"]/(halo["Mdm200"]+halo["Mgas200"])
    print "r200                   = {0:3.1f}".format(halo["r200"] * cm2kpc)
    print "rho_avg(r200)/rho_crit = {0:.1f}".format(halo["rho200_over_rhocrit"])
    print "bf200                  = {0:1.4f}".format(bf_200)
    print "rho0                   = {0:1.4e}".format(halo["rho0"])
    print "rc                     = {0:.3f}".format(halo["rc"] * cm2kpc)
    print "Mgas200                = {0:1.4e}".format(halo["Mgas200"] * g2msun)
    print "Mdm200                 = {0:1.4e}".format(halo["Mdm200"] * g2msun)
    print "M200                   = {0:1.4e}".format(halo["M200"] * g2msun)
    print "cNFW                   = {0:1.4f}".format(halo["cNFW"])
    print "rs                     = {0:3.1f}".format(halo["rs"] * cm2kpc)
    print "a                      = {0:3.1f}".format(halo["a"] * cm2kpc)
    print "Mdm                    = {0:1.4e}".format(halo["Mdm"] * g2msun)
    print


def plot_profiles(cygA, cygB, ratio=False):
    r = numpy.arange(1, 1e4, 0.1)  # kpc!

    cygA_gas_rhom = gas_density_beta(r, cygA['rho0'], cygA['rc']*cm2kpc)
    cygA_gas_mass = M_gas_below_r(r, cygA['rho0']*g2msun/p3(cm2kpc), cygA['rc']*cm2kpc)
    cygB_gas_rhom = gas_density_beta(r, cygB['rho0'], cygB['rc']*cm2kpc)
    cygB_gas_mass = M_gas_below_r(r, cygB['rho0']*g2msun/p3(cm2kpc), cygB['rc']*cm2kpc)

    cygA_dm_rhom = dm_density_hernquist(r*kpc2cm, cygA['Mdm'], cygA['a'])
    cygA_dm_mass = M_dm_below_r(r, cygA['Mdm']*g2msun, cygA['a']*cm2kpc)
    cygB_dm_rhom = dm_density_hernquist(r*kpc2cm, cygB['Mdm'], cygB['a'])
    cygB_dm_mass = M_dm_below_r(r, cygB['Mdm']*g2msun, cygB['a']*cm2kpc)

    if ratio:
        fig = pyplot.figure(figsize=(12, 9))
        cygA_mass = cygA_gas_mass + cygA_dm_rhom
        cygB_mass = cygB_gas_mass + cygB_dm_rhom

        ratio = numpy.array(cygA_mass/cygB_mass)
        pyplot.plot(r, ratio, c="b")
        pyplot.fill_between(r, 0.9*ratio, 1.1*ratio, facecolor='green', alpha=0.2)

        pyplot.gca().set_xscale("log")
        pyplot.axvline(500, c="k", ls="dotted")
        pyplot.axvline(cygA['r200']*cm2kpc, c="r")
        pyplot.xlabel(r"$r$ [kpc]")
        pyplot.ylabel(r"Mass Ratio [Cyg$_{\rm A}$/Cyg$_{\rm NW}$]")
        pyplot.xlim(100, 1000)
        pyplot.ylim(0, 10)
        pyplot.savefig("out/cygA_cygNW_massRatio.png", dpi=300)
        #pyplot.show()
    else:
        fig, (ax0, ax1) = pyplot.subplots(1, 2, sharex=True, figsize=(16, 8))

        pyplot.sca(ax0)
        pyplot.loglog(r, cygA_gas_rhom, label="cygA gas", c="r", ls="dotted")
        pyplot.loglog(r, cygA_dm_rhom, label="cygA dm", c="r")
        pyplot.loglog(r, cygB_gas_rhom, label="cygNW gas", c="g", ls="dashed")
        pyplot.loglog(r, cygB_dm_rhom, label="cygNW dm", c="g")

        pyplot.sca(ax1)
        pyplot.loglog(r, cygA_gas_mass, label="cygA gas", c="r", ls="dotted")
        pyplot.loglog(r, cygA_dm_mass, label="cygA dm", c="r")
        pyplot.loglog(r, cygB_gas_mass, label="cygNW gas", c="g", ls="dashed")
        pyplot.loglog(r, cygB_dm_mass, label="cygNW dm", c="g")

        for ax in [ax0, ax1]:
            pyplot.sca(ax)
            pyplot.xlabel(r"$r$ [kpc]")
        ax0.axvline(cygA['r200']*cm2kpc, c="r")
        ax1.axvline(cygB['r200']*cm2kpc, c="g")
        ax0.set_ylabel(r"$\rho$ [g/cm$^3$]")
        ax1.set_ylabel(r"$M(<r)$ [$M_{\odot}$]")

        ax0.legend(loc=3)
        ax1.legend(loc=4)
        pyplot.tight_layout()
        pyplot.savefig("out/cygA_cygNW_massAndDensity.png", dpi=300)


if __name__ == "__main__":
    cygA_observed = ObservedCluster("cygA")
    cygB_observed = ObservedCluster("cygB")

    print "CygnusA"
    # obtained from beta (2/3) model to Chandra data
    cygA_rc = 27.645 * kpc2cm
    cygA_ne0 = 0.136  # 1/cm**3
    cygA_rho0 = ne_to_rho(cygA_ne0)

    # convert -delay 100 -loop 0 out/findmass_cygA_*.png out/bisection_cygA.gif
    cygA = obtain_M200_bisection(cygA_rc, cygA_rho0, verbose=True,
                                 visualise=True, observed=cygA_observed)
    print_inferred_values(cygA)
    print

    print "CygnusB"
    # obtained from beta (2/3) model to Chandra data
    cygB_rc = 290.909 * kpc2cm
    cygB_ne0 = 1.9397e-03  # 1/cm**3
    cygB_rho0 = ne_to_rho(cygB_ne0)

    cygB = obtain_M200_bisection(cygB_rc, cygB_rho0, verbose=True,
                                 visualise=True, observed=cygB_observed)
    print_inferred_values(cygB)
    print

    print "cygA M200              = {0:1.4e} MSun".format(cygA['M200'] * g2msun)
    print "cygB_M200              = {0:1.4e} MSun".format(cygB['M200'] * g2msun)
    print
    print "M_CygA/M_CygB          = {0:1.4f}".format(cygA['M200']/cygB['M200'])
    print "M_CygB/M_CygA          = {0:1.4f}".format(cygB['M200']/cygA['M200'])
    print "Mtotal                 = {0:1.4e} MSun".format((cygA['M200'] + cygB['M200']) * g2msun)

    plot_profiles(cygA, cygB)
    plot_profiles(cygA, cygB, ratio=True)
