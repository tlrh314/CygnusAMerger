import numpy
import pandas
import csv
import scipy
from scipy import stats

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
import globals
from cosmology import CosmologyCalculator
from ioparser import parse_toycluster_parms


class ObservedCluster(object):
    """ Observed Cygnus A-Cygnus B situation """
    def __init__(self, name):
        """ Generate cluster instance with volume, pressure, density,
            Compton Y for bins with given inner and outer radii.
            Data obtained by M. N. de Vries from 800 ksec Chandra data """

        self.name = name
        density_file = "data/{0}_1T_fixnH_pressureprofile.dat".format(name)
        radius_file = "data/{0}_sn100_sbprofile.dat".format(name)

        self.parse_data(density_file, radius_file)

        self.cc = CosmologyCalculator()
        arcsec2kpc = self.cc.kpc_DA # | units.kpc

        self.radius = (self.inner_radius+self.outer_radius)/2 * arcsec2kpc
        self.binsize = (self.outer_radius - self.inner_radius) * arcsec2kpc

    def __str__(self):
        tmp = ""

        tmp += "\nbin_number\n" + str(self.bin_number)
        tmp += "\n\nbin_volume\n" + str(self.bin_volume)
        tmp += "\n\ndensity\n" + str(self.density)
        tmp += "\n\ndensity_std\n" + str(self.density_std)
        tmp += "\n\npressure\n" + str(self.pressure)
        tmp += "\n\npressure_std\n" + str(self.pressure_std)
        tmp += "\n\ncompton_y\n" + str(self.compton_y)

        tmp += "\n\nbin_number_sn100\n" + str(self.bin_number_sn100)
        tmp += "\n\ninner_radius\n" + str(self.inner_radius)
        tmp += "\n\nouter_radius\n" + str(self.outer_radius)
        tmp += "\n\nsource_sb\n" + str(self.source_sb)
        tmp += "\n\nsource_sb_std\n" + str(self.source_sb_std)
        tmp += "\n\nbackground_sb\n" + str(self.background_sb)
        tmp += "\n\nbackground_sb_std\n" + str(self.background_sb_std)

        return tmp

    def parse_data(self, density_file, radius_file):
        # Read Martijn Data
        raw = pandas.read_csv(density_file, delimiter="|")
        # print raw.keys()

        self.bin_number = raw[" Bin number "].as_matrix()
        self.bin_volume = raw[" Bin volume (cm^3) "].as_matrix()
        self.density = raw["   Density (cm^-3) "].as_matrix()
        self.density_std = raw["     Sigma density "].as_matrix()
        self.pressure = raw[" Pressure (erg cm^-3) "].as_matrix()
        self.pressure_std = raw["    Sigma Pressure "].as_matrix()
        self.compton_y = raw[" Compton Y parameter "].as_matrix()

        raw_sn100 = pandas.read_csv(radius_file, delimiter="|")
        # print raw_sn100.keys()

        self.inner_radius = raw_sn100[" Inner radius (arcsec) "].as_matrix()
        self.outer_radius = raw_sn100[" Outer radius (arcsec) "].as_matrix()

        if " Source SB (counts/cm^2/arcsec^2/s) " in raw_sn100.keys():  # CygA
            self.source_sb = raw_sn100[" Source SB (counts/cm^2/arcsec^2/s) "].as_matrix()
        elif " SB (cnts/cm^2/arcsec^2/s) " in raw_sn100.keys():  # CygB
            self.source_sb = raw_sn100[" SB (cnts/cm^2/arcsec^2/s) "].as_matrix()

        self.source_sb_std = raw_sn100["   Sigma source SB "].as_matrix()

        if " Background SB(counts/cm^2/arcsec^2/s) " in raw_sn100.keys():  # CygA
            self.background_sb = raw_sn100[" Background SB(counts/cm^2/arcsec^2/s) "].as_matrix()
        elif " Bkg SB(cnts/cm^2/arcsec^2/s) " in raw_sn100.keys():  # CygB
            self.background_sb = raw_sn100[" Bkg SB(cnts/cm^2/arcsec^2/s) "].as_matrix()

        if " Sigma Background SB " in raw_sn100.keys():
            self.background_sb_std = raw_sn100[" Sigma Background SB "].as_matrix()
        elif "      Sigma Bkg SB " in raw_sn100.keys():
            self.background_sb_std = raw_sn100["      Sigma Bkg SB "].as_matrix()

        self.bin_number_sn100 = raw_sn100[" Bin number "].as_matrix()

        if not numpy.array_equal(self.bin_number, self.bin_number_sn100):
# Yeah, I know way too generic to raise, but w/e dont wanna define ReadingTwoDataFilesIsProneToMakeMistakesExeption *_*
            raise Exception("Mismatch in the bin numbers of the two data input files")


def simple_plot(cluster, yparm="density"):
    """ Possible properties to plot on y-axis are:
            bin_volume, density*, pressure*, compton_y, source_sb*, background_sb*
        * means std is available --> errorbar plot
        x-axis is computer from inner_radius and/or outer_radius
    """

    density_as_mass_density = False
    if density_as_mass_density and yparm == "density":
        density = globals.m_p * globals.mu * cluster.density
    else:
        density = cluster.density

    pyplot.figure(figsize=(12, 9))
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    if yparm not in ["bin_volume", "compton_y"]:
        pyplot.errorbar(cluster.radius, getattr(cluster, yparm),
            yerr=getattr(cluster, yparm+"_std"))
    else:
        pyplot.plot(cluster.radius, getattr(cluster, yparm))

    pyplot.xlabel("Radius")
    pyplot.ylabel(yparm)
    # pyplot.legend()


def Gas_density_profile(parm, r):
    """ Beta model if len(parm)==2 (Donnert 2014);
        Cut-off beta profile at rc and rcut if len(parm)==3 (Donnert et al. 2016, in prep)
        Double cut-off beta profile if len(parm)==5. See Toycluster/setup.c """
    DOUBLE_BETA_COOL_CORES = False
    rho0 = parm[0]
    rc = parm[1]

    rho = rho0 / (1 + p2(r/rc))   # not cut-off (Donnert 2014)
    if len(parm) == 3:
        rcut = parm[2]
        rho /= (1 + p3(r/rcut) * (r/rcut))

    if len(parm) == 5:  # Double Beta model (for Coolcore clusters)
        Rho0_Fac = parm[3]
        Rc_Fac = parm[4]

        rho0_cc = rho0 * Rho0_Fac
        rc_cc = rc / Rc_Fac

        rho += rho0_cc / (1 + p2(r/rc_cc)) / (1 + p3(r/rcut) * (r/rcut));

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


def fit_betamodel_to_chandra(cluster, parm=[1., 1., 1.], double=False):
    """ Fit betamodel to cluster observation, or
        double betamodel when len(parm)==5. """

    # Fit to data
    if cluster.name == "cygA":
        if len(parm) == 3:
            bounds = [(None, None), (None, None), (800, 1400)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1400), (None, None), (None, None)]
        else:
            bounds = [(None, None), (100, 400)]
    elif cluster.name == "cygB":
        if len(parm) == 3:
            bounds = [(None, None), (None, None), (600, 800)]
        elif len(parm) == 5:
            bounds = [(None, None), (None, None), (800, 1400), (None, None), (None, None)]
        else:
            bounds = [(None, None), (None, None)]
    result = scipy.optimize.minimize(stat, parm,
            args=(cluster.radius, cluster.density, cluster.density_std),
            method='L-BFGS-B', bounds=bounds)

    # Obtain and print MLEs
    obtain_mles(cluster, result)

    return result


def obtain_mles(cluster, result):
    ml_vals = result["x"]
    ml_func = result["fun"]

    moddof = len(ml_vals)  # Model degrees of freedom; nr of fit parameters
    dof = len(cluster.density) - moddof  # degrees of freedom

    ch = scipy.stats.chi2(dof)
    pval = 1.0 - ch.cdf(ml_func)

    print cluster.name
    print "Results for the '{0}' model:".format("Cut-off beta")
    print "  Using scipy.optimize.minimze to minimize chi^2 yields:"
    print "    n_e,0       = {0:.3f}".format(ml_vals[0])
    print "    r_c         = {0:.3f}".format(ml_vals[1])
    if len(ml_vals) == 3:
        print "    r_cut       = {0:.3f}".format(ml_vals[2])
    print "    chisq/dof   = {0:.3f}".format(ml_func/dof)
    print "    p-value     = {0:.5f}".format(pval)

    ml_vals, ml_covar = scipy.optimize.curve_fit(
            Gas_density_profile_wrapper,
            cluster.radius, cluster.density, p0=ml_vals,
            sigma=cluster.density_std,)
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


def plot_fit_results(cluster, result, save=False):
    poster_style = True
    if poster_style:
        pyplot.style.use(["dark_background"])
        data_colour = (255./255, 64./255, 255./255)
        fit_colour = "white"
    else:
        data_colour = "r"
        fit_colour = "k"

    fig, (ax, ax_r) = pyplot.subplots(2, 2, sharex=True, figsize=(16, 12))
    gs1 = matplotlib.gridspec.GridSpec(3, 3)
    gs1.update(hspace=0)
    ax = pyplot.subplot(gs1[:-1,:])
    ax_r = pyplot.subplot(gs1[-1,:])  # residuals

    # Plot data
    pyplot.sca(ax)
    pyplot.title(cluster.name)
    pyplot.errorbar(cluster.radius+cluster.binsize/2, cluster.density, xerr=cluster.binsize/2,
            yerr=cluster.density_std, marker='o', ms=6, ls='', c=data_colour,
            label="800 ks Chandra\n(De Vries, 2016)",)

    # Plot fit results
    rho0_fit = result["x"][0]  # | units.g/units.cm**3
    rc_fit = result["x"][1]  # | units.kpc
    if len(result["x"]) == 3:
        rcut_fit = result["x"][2]  # | units.kpc
        analytical_density = Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
        label = r"cut-off $\beta$-model:"+"\n\t"+r"$n_{{e,0}} \,$ = {0:.2e}".format(rho0_fit)\
            +"\n\t"+r"$r_c \,$ = {0:.2f}".format(rc_fit)+"\n\t"+\
            r"$r_{{\rm cut}}$ = {0:.2f}".format(rcut_fit)
        residual_density = cluster.density - \
            Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
    elif len(result["x"]) == 5:
        Rho0_Fac_fit = result["x"][3] | units.kpc
        Rc_Fac_fit = result["x"][4] | units.kpc
        # TODO
        analytical_density = Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
        label = r"double $\beta$-model:"+"\n\t"+r"$n_{{e,0}} \,$ = {0:.2e}".format(rho0_fit)\
            +"\n\t"+r"$r_c \,$ = {0:.2f}".format(rc_fit)+"\n\t"+\
            r"$r_{{\rm cut}}$ = {0:.2f}".format(rcut_fit)
        residual_density = cluster.density - \
            Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
    else:
        analytical_density = Gas_density_profile((rho0_fit, rc_fit), cluster.radius)
        label = r"$\beta$-model:"+"\n\t"+r"$n_{{e,0}} \,$ = {0:.2e}".format(rho0_fit)\
            +"\n\t"+r"$r_c \,$ = {0:.2f}".format(rc_fit)
        residual_density = cluster.density - \
            Gas_density_profile((rho0_fit, rc_fit), cluster.radius)

    pyplot.plot(cluster.radius+cluster.binsize/2, analytical_density,
            c=fit_colour, lw=4, label=label)  #, drawstyle="steps-mid")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$ [kpc]")
    ax.set_ylabel(r"$n_e$ [cm$^{-3}$]")
    ax.legend(loc=3, prop={'size':20})

    # Plot Residuals
    pyplot.sca(ax_r)

    pyplot.errorbar(cluster.radius+cluster.binsize/2, residual_density,
            yerr=cluster.density_std, c=fit_colour, drawstyle="steps-mid")
    ax_r.axhline(y=0, lw=2, ls="dashed", c="white")

    ax_r.set_xscale("log")
    # ax_r.set_yscale("log")

    ax.set_xlim(min(cluster.radius)-0.3, max(cluster.radius)+2000)
    ax_r.set_xlim(min(cluster.radius)-0.3, max(cluster.radius)+2000)

    if cluster.name == "cygA":
        ax_r.set_ylim(-0.03, 0.06)
    elif cluster.name == "cygB":
        ax_r.set_ylim(-0.0005, 0.002)

    # Fix for overlapping y-axis markers
    from matplotlib.ticker import MaxNLocator
    ax.tick_params(labelbottom='off')
    nbins = len(ax_r.get_yticklabels())
    ax_r.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

    if save:
        pyplot.savefig("out/density_betamodel_fit_{0}.png".format(cluster.name), dpi=600)
    # pyplot.show()


def cummulative_mass_profile(r, r_c, rho_0):
    return 4*numpy.pi*r_c**3*rho_0 * (r/r_c - numpy.arctan(r/r_c))


def get_mass_profile(cluster, result):
    ml_vals = result["x"]

    rho0_fit = ml_vals[0]  # | 1/units.cm**3
    rc_fit = ml_vals[1]  # | units.kpc
    if len(ml_vals) == 3:
        rcut_fit = ml_vals[2]  # | units.kpc
    if len(ml_vals) == 5:
        Rho0_Fac = ml_vals[3]
        Rc_Fac = ml_vals[4]

    kpc2cm = units.kpc(1).value_in(units.cm)
    rho0_analytical = rho0_fit * globals.mu * globals.m_p  # g/cm**3
    rc_analytical = rc_fit * kpc2cm  # cm

    analytical_radius = numpy.arange(min(cluster.radius), max(cluster.radius), 0.01) * kpc2cm
    analytical_mass = cummulative_mass_profile(analytical_radius, rc_analytical,  rho0_analytical)

    mass = globals.mu * globals.m_p * cluster.density * cluster.bin_volume
    volume = 4./3*numpy.pi*(cluster.radius*kpc2cm)**3  # cm cubed
    mass2 = globals.mu * globals.m_p * cluster.density * volume
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
    # amuse_plot.plot((cluster.radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
    #     (mass_summed | units.g).as_quantity_in(units.MSun),
    #     label="Martijn selfsummed")
    amuse_plot.plot((cluster.radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
        (mass.cumsum() | units.g).as_quantity_in(units.MSun),
        label="Martijn")
    # amuse_plot.plot((cluster.outer_radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
    #     (mass2_summed | units.g).as_quantity_in(units.MSun), label="Timo selfsummed")
    amuse_plot.plot((cluster.outer_radius*kpc2cm | units.cm).as_quantity_in(units.kpc),
        (mass2.cumsum() | units.g).as_quantity_in(units.MSun), label="Timo")
    pyplot.legend(loc=2)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    # pyplot.show()


if __name__ == "__main__":
    cygA = ObservedCluster("cygA")
    cygB = ObservedCluster("cygB")

    debug = False
    if debug:
        print "Debug information after parsing Martijn's files"
        print 80*"-"
        print cygA
        print cygB

        for property in ["bin_volume", "density", "pressure",
                "compton_y", "source_sb", "background_sb"]:
            print "property:", property
            simple_plot(cygA, property)
            simple_plot(cygB, property)
        print 80*"-"

    fit = True
    if fit:
        print "Obtaining parameters from observed density profiles"
        print 80*"-"
        discard_firstbins = False
        if discard_firstbins:
            cygA_fit = fit_betamodel_to_chandra(cygA, parm=[0.135, 27, 1.])
            cygA.radius = cygA.radius[4:]
            cygA.binsize = cygA.binsize[4:]
            cygA.density= cygA.density[4:]
            cygA.density_std = cygA.density_std[4:]
        cygA_fit = fit_betamodel_to_chandra(cygA, parm=[0.135, 27, 1.])
        cygB_fit = fit_betamodel_to_chandra(cygB, parm=[1., 1., 1.])
        cygB_fit = fit_betamodel_to_chandra(cygB, parm=[0.002, 200, 700])
        cygB_fit = fit_betamodel_to_chandra(cygB, parm=[0.002, 1.])
        print 80*"-"

    plot_fit = False
    if plot_fit:
        plot_fit_results(cygA, cygA_fit)
        plot_fit_results(cygB, cygB_fit)

    obtain_mass = False
    if obtain_mass:
        print "Obtaining observed mass profiles"
        print 80*"-"
        cygA_mass = get_mass_profile(cygA, cygA_fit)
        # cygB_mass = get_mass_profile(cygB, cygA_fit)
        print 80*"-"

    plot_mass = False
    if plot_mass:
        # get_mass_profile(cygA, cygA_mass)
        # get_mass_profile(cygB, cygB_mass)
        pass

    # pyplot.show()


    print "Obtaining r_200"
    print 80*"-"
    # cc = CosmologyCalculator(0)
    # print "rho_crit(z={0}) = {1:.4e}".format(cc.z, cc.rho_crit())
    # print "rho_crit(z={0}) = {1:.4e}".format(cluster.cc.z, cluster.cc.rho_crit())
    rhocrit200 = 200*cygA.cc.rho_crit()

    print "200*rhocrit =", rhocrit200
    mass_density = cygA.density * globals.mu * globals.m_p
    rho_average = (cygA.density / (4./3*numpy.pi*cygA.radius**3))
    pyplot.figure(figsize=(12, 9))
    pyplot.loglog(cygA.radius, rho_average)
    pyplot.show()

    print "Indicies of rho_average where rho_average > 200*rhocrit"
    print numpy.where(rho_average > rhocrit200)
    r_200 = cygA.radius[(numpy.abs(rho_average-rhocrit200)).argmin()]
    print "r200 =", r_200
    print 80*"-"







    import sys; sys.exit(0)


    # M_dm_below_r = mu * m_p * int n_e(r) dV
    analytical_radius = numpy.arange(min(cluster.radius), max(cluster.radius), 0.01)


    kpc2cm = 3.08568025e21
    A = 4./3*numpy.pi
    B = cluster.cc.kpc_DA*kpc2cm

    volume1, volume2, volume3, volume_martijn = [], [], [], []

    rho0_fit = result["x"][0] * globals.mu * globals.m_p
    r_c = result["x"][1]

    for i in range(len(cluster.bin_number)):
        volume1.append(A*(B*cluster.inner_radius[i])**3)
        volume2.append(A*(B*cluster.outer_radius[i])**3)
        volume3.append(A*(B*(cluster.outer_radius[i]+cluster.inner_radius[i])/2)**3)
        # print "Timo volume    = {0:.2e}".format(volume1)
        # print "Timo volume    = {0:.2e}".format(volume2)
        # print "Timo volume    = {0:.2e}".format(volume3)
        # print "Martijn volume = {0:.2e}\n".format(volume_martijn)
    # pyplot.figure(figsize=(12, 9))
    # pyplot.plot(cluster.radius, volume1, label="inner radius")
    # pyplot.plot(cluster.radius, volume2, label="outer radius")
    # pyplot.plot(cluster.radius, volume3, label="(inner+outer)/2")
    # pyplot.plot(cluster.radius, cluster.bin_volume, label="Martijn")
    # pyplot.plot(cluster.radius, cluster.bin_volume.cumsum(), label="Martijn, cumsum")
    # pyplot.legend(loc=2)
    # pyplot.gca().set_xscale("log")
    # pyplot.gca().set_yscale("log")

    print cluster.radius*kpc2cm
    import sys;sys.exit(0)

    M_gas_below_r = cummulative_mass_profile(cluster.radius, r_c, globals.mu*globals.m_p*rho0_fit)
    M_gas_below_r_1 = globals.mu * globals.m_p * cluster.density * cluster.bin_volume.cumsum() / 1.98892e33

    pyplot.figure(figsize=(12, 9))
    pyplot.plot(cluster.radius, M_gas_below_r, label="analytical")
    pyplot.plot(cluster.radius, M_gas_below_r_1, label="martijn")
    pyplot.legend(loc=6)
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    pyplot.show()

    mass = globals.mu*globals.m_p*cluster.density*cluster.bin_volume / 1.98892e33
    pyplot.loglog(cluster.radius, mass)
    pyplot.loglog(cluster.radius, mass.cumsum())
    pyplot.show()

    rho0_fit = result["x"][0]  # | units.g/units.cm**3
    rc_fit = result["x"][1]  # | units.kpc
    if len(result["x"]) == 3:
        rcut_fit = result["x"][2]  # | units.kpc
    if len(result["x"]) == 5:
        Rho0_Fac = result["x"][3]
        Rc_Fac = result["x"][4]

    rho_fit = ((rho0_fit | units.cm**-3) * globals.mu * (globals.m_p | units.g)).value_in(units.MSun/units.kpc**3)
    mass = Mass_profile(analytical_radius, rho0_fit, rc_fit, rcut_fit) * 1e10

    m200 = Mass_profile(rcut_fit, rho0_fit, rc_fit, rcut_fit)
    print m200
    print m200 *(1+0.17)
    print concentration_parameter(cluster, m200)

    pyplot.figure(figsize=(12, 9))
    mass2 = globals.mu*globals.m_p*cluster.density*cluster.bin_volume / 1.98892e33
    amuse_plot.loglog(cluster.radius | units.kpc, mass2 | units.MSun)
    amuse_plot.loglog(cluster.radius | units.kpc, mass2.cumsum() | units.MSun)
    amuse_plot.loglog(analytical_radius | units.kpc, mass | units.MSun)
    amuse_plot.xlabel(r"$r$")
    amuse_plot.ylabel(r"$M(<r)$")
    pyplot.xlim(10, 4500)
    # pyplot.ylim(9e9, 5e15)
    pyplot.show()
