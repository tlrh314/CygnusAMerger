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
from cosmology import CosmologyCalculator
from toyclusterparamsparser import parse_toycluster_parms

class AnalyticalCluster(object):
    def __init__(self):
        pass

    def set_dm_density(self):
        """ Dark Matter density radial profile rho(r). Given by Hernquist (1990) """
        rho_dm = self.M_dm/(2*numpy.pi) * self.a / (self.dm.r*(self.dm.r + self.a)**3)
        self.dm.rho = rho_dm.as_quantity_in(units.g/units.cm**3)

    def dm_density(self, r):
        """ Dark Matter density radial profile rho(r). Given by Hernquist (1990) """
        rho_dm = self.M_dm/(2*numpy.pi) * self.a / (r*(r + self.a)**3)
        return rho_dm.as_quantity_in(units.g/units.cm**3)

    def dm_cummulative_mass(self, r):
        """ Dark Matter cummulative mass profile M_DM(<r). Given by Hernquist (1990) """
        M_dm_below_r = self.M_dm * r**2 / (r + self.a)**2
        return M_dm_below_r.as_quantity_in(units.MSun)

    def average_dm_density(self, r):
        rho_average = (self.dm_cummulative_mass(r) / (4./3*numpy.pi*r**3))
        return rho_average.as_quantity_in(units.g/units.cm**3)

    def gas_density_beta(self, r):
        # Gas profile
        # "Should any other value of beta be used, this step would involve the confluent
        # hypergeometric function, which usually makes the analytical solution of the
        # hydrostatic equation impossible. Another elegant solution would be beta=1,
        # which however is unfavoured from observations" - Donnert in MNRAS 438 2014
        # self.beta = 2./3  # Mastropietro & Burkert (2008)
        rho_gas_beta = self.rho0gas * (1 + r**2/self.rc**2)**(-3*self.beta/2.)
        return rho_gas_beta.as_quantity_in(units.g/units.cm**3)

    def gas_density_double_beta(self, r):
        """ Double beta model. The additional cutoff ensures rho_gas < rho_dm
            for large r (r > r_200?). See Donnert et al. (2016) in prep. """
        # TODO: understand exactly why double beta profile!
        rho_gas = self.rho0gas / (1 + (r/self.rc)**2) / (1 + (r/self.rcut)**3 * (r/self.rcut))
        return rho_gas.as_quantity_in(units.g/units.cm**3)

    def gas_cummulative_mass_beta(self, r):
        """ Analytical (single) beta model: cummulative mass profile """
        M_gas_below_r_beta = 4*numpy.pi*self.rc**3*self.rho0gas * (r/self.rc - numpy.arctan(r/self.rc))
        return M_gas_below_r_beta.as_quantity_in(units.MSun)

    def gas_cummulative_mass_double_beta(self, r):
        """ Analytical double beta model: cummulative mass profile
            See Donnert et al. (2016) in prep for details. """
        # Copied from setup.c. Is this faster?
        r2 = r**2
        rc = self.rc
        rc2 = self.rc**2
        rcut = self.rcut
        rcut2 = self.rcut**2
        sqrt2 = numpy.sqrt(2)

        A = (rc2 - rcut2)*(numpy.log((rcut2 - sqrt2*rcut*r + r2).value_in(units.kpc**2)) \
            - numpy.log((rcut2 + sqrt2*rcut*r + r2).value_in(units.kpc**2)))
        Bplus = 2*(rc2 + rcut2)*numpy.arctan(1 + sqrt2*r/rcut)
        Bmin = 2*(rc2 + rcut2)*numpy.arctan(1 - sqrt2*r/rcut)

        # NB the paper is slightly different: equation 2 does not contain 4*pi*rho
        M_gas_below_r = (4 * numpy.pi * self.rho0gas) *\
            rc2*rcut**3/(8*(rcut**4+rc**4)) *\
            (sqrt2*(A - Bmin + Bplus) - 8*rc*rcut*numpy.arctan(r/rc))
        return M_gas_below_r.as_quantity_in(units.MSun)

    def gas_temperature(self, r):
        F_1 = numpy.pi**2/(8*self.rc) - numpy.arctan(r/self.rc)**2/(2*self.rc) -\
            numpy.arctan(r/self.rc)/r

        T_r_gas = constants.G*self.mu*constants.proton_mass/constants.kB *\
            (1 + r**2/self.rc**2) * (4*numpy.pi*self.rc**3*self.rho0gas*F_1)
        return T_r_gas.as_quantity_in(units.K)

    def dm_temperature(self, r):
        F_0 = self.rc/(self.a**2 + self.rc**2)**2 * (numpy.pi/2*(self.a**2 - self.rc**2) +\
            self.r_c*(self.a**2 + self.rc**2)/(self.a + r) - (self.a**2 - self.rc**2)*\
            numpy.arctan(r/self.rc) - self.rc*self.a*numpy.log((self.a + r)**2/(r**2 + self.rc**2)))
        T_r_dm = constants.G*self.mu*constants.proton_mass/constants.kB * (1 + r**2/self.rc**2) *\
            (self.M_dm*F_0)
        return T_r_dm.as_quantity_in(units.K)

    def temperature(self, r):
        return (self.gas_temperature(r) + self.dm_temperature(r)).as_quantity_in(units.K)

    def characteristic_temperature_analytically(self):
        # Charecteristic temperature
        T_c = 2*constants.G*self.mu*constants.proton_mass/constants.kB * \
            (self.M_dm * self.rc**2/(self.a**2 + self.rc**2)**2 * (numpy.pi/(4*self.rc) * \
            (self.a**2 - self.rc**2) + (self.a**2 + self.rc**2)/(self.a + self.rc) - \
            self.a*numpy.log((self.a + self.rc)**2/(2*self.rc**2))) \
            + numpy.pi**2*self.rc**2*self.rho_0*(3*numpy.pi/8 - 1))
        return T_c.as_quantity_in(units.K)

    def gas_pressure(self, r):
        P_gas = constants.kB*self.temperature(r)*self.gas_density(r)/\
            (self.mu*constants.proton_mass)
        return P_gas.as_quantity_in(units.Pa)

    def dm_pressure(self, r):
        P_dm = constants.kB*self.temperature(r)*self.dm_density(r)/\
            (self.mu*constants.proton_mass)
        return P_dm.as_quantity_in(units.Pa)

    def find_characteristic_r(self, mass=None, rscale=200):
        """ We assume a spherically symmetric ball of gas/dark matter """
        if not mass:
            print "Assuming r_200 should be found "
            mass = self.M_200
        r_cubed = (0.75/numpy.pi * mass/(rscale*self.rho_crit))
        return amuse_nth_root(r_cubed, 3).as_quantity_in(units.kpc)

    def find_p500_analytically(self):
        """ Calculate characteristic pressure at r_500 using
        the analytical formula given by Arnaud et al. (2010) in Appendix A """
        P_500 = 3/(8*numpy.pi) * ((500./2)*constants.G**(-1./4)*self.Hubble_of_z**2)**(4./3) *\
            self.mu/self.mu_e * self.f_b * self.M_500**(2./3)

        return P_500.as_quantity_in(units.Pa)

# TODO: end tocheck




class Toycluster(object):
    """ Numerical Representation of Clusters """
    def __init__(self):
        self.parms = parse_toycluster_parms("toycluster.par")

        M200 = self.parms["Mtotal"]
        Xm = self.parms["Mass_Ratio"]
        redshift = self.parms["Redshift"]
        h_100 = self.parms["h_100"]
        bf = self.parms["bf"]

        M200_A = M200 / (1+Xm)
        M200_B = M200 - M200_A

        M200_A_gas = M200_A / (1+bf)
        M200_B_gas = M200_B / (1+bf)

        M200_A_dm = M200 - M200_A_gas
        M200_B_dm = M200 - M200_B_gas

        # print M200_A
        # print M200_A_dm
        # print M200_A_gas
        # print concentration_parameter(M200_A)
        # print M200_B
        # print M200_B_dm
        # print M200_B_gas
        # print concentration_parameter(M200_B)


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

        m_p = constants.proton_mass.value_in(units.g)
        mu = 0.17
        cc = CosmologyCalculator()
        arcsec2kpc = cc.kpc_DA # | units.kpc

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
        m_p = constants.proton_mass.value_in(units.g)
        mu = 0.17
        density = m_p * mu * cluster.density
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
    pyplot.show()


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


def concentration_parameter(gas_mass):
    cc = CosmologyCalculator()
    A = 5.74
    B = -0.097
    C = -0.47
    Mpivot = 2e12 / cc.h

    UnitMass_in_g = 1.989e43  # 1e10 Msun
    Msol2cgs = 1.98892e33

    bf = 0.17
    total_mass = gas_mass * (1+1/bf)  #Mgas+Mdm, Mdm = Mgas/bf

    # total mass has to be in 1e10 Msun

    mass = total_mass*UnitMass_in_g/Msol2cgs

    c_NFW = A * pow(mass/Mpivot, B) * pow(1+cc.z, C)

    return c_NFW


# Define the statistical model, in this case we shall use a chi-squared distribution, assuming normality in the errors
def stat(parm, x, y, dy):
    ymodel = Gas_density_profile(parm, x)
    chisq = numpy.sum((y - ymodel)**2 / dy**2)
    return(chisq)


def fit_betamodel_to_chandra(cluster, parm=[1., 1., 1.], double=False):
    """ Fit betamodel to cluster observation, or
        double betamodel when len(parm)==5. """
    analytical_radius = numpy.arange(min(cluster.radius), max(cluster.radius), 0.01)

    # set_printing_strategy("custom", preferred_units = [units.kpc, units.cm**-3*units.g])

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

    # Plot fit results
    rho0_fit = result["x"][0]  # | units.g/units.cm**3
    rc_fit = result["x"][1]  # | units.kpc
    if len(parm) == 3:
        rcut_fit = result["x"][2]  # | units.kpc
        analytical_density = Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
        label = r"cut-off $\beta$-model:"+"\n\t"+r"$n_{{e,0}} \,$ = {0:.2e}".format(rho0_fit)\
            +"\n\t"+r"$r_c \,$ = {0:.2f}".format(rc_fit)+"\n\t"+\
            r"$r_{{\rm cut}}$ = {0:.2f}".format(rcut_fit)
        residual_density = cluster.density - \
            Gas_density_profile((rho0_fit, rc_fit, rcut_fit), cluster.radius)
    elif len(parm) == 5:
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


    pyplot.savefig("out/density_betamodel_fit_{0}.png".format(cluster.name), dpi=600)
    # pyplot.show()

    # set_printing_strategy("cgs")

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


def get_mass_profile(cluster):
    # M_dm_below_r = mu * m_p * int n_e(r) dV
    m_p = constants.proton_mass.value_in(units.g)
    mu = 0.17

    cc = CosmologyCalculator(0)
    print "rho_crit(z={0}) = {1:.4e}".format(cc.z, cc.rho_crit())
    cc = CosmologyCalculator()
    print "rho_crit(z={0}) = {1:.4e}".format(cc.z, cc.rho_crit())

    rhocrit200 = 200*cc.rho_crit()

    return

    for i in range(len(cluster.bin_number)):
        volume = 4./3 * numpy.pi * (cluster.outer_radius[i] * arcsec2kpc * 3.08568025e21)**3
        volume_martijn = cluster.bin_volume[i]

        print "Timo volume    = {0:.2e}".format(volume)
        print "Martijn volume = {0:.2e}\n".format(volume_martijn)

    mass = mu*m_p*cluster.density*cluster.bin_volume / 1.98892e33
    pyplot.loglog(cluster.radius, mass)
    pyplot.loglog(cluster.radius, mass.cumsum())
    pyplot.show()

    m_p = constants.proton_mass
    mu = 0.17
    rho_fit = ((rho0_fit | units.cm**-3) * mu * m_p).value_in(units.MSun/units.kpc**3)
    mass = Mass_profile(analytical_radius, rho0_fit, rc_fit, rcut_fit) * 1e10

    m200 = Mass_profile(rcut_fit, rho0_fit, rc_fit, rcut_fit)
    print m200
    print m200 *(1+0.17)
    print concentration_parameter(m200)

    pyplot.figure(figsize=(12, 9))
    mass2 = mu*m_p.value_in(units.g)*cluster.density*cluster.bin_volume / 1.98892e33
    amuse_plot.loglog(cluster.radius | units.kpc, mass2 | units.MSun)
    amuse_plot.loglog(cluster.radius | units.kpc, mass2.cumsum() | units.MSun)
    amuse_plot.loglog(analytical_radius | units.kpc, mass | units.MSun)
    amuse_plot.xlabel(r"$r$")
    amuse_plot.ylabel(r"$M(<r)$")
    pyplot.xlim(10, 4500)
    # pyplot.ylim(9e9, 5e15)
    pyplot.show()


    # r_200 = ClusterA.r[(numpy.abs(rho_average-200*ClusterA.rho_crit)).argmin()]


if __name__ == "__main__":
    cygA = ObservedCluster("cygA")
    cygB = ObservedCluster("cygB")

    debug = False
    if debug:
        print cygA
        print cygB

        for property in ["bin_volume", "density", "pressure",
                "compton_y", "source_sb", "background_sb"]:
            print "property:", property
            simple_plot(cygA, property)
            simple_plot(cygB, property)

    fit = True
    if fit:
        # cygA_fit = fit_betamodel_to_chandra(cygA, parm=[0.135, 27, 1.])
        # cygA.radius = cygA.radius[4:]
        # cygA.binsize = cygA.binsize[4:]
        # cygA.density= cygA.density[4:]
        # cygA.density_std = cygA.density_std[4:]
        # cygA_fit = fit_betamodel_to_chandra(cygA, parm=[0.135, 27, 1.])
        # cygB_fit = fit_betamodel_to_chandra(cygB, parm=[1., 1., 1.])
        # cygB_fit = fit_betamodel_to_chandra(cygB, parm=[0.002, 200, 700])
        cygB_fit = fit_betamodel_to_chandra(cygB, parm=[0.002, 1.])

    obtain_mass = False
    if obtain_mass:
        get_mass_profile(cygA)

    pyplot.show()
