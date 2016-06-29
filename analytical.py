"""
File: analytical.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Mon May 16, 2016 05:23 pm
Last modified: Fri Jun 24, 2016 02:35 pm

Set up analytical cluster following Donnert (2014) and Donnert et al. (2016, in prep)
"""

import numpy
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import new_quantity
from amuse.units.quantities import VectorQuantity
from amuse import datamodel
import amuse.plot as amuse_plot

from macro import amuse_nth_root
from cosmology import CosmologyCalculator
from parser import parse_toycluster_parms

class AnalyticalCluster(object):
    def __init__(self):
        self.cc = CosmologyCalculator()

        self.set_density()
        self.set_mass()
        self.set_temperature()
        self.set_density()

    def set_dm_density(self):
        """ Dark Matter density radial profile rho(r). Given by Hernquist (1990) """
        rho_dm = self.M_dm/(2*numpy.pi) * self.a / (self.dm.r*(self.dm.r + self.a)**3)
        self.dm.rho = rho_dm.as_quantity_in(units.g/units.cm**3)

    def set_mass(self):
        pass

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



if __name__ == "__main__":
    print "Parsing Toycluster Parameter file"
    print 80*'-'
    parms = parse_toycluster_parms("toycluster.par")
    for key, value in parms.iteritems():
        print key, "=", value

    Xm = parms["Mass_Ratio"]
    Xm = min(Xm, 1/Xm)

    massA = parms["Mtotal"] * Xm
    rcoreA = parms["rc_1"]
    c_nfwA = parms["c_nfw_1"]
    v_comA = parms["v_com_0"]
    print "Cluster A"
    print "Mass       =", massA
    print "Rcore      =", rcoreA
    print "c_nfw      =", c_nfwA
    print "v_com      =", v_comA

    print "Cluster B"
    print "Mass       =", massB
    print "Rcore      =", rcoreB
    print "c_nfw      =", c_nfwB
    print "v_com      =", v_comB
    print 80*'-'
