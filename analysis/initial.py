"""
File: initial.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Mon Apr 18, 2016 02:25 pm
Last modified: Tue May 03, 2016 02:50 pm

Set up Galaxy Cluster initial conditions from the parsed Toycluster 2.0 output

Place particles in AMUSE datamodel to enable feeding it to Gadget-2 later on.
"""

import numpy
# Bigger fontsize is better *_*
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import new_quantity
from amuse.units.quantities import VectorQuantity
from amuse import datamodel
import amuse.plot as amuse_plot

from parser import Gadget2BinaryF77UnformattedType2Parser
from parser import Toycluster2RuntimeOutputParser


def amuse_nth_root(quant, n):
    """ Simply telling AMUSE e.q. quant**(1./3) breaks the units :-( """
    return new_quantity((quant.number)**(1./n), (quant.unit ** (1./n)).to_simple_form())


class Cluster(object):
    def __init__(self, icdir, snapdir, icfile=None):
        # Output of runToycluster writes files with these filename
        logfile="runToycluster.log"
        if icfile is None:
            icfile="IC_single_0"

        # Read runtime output of Toycluster 2.0. Parse runtime output
        self.toyclusterlog = Toycluster2RuntimeOutputParser(filename=icdir+logfile)
        self.set_toycluster2_values()


        self.raw_data = Gadget2BinaryF77UnformattedType2Parser(snapdir+icfile)
        # 1e10 because all masses are given in code units in cluster.par, which is set to 1e10 Msun
        self.M_gas = self.raw_data.Ngas * self.raw_data.massarr[0] * 1e10 | units.MSun
        self.M_dm = self.raw_data.Ndm * self.raw_data.massarr[1] * 1e10 | units.MSun
        self.gas, self.dm = self.place_ic_data_in_datamodel(self.raw_data)
        self.set_dm_density()
        self.set_mass()

    def set_dm_density(self):
        """ Dark Matter density radial profile rho(r). Given by Hernquist (1990) """
        rho_dm = self.M_dm/(2*numpy.pi) * self.a / (self.dm.r*(self.dm.r + self.a)**3)
        self.dm.rho = rho_dm.as_quantity_in(units.g/units.cm**3)

    def set_mass(self):
        # TODO: check with Julius that we can do this
        # Gadget-2 requires 'mass' attribute to be set for gas/dm datamodel
        # particleset when it is added to the Gadget2 instance.
        # In the IC file, however, no individual mass is set.
        # Also, the particles should have the same mass. This method is
        # probably wrong and should not be used. Hmmm...

        # self.dm.mass = self.dm_cummulative_mass(self.dm.r)
        # self.gas.mass = self.gas_cummulative_mass_double_beta(self.gas.r)

        # TODO: check if this method is correct
        self.gas.mass = self.raw_data.massarr[0] * 1e10 | units.MSun
        self.dm.mass = self.raw_data.massarr[1] * 1e10 | units.MSun

# TODO: check the methods below
    @property
    def rho_crit(self):
        """ Critical density of the Universe as a function of redshift """
        # rho_crit = 3 * (70 | units.km / units.s / units.Mpc)**2 / (8 * numpy.pi * constants.G)
        rho_crit = 3 * self.Hubble_of_z**2 / (8 * numpy.pi * constants.G)
        return rho_crit.as_quantity_in(units.g/units.cm**3)

    @property
    def Hubble_of_z(self):
        """ Hubble constant as a function of redshift """
        units.H0 = units.named('H0', 'km/s/Mpc', (1 | units.km / units.s / units.Mpc).to_unit())
        H_0 = self.h * 100 | units.H0

        return H_0 * numpy.sqrt(self.Omega_M*(1+self.z)**3 + self.Omega_Lambda)

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

    def set_toycluster2_values(self):
        self.model = self.toyclusterlog.halo0_model
        self.rgas = self.toyclusterlog.halo0_rgas
        self.rdm  = self.toyclusterlog.halo0_rdm
        self.qmax = self.toyclusterlog.halo0_qmax
        self.Mass = self.toyclusterlog.halo0_Mass
        self.Mass_in_DM = self.toyclusterlog.halo0_Mass_in_DM
        self.Mass_in_gas = self.toyclusterlog.halo0_Mass_in_gas
        self.Mass_in_R200 = self.toyclusterlog.halo0_Mass_in_R200
        self.c_nfw = self.toyclusterlog.halo0_c_nfw
        self.R200 = self.toyclusterlog.halo0_R200
        self.a = self.toyclusterlog.halo0_a_hernquist
        self.rho0gas = self.toyclusterlog.halo0_rho0gas_cgs
        self.rho0gas_gadget = self.toyclusterlog.halo0_rho0gas_gadget
        self.beta = self.toyclusterlog.halo0_beta
        self.rc = self.toyclusterlog.halo0_rc
        self.rcut = self.toyclusterlog.halo0_R200
        self.R500  = self.toyclusterlog.halo0_R500
        self.bf_200 = self.toyclusterlog.halo0_bf_200
        self.bf_500 = self.toyclusterlog.halo0_bf_500

        self.Omega_M = self.toyclusterlog.system_Omega_M
        self.Omega_Lambda = 1 - self.Omega_M
        self.h = self.toyclusterlog.system_H_over_100
        self.z = self.toyclusterlog.system_z

    def place_ic_data_in_datamodel(self, ic_data):
        # Split gas and DM in the IC position arrays.
        gaspos = ic_data.pos[0:ic_data.Ngas]
        dmpos = ic_data.pos[ic_data.Ngas:ic_data.N]

        # Split positions into x, y, z; calculate and set radius.
        xgas = gaspos[:,0] - ic_data.boxSize/2
        ygas = gaspos[:,1] - ic_data.boxSize/2
        zgas = gaspos[:,2] - ic_data.boxSize/2
        rgas = numpy.sqrt(xgas**2+ygas**2+zgas**2)
        xdm = dmpos[:,0] - ic_data.boxSize/2
        ydm = dmpos[:,1] - ic_data.boxSize/2
        zdm = dmpos[:,2] - ic_data.boxSize/2
        rdm = numpy.sqrt(xdm**2+ydm**2+zdm**2)

        # Split gas and DM in the IC velocity arrays.
        gasvel = ic_data.vel[0:ic_data.Ngas]
        dmvel = ic_data.vel[ic_data.Ngas:ic_data.N]

        # Split velocities into vx, vy, vz.
        vxgas = gasvel[:,0]
        vygas = gasvel[:,1]
        vzgas = gasvel[:,2]
        vxdm = dmvel[:,0]
        vydm = dmvel[:,1]
        vzdm = dmvel[:,2]

        # Split gas and DM in the IC id arrays.
        gasids = ic_data.ids[0:ic_data.Ngas]
        dmids = ic_data.ids[ic_data.Ngas:ic_data.N]

        # Set RHO and RHOm
        rhogas = ic_data.rho#[0:ic_data.Ngas]
        # rhomgas = ic_data.rhom#[0:ic_data.Ngas]

        # Set up datamodel
        gas = datamodel.Particles(ic_data.Ngas)
        gas.x = xgas | units.kpc
        gas.y = ygas | units.kpc
        gas.z = zgas | units.kpc
        gas.r = rgas | units.kpc
        gas.vx = vxgas | units.kms
        gas.vy = vygas | units.kms
        gas.vz = vzgas | units.kms
        gas.rho = (rhogas * 1e10 * ic_data.hubbleParam**(2) | (units.MSun / units.kpc**3))\
            .as_quantity_in(units.g/units.cm**3)
        # gas.rhom = (rhomgas * 1e10 * ic_data.hubbleParam**(2) | (units.MSun / units.kpc**3))\
        #    .as_quantity_in(units.g/units.cm**3)

        dm = datamodel.Particles(ic_data.Ndm)
        dm.x = xdm | units.kpc
        dm.y = ydm | units.kpc
        dm.z = zdm | units.kpc
        dm.r = rdm | units.kpc
        dm.vx = vxdm | units.kms
        dm.vy = vydm | units.kms
        dm.vz = vzdm | units.kms

        return gas, dm

    def perform_sanity_checks(self):

        # Let's see what the radius of the gas and dm looks like

        # Donnert (2014), Sec. 3: "The gas profile of one cluster is sampled to a maximal radius $R_{\rm max}$, which is half the box size"
        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.gas.r, bins=int(numpy.sqrt(self.raw_data.Ngas)), label="Gas")
        amuse_plot.ylabel(r'$N$')
        amuse_plot.xlabel(r'$r$')
        pyplot.legend()

        # Donnert (2014), Sec. 3: "The clusters DM profile is sampled up to infinity"
        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.dm.r, bins=int(numpy.sqrt(self.raw_data.Ndm)), label="DM")
        pyplot.legend()
        amuse_plot.ylabel(r'$N$')
        amuse_plot.xlabel(r'$r$')

        # Check the velocities of the gas (should only be zero) and of the dm (should follow Hernquist DF)
        for attr in ["vx", "vy", "vz"]:
            # NB, turn into numpy array using value_in method! VectorQuantity has no nonzero method
            v = getattr(self.gas, attr).value_in(units.kms)
            if len(v.nonzero()[0]) != 0:
                print "Error: gas {0} has nonzero values".format(attr)
            else:
                print "Passed: gas {0} only has nonzero values".format(attr)
            v = getattr(self.dm, attr).value_in(units.kms)
            if len(v.nonzero()[0]) == self.raw_data.Ndm:
                print "Passed: dm {0} has no nonzero values".format(attr)
            else:
                print "Error: dm {0} has nonzero values".format(attr)

        # Look at the gas velocity distribution (should only be zero)
        # Donnert (2014), Sec. 3: "The initial velocity of the gas particles is set to zero"
        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.gas.vx, bins=int(numpy.sqrt(self.raw_data.Ngas)))
        pyplot.xlim(-1, 1)

        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.gas.vy, bins=int(numpy.sqrt(self.raw_data.Ngas)))
        pyplot.xlim(-1, 1)

        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.gas.vz, bins=int(numpy.sqrt(self.raw_data.Ngas)))
        pyplot.xlim(-1, 1)

        # Look at the dm velocity distribution (should follow Hernquist DF)
        # Donnert (2014), Sec. 3 "The velocity of the DM component is satisfied so the Hernquist DF, equation (5), is satisfied"
        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.dm.vx, bins=int(numpy.sqrt(self.raw_data.Ndm)))

        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.dm.vy, bins=int(numpy.sqrt(self.raw_data.Ndm)))

        pyplot.figure(figsize=(9, 9))
        amuse_plot.hist(self.dm.vz, bins=int(numpy.sqrt(self.raw_data.Ndm)))
        # pyplot.show()
