import numpy
import pandas
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import VectorQuantity
from amuse import datamodel
import amuse.plot as amuse_plot

from macro import *
from cosmology import CosmologyCalculator
from ioparser import parse_toycluster_parms
from ioparser import Toycluster2RuntimeOutputParser
from ioparser import Gadget2BinaryF77UnformattedType2Parser
import convert


class ObservedCluster(object):
    """ Observed Cygnus A-Cygnus B situation """
    def __init__(self, name, oldICs=False):
        """ Generate cluster instance with volume, pressure, density,
            Compton Y for bins with given inner and outer radii.

            Data obtained by M. N. de Vries from 800 ksec Chandra data

            @param oldICs: True --> data prior to background (thus rho)

            TODO: run with 900 ksec including pointing at CygB :-)! """

        self.name = name
        self.oldICs = oldICs

        if self.oldICs:
            density_file = "data/{0}_1T_fixnH_pressureprofile_OLD.dat".format(name)
            radius_file = "data/{0}_sn100_sbprofile.dat".format(name)
        else:
            density_file = "data/{0}_1T_fixnH_pressureprofile.dat".format(name)
            radius_file = "data/{0}_sn100_sbprofile.dat".format(name)

        if self.name == "cygA":
            z = 0.0562
        if self.name == "cygB":
            z = 0.070
        self.cc = CosmologyCalculator(z)

        self.parse_data(density_file, radius_file)

        arcsec2kpc = self.cc.kpc_DA # | units.kpc

        self.radius = (self.inner_radius+self.outer_radius)/2 * arcsec2kpc
        self.binsize = (self.outer_radius - self.inner_radius) * arcsec2kpc

    def __str__(self):
        tmp = "" + str(self.name)

        tmp += "\nbin_number\n" + str(self.bin_number)
        tmp += "\n\nbin_volume\n" + str(self.bin_volume)
        tmp += "\n\nnumber_density\n" + str(self.number_density)
        tmp += "\n\nnumber_density_std\n" + str(self.number_density_std)
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
        # TODO: old datafile also has column Analytical bin volume. What this?
        self.bin_volume = raw[" Bin volume (cm^3) "].as_matrix()
        self.number_density = raw["   Density (cm^-3) "].as_matrix()
        self.number_density_std = raw["     Sigma density "].as_matrix()
        self.density = convert.ne_to_rho(self.number_density, self.cc.z)
        self.density_std = convert.ne_to_rho(self.number_density_std, self.cc.z)
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


class NumericalCluster(object):
    def __init__(self, icdir, snapdir, logfile=None, icfile=None):
        # Output of runToycluster writes files with these filename
        if logfile is None:
            logfile="runToycluster.log"
        if icfile is None:
            icfile="IC_single_0"

        # Read runtime output of Toycluster 2.0. Parse runtime output
        self.toyclusterlog = Toycluster2RuntimeOutputParser(filename=icdir+logfile)
        self.raw_data = Gadget2BinaryF77UnformattedType2Parser(snapdir+icfile)

        # if the mass ratio is 0.0 we only have one cluster in the box
        if (-2**-14 < self.toyclusterlog.systemsetup['Mass_Ratio'] < 2**-14):
            self.set_toycluster2_values(i=0)
            # 1e10 because all masses are given in code units in cluster.par, which is set to 1e10 Msun
            # Only makes sense if Mass_Ratio = 0.0
            self.M_gas = self.raw_data.Ngas * self.raw_data.massarr[0] * 1e10 | units.MSun
            self.M_dm = self.raw_data.Ndm * self.raw_data.massarr[1] * 1e10 | units.MSun

            """ This actually also works if Mass_Ratio != 0.0
            NB then the gas and dm contain datamodels with all
            gas/dm particles.

            Elsewhere we explicitly split both particles,
            then use this class for a single cluster """

            # TODO: this step is already rather slow fow 2e6 particles...
            self.gas, self.dm = self.place_ic_data_in_datamodel(self.raw_data)

            # Only makes sense if Mass_Ratio = 0.0
            # Use these functions in density_stability
            # self.get_gas_mass_via_density()
            # self.get_dm_mass_via_number_density()
            # self.set_dm_density()

    def set_toycluster2_values(self, i=0):
        self.model = self.toyclusterlog.halosetup[i]['model']
        self.rgas = self.toyclusterlog.halosetup[i]['rgas']
        self.rdm  = self.toyclusterlog.halosetup[i]['rdm']
        self.qmax = self.toyclusterlog.halosetup[i]['qmax']
        self.Mass = self.toyclusterlog.halosetup[i]['Mass']
        self.Mass_in_DM = self.toyclusterlog.halosetup[i]['Mass_in_DM']
        self.Mass_in_gas = self.toyclusterlog.halosetup[i]['Mass_in_gas']
        self.Mass_in_R200 = self.toyclusterlog.halosetup[i]['Mass_in_R200']
        self.c_nfw = self.toyclusterlog.halosetup[i]['c_nfw']
        self.R200 = self.toyclusterlog.halosetup[i]['R200']
        self.a = self.toyclusterlog.halosetup[i]['a_hernquist']
        self.rho0gas = self.toyclusterlog.halosetup[i]['rho0gas_cgs']
        self.rho0gas_gadget = self.toyclusterlog.halosetup[i]['rho0gas_gadget']
        self.beta = self.toyclusterlog.halosetup[i]['beta']
        self.rc = self.toyclusterlog.halosetup[i]['rc']
        self.rcut = self.toyclusterlog.halosetup[i]['R200']
        self.R500  = self.toyclusterlog.halosetup[i]['R500']
        self.bf_200 = self.toyclusterlog.halosetup[i]['bf_200']
        self.bf_500 = self.toyclusterlog.halosetup[i]['bf_500']

        self.Omega_M = self.toyclusterlog.systemat['Omega_M']
        self.Omega_Lambda = 1 - self.Omega_M
        self.h = self.toyclusterlog.systemat['H_over_100']
        self.z = self.toyclusterlog.systemat['z']

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
        gas.h = ic_data.hsml
        gas.u = ic_data.u

        # Bug in Toycluster's unit.c line 35: "* p2(0.7)"
        # Removed in commit e88b863319e969f4a15765baad447de9e3571d7e
        #gas.rho = (rhogas * 1e10 * ic_data.hubbleParam**(2) | (units.MSun / units.kpc**3))\
        #    .as_quantity_in(units.g/units.cm**3)

        gas.rho = (rhogas * 1e10  | (units.MSun / units.kpc**3))\
            .as_quantity_in(units.g/units.cm**3)

        dm = datamodel.Particles(ic_data.Ndm)
        dm.x = xdm | units.kpc
        dm.y = ydm | units.kpc
        dm.z = zdm | units.kpc
        dm.r = rdm | units.kpc
        dm.vx = vxdm | units.kms
        dm.vy = vydm | units.kms
        dm.vz = vzdm | units.kms

        return gas, dm

    def get_gas_mass_via_density(self, DESNNGB=50):
        """ Kernel weigthed sph density beautifully fits analytical solution <3
            DM has no density, so we obtain it trough counting particles.

            @param DESNNGB: 50 for Gadget-2 B-spline, 295 for toycluster WC6"""

        gas_i = numpy.argsort(self.gas.r.value_in(units.kpc))
        gas_r = self.gas.r[gas_i].value_in(units.kpc)
        gas_h = self.gas.h[gas_i]
        gas_rho = (self.gas.rho[gas_i]).value_in(units.MSun/units.kpc**3)

        # Price (2012) equation 11: Mtot = 4/3 pi R_kern**3 rho.
        # Rkern**3 = hsml**3/295 (Toycluster: Wendland C6 with 295 neighbours)
        # For ICs!
        # gas_mass = (4./3*numpy.pi*(gas_h**3/295)*gas_rho)
        # For Gadget-2 (TODO: check Kernel, also check Nngb)
        gas_mass = (4./3*numpy.pi*(gas_h**3/DESNNGB)*gas_rho)

        self.gas_radii = gas_r | units.kpc
        self.M_gas_below_r = gas_mass.cumsum() | units.MSun

    def get_dm_mass_via_number_density(self):
        """ Count particles <r (= number density). Obtain DM mass from it """

        # TODO: find out if this method can also be used for sph mass
        # should we take into account the kernel, and if so how?

        print "Counting particles for which radii < r to obtain M(<r)"

        radii = numpy.arange(0, 1e4, 10)
        N = len(radii)

        particles = numpy.zeros(N)
        dr = radii[1] - radii[0]  # Because evenly spaced
        for i, r in enumerate(radii):
            particles[i] = ((numpy.where(self.dm.r.value_in(units.kpc) < r)[0]).size)
            if i==(N-1) or i%100 == 0:
                print_progressbar(i, N)

            # print i, r, particles[i]

        particles_in_shell = numpy.zeros(len(particles))
        for i in range(1, len(particles)):
            particles_in_shell[i-1] = particles[i] - particles[i-1]

        self.dm_volume = 4 * numpy.pi * radii**2 * dr
        self.dm_particles_in_shell = particles_in_shell
        # density = particles_in_shell/(self.raw_data.Ndm*volume*self.M_dm.number)
        # density *= globals.mu  # TODO: Why ??
        # Using volume = 4/3 pi r**3
        # M_dm_below_r = numpy.zeros(len(radii))
        # N = len(radii)
        # for i, r in enumerate(radii):
        #     M_dm_below_r[i] = ((numpy.where(self.dm.r < r)[0]).size)
        #     if i==(N-1) or i%100 == 0:
        #         globals.print_progressbar(i, N)

        # M_dm_below_r *= (self.M_dm/self.raw_data.Ndm)

        # Works
        M_dm_below_r = particles * self.M_dm/self.raw_data.Ndm

        # Also works
        # DM_Part_mass = self.toyclusterlog.systemsetup['DM_Part_mass']
        # M_dm_below_r = particles * DM_Part_mass
        self.dm_radii = radii | units.kpc
        self.M_dm_below_r = M_dm_below_r

    def set_dm_density(self):
        # Works, but what is up with the 0.1 ? O_o
        # self.rho_dm_below_r = 0.1*(self.dm_particles_in_shell/(self.raw_data.Ndm*self.dm_volume*self.M_dm.number)) | units.g/units.cm**3

        # Does not work
        # self.rho_dm_below_r = self.M_dm_below_r / (self.dm_radii**3)

        # Does not work
        # volume = self.dm_volume.cumsum()
        # self.rho_dm_below_r = self.M_dm_below_r / (volume | units.kpc**3)

        # Works, seems good enough alright
        # DM_Part_mass = self.toyclusterlog.systemsetup['DM_Part_mass']
        # self.rho_dm_below_r = DM_Part_mass * self.dm_particles_in_shell / (self.dm_volume | units.kpc**3)

        self.rho_dm_below_r = self.M_dm * (self.dm_particles_in_shell/self.raw_data.Ndm) / (self.dm_volume | units.kpc**3)

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

class AnalyticalCluster(object):
    """ Set up an analytical cluster based on the Chandra observed density """
    def __init__(self, parms, dm_parms=None, radius=None, z=0):
        """ parms = ne0, rc, [rcut], [ne0_fac, rc_fac]
            dm_parms = M_DM, a  """

        # TODO: implement without need for dm_parms

        if radius is None:
            self.radius = VectorQuantity.arange(units.kpc(1), units.kpc(1e5), units.parsec(100))
        else:
            if type(radius) == numpy.ndarray:
                self.radius = radius | units.kpc
            else:
                self.radius = radius
        self.z = z
        self.ne0 = parms[0]
        self.rho0 = convert.ne_to_rho(self.ne0, self.z) | units.g/units.cm**3
        self.ne0 = self.ne0 | units.cm**-3
        self.rc = parms[1] | units.kpc
        self.model = 0
        self.rcut = None
        self.ne0_cc = None
        self.rho0_cc = None
        self.rc_cc = None
        if len(parms) == 3:
            self.model = 1
            self.rcut = parms[2] | units.kpc
        if len(parms) == 5:
            self.model = 2
            ne0_fac = parms[3]
            rc_fac = parms[4]
            self.ne0_cc = ne0 * ne0_fac
            self.rho0_cc = self.ne0_cc * globals.mu * (globals.m_p | units.g)
            self.rc_cc = rc / rc_fac

        modelnames = {0: r"$\beta (2/3)$-model",
                      1: r"cut-off $\beta (2/3)$-model",
                      2: r"cut-off double $\beta (2/3)$-model"}
        self.modelname = modelnames[self.model]

        if dm_parms:
            self.M_dm = dm_parms[0]
            self.a = dm_parms[1]

        # self.get_r200()

    def gas_density(self, r=None):
        """ For details of beta-model see Cavaliere & Fusco-Femiano (1978).

            model    name                   reference
            0        beta-model             Donnert 2014)
            1        cut-off beta           Donnert et al. (2016, in prep)
            2        cut-off double beta    Donnert et al. (2016, in prep)
        """

        if r is None:
            r = self.radius

        beta = 2./3  # Mastropietro & Burkert (2008)
        rho_gas = self.rho0 * (1 + p2(r/self.rc))**(-3*beta/2.)  # model >= 0
        if self.model >= 1:
            rho_gas /= (1 + p3(r/self.rcut) * (r/self.rcut))
        if self.model == 2:
            rho_gas += self.rho0_cc / (1 + p2(r/rc_cc)) / (1 + p3(r/rcut) * (r/rcut))

        return rho_gas.as_quantity_in(units.g/units.cm**3)

    def gas_number_density(self, r=None):
        return (convert.rho_to_ne(self.gas_density(r).value_in(units.g/units.cm**3),
            self.z) | (1/units.cm**3))

    def gas_mass(self, r=None):
        """ return M(<r) """

        if not r:
            r = self.radius

        if self.model == 0:  # beta-model
            M_gas_below_r = 4*numpy.pi*self.rc**3*self.rho0 * (r/self.rc - numpy.arctan(r/self.rc))

        if self.model == 1:  # cut-off beta-model
            rho0 = self.rho0.value_in(units.MSun/units.kpc**3)
            r = r.value_in(units.kpc)
            r2 = p2(r)
            rc = self.rc.value_in(units.kpc)
            rc2 = p2(rc)
            rcut = self.rcut.value_in(units.kpc)
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
            M_gas_below_r = M_gas_below_r | units.MSun

        if self.model == 2:  # cut-off double beta-model
            raise Exception("Error: analytical double beta model not implemented")

        return M_gas_below_r.as_quantity_in(units.MSun)

    def get_r200(self):
        # TODO: implement this function from analytical solution
        cc = CosmologyCalculator()
        rhocrit200 = 200*cc.rho_crit()
        M_gas_below_r = self.gas_mass()
        M_dm_below_r = M_gas_below_r/0.17
        rho_average = (M_dm_below_r / (4./3*numpy.pi*self.radius**3)).value_in(units.g/units.cm**3)
        r200 = self.radius[(numpy.abs(rho_average-rhocrit200)).argmin()]
        print rhocrit200
        print r200

        pyplot.figure()
        amuse_plot.plot(self.radius, self.gas_density())
        amuse_plot.plot(self.radius, rho_average)
        pyplot.gca().set_xscale("log")
        pyplot.gca().set_yscale("log")
        pyplot.axhline(rhocrit200)
        pyplot.show()

    def dm_density(self, r=None):
        """ Dark Matter density radial profile rho(r). Given by Hernquist (1990) """
        if not r:
            r = self.radius

        rho_dm = self.M_dm/(2*numpy.pi) * self.a / (r*p3(r + self.a))
        return rho_dm.as_quantity_in(units.g/units.cm**3)

    def dm_cummulative_mass(self, r=None):
        """ Dark Matter cummulative mass profile M_DM(<r). Given by Hernquist (1990) """
        if not r:
            r = self.radius

        M_dm_below_r = self.M_dm * r**2 / (r + self.a)**2
        return M_dm_below_r.as_quantity_in(units.MSun)

    def average_dm_density(self, r):
        rho_average = (self.dm_cummulative_mass(r) / (4./3*numpy.pi*r**3))
        return rho_average.as_quantity_in(units.g/units.cm**3)

    def gas_temperature(self, r):
        F_1 = numpy.pi**2/(8*self.rc) - numpy.arctan(r/self.rc)**2/(2*self.rc) -\
            numpy.arctan(r/self.rc)/r

        # TODO: use mu from convert.py?
        # T_r_gas = constants.G*self.mu*constants.proton_mass/constants.kB *\
        #    (1 + r**2/self.rc**2) * (4*numpy.pi*self.rc**3*self.rho0gas*F_1)
        return T_r_gas.as_quantity_in(units.K)

    def dm_temperature(self, r):
        F_0 = self.rc/(self.a**2 + self.rc**2)**2 * (numpy.pi/2*(self.a**2 - self.rc**2) +\
            self.r_c*(self.a**2 + self.rc**2)/(self.a + r) - (self.a**2 - self.rc**2)*\
            numpy.arctan(r/self.rc) - self.rc*self.a*numpy.log((self.a + r)**2/(r**2 + self.rc**2)))
        # TODO: use mu from convert.py?
        #T_r_dm = constants.G*self.mu*constants.proton_mass/constants.kB * (1 + r**2/self.rc**2) *\
        #    (self.M_dm*F_0)
        return T_r_dm.as_quantity_in(units.K)

    def temperature(self, r):
        return (self.gas_temperature(r) + self.dm_temperature(r)).as_quantity_in(units.K)

    def characteristic_temperature_analytically(self):
        # Charecteristic temperature
        # TODO: use mu from convert.py?
        #T_c = 2*constants.G*self.mu*constants.proton_mass/constants.kB * \
        #    (self.M_dm * self.rc**2/(self.a**2 + self.rc**2)**2 * (numpy.pi/(4*self.rc) * \
        #    (self.a**2 - self.rc**2) + (self.a**2 + self.rc**2)/(self.a + self.rc) - \
        #    self.a*numpy.log((self.a + self.rc)**2/(2*self.rc**2))) \
        #    + numpy.pi**2*self.rc**2*self.rho_0*(3*numpy.pi/8 - 1))
        return T_c.as_quantity_in(units.K)

    def gas_pressure(self, r):
        # TODO: use mu from convert.py?
        #P_gas = constants.kB*self.temperature(r)*self.gas_density(r)/\
        #    (self.mu*constants.proton_mass)
        return P_gas.as_quantity_in(units.Pa)

    def dm_pressure(self, r):
        # TODO: use mu from convert.py?
        #P_dm = constants.kB*self.temperature(r)*self.dm_density(r)/\
        #    (self.mu*constants.proton_mass)
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
        # TODO: use mu from convert.py?
        #P_500 = 3/(8*numpy.pi) * ((500./2)*constants.G**(-1./4)*self.Hubble_of_z**2)**(4./3) *\
        #    self.mu/self.mu_e * self.f_b * self.M_500**(2./3)

        return P_500.as_quantity_in(units.Pa)


class SampledBox(object):
    """ The Toycluster output contains a box with two haloes living
    inside it. A halo consists of gas and DM particles, but we do
    not know the order of the particles.

    We do know the x-axis is the merger axis, so we split the two
    haloes by taking all particles on the left (negative x), and the
    particles on the right (positive x).

    Here we explicitly split the particles of both haloes, but this
    only works for Toycluster output (i.e. when we explicitly know the
    positions of the centers of mass for the haloes).

    We then give the SampledBox class an AnalyticalCluster instance
    for both haloes, and place the sampled particles for the halo
    in a NumericalCluster instance which we again place at (0, 0, 0).

    NB, we initially do use the NumericalCluster instance because
    it already parses the F77Unformatted data rather neatly and places
    it in an AMUSE datamodel, which we then use take subsets of the
    particles in the datamodel :-)... """

    def __init__(self, timestamp):
        debug = True

        """ First read IC output where two clusters live in the box.
        NB simulation to make distinction between `numerical' (has 1 cluster),
        and `simulation' where two clusters live in the box. """
        simulation = NumericalCluster(
            icdir="../runs/{0}/ICs/".format(timestamp),
            snapdir="../runs/{0}/ICs/".format(timestamp),
            logfile="runToycluster.log",
            icfile="IC_single_0")

        # NB simulation contains particles of both haloes
        gas, dm = simulation.place_ic_data_in_datamodel(simulation.raw_data)

        # Placeholders for now, will be filled with single-cluster data later
        # NB this is rather stupid but we cannot say halo0_numerical = simulation
        # because then id(halo0_numerical) = id(halo1_numerical) = id(simulation)
        halo0_numerical = NumericalCluster(
            icdir="../runs/{0}/ICs/".format(timestamp),
            snapdir="../runs/{0}/ICs/".format(timestamp),
            logfile="runToycluster.log",
            icfile="IC_single_0")

        halo1_numerical = NumericalCluster(
            icdir="../runs/{0}/ICs/".format(timestamp),
            snapdir="../runs/{0}/ICs/".format(timestamp),
            logfile="runToycluster.log",
            icfile="IC_single_0")

        """ Second, set up analytical models of both sampled clusters. """
        halo0_analytical = setup_analytical_cluster(simulation, i=0)
        halo1_analytical = setup_analytical_cluster(simulation, i=1)

        # Histogram of gas x-values is beautifully bimodal (y, z aint)
        if debug:
            pyplot.figure(figsize=(12, 12))
            amuse_plot.hist(gas.x, bins=int(numpy.sqrt(len(gas.x))))
            pyplot.show()


        """ Third, split up particles in two haloes. Shift back to (0,0,0) """
        print "Splitting up haloes, halo0: x<0; halo1: x>0."
        # The merger-axis is x. We know halo center x position, is D_CoM_0/1
        # TODO: can we assume the box is split in two at x=0?
        # halo0 lives on the left-hand side of the box (negative x)
        halo0gas = gas.select_array(lambda l : l < 0.0 | units.kpc, ["x"])
        halo0dm = dm.select_array(lambda l : l < 0.0 | units.kpc, ["x"])
        # halo1 lives on the left-hand side of the box (negative x)
        halo1gas = gas.select_array(lambda l : l > 0.0 | units.kpc, ["x"])
        halo1dm = dm.select_array(lambda l : l > 0.0 | units.kpc, ["x"])

        # TODO: check if this routine is valid only for -DCOMET ?
        # TODO: boxhalf is added in Shift_Origin, but then subtracted in
        # Apply_kinematics ?
        # NB, already corrected for in place_ic_data_in_datamodel
        boxhalf = simulation.toyclusterlog.systemsetup['Boxsize']/2

        # The x-position is shifted back from the CoM to the center (x=0)
        # TODO: properly shift back
        # d_clusters = 0.9 * (halo[0].R200 + halo[1].R200)
        # D_CoM_0 = -1 * d_clusters * halo[1].Mtotal200/param.Mtot200 (why?)
        # D_CoM_1 = d_clusters + D_CoM_0
        # dx =
        print "Shifting haloes, back to origin."
        halo0gas.x -= simulation.toyclusterlog.kinematics['D_CoM_0']
        halo0dm.x -= simulation.toyclusterlog.kinematics['D_CoM_0']
        halo1gas.x -= simulation.toyclusterlog.kinematics['D_CoM_1']
        halo1dm.x -= simulation.toyclusterlog.kinematics['D_CoM_1']

        # TODO: when an impact parameter is given, then y is also shifted!
        if not (-2**-14 < simulation.toyclusterlog.kinematics['Impact_Parameter'].value_in(units.kpc) < 2**-14):
            # TODO: properly shift back
            # b_CoM_0 = -1 * param.Impact_Param * halo[0].Mtotal200/param.Mtot200 (why?)
            # b_CoM_1 = param.Impact_Param + b_CoM_0
            #
            print "WARNING: correcting for impact parameter is untested"
            halo0gas.y -= simulation.toyclusterlog.kinematics['b_CoM_0']
            halo0dm.y -= simulation.toyclusterlog.kinematics['b_CoM_0']
            halo1gas.y -= simulation.toyclusterlog.kinematics['b_CoM_1']
            halo1dm.y -= simulation.toyclusterlog.kinematics['b_CoM_1']

        if debug:
            # Bimodality is gone; cluster now centered around x=0
            pyplot.figure(figsize=(12, 12))
            amuse_plot.hist(halo0gas.x, bins=int(numpy.sqrt(len(halo0gas.x))))
            pyplot.show()

            pyplot.figure(figsize=(12, 12))
            amuse_plot.hist(halo1gas.x, bins=int(numpy.sqrt(len(halo1gas.x))))
            pyplot.show()

        import sys; sys.exit(0)

        # recalculate r because now it is calculated with uncentered x, y values
        print "Recalculating halo radii."
        halo0gas.r = numpy.sqrt(p2(halo0gas.x.value_in(units.kpc))+
                                p2(halo0gas.y.value_in(units.kpc))+
                                p2(halo0gas.z.value_in(units.kpc))) | units.kpc
        halo0dm.r = numpy.sqrt(p2(halo0dm.x.value_in(units.kpc))+
                               p2(halo0dm.y.value_in(units.kpc))+
                               p2(halo0dm.z.value_in(units.kpc))) | units.kpc
        halo1gas.r = numpy.sqrt(p2(halo1gas.x.value_in(units.kpc))+
                                p2(halo1gas.y.value_in(units.kpc))+
                                p2(halo1gas.z.value_in(units.kpc))) | units.kpc
        halo1dm.r = numpy.sqrt(p2(halo1dm.x.value_in(units.kpc))+
                               p2(halo1dm.y.value_in(units.kpc))+
                               p2(halo1dm.z.value_in(units.kpc))) | units.kpc

        # Now fill the placeholders created earlier with the split up, reshifted haloes
        print "Filling AMUSE datamodel with particle properties."
        halo0_numerical.set_toycluster2_values(i=0)
        halo0_numerical.gas = halo0gas
        halo0_numerical.dm = halo0dm
        halo0_numerical.M_dm = halo0_numerical.Mass_in_DM
        halo0_numerical.M_gas = halo0_numerical.Mass_in_gas
        # Rather ugly hack, but ensures set_dm_density() method works
        halo0_numerical.raw_data.Ndm = len(halo0_numerical.dm)

        halo1_numerical.set_toycluster2_values(i=1)
        halo1_numerical.gas = halo1gas
        halo1_numerical.dm = halo1dm
        halo1_numerical.M_dm = halo1_numerical.Mass_in_DM
        halo1_numerical.M_gas = halo1_numerical.Mass_in_gas
        # Rather ugly hack, but ensures set_dm_density() method works
        halo1_numerical.raw_data.Ndm = len(halo1_numerical.dm)

        self.halo0_analytical = halo0_analytical
        self.halo0_numerical = halo0_numerical
        self.halo1_analytical = halo1_analytical
        self.halo1_numerical = halo1_numerical

        # To plot density profiles :-)...
        # self.halo0_numerical.get_gas_mass_via_density()
        # self.halo0_numerical.get_dm_mass_via_number_density()
        # self.halo0_numerical.set_dm_density()

        # self.halo1_numerical.get_gas_mass_via_density()
        # self.halo1_numerical.get_dm_mass_via_number_density()
        # self.halo1_numerical.set_dm_density()

        return

        # TODO: this can be used in plot_individual_cluster_density to sort
        # Sort subset of the particles by the radius such that we can
        # integrate in spherical shells to obtain the radial density profile
        sorted0gas = halo0gas.sorted_by_attribute('r')
        sorted1gas = halo1gas.sorted_by_attribute('r')
        sorted0dm = halo0dm.sorted_by_attribute('r')
        sorted1dm = halo1dm.sorted_by_attribute('r')


def setup_analytical_cluster(simulation, i=0):
    """ Set up an analytical cluster given sampled parameters of the halo.
        In Toycluster we have two haloes: halo[0] and halo[1].

        Toycluster runtime output is piped to runToycluster.log.
        This output is then parsed by `Toycluster2RuntimeOutputParser' class
        and saved in the `NumericalCluster' instance `simulation' class
        variable `toyclusterlog'.

        toyclusterlog has a 2D dict halosetup, where
        halosetup[0][...] gives the specific halo properties of halo0
    """

    # Caution: parms[0] is number density! Caution: use unitsless numbers!
    rho0 = simulation.toyclusterlog.halosetup[i]['rho0gas_cgs']
    ne0 = convert.rho_to_ne(rho0.value_in(units.g/units.cm**3),
                            simulation.toyclusterlog.systemat['z'])
    rc = simulation.toyclusterlog.halosetup[i]['rc']
    R200 = simulation.toyclusterlog.halosetup[i]['R200']
    Mass_in_DM = simulation.toyclusterlog.halosetup[i]['Mass_in_DM']
    a = simulation.toyclusterlog.halosetup[i]['a_hernquist']

    # print "rho0          :", rho0
    # print "ne0           :", ne0
    # print "rc            :", rc
    # print "R200          :", R200
    # print "Mass_in_DM    :", Mass_in_DM
    # print "a             :", a

    parms = (ne0, rc.value_in(units.kpc),
             R200.value_in(units.kpc))
    dm_parms = (Mass_in_DM, a)

    return AnalyticalCluster(parms, dm_parms,
        z=simulation.toyclusterlog.systemat['z'])


def simple_plot(observed_cluster, yparm="density"):
    """ Possible properties to plot on y-axis are:
            bin_volume, density*, pressure*, compton_y, source_sb*, background_sb*
        * means std is available --> errorbar plot
        x-axis is computer from inner_radius and/or outer_radius
    """

    pyplot.figure(figsize=(12, 9))
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    if yparm not in ["bin_volume", "compton_y"]:
        pyplot.errorbar(observed_cluster.radius, getattr(observed_cluster, yparm),
            yerr=getattr(observed_cluster, yparm+"_std"))
    else:
        pyplot.plot(observed_cluster.radius, getattr(observed_cluster, yparm))

    pyplot.xlabel("Radius")
    pyplot.ylabel(yparm)
    # pyplot.legend()


def compare_old_and_new_data(old, new, yparm="density"):
    pyplot.figure(figsize=(12, 9))
    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")

    if yparm not in ["bin_volume", "compton_y"]:
        pyplot.errorbar(old.radius, getattr(old, yparm),
            yerr=getattr(old, yparm+"_std"), label="old")
        pyplot.errorbar(new.radius, getattr(new, yparm),
            yerr=getattr(new, yparm+"_std"), label="new")
    else:
        pyplot.plot(old.radius, getattr(old, yparm), label="old")
        pyplot.plot(new.radius, getattr(new, yparm), label="new")

    pyplot.xlabel("Radius")
    pyplot.ylabel(yparm)
    pyplot.legend()


if __name__ == "__main__":
    print "Reading Observed Cluster"
    print 80*'-'
    cygA_observed_800ksecOLD = ObservedCluster("cygA", oldICs=True)
    cygB_observed_800ksecOLD = ObservedCluster("cygB", oldICs=True)

    cygA_observed_800ksec = ObservedCluster("cygA")
    cygB_observed_800ksec = ObservedCluster("cygB")

    compare_old_and_new_data(cygA_observed_800ksecOLD, cygA_observed_800ksec)
    compare_old_and_new_data(cygB_observed_800ksecOLD, cygB_observed_800ksec)

    pyplot.show()
    # TODO: Plot difference between 800 and 900 ksec observation
    import sys; sys.exit(0)

    debug = False
    if debug:
        print "Debug information after parsing Martijn's Chandra observation data"
        print 80*"-"
        print cygA_observed
        print cygB_observed
        print

        for property in ["bin_volume", "density", "number_density", "pressure",
                "compton_y", "source_sb", "background_sb"]:
            print "Plotting property:", property
            simple_plot(cygA_observed, property)
            simple_plot(cygB_observed, property)
        print 80*"-"

    print "Reading Toycluster Run without WVT relax"
    print 80*'-'
    numerical_cluster = NumericalCluster(
        icdir="../runs/20160623T1755/ICs/",
        snapdir="../runs/20160623T1755/ICs/",
        logfile="runToycluster.log",
        icfile="IC_single_0")
    # numerical_cluster.perform_sanity_checks()
    print numerical_cluster.gas.u
    print 80*'-'

    pyplot.show()
