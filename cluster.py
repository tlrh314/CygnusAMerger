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

            Data obtained by M. N. de Vries from 900 ksec Chandra data

            @param oldICs: True -> old (800 ksec) data, else latest (900 ksec)
        """

        self.name = name
        self.oldICs = oldICs

        if self.oldICs:
            density_file = "data/{0}_1T_fixnH_pressureprofile_800ksec.dat".format(name)
            radius_file = "data/{0}_sn100_sbprofile_800ksec.dat".format(name)
        else:
            density_file = "data/{0}_1T_pressureprofile_900ksec.dat".format(name)
            radius_file = "data/{0}_sb_sn100_900ksec.dat".format(name)

        if self.name == "cygA":
            z = 0.0562
        if self.name == "cygB":
            z = 0.0562
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

        raw_sn100 = pandas.read_csv(radius_file,
            delimiter="|" if self.oldICs else None,
            delim_whitespace=None if self.oldICs else True,
            header=0 if self.oldICs else 18)

        if self.oldICs:
            self.inner_radius = \
                raw_sn100[" Inner radius (arcsec) "].as_matrix()
            self.outer_radius = \
                raw_sn100[" Outer radius (arcsec) "].as_matrix()

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

        else:
            self.inner_radius = raw_sn100["Radius1"].as_matrix()
            self.outer_radius = raw_sn100["Radius2"].as_matrix()
            self.source_sb = raw_sn100["SB"].as_matrix()
            self.source_sb_std = raw_sn100["SBError"].as_matrix()
            self.background_sb = raw_sn100["BGRD"].as_matrix()
            self.background_sb_std = raw_sn100["BGError"].as_matrix()
            self.bin_number_sn100 = range(len(raw_sn100))

        if not numpy.array_equal(self.bin_number, self.bin_number_sn100):
# Yeah, I know way too generic to raise, but w/e dont wanna define ReadingTwoDataFilesIsProneToMakeMistakesExeption *_*
            raise Exception("Mismatch in the bin numbers of the two data input files")


class NumericalCluster(object):
    def __init__(self, icdir, snapdir, logfile=None, icfile=None, verbose=True):
        # Output of runToycluster writes files with these filename
        if logfile is None:
            logfile="runToycluster.log"
        if icfile is None:
            icfile="IC_single_0"
        self.verbose = verbose

        # Read runtime output of Toycluster 2.0. Parse runtime output
        self.toyclusterlog = Toycluster2RuntimeOutputParser(filename=icdir+logfile)
        self.raw_data = Gadget2BinaryF77UnformattedType2Parser(snapdir+icfile, verbose=self.verbose)

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

    def get_dm_mass_via_number_density(self, log_binning=True):
        """ Count particles <r (= number density). Obtain DM mass from it """

        # TODO: find out if this method can also be used for sph mass
        # should we take into account the kernel, and if so how?

        if self.verbose:
            print "Counting particles for which radii < r to obtain M(<r)"

        if log_binning:
            radii = numpy.power(10, numpy.linspace(numpy.log(1), numpy.log(1e5), 1001))
            dr = radii[1:] - radii[:-1]
            radii = radii[:-1]
        else:
            radii = numpy.arange(0, 1e5, 10)
            dr = radii[1] - radii[0]  # Because evenly spaced
        N = len(radii)

        particles = numpy.zeros(N)
        for i, r in enumerate(radii):
            particles[i] = ((numpy.where(self.dm.r.value_in(units.kpc) < r)[0]).size)
            if self.verbose and (i==(N-1) or i%100 == 0):
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

    def set_gas_temperature(self):
        """ Thermodynamics in Gadget-2:
            - Ideal gas, close to pure hydrogen.
            - SPH is energy conserving, expansion is always adiabatic
            - Eint = (gamma-1)^-1 N kB T; where gamma = 5/3
            - PV = N kB T ==> P = rho kB/(mu m_p) T

            - U = (gamma-1)^-1 kB/(mu m_p) T [erg/g]

            Here we set the gas temperature from the internal energy
            """

        # T = (gamma-1)*(mu*m_p)/kB U
        gamma = 5.0/3
        kB = constants.kB.value_in(units.erg/units.K)
        m_p = constants.proton_mass.value_in(units.g)
        factor = (gamma-1)*convert.umu*m_p/kB
        # NB internal energy u is in code units energy/mass
        # unit.mass = 1e10 MSun = 1e10 * 1.9889e33 gram
        # unit.velocity = km/sec = 1e5 cm/sec
        # unit.energy = unit.mass*unit.velocity**2 = 1.9889e53 erg
        # unit.energy/unit.mass = 1.9889e53/1.9889e43 = 1e10
        # So code unit energy per unit mass yields a factor 1e10
        self.gas.T = factor * self.gas.u * 1e10 | units.K

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
    def __init__(self, parms, dm_parms=None, radius=None, z=0, free_beta=False):
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
        self.free_beta = free_beta
        if len(parms) == 3 and not free_beta:
            self.model = 1
            self.rcut = parms[2] | units.kpc
        if len(parms) == 3 and free_beta:
            self.model = 3
            self.beta = parms[2]
        if len(parms) == 5:
            self.model = 2
            ne0_fac = parms[3]
            rc_fac = parms[4]
            self.ne0_cc = ne0 * ne0_fac
            self.rho0_cc = self.ne0_cc * globals.mu * (globals.m_p | units.g)
            self.rc_cc = rc / rc_fac

        modelnames = {0: r"$\beta=2/3$",
                      1: r"cut-off $\beta=2/3$",
                      2: r"cut-off double $\beta (2/3)$",
                      3: r"free $\beta$"}
        self.modelname = modelnames[self.model]

        if dm_parms:
            self.M_dm = dm_parms[0]
            self.a = dm_parms[1]

        # self.get_r200()

    def gas_density(self, r=None):
        """ For details of beta-model see Cavaliere & Fusco-Femiano (1978).

            model    name                   reference
            0        beta-model             Donnert (2014)
            1        cut-off beta           Donnert et al. (2016, in prep)
            2        cut-off double beta    Donnert et al. (2016, in prep)
            3        beta-model             Donnert (2014), beta free param
        """

        if r is None:
            r = self.radius

        beta = 2./3  # Mastropietro & Burkert (2008)
        if self.free_beta:
            beta = self.beta
        rho_gas = self.rho0 * (1 + p2(r/self.rc))**(-3*beta/2.)  # model >= 0
        if self.model >= 1 and not self.free_beta:
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

    def F0(self, r=None):
        """ Helper function for temperature calculation -- dm part
                Donnert (2014) equation 11 """
        if not r:
            r = self.radius

        rc2 = p2(self.rc)
        a2 = p2(self.a)
        a = self.a
        rc = self.rc

        result = rc / p2(a2 + rc2)
        result *= (numpy.pi/2*(a2 - rc2) + rc*(a2 + rc2)/(a + r) \
            - (a2 - rc2)*numpy.arctan(r/rc) -\
            rc*a*numpy.log(p2(a + r)/(p2(r) + rc2)))

        return result

    def F1(self, r=None):
        """ Helper function for temperature calculation -- gas part
                Donnert (2014) equation 12 """
        if not r:
            r = self.radius

        rc = self.rc

        return p2(numpy.pi)/(8*rc) - p2(numpy.arctan(r/rc))/(2*rc) -\
            numpy.arctan(r/rc)/r

    def dm_temperature(self, r=None):
        """ Donnert (2014) equation 10 -- dm part: ...*[ M_DM F0(r) ]
        TODO: what is DM temperature? """
        if not r:
            r = self.radius

        rc = self.rc
        m_p = constants.proton_mass
        kB = constants.kB

        T_r_dm = constants.G*convert.umu*m_p/kB *\
            (1 + p2(r/rc)) * (self.M_dm*self.F0(r))
        return T_r_dm.as_quantity_in(units.K)

    def gas_temperature(self, r=None):
        """ Donnert (2014) equation 10 -- gas part: ...*[ ...*F1(r) ] """
        if not r:
            r = self.radius

        rc = self.rc
        m_p = constants.proton_mass
        kB = constants.kB

        T_r_gas = constants.G*convert.umu*m_p/kB *\
            (1 + p2(r/rc)) * (4*numpy.pi*p3(rc)*self.rho0*self.F1(r))
        return T_r_gas.as_quantity_in(units.K)

    def temperature(self, r=None):
        """ Donnert (2014) equation 10 -- gas + dm (full equation) """
        if not r:
            r = self.radius

        return (self.gas_temperature(r) + \
                self.dm_temperature(r)).as_quantity_in(units.K)

    def characteristic_temperature(self):
        """ Donnert (2014) equation 13 -- Tc = T(rc) """

        rc2 = p2(self.rc)
        a2 = p2(self.a)
        a = self.a
        rc = self.rc
        m_p = constants.proton_mass
        kB = constants.kB

        T_c = 2*constants.G*convert.umu*m_p/kB * \
            (self.M_dm * rc2/p2(a2 + rc2) * (numpy.pi/(4*rc) * \
            (a2 - rc2) + (a2 + rc2)/(a + rc) - \
            a*numpy.log(p2(a + rc)/(2*rc2))) \
            + p2(numpy.pi)*rc2*self.rho0*(3*numpy.pi/8 - 1))
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

        icdir = "../runs/{0}/ICs/".format(timestamp)
        simulation = NumericalCluster(
            icdir=icdir,
            snapdir=icdir,
            logfile="runToycluster.log",
            icfile="IC_single_0")

        boxsize = simulation.toyclusterlog.systemsetup['Boxsize']\
            .value_in(units.kpc)

        # import os
        # tc_par_name = [file for file in os.listdir(icdir) if "par" in file][0]
        # NB simulation contains particles of both haloes
        gas, dm = simulation.place_ic_data_in_datamodel(simulation.raw_data)

        # Placeholders for now, will be filled with single-cluster data later
        # NB this is rather stupid but we cannot say halo0_numerical = simulation
        # because then id(halo0_numerical) = id(halo1_numerical) = id(simulation)
        halo0_numerical = NumericalCluster(
            icdir=icdir,
            snapdir=icdir,
            logfile="runToycluster.log",
            icfile="IC_single_0")

        halo1_numerical = NumericalCluster(
            icdir=icdir,
            snapdir=icdir,
            logfile="runToycluster.log",
            icfile="IC_single_0")

        """ Second, set up analytical models of both sampled clusters. """
        halo0_analytical = setup_analytical_cluster(simulation, i=0)
        halo1_analytical = setup_analytical_cluster(simulation, i=1)

        hist, edges = numpy.histogram(gas.x.value_in(units.kpc),
            bins=int(numpy.sqrt(len(gas.x))))
        if debug:
            pyplot.figure()
            pyplot.plot((edges[1:]+edges[:-1])/2, hist)
            pyplot.show()
        # Domain contains indices of x values between -750 and 750
        # somewhere in this range there is a minimum x-value, which
        # is the center that we need to shift back the haloes.
        domain = numpy.where(numpy.logical_and(edges>=-1000, edges<=1000))
        ymin_index = numpy.argmin(hist[domain])
        center = edges[domain][ymin_index]
        if debug:
            # numpy.amin(hist[domain])
            pyplot.figure()
            pyplot.plot(edges[domain], hist[domain])
            ymin = hist[domain][ymin_index]
            pyplot.axhline(ymin)
            pyplot.axvline(center)
            pyplot.show()

        # Histogram of gas x-values is beautifully bimodal (y, z aint)
        if debug:
            pyplot.figure(figsize=(12, 12))
            amuse_plot.hist(gas.x, bins=int(numpy.sqrt(len(gas.x))))
            pyplot.show()


        """ Third, split up particles in two haloes. Shift back to (0,0,0) """
        print "Splitting up haloes, halo0: x<{0}; halo1: x>{0}.".format(center)
        # The merger-axis is x. We know halo center x position, is D_CoM_0/1
        # TODO: can we assume the box is split in two at x=0?
        # halo0 lives on the left-hand side of the box (negative x)
        halo0gas = gas.select_array(lambda l : l < center | units.kpc, ["x"])
        halo0dm = dm.select_array(lambda l : l < center | units.kpc, ["x"])
        # halo1 lives on the left-hand side of the box (negative x)
        halo1gas = gas.select_array(lambda l : l > center | units.kpc, ["x"])
        halo1dm = dm.select_array(lambda l : l > center | units.kpc, ["x"])

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
    pyplot.title(new.name)
    pyplot.legend()


if __name__ == "__main__":
    print "Reading Observed Cluster"
    print 80*'-'

    cygA_observed_800ksec = ObservedCluster("cygA", oldICs=True)
    cygA_observed_900ksec = ObservedCluster("cygA")

    cygB_observed_800ksec = ObservedCluster("cygB", oldICs=True)
    cygB_observed_900ksec = ObservedCluster("cygB")
    print "Done reading ObservedCluster"
    print 80*'-'

    debug = False
    if debug:
        print "Printing attributes of old- and new ICs"
        print 80*'-'
        for attr in dir(cygA_observed_800ksec):
            if not attr.startswith("__"):
                print "Comparing", attr
                old = getattr(cygA_observed_800ksec, attr)
                new = getattr(cygA_observed_900ksec, attr)
                if type(old) == numpy.ndarray:
                    print "old head\n", old[0:10]
                    print "new head\n", new[0:10]
                    print "old tail\n", old[-9:]
                    print "new tail\n", new[-9:]
                else:
                    print "old\n", old, "\nnew\n", new
                # raw_input("Press enter to continue\n")
                print "\n"
        print 80*'-'

    compare_clusters = True
    if compare_clusters:
        # print cygA_observed_800ksec
        # print cygA_observed_900ksec
        print "Generating plots for comparison of attributes"
        print 80*'-'
        for property in ["bin_volume", "density", "number_density", "pressure",
                "compton_y", "source_sb", "background_sb"]:
            compare_old_and_new_data(cygA_observed_800ksec,
                cygA_observed_900ksec, yparm=property)
            compare_old_and_new_data(cygB_observed_800ksec,
                cygB_observed_900ksec, yparm=property)
        print "Done generating plots for comparison of attributes"
        print 80*'-'

    if debug:
        print "Debug information after parsing Martijn's Chandra observation data"
        print 80*"-"
        print cygA_observed_900ksec
        print cygB_observed_900ksec
        print

        for property in ["bin_volume", "density", "number_density", "pressure",
                "compton_y", "source_sb", "background_sb"]:
            print "Plotting property:", property
            simple_plot(cygA_observed_900ksec, property)
            simple_plot(cygB_observed_900ksec, property)
        print 80*"-"

    numerical_sanity_check = False
    if numerical_sanity_check:
        print "Reading Toycluster Run without WVT relax"
        print 80*'-'
        numerical_cluster = NumericalCluster(
            icdir="../runs/20160623T1755/ICs/",
            snapdir="../runs/20160623T1755/ICs/",
            logfile="runToycluster.log",
            icfile="IC_single_0")
        numerical_cluster.perform_sanity_checks()
        print numerical_cluster.gas.u
        print 80*'-'

    pyplot.show()
