""" Number density to mass density; gadget <--> cgs density
    Based on Julius' IDL script gadget_density

    NB we assume lambda CDM concordance cosmology /w h=0.7
"""

from amuse.units import units
from amuse.units import constants

uMass = 1.989e43       # in gram (1e10 MSun)
uLength = 3.08568e21   # cm (kpc)

# Hydrogen fraction
xH = 0.76              # Gadget-2 uses 0.76 in allvars.h
# Mean molecular weight
umu = 4./(5.*xH+3.)    # Gadget-2 read_ic.c:129, assuming full ionisation (T>1e4)

h = 0.7   # Assuming Lambda CDM, concordance, H0(z=0) = 70 km/s/Mpc

# Ununed here, but present in Julius' IDL script. What does this do?
# total number density to electron number density only?
# conversion n_pat -> n_electrons
# n2ne = (xH+0.5*(1-xH))/(2*xH+0.75*(1-xH))

def gadget_units_to_cgs(rho):
    """ convert mass density in gadget units to cgs units """
    return uMass/uLength**3 * rho * h**2

def cgs_to_gadget_units(rho):
    """ convert mass density in cgs units to gadget units """
    return 1./(uMass/uLength**3) * rho / h**2

def rho_to_ne(rho, z=None):
    """ convert mass density to electron number density """

    # TODO: if z: comoving?
    # Julius uses (+z)**3*h**2 too, but this is for comoving density?

    ne = rho/(umu*constants.proton_mass.value_in(units.g))
    return ne

def ne_to_rho(ne, z):
    """ convert electron number density to mass density """

    # TODO: if z: comoving?

    rho = umu*constants.proton_mass.value_in(units.g)*ne
    return rho

if __name__ == "__main__":
    print '-'*80
    print "Example Toycluster output"
    print "rho0_gas            = 1.97444e-26 g/cm^3"
    print "rho0_gas            = 5.95204e-05 [gadget]"
    example_mass_density = 1.97444e-26
    mass_density_to_gadget_density = cgs_to_gadget_units(example_mass_density)
    print
    print "Given mass density  : {0:1.5e}".format(example_mass_density)
    print "Gadget density      : {0:1.5e}".format(mass_density_to_gadget_density)
    example_gadget_density = 5.95204e-05
    gadget_density_to_mass_density = gadget_units_to_cgs(example_gadget_density)
    print
    print "Given gadget density: {0:1.5e}".format(example_gadget_density)
    print "Mass density        : {0:1.5e}".format(gadget_density_to_mass_density)
    print '-'*80

    print  "Converting CygA/B electron number density to mass density"
    print '-'*80
    # best fit value for CygA n_e0
    cygA_electron_number_density = 1.35e-1  # cm**-3, from fit to Chandra data
    print "Number density CygA = {0:1.3e} 1/cm**3".format(cygA_electron_number_density)
    print "Mass density CygA   = {0:1.3e} g/cm**3".format(
        ne_to_rho(cygA_electron_number_density, 0.0562))
    print
    # best fit value for CygB n_e0
    cygB_electron_number_density = 1.94e-3  # cm**-3, from fit to Chandra data
    print "Number density CygB = {0:1.3e} 1/cm**3".format(cygB_electron_number_density)
    print "Mass density CygB   = {0:1.3e} g/cm**3".format(
        ne_to_rho(cygB_electron_number_density, z=0.070))
    print '-'*80
