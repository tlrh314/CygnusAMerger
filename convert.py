""" Number density to mass density; gadget <--> cgs density
    Based on Julius' IDL script gadget_density

    NB we assume lambda CDM concordance cosmology /w h=0.7
"""

import numpy

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

def ne_to_rho(ne, z=None):
    """ convert electron number density to mass density """

    # TODO: if z: comoving?

    rho = umu*constants.proton_mass.value_in(units.g)*ne
    return rho


# --------------------------------------------------------------------------- #
# Stolen from https://github.com/aplpy/aplpy/blob/master/aplpy/wcs_util.py
def precession_matrix(equinox1, equinox2, fk4=False):
    "Adapted from the IDL astronomy library"

    deg_to_rad = numpy.pi / 180.
    sec_to_rad = deg_to_rad / 3600.

    t = 0.001 * (equinox2 - equinox1)

    if not fk4:

        st = 0.001 * (equinox1 - 2000.)

        # Compute 3 rotation angles
        a = sec_to_rad * t * (23062.181 + st * (139.656 + 0.0139 * st) + t * (30.188 - 0.344 * st + 17.998 * t))
        b = sec_to_rad * t * t * (79.280 + 0.410 * st + 0.205 * t) + a
        c = sec_to_rad * t * (20043.109 - st * (85.33 + 0.217 * st) + t * (- 42.665 - 0.217 * st - 41.833 * t))

    else:

        st = 0.001 * (equinox1 - 1900.)

        # Compute 3 rotation angles
        a = sec_to_rad * t * (23042.53 + st * (139.75 + 0.06 * st) + t * (30.23 - 0.27 * st + 18.0 * t))
        b = sec_to_rad * t * t * (79.27 + 0.66 * st + 0.32 * t) + a
        c = sec_to_rad * t * (20046.85 - st * (85.33 + 0.37 * st) + t * (- 42.67 - 0.37 * st - 41.8 * t))

    sina = numpy.sin(a)
    sinb = numpy.sin(b)
    sinc = numpy.sin(c)
    cosa = numpy.cos(a)
    cosb = numpy.cos(b)
    cosc = numpy.cos(c)

    r = numpy.matrix([[cosa * cosb * cosc - sina * sinb, sina * cosb + cosa * sinb * cosc,  cosa * sinc],
                   [- cosa * sinb - sina * cosb * cosc, cosa * cosb - sina * sinb * cosc, - sina * sinc],
                   [- cosb * sinc, - sinb * sinc, cosc]])

    return r

P1 = precession_matrix(1950., 2000.)
P2 = precession_matrix(2000., 1950.)


def b1950toj2000(ra, dec):
    '''
    Convert B1950 to J2000 coordinates.
    This routine is based on the technique described at
    http://www.stargazing.net/kepler/b1950.html
    '''

    # Convert to radians
    ra = numpy.radians(ra)
    dec = numpy.radians(dec)

    # Convert RA, Dec to rectangular coordinates
    x = numpy.cos(ra) * numpy.cos(dec)
    y = numpy.sin(ra) * numpy.cos(dec)
    z = numpy.sin(dec)

    # Apply the precession matrix
    x2 = P1[0, 0] * x + P1[1, 0] * y + P1[2, 0] * z
    y2 = P1[0, 1] * x + P1[1, 1] * y + P1[2, 1] * z
    z2 = P1[0, 2] * x + P1[1, 2] * y + P1[2, 2] * z

    # Convert the new rectangular coordinates back to RA, Dec
    ra = numpy.arctan2(y2, x2)
    dec = numpy.arcsin(z2)

    # Convert to degrees
    ra = numpy.degrees(ra)
    dec = numpy.degrees(dec)

    # Make sure ra is between 0. and 360.
    ra = numpy.mod(ra, 360.)
    dec = numpy.mod(dec + 90., 180.) - 90.

    return ra, dec


def j2000tob1950(ra, dec):
    '''
    Convert J2000 to B1950 coordinates.
    This routine was derived by taking the inverse of the b1950toj2000 routine
    '''

    # Convert to radians
    ra = numpy.radians(ra)
    dec = numpy.radians(dec)

    # Convert RA, Dec to rectangular coordinates
    x = numpy.cos(ra) * numpy.cos(dec)
    y = numpy.sin(ra) * numpy.cos(dec)
    z = numpy.sin(dec)

    # Apply the precession matrix
    x2 = P2[0, 0] * x + P2[1, 0] * y + P2[2, 0] * z
    y2 = P2[0, 1] * x + P2[1, 1] * y + P2[2, 1] * z
    z2 = P2[0, 2] * x + P2[1, 2] * y + P2[2, 2] * z

    # Convert the new rectangular coordinates back to RA, Dec
    ra = numpy.arctan2(y2, x2)
    dec = numpy.arcsin(z2)

    # Convert to degrees
    ra = numpy.degrees(ra)
    dec = numpy.degrees(dec)

    # Make sure ra is between 0. and 360.
    ra = numpy.mod(ra, 360.)
    dec = numpy.mod(dec + 90., 180.) - 90.

    return ra, dec
# --------------------------------------------------------------------------- #


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
