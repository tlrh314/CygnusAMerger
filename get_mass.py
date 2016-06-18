import numpy

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


def obtain_M200_bisection(rc, rho0, verbose=False):
    """ We follow the Toycluster (Donnert 2014) setup.c method in reverse.
    If we assume a value for r200 and we assume the baryon fraction at r200
    is equal to 0.17 we are able to obtain a total mass M200 for which
    rho_average(r200) == 200 rho_crit."""

    # Find r200 such that rho200 / rho_crit == 200
    lower = 10 * kpc2cm
    upper = 2000 * kpc2cm

    # bisection method
    epsilon = 0.01
    while upper/lower > 1+epsilon:
        r200 = (lower+upper)/2.

        bf = 0.17  # Mission critical assumption that bf = 0.17 at r200!

        Mgas200 = M_gas_below_r(r200, rho0, rc)
        Mdm200 = Mgas200 * (1/bf - 1)  # TODO: check this step

        M200 = Mgas200 + Mdm200

        cNFW = concentration_parameter(M200)
        rs = r200 / cNFW
        a = hernquist_a(rs, cNFW)

        Mdm = Mdm200 / M_dm_below_r(r200, 1, a)

        # Now r200/rhocrit should of course equal 200. If not? Try different r200
        rho200_over_rhocrit = ( M200 / (4./3 * numpy.pi * p3(r200))) / rho_crit()
        if verbose:
            print "Lower                  = {0:3.1f}".format(lower * cm2kpc)
            print "r200                   = {0:3.1f}".format(r200 * cm2kpc)
            print "Upper                  = {0:3.1f}".format(upper * cm2kpc)
            print "rho_avg(r200)/rho_crit = {0:.1f}".format(rho200_over_rhocrit)
            print "Ratio                  = {0:.1f}".format(rho200_over_rhocrit/200)
            print

        if rho200_over_rhocrit < 200:
            upper = r200
        if rho200_over_rhocrit > 200:
            lower = r200

    print "r200                   = {0:3.1f}".format(r200 * cm2kpc)
    print "rho_avg(r200)/rho_crit = {0:.1f}".format(rho200_over_rhocrit)
    print "bf200                  = {0:1.4f}".format(Mgas200/(Mdm200+Mgas200))
    print "rho0                   = {0:1.4e}".format(rho0)
    print "rc                     = {0:.3f}".format(rc * cm2kpc)
    print "Mgas200                = {0:1.4e}".format(Mgas200 * g2msun)
    print "Mdm200                 = {0:1.4e}".format(Mdm200 * g2msun)
    print "M200                   = {0:1.4e}".format(M200 * g2msun)
    print "cNFW                   = {0:1.4f}".format(cNFW)
    print "rs                     = {0:3.1f}".format(rs * cm2kpc)
    print "a                      = {0:3.1f}".format(a * cm2kpc)
    print "Mdm                    = {0:1.4e}".format(Mdm * g2msun)
    print

    return M200


if __name__ == "__main__":
    print "CygnusA"
    # obtained from beta (2/3) model to Chandra data
    cygA_rc = 27.645 * kpc2cm
    cygA_ne0 = 0.136  # 1/cm**3
    cygA_rho0 = ne_to_rho(cygA_ne0)

    cygA_M200 = obtain_M200_bisection(cygA_rc, cygA_rho0)
    print

    print "CygnusB"
    # obtained from beta (2/3) model to Chandra data
    cygB_rc = 290.909 * kpc2cm
    cygB_ne0 = 1.9397e-03  # 1/cm**3
    cygB_rho0 = ne_to_rho(cygB_ne0)

    cygB_M200 = obtain_M200_bisection(cygB_rc, cygB_rho0)
    print

    print "cygA M200              = {0:1.4e} MSun".format(cygA_M200 * g2msun)
    print "cygB_M200              = {0:1.4e} MSun".format(cygB_M200 * g2msun)
    print
    print "M_CygA/M_CygB          = {0:1.4f}".format(cygA_M200/cygB_M200)
    print "M_CygB/M_CygA          = {0:1.4f}".format(cygB_M200/cygA_M200)
    print "Mtotal                 = {0:1.4e} MSun".format((cygA_M200 + cygB_M200) * g2msun)
