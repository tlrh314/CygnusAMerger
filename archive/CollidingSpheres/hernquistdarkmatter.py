"""
Hernquist dark matter halo generator for clusters of galaxies

This module contains a function used to create Hernquist (1990) dark matter
halos for clusters of galaxies.
Hernquist models follow a spherically symmetric density profile in the form:
rho_h(r) = M_h/(2pi) * r_h/(r*(r+r_h)**3),
    where M_h is the dark matter halo mass, and r_h is a scale length.

This profile is similar to the NFW profile, except in the outermost parts
beyond the r_200 radius (also understood as the virial radius)
"""

import numpy

from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel

__all__ = ["new_hernquist_dm_model"]

class MakeHernquistDarkMatterModel(object):

    def __init__(self, targetN, convert_nbody=None, base_grid=None, rscale=1/1.695,
                 mass_cutoff=0.999, do_scale=False):
        self.targetN = targetN
        self.convert_nbody = convert_nbody
        self.rscale = rscale
        self.mass_frac = mass_cutoff
        # self.do_scale = do_scale
        # self.internal_energy = 0.25 / self.rscale
        self.base_sphere = uniform_unit_sphere(targetN, base_grid)

    def new_model(self):
        x, y, z = self.base_sphere.make_xyz()
        self.actualN = len(x)
        r = numpy.sqrt(x**2+y**2+z**2) * self.mass_frac**(1/3.)
        rtarget = self.rscale*(r**2

def new_hernquist_dm_model(number_of_particles, *list_arguments, **keyword_arguments):
    """
    Create a Hernquist dark matter model with the given numner of particles.
    Returns a set of DM particles with equal masses and positions distributed
    to fit a Hernquist distribution model. The model is centered around the
    origin. Velocities are set to zero.

    :argument number_of_particles: Number of particles to include in the Hernquist sphere
    :argument convert_nbody: When given will convert the resulting set to SI units
    :argument mass_cutoff: Mass percentage inside radius of 1
    :argument do_scale: scale the result, similar to true nbody units (M=1, Q=0.25, U=-0.5)
    """
    uc = MakeHernquistDarkMatterModel(number_of_particles, *list_arguments, **keyword_arguments)
    return uc.result
