"""
File: compare_kernels.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Date created: Sun Sep 04, 2016 01:38 pm
Last modified: Sun Sep 04, 2016 07:08 pm

Compare Wendland C6, Cubic Spline M4, and Gaussian kernel

"""

# For kernels it is important that:
# i) 'the top is flat' (around q=0)
# ii) shit is finite (truncated at a certain q). Is this compact support?
# iii) symmetry W(|r'-r|,h) = W(|r-r'|,h) (can we check this explicitly?)
# iv) looks smooth and has smooth derivatives (first and second order)


import numpy
import matplotlib
from matplotlib import pyplot
from plotsettings import PlotSettings
style = PlotSettings()
from macro import *

M_1_PI = 1.0/numpy.pi
M_2_SQRTPI = 2.0/numpy.sqrt(numpy.pi)


def gaussian_kernel(r=1.0, h=1.0, dim=3, norm=False):
    sigma =  0.5*M_2_SQRTPI
    if dim > 1:
        sigma *=  0.5*M_2_SQRTPI
    if dim > 2:
        sigma *=  0.5*M_2_SQRTPI

    h1 = 1.0/h
    q = r*h1

    if dim == 1:
        sigma = sigma * h1
    elif dim == 2:
        sigma = sigma * p2(h1)
    elif dim == 3:
        sigma = sigma * p3(h1)

    if norm:
        sigma = 1.0

    val = 0.0
    if ( q < 3.0):
        val = sigma * numpy.exp(-1.0*p2(q))

    return val


def cubic_spline(r=1.0, h=1.0, dim=3, norm=False):
    h1 = 1.0/h
    q = r*h1
    t = 2. - q

    if dim == 3:
        sigma = M_1_PI * p3(h1)
    elif dim == 2:
        sigma = 10*M_1_PI/7.0 * p2(h1)
    else:
        sigma = 2.0/3.0 * h1

    if q > 2.0:
        val = 0.0
    elif q > 1.0:
        val = 0.25 * p3(t)
    else:
        val = 1 - 1.5 * p2(q) * (1 - 0.5 * q)

    if norm:
        sigma = 1.0
    return sigma * val


def gadget_kernel(r, h):
    """ copy-paste from Gadget2 density.c 467 onwards """
    KERNEL_COEFF_1 = 2.546479089470    # 1 * 1 * 8/pi
    KERNEL_COEFF_2 = 15.278874536822   # 1 * 6 * 8/pi
    KERNEL_COEFF_3 = 45.836623610466   # 3 * 6 * 8/pi
    KERNEL_COEFF_4 = 30.557749073644   # 3 * 6 * 8/pi * 2/3
    KERNEL_COEFF_5 = 5.092958178941    # 2 * 1 * 8/pi
    KERNEL_COEFF_6 = (-15.278874536822)
    hinv = 1.0/h
    hinv3 = hinv*hinv*hinv
    hinv4 = hinv3 * hinv
    u = r * hinv
    if (u > 1.0):
        wk = 0
    elif(u < 0.5):
        wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u)
        # dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4)
    else:
        wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u)
        # dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u)
    return wk


def wendland_c6(r=1.0, h=1.0, dim=3, norm=False):
    """ Toycluster sph.c 426 (Donnert 2014, 2016 in prep),
        also see Dehnen & Aly (2012) """
    h1 = 1.0/h
    q = r*h1
    t = 1.0 - q

    if dim == 1:
        print "Wendland C6 not implemented for 1 dimension"
    elif dim == 2:
        print "Wendland C6 not implemented for 2 dimensions"
    elif dim == 3:
        sigma = 1365.0/(64*numpy.pi) * p3(h1)

    val = 0.0

    if q < 1.0:
         val = t*t*t*t*t*t*t*t*(1+8*q + 25*p2(q) + 32*p3(q))

    if norm:
        sigma = 1.0
    return sigma * val


def compare_kernels(norm=False):
    h = 1.0
    radii = numpy.arange(0, 3*h, 0.001)

    pyplot.figure(figsize=(12, 9))

    pysph = False
    if pysph:
        import pysph
        from pysph.base import kernels
        m4 = pysph.base.kernels.CubicSpline(dim=3)
        gauss = pysph.base.kernels.Gaussian(dim=3)
        wc2 = pysph.base.kernels.WendlandQuintic(dim=3)
        pyplot.plot(radii/h, [gauss.kernel(rij=r, h=h) for r in radii], label="G")
        pyplot.plot(radii/h, [m4.kernel(rij=r, h=h) for r in radii], label="M4")
        pyplot.plot(radii/h, [wc2.kernel(rij=r, h=h) for r in radii], label="WC2")

    gauss = numpy.zeros(len(radii))
    spline = numpy.zeros(len(radii))
    gadget = numpy.zeros(len(radii))
    wendland = numpy.zeros(len(radii))
    for i, r in enumerate(radii):
        # NB domain 0 <= q <= 2
        gauss[i] = gaussian_kernel(r, h, norm=norm)
        spline[i] = cubic_spline(r, h, norm=norm)
        # NB domain 0 <= q <= 1
        # gadget[i] = gadget_kernel(r/2, h)/8
        if norm:
            wendland[i] = wendland_c6(r/2, h, norm=norm)
        else:
            wendland[i] = wendland_c6(r/2, h, norm=norm)/8

    pyplot.plot(radii/h, gauss, c="k", ls=":", label="Gaussian")
    pyplot.plot(radii/h, spline, c="k", ls="--", label="Cubic Spline")
    # pyplot.plot(radii/h, gadget, c="magenta", label="Gadget-2")
    pyplot.plot(radii/h, wendland, c="k", ls="-", label="Wendland C6")
    pyplot.xlabel(r"r / h", fontsize=28)
    pyplot.ylabel(r"w (q)", fontsize=28)
    # pyplot.yscale("log")  # makes compact support clearly visible :)
    pyplot.legend(fontsize=22)
    pyplot.show()


if __name__ == "__main__":
    compare_kernels()
    compare_kernels(norm=True)
