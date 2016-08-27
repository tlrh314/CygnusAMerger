import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('text', usetex=True)
import numpy

# For kernels it is important that:
# i) 'the top is flat' (around q=0)
# ii) shit is finite (truncated at a certain q). Is this compact support?
# iii) symmetry W(|r'-r|,h) = W(|r-r'|,h) (can we check this explicitly?)
# iv) looks smooth and has smooth derivatives (first and second order)

# TODO: plot first and second derivate to see smoothness of derivates
# See Price (2012) Figure 2

def gaussian_kernel(r, h):
    sigma1d, sigma2d, sigma3d = [1/numpy.sqrt(numpy.pi), 1/numpy.pi, 1/(numpy.pi*numpy.sqrt(numpy.pi))]

    d = 3  # dimensions
    return sigma3d / h**d * numpy.exp(-1.0*(r/h)**2)


def cubic_spline(q):
    """ W(|r-r'|, h) = 1/h**d w(q). This function implements w(q)
        Schoenberg B-spline, M_4, cubic spline """
    sigma1d, sigma2d, sigma3d = [2./3, 10/(7*numpy.pi), 1/numpy.pi]
    if 0 <= q < 1:
        return sigma3d * (1./4 * (2-q)**3 - (1-q)**3)
    elif 1 <= q < 2:
        return sigma3d * (1./4 * (2-q)**3)
    elif q >= 2:
        return 0

def quartic_spline(q):
    sigma1d, sigma2d, sigma3d = [1./24, 96/(1199*numpy.pi), 1/(20*numpy.pi)]
    if 0 <= q < 1:
        return sigma3d * ((3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5)
    elif 1 <= q < 2:
        return sigma3d * ((3-q)**5 - 6*(2-q)**5)
    elif 2 <= q < 3:
        return sigma3d * ((3-q)**5)
    elif q >= 3:
        return 0

def quintic_spline(q):
    sigma1d, sigma2d, sigma3d = [1./120, 7/(478*numpy.pi), 1/(120*numpy.pi)]
    if 0 <= q < 1./2:
        return sigma3d * ((5./2-q)**4 - 5*(3./2-q)**4 + 10*(1./2-q)**4)
    elif 1./2 <= q < 3./2:
        return sigma3d * ((5./2-q)**4 - 5*(3./2-q)**4)
    elif 3./2 <= q < 5./2:
        return sigma3d * ((5./2-q)**4)
    elif q >= 2:
        return 0


def spline_kernel(r, h):
    """ Used in Gadget-2. See Monaghan & Lattanzio (1985) """
    sigma3d = 8/(numpy.pi*h)
    if 0 <= r/h < 1./2:
        return sigma3d * (1 - 6*(r/h)**2 + 6*(r/h)**3)
    elif 1./2 <= r/h <= 1:
        return sigma3d * (2*(1 - r/h)**3)
    elif r/h >= 1:
        return 0


def gadget_kernel(r, h):
    """ copy-paste from Gadget2 density.c 467 onwards """
    KERNEL_COEFF_1 = 2.546479089470
    KERNEL_COEFF_2 = 15.278874536822
    KERNEL_COEFF_3 = 45.836623610466
    KERNEL_COEFF_4 = 30.557749073644
    KERNEL_COEFF_5 = 5.092958178941
    KERNEL_COEFF_6 = (-15.278874536822)
    h2 = h*h
    hinv = 1.0/h
    hinv3 = hinv*hinv*hinv
    hinv4 = hinv3 * hinv
    u = r * hinv
    if (u > 1):
        wk = 0
    elif(u < 0.5):
        wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u)
        # dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4)
    else:
        wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u)
        # dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u)
    return wk



def kernel_spline_3d(r, h):
    """ Copy-paste from P-Smac2 project_sph.c line 336 """
    r /= h

    u = 1 - r

    if (r > 1):
        return 0
    elif (r > 0.5):
        return 8 / numpy.pi * 2 * u * u * u
    else:
        return 8 / numpy.pi * (1 - 6 * r * r * u)


def psmac2_kernel_vs_gadget_kernel():
    from numpy import vectorize
    vkernel_spline_3d = vectorize(kernel_spline_3d)
    epsilon = 1
    h = 2.8 * epsilon
    r = numpy.arange(0, 2*h, 0.001)
    psmac2 = vkernel_spline_3d(r, h)
    psmac2_hank = numpy.zeros(len(r))
    gadget = numpy.zeros(len(r))
    price = numpy.zeros(len(r))
    for i, hank in enumerate(r):
        # all three differ, yay
        psmac2_hank[i] = kernel_spline_3d(hank, h)
        gadget[i] = gadget_kernel(hank, h)
        price[i] = spline_kernel(hank, h)

    pyplot.figure(figsize=(12, 9))
    pyplot.plot(r, psmac2, label="psmac")
    pyplot.plot(r, gadget, label="gadget")
    pyplot.plot(r, price, label="price")
    # pyplot.plot(r, psmac2_hank)
    pyplot.legend()
    pyplot.show()


if __name__ == "__main__":
    psmac2_kernel_vs_gadget_kernel()
    import sys; sys.exit(0)


    r = numpy.arange(0, 2.8, 0.001)
    W_of_r_cubic = numpy.zeros(len(r))
    W_of_r_quartic = numpy.zeros(len(r))
    W_of_r_quintic = numpy.zeros(len(r))
    for n, i in enumerate(r):
        W_of_r_cubic[n] = cubic_spline(i)
        W_of_r_quartic[n] = quartic_spline(i)
        W_of_r_quintic[n] = quintic_spline(i)

    pyplot.figure(figsize=(12, 9))
    pyplot.style.use(["dark_background"])
    # data_colour = (255./255, 64./255, 255./255)
    pyplot.plot(r, W_of_r_cubic, label="Cubic", c="white", lw=4)
    pyplot.plot(r, W_of_r_quartic, label="Quartic", c="white", lw=4, ls="dashed")
    pyplot.plot(r, W_of_r_quintic, label="Quintic", c="white", lw=4, ls="dotted")
    # pyplot.legend()
    # pyplot.xlabel(r"$q$")
    # pyplot.ylabel(r"$W(q)$")
    # pyplot.show()
    # pyplot.savefig("out/CubicSpline.png", dpi=600)

    epsilon = 1
    h = 2.8 * epsilon
    W_of_r = numpy.zeros(len(r))
    for n, i in enumerate(r):
        W_of_r[n] = spline_kernel(i, h)

    # pyplot.figure(figsize=(12, 9))
    # pyplot.style.use(["dark_background"])
    data_colour = (255./255, 64./255, 255./255)
    pyplot.plot(r, W_of_r, label="Gadget2 spline", c=data_colour, lw=4)
    # pyplot.legend()
    # pyplot.xlabel(r"$q$")
    # pyplot.ylabel(r"$W(q)$")
    # pyplot.show()
    # pyplot.savefig("out/Gadget2Spline.png", dpi=600)

    # pyplot.figure(figsize=(12, 9))
    # pyplot.style.use(["dark_background"])
    W_of_r = gaussian_kernel(r, 1)
    pyplot.plot(r, W_of_r, label="Gaussian", c="g", lw=4)
    pyplot.legend()
    pyplot.xlabel(r"$q$")
    pyplot.ylabel(r"$W(q)$")
    pyplot.show()
    # pyplot.savefig("out/GaussianKernel.png", dpi=600)
