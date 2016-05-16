import matplotlib
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('text', usetex=True)
from matplotlib import pyplot
import numpy


def cubic_spline(q):
    sigma1d, sigma2d, sigma3d = [2./3, 10/(7*numpy.pi), 1/numpy.pi]
    if 0 <= q < 1:
        return sigma3d * (1./4 * (2-q)**3 - (1-q)**3)
    elif 1 <= q < 2:
        return sigma3d * (1./4 * (2-q)**3)
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

if __name__ == "__main__":
    r = numpy.arange(0, 2.8, 0.001)
    W_of_r = numpy.zeros(len(r))
    for n, i in enumerate(r):
        W_of_r[n] = cubic_spline(i)

    pyplot.figure(figsize=(12, 9))
    pyplot.style.use(["dark_background"])
    # data_colour = (255./255, 64./255, 255./255)
    pyplot.plot(r, W_of_r, label="Cubic spline", c="white", lw=4)
    pyplot.legend()
    pyplot.xlabel(r"$q$")
    pyplot.ylabel(r"$W(q)$")
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
    pyplot.legend()
    pyplot.xlabel(r"$q$")
    pyplot.ylabel(r"$W(q)$")
    pyplot.show()
    # pyplot.savefig("out/Gadget2Spline.png", dpi=600)
