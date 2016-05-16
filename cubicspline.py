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


r = numpy.arange(0, 2, 0.001)
W_of_r = numpy.zeros(len(r))
for n, i in enumerate(r):
    W_of_r[n] = cubic_spline(i)

pyplot.figure(figsize=(12,9))
pyplot.style.use(["dark_background"])
# data_colour = (255./255, 64./255, 255./255)
pyplot.plot(r, W_of_r, label="Cubic spline", c="white", lw=4)
pyplot.legend()
pyplot.xlabel(r"$q$")
pyplot.ylabel(r"$W(q)$")
# pyplot.show()
pyplot.savefig("out/CubicSpline.png", dpi=600)
