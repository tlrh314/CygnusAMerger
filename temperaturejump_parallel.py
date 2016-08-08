import os
import argparse
import numpy
import scipy
from scipy import ndimage
from scipy.ndimage import filters

from deco import concurrent, synchronized

import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})
pyplot.rcParams.update({"xtick.major.size": 8})
pyplot.rcParams.update({"xtick.minor.size": 4})
pyplot.rcParams.update({"ytick.major.size": 8})
pyplot.rcParams.update({"ytick.minor.size": 4})
pyplot.rcParams.update({"xtick.major.width": 4})
pyplot.rcParams.update({"xtick.minor.width": 2})
pyplot.rcParams.update({"ytick.major.width": 4})
pyplot.rcParams.update({"ytick.minor.width": 2})
pyplot.rcParams.update({"xtick.major.pad": 8})
pyplot.rcParams.update({"xtick.minor.pad": 8})
pyplot.rcParams.update({"ytick.major.pad": 8})
pyplot.rcParams.update({"ytick.minor.pad": 8})
pyplot.rcParams.update({"legend.loc": "best"})
pyplot.rcParams.update({"figure.autolayout": True})

from astropy.io import fits

from ioparser import parse_gadget_parms
from rotate import helix_tables
from macro import *


@concurrent(processes=4)
def plot_snapshot(rundir, data, n, scale, cmap, TimeBetSnapshot, plot=True):
    if plot:
        print "Plotting xray surface brightness for snapshot", n

    xlen, ylen = data.shape
    pixelscale = float(scale)/int(xlen)
    if plot:
        scale_text = "[{0:.1f} Mpc]^2".format(float(scale)/1000)

        fig = pyplot.figure(figsize=(12, 9))
        fig.gca().set_xlim(0, 0.75*xlen)
        fig.gca().set_ylim(0, 0.75*ylen)

        pad = 4
        pyplot.text(0.75*xlen-pad, pad, scale_text, color="white", size=18,
                    horizontalalignment="right", verticalalignment="bottom")
        pyplot.text(pad, pad, "T = {0:04.2f} Gyr".format(TimeBetSnapshot*n),
                    color="white", size=18, ha="left", va="bottom")

    # https://stackoverflow.com/questions/9111711
    neighborhood_size = int(2**(-5)*xlen)
    threshold = 0

    snapshot = numpy.log10(data.clip(min=1e-12))
    # Cut borders because the algorithm finds peaks at the borders
    snapshot = snapshot[int(0.125*xlen):int(0.875*xlen),
                        int(0.125*ylen):int(0.875*ylen)]

    data_max = scipy.ndimage.filters.maximum_filter(snapshot, neighborhood_size)
    maxima = (snapshot == data_max)
    data_min = scipy.ndimage.filters.minimum_filter(snapshot, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = scipy.ndimage.label(maxima)
    slices = scipy.ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2
        y.append(y_center)

    # for a,b in zip(x,y):
    #     print "  (x, y) = ({0:3d}, {1:3d})".format(a,b)
    if len(x) == 2:
        distance = numpy.sqrt(p2(x[0]-x[1]) + p2(y[0]-y[1]))
        # print "  distance = {0:<5.2f}".format(distance)
        # print "  distance = {0:<5.2f} kpc".format(distance*pixelscale)
    else:
        distance = numpy.nan

    if plot:
        pyplot.suptitle("distance = {0:<9.1f} kpc".format(distance*pixelscale),
                        color="white", size=32, y=0.95)

        pyplot.imshow(snapshot, cmap=cmap)
        pyplot.plot(x, y, "ro")

        pyplot.savefig(rundir+"out/distance_{0:03d}.png".format(n))
        pyplot.close()

    return distance*pixelscale


@synchronized
def generate_distance_plots(rundir):
    # For time counter
    gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
    TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

    # toyclusterparms = parse_toycluster_parms(rundir+"ICs/ic_both_hybrid.par")
    # particles = toyclusterparms['Ntotal']

    projection = "z"
    projection = "_projection-{0}".format(projection)
    option = "xray-surface-brightness"

    xraycube = rundir+"analysis/"+option+projection+".fits.fz"
    print "Simulation surface brightness:", xraycube
    with fits.open(xraycube) as f:
        header = f[0].header
        data = f[0].data

        # Set colourtable based on P-Smac module/flag (in fits header)
        for line in repr(header).split("\n"):
            if "Effect_Module" in line:
                module = line.strip().split("=")[-1].strip().split("/")[0]
            if "Effect_Flag" in line:
                flag = line.strip().split("=")[-1].strip().split("/")[0]
            if "XYSize" in line:  # Could also be obtained from gadgetparms
                scale = line.strip().split("=")[-1].strip().split("/")[0]
        cmap = helix_tables(int(module.strip()), int(flag.strip()))

    number_of_snapshots = header['NAXIS3']
    distance = numpy.zeros(number_of_snapshots)

    for n in range(number_of_snapshots):
        distance[n] = plot_snapshot(rundir, data[n], n, scale, cmap, TimeBetSnapshot)

    return distance, scale


@concurrent(processes=4)
def plot_temperature_along_merger_axis(bestsnap, distance, sigma=None):
    if sigma:
        print "Plotting 2D Gaussian convolved (sigma = {0}) temperature for snapshot {1}"\
            .format(sigma, bestsnap)
    else:
        print "Plotting temperature for snapshot", bestsnap

    option = "temperature-emission-weighted"
    temcube = rundir+"analysis/"+option+"_projection-z.fits.fz"
    option = "temperature-spectroscopic"
    tspeccube = rundir+"analysis/"+option+"_projection-z.fits.fz"
    with fits.open(temcube) as f:
        temheader = f[0].header
        temdata = f[0].data[bestsnap]

    with fits.open(tspeccube) as f:
        tspecheader = f[0].header
        tspecdata = f[0].data[bestsnap]

    kB_kev_per_kelvin = 8.6173427909e-08
    xlen, ylen = temdata.shape

    pixelscale = float(scale)/int(xlen)
    if sigma:
        sigma = 2*sigma/pixelscale
        temdata = scipy.ndimage.filters.gaussian_filter(
            temdata, order=0, sigma=sigma)
        tspecdata = scipy.ndimage.filters.gaussian_filter(
            tspecdata, order=0, sigma=sigma)

    fig, (ax0, ax1) = pyplot.subplots(1, 2, figsize=(16, 8))
    pyplot.suptitle("core separation = {0:<9.1f} kpc".format(distance))
    ax0.set_title(r"$T_{\rm em}$")
    tem_mergeraxis = temdata[ylen/2-1:ylen/2+1,
                             int((2**-3)*ylen):int((1-(2**-3))*ylen)]
    ax0.imshow((temdata*kB_kev_per_kelvin))
    #ax0.imshow(tem_mergeraxis.)
    # ax0.set_xlim(0, xlen)
    # ax0.set_ylim(0, ylen)
    ax1.set_title(r"$T_{\rm spec}$")
    ax1.imshow((tspecdata*kB_kev_per_kelvin))
    #pyplot.show()
    pyplot.close()

    # temdata*kB_kev_per_kelvin
    fig, (ax0, ax1) = pyplot.subplots(1, 2, figsize=(16, 8))
    pyplot.suptitle("core separation = {0:<9.1f} kpc".format(distance))
    ax0.set_xlabel("Pixel")
    ax0.set_ylabel(r"$T_{\rm em}$ (keV)")
    ax0.set_ylim(0, 4)
    ax0.plot(tem_mergeraxis[0]*kB_kev_per_kelvin)
    tspec_mergeraxis = tspecdata[ylen/2-1:ylen/2+1,
                                 int((2**-3)*ylen):int((1-(2**-3))*ylen)]
    ax1.set_xlabel("Pixel")
    ax1.set_ylabel(r"$T_{\rm spec}$ (keV)")
    ax1.set_ylim(0, 4)
    ax1.plot(tspec_mergeraxis[0]*kB_kev_per_kelvin)
    pyplot.savefig(rundir+"out/temperaturejump{1}_{0:03d}.png".format(bestsnap, "_blur" if sigma else ""))
    pyplot.savefig(rundir+"out/distance_{0:03d}.png".format(n))
    pyplot.close()
    #pyplot.show()

    # Sarazin 2012 Figure with Suzaku temperature jump
    # arcmin, keV
    # x, y         error x    error y
    # -2, 5        -1 +1      -1 +1
    # 0, 5.6       -1 +1      -0.7 +0.3
    # 2, 6.2       -1 +1      -1 +1
    # 4, 7         -1 +1      -2 +2
    # 6, 8.2       -1 +1      -1, 1.5
    # 8, 8.5       -1 +1      -1.3 +1.6
    # 10, 7.2      -1 +1      - 1.3 +1.1
    # 12, 6.4      -1 +1      -0.8 +0.7

@synchronized
def generate_temperature_plot(rundir, distance, scale, sigma=None):
    print "Distance matrix:\n", distance
    bestsnap = numpy.intersect1d(numpy.where(distance > 650), numpy.where(distance < 750))[0]

    for bestsnap in range(bestsnap+3):
        plot_temperature_along_merger_axis(bestsnap, distance[bestsnap], sigma)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        description="Calculate distances between haloes.")
    parser.add_argument("-s", "--sigma", dest="sigma", nargs=1, default=None,
        help="2D Gaussian to convolve Tem/Tspec with to mimic Suzaku SPF.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="simulationID", nargs=1,
        help="string of the Simulation ID")

    return parser


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    simulationID = arguments.simulationID[0]
    sigma = int(arguments.sigma[0]) if arguments.sigma else None

    rundir = "../runs/{0}/".format(simulationID)

    if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
        os.mkdir(rundir+"out")

    distance, scale = generate_distance_plots(rundir)
    generate_temperature_plot(rundir, distance, scale, sigma)
