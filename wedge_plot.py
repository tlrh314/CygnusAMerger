import os
import argparse
import numpy
import scipy
from scipy import ndimage
from scipy.ndimage import filters

from deco import concurrent, synchronized

import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
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


@concurrent(processes=8)
def plot_snapshot(rundir, data, n, scale, cmap, TimeBetSnapshot, plot=True):
    if plot:
        print "Plotting xray surface brightness for snapshot", n

    xlen, ylen = data.shape
    pixelscale = float(scale)/int(xlen)
    if plot:
        scale_text = "[{0:.1f} Mpc]^2".format(float(scale)/1000)

        fig = pyplot.figure(figsize=(12, 9))
        #fig.gca().set_xlim(0, 0.75*xlen)
        #fig.gca().set_ylim(0, 0.75*ylen)
        fig.gca().set_xlim(0, xlen)
        fig.gca().set_ylim(0, ylen)

        pad = 4
        pyplot.text(0.75*xlen-pad, pad, scale_text, color="white", size=18,
                    horizontalalignment="right", verticalalignment="bottom")
        pyplot.text(pad, pad, "T = {0:04.2f} Gyr".format(TimeBetSnapshot*n),
                    color="white", size=18, ha="left", va="bottom")

    # https://stackoverflow.com/questions/9111711
    neighborhood_size = int(2**(-4)*xlen)
    threshold = 0.75

    snapshot = numpy.log10(data.clip(min=1e-12, max=1e-5))
    # Cut borders because the algorithm finds peaks at the borders
    #snapshot = snapshot[int(0.125*xlen):int(0.875*xlen),
    #                    int(0.125*ylen):int(0.875*ylen)]

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
        y_center = (dy.start + dy.stop - 1)/2
        if x_center > int(0.125*xlen) and x_center < int(0.875*xlen) and \
                y_center > int(0.125*ylen) and y_center < int(0.875*ylen):
            x.append(x_center)
            y.append(y_center)

    # for a,b in zip(x,y):
    #     print "  (x, y) = ({0:3d}, {1:3d})".format(a,b)
    if len(x) == 2:
        distance = numpy.sqrt(p2(x[0]-x[1]) + p2(y[0]-y[1]))
        print "  distance = {0:<5.2f}".format(distance)
        print "  distance = {0:<5.2f} kpc".format(distance*pixelscale)
    elif len(x) == 3 and n < 15:
        print "Found three points. Taking max distance. Please verify!"
        # Assuming one point lies between both clusters, in shock region
        max_dist = 0
        index0 = 0
        index1 = 0
        for i in range(3):
            for j in range(i):
                distance = numpy.sqrt(p2(x[i]-x[j]) + p2(y[i]-y[j]))
                if distance > max_dist:
                    max_dist = distance
                    index0 = i
                    index1 = j

        xcenter0 = x[index0]
        ycenter0 = y[index0]
        xcenter1 = x[index1]
        ycenter1 = y[index1]
        x, y = [], []
        x.append(xcenter0)
        x.append(xcenter1)
        y.append(ycenter0)
        y.append(ycenter1)
        distance = numpy.sqrt(p2(x[0]-x[1]) + p2(y[0]-y[1]))
        print "  distance = {0:<5.2f}".format(distance)
        print "  distance = {0:<5.2f} kpc".format(distance*pixelscale)
    else:
        distance = numpy.nan

    if plot:
        pyplot.suptitle("distance = {0:<9.1f} kpc".format(distance*pixelscale),
                        color="white", size=32, y=0.95)

        pyplot.imshow(snapshot, cmap=cmap)
        pyplot.plot(x, y, "ro")

        pyplot.savefig(rundir+"out/{0:03d}_distance.png".format(n))
        pyplot.close()

    pyplot.figure(figsize=(12,9))
    pyplot.suptitle("raw")
    pyplot.gca().set_xlim(0, xlen)
    pyplot.gca().set_ylim(0, ylen)
    pyplot.plot(x, y, "o", ms=neighborhood_size, fillstyle="none", mfc=None, mec="k", mew=2.0)
    pyplot.imshow(data)
    pyplot.savefig(rundir+"out/{0:03d}_raw.png".format(n))
    pyplot.close()

    pyplot.figure(figsize=(12,9))
    pyplot.suptitle("log10, clipped")
    pyplot.gca().set_xlim(0, xlen)
    pyplot.gca().set_ylim(0, ylen)
    pyplot.plot(x, y, "o", ms=neighborhood_size, fillstyle="none", mfc=None, mec="k", mew=2.0)
    pyplot.imshow(snapshot)
    pyplot.savefig(rundir+"out/{0:03d}_log10clipped.png".format(n))
    pyplot.close()

    pyplot.figure(figsize=(12,9))
    pyplot.suptitle("max")
    pyplot.gca().set_xlim(0, xlen)
    pyplot.gca().set_ylim(0, ylen)
    pyplot.plot(x, y, "o", ms=neighborhood_size, fillstyle="none", mfc=None, mec="k", mew=2.0)
    pyplot.imshow(data_max)
    pyplot.savefig(rundir+"out/{0:03d}_max.png".format(n))
    pyplot.close()

    pyplot.figure(figsize=(12,9))
    pyplot.suptitle("min")
    pyplot.gca().set_xlim(0, xlen)
    pyplot.gca().set_ylim(0, ylen)
    pyplot.plot(x, y, "o", ms=neighborhood_size, fillstyle="none", mfc=None, mec="k", mew=2.0)
    pyplot.imshow(data_min)
    pyplot.savefig(rundir+"out/{0:03d}_min.png".format(n))
    pyplot.close()

    pyplot.figure(figsize=(12,9))
    pyplot.suptitle("max-min")
    pyplot.gca().set_xlim(0, xlen)
    pyplot.gca().set_ylim(0, ylen)
    pyplot.imshow(data_max - data_min)
    pyplot.plot(x, y, "o", ms=neighborhood_size, fillstyle="none", mfc=None, mec="k", mew=2.0)
    pyplot.savefig(rundir+"out/{0:03d}_minmax.png".format(n))
    pyplot.close()

    return distance*pixelscale


@synchronized
def find_xray_distance(rundir):
    print "-"*80
    print "Obtaining distance between subclusters in X-ray surface brightness."
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


#@concurrent(processes=8)
def temperature_wedge(rundir, n, distance, xcenter, ycenter,
        scale, TimeBetSnapshot, vel, fig, ax0, ax1):
    print "Plotting temperature for snapshot", n

    option = "temperature-emission-weighted"
    temcube = rundir+"analysis/"+option+"_projection-z.fits.fz"
    with fits.open(temcube) as f:
        temheader = f[0].header
        temdata = f[0].data[n]

    kB_kev_per_kelvin = 8.6173427909e-08
    xlen, ylen = temdata.shape
    pixelscale = float(scale)/int(xlen)

    tem_mergeraxis = temdata[ylen/2-1:ylen/2+1,
                             int((2**-3)*ylen):int((1-(2**-3))*ylen)]
    ax0.plot(tem_mergeraxis[0]*kB_kev_per_kelvin,
            #label="{0:3.1f}".format(vel))
            label="{0:3.1f} -> t={1:04.2f}; r={2:03.1f}"\
                  .format(vel, n*TimeBetSnapshot, distance))

    sigma=120
    sigma = sigma/pixelscale
    temdata = scipy.ndimage.filters.gaussian_filter(
        temdata, order=0, sigma=sigma)

    # tem_mergeraxis = temdata[ycut-1:ycut+1,xcenter[0]-20:xcenter[1]+20]
    tem_mergeraxis = temdata[ylen/2-1:ylen/2+1,
                             int((2**-3)*ylen):int((1-(2**-3))*ylen)]
    ax1.plot(tem_mergeraxis[0]*kB_kev_per_kelvin)
    #        label="{0:3.1f} -> t={1:4.1f} Gyr; r={2:3.1f} kpc"\
    #              .format(vel, n*TimeBetSnapshot, distance))

    return pixelscale

#@synchronized
def generate_temperature_plot(rundir, distance, scale, sigma=None):
    gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
    TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

    print "Distance matrix:\n", distance
    # bestsnap = numpy.intersect1d(numpy.where(distance > 650), numpy.where(distance < 750))[0]
    bestsnap = numpy.nanargmin(numpy.abs(distance-701))
    print "Bestsnap, distance =", bestsnap, distance[bestsnap]

    # snaprange = len(distance)
    snaprange = bestsnap+3 if bestsnap+3 < len(distance) else bestsnap
    for snap in range(snaprange):
        plot_temperature_along_merger_axis(rundir, snap, distance[snap], scale,
                                           TimeBetSnapshot, sigma)


def new_argument_parser():
    description="Plot T-jump (keV) along 90 deg wedge around merger axis."
    parser = argparse.ArgumentParser(description=description)
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

    print "-"*80
    print "Creating temperature jump for SimulationID =", simulationID
    if sigma: print "Convolving temperature map with 2D Gaussian"
    rundir = "../runs/{0}/".format(simulationID)
    print "Rundir =", rundir
    print "-"*80

    if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
        print "Creating output directory:", rundir+"out"
        os.mkdir(rundir+"out")


    distances, pixelscale = find_xray_distance(rundir)

    print "-"*80
    print "Distance matrix:\n", distances
    # bestsnap = numpy.intersect1d(numpy.where(distance > 650), numpy.where(distance < 750))[0]
    bestsnap = numpy.nanargmin(numpy.abs(distances-701))
    print "Bestsnap, distance =", bestsnap, distances[bestsnap]
    print "-"*80
