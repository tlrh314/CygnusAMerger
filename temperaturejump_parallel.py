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
from suzaku_temperaturejump import plot_sarazin_suzaku


#@concurrent(processes=8)
def plot_snapshot(rundir, data, n, scale, cmap, TimeBetSnapshot, plot=False):
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

    return distance*pixelscale, x, y


#@synchronized
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
    x = numpy.zeros(number_of_snapshots, dtype=object)
    y = numpy.zeros(number_of_snapshots, dtype=object)

    for n in range(number_of_snapshots):
        distance[n], x[n], y[n] = plot_snapshot(rundir, data[n], n, scale, cmap, TimeBetSnapshot)

    return distance, x, y, scale


#@concurrent(processes=8)
def plot_temperature_along_merger_axis(rundir, bestsnap, distance,
                                       TimeBetSnapshot, sigma=None):
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


#@concurrent(processes=8)
def plot_temperature_along_merger_axis_2(rundir, n, distance, xcenter, ycenter,
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

    # temdata = numpy.log10(temdata.clip(min=1e-12))
    # temdata = temdata[int(0.125*xlen):int(0.875*xlen),
    #                    int(0.125*ylen):int(0.875*ylen)]

    # print xcenter
    # print ycenter

    # ycut = int((ycenter[0]+ycenter[1])/2)
    # tem_mergeraxis = temdata[ycut-1:ycut+1,xcenter[0]-20:xcenter[1]+20]
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
def plot_temperaturejump_bestsnap_vs_velocity():
    fig, (ax0, ax1) = pyplot.subplots(1, 2, figsize=(16, 8))

    velocity = numpy.arange(0, 1.4, 0.1)
    for i,run in enumerate(["20160727T1105", "20160727T1108", "20160727T1112"]):#,
                            #"20160727T1116", "20160727T1120", "20160727T1124",
                            #"20160727T1128"]):#, "20160727T1132", "20160727T1303",
                            #"20160727T1307", "20160727T1312", "20160727T1316",
                            #"20160727T1321", "20160727T1325"]):
        print velocity[i]

        rundir = "../runs/{0}/".format(run)

        if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
            os.mkdir(rundir+"out")

        # gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
        # TimeBetSnapshot = gadgetparms['TimeBetSnapshot']
        TimeBetSnapshot = 0.1

        distance, x, y, scale = generate_distance_plots(rundir)
        print "Distance matrix:\n", distance
        bestsnap = numpy.nanargmin(numpy.abs(distance-701))
        print "Bestsnap, distance =", bestsnap, distance[bestsnap]
        pixelscale = plot_temperature_along_merger_axis_2(
            rundir, bestsnap, distance[bestsnap], x[bestsnap], y[bestsnap], scale,
            TimeBetSnapshot, velocity[i], fig, ax0, ax1)

    pyplot.show()


    pyplot.sca(ax0)
    ax0.set_title(r"unconvolved")
    ax0.set_xlabel("Mpc")
    ax0.set_ylabel(r"$T_{\rm em}$ (keV)")
    ax0.set_ylim(0, 10)
    pyplot.xticks([x for x in range(0, 800, 100)],
        ["{0:.1f}".format(x*pixelscale/1000) for x in range(0, 800, 100)])
    ax0.legend(loc="upper left", fontsize=10)

    pyplot.sca(ax1)
    ax1.set_title(r"120 arcsec convolved")
    ax1.set_xlabel("Mpc")
    ax1.set_ylabel(r"$T_{\rm em}$ (keV)")
    ax1.set_ylim(0, 10)
    ax1.legend(loc="upper right")
    pyplot.xticks([x for x in range(0, 800, 100)],
        ["{0:.1f}".format(x*pixelscale/1000) for x in range(0, 800, 100)])
    pyplot.tight_layout()
    pyplot.setp(ax0.xaxis.get_majorticklabels(), rotation=45)
    pyplot.setp(ax1.xaxis.get_majorticklabels(), rotation=45)

    pyplot.sca(ax1)
    plot_sarazin_suzaku(convert_to_mpc=True, pixelscale=pixelscale/800)

    pyplot.savefig("out/jump_vs_ZeroEOrbitFrac.png")
    # pyplot.close()
    #pyplot.show()


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
    # arguments = new_argument_parser().parse_args()
    # simulationID = arguments.simulationID[0]
    # sigma = int(arguments.sigma[0]) if arguments.sigma else None

    # rundir = "../runs/{0}/".format(simulationID)

    # if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
    #     os.mkdir(rundir+"out")

    #distance, x, y scale = generate_distance_plots(rundir)
    #generate_temperature_plot(rundir, distance, scale, sigma)
    plot_temperaturejump_bestsnap_vs_velocity()
