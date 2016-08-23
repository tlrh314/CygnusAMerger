import os
import argparse
from collections import OrderedDict
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
from suzaku_temperaturejump import plot_sarazin_suzaku
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


# -----------------------------------------------------------------------------
def new_argument_parser():
    description="Plot T-jump (keV) along 90 deg wedge around merger axis."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s", "--sigma", dest="sigma", nargs=1, default=None,
        help="2D Gaussian to convolve Tem/Tspec with to mimic Suzaku SPF.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="simulationID", nargs=1,
        help="string of the Simulation ID")

    return parser


def create_T_jump(rundir, cube_number, (cygA_center_x, cygA_center_y),
                  (cygB_center_x, cygB_center_y)):
    x, y = [cygA_center_x, cygB_center_x], [cygA_center_y, cygB_center_y]

    # First show the X-ray surface brightness, and indicate cluster centers
    option = "xray-surface-brightness"
    xraycube = rundir+"analysis/"+option+"_projection-z.fits.fz"
    with fits.open(xraycube) as f:
        xrayheader = f[0].header
        xraydata = f[0].data[cube_number]

        # Set colourtable based on P-Smac module/flag (in fits header)
        for line in repr(xrayheader).split("\n"):
            if "Effect_Module" in line:
                module = line.strip().split("=")[-1].strip().split("/")[0]
            if "Effect_Flag" in line:
                flag = line.strip().split("=")[-1].strip().split("/")[0]
            if "XYSize" in line:  # Could also be obtained from gadgetparms
                scale = line.strip().split("=")[-1].strip().split("/")[0]
        cmap = helix_tables(int(module.strip()), int(flag.strip()))

    xlen, ylen = xraydata.shape
    pixelscale = float(scale)/int(xlen)

    distance = numpy.sqrt(p2(x[0]-x[1]) + p2(y[0]-y[1]))
    print "Distance      = {0:<5.2f}".format(distance)
    print "Distance      = {0:<5.2f} kpc".format(distance*pixelscale)

    pyplot.figure(figsize=(12, 9))
    pyplot.gca().set_xlim(xlen*2**-2, xlen*(1-2**-2))
    pyplot.gca().set_ylim(ylen*2**-2, ylen*(1-2**-2))
    pyplot.imshow(numpy.log10(xraydata.clip(min=1e-12)))
    pyplot.plot(x, y, "o", ms=int(2**(-5)*xlen), fillstyle="none",
                mfc=None, mec="k", mew=2.0)
    pyplot.savefig(rundir+"out/xray_{0:03d}.png".format(cube_number))

    # Move to temperature
    kB_kev_per_kelvin = 8.6173427909e-08
    option = "temperature-emission-weighted"
    temcube = rundir+"analysis/"+option+"_projection-z.fits.fz"
    with fits.open(temcube) as f:
        temheader = f[0].header
        temdata = f[0].data[cube_number]

    # Convolve with 2D Gaussian
    # sigma=120
    # sigma = sigma/pixelscale
    # temdata = scipy.ndimage.filters.gaussian_filter(
    #     temdata, order=0, sigma=sigma)


    pyplot.figure(figsize=(12,9))
    pyplot.gca().set_xlim(xlen*2**-2, xlen*(1-2**-2))
    pyplot.gca().set_ylim(ylen*2**-2, ylen*(1-2**-2))
    pyplot.imshow(temdata.clip(min=1e-12)*kB_kev_per_kelvin)
    pyplot.plot(x, y, "o", ms=int(2**(-5)*xlen), fillstyle="none",
                mfc=None, mec="k", mew=2.0)

    ycenter = int((y[0]+y[1])/2)
    mergeraxis_pixels = numpy.arange(x[0]-50, x[1]+1+50, 4, dtype=numpy.int)
    #mergeraxis_pixels = numpy.arange(x[1]+1+5, x[0]-50, -3, dtype=numpy.int)
    temperature_along_mergeraxis = numpy.zeros(len(mergeraxis_pixels))
    temperature_along_mergeraxis_std = numpy.zeros(len(mergeraxis_pixels))
    for j, i in enumerate(mergeraxis_pixels[:-1], 1):
        pyplot.plot(i, ycenter, "wo", mfc="w", mew=0.0)
        pyplot.plot(i, ycenter+j, "ko", mfc="k", mew=0.0)
        pyplot.plot(i, ycenter-j, "ko", mfc="k", mew=0.0)
        T_at_x = temdata.clip(min=1e-12)[ycenter-j:ycenter+j, i]
        temperature_along_mergeraxis[j] = numpy.mean(T_at_x)*kB_kev_per_kelvin
        temperature_along_mergeraxis_std[j] = numpy.std(T_at_x)*kB_kev_per_kelvin

    temperature_along_mergeraxis[0] = numpy.nan  # because j=0 slice is empty
    temperature_along_mergeraxis_std[0] = numpy.nan  # because j=0 slice is empty
    pyplot.savefig(rundir+"out/wedge_{0:03d}.png".format(cube_number))

    # Get temperature plot along the merger axis
    pyplot.figure(figsize=(12,9))
    pyplot.errorbar((mergeraxis_pixels-ycenter)*pixelscale,
    #pyplot.errorbar((mergeraxis_pixels[::-1]-ycenter)*pixelscale,
                    temperature_along_mergeraxis,
                    yerr=temperature_along_mergeraxis_std,
                    marker="o", ls="", ms=4, elinewidth=2)

    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Black points.
    r_arcmin_Suzaku = numpy.array([-2, -0.1, 2, 4, 5.9, 8, 10, 12.3])
    r_error_plus = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    r_error_min = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    T_kev_Suzaku = numpy.array([5, 5.5, 6.2, 7, 8.2, 8.5, 7.2, 6.4])
    T_error_plus = numpy.array([1, 0.3, 1, 2, 1.5, 1.6, 1.1, 0.7])
    T_error_min = numpy.array([1, 0.7, 1, 2, 1, 1.3, 1.3, 0.8])
    pyplot.errorbar(r_arcmin_Suzaku*60, T_kev_Suzaku,
        xerr=[r_error_min*60, r_error_plus*60],
        yerr=[T_error_min, T_error_plus],
        c="k", fmt="+", elinewidth=4, capsize=0)

    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Grey points.
    r_arcmin_Suzaku = numpy.array([0.1, 6.1, 11.1])
    r_error_plus = numpy.array([3.5, 1.8, 3])
    r_error_min = numpy.array([3.1, 2.8, 3.2])
    T_kev_Suzaku = numpy.array([5.6, 8.8, 7])
    T_error_plus = numpy.array([0.3, 0.7, 0.5])
    T_error_min = numpy.array([0.4, 0.7, 0.4])
    pyplot.errorbar(r_arcmin_Suzaku*60, T_kev_Suzaku,
        xerr=[r_error_min*60, r_error_plus*60],
        yerr=[T_error_min, T_error_plus],
        c="grey", fmt="+", elinewidth=4, capsize=0)
    pyplot.xlabel("Distance Along Merger Axis (kpc)")
    pyplot.ylabel(r"$T$ (keV)")
    pyplot.savefig(rundir+"out/tjump_{0:03d}.png".format(cube_number))
    pyplot.show()


if __name__ == "__main__":
    # arguments = new_argument_parser().parse_args()
    # simulationID = arguments.simulationID[0]
    # sigma = int(arguments.sigma[0]) if arguments.sigma else None

    # Visually selected from ds9 in x-ray surface brightness cube fits file
    parameters =  {
        #"20160727T1105": {"n": 23, "cygA": (505, 512), "cygB": (568, 512)},
        #"20160727T1108": {"n": 20, "cygA": (506, 512), "cygB": (562, 512)},
        #"20160727T1112": {"n": 18, "cygA": (503, 512), "cygB": (563, 512)},
        #"20160727T1116": {"n": 15, "cygA": (500, 512), "cygB": (582, 512)},
        #"20160727T1120": {"n": 14, "cygA": (501, 512), "cygB": (572, 512)},
        #"20160727T1124": {"n": 13, "cygA": (503, 512), "cygB": (568, 512)},
        #"20160727T1128": {"n": 12, "cygA": (503, 512), "cygB": (567, 512)},
        #"20160727T1132": {"n": 10, "cygA": (499, 512), "cygB": (579, 512)},
        #"20160727T1303": {"n": 9, "cygA": (502, 512), "cygB": (564, 512)},
        #"20160727T1307": {"n": 9, "cygA": (501, 512), "cygB": (570, 512)},
        #"20160727T1312": {"n": 8, "cygA": (499, 512), "cygB": (571, 512)},
        #"20160727T1316": {"n": 8, "cygA": (501, 512), "cygB": (565, 512)},
        #"20160727T1321": {"n": 7, "cygA": (499, 512), "cygB": (578, 512)},
        "20160727T1325": {"n": 6, "cygA": (500, 512), "cygB": (576, 512)}
    }
    parameters = OrderedDict(sorted(parameters.items(), key=lambda t: t[0]))

    for simulationID in parameters.keys():
        print "-"*80
        print "Creating temperature jump for SimulationID =", simulationID
        # if sigma: print "Convolving temperature map with 2D Gaussian"
        rundir = "/Volumes/Taurus/runs/{0}/".format(simulationID)
        print "Rundir =", rundir

        if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
            print "Creating output directory:", rundir+"out"
            os.mkdir(rundir+"out")

        print "Best snapshot =", parameters[simulationID]["n"]
        print "Cygnus A (x,y)=", parameters[simulationID]["cygA"]
        print "Cygnus B (x,y)=", parameters[simulationID]["cygB"]

        create_T_jump(rundir, parameters[simulationID]["n"],
                      parameters[simulationID]["cygA"],
                      parameters[simulationID]["cygB"])
        print "-"*80
