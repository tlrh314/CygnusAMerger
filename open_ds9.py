import pyds9
import numpy
import scipy
from scipy import ndimage
import peakutils

import matplotlib
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
#from matplotlib import rc
#matplotlib.rc("font", **{"family":"serif", "serif":["Computer Modern Roman"],
#                         "size":28, "weight":"bold"})
params = { "backend": "pdf",
           # "axes.labelsize": 12,
           # "text.fontsize": 12,
           # "legend.fontsize": 12,
           # "xtick.labelsize": 10,
           # "ytick.labelsize": 10,
           # "text.usetex": True,
           # "text.latex.preamble": r"\boldmath",
           # "axes.unicode_minus": True,
           "font.size": 28,
           "figure.figsize": (12,9)}
matplotlib.rcParams.update(params)

from astropy.io import fits

from deco import concurrent, synchronized
import amuse.plot as amuse_plot


#@concurrent(processes=4)
def plot_density(snapnr, snapshot, pixelscale):
    # Sum the data in x and in y direction

    print "Checking snapshot:", snapnr
    xlen, ylen = snapshot.shape
    xpix = numpy.arange(0, xlen, 1, dtype=numpy.int)
    ypix = numpy.arange(0, ylen, 1, dtype=numpy.int)
    xsum = numpy.sum(snapshot, axis=0)
    ysum = numpy.sum(snapshot, axis=1)

    # The centroids will appear as peaks in the density sum
    xpeaks = peakutils.indexes(xsum, min_dist=300/pixelscale)
    ypeaks = peakutils.indexes(ysum, min_dist=300/pixelscale)

    found_peaks = False
    if len(xpeaks) == 2 and len(ypeaks) == 1:
        # Further optimize peakfinding by interpolating
        xpeaks = peakutils.interpolate(range(0, len(xsum)), xsum, ind=xpeaks)
        ypeaks = peakutils.interpolate(range(0, len(ysum)), ysum, ind=ypeaks)
        found_peaks = True
        CygA = xpeaks[0], ypeaks[0]
        CygB = xpeaks[1], ypeaks[0]
        distance = numpy.sqrt((CygA[0]-CygB[0])**2+(CygA[1]-CygB[1])**2)*pixelscale
        print "Success: found 2 xpeaks, and 1 ypeak!"
        print "  CygA: (x, y) = {0}".format(CygA)
        print "  CygB: (x, y) = {0}".format(CygB)
        print "  distance     = {0:.2f}".format(distance)
        print
    elif len(xpeaks) == 3 and len(ypeaks) == 1:
        found_peaks = True
        CygA = xpeaks[0], ypeaks[0]
        CygBx = numpy.argmax(xpeaks)
        CygB = xpeaks[CygBx], ypeaks[0]
        distance = numpy.sqrt((CygA[0]-CygB[0])**2+(CygA[1]-CygB[1])**2)*pixelscale
        print "Problem: found 3 xpeaks, and 1 ypeak!"
        print "Assuming CygB is at xmax"
        print "  CygA: (x, y) = {0}".format(CygA)
        print "  CygB: (x, y) = {0}".format(CygB)
        print "  distance     = {0:.2f}".format(distance)
        print
    else:
        distance = numpy.nan

    fig = pyplot.figure(figsize=(12, 12))

    gs = gridspec.GridSpec(3, 3)
    axx = fig.add_subplot(gs.new_subplotspec((0, 0), colspan=2))
    axy = fig.add_subplot(gs.new_subplotspec((1, 2), rowspan=2))
    axd = fig.add_subplot(gs.new_subplotspec((1, 0), colspan=2, rowspan=2),
        sharex=axx, sharey=axy)
    axt = fig.add_subplot(gs.new_subplotspec((0, 2)))
    axt.axis("off")
    gs.update(wspace=0, hspace=0, top=0.94, left=0.15)

    axx.text(0.5, 1.01, "Summed Physical Density X", ha="center", va="bottom",
            transform=axx.transAxes, fontsize=16)
    axx.plot(xpix, xsum)
    axx.plot(xpeaks, xsum[(xpeaks+0.5).astype(int)], "ro")
    #axx.set_xlim(CygA[0]-1000/pixelscale, CygB[0]+750/pixelscale)
    axx.set_xlim(800, 1400)
    axx.set_ylim(0, 5)
    tickrange=range(800, 1401, 200)
    axx.set_xticks(tickrange)
    axx.set_xticklabels(["{0}".format(t) for t in tickrange])

    axy.text(1.01, 0.5, "Summed Physical Density Y", ha="left", va="center",
            transform=axy.transAxes, fontsize=16, rotation=-90)
    axy.plot(ysum, ypix)
    axy.plot(ysum[(ypeaks+0.5).astype(int)], ypeaks, "ro")
    axy.set_xlim(0, 5)
    #axy.set_ylim(CygA[1]-1000/pixelscale, CygA[1]+1000/pixelscale)
    axy.set_ylim(824, 1224)

    # axd.imshow(snapshot, cmap="cubehelix", origin="lower")

    if found_peaks:
        axd.plot(CygA[0], CygA[1], "ro", lw=20)
        axd.plot(CygB[0], CygB[1], "ro", lw=10)
    axd.set_axis_bgcolor("k")
    axd.set_xlabel("x [pixel]")
    axd.set_ylabel("y [pixel]")

    # info = r"\begin{tabular}{llll}"
    # info += " T & = & {0:02.2f} & [Gyr] \\\\".format(snapnr*0.025)
    # info += " R & = & {0:5.2f}  & [kpc]  \\\\".format(distance)
    # info += (" \end{tabular}")
    # axt.text(0.02, 1, info, ha="left", va="top")

    axt.text(0.02, 1, "T = {0:2.2f} [Gyr]".format(snapnr*0.025),
             ha="left", va="top")
    axt.text(0.02, 0.8, "R = {0:5.2f} [kpc]".format(distance),
             ha="left", va="top")
    if 625 < distance < 725:
        axt.text(0.02, 0.6, "CANDIDATE FOUND", ha="left", va="top", color="red")

    for ax in [axx, axy, axt]:
        for tl in ax.get_xticklabels() + ax.get_yticklabels()\
                + ax.get_xticklines()  + ax.get_yticklines():
            tl.set_visible(False)
    #pyplot.show()
    #pyplot.savefig("out/xray_{0:03d}.png".format(snapnr))
    #pyplot.close()
    return axd


def plot_xray(snapnr, snapshot, pixelscale, ax, smooth=False):
    ax.text(0.5, 0.95, "X-ray Surface Brightness", ha="center", va="top",
            color="white", transform=ax.transAxes)
    xray = numpy.log10(snapshot.clip(min=2e-8, max=0.02))
    if smooth:
        xray = scipy.ndimage.filters.gaussian_filter(
            xray, order=0, sigma=9*pixelscale)
    ax.imshow(xray, cmap="spectral", origin="lower")

#@synchronized
def read_fitsfile(rhocube, xraycube):
    # Find center pixelvalues in a histogram of the physical density
    with fits.open(rhocube) as f:
        rhoheader = f[0].header
        rhodata = f[0].data
        for line in repr(rhoheader).split("\n"):
            # if "Effect_Module" in line:
            #     module = line.strip().split("=")[-1].strip().split("/")[0]
            # if "Effect_Flag" in line:
            #     flag = line.strip().split("=")[-1].strip().split("/")[0]
            if "XYSize" in line:  # Could also be obtained from gadgetparms
                scale = line.strip().split("=")[-1].strip().split("/")[0]
        # cmap = helix_tables(int(module.strip()), int(flag.strip()))

    xlen, ylen = rhodata[0].shape
    pixelscale = float(scale)/int(xlen)
    number_of_snapshots = len(rhodata)

    with fits.open(xraycube) as f:
        xrayheader = f[0].header
        xraydata = f[0].data

    for snapnr in range(number_of_snapshots):
        ax = plot_density(snapnr, rhodata[snapnr], pixelscale)
        plot_xray(snapnr, xraydata[snapnr], pixelscale, ax)
        pyplot.savefig("out/xray-density_{0:03d}.png".format(snapnr))
        pyplot.close()


if __name__ == "__main__":
    timestamp = "20160820T0438"
    simdir = "/Volumes/SURFlisa/runs/{0}/analysis/".format(timestamp)
    rhocube = simdir+"physical-density_projection-z.fits.fz"
    xraycube = simdir+"xray-surface-brightness_projection-z.fits.fz"
    temcube = simdir+"temperature-emission-weighted_projection-z.fits.fz"
    projection = "temperature_700_EA1.fits.fz"
    projection = "xray_700_EA1.fits.fz"

    # read_fitsfile(rhocube, xraycube)
    # import sys; sys.exit(0)

    view = simdir+projection


    # ds9_targets returns None if no ds9 window is opened
    ds9_running = pyds9.ds9_targets()
    if ds9_running:
        print "ds9 is already running"
        for instance in ds9_running:
            targetid = instance.replace("DS9:", "")
            d = pyds9.DS9(targetid)
            print "  instance =", targetid
            while True:
                close = raw_input("    Close [y/n]: ")
                if close.lower() in "yn":
                    if close == "y":
                        d.set("exit")
                    break
                else:
                    print "    Invalid character"
    else:
        # List of commands: http://ds9.si.edu/doc/ref/command.html
        d = pyds9.DS9()
        print d.id
        print d.target
        d.set("file '{0}'".format(view))

        if "temperature" in view:
            d.set("scale linear")
            d.set("scale limits 2e6 2e8")
            d.set("cmap bb")
        if "xray" in view:
            d.set("scale log")
            d.set("scale limits 2e-16 0.02")
            d.set("cmap sls")
        if "density" in view:
            d.set("scale log")
            d.set("scale limits 2e-16 0.02")
            d.set("cmap sls")

        d.set("scale mode User")
        # d.set("scale scope global")  # takes forever
        # d.set("scale open")

        # seems slow until al frames have been loaded at least once
        d.set("cube interval 0.125")
        d.set("cube play")
        # d.set("cube close")
        # d.set("cube 35")

        # reg = "/Users/timohalbesma/Desktop/TemperatureJump/Simulation/"+\
        #    "20160819T2357.reg"
        # d.set("regions load '{0}'".format(reg))

        # This takes forever
        # fits = d.get_arr2np()
        # print fits.shape
        # print fits.dtype

        #pyplot.figure(figsize=(12,9))
        #pyplot.imshow(fits[0], origin="lower")
        #pyplot.show()

        # Close the ds9 window if it is still open
        raw_input("Press any key to exit")
        if pyds9.ds9_targets() and "DS9:{0} {1}".format(d.target, d.id) in pyds9.ds9_targets():
            d.set("exit")
