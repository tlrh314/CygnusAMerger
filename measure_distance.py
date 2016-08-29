import os
import time
import pyds9
import numpy
import scipy
from scipy import ndimage
import peakutils

from plotsettings import PlotSettings
style = PlotSettings()

from astropy.io import fits
from deco import concurrent, synchronized

from cluster import NumericalCluster
from ioparser import SimulationOutputParser


#@concurrent(processes=4)
def plot_xray(sim, snapnr, smooth=False):
    """ Create plot of X-ray Surface Brightness """

    print "Checking snapshot:", snapnr
    # xpix = numpy.arange(0, sim.xlen, 1, dtype=numpy.int)
    # ypix = numpy.arange(0, sim.ylen, 1, dtype=numpy.int)

    # The sum of the snapshot in the fits cube of the physical/darkmatter
    # density along the x-axis/y-axis also gives cluster centroids
    # xsum = numpy.sum(snapshot, axis=0)
    # ysum = numpy.sum(snapshot, axis=1)
    # pyplot.figure()
    # pyplot.plot(xpix, xsum, label="x")
    # pyplot.plot(ysum, ypix, label="y")
    # pyplot.legend()
    # pyplot.show()
    # exit(0)

    # Find the cluster centroids by taking the median values of the DM positions
    # along the merger-axis x, and along the impact-param axis x
    sim.read_gadget(snapnr)
    dmpos = sim.numerical.raw_data.pos[sim.numerical.raw_data.Ngas:
                                       sim.numerical.raw_data.N]
    xdm = dmpos[:,0]
    ydm = dmpos[:,1]
    zdm = dmpos[:,2]

    xhist, xbins = numpy.histogram(xdm, bins=sim.xlen)
    xcenters = xbins[:-1] + 0.5 * (xbins[1:] - xbins[:-1])
    yhist, ybins = numpy.histogram(ydm, bins=sim.ylen)
    ycenters = ybins[:-1] + 0.5 * (ybins[1:] - ybins[:-1])
    fig = pyplot.figure(figsize=(16,6))
    ax = fig.add_subplot(121)
    ax.plot(xcenters, xhist)
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(ybins[0], ybins[-1])
    ax.set_aspect("equal")
    ax = fig.add_subplot(122)
    ax.plot(yhist, ycenters)
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(ybins[0], ybins[-1])
    ax.set_aspect("equal")
    pyplot.show()

    H, xbins, ybins = numpy.histogram2d(ydm, xdm, bins=(sim.xlen, sim.ylen))
    fig = pyplot.figure(figsize=(14, 6))
    ax = fig.add_subplot(131)
    im = pyplot.imshow(H, interpolation="nearest", origin="low",
                extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]])
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(ybins[0], ybins[-1])
    ax = fig.add_subplot(132)
    X, Y = numpy.meshgrid(xbins, ybins)
    ax.pcolormesh(X, Y, H)
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(ybins[0], ybins[-1])
    ax.set_aspect("equal")
    ax = fig.add_subplot(133)
    im = matplotlib.image.NonUniformImage(ax, interpolation="bilinear")
    xcenters = xbins[:-1] + 0.5 * (xbins[1:] - xbins[:-1])
    ycenters = ybins[:-1] + 0.5 * (ybins[1:] - ybins[:-1])
    im.set_data(xcenters, ycenters, H)
    ax.images.append(im)
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(ybins[0], ybins[-1])
    ax.set_aspect("equal")

    pyplot.show()
    exit(0)

    boxhalf = sim.pixelscale*sim.xlen/2

    left = xdm[numpy.where(xdm < boxhalf)]
    right = xdm[numpy.where(xdm > boxhalf)]
    bottom = ydm[numpy.where(ydm < boxhalf)]
    top = ydm[numpy.where(ydm > boxhalf)]

    x1 = numpy.median(left)/sim.pixelscale
    x2 = numpy.median(right)/sim.pixelscale
    y = numpy.median(numpy.concatenate((bottom,top)))/sim.pixelscale
    print "{0}, {1}, {2}".format(snapnr, (x1,y), (x2,y))

    pyplot.figure()
    xray = numpy.log10(snapshot.clip(min=2e-8, max=0.02))
    if smooth:
        xray = scipy.ndimage.filters.gaussian_filter(
            xray, order=0, sigma=9*sim.pixelscale)
    pyplot.imshow(xray, cmap="gray", origin="lower")
    pyplot.plot(x1, y, "ro", lw=20)
    pyplot.plot(x2, y, "bo", lw=20)
    pyplot.xlim(350, 750)
    pyplot.ylim(412, 612)
    pyplot.savefig(sim.outdir+"xray_dmmedian_{0:03d}.png".format(snapnr))
    pyplot.close()
    #pyplot.show()
    return

    # Now find the peaks in the median of the dark matter positions
    xpeaks = peakutils.indexes(xsum, min_dist=100/sim.pixelscale)
    ypeaks = peakutils.indexes(ysum, min_dist=100/sim.pixelscale)

    found_peaks = False
    if len(xpeaks) == 2 and len(ypeaks) == 1:
        # Further optimize peakfinding by interpolating
        try:
            xpeaks = peakutils.interpolate(range(0, len(xsum)), xsum, ind=xpeaks)
            ypeaks = peakutils.interpolate(range(0, len(ysum)), ysum, ind=ypeaks)
        except RuntimeError as err:
            if "Optimal parameters not found: Number of calls to function has reached" in str(err):
                print "peakutils.interpolate broke"
            else:
                raise
        found_peaks = True
        CygA = xpeaks[0], ypeaks[0]
        CygB = xpeaks[1], ypeaks[0]
        distance = numpy.sqrt((CygA[0]-CygB[0])**2+(CygA[1]-CygB[1])**2)*sim.pixelscale
        print "Success: found 2 xpeaks, and 1 ypeak!"
        print "  CygA: (x, y) = {0}".format(CygA)
        print "  CygB: (x, y) = {0}".format(CygB)
        print "  distance     = {0:.2f}".format(distance)
        print
    #elif len(xpeaks) == 3 and len(ypeaks) == 1:
    #    found_peaks = True
    #    CygA = xpeaks[0], ypeaks[0]
    #    CygBx = numpy.argmax(xpeaks)
    #    CygB = xpeaks[CygBx], ypeaks[0]
    #    distance = numpy.sqrt((CygA[0]-CygB[0])**2+(CygA[1]-CygB[1])**2)*sim.pixelscale
    #    print "Problem: found 3 xpeaks, and 1 ypeak!"
    #    print "Assuming CygB is at xmax"
    #    print "  CygA: (x, y) = {0}".format(CygA)
    #    print "  CygB: (x, y) = {0}".format(CygB)
    #    print "  distance     = {0:.2f}".format(distance)
    #    print
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

    axx.text(0.5, 1.01, "Summed Dark Matter Density X", ha="center", va="bottom",
            transform=axx.transAxes, fontsize=16)
    axx.plot(xpix, xsum)
    axx.plot(xpeaks, xsum[(xpeaks+0.5).astype(int)], "ro")
    #axx.set_xlim(CygA[0]-1000/sim.pixelscale, CygB[0]+750/sim.pixelscale)
    #axx.set_xlim(800, 1400)
    #axx.set_ylim(0, 5)
    #tickrange=range(800, 1401, 200)
    #axx.set_xticks(tickrange)
    #axx.set_xticklabels(["{0}".format(t) for t in tickrange])

    axy.text(1.01, 0.5, "Summed Dark Matter Density Y", ha="left", va="center",
            transform=axy.transAxes, fontsize=16, rotation=-90)
    axy.plot(ysum, ypix)
    axy.plot(ysum[(ypeaks+0.5).astype(int)], ypeaks, "ro")
    #axy.set_xlim(0, 5)
    #axy.set_ylim(CygA[1]-1000/sim.pixelscale, CygA[1]+1000/sim.pixelscale)
    #axy.set_ylim(824, 1224)

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
    ax.text(0.5, 0.95, "X-ray Surface Brightness", ha="center", va="top",
            color="white", transform=ax.transAxes)
    xray = numpy.log10(snapshot.clip(min=2e-8, max=0.02))
    if smooth:
        xray = scipy.ndimage.filters.gaussian_filter(
            xray, order=0, sigma=9*sim.pixelscale)
    ax.imshow(xray, cmap="spectral", origin="lower")
    #pyplot.show()
    pyplot.savefig(sim.outdir+"xray-density_{0:03d}.png".format(snapnr))
    pyplot.close()


#@synchronized
def parallel_plot(sim, xraycube ):
    sim.read_smac(xraycube)
    for snapnr in range(sim.nsnaps):
        plot_xray(sim, snapnr)
        break


if __name__ == "__main__":
    timestamp = "20160727T1105"
    sim = SimulationOutputParser("/Volumes/SURFlisa", timestamp)

    xraycube = "xray-surface-brightness_projection-z.fits.fz"
    tspeccube = "analysis/temperature-emission-weighted_projection-z.fits.fz"

    parallel_plot(sim, xraycube)

    os.chdir(outdir)
    os.system('ffmpeg -y -r 8 -i "xray-density_%3d.png" -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s "2000:2000" -an "xray-density.mp4"')
