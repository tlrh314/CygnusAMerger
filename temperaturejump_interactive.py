import os
import time
import aplpy
import numpy
import scipy
from scipy import ndimage
from scipy.ndimage import filters

import matplotlib
from matplotlib import pyplot
from matplotlib import animation
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
from ioparser import parse_toycluster_parms
from macro import print_progressbar
from stitch import helix_tables

from macro import *

def find_snapshot(rundir):
    """ Find the most likely snapshot number based on the X-ray observation.

    From the X-ray observation we find the two peaks of emission are separated
    by 700 seconds of arc. """

    # For time counter
    gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
    TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

    # toyclusterparms = parse_toycluster_parms(rundir+"ICs/ic_both_hybrid.par")
    # particles = toyclusterparms['Ntotal']

    projection = "z"
    projection = "_projection-{0}".format(projection)
    option = "xray-surface-brightness"

    xraycube = rundir+"analysis/"+option+projection+".fits.fz"
    print xraycube
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
        cmap = helix_tables(module.strip(), flag.strip())

        # To see entire header, including comments starting with "/"
        # print line

    print "Effect_Module", module
    print "Effect_Flag", flag

    number_of_snapshots = header['NAXIS3']
    xlen = int(header['NAXIS1'])
    ylen = int(header['NAXIS2'])
    pixelscale = float(scale)/int(xlen)

    print "# snapshots =", number_of_snapshots
    print "BoxSize     =", scale
    print "xlen        =", xlen
    print "ylen        =", ylen
    print "kpc/pixel   =", pixelscale

    # https://stackoverflow.com/questions/9111711
    neighborhood_size = 2**(-5)*xlen
    threshold = 0

    fig = pyplot.figure(figsize=(12, 9))
    fig.gca().set_xlim(0, 0.75*xlen)
    fig.gca().set_ylim(0, 0.75*ylen)
    im = pyplot.imshow(numpy.log(data[0][0.125*xlen:0.875*xlen, 0.125*ylen:0.875*ylen]), cmap=cmap)
    hank = fig.gca().plot(0, 0, "ro")
    scale_text = "[{0:.1f} Mpc]^2".format(float(scale)/1000)
    pad=4
    pyplot.text(xlen-pad, pad, scale_text, color="white", size=42,
                horizontalalignment="right", verticalalignment="bottom")

    def init():
        print "hank"
        im.set_data(numpy.zeros((xlen, ylen)))
        global hank
        hank = fig.gca().plot(0, 0, "ro")
        for n in range(len(hank)):
            hank.pop(n).remove()
        return [im]

    def animate(i):
        global hank
        for n in range(len(hank)):
            hank.pop(n).remove()
        snapshot = numpy.log(data[i])
        #borderless_snapshot = snapshot[]
        snapshot = snapshot[0.125*xlen:0.875*xlen, 0.125*ylen:0.875*ylen]
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

        print "Snapshot:", i
        for a,b in zip(x,y):
            print "  (x, y) = ({0:3d}, {1:3d})".format(a,b)
        if len(x) == 2:
            distance = numpy.sqrt(p2(x[0]-x[1]) + p2(y[0]-y[1]))
            print "  distance = {0:<5.2f}".format(distance)
            print "  distance = {0:<5.2f} kpc".format(distance*pixelscale)
        im.set_data(snapshot)
        if i < number_of_snapshots-1:
            hank = fig.gca().plot(x, y, "ro")
        return [im]

    anim = animation.FuncAnimation(fig, animate, init_func=init,
        frames=number_of_snapshots, interval=500, repeat_delay=5000)
    pyplot.show()

    return

    pad = 4  # number of pixels padding for text placement (of tile titles)
    for n in range(number_of_snapshots):
        # Set up four-panel plot, stitched together
        pyplot.figure(figsize=(16,16))
        pyplot.style.use(["dark_background"])
        gs1 = gridspec.GridSpec(2, 2)
        #gs1 = gridspec.GridSpec(1, 1)
        gs1.update(wspace=0, hspace=0) # remove spacing between axes.

        for i in range(4):  # gridspec indexes start at 0
        #for i in range(1):  # gridspec indexes start at 0
            ax = pyplot.subplot(gs1[i])
            pyplot.axis('on')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            # NB: breaks when using pyplot.show(), save only!
            ax.set_aspect('equal')

            # Plot every panel. TODO: change way log/sqrt scheme is chosen
            # Now all properties look good in log stretch, but DM density
            # seems to look better using sqrt scaling!
            if "dm" not in chosen[i]:
                ax.imshow(numpy.log(data[n]), cmap=cmap[i])
            else:
                ax.imshow(numpy.sqrt(data[n]), cmap=cmap[i])
            # Tile text: name of physical property
            pyplot.text(pad if i%2==0 else xlen-pad, pad,
                        names[i], color="white", size=42,
                        horizontalalignment="left" if i%2==0 else "right",
                        verticalalignment="top")
            # pyplot.text(pad if i%2==0 else xlen-pad, ylen-pad,
            #             names[i], color="white", size=42,
            #             horizontalalignment="left" if i%2==0 else "right",
            #             verticalalignment="bottom")
        # Image scale (lives in lower right corner of tile 3)
        pyplot.text(xlen-pad, ylen-pad, scale, color="white", size=42,
                    horizontalalignment="right", verticalalignment="bottom")

        # pyplot.suptitle("T = {0:05.2f} Myr".format(0.05*n),
        pyplot.suptitle("T = {0:04.2f} Gyr".format(TimeBetSnapshot*n),
            color="white", size=52, y=0.95)
        pyplot.tight_layout()
        pyplot.savefig(rundir+"out/snapshot{0}_{1:03d}.png".format(projection, n))
        pyplot.close()

        print_progressbar(n, number_of_snapshots)


if __name__ == "__main__":
    timestamp = "20160707T0034"
    rundir = "../runs/{0}/".format(timestamp)

    if not (os.path.isdir(rundir+"out") or os.path.exists(rundir+"out")):
        os.mkdir(rundir+"out")

    find_snapshot(rundir)
