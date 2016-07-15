import os
import time
import glob
import aplpy
import numpy

from deco import concurrent, synchronized

import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
from astropy.io import fits

# http://www.ifweassume.com/2014/04/cubehelix-colormap-for-python.html
# https://github.com/jradavenport/cubehelix
import cubehelix

from ioparser import parse_gadget_parms
from ioparser import parse_toycluster_parms
from macro import print_progressbar

def helix_tables(module, flag, inv=False, values=None):
    """ Port of P-Smac2/lib/idl/helix_tables.pro of P-Smac colormaps.

        @param module --> int, corresponds to P-Smac2 module
        @param flag   --> int, corresponds to P-Smac2 flag
        @param inv    --> bool, invert colormap
        @param values --> iterable of length 4, use given values instead

        Cubehelix parameters: [ start, rots, hue, gamma]
        returns a cubehelix colormap ready to use in e.g. pyplot :-) """

    nModule = 15
    nFlag = 15
    setup = numpy.zeros((4, nModule, nFlag), dtype=numpy.float)

    # std values
    for i in range(0, nModule-1):
        for j in range(0, nFlag-1):
            setup[:,i,j] = [ -0, 0, 2.5, 1 ]

    # density
    setup[:, 0, 0] = [ 3, 3, 3, 2 ]

    # velocities
    setup[:, 1, 0] = [ 1, 1.5, 2, 1 ]

    # x-rays
    setup[:, 2, 0] = [ 2, -2, 3, 1 ]

    # Temperature
    setup[:, 4, 0] = [ 0.9, 0, 1, 1 ]     # Mass Weighted
    setup[:, 4, 1] = [ 1, 0, 1, 1 ]       # Sound Speed
    setup[:, 4, 2] = [ 1, 0, 1, 1 ]       # Emission Weighted
    setup[:, 4, 3] = [ 1.3, -0.5, 2, 2 ]  # Spectroscopic

    # Pressure
    setup[:, 5, 0] = [ 0, 3, 1, 1 ]

    # magnetic field
    setup[:, 6, 0] = [ 3, 0, 3, 2]

    # Compton-Y / SZ
    setup[:, 7, 0] = [ 3, -1, 4, 1 ]
    setup[:, 7, 1] = [ 2, 1, 4, 1 ]

    # Dm density
    setup[:, 10, 0] = [ 3, -1, 4, 1 ]
    setup[:, 11, 0] = [ 1, 0, 1.4, 1 ]

    if values:
        setup[:, module, flag] = values

    # set
    start, rots, hue, gamma = setup[:, module, flag]
    # Cubehelix guy changed "hue" to "sat".
    cx = cubehelix.cmap(start=start, rot=rots, sat=hue, gamma=gamma, reverse=inv)

    return cx


@concurrent(processes=8)
def make_plot(rotation, cmap, scale, xlen, ylen, name):
    replot = True
    print "Generating plot of", rotation
    rundir = "/".join(rotation.split("/")[0:-2])


    with fits.open(rotation) as f:
        header = f[0].header
        data = f[0].data
        filename = rotation.split("/")[-1].split(".")[0]

    pad = 4  # number of pixels padding for text placement (of tile titles)
    imshowfilename = rundir+"/out/{0}.png".format(filename)
    if (os.path.isfile(imshowfilename) and os.path.exists(imshowfilename)
            and not replot):
        print "Plots exist. Skipping imshow of", filename
    else:
        pyplot.figure(figsize=(16, 16))
        ax = pyplot.gca()
        pyplot.style.use(["dark_background"])
        pyplot.axis('on')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')

        # Clip to avoid "RuntimeWarning: divide by zero encountered in log"
        ax.imshow(numpy.log10(data[0].clip(min=1e-10)), cmap=cmap)
        # Tile text: name of physical property
        pyplot.text(pad, pad,
                    name, color="white", size=42,
                    horizontalalignment="left",
                    verticalalignment="top")

        pyplot.text(xlen-pad, ylen-pad, scale, color="white", size=42,
                    horizontalalignment="right", verticalalignment="bottom")

        #pyplot.suptitle("T = {0:04.2f} Gyr".format(TimeBetSnapshot*n),
        #    color="white", size=52, y=0.95)
        pyplot.tight_layout()
        pyplot.savefig(imshowfilename)
        pyplot.close()
    return

    aplfilename = rundir+"/out/apl_{0}.pdf".format(filename)
    if (os.path.isfile(aplfilename) and os.path.exists(aplfilename)
            and not replot):
        print "Plots exist. Skipping apl show_colorscale of", filename
    else:

        f = aplpy.FITSFigure(rundir+"/rotation/{0}.fits.fz".format(filename))
        f.show_colorscale(vmin=1.0e-12,vmax=1.0e-4,stretch="log",cmap="spectral",smooth=9)
        f.save(aplfilename, dpi=300)
        # raw_input()

    print "Done generating plots of", rotation


@synchronized
def plot_all(rundir):
    # For time counter
    gadgetparms = parse_gadget_parms(rundir+"snaps/gadget2.par")
    TimeBetSnapshot = gadgetparms['TimeBetSnapshot']

    toyclusterparms = [file for file in os.listdir(rundir+"ICs") if "par" in file][0]
    toyclusterparms = parse_toycluster_parms(toyclusterparms)
    particles = toyclusterparms['Ntotal']

    if particles > 2e6:
        print "Strap in, this could take a while :-) ..."

    files = sorted(glob.glob(rundir+"rotation/*.fits.fz"), key=os.path.getmtime)

    with fits.open(files[0]) as f:
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
            if "Description" in line:  # For tile title: physical property
                name = line.strip().split("/")[-1].strip()
            # To see entire header, including comments starting with "/"
            # print line

        cmap = helix_tables(int(module.strip()), int(flag.strip()))

        scale = "[{0:.1f} Mpc]^2".format(float(scale)/1000)

        xlen = header['NAXIS1']
        ylen = header['NAXIS2']

    for filename in files:
        make_plot(filename, cmap, scale, xlen, ylen, name)


if __name__ == "__main__":
    timestamp = "20160707T0034"
    rundir = "../runs/{0}/".format(timestamp)

    if not (os.path.isdir(rundir+"/out") or os.path.exists(rundir+"/out")):
        os.mkdir(rundir+"/out")

    plot_all(rundir)
