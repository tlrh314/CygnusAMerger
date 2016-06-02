import os
import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
# pyplot.rcParams.update({'font.size': 22})
# pyplot.rcParams.update({"text.usetex": True})
from astropy.io import fits

# http://www.ifweassume.com/2014/04/cubehelix-colormap-for-python.html
# https://github.com/jradavenport/cubehelix
import cubehelix

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


def make_video(analysis_dir):
    """ Make a video of four physical properties (options) """
    def plot_panels(chosen):
        """ Fill four panels with one snapshot, save png """
        projection = "_projection-z"

        headers = []
        data = []
        cmap = []
        names = []
        for option in chosen:
            with fits.open(analysis_dir+option+projection+".fits.fz") as f:
                headers.append(f[0].header)
                data.append(f[0].data)

                # Set colourtable based on what the fits header says it is
                for line in repr(f[0].header).split("\n"):
                    if "Effect_Module" in line:
                        module = line.strip().split("=")[-1].strip().split("/")[0]
                    if "Effect_Flag" in line:
                        flag = line.strip().split("=")[-1].strip().split("/")[0]
                    if "XYSize" in line:
                        scale = line.strip().split("=")[-1].strip().split("/")[0]
                    if "Description" in line:
                        names.append(line.strip().split("/")[-1])
                cmap.append(helix_tables(module, flag))
        scale = "[{0:.1f} Mpc]^2".format(float(scale)/1000)

        number_of_snapshots = headers[0]['NAXIS3']
        xlen = headers[0]['NAXIS1']
        ylen = number_of_pixels_y = headers[0]['NAXIS2']
        # See entire header, including comments starting with "/"

        pad = 4  # number of pixels padding for text placement
        for n in range(number_of_snapshots):
            # Set up four-panel plot, stitched together
            pyplot.figure(figsize=(16,16))
            gs1 = gridspec.GridSpec(2, 2)
            gs1.update(wspace=0, hspace=0) # remove spacing between axes.

            for i in range(4):  # gridspec indexes start at 0
                ax = pyplot.subplot(gs1[i])
                pyplot.axis('on')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_aspect('equal')

                # Plot every panel
                if i != 0:
                    ax.imshow(numpy.log(data[i][n]), cmap=cmap[i])
                else:
                    ax.imshow(data[i][n], cmap=cmap[i])
                pyplot.text(pad if i%2==0 else xlen-pad, pad,
                            names[i], color="white", size=22,
                            horizontalalignment="left" if i%2==0 else "right",
                            verticalalignment="top")
            pyplot.text(xlen-pad, ylen-pad, scale, color="white", size=22,
                        horizontalalignment="right", verticalalignment="bottom")

            # pyplot.suptitle("T = {0:05.2f} Myr".format(0.05*n),
            pyplot.suptitle("T = {0:04.2f} Myr".format(0.05*n),
                color="white", size=30)
            pyplot.tight_layout()
            pyplot.savefig("out/snapshot_{0:03d}.png".format(n))
            pyplot.close()

    options = {0: "dm-density",
               1: "dm-annihilation",
               2: "physical-density",
               3: "pressure",
               4: "sz-compton-y",
               5: "sz-thermal-dt-over-t",
               6: "temperature-emission-weighted",
               7: "temperature-spectroscopic",
               8: "xray-surface-brightness",
               }

    chosen = (options[0], options[7], options[8], options[4])
    number_of_snapshots = plot_panels(chosen)
    os.system("./make_movie.sh")


if __name__ == "__main__":
    make_video(analysis_dir="../runs/20160526T1354/analysis/")
