import aplpy
import numpy
import scipy
from scipy import ndimage

from deco import concurrent, synchronized

import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
from matplotlib import rc
matplotlib.rc("font", **{"size":16})

from astropy.io import fits


# Experimental
@concurrent(processes=4)
def plot_usm(original, sigma, amount, threshold):
    """ Unsharp mask: sharpened = original + (original - blurred) * amount """


    print "USM. Sigma =", sigma, "; amount =", amount,
    print "; threshold =", threshold

    blurred = scipy.ndimage.filters.gaussian_filter(
        original, order=0, sigma=sigma)

    toadd = ((original-blurred)*amount).clip(min=threshold)

    # pyplot.figure(figsize=(16, 16))
    # ax = pyplot.gca()
    # # order=0: convolve with Gaussian kernel; sigma = kernel stdev
    # # Imshow is "flipped" with respect to ds9 image. Origin="lower" fixed it.
    # ax.imshow(blurred, origin="lower")
    # ax.set_title("Blurred. Sigma = {0}".format(sigma))
    # pyplot.savefig("out/core_blurred_sigma_{0}.png"
    #                .format(sigma))
    # pyplot.close()

    # pyplot.figure(figsize=(16, 16))
    # ax = pyplot.gca()
    # ax.imshow(toadd, origin="lower")
    # ax.set_title(
    #     "Org-Blr. Sigma = {0:03d}. Amount = {1:1.03f}. Threshold = {2:1.03f}"
    #         .format(sigma, amount, threshold))
    # pyplot.savefig(
    #     "out/core_added_sigma_{0:03d}_amount_{1:3.03f}_{2:3.03f}.png"
    #         .format(sigma, amount, threshold))
    # pyplot.close()

    # pyplot.figure(figsize=(16, 16))
    # ax = pyplot.gca()
    # ax.imshow(original+toadd, origin="lower")
    # ax.set_title(
    #     "USM. Sigma = {0:03d}. Amount = {1:1.03f}. Threshold = {2:1.03f}"
    #     .format(sigma, amount, threshold))
    # pyplot.savefig(
    #     "out/core_USM_sigma_{0:03d}_amount_{1:3.03f}_{2:3.03f}.png"
    #         .format(sigma, amount, threshold))
    # pyplot.close()

    pyplot.figure(figsize=(16, 16))
    ax = pyplot.gca()
    ax.imshow(original-blurred+toadd, origin="lower")
    ax.set_title(
        "USM. Sigma = {0:03d}. Amount = {1:1.03f}. Threshold = {2:1.03f}"
        .format(sigma, amount, threshold))
    pyplot.savefig(
        "out/core_experiment_sigma_{0:03d}_amount_{1:3.03f}_{2:3.03f}.png"
            .format(sigma, amount, threshold))
    pyplot.close()


    return


# Experimental
@synchronized
def unsharp_mask(mosaic):
    with fits.open(mosaic) as f:
        header = f[0].header
        for line in repr(header).split("\n"):
            print line
        data = f[0].data

        xlen = int(header["NAXIS1"])
        ylen = int(header["NAXIS2"])

        xcenter = int(header["CRPIX1"])
        ycenter = int(header["CRPIX2"])

        original = numpy.log10(data[xcenter-150:xcenter+250,
            ycenter-350:ycenter+50].clip(min=8e-10, max=2.0e-6))

    """
    The parallel library seems to horibly break when I first plot in the
    synchronized function. The error thrown is:

    The process has forked and you cannot use this CoreFoundation functionality safely. You MUST exec(). Break on __THE_PROCESS_HAS_FORKED_AND_YOU_CANNOT_USE_THIS_COREFOUNDATION_FUNCTIONALITY___YOU_MUST_EXEC__() to debug.
    """

    # Imshow is "flipped" with respect to ds9 image. Origin="lower" fixed it.
    # print "Showing original"
    # pyplot.figure()
    # pyplot.imshow(original, origin="lower")
    # pyplot.savefig("out/core.pdf", dpi=300)
    # pyplot.close()

    for sigma in numpy.arange(30, 110, 10):
        for amount in numpy.arange(0.05, 0.25, 0.05):
            for threshold in range(11):
                plot_usm(original, sigma, amount, 0.005*threshold)

    # This really does not seem to work well
    # for sigma in numpy.arange(0.5, 2, 0.5):
    #     for amount in numpy.arange(0.5, 1.5, 0.5):
    #         plot_usm(original, sigma, amount)


def plot_zoomin_of_core(mosaic, radio, cygA):
    """ Zoomin of the AGN feedback core region where FR-II interacts with gas
        @param mosaic: path to the Chandra x-ray mosaic fits file
        @param cygA  : tuple with RA, dec of CygA centroid
    """
    gc = aplpy.FITSFigure(mosaic)

    # Show xray observation with stretch to highlight the core
    gc.show_colorscale(vmin=8.0e-10, vmax=2.0e-6,
                       stretch="log", cmap="spectral")

    # Add the famous 5GHz radio contours
    gc.show_contour(radio, vmin=0.002, vmax=0.1, levels=15, smooth=1,
                    colors="black", lw=8)

    # Show a scale bar of 30 kpc after unit conversions
    from astropy import units
    from cosmology import CosmologyCalculator
    cc = CosmologyCalculator(0.0562)
    kpc2arcsec = 1/cc.kpc_DA
    gc.add_scalebar(30 * kpc2arcsec * units.arcsecond)
    gc.scalebar.set_corner("bottom right")
    gc.scalebar.set_linewidth(4)
    gc.scalebar.set_label("30 kpc")
    gc.scalebar.set_color("black")

    # Zoom in on the central region
    gc.recenter(cygA[0], cygA[1], width=0.037, height=0.018)
    pyplot.gca().tick_params(axis="both", which="both", colors="k", reset=True)

    # Pretty notation on the axes
    gc.tick_labels.set_xformat("hh:mm:ss")
    gc.tick_labels.set_yformat("dd:mm:ss")

    ax = pyplot.gca()
    ax.tick_params(axis="both", which="minor", colors="k",
                   pad=8, width=2, size=4, reset=True)
    ax.tick_params(axis="both", which="major", colors="k",
                   pad=8, width=2, size=8, reset=True)

    # gc.add_colorbar()
    # gc.colorbar.set_pad(0.1)

    pyplot.tight_layout()
    gc.save("out/CygA_Radio_5GHz.pdf", dpi=300)


def plot_mosaic_with_ruler(mosaic, cygA, cygB):
    gc = aplpy.FITSFigure(mosaic)

    # Add smoothed log-stretch of the entire mosaic
    gc.show_colorscale(vmin=7.0e-10, vmax=1.0e-6, stretch="log",
                       cmap="spectral", smooth=9)

    # Add scale. Length is 500 kpc after unit conversions
    gc.add_scalebar(0.13227513)
    gc.scalebar.set_corner("bottom right")
    gc.scalebar.set_length(0.1)
    gc.scalebar.set_linewidth(4)
    gc.scalebar.set_label("500 kpc")
    gc.scalebar.set_color("white")

    # Find the pixels of the centroids
    cygA_x, cygA_y = gc.world2pixel(cygA[0], cygA[1])
    cygB_x, cygB_y = gc.world2pixel(cygB[0], cygB[1])

    ax = pyplot.gca()
    ax.plot([cygA_x, cygB_x], [cygA_y, cygB_y], c="w", lw=1)

    # Eyeballed coordinates in ds9 :-) ...
    text_x, text_y = gc.world2pixel( 299.78952, 40.816273 )
    ax.text(text_x, text_y, '700.621"', ha="center", va="center", color="white",
            rotation=51, weight="bold", fontsize=22)

    # Pretty notation on the axes
    gc.tick_labels.set_xformat("hh:mm:ss")
    gc.tick_labels.set_yformat("dd:mm:ss")

    ax.tick_params(axis="both", which="minor", colors="k",
                   pad=8, width=2, size=4, reset=True)
    ax.tick_params(axis="both", which="major", colors="k",
                   pad=8, width=2, size=8, reset=True)

    # Zoom in a bit more on the merger region
    gc.recenter(299.78952, 40.81, width=0.185, height=0.185)
    pyplot.tight_layout()
    pyplot.show()
    gc.save("out/mosaic_xray_ruler.pdf", dpi=300)


def plot_mosaic_with_wedges(mosaic, cygA):
    """ Plot the merger, hot, cold regions in smoothed mosaic
        @param mosaic: path to the Chandra x-ray mosaic fits file
        @param cygA  : tuple with RA, dec of CygA centroid
    """
    gc = aplpy.FITSFigure(mosaic)

    # Add smoothed log-stretch of the entire mosaic
    gc.show_colorscale(vmin=7.0e-10, vmax=4.0e-8, stretch="log",
                       cmap="spectral", smooth=9)

    # Add scale. Length is 500 kpc after unit conversions
    gc.add_scalebar(0.13227513)
    gc.scalebar.set_corner("bottom right")
    gc.scalebar.set_length(0.1)
    gc.scalebar.set_linewidth(4)
    gc.scalebar.set_label("500 kpc")
    gc.scalebar.set_color("white")

    # Find the pixels of the centroids
    x_pix, y_pix = gc.world2pixel(cygA[0], cygA[1])

    # Cut-out angles: 6, 96, 225 and 315 degrees.
    radii = numpy.linspace(0, 4500, 100)
    x6 = numpy.zeros(len(radii))
    x96 = numpy.zeros(len(radii))
    x225 = numpy.zeros(len(radii))
    x315 = numpy.zeros(len(radii))
    y6 = numpy.zeros(len(radii))
    y96 = numpy.zeros(len(radii))
    y225 = numpy.zeros(len(radii))
    y315 = numpy.zeros(len(radii))
    for i, r in enumerate(radii):
        x6[i] = r*numpy.cos(9*numpy.pi/180)
        y6[i] = r*numpy.sin(6*numpy.pi/180)

        x96[i] = r*numpy.cos(96*numpy.pi/180)
        y96[i] = r*numpy.sin(96*numpy.pi/180)

        x225[i] = r*numpy.cos(225*numpy.pi/180)
        y225[i] = r*numpy.sin(225*numpy.pi/180)

        x315[i] = r*numpy.cos(315*numpy.pi/180)
        y315[i] = r*numpy.sin(315*numpy.pi/180)

    ax = pyplot.gca()
    ax.plot(x6+x_pix, y6+y_pix, c="w", lw=2)
    ax.plot(x96+x_pix, y96+y_pix, c="w", lw=2)
    ax.plot(x225+x_pix, y225+y_pix, c="w", lw=2)
    ax.plot(x315+x_pix, y315+y_pix, c="w", lw=2)
    ax.text(0.65, 0.85, "MERGER", ha="center", va="center", color="white",
            bbox=dict(facecolor="green", edgecolor="green", pad=6),
            weight="bold", transform=pyplot.gca().transAxes, fontsize=16)
    ax.text(0.1, 0.5, "HOT", ha="center", va="center", color="white",
            bbox=dict(facecolor="red", edgecolor="red", pad=6),
            weight="bold", transform=pyplot.gca().transAxes, fontsize=16)
    ax.text(0.65, 0.4, "HOT", ha="center", va="center", color="white",
            bbox=dict(facecolor="red", edgecolor="red", pad=6),
            weight="bold", transform=pyplot.gca().transAxes, fontsize=16)
    ax.text(0.3, 0.05, "COLD", ha="center", va="center", color="white",
            bbox=dict(facecolor="purple", edgecolor="purple", pad=6),
            weight="bold", transform=pyplot.gca().transAxes, fontsize=16)

    # Pretty notation on the axes
    gc.tick_labels.set_xformat("hh:mm:ss")
    gc.tick_labels.set_yformat("dd:mm:ss")

    ax.tick_params(axis="both", which="minor", colors="k",
                   pad=8, width=2, size=4, reset=True)
    ax.tick_params(axis="both", which="major", colors="k",
                   pad=8, width=2, size=8, reset=True)

    pyplot.tight_layout()
    gc.save("out/mosaic_xray_wedges.pdf", dpi=300)


if __name__ == "__main__":
    # Coordinates of the CygA centroid
    cygA = ( 299.8669, 40.734496 )
    cygB = ( 299.7055, 40.884849 )

    # Data directory and radio/xray observation fits files
    obsdir = "../runs/ChandraObservation/"
    radio = obsdir+"StruisMosaics/radio/radio5GHz.fits"
    xray = obsdir+"StruisMosaics/mosaic/cygnus_tot_flux.fits"

    # plot_zoomin_of_core(xray, radio, cygA)
    # plot_mosaic_with_wedges(xray, cygA)
    plot_mosaic_with_ruler(xray, cygA, cygB)

    # Experimental
    # unsharp_mask(mosaic)
