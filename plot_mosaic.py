import aplpy
import numpy
import scipy
from scipy import ndimage

from deco import concurrent, synchronized

import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 28})
from astropy.io import fits


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

core = True
usm = False

obsdir = "../runs/ChandraObservation/"
# radio = obsdir+"radio5GHz.fits"
radio = obsdir+"StruisMosaics/radio/radio5GHz.fits"
mosaic = obsdir+"StruisMosaics/mosaic/cygnus_tot_flux.fits"

if usm:
    unsharp_mask(mosaic)
    import sys
    sys.exit(0)


gc = aplpy.FITSFigure(mosaic)

# aplpy.wcs_util(equinox="J2000")
# gc.pixel2world(BoxSize/XYPix)  for simulations?

if core:
    gc.show_colorscale(vmin=8.0e-10, vmax=2.0e-6,
                       stretch="log", cmap="spectral")
else:
    # gc.show_colorscale(vmin=9.0e-10, vmax=4.0e-8,
    #     stretch="log", cmap="spectral", smooth=9)
    gc.show_colorscale(vmin=7.0e-10, vmax=4.0e-8, stretch="log",
                       cmap="spectral", smooth=9)

gc.tick_labels.set_font(size=11, family="Times New Roman")
gc.axis_labels.set_font(size=14, family="Times New Roman")
gc.tick_labels.set_xformat("hh:mm:ss")
gc.tick_labels.set_yformat("dd:mm:ss")

# gc.add_colorbar()
# gc.colorbar.set_pad(0.1)

if core:
    from astropy import units
    from cosmology import CosmologyCalculator
    cc = CosmologyCalculator(0.0562)
    kpc2arcsec = 1/cc.kpc_DA
    gc.add_scalebar(30 * kpc2arcsec * units.arcsecond)
    gc.scalebar.set_corner("bottom right")
    gc.scalebar.set_label("30 kpc")
    gc.scalebar.set_color("black")
    gc.scalebar.set_font(size=24)
else:
    gc.add_scalebar(0.13227513)
    gc.scalebar.set_corner("bottom right")
    gc.scalebar.set_length(0.1)
    gc.scalebar.set_label("500 kpc")
    gc.scalebar.set_color("white")
    gc.scalebar.set_font(family="Times New Roman", size=24)

if core:
    gc.show_contour(radio, vmin=0.002, vmax=0.1, levels=15, smooth=1,
                    colors="black", lw=8)
    gc.recenter(299.86652, 40.734496,
                width=0.035694445, height=0.018194445)

gc.save("out/cygnus{0}_5GHz_xray.pdf".format("_core" if core else ""), dpi=300)
