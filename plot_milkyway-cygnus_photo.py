import aplpy
from matplotlib import pyplot

def plot_optical(optical, cygA):
    # Started with photo of E. Hanko (png). Split up RGB in separate channels.
    # Generated fits for each channel using nova.astrometry.net. Now combining
    # the fits with WCS to overplot Xray
    #aplpy.make_rgb_cube(["../runs/MWCyg_red.fits", "../runs/MWCyg_green.fits",
    #                 "../runs/MWCyg_blue.fits"], "out/MWCyg_rgb.fits", north=True)
    #aplpy.make_rgb_image("out/MWCyg_rgb.fits", "out/MWCyg_rgb.png")
    gc = aplpy.FITSFigure("out/MWCyg_rgb_2d.fits")
    gc.show_rgb("out/MWCyg_rgb.png")

    # Find the pixels of the centroids
    cygA_x, cygA_y = gc.world2pixel(cygA[0], cygA[1])
    print cygA_x, cygA_y
    ax = pyplot.gca()
    ax.plot(cygA_x, cygA_y, c="r", lw=10)

    # Pretty notation on the axes
    gc.tick_labels.set_xformat("hh")
    gc.tick_labels.set_yformat("dd")

    ax.tick_params(axis="both", which="minor", colors="k",
                   pad=8, width=2, size=4, reset=True)
    ax.tick_params(axis="both", which="major", colors="k",
                   pad=8, width=2, size=8, reset=True)

    # Zoom in a bit more on the merger region
    gc.recenter(299.78952, 40.81, radius=30)
    # pyplot.tight_layout()
    pyplot.show()
    # gc.save("out/optical.pdf", dpi=300)
    return

    from astropy.io import fits
    with fits.open(optical) as f:
        print f.info()
        header = f[0].header
        data = f[0].data
    for line in repr(header).split("\n"):
        print line

    print data.shape
    red = data[0]
    green = data[1]
    blue = data[2]
    rgb = numpy.dstack((red, green, blue))

    # titles = ['R', 'G', 'B', 'RGB']
    # cmaps = [pyplot.cm.Reds_r, pyplot.cm.Greens_r, pyplot.cm.Blues_r, None]

    # fig, axes = pyplot.subplots(1, 4, figsize=(13,3))
    # objs = zip(axes, [red, green, blue, rgb], titles, cmaps)

    # for ax, channel, title, cmap in objs:
    #     ax.imshow(channel, cmap=cmap)
    #     ax.set_title(title)
    #     # ax.set_xticks(())
    #     # ax.set_yticks(())

    pyplot.figure()
    pyplot.imshow(rgb)
    pyplot.show()

    return


if __name__ == "__main__":
    # Coordinates of the CygA centroid
    cygA = ( 299.8669, 40.734496 )
    cygB = ( 299.7055, 40.884849 )

    # Data directory and radio/xray/optical observation fits files
    obsdir = "../runs/ChandraObservation/"
    optical = "../runs/MWCyg.fits"

    plot_optical(optical, cygA)

