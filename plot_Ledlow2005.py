import numpy
import scipy
from scipy.stats import kde
import astropy
import astropy.units as u
from matplotlib import pyplot
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = None

from plotsettings import PlotSettings
style = PlotSettings()


def plot_velocity_histogram(t):
    cygA = t[49]  # eyeballed

    binwidth = 500  # km /s
    xmin, xmax = 14000, 22000  # km/s

    pyplot.figure(figsize=(12, 12))
    pyplot.hist(t["HV"], bins=range(xmin, xmax+binwidth, binwidth),
                facecolor="none", edgecolor="k", histtype="step")
    pyplot.xlim(xmin-binwidth, xmax+binwidth)
    pyplot.xlabel(r"$V_H$ (km s$^{-1}$)")
    pyplot.ylabel("Number of Galaxies")
    pyplot.tight_layout()


def plot_galaxy_positions(t):
    ra = astropy.coordinates.Angle(t['RAJ2000'].filled(numpy.nan), unit=u.degree)
    dec = astropy.coordinates.Angle(t['DEJ2000'].filled(numpy.nan), unit=u.degree)

    fig = pyplot.figure(figsize=(12, 12))
    pyplot.gca().invert_xaxis()

    # Filled symbols show velocities within 1 sigma of the mean
    # Squares indicate velocity less than the mean, triangles greater

    # caption says biweight mean, but median seems to give matching results

    mean = numpy.mean(t["HV"])
    median = numpy.median(t["HV"])
    print mean
    print median
    for ra_i, dec_i, v_i, ve_i in zip(ra, dec, t["HV"], t["e_HV"]):
        if v_i < median:
            if v_i-median < ve_i:
                pyplot.plot(ra_i, dec_i, "s", mec="k", mfc="k", ms=8)
            else:
                pyplot.plot(ra_i, dec_i, "s", mec="k", mfc="none", ms=8)
        else:
            if v_i+mean > ve_i:
                pyplot.plot(ra_i, dec_i, "^", mec="k", mfc="k", ms=8)
            else:
                pyplot.plot(ra_i, dec_i, "^", mec="k", mfc="none", ms=8)
    pyplot.gca().ticklabel_format(useOffset=False)# no offset
    pyplot.xlim(20.025, 19.950)
    pyplot.xlabel("RA (hours)")
    pyplot.ylim(40.5, 41.2)
    pyplot.ylabel("DEC (degrees)")
    # pyplot.tight_layout()
    # pyplot.show()

def plot_galaxy_contours(t):
    ra = astropy.coordinates.Angle(t['RAJ2000'].filled(numpy.nan), unit=u.degree)
    dec = astropy.coordinates.Angle(t['DEJ2000'].filled(numpy.nan), unit=u.degree)

    # https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.stats.gaussian_kde.html
    # Perform kernel density estimate
    xmin, xmax = ra.min().value, ra.max().value
    ymin, ymax = dec.min().value, dec.max().value
    X, Y = numpy.mgrid[xmax:xmin:100j, ymin:ymax:100j]
    positions = numpy.vstack([X.ravel(), Y.ravel()])
    values = numpy.vstack([ra, dec])
    kernel = kde.gaussian_kde(values)
    Z = numpy.reshape(kernel(positions).T, X.shape)
    Z /= ((41.1-40.6)*60*(20.025-19.50)*60)

    print numpy.max(Z)

    # fig = pyplot.figure(figsize=(12, 12))

    #pyplot.imshow(numpy.rot90(Z), cmap=pyplot.cm.gist_earth_r,
    #              extent=[xmax, xmin, ymin, ymax])

    cset = pyplot.contour(X, Y, Z, colors="k",
            levels=numpy.array([1,2,3,5,8,10])*Z.max()/10)
    #cset = pyplot.contour(X, Y, Z, colors="k",
    #        levels=numpy.linspace(Z.min(), Z.max(), 15))


    pyplot.clabel(cset, inline=1, fontsize=10)
    pyplot.xscale("log")
    #pyplot.gca().ticklabel_format(useOffset=False)# no offset
    # pyplot.xlim(20.025, 19.950)
    pyplot.xlabel("RA (hours)")
    # pyplot.ylim(40.5, 41.2)
    pyplot.ylabel("DEC (degrees)")
    pyplot.xticks([20.02, 20, 19.98, 19.96], ["20.02", "20", "19.98", "19.96"])
    # pyplot.tight_layout()
    pyplot.show()




if __name__ == "__main__":
    # print Vizier.query_object("Cygnus A")
    # print Vizier.find_catalogs("Ledlow Cygnus").keys()
    ledlow2005 = Vizier.get_catalogs("J/AJ/130/47")[0]  # table1
    ledlow2005_not = Vizier.get_catalogs("J/AJ/130/47")[1]  # table1
    # plot_galaxy_positions(ledlow2005_not)
    # plot_galaxy_contours(ledlow2005_not)

    # plot_velocity_histogram(ledlow2005)
    plot_galaxy_positions(ledlow2005)
    plot_galaxy_contours(ledlow2005)

    pyplot.show()
