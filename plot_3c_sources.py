import numpy
import astropy.units as u
from matplotlib import pyplot
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = None

from cosmology import CosmologyCalculator
# from macro import print_progressbar
# from plotsettings import PlotSettings
# style = PlotSettings()

def plot_3cr(t):
    """ Reproduction of Stockton & Ridgway 1996 """
    jy = t["S_178_"]  # already in units of Janksy b/c Vizier astroquery
    # Convert Janksy to SI equivalent
    power = jy.to(u.W / u.m**2 / u.Hz,
                  equivalencies=u.spectral_density(178 * u.MHz))
    # Stockton & Ridgway assume H0 = 50 km/s/Mpc, we take 70
    distance = [CosmologyCalculator(z, 70).DL_Mpc for z in t["z"]]
    t["P"] = power  # already has units
    t["d"] = distance * u.Mpc

    # Sketchy factor pi. 4*pi would be understandable, but now plot matches
    # when H_0 = 50
    p178 = numpy.pi*t["d"].to(u.m)**2 * t["P"]

    pyplot.figure(figsize=(12,12))
    pyplot.scatter(t["z"], p178, marker="+", s=50, c="k")
    for i in ["405.0", "348.0", "123.0", "20.0", "427.1",
                "295.0", "265.0", "237.0", "268.1", "280.0"]:
        index = numpy.where(t["_3CR"] == i)[0][0]
        pyplot.text(t[index]["z"], p178[index].value,
                    i[:-2] if i[-1]=="0" else i, fontsize=22)
    pyplot.xlim(-0.05, 1.1)
    pyplot.ylim(-1.5e27, 3.2e28)
    pyplot.xlabel(r"$z$")
    pyplot.ylabel(r"$P_{178}$ (W Hz$^{-1}$)")
    pyplot.savefig("out/3CR.pdf", dpi=300)


if __name__ == "__main__":
    """ References
        3C : Edge+ 1959-1962 -- 1959MmRAS..68...37E
             Benett 1962 -- 1962MmRAS..68..163B  (VizieR)
        3CR: Spinrad+ 1985 -- 1985PASP...97..932S
                           -- 1986PASP...98..628S (Errratum)
    """
    c_3cr = Vizier.get_catalogs("J/PASP/97/932/3cr")[0]
    plot_3cr(c_3cr)
