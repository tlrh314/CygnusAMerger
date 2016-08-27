import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
pyplot.rcParams.update({'font.size': 22})
from scipy import special
import numpy

kpc2cm = 3.08568e+21
g2msun = 5.02785e-34

def p2(a):
    return ((a)*(a))
def p3(a):
    return ((a)*(a)*(a))

def analytical_mass(r, rho0, rc, rcut, beta):
    r2 = p2(r)
    rc2 = p2(rc)
    rcut2 = p2(rcut)
    sqrt2 = numpy.sqrt(2)

    A = (rc2 - rcut2)*(numpy.log(rcut2 - sqrt2*rcut*r + r2) \
        - numpy.log(rcut2 + sqrt2*rcut*r + r2))
    Bplus = 2*(rc2 + rcut2)*numpy.arctan(1 + sqrt2*r/rcut)
    Bmin = 2*(rc2 + rcut2)*numpy.arctan(1 - sqrt2*r/rcut)

    # NB the paper is slightly different: equation 2 does not contain 4*pi*rho
    M_gas_below_r = (4 * numpy.pi * rho0) *\
        rc2*p3(rcut)/(8*(p2(rcut2)+p2(rc2))) *\
        (sqrt2*(A - Bmin + Bplus) - 8*rc*rcut*numpy.arctan(r/rc))
    return M_gas_below_r

def analytical_density(r, rho0, rc, rcut, beta):
    rho_gas = rho0 * (1 + p2(r/rc))**(-1.5*beta)
    rho_gas /= (1 + p3(r/rcut) * (r/rcut))
    return 4*numpy.pi*p2(r)*rho_gas



r_list = []
rho_list = []
r_spline_list = []
rho_spline_list = []
Mr_numerical_list = []
Mr_twothirds_list = []
found_spline = False

with open("test.txt", "r") as f:
    for i, line in enumerate(f.readlines()):
        if "END OF TABLE" in line:
            found_spline = True
            continue
        columns = line.split(",")
        r = columns[0].split("=")[-1]
        rho = columns[1].split("=")[-1]
        if not found_spline:
            r_list.append(float(r))
            rho_list.append(float(rho))
        else:
            r_spline_list.append(float(r))
            rho_spline_list.append(float(rho))

pyplot.figure(figsize=(12, 9))
pyplot.loglog(r_list, rho_list, label="rho spline")
pyplot.loglog(r_spline_list, rho_spline_list, label="integrated spline mass")
pyplot.loglog(r_spline_list, analytical_mass(numpy.array(r_spline_list), 1, 25, 1200, 1/1.5), label="analytical 2/3 mass")
pyplot.loglog(r_spline_list, analytical_density(numpy.array(r_spline_list), 1, 25, 1200, 1/1.5), label="analytical 2/3 density")
pyplot.legend(loc="best")
pyplot.show()
