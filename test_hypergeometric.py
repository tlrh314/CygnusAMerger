import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({'font.size': 22})
from scipy import special
import numpy

from macro import *
import convert

kpc2cm = 3.08568e+21
g2msun = 5.02785e-34


r_list = []
Mr_list = []
with open("gsl_hyper.txt", "r") as f:
    for i, line in enumerate(f.readlines()):
        columns = line.split()
        r = columns[0].split("=")[-1]
        Mr = columns[-1].split("=")[-1]
        r_list.append(r)
        Mr_list.append(float(Mr))

print Mr_list[0:20]
print Mr_list[-19:]
pyplot.figure(figsize=(16,12))
pyplot.plot(r_list, Mr_list, label="GSL with transformations")

r = numpy.arange(1, 1000, 1) * kpc2cm
#r = 1346.84 * kpc2cm
rho0 = convert.ne_to_rho(5.895e-2)
#rho0 = 1
rc = 26.141 * kpc2cm
beta = 0.53869
# beta = 2.0/3
# print -1.0*p2(r/rc)
M =  special.hyp2f1(1.5, 1.5*beta, 2.5, -p2(r/rc))
# print M
M *= rho0/3.0 * 4 * numpy.pi * p3(r)
# print M
#import sys; sys.exit(0)
M_analytical =  4*numpy.pi*p3(rc)*rho0 * (r/rc - numpy.arctan(r/rc))
pyplot.loglog(r/kpc2cm, M*g2msun, label="scipy.special.hyp2f1, beta={0}".format(beta))
pyplot.loglog(r/kpc2cm, M_analytical*g2msun, label="Analytical (beta 2/3)")
pyplot.legend(loc=4)
pyplot.xlabel(r"$r$ [kpc]")
pyplot.ylabel(r"$M(<r)$[MSun]")
pyplot.show()
