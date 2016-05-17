import sys

from amuse.units import units
from amuse.units import constants

m_p = constants.proton_mass.value_in(units.g)
mu = 0.17

def print_progressbar(i, N):
    # progressbar
    pbwidth = 42

    progress = float(i)/N
    block = int(round(pbwidth*progress))
    text = "\rProgress: [{0}] {1:.1f}%".format( "#"*block + "-"*(pbwidth-block), progress*100)
    sys.stdout.write(text)
    sys.stdout.flush()

    if i == (N-1):
        print " .. done"
