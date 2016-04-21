"""
File: helper_functions.py
Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
Version: 0.01 (Initial)
Date created: Fri Mar 04, 2016 03:46 pm
Last modified: Mon Mar 07, 2016 07:09 pm

Description: Helper functions for Cygnus A merger

"""

# From AMUSE code
def smart_length_units_for_vector_quantity(quantity):
    from amuse.units import units
    length_units = [units.Mpc, units.kpc, units.parsec, units.AU, units.RSun, units.km]
    total_size = max(quantity) - min(quantity)
    for length_unit in length_units:
        if total_size > (1 | length_unit):
            return length_unit
    return units.m


def print_progressbar(i, tot):
    import sys
    bar_width = 42  # obviously
    progress = float(i)/tot
    block = int(round(bar_width * progress))
    sys.stdout.write(
        "\r[{0}{1}] {2:.2f}% \t{3}/{4}"
        .format('#'*block, ' '*(bar_width - block),
                progress*100, i, tot))
    sys.stdout.flush()


def load_merger(timestamp):
    import pickle

    filename = "out/{0}/data/merger.dat".format(timestamp)
    print "Loading dumped ClusterMerger instance from", filename, "\n"
    return pickle.load(open(filename, "rb"))
