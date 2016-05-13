"""
File: toyclusterparamsparser.py
Author: Timo L. R. Halbesma <timohalbesma@gmail.com>
Date created: Fri May 13, 2016 05:28 pm
Last modified: Fri May 13, 2016 06:00 pm


"""

def parse_toycluster_parms(filename):
    parameters = dict()

    with open(filename, "r") as f:
        for line in f:
            # Ignore commented lines
            if len(line) > 1 and not line.strip().startswith("%"):
                line = line.strip().split("%")[0]  # Ignore comments in lines
                keyvaluepair = line.split()
                if keyvaluepair[0] != "Output_file":
                    parameters[keyvaluepair[0]] = float(keyvaluepair[1])
                else:
                    parameters[keyvaluepair[0]] = keyvaluepair[1]

    return parameters


if __name__ == "__main__":
    parms = parse_toycluster_parms("toycluster.par")

    for key, value in parms.items():
        print key, "=", value
