import sys

def p2(a):
    return ((a)*(a))

def p3(a):
    return ((a)*(a)*(a))

def amuse_nth_root(quant, n):
    """ Simply telling AMUSE e.q. quant**(1./3) breaks the units :-( """
    return new_quantity((quant.number)**(1./n), (quant.unit ** (1./n)).to_simple_form())

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
