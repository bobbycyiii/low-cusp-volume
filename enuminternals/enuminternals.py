import regina
import sys
from .faultfinding import *
from .sanity import *
from .bead import *
from .ne import *

def enumerate_internal_necklaces(bead, verbose=False):
    neckl_sigs = {}
    if verbose:
        print("We abbreviate enumerate_internal_necklaces to e_i_n.")
        print("e_i_n, {0} beads: Enumerating isosigs...".format(bead))
    match_sigs = enumerate_isosigs(bead)
    if verbose:
        print("e_i_n, {0} beads: Running basic sanity checks...".format(bead))
    for match_sig in match_sigs:
        mfld = regina.Triangulation3(match_sig)
        mfld.intelligentSimplify()
        sig = mfld.isoSig()
        if sig in neckl_sigs:
            continue
        if not is_nontrivial_link_exterior(mfld):
            continue
        # full necklace structures always have 2 or 3 cusps
        if mfld.countBoundaryComponents() not in (2, 3):
            continue
        neckl_sigs[sig] = match_sigs[match_sig]
    if verbose:
        print("e_i_n, {0} beads: Finding hyp. ancestral set...".format(bead))
    ne_sigs = {}
    for sig in neckl_sigs:
        mfld = regina.Triangulation3(sig)
        X = solve_problem_ne(mfld, verbose)
        for ne_sig in X:
            if not ne_sig in ne_sigs:
                if verbose:
                    print("e_i_n, {0} beads: Found new sig ".format(bead)
                          + ne_sig)
                ne_sigs[ne_sig] = neckl_sigs[sig]
    names = {}
    print("e_i_n, {0} beads: Recognizing signatures...".format(bead))
    for ne_sig in ne_sigs:
        mfld = regina.Triangulation3(ne_sig)
        hits = regina.Census.lookup(mfld)
        if len(hits) == 0:
            # hLLAMkbeddfggghhbgahha is s441.
            # Simplifying harder shows this is true.
            mfld.simplifyExhaustive(2)
            hits = regina.Census.lookup(mfld)
        name_with_stuff = hits[0].name()
        if verbose:
            print(name_with_stuff)
        name = name_with_stuff.split(' ')[0]
        # we have an L(3,1) compact piece which we eliminate by hand
        if name != "L(3,1)":
          names[name] = ne_sigs[ne_sig]
    return names

if __name__ == "__main__":
    neckls = set()
    for bead in [4,5,6,7]:
        neckls = neckls.union(enumerate_internal_necklaces(bead, verbose=True))
    print(neckls)
