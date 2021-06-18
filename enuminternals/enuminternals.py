import regina
import sys
from enuminternals.faultfinding import *
from enuminternals.sanity import *
from enuminternals.bead import *
from enuminternals.ne import *

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
        if hits.empty():
            # hLLAMkbeddfggghhbgahha is s441.
            # Simplifying harder shows this is true.
            mfld.simplifyExhaustive(2)
            hits = regina.Census.lookup(mfld)
        name_with_stuff = hits.first().name()
        name = name_with_stuff.split(' ')[0]
        names[name] = ne_sigs[ne_sig]
    return names

if __name__ == "__main__":
    print(enumerate_internal_necklaces(7, verbose=True))
