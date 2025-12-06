from enuminternals.bead import enumerate_isosigs
from enuminternals.quick_checks import is_nontrivial_link_exterior
from gordon.hyperbolicity import hyp_info
import regina
import snappy

def get_necklace_sigs(bead, verbose=True):
    neckl_sigs = {}
    if verbose:
        print("We abbreviate enumerate_internal_necklaces to e_i_n.")
        print("e_i_n, {0} beads: Enumerating isosigs...".format(bead))
    match_sigs = enumerate_isosigs(bead)
    if verbose:
        print("e_i_n, {0} beads: Getting isosigs for generalized link exteriors with 2 or 3 cusps...".format(bead))
    for match_sig in match_sigs:
        mfld = regina.Triangulation3(match_sig)
        # mfld.intelligentSimplify()
        sig = mfld.isoSig()
        if sig in neckl_sigs:
            continue
        if not is_nontrivial_link_exterior(mfld):
            continue
        # full necklace structures always have 2 or 3 cusps
        if mfld.countBoundaryComponents() not in (2, 3):
            continue
        neckl_sigs[sig] = match_sigs[match_sig]
    return neckl_sigs 

if __name__ == "__main__":
    n_s = [get_necklace_sigs(i) for i in [4,5,6,7]]
    num_neckl_sigs = sum(len(n_sigs) for n_sigs in n_s)
    print(f"There are {num_neckl_sigs} cusped necklace structures with 4 to 7 beads.\n")

    mfld_sigs = set() 
   
    for n_sigs in n_s:
        for sig in n_sigs:
            mfld = regina.Triangulation3(sig)
            mfld.intelligentSimplify()
            mfld_sigs.add(mfld.isoSig())
    print(f"These simplify to {len(mfld_sigs)} isomorphism signatures.")

    hyperbolic = set()
    common_axis = set()
    non_sep = set()
    other = set()

    for sig in mfld_sigs:
        M = snappy.ManifoldHP(sig)
        x = hyp_info(M, 0, (0,0))
        if x[0]:
            hyperbolic.add(sig)
        elif 'common axis' in x[1]:
            common_axis.add(sig)
        elif 'nonsep. fault' in x[1]:
            non_sep.add(sig)
        else:
            other.add(sig)

    print(f"There were {len(hyperbolic)} sigs: {hyperbolic}\n\n")
    print(f"There were {len(common_axis)} sigs: {common_axis}\n\n")
    print(f"There were {len(non_sep)} sigs: {non_sep}\n\n")
    print(f"There were {len(other)} sigs: {other}\n\n")
