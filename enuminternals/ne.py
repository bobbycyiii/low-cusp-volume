import regina
from .faultfinding import *
from .sanity import is_nontrivial_link_exterior, is_closed_oriented, has_common_axis_obstruction

def find_from(predicate, F):
    n = F.size()
    for i in range(n):
        surf = F.surface(i)
        if predicate(surf):
            return i
    else:
        return None

def solve_problem_ne(mfld, verbose=False):
    """Returns a set containing sigs whose manifolds constitute
    an ancestral set for the finite-volume hyperbolic manifolds into
    which mfld embeds nonelementarily."""
    assert is_nontrivial_link_exterior(mfld) or is_closed_oriented(mfld)
    sig = mfld.isoSig()
    if has_common_axis_obstruction(mfld):
        # No nonelementary embeddings.
        if verbose:
            print("ne: {0}: common axis obstruction".format(sig))
        return set()
    mu = regina.Triangulation3(sig)
    mu.finiteToIdeal()
    mu.intelligentSimplify()
    if mu.hasStrictAngleStructure():
        if verbose:
            print("ne: {0}: strict angle structure".format(sig))
        return set([mu.isoSig()])
    # Regina's simplification algorithm is
    # randomized, and does not always return the
    # same isomorphism class of triangulation.
    # In hopes of ensuring a more reliably
    # consistent computation, we appeal to
    # this procedure several times before
    # selecting a triangulation to work with.
    #
    # We found it most economical to place
    # this simplification here, rather than at the
    # end of every cutting operation.
    musigs = []
    for i in range(20):
        mu = regina.Triangulation3(sig)
        mu.idealToFinite()
        mu.intelligentSimplify()
        musigs.append((mu.countTetrahedra(),mu.isoSig()))
    musigs.sort()
    M = regina.Triangulation3(musigs[0][1])
    material_sig = M.isoSig()
    nsl = regina.NormalSurfaces.enumerate
    F = nsl(M, regina.NS_QUAD, regina.NS_VERTEX)
    idx = find_from(is_nonseparating_closed_fault, F)
    if idx != None:
        if verbose:
            s = "ne: {0} homeo. {1}: nonsep. at index {2} in {1}"
            print(s.format(sig, material_sig,idx))
        return set()
    if verbose:
        print("no nonsep in {}".format(material_sig))
    idx = find_from(is_essential_sphere, F)
    if idx != None:
        S = F.surface(idx)
        # The appeal to signatures ensures there are
        # no references to objects deleted after the call to unSum.
        Lsig,Rsig = unsum(S)
        L = regina.Triangulation3(Lsig)
        R = regina.Triangulation3(Rsig)
        if verbose:
            s = "ne: {0} homeo. {1}: reducing S2 at index {2} in {1}"
            print(s.format(sig, material_sig, idx))
        return solve_problem_ne(L).union(solve_problem_ne(R))
    if verbose:
        print("{} is irreducible".format(material_sig))
    idx = find_from(is_torus_fault, F)
    if idx != None:
        S = F.surface(idx)
        LR = S.cutAlong()
        LR.intelligentSimplify()
        L, R = LR.triangulateComponents()
        if verbose:
            s = "ne: {0} homeo. {1}: essential T2 at index {2} in {1}"
            print(s.format(sig, material_sig, idx))
        return solve_problem_ne(L).union(solve_problem_ne(R))
    if verbose:
        print("{} is atoroidal".format(material_sig))
    idx = find_from(is_solid_torus_annulus, F)
    if idx != None:
        if verbose:
            s = "ne: {0} homeo. {1}: s.t.s. A2 at index {2} in {1}"
            print(s.format(sig, material_sig, idx))
        return set()
    if verbose:
        print("{} is not two solid tori".format(material_sig))
    # If we are here, then the manifold is either faultness or
    # only contains non-separating faults with boundary.
    # In particular the manifold must be either
    # faultless or Seifert fibred
    idx = find_from(is_nonseparating_nonclosed_fault, F)
    if idx != None:
        if verbose:
            s = "ne: {0} homeo. {1}: only nonsep. nonclosed, e.g. index {2} in {1}"
            print(s.format(sig, material_sig,idx))
        return set()
    # Manifold must be faultless
    if verbose:
        print("ne: {0} homeo. {1}: faultless".format(sig, material_sig))
    return set([mfld.isoSig()])

