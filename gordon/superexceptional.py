# This module generates a list of candidate superexceptional manifolds,
# then determines which among these manifolds is in fact superexceptional.

import snappy
from .hyperbolicity import *

def find_superexceptional_fillings(name, verbose=False):
    supers = set()
    M = snappy.ManifoldHP(name)
    psf = {}
    if verbose:
        s = "fsx: Finding potential superexceptional fillings for {0}"
        print(s.format(name))
    find_potential_superexceptional_fillings(M, psf, verbose=verbose)
    for cusp in range(M.num_cusps()):
        for filling in psf[cusp]:
            mu = M.copy()
            mu.dehn_fill(filling, cusp)
            if verbose:
                print("fsx: Determining if {0} is superexceptional".format(mu))
            mu.canonize()
            if verbose:
                print("fsx: Canonized {0}".format(mu))
            N = snappy.ManifoldHP(mu.filled_triangulation().canonical_retriangulation(verified=True))
            if verbose:
                print("fsx: Filled {0}: now running is_superexceptional".format(mu))
            if is_superexceptional(N, verbose=verbose):
                supers.add((name,cusp,filling,N.identify()[0].name()))
    return supers

def find_potential_superexceptional_fillings(M, outdict, verbose=False):
    assert M.num_cusps() == 2
    cusps = range(2)

    # Our first order of business is constructing exc and hyp for cusps.
    # These are dicts of type \{0,1\}:\{(a,b)|a,b in N\}.
    # 0 and 1 representing the two cusps.
    cusp_hyp, cusp_exc = {}, {} 
    for tau in cusps:
        if verbose:
            print("fpsx_1: Classifying slopes on cusp {0}".format(tau))
        cusp_hyp[tau], cusp_exc[tau] = set(), set()
        classify_short_slopes(M, tau, cusp_hyp[tau], cusp_exc[tau], verbose=verbose)

    # Next we construct hyp and exc for fillings along short hyperbolic slopes.
    slope_hyp, slope_exc = {}, {}
    for tau in cusps:
        if verbose:
            print("fpsx_2: Constructing exc and hyp for cusp {0}".format(tau))
        for slope in cusp_hyp[tau]:
            if verbose:
                print("fpsx_2: Looking at slope {0}".format(slope))
            mfld = M.copy()
            mfld.dehn_fill(slope,tau)
            N = mfld.filled_triangulation()
            try:
                N.canonize()
            except:
                s = "fpsx_2: {0} {1} failed to canonize"
                raise Exception(s.format(M.name(), filling))
            if verbose:
                s = "fpsx_2: Classifying short slopes for {0}..."
                print(s.format((slope,tau)))
            slope_hyp[(slope,tau)], slope_exc[(slope,tau)] = set(), set()
            classify_short_slopes(N, 0,
                                  slope_hyp[(slope,tau)],
                                  slope_exc[(slope,tau)], verbose=verbose)

    # Next we verify that $exc(\tau)$ is always at most 8.
    if verbose:
        print("fpsx_3: Verifying exc(tau) is at most 8")
    for tau in cusps:
        if len(cusp_exc[tau]) > 8:
            s = "fpsx_3: {0} has more than 8 exceptional fillings on cusp {1}!"
            raise Exception(s.format(M.name(), tau))

    # Now we want all long slopes that might be superexceptional fillings.
    # For each cusp $\tau$, for each $s$ in $\tau$,
    # let $X(s) = \{s' \in hyp(\tau')\ |\ s \in exc(s')\}$.
    # Obviously we can't do this by iterating over all slopes $s$ in $\tau$.
    # Instead, we iterate over all $s'$ in $hyp(\tau')$,
    # and add $s'$ to each set $X(s)$ for all $s \in exc(s')$.
    long_maybe_superexc = {}
    for tau in cusps:
        if verbose:
            s = "Are there long possibly superexceptional fillings on cusp {0}?"
            print("fpsx_4: " + s.format(tau))
        long_maybe_superexc[tau] = set()
        X_dict = {}
        taup = 1 - tau
        if verbose:
            print("fpsx_4: Constructing X_dict")
        for sp in cusp_hyp[taup]:
            for s in slope_exc[(sp,taup)]:
                if not s in X_dict:
                    X_dict[s] = set()
                X_dict[s].add(sp)
        # Now we check whether or not equation (7--2) from Lemma 7.16 gives us the bound we want.
        if verbose:
            print("fpsx_4: Checking slopes in X_dict")
        for s in cusp_hyp[tau]:
            if verbose:
                print("fpsx_4: Checking slope {0}...".format(s))
            if not s in X_dict:
                X_dict[s] = set()
            rhs = len(cusp_exc[tau]) + len(X_dict[s])
            if rhs > 9:
                if verbose:
                    print("fpsx_4: Maybe superexceptional!")
                long_maybe_superexc[tau].add(s)
    for tau in cusps:
        outdict[tau] = set()
        outdict[tau].update(cusp_hyp[tau])
        outdict[tau].update(long_maybe_superexc[tau])

def classify_short_slopes(mfld, cusp, hyps, excs, verbose=False):
    """Classifies short slopes on cusp in mfld into hyp. and excs.
    mfld with whatever Dehn fillings it has attached should be verifiably hyperbolic.
    Otherwise short_slopes may raise an error or even dump core."""
    M = snappy.ManifoldHP(mfld)
    try:
        M.canonize()
    except:
        pass
    matrices = M.set_peripheral_curves('shortest', return_matrices=True)
    [(a,c),(b,d)] = matrices[cusp]
    if verbose:
        print("classify_short_slopes: enumerating short slopes on {0}".format(mfld))
    short_shorts = M.short_slopes(policy='greedy',
                                  method='maximal',
                                  verified=True,
                                  first_cusps=[cusp,])
    original_shorts = [(d*p-b*q, -c*p+a*q) for (p,q) in short_shorts[cusp]]
    for slope in original_shorts:
        if hyp_info(mfld, cusp, slope)[0]:
            hyps.add(slope)
        else:
            excs.add(slope)

def is_superexceptional(mfld,verbose=False):
    assert mfld.num_cusps() == 1
    assert mfld.is_orientable()
    if mfld.num_tetrahedra() > 9:
        print(f"more than 9 tets: {mfld}")
    hyps, excs = set(), set()
    classify_short_slopes(mfld, 0, hyps, excs, verbose=verbose)
    return len(excs) > 8

