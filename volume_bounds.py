# Uses Snappy 2.8 in Sage 9.1
from snappy import ManifoldHP
from snappy.verify.exceptions import ShapePositiveImaginaryPartNumericalVerifyError
from gordon.hyperbolicity import hyp_info
from regina import Triangulation3
from enuminternals.ne import solve_problem_ne
import sys

def length_upper_bound_Futer_Kalfagianni_Purcell(mfld, vol_upper_bound):
    V = mfld.volume(verified=True)
    RR = V.base_ring()
    pi = RR.pi()
    return 2*pi*1/(1-((RR(2)/RR(3)*(RR(vol_upper_bound).log()-V.log())).exp())).square_root()

fkp = length_upper_bound_Futer_Kalfagianni_Purcell

def fkp_slopes(name, cusp, volume_bound, verbose=False):
    parent = ManifoldHP(name)
    ell_bd = fkp(parent, volume_bound)
    all_short_cusp_slopes =             \
        parent.short_slopes(length=ell_bd,   \
                            policy='greedy', \
                            method='maximal',\
                            verified=True,   \
                            first_cusps=[cusp,])
    return all_short_cusp_slopes[cusp]
        
if __name__ == "__main__":
    verbose = True

    # One cusped manifolds with volume at most 3.07
    one_cusped_volume_bound = 3.07
    one_cusped_list = ['m003', 'm004', 'm006', 'm007', 'm009', 'm010', 'm011', \
                       'm015', 'm016', 'm017', 'm019', 'm022', 'm023', 'm026']

    if "one_cusp" in sys.argv: 
        # s776 has an isometry permuting the three cusps in a cycle.
        # So it will suffice to investigate one cusp, here cusp 0.
        print("Working on two-cusped fillings of s776...\n")
        magic_filling_names = set()
        for slope in fkp_slopes('s776', 0, one_cusped_volume_bound, verbose):
            magic = ManifoldHP('s776')
            magic.dehn_fill(slope, 0)
            N_snappy = magic.filled_triangulation()
            N_regina = Triangulation3(N_snappy)
            N_regina_sigs = solve_problem_ne(N_regina)
            reg_to_snap = lambda sig: ManifoldHP(Triangulation3(sig).snapPea())
            N_snappy_pieces = [reg_to_snap(sig) for sig in N_regina_sigs]
            for piece in N_snappy_pieces:
                if piece.num_cusps() < 2:
                    raise Exception("s776{0} has an interesting piece".format(slope))
                assert piece.num_cusps() == 2
                piece.canonize()
                names = [X.name() for X in piece.identify()]
                if names == []:
                    raise Exception("Didn't identify s776{0}".format(slope))
                magic_filling_names.add(names[0])

        print("\nThe two-cusped fillings of s776 of interest are as follows:")
        L = list(magic_filling_names); L.sort()
        print(L)

        reduced_necklace_names = {'s596', 's647', 's774', #'s776',
                                  's780', 's782', 's785', 'v2124',
                                  'v2355', 'v2533', 'v2644', 'v2731',
                                  'v3108', 'v3127', 'v3211', 'v3376'}

        two_cusped_parent_names = list(reduced_necklace_names.union(magic_filling_names))
        two_cusped_parent_names.sort()

        one_cusped_names = set()
        for name in two_cusped_parent_names:
            for cusp in [0,1]:
                for slope in fkp_slopes(name, cusp, one_cusped_volume_bound, verbose):
                    M = ManifoldHP(name)
                    (is_hyp, reason) = hyp_info(M, cusp, slope, verbose)
                    if verbose and not is_hyp:
                        print("{0}, cusp {1}, slope {2}: {3}".format(name, cusp, slope, reason))
                        continue
                    M.dehn_fill(slope, cusp)
                    N = M.filled_triangulation()
                    N.canonize()
                    if N.volume(verified=True) > one_cusped_volume_bound:
                        print("{0}, cusp {1}, slope {2}: volume > {3}".format(name, cusp, slope, one_cusped_volume_bound))
                        continue
                    print("{0}, cusp {1}, slope {2}: volume <~ {3} ".format(name, cusp, slope, one_cusped_volume_bound))
                    one_cusped_names.add(N.identify()[0].name())
                
        L = list(one_cusped_names); L.sort()
        if not L == one_cusped_list:
            raise Exception("Names differ from Table in paper!")
        print("\nOne-cusped manifolds with volume < 3.07")
        print(one_cusped_list)
        
    if "closed" in sys.argv: # not elif
    # Closed manifolds with volume at most 3.07/3.02
        drilling_constant = 3.02
        # From Lemma 3.1 of "Dehn surgery, homology and hyperbolic volume,"
        # by I. Agol, M. Culler, and P. Shalen,
        # at [AGT 6 (2006) 2297–2312](https://msp.org/agt/2006/6-5/p10.xhtml),
        # with [DOI link here](http://dx.doi.org/10.2140/agt.2006.6.2297)
        # and [arXiv link here](http://arxiv.org/abs/math.GT/0508208).
        #
        # If M is a closed orientable hyperbolic 3-manifold,
        # and M admits a (shortest) geodesic C with tuberad(C) >= log(3)/2,
        # then vol(M - C) < 3.02 * vol(M).
        #
        # A careful examination of this lemma's proof
        # and the statement and proof of Proposition 10.1 from
        # Agol-Storm-Thurston-Dunfield shows that the assumption
        # of C being a shortest geodesic can be relaxed to being an embedded geodesic.
        #
        # By Gabai & Trnkova, a closed orientable hyperbolic 3-manifold
        # either is Vol3, or admits an embedded geodesic C with tuberad(C) >= log(3)/2.
        #
        # So if vol(M) < 3.07/3.02, then either M is Vol3,
        # or M is a Dehn filling of a one-cusped manifold N with vol(N) <= 3.07,
        # i.e. one of the above one-cusped manifolds.
        
        closed_volume_bound = one_cusped_volume_bound/drilling_constant
        print("\n*****************")
        print("Determining closed orientable hyperbolic 3-manifolds")
        print("of volume at most {0}".format(closed_volume_bound))

        closed_fillings = []
        vol_failed = []
        for N_name in one_cusped_list:
            cusp = 0
            for slope in fkp_slopes(N_name, cusp, closed_volume_bound, verbose):
                N = ManifoldHP(N_name)
                (is_hyp, reason) = hyp_info(N, cusp, slope, verbose)
                if is_hyp == None:
                    raise Exception("{0}{1} has unknown hyperbolicity".format(N_name, slope))
                if not is_hyp:
                    if verbose:
                        print("{0}{1}: {2}".format(N_name, slope, reason))
                    continue
                N.dehn_fill(slope, cusp)
                try:
                    vol = N.volume(verified=True)
                    closed_fillings.append((vol, N_name, slope))
                    continue
                except ShapePositiveImaginaryPartNumericalVerifyError:
                    vol_failed.append((N.volume(), N_name, slope))

        closed_fillings.sort()
        print("\nClosed fillings with volume verifiably at most {0}:".format(closed_volume_bound))
        for datum in closed_fillings:
            if not (datum[0] > closed_volume_bound):
                print(datum)

        verified = [(N_name, N_slope) for (v, N_name, N_slope) in closed_fillings if v.upper() < closed_volume_bound]
        print(verified)
        assert verified == [('m003', (-3,1)), ('m003', (2,1)), ('m003', (-2,3)), ('m003', (-1,3)), ('m004', (-5,1)), ('m004', (5,1))]
        fmw = ManifoldHP('m003(-3,1)')
        assert fmw.is_isometric_to(ManifoldHP('m003(2,1)'))
        print("\nm003(-3,1) and m003(2,1) are isometric, with volume at most")
        print(fmw.volume(verified=True).upper())

        vol2 = ManifoldHP('m003(-2,3)')
        assert vol2.is_isometric_to(ManifoldHP('m003(-1,3)'))
        assert vol2.is_isometric_to(ManifoldHP('m004(-5,1)'))
        assert vol2.is_isometric_to(ManifoldHP('m004(5,1)'))
        print("\nm003(-2,3), m003(-1,3), m004(-5,1), and m004(5,1) are isometric,")
        print("with volume at most")
        print(vol2.volume(verified=True).upper())

        vol_failed.sort()
        still_failed = []
        print("\nClosed fillings with nonverified volume:")
        for (N_vol, N_name, N_slope) in vol_failed:
            N = ManifoldHP(N_name)
            N.dehn_fill(N_slope, 0)
            for (X_vol, X_name, X_slope) in closed_fillings:
                if abs(X_vol-N_vol) < 0.0001:
                    X = ManifoldHP(X_name)
                    X.dehn_fill(X_slope, 0)
                    try:
                        if X.is_isometric_to(N):
                            if verbose:
                                print("{0}{1} is isometric to {2}{3}".format(N_name, N_slope, X_name, X_slope))
                            break
                    except RuntimeError:
                        pass
            else:
                still_failed.append((N_vol, N_name, N_slope))
        still_failed.sort()
        print("\nClosed fillings still with unverified volume:")
        for (N_vol, N_name, N_slope) in still_failed:
            print("{0}{1} still failed".format(N_name, N_slope))

        failures = [('m007', (3,1)), ('m010', (-1,2)), ('m007', (3,2)), ('m006', (4,1)), ('m022', (-1,2))]
        assert [(N_name, N_slope) for (N_vol, N_name, N_slope) in still_failed] == failures

        print("\nm007(3,1) is Vol3 by definition.")
        Vol3 = ManifoldHP('m007(3,1)')

        assert ManifoldHP('m010(-1,2)').is_isometric_to(Vol3)
        print("\nm010(-1,2) is isometric to Vol3.")

        M = ManifoldHP('m036(3,-1)')
        assert ManifoldHP('m007(3,2)').is_isometric_to(M)
        assert M.verify_hyperbolicity()[0]
        v = M.volume(verified=True)
        print("\nm007(3,2) is isometric to m036(3,-1), with volume at least")
        print(v.lower())

        M = ManifoldHP('m027(3,1)')
        assert ManifoldHP('m006(4,1)').is_isometric_to(M)
        assert M.verify_hyperbolicity()[0]
        v = M.volume(verified=True)
        print("\nm006(4,1) is isometric to m027(3,1), with volume at least")
        print(v.lower())

        assert ManifoldHP('m022(-1,2)').is_isometric_to(M)
        print("\nm022(-1,2) coincidentally is also isometric to m027(3,1).")
        
                
