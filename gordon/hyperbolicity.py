import regina
import snappy
from enuminternals.faultfinding import *
from enuminternals.quick_checks import *
from enuminternals.ne import find_from

def hyp_regina(given_sig):
    M = regina.Triangulation3(given_sig)
    M.intelligentSimplify()
    # Essential vtx surfaces are only guaranteed for material triangulations.
    # We therefore truncating all ideal vertices, then simplify.
    M.idealToFinite()
    M.intelligentSimplify()
    material_sig = M.isoSig()
    # We want our surface labellings to be consistent.
    # To do so, we fix the triangulation's labelling.
    M = regina.Triangulation3(material_sig)
    
    # pi1 obstructions
    if has_common_axis_obstruction(M):
        s = "hyp_regina: {0} homeo. {1}: common axis"
        s = s.format(given_sig,material_sig)
        print(s)
        return (False, s.format(given_sig,material_sig))

    # Strict angle structures only exist for ideal triangulations.
    mu = regina.Triangulation3(given_sig)
    mu.finiteToIdeal()
    # N.B. finiteToIdeal does not necessarily yield an ideal triangulation.
    #      It may fail to crush all finite vertices into ideal vertices.
    #      In this case we will not get a strict angle structure.
    mu.intelligentSimplify()
    if mu.hasStrictAngleStructure():
        s = "hyp_regina: {0} homeo. {1}: strict angle structure"
        s = s.format(given_sig, mu.isoSig())
        print(s)
        return (True, s)

    # We use the material triangulation M for the surface enumeration.
    nsl = regina.NormalSurfaces
    qvs = nsl(M, regina.NS_QUAD, regina.NS_VERTEX)
    # The literature only guarantees existence results for fundamental surfaces.
    # Happily, for our data set, vertex surfaces alone turn out to suffice.
    
    idx = find_from(is_nonseparating_fault, qvs)
    if idx != None:
        s = "hyp_regina: {0} homeo. {1}: nonsep. fault at {2} in {1}"
        s = s.format(given_sig, material_sig, idx)
        print(s)
        return (False, s)

    idx = find_from(is_essential_sphere, qvs)
    if idx != None:
        s = "hyp_regina: {0} homeo. {1}: essential S2 at {2} in {1}"
        s = s.format(given_sig, material_sig, idx)
        print(s)
        return (False, s)

    # Solid tori were already caught by the search for a nonseparating fault.
    # if M.isSolidTorus():
    #    return (False, "Solid torus")

    idx = find_from(is_torus_fault, qvs)
    if idx != None:
        s = "hyp_regina: {0} homeo. {1}: torus fault at {2} in {1}"
        s = s.format(given_sig, material_sig, idx)
        print(s)
        return (False, s)

    idx = find_from(is_solid_torus_annulus, qvs)
    if idx != None:
        s = "hyp_regina: {0} homeo. {1}: solid torus A2 fault at {2} in {1}"
        s = s.format(given_sig, material_sig, idx)
        print(s)
        return (False, s)

    if M.hasBoundaryFacets():
        s = "hyp_regina: {0} homeo. {1}: faultless with nonempty boundary"
        s = s.format(given_sig,material_sig)
        print(s)
        return (True, s.format(given_sig,material_sig))

    if M.isHaken():
        s = "hyp_regina: {0} homeo. {1}: faultless and Haken"
        s = s.format(given_sig,material_sig)
        print(s)
        return (True, s.format(given_sig,material_sig))
    else:
        s = "hyp_regina: {0} homeo. {1}: faultless and not Haken!"
        s = s.format(given_sig,material_sig)
        print(s)
        return (None, s.format(given_sig,material_sig))

def hyp_snappy(snappy_mfld, verbose=False):
    M = snappy.Manifold(snappy_mfld)
    # As of 2019/11/21, one cannot use ManifoldHP in this step.
    # That causes a segfault for v1060(-2,1)(0,0)!
    # print "hyp_snappy: {0}".format(M)
    try:
        if M.verify_hyperbolicity()[0]:
            s = 'hyp_snappy: {0}: immediate'
            s = s.format(M.triangulation_isosig(decorated=False))
            print(s)
            return (True,s)
    except ZeroDivisionError:
        pass
    
    for cusp in range(M.num_cusps()):
        if M.cusp_info('is_complete')[cusp]:
            continue
        X = M.filled_triangulation([cusp,])
        if X.num_cusps() > 0:
            try:
                if X.verify_hyperbolicity()[0]:
                    s = 'hyp_snappy: {0}: after filling cusp {1} of {2}'
                    s = s.format(M.triangulation_isosig(decorated=False), cusp, M.num_cusps())
                    print(s)
                    return (True, s)
            except ZeroDivisionError:
                pass

    for curve in M.dual_curves():
        X = M.drill(curve)
        X.dehn_fill([(1,0)])
        try:
            X.canonize()
            assert X.verify_hyperbolicity()[0]
            s = 'hyp_snappy: {0}: after drilling along {1}'
            s = s.format(M.triangulation_isosig(decorated=False), curve)
            print(s)
            return (True, s)
        except:
            pass
    return (None, '')

def hyp_census(sig):
    # The following code works for Regina version 7.3.
    # The user will have to modify it appropriately for different versions.
    # It may be that census naming conventions could change.
    
    vs = regina.versionString()
    # String literals are automatically in Unicode in Python3.
    if not (vs == '7.3'): # As of 2023-07-31 sagedocker has version 7.3 installed.
        raise Exception("Unknown version of Regina")
    hits = regina.Census.lookup(sig)
    for hit in hits: 
        name = hit.name()
        print("hyp_census: {0} {1}".format(sig, name))
        if name[0:3] == 'Hyp':
            # A hyperbolic manifold from the closed orientable census.
            return (True, "hyp_census: {0}: {1}".format(sig,name))
        if name[0] in "0123456789" and name[1] == '.':
            # This is a manifold from the closed hyperbolic census.
            # Its name begins with its volume estimate.
            return (True, "hyp_census: {0}: {1}".format(sig,name))
        if name[0] in "msvt" or name[0:2] == "o9":
            # This is a manifold from the cusped orientable hyperbolic census.
            return (True, "hyp_census: {0}: {1}".format(sig,name))
        if name[0] in "0123456789" and name.split('a')[1][0] == 'h': 
            # This is a hyperbolic manifold from the prime knot census.
            return (True, "hyp_census: {0}: {1}".format(sig,name))
            
        if name[0:2] == 'L(':
            # A lens space.
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0:3] == 'SFS':
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0:5] == "T x I":
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0:2] == "S3":
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0:3] == "RP3":
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0:7] == "S2 x S1":
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0] in "0123456789" and name.split('a')[1][0] != 'h':
            # This is a nonhyperbolic manifold from the prime knot census.
            return (False, "hyp_census: {0}: {1}".format(sig,name))
        if name[0] == 'L' and name[1] in "123456789":
            # In the Christy knot and link census.
            # This census has both nonhyperbolic and hyperbolic elements.
            # For instance, L108019 is a Seifert-fibered space.
            continue 

        raise Exception("hyp_census: {0}: new name type: " + name)

    return (None, "hyp_census: " + sig)

def generate_sigs(snappy_mfld, fuel, verbose=False):
    # Regina does this more methodically.
    # But this function is not available in its Python bindings.
    # So we use the following jury-rigged alternative in SnapPy.
    L = fuel
    x = 0
    N = snappy_mfld.filled_triangulation()
    sig = N.triangulation_isosig(decorated=False)
    S = {sig}
    while L > 0:
        x = x+1
        for i in range(fuel):
            M = N.copy()
            for j in range(x):
                M.randomize()
            if M.num_cusps() > 0:
                # SnapPy 2.7 sends an unhandled abort signal
                # if one tries to fill a closed Triangulation.
                N = M.filled_triangulation()
            else:
                N = M
            N.simplify()
            this_sig = N.triangulation_isosig(decorated=False)
            if (this_sig not in S) and verbose:
                print("new sig at ({0}.{1}.{2})".format(L,x,i))
            S.add(this_sig)
        L = L // 2
    X = [(len(sig),sig) for sig in S]
    X.sort()
    return [sig for (n,sig) in X]

def hyp_info(mfld, cusp, slope, verbose=False):
    M = mfld.copy()
    M.dehn_fill(slope, cusp)
    sigs = generate_sigs(M, 4)
    # 32 seems to be enough to get > 1 choice for v1060 consistently for Manifolds.
    # For instances of ManifoldHP, 4 suffices.

    for sig in sigs:
        is_hyp = hyp_census(sig)
        if not is_hyp[0] == None:
            return is_hyp

    for sig in sigs:
        reg_N = regina.Triangulation3(sig)
        snp_N = snappy.ManifoldHP(reg_N.snapPea())
        is_hyp = hyp_snappy(snp_N, verbose)
        if not is_hyp[0] == None:
            return is_hyp

    sig = sigs[0]
    N = regina.Triangulation3(sig)
    x = hyp_regina(N)
    if x[0] != None:
        return x

    # Take care of jLLLAAQacfggghiiinkgokohqkv.
    N.simplifyExhaustive(2)
    x = hyp_census(N.isoSig())
    if x[0] != None:
        print("hyp_info: {0} simplifies to {1}".format(sig,N.isoSig()))
        return x
    
    raise Exception("Can't determine hyperbolicity of {0}".format(sig))

