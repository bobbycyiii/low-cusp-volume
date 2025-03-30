import regina

def is_nontrivial_link_exterior(mfld):
    M = regina.Triangulation3(mfld)
    M.idealToFinite()
    M.intelligentSimplify()
    if not M.isConnected():
        return False
    if not M.isOrientable():
        return False
    if M.countBoundaryComponents() == 0:
        return False
    for cpt in M.boundaryComponents():
        t = cpt.build()
        if not t.eulerChar() == 0:
            return False
    else:
        return True

def is_closed_oriented(mfld):
    M = regina.Triangulation3(mfld)
    M.idealToFinite()
    M.intelligentSimplify()
    if not M.isConnected():
        return False
    if not M.isOrientable():
        return False
    if M.countBoundaryComponents() > 0:
        return False
    return True

def is_common_axis_commutator(expr):
    l = expr.terms()
    if (len(l) != 4
        or l[0] != l[2].inverse()
        or l[1] != l[3].inverse()):
        return False
    else:
        return True

def is_common_axis_equation(expr):
    l = expr.terms()
    return len(l) == 2

def has_common_axis_obstruction(mfld):
    G = mfld.fundamentalGroup()
    G.intelligentSimplify()
    n = G.countGenerators()
    CA = {i:set() for i in range(n)}
    for r in G.relations():
        if is_common_axis_commutator(r) or is_common_axis_equation(r):
            gx,gy = r.generator(0), r.generator(1)
            CA[gx].add(gy)
            CA[gy].add(gx)
    # CA_P is now in CA represented by adjacency sets.
    # Now check if CA_P is connected via depth-first search.
    novel = [True for i in range(n)]
    def depth_first(x):
        if novel[x]:
            novel[x] = False
            for y in CA[x]:
                depth_first(y)
    depth_first(0)
    return all(not x for x in novel)
