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
    if G.countGenerators() == 2:
        for i in range(G.countRelations()):
            if is_common_axis_commutator(G.relation(i)):
                return True
            if is_common_axis_equation(G.relation(i)):
                return True
    return False
    
