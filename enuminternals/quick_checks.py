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
    """
    
    if G.countGenerators() == 2:
        for i in range(G.countRelations()):
            if is_common_axis_commutator(G.relation(i)):
                return True
            if is_common_axis_equation(G.relation(i)):
                return True
    return False
    """    
    G = mfld.fundamentalGroup()
    G.intelligentSimplify()
    n = G.countGenerators()
    gen_commute_graph_components = [[i,1] for i in range(n)]
    for r in G.relations():
        if is_common_axis_commutator(r) or is_common_axis_equation(r):
            gx = r.generator(0)
            gy = r.generator(1)
            unify(gen_commute_graph_components, gx, gy)
    cpts = [find(gen_commute_graph_components, i) for i in range(n)]
    return all(cpt == cpts[0] for cpt in cpts)
    
def unify(L, u, v):
    x, y = find(L, u), find(L, v)
    if x != y:
        if L[x][1] > L[y][1]:
            x,y = y,x
        L[x][0] = y
        L[y][1] += L[x][1]
        
def find(L, u):
    if L[u][0] != u:
        L[u][0] = find(L, L[u][0])
    return L[u][0]

