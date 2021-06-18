import regina

def is_homology_T2xI(M):
    mfld = regina.Triangulation3(M)
    mfld.idealToFinite()
    mfld.intelligentSimplify()
    if not mfld.isConnected():
        return False
    x = mfld.countBoundaryComponents()
    if x != 2:
        return False
    h1 = mfld.homology()
    if h1.str() != '2 Z':
        return False
    h1r = mfld.homologyRel()
    if h1r.str() != 'Z':
        return False
    bh1 = mfld.homologyBdry()
    if bh1.str() != '4 Z':
        return False
    return True

def is_embedded(edge):
    source = edge.face(0,0)
    sink = edge.face(0,1)
    return source.index() != sink.index()

def get_embedded_boundary_edge(mfld):
    bdycpts = mfld.boundaryComponents()
    for k in bdycpts:
        edges = k.faces(1)
        for e in edges:
            if is_embedded(e):
                return e.index()
    else:
        return None

def get_coembedded_boundary_edge(mfld):
    bdycpts = mfld.boundaryComponents()
    for k in bdycpts:
        edges = k.faces(1)
        for e in edges:
            if mfld.closeBook(e, perform=False):
                return e.index()
    else:
        return None

def simplify_boundary(mfld):
    mfld.intelligentSimplify()
    em = get_embedded_boundary_edge(mfld)
    while em != None:
        mfld.layerOn(em)
        coem = get_coembedded_boundary_edge(mfld)
        while coem != None:
            mfld.closeBook(coem)
            coem = get_coembedded_boundary_edge(mfld)
        em = get_embedded_boundary_edge(mfld)

def is_T2xI(mfld):
    T = regina.Triangulation3(mfld)
    T.idealToFinite()
    T.intelligentSimplify()
    if not is_homology_T2xI(T):
        return False
    simplify_boundary(T)
    k = T.boundaryComponent(0)
    for e in k.faces(1):
        S = regina.Triangulation3(T)
        f = S.face(1,e.index())
        S.closeBook(f, False, True)
        if not S.isSolidTorus():
            return False
    return True

if __name__ == "__main__":
    print("non T2xIs")
    nont2is = ['hLLLQkbeegefgghhhhhhgb','bGaj']
    for nont2i in nont2is:
        m = regina.Triangulation3(nont2i)
        print("Working on {0}".format(nont2i))
        x = is_T2xI(m)
        print(str(x))
        if x:
            print("is homology t2i: {0}".format(is_homology_T2xI(m)))
            raise Exception("is_T2xI is wrong: {0}".format(nont2i))

    print("T2xIs")
    t2is = ['gfLPIadfdefrdwun', 'mfLLjzLOQcdffhkjiljllrtariqtgvfg']
    for t2i in t2is:
        m = regina.Triangulation3(t2i)
        if not is_T2xI(m):
            raise Exception("is_T2xI is wrong: {0}".format(t2i))

