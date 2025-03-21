import regina

def make_dipyr(n):
    newt = regina.Triangulation3()
    for i in range(0,n):
        newt.newTetrahedron()
    for i in range(0,n):
        me = newt.simplex(i)
        you = newt.simplex((i+1)%n)
        me.join(2,you,regina.Perm4(2,3))
    return newt
        
def which_tet(face):
    return face // 2
    
def make_necklace(match):
    n = len(match)
    ndipyr = make_dipyr(n)
    local = list(match)
    while local != []:
        (i,j) = local.pop()
        me = ndipyr.tetrahedron(which_tet(i))
        you = ndipyr.tetrahedron(which_tet(j))    
        me_face = i % 2
        if (i % 2 == j % 2):
            me.join(me_face,you,regina.Perm4(2,3))
        else:
            me.join(me_face,you,regina.Perm4(0,1))
    return ndipyr

def smush(l,m):
    L = list(l)
    M = list(m)
    while M != []:
        L.append(M.pop())
    return L

def old_for_all_perfect_matchings(p,l,m,f):
    pfm = old_for_all_perfect_matchings
    if l == []:
        if m == []:
            f(p)
    elif l[1:] != []:
        # The head of a Python list is its last element.
        x0, x1, l2 = l[-1], l[-2], l[:-2]
        pfm(p+[(x0,x1)], smush(l2,m), [], f)
        pfm(p, l2+[x0], m+[x1], f)

def for_all_perfect_matchings(p,l,m,f):
    pfm = for_all_perfect_matchings
    rv = lambda x: list(reversed(x))
    x,y = l[-1], l[-2]
    if len(l) == 2:
        if m == []:
            f(p + [(x,y)])
        else:
            pfm(p + [(x,y)], rv(m), [], f)
    else:
        pfm(p, l[:-2] + [x], m + [y], f)
        pfm(p + [(x,y)], rv(m) + l[:-2], [], f)
            
    
def enumerate_isosigs(n):
    labels = list(range(2*n))
    sigs = {}
    def put_in_sigs(p):
        mfld = make_necklace(p)
        sigs[mfld.isoSig()] = p
    for_all_perfect_matchings([],labels,[],put_in_sigs)
    return sigs

import sys
import json
if __name__ == """__main__""":
    n = int(sys.argv[1])
    sigs = enumerate_isosigs(n)
    l = list(sigs)
    l.sort()
    fn = sys.argv[2]
    with open(fn, 'w', encoding='utf-8') as f:
        for sig in l:
            f.write(sig + "\n")


