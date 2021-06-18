def verify(Q):
    a,b,c,d = Q(1,0), Q(0,1), Q(1,1), Q(-1,1)
    e,f,g,h = b+c-a, b+d-a, a+c-b, a+d-b
    return e > 0 and f > 0 and g > 0 and h > 0

def fold_until(op,e,stop,M):
    if stop(M):
        return e
    (a,b,c,d) = M
    L = fold_until(op,e,stop,(a,a+b,c,c+d))
    R = fold_until(op,e,stop,(a+b,b,c+d,d))
    return op(M,L,R)

def union(l):
    a = []
    for s in l:
        a = a + s
    return a

def short_slopes(Q, bound):
    assert verify(Q)
    J = [(1,0),(0,1),(1,1),(-1,1)]
    Jp = [v for v in J if not (Q(v[0],v[1]) > bound)]
    A = (1,1,0,1)
    B = (1,0,1,1)
    C = (-1,-1,1,0)
    D = (0,-1,1,1)
    too_big = lambda M: Q(M[0]+M[1], M[2]+M[3]) > bound
    branch_set = lambda M,L,R: L + [M] + R
    short_descendants = lambda M: fold_until(branch_set,[],too_big,M)
    matrices = union(list(map(short_descendants, [A,B,C,D])))
    return Jp + [(M[0]+M[1], M[2]+M[3]) for M in matrices]
