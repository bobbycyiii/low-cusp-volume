import regina
from .t2i import is_T2xI
       
def is_nonseparating_fault(surf):
    T = surf.cutAlong()
    isFault = surf.eulerChar() >= 0
    return (T.isConnected() and
            surf.isCompact() and
            isFault)

def is_nonseparating_closed_fault(surf):
    isClosed = not surf.hasRealBoundary()
    return isClosed and is_nonseparating_fault(surf)

def is_nonseparating_nonclosed_fault(surf):
    isNotClosed = surf.hasRealBoundary()
    return isNotClosed and is_nonseparating_fault(surf)

def is_sphere(surf):
    return (surf.isConnected()
            and surf.isCompact()
            and surf.eulerChar() == 2)

def unsum(sphere):
    T = sphere.cutAlong()
    L, R = T.triangulateComponents()
    L.finiteToIdeal()
    L.intelligentSimplify()
    L.idealToFinite()
    L.intelligentSimplify()
    R.finiteToIdeal()
    R.intelligentSimplify()
    R.idealToFinite()
    R.intelligentSimplify()
    return (L.isoSig(), R.isoSig())

def is_essential_sphere(surf):
    if not is_sphere(surf):
        return False
    T = surf.cutAlong()
    if T.isConnected():
        return True
    Lsig, Rsig = unsum(surf)
    L = regina.Triangulation3(Lsig)
    R = regina.Triangulation3(Rsig)
    if L.isSphere() or R.isSphere():
        return False
    return True

def is_torus_fault(surf):
    is_torus = (surf.isConnected()
                and surf.isCompact()
                and surf.isOrientable()
                and not surf.hasRealBoundary()
                and surf.eulerChar() == 0)
    if not is_torus:
        return False
    T = surf.cutAlong()
    T.intelligentSimplify()
    if T.hasCompressingDisc():
        return False
    if T.isConnected():
        return True
    L, R = T.triangulateComponents()
    if is_T2xI(L) or is_T2xI(R):
        return False
    return True

def is_mobius_band(surf):
    return (surf.isConnected()
            and surf.isCompact()
            and not surf.isOrientable()
            and surf.hasRealBoundary()
            and surf.eulerChar() == 0)

def is_solid_torus_annulus(surf):
    is_annulus = (surf.isConnected()
                  and surf.isCompact()
                  and surf.isOrientable()
                  and surf.hasRealBoundary()
                  and surf.eulerChar() == 0)
    if not is_annulus:
        return False
    T = surf.cutAlong()
    T.intelligentSimplify()
    if T.isConnected():
        return False
    L, R = T.triangulateComponents()
    if L.isSolidTorus() and R.isSolidTorus():
        return True
    return False
