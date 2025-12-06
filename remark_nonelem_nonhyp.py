from enuminternals.bead import enumerate_isosigs 
from enuminternals.quick_checks import is_nontrivial_link_exterior
from enuminternals.ne import solve_problem_ne
from gordon.hyperbolicity import hyp_info
import regina
import snappy


if __name__ == "__main__":
    sigs = {}

    for i in range(4,8):
        sigs = {**sigs, **enumerate_isosigs(i)}

    n = len(sigs)
    i = 0
    for sig in sigs:
        i += 1
        print(f"{i} of {n}")
        rmfld = regina.Triangulation3(sig)
        if not is_nontrivial_link_exterior(rmfld):
            continue
        x = solve_problem_ne(rmfld)
        if x != set(): 
            smfld = snappy.ManifoldHP(sig) 
            x = hyp_info(smfld,0,(0,0))
            if not x[0] == True:
                print(f"{sig} is not hyperbolic, but fills to hyperbolic\nmatching: {sigs[sig]}\n")
            
