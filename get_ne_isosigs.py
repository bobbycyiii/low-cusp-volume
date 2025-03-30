from enuminternals.ne import solve_problem_ne
from enuminternals.quick_checks import *
import sys
import regina

if __name__ == "__main__":
    sig = sys.argv[1]
    print(f"{sig}: start")
    mfld = regina.Triangulation3(sig)
    if is_nontrivial_link_exterior(mfld) and mfld.countBoundaryComponents() in [2,3]:
        X = ""
        try:
            ne_sigs = solve_problem_ne(mfld)
            for ne_sig in ne_sigs:
                X += ne_sig + "\n"
        except Exception as x:
            X = x
        if len(sys.argv) == 2:
            print(X)
        else:
            with open(sys.argv[2], 'w') as fl:
                fl.write(X)
    print(f"{sig}: finished")
