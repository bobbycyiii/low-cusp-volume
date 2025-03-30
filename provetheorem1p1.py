from enuminternals import enumerate_internal_necklaces
from gordon.superexceptional import *
from snappy import *
from time import process_time

if __name__ == "__main__":
    beginning = process_time()
    print("*** Part 0: Enumerating the internal necklace structures ***")
    part0_start = process_time()
    neckls = set()
    for bead in [4,5,6,7]:
        neckls = neckls.union(enumerate_internal_necklaces(bead, verbose=True))
    print("------------------------------------")
    print("| The ancestral set is as follows: |")
    print("------------------------------------\n")
    names = list(neckls)
    names.sort()
    for neckl in names:
        print(neckl)
    print("\n*** Part 0: Finished in {0} seconds. ***\n".format(process_time()
                                                                - part0_start))

    print("\n*** Part 1: Reducing the ancestral set. ***\n")
    part1_start = process_time()
    reduced_neckls = set()
    for neckl in names:
        M = snappy.Manifold(neckl)
        n = M.dual_curves()
        matched = False
        for i in range(len(n)):
            x = M.drill(i)
            for y in x.identify():
                if y.name() in names:
                    print("{0} is a Dehn filling of {1}".format(neckl, y.name()))
                    matched = True
        if not matched:
            reduced_neckls.add(neckl)
    print("\n*** Part 1: Finished in {0} seconds. ***\n".format(process_time()
                                                                - part1_start))

