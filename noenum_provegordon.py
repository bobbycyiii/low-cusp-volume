# from enuminternals.enuminternals import enumerate_internal_necklaces
from gordon.superexceptional import *
from snappy import *
from time import process_time

reduced_neckls = ['s596', 's647', 's774', 's776', 's780', 's782', 's785', 'v2124', 'v2355', 'v2533', 'v2644', 'v2731', 'v3108', 'v3127', 'v3211', 'v3376']
beginning = process_time()
part2_start = process_time()
sups = set()
for neckl in reduced_neckls:
    if neckl == 's776':
        print("*** The fillings of s776 are treated elsewhere. ***\n")
        continue
    print("*** Working on {0} ***".format(neckl))
    neckl_start = process_time()
    sups = sups.union(find_superexceptional_fillings(neckl, verbose=True))
    print("*** Done with {0} in {1} seconds. ***\n".format(neckl,
                                                           process_time()
                                                           - neckl_start))

print("*** Part 2: Finished in {0} seconds. ***\n".format(process_time()
                                                              - part2_start))
print("The superexceptional manifolds are as follows: {0}".format(sups))
print("\nAll in all, that took {0} seconds.".format(process_time()
                                                        - beginning))
