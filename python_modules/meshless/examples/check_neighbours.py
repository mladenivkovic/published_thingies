#!/usr/bin/env python3


#===============================================================
# Get neighbour lists in naive way by looping over all
# particles and by sorting them out
#===============================================================


import numpy as np
import matplotlib.pyplot as plt
import time

try:
    import meshless as ms
except ImportError:
    print("Didn't find 'meshless' module... be sure to add it to your pythonpath!")
    quit(2)





#---------------------------
# initialize variables
#---------------------------


# temp during rewriting
srcfile = './snapshot_sodshock.hdf5'    # swift output file
ptype = 'PartType0'                 # for which particle type to look for





#========================
def main():
#========================
    
    # get data
    x, y, h, rho, m, ids, npart = ms.read_file(srcfile, ptype)
    H = ms.get_H(h)

    # get neighbours
    start_n = time.time()
    ndata_naive = ms.get_neighbour_data_for_all_naive(x, y, H, fact=1.0, L=1.0, periodic=True)
    stop_n = time.time()

    start_s = time.time()
    ndata_smart = ms.get_neighbour_data_for_all(x, y, H, fact=1.0, L=1.0, periodic=True)
    stop_s = time.time()


    # compare results
    mx_n = ndata_naive.maxneigh
    mx_s = ndata_smart.maxneigh

    if mx_n != mx_s:
        print("Max number of neighbours is different!", mx_n, mx_s)

    neigh_n = ndata_naive.neighbours
    nneigh_n = ndata_naive.nneigh
    neigh_s = ndata_smart.neighbours
    nneigh_s = ndata_smart.nneigh

    found_difference = False

    for p in range(npart):
        nn = nneigh_n[p]
        ns = nneigh_s[p]
        if nn != ns:
            print("Got different number of neighbours for particle ID", ids[p], ":", nn, ns)
            print(neigh_n[p])
            print(neigh_s[p])
            quit()

        for n in range(nneigh_n[p]):
            nn = neigh_n[p][n]
            ns = neigh_s[p][n]
            if nn != ns:
                print("Got different neighbour IDs:", nn, ns)
                print(neigh_n[p])
                print(neigh_s[p])
                found_difference = True

    if not found_difference:
        print("Found no difference.")


    print()
    print('{0:18} {1:18} {2:18}'.format("time naive", "time better", "|1 - naive/better|"))
    tn = stop_n - start_n
    ts = stop_s - start_s
    print('{0:18.4f} {1:18.4f} {2:18.4f}'.format(tn, ts, abs(1 - tn/ts)))

    

if __name__ == '__main__':
    main()
