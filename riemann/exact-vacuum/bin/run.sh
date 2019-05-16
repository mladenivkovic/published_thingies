#!/bin/bash

make
if [[ $? -ne 0 ]]; then
    exit 1
fi
echo "############################################################################"


genparamfile() {
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax

    f=paramfile.txt
    echo "// parameter file for riemann program" > $f
    echo "" >> $f
    echo "verbose = 1" >> $f
    echo "nx = 1000" >> $f
    echo "gamma = 1.4" >> $f
    echo "nsteps = $1" >> $f
    echo "tmax = $2" >> $f


}

genparamfile 8 0.2
./riemann-vacuum paramfile.txt ../ic/left_vacuum.dat
../plot_riemann_solution.py left_vacuum

genparamfile 8 0.2
./riemann-vacuum paramfile.txt ../ic/right_vacuum.dat
../plot_riemann_solution.py right_vacuum

genparamfile 8 0.2
./riemann-vacuum paramfile.txt ../ic/vacuum_generating.dat
../plot_riemann_solution.py vacuum_generating
