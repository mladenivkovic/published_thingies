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

genparamfile 8 1
./riemann paramfile.txt ../ic/sod_test.dat
../plot_riemann_solution.py sod_test

genparamfile 8 0.1
./riemann paramfile.txt ../ic/sod_test_modified.dat
../plot_riemann_solution.py sod_test_modified

genparamfile 8 0.2
./riemann paramfile.txt ../ic/123problem.dat
../plot_riemann_solution.py 123problem

genparamfile 8 0.025
./riemann paramfile.txt ../ic/left_blast_wave.dat
../plot_riemann_solution.py left_blast_wave

genparamfile 8 0.1
./riemann paramfile.txt ../ic/right_blast_wave.dat
../plot_riemann_solution.py right_blast_wave

genparamfile 8 0.035
./riemann paramfile.txt ../ic/two_shocks.dat
../plot_riemann_solution.py two_shocks

# genparamfile 8 0.2
# ./riemann paramfile.txt ../ic/left_vacuum.dat
# ../plot_riemann_solution.py left_vacuum
#
# genparamfile 8 0.2
# ./riemann paramfile.txt ../ic/right_vacuum.dat
# ../plot_riemann_solution.py right_vacuum
#
# genparamfile 8 0.2
# ./riemann paramfile.txt ../ic/vacuum_generating.dat
# ../plot_riemann_solution.py vacuum_generating
#
