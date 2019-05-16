#!/bin/bash

#-------------------------------------------
# Automatically compile and run the godunov
# program with all solvers and various IC
# files.
# Unless you specify noplot, it will also
# plot the solutions in the directory of 
# each output. To skip plotting, use any of the
# cmd line args
# -n | -np | --np | --noplot | noplot | no | n | np)
#----------------------------------------




#===================================
genparamfile() {
#===================================
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax
    # $3 = output frequency: after how many steps to write

    f=paramfile.txt
    echo "// parameter file for godunov program" > $f
    echo ""             >> $f
    echo "verbose = 1"  >> $f
    echo "nx = 1000"     >> $f
    echo "gamma = 1.4"  >> $f
    echo "ccfl = 0.9"   >> $f
    echo "nsteps = $1"  >> $f
    echo "tmax = $2"    >> $f
    echo "foutput = $3" >> $f
}



#=====================================
# Check for noplot cmd line arg
#=====================================
noplot=false
if [[ $# > 0 ]]; then
    case $1 in

        -n | -np | --np | --noplot | noplot | no | n | np)
            noplot=true
        ;;
        
        *)
            echo "unrecognized cmd line param '" $1 "'" 
            exit
        ;;

    esac
fi




#========================================
# MAIN LOOP
#========================================


#---------------------------------------------
# for SOLVER in HLLC; do
# for SOLVER in HLL; do
# for SOLVER in TSRS; do
# for SOLVER in TRRS; do
# for SOLVER in EXACT; do
for SOLVER in HLL HLLC; do
# for SOLVER in EXACT TRRS TSRS HLL HLLC; do
#---------------------------------------------

    sed -i "s/^RIEMANN=.*/RIEMANN=${SOLVER}/" Makefile

    make clean
    make
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    echo "############################################################################"



    # genparamfile 20 1 1
    # rm -r $SOLVER/sod_test_modified
    # ./godunov paramfile.txt ../ic/sod_test_modified.dat
    # # gdb --args ./godunov paramfile.txt ../ic/sod_test_modified.dat
    # valgrind ./godunov paramfile.txt ../ic/sod_test_modified.dat
    # ../plot_godunov_solution.py $SOLVER/sod_test_modified

    # 30 lines down to (un)comment
    genparamfile 0 1 0
    rm -r $SOLVER/sod_test
    ./godunov paramfile.txt ../ic/sod_test.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/sod_test ; fi;

    genparamfile 0 0.1 0
    rm -r $SOLVER/sod_test_modified
    ./godunov paramfile.txt ../ic/sod_test_modified.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/sod_test_modified ; fi;

    genparamfile 0 0.15 0
    rm -r $SOLVER/123problem
    ./godunov paramfile.txt ../ic/123problem.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/123problem ; fi;

    genparamfile 0 0.025 0
    rm -r $SOLVER/left_blast_wave
    ./godunov paramfile.txt ../ic/left_blast_wave.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/left_blast_wave ; fi;

    genparamfile 0 0.1 0
    rm -r $SOLVER/right_blast_wave
    ./godunov paramfile.txt ../ic/right_blast_wave.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/right_blast_wave ; fi;

    genparamfile 0 0.035 0
    rm -r $SOLVER/two_shocks
    ./godunov paramfile.txt ../ic/two_shocks.dat
    if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/two_shocks ; fi;
    #
    # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/left_vacuum
    # # ./godunov paramfile.txt ../ic/left_vacuum.dat
    # # if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/left_vacuum ; fi;
    #
    # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/right_vacuum
    # # ./godunov paramfile.txt ../ic/right_vacuum.dat
    # # if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/right_vacuum ; fi;
    #
    # # genparamfile 0 0.2 0
    # # rm -r $SOLVER/vacuum_generating
    # # ./godunov paramfile.txt ../ic/vacuum_generating.dat
    # # if [ "$noplot" = false ]; then ../plot_godunov_solution.py $SOLVER/vacuum_generating ; fi;
    #
    #
done;
