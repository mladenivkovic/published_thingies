Godunov Solver
========================

Simple 1D upwind-differencing godunov solver using different Riemann solvers.
Use `presentation_compare_all.py` or `presentation_compare_exact.py` to 
compare results to exact results; The scripts will however assume that you
have a ../riemann/exact-full directory relative to the scripts paths; And 
that that the results of the Riemann solver are present. (Just run run.sh in
the Riemann/bin directory)


Usage:
./godunov paramfile ic-file

paramfile can contain the following values (default values are given):
```
verbose = 0      // How talkative the code should be. 1 = true, 0 = false
nx      = 100    // Number of cells to use
gamma   = 1.4    // adiabatic index
ccfl    = 0.9    // courant factor; dt = ccfl * dx / vmax
nsteps  = 0      // Up to how many steps to do. If = 0, run until t >= tmax
tmax    = 0.1    // Up to which time to simulate
foutput = 0      // Frequency of writing outputs. If = 0, will only write initial and final steps.    
twopstate_ic = 1 // If initial conditions contain Riemann-style only two neighbouring states.
                 // currently this is the only implemented option.

```


NOTE: The methods used can be disappointingly inaccurate.
For example, if you let most of the simulations run to an end time longer then
currently specified in the ./bin/run.sh file, some weird stuff happens. Especially
for the 123 problem and the two shocks problem.
Also, for the HLL solver, the wave speed estimates aren't good for the 123problem.
It'll produce NANs because it produces negative densities :( 
