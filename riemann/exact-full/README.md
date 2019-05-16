Exact Riemann Solver
========================

Computes all primitive variables for an ideal gas of a Riemann problem using an
exact Riemann solver for given initial conditions.
This is the full Riemann solver that can handle both vacuum and non-vacuum ICs.
If parameters 'tmax' and 'nsteps' are given in the parameter file, it will
create 'nsteps' output in linear spacings from 0 to 'tmax' in arbitrary units.
