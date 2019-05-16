Exact Riemann Solver
========================

Computes all primitive variables for an ideal gas of a Riemann problem using an
exact Riemann solver for given initial conditions.
This is a solver for vacuum cases only: Either left state vacuum, right state
vacuum or vacuum generating conditions.
If parameters 'tmax' and 'nsteps' are given in the parameter file, it will
create 'nsteps' output in linear spacings from 0 to 'tmax' in arbitrary units.
