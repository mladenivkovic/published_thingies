/* all stuff concerning Riemann problem */

#ifndef RIEMANN_H
#define RIEMANN_H
#include "gas.h"


extern void check_ic(state *left, state *right);
extern void compute_solution_vacuum(double t, double* x, double* rho, double* u, double* p, \
  params* par, state* left, state* right);
extern double rho_fanL(double S, state* s);
extern double u_fanL(double S, state* s);
extern double p_fanL(double S, state* s);
extern double rho_fanR(double S, state* s);
extern double u_fanR(double S, state* s);
extern double p_fanR(double S, state* s);
#endif
