/* IO routines */

#ifndef IO_H
#define IO_H
#include "gas.h"
#include "params.h"

extern void read_cmdlineargs(int argc, char* argv[], params* p);
extern void read_ic(params* p, state* left, state* right);
extern void read_paramfile(params* p);
extern void write_output(int step, double t, double* x, double* rho, double* u, double* p, params* par);

#endif
