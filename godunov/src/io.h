/* IO routines */

#ifndef IO_H
#define IO_H
#include "gas.h"
#include "params.h"

extern void read_cmdlineargs(int argc, char* argv[]);
extern void read_ic(pstate* left, pstate* right);
extern void read_paramfile();
void write_output(int step, double t, double* x, pstate* s);

#endif
