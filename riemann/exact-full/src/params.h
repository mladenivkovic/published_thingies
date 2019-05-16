/* All around parameters used in the simulation. */

#ifndef PARAMS_H
#define PARAMS_H
#include "gas.h"

#define MAX_FNAME_SIZE 100

/* GLOBAL PARAMETERS */
typedef struct {
  int verbose;
  int nsteps;
  double tmax;
  int nx;
  char paramfilename[MAX_FNAME_SIZE];
  char datafilename[MAX_FNAME_SIZE];
} params;



// extern void check_params();
extern void init_params(params *p);
extern void print_params(params *p);


#endif
