#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "params.h"


extern double gamma;


/* ========================================== */
void init_params(params * p){
/* ========================================== */
  /* Initialize parameters to default values  */
  /*------------------------------------------*/

  p->verbose = 0;
  for (int i = 0; i<MAX_FNAME_SIZE; i++){
    p->paramfilename[i] = 0;
    p->datafilename[i] = 0;
  }
  p->nsteps = 0;
  p->tmax = 1;
  p->nx = 100;

}




/* ========================================== */
void print_params(params *p){
/* ========================================== */
  /* Print out current parameters             */
  /*------------------------------------------*/

  printf("==================================================\n");
  printf("Starting calculation. Parameters are:\n");

  printf("Verbose?              ");
  if (p->verbose) {
    printf("True\n");
  } else {
    printf("False\n");
  }
  printf("tmax:                 %g\n", p->tmax);
  printf("nsteps:               %d\n", p->nsteps);
  printf("nx:                   %d\n", p->nx);

  printf("gamma:                %g\n", gamma);

  printf("==================================================\n");

}
