/* ======================================== 
 *  
 * An exact Riemann solver.
 *
 * Usage: ./riemann paramfile.txt ic.dat
 *  
 * ======================================== */



#include <stdlib.h>
#include <stdio.h>

#include "gas.h"
#include "io.h"
#include "params.h"
#include "riemann.h"


double gamma = 1.4;   /* Ratio of specific heats; Adiabatic exponent... */
int vacuum = 0;       /* whether we have a vacuum case or not */


/*==========================================================================*/
void print_results(state* left, state* right, state* starL, state* starR){
/*==========================================================================*/
  /* print the 2 resulting star states between the initial states.          */
  /*------------------------------------------------------------------------*/
    
  printf("        %12s %12s %12s %12s\n", "Left", "Star left", "Star right", "Right");
  printf(" rho =  %12.5lf %12.5lf %12.5lf %12.5lf \n", left->rho, starL->rho, starR->rho, right->rho);
  printf(" u   =  %12.5lf %12.5lf %12.5lf %12.5lf \n", left->u,   starL->u,   starR->u,   right->u);
  printf(" p   =  %12.5lf %12.5lf %12.5lf %12.5lf \n", left->p,   starL->p,   starR->p,   right->p);
  printf("=======================================================================\n");

}




/* ====================================== */
int main(int argc, char* argv[]){
/* ====================================== */

  params p;
  init_params(&p);

  /* constants c; */
  /* init_constants(&c); */

  state left, right, starL, starR;
  init_states(&left, &right, &starL, &starR);

  read_cmdlineargs(argc, argv, &p);
  read_paramfile(&p);
  if (p.verbose) print_params(&p);

  read_ic(&p, &left, &right);
  check_ic(&left, &right);


  if (!vacuum) {
    compute_star_state(&left, &right, &starL, &starR);
    print_results(&left, &right, &starL, &starR);
  }

  if (p.nsteps > 0) {
  
    double *X = malloc(p.nx*sizeof(double));
    double *RHO = malloc(p.nx*sizeof(double));
    double *U = malloc(p.nx*sizeof(double));
    double *P = malloc(p.nx*sizeof(double));

    /* initialize x values */
    double dx = 2./(p.nx-1);
    for (int i=0; i<p.nx; i++){
      X[i] = i*dx - 1;
    }
    
    for (int istep=0; istep<p.nsteps+1; istep++){
      if (p.verbose) printf("Computing and writing output for step=%d\n", istep);
      double t = ((double) istep)/((double) p.nsteps)*p.tmax;
      if (vacuum) {
        compute_solution_vacuum(t, X, RHO, U, P, &p, &left, &right);
      }
      else{
        compute_solution(t, X, RHO, U, P, &p, &left, &right, &starL, &starR);
      }
      write_output(istep, t, X, RHO, U, P, &p);
    }


  }


  return(0);

}




