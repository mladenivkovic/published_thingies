#include <stdio.h>
#include <math.h>

#include "gas.h"

extern double gamma;

/* ========================================= */
/* void init_constants(constants* c){ */
/* ========================================= */
  /* Initialize default values for constants */
  /*-----------------------------------------*/
/*  */
/*   c->gamma = 1.4;    [> 7/5 for diatomic ideal gas <] */
/*  */
/* } */




/* ===================================================================== */
void init_states(state *left, state *right, state* starL, state* starR){
/* ===================================================================== */
  /* This function initializes the states.                               */
  /*---------------------------------------------------------------------*/

  left->rho  = 0;
  left->u    = 0;
  left->p    = 0;

  right->rho = 0;
  right->u   = 0;
  right->p   = 0;

  starL->rho  = 0;
  starL->u    = 0;
  starL->p    = 0;

  starR->rho  = 0;
  starR->u    = 0;
  starR->p    = 0;
}







/* ========================================= */
double soundspeed(state* s){
/* ========================================= */
  /* compute sound speed of ideal gas        */
  /*-----------------------------------------*/
  return(sqrt(gamma * s->p / s->rho));
}
