/* Relations, constants etc concerning ideal gasses */

#ifndef GAS_H
#define GAS_H

// typedef struct {
//   double gamma;
// } constants;


typedef struct {
  double rho;
  double u;
  double p;
} pstate; /*primitive state*/

typedef struct {
  double rho;
  double rhou;
  double E;
} cstate; /*conserved state*/





// extern void init_constants(constants* c);
extern void init_pstates(pstate *left, pstate *right, pstate* starL, pstate* starR);
extern double soundspeed(pstate* s);
extern double energy(pstate* s);


#endif
