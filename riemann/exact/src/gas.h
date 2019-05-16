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
} state;


// extern void init_constants(constants* c);
extern void init_states(state *left, state *right, state* starL, state* starR);
extern double soundspeed(state* s);


#endif
