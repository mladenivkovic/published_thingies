#ifndef GODUNOV_H
#define GODUNOV_H

/* number of boundary cells */
#define NBC 1
#define NBCT 2*NBC

extern void compute_new_states();
extern void set_boundaries();
extern void compute_conserved_states();
extern double compute_dt();

#endif
