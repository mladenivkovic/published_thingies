#include <stdio.h>

#include "gas.h"
#include "godunov.h"
#include "params.h"
#include <math.h>

#ifdef RIEMANN_EXACT
#include "riemann-exact.h"
#elif defined RIEMANN_TRRS
#include "riemann-trrs.h"
#elif defined RIEMANN_TSRS
#include "riemann-tsrs.h"
#elif defined RIEMANN_HLL
#include "riemann-hll.h"
#elif defined RIEMANN_HLLC
#include "riemann-hllc.h"
#endif

extern pstate* w_old;
extern pstate* w_new;
extern pstate* w_intercell;
extern cstate* u_old;
extern cstate* u_new;
extern cstate* flux;
extern params pars;
extern double dx;
extern double dt;
extern double gamma;

extern double *x;





/* ================================ */
void compute_conserved_states(){
/* ================================ */
  /* Computes conserved states from */
  /* primitive states               */
  /*--------------------------------*/

  for (int i=NBC; i<pars.nx+NBCT; i++){
    u_old[i].rho = w_old[i].rho;
    u_old[i].rhou = w_old[i].rho * w_old[i].u;
    u_old[i].E = energy(&w_old[i]);
  }


}

/* =================================== */
void compute_new_states(){
/* =================================== */
  /* Compute the new pstates           */
  /*-----------------------------------*/

  double dtdx = dt/dx;

/* TODO remove */
/* printf("dt= %g dx=%g dtdx=%g\n", dt, dx, dtdx); */
/* printf("GODUNOV CHECK\n"); */
  for (int i=NBC; i<pars.nx+NBC; i++){
    u_new[i].rho = u_old[i].rho + dtdx*(flux[i].rho - flux[i+1].rho);
    u_new[i].rhou = u_old[i].rhou + dtdx*(flux[i].rhou - flux[i+1].rhou);
    u_new[i].E = u_old[i].E + dtdx*(flux[i].E - flux[i+1].E);
/* TODO remove */
/* printf("%lf \t %lf \t %lf \n",u_old[i].rho, u_old[i].rhou, u_old[i].E); */
/* printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",flux[i].rho, flux[i].rhou, flux[i].E, flux[i+1].rho, flux[i+1].rhou, flux[i+1].E); */
/* printf("%lf \t %lf \t %lf \n",u_new[i].rho, u_new[i].rhou, u_new[i].E); */
/* printf("\n"); */
/*  */
    w_new[i].rho = u_new[i].rho;
    w_new[i].u = u_new[i].rhou/u_new[i].rho;
    w_new[i].p = (gamma-1)*(u_new[i].E - 0.5*u_new[i].rhou*w_new[i].u);
  }

/* TODO remove */
/* printf("GODUNOV CHECK END=============================================================\n"); */
}



/* ================================= */
double compute_dt(){
/* ================================= */
  /* Compute new dt                  */
  /*---------------------------------*/

  double vmax = 0;
  
  for (int i = NBC; i<pars.nx+NBC; i++){
#ifdef HLL
    double SL, SR;
    compute_wave_speeds(u_old[i-1], u_old[i], &SL, &SR);
    if (fabs(SL)>vmax) vmax = fabs(SL);
    if (fabs(SR)>vmax) vmax = fabs(SR);
#else
    double u, a, S;
    u = w_old[i].u;
    a = soundspeed(&w_old[i]);
    S = fabs(u)+a;
    if (S>vmax) vmax = S;
#endif
  }
  
  return(pars.ccfl * dx/vmax);


}


/* ==========================*/
void set_boundaries(){
/* ==========================*/
  /* Set boundary conditions */
  /*-------------------------*/

  /* ----------------------- */
  /* transmissive boundaries */
  /* ----------------------- */

  int re = pars.nx+NBC; /* right edge of actual sim */

  if (NBC==0) return;
  for (int i=0; i<NBC; i++){
    /* left boundary    */
    w_old[NBC-i-1].rho = w_old[NBC+i].rho;
    w_old[NBC-i-1].u = w_old[NBC+i].u;
    w_old[NBC-i-1].p = w_old[NBC+i].p;

    /* right boundary   */
    w_old[re+i].rho = w_old[re-i-1].rho;
    w_old[re+i].u = w_old[re-i-1].u;
    w_old[re+i].p = w_old[re-i-1].p;

    /* left boundary    */
    u_old[NBC-i-1].rho = u_old[NBC+i].rho;
    u_old[NBC-i-1].rhou = u_old[NBC+i].rhou;
    u_old[NBC-i-1].E = u_old[NBC+i].E;

    /* right boundary   */
    u_old[re+i].rho = u_old[re-i-1].rho;
    u_old[re+i].rhou = u_old[re-i-1].rhou;
    u_old[re+i].E = u_old[re-i-1].E;
  }
 
}
