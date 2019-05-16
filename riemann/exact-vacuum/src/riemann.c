#include "params.h"
#include "gas.h"
#include "riemann.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern double gamma;


/* ====================================================== */
void check_ic(state *left, state *right){
/* ====================================================== */
  /* Check that the initial conditions are proper         */
  /*------------------------------------------------------*/

  if (left->rho == 0 && left->p == 0){
    printf("Found left vacuum state.\n");
    return;
  }
  if (right->rho == 0 && right->p == 0){
    printf("Found right vacuum state.\n");
    return;
  }

  double delta_u = right->u - left->u;
  double u_crit = 2./(gamma - 1) * (soundspeed(left) + soundspeed(right));

  if (delta_u < u_crit){
    printf("delta_u < u_crit: This is not a vacuum generating IC. Quitting.\n");
    exit(0);
  }
  else {
    printf("delta_u >= u_crit: We have vacuum generating conditions.\n");
  }

}





/* ================================================================================ */
void compute_solution_vacuum(double t, double* x, double* rho, double* u, double* p, \
  params* par, state* left, state* right){
/* ================================================================================ */
  /* Compute the solution of the riemann problem at given time t for all x          */
  /*--------------------------------------------------------------------------------*/

  if(t==0){
    for(int i=0; i<par->nx; i++){
      if(x[i]<=0){
        rho[i] = left->rho;
        u[i] = left->u;
        p[i] = left->p;
      }
      else {
        rho[i] = right->rho;
        u[i] = right->u;
        p[i] = right->p;
      }
    }

    return;
  }


  if (left->rho==0){
    /*------------------------*/
    /* Left vacuum state      */
    /*------------------------*/
    double ar = soundspeed(right);
    double SR = right->u - 2*ar/(gamma-1);

    for (int i=0; i<par->nx; i++){

      double S = x[i]/t;

      if (S <= SR){
        /* left vacuum */
        rho[i] = left->rho;
        u[i] = SR;
        p[i] = left->p;
      }
      else if (S < right->u + ar){
        /* inside rarefaction */
        rho[i] = rho_fanR(S, right);
        u[i] = u_fanR(S, right);
        p[i] = p_fanR(S, right);
      }
      else{
        /* original right state */
        rho[i] = right->rho;
        u[i] = right->u;
        p[i] = right->p;
      }
    }
  }

  else if (right->rho==0){
    /*------------------------*/
    /* Right vacuum state     */
    /*------------------------*/

    double al = soundspeed(left);
    double SL = left->u + 2*al/(gamma-1);

    for (int i=0; i<par->nx; i++){

      double S = x[i]/t;

      if (S >= SL){
        /* right vacuum */
        rho[i] = right->rho;
        u[i] = SL;
        p[i] = right->p;
      }
      else if (S > left->u - al){
        /* inside rarefaction */
        rho[i] = rho_fanL(S, left);
        u[i] = u_fanL(S, left);
        p[i] = p_fanL(S, left);
      }
      else{
        /* original left state */
        rho[i] = left->rho;
        u[i] = left->u;
        p[i] = left->p;
      }
    }
  }
  else {
    /*------------------------*/
    /* Vacuum generating case */
    /*------------------------*/

    double al = soundspeed(left);
    double ar = soundspeed(right);
    double SL = left->u + 2*al/(gamma-1);
    double SR = right->u - 2*ar/(gamma-1);

    for (int i=0; i<par->nx; i++){

      double S = x[i]/t;

      if (S <= left->u-al){
        /* left original state*/
        rho[i] = left->rho;
        u[i] = left->u;
        p[i] = left->p;
      }
      else if (S < SL){
        /* rarefaction fan from right to left */
        rho[i] = rho_fanL(S, left);
        u[i] = u_fanL(S, left);
        p[i] = p_fanL(S, left);
      }
      else if (S < SR) {
        /* vacuum region */
        rho[i] = 0;
        u[i] = 0.5*(SL+SR); /* just made something up here */
        p[i] = 0;
      }
      else if (S < right->u + ar){
        /* rarefaction fan from left to right */
        rho[i] = rho_fanR(S, right);
        u[i] = u_fanR(S, right);
        p[i] = p_fanR(S, right);
      }
      else{
        /* right original state */
        rho[i] = right->rho;
        u[i] = right->u;
        p[i] = right->p;
      }
    }
  }
}




/* ==================================== */
double rho_fanL(double S, state* s){
/* ==================================== */
  /* Compute rho inside rarefaction fan */
  /*------------------------------------*/
  return s->rho * pow((2./(gamma+1) + (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2/(gamma-1)));
}

/* ==================================== */
double u_fanL(double S, state* s){
/* ==================================== */
  /* Compute u inside rarefaction fan   */
  /*------------------------------------*/
  return 2/(gamma+1) * (soundspeed(s) + 0.5*(gamma-1)*s->u + S);
}

/* ==================================== */
double p_fanL(double S, state* s){
/* ==================================== */
  /* Compute p inside rarefaction fan   */
  /*------------------------------------*/
  return s->p * pow((2/(gamma+1) + (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2*gamma/(gamma-1)));
}





/* ==================================== */
double rho_fanR(double S, state* s){
/* ==================================== */
  /* Compute rho inside rarefaction fan */
  /*------------------------------------*/
  return s->rho * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2./(gamma-1)));
}

/* ==================================== */
double u_fanR(double S, state* s){
/* ==================================== */
  /* Compute u inside rarefaction fan   */
  /*------------------------------------*/
  return 2/(gamma+1) * (-soundspeed(s) + 0.5*(gamma-1)*s->u + S);
}

/* ==================================== */
double p_fanR(double S, state* s){
/* ==================================== */
  /* Compute p inside rarefaction fan   */
  /*------------------------------------*/
  return s->p * pow((2./(gamma+1) - (gamma-1)/(gamma+1)/soundspeed(s) * (s->u - S)), (2*gamma/(gamma-1)));
}
