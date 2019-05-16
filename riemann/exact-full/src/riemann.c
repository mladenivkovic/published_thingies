#include "params.h"
#include "gas.h"
#include "riemann.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern double gamma;
extern int vacuum;


/* ====================================================== */
void check_ic(state *left, state *right){
/* ====================================================== */
  /* Check whether we work with vacuum                    */
  /*------------------------------------------------------*/

  if (left->rho == 0 && left->p == 0){
    printf("Found left vacuum state.\n");
    vacuum = 1;
    return;
  }
  if (right->rho == 0 && right->p == 0){
    printf("Found right vacuum state.\n");
    vacuum = 1;
    return;
  }

  double delta_u = right->u - left->u;
  double u_crit = 2./(gamma - 1) * (soundspeed(left) + soundspeed(right));

  if (delta_u < u_crit){
    printf("no vacuum detected.\n");
  }
  else {
    printf("delta_u >= u_crit: We have vacuum generating conditions.\n");
    vacuum = 1;
  }
}




/* ========================================================================================== */
void compute_star_state( state *left, state *right, state* starL, state* starR){
/* ========================================================================================== */
  /* computes the star state given the left and right states.                                 */
  /*------------------------------------------------------------------------------------------*/

  double AL = 2. / ((gamma + 1) * left->rho);
  double AR = 2. / ((gamma + 1) * right->rho);
  double BL = (gamma - 1)/(gamma + 1) * left->p;
  double BR = (gamma - 1)/(gamma + 1) * right->p;
  double aL = soundspeed(left);
  double aR = soundspeed(right);

  double delta_u = right->u - left->u;

  double tolerance = 1e-6;

  double pguess, pold;
  double fL, fR, dfpdpL, dfpdpR;


  pguess = 0.5*(left->p + right->p);

  int niter = 0;

  do {
    niter += 1;
    fL = fp(pguess, left,  gamma, AL, BL, aL);
    fR = fp(pguess, right, gamma, AR, BR, aR);
    dfpdpL = dfpdp(pguess, left,  gamma, AL, BL, aL);
    dfpdpR = dfpdp(pguess, right, gamma, AR, BR, aR);
    pold = pguess;
    pguess -= fpfull(fL, fR, delta_u)/dfpdpfull(dfpdpL, dfpdpR);
    if (pguess<tolerance) pguess=tolerance;
    /* printf("pstar iter %d pguess=%10.6lf, pold=%10.6lf\n", niter, pguess, pold); */
  }
  while (2*fabs((pguess-pold)/(pguess+pold)) >= tolerance);
  printf("p* found after %d iterations.\n", niter);


  starL->p = pguess;
  starR->p = pguess;

  fL = fp(pguess, left,  gamma, AL, BL, aL);
  fR = fp(pguess, right, gamma, AR, BR, aR);
  double u = 0.5*((left->u + right->u) + fR - fL);
  starL->u = u;
  starR->u = u;

  starL->rho = rho_star(left,  starL);
  starR->rho = rho_star(right, starR);
  

}



/* ==================================================== */
double rho_star( state *s, state *star){
/* ==================================================== */
  /* Compute the density in the star region             */
  /*----------------------------------------------------*/

  double pdiv = (star->p/s->p);
  if (star->p > s->p) {
    /* shocking matters */
    double gamfact = (gamma - 1)/(gamma + 1);
    return s->rho * (gamfact + pdiv) / (gamfact*pdiv + 1);
  } 
  else{
    /* rare occasions */
    return s->rho * pow(pdiv, (1./gamma));
  }

}





/* ========================================================== */
double fpfull(double fpL, double fpR, double delta_u){
/* ========================================================== */
  /* Full pressure function                                   */
  /*----------------------------------------------------------*/
  return fpL + fpR + delta_u;
}


/* ================================================== */
double dfpdpfull(double dfpdpL, double dfpdpR){
/* ================================================== */
  /* Full derivative of pressure function             */
  /*--------------------------------------------------*/
  return dfpdpL + dfpdpR;
}






/* =============================================================================== */
double fp(double pguess, state *s, double gamma, double A, double B, double a){
/* =============================================================================== */
  /* Left/Right part of the pressure function                                      */
  /*-------------------------------------------------------------------------------*/

  if (pguess > s->p){
    /* we have a shock situation */
    return (pguess - s->p)*sqrt(A/(pguess + B));
  }
  else{
    /* we have a rarefaction situation */
    return 2 * a / (gamma - 1) * ( pow(pguess/s->p, 0.5*(gamma-1)/gamma) - 1 );
  }
}


/* =============================================================================== */
double dfpdp(double pguess, state *s, double gamma, double A, double B, double a){
/* =============================================================================== */
  /* First derivative of Left/Right part of the pressure function                  */
  /*-------------------------------------------------------------------------------*/

  if (pguess > s->p){
    /* we have a shock situation */
    return sqrt(A/(pguess + B)) * (1 + 0.5 * (pguess - s->p)/(pguess + B));
  }
  else{
    /* we have a rarefaction situation */
    return 1./(s->rho * a) * ( pow(pguess/s->p, -0.5*(gamma+1)/gamma) );
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






/* ================================================================================ */
void compute_solution(double t, double* x, double* rho, double* u, double* p, \
  params* par, state* left, state* right, state* starL, state* starR){
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


  double al =  soundspeed(left);
  double ar =  soundspeed(right);
  double asl = soundspeed(starL);
  double asr = soundspeed(starR);

  double SHL = left->u - al;    /* speed of head of left rarefaction fan */
  double STL = starL->u - asl;  /* speed of tail of left rarefaction fan */
  double SL  = left->u  - al*sqrt(0.5*(gamma+1)/gamma * starL->p/left->p  + 0.5*(gamma-1)/gamma); /* left shock speed */
  double SHR = right->u + ar;   /* speed of head of right rarefaction fan */
  double STR = starR->u + asr;  /* speed of tail of right rarefaction fan */
  double SR  = right->u + ar*sqrt(0.5*(gamma+1)/gamma * starR->p/right->p + 0.5*(gamma-1)/gamma); /* right shock speed */

  for (int i=0; i<par->nx; i++){
    double S = x[i]/t;

    if (S < starL->u){
      /*------------------------*/
      /* We're on the left side */
      /*------------------------*/
      if (starL->p <= left->p){
        /*------------------*/
        /* left rarefaction */
        /*------------------*/
        if (S < SHL) {
          /* we're outside the rarefaction fan */
          rho[i] = left->rho;
          u[i] = left->u;
          p[i] = left->p;
        }
        else {
          if (S < STL){
            /* we're inside the fan */
            rho[i] = rho_fanL(S, left);
            u[i] = u_fanL(S, left);
            p[i] = p_fanL(S, left);
          }
          else{
            /* we're in the star region */
            rho[i] = starL->rho;
            u[i] = starL->u;
            p[i] = starL->p;
          }
        }
      }
      else{
        /*------------------*/
        /* left shock       */
        /*------------------*/
        if (S<SL){
          /* we're outside the shock */
          rho[i] = left->rho;
          u[i] = left->u;
          p[i] = left->p;
        }
        else{
          /* we're in the star region */
          rho[i] = starL->rho;
          u[i] = starL->u;
          p[i] = starL->p;
        }
      }
    }
    else{
      /*-------------------------*/
      /* We're on the right side */
      /*-------------------------*/
      if (starR->p <= right->p){
        /*-------------------*/
        /* right rarefaction */
        /*-------------------*/
        if (S > SHR) {
          /* we're outside the rarefaction fan */
          rho[i] = right->rho;
          u[i] = right->u;
          p[i] = right->p;
        }
        else {
          if (S > STR){
            /* we're inside the fan */
            rho[i] = rho_fanR(S, right);
            u[i] = u_fanR(S, right);
            p[i] = p_fanR(S, right);
          }
          else{
            /* we're in the star region */
            rho[i] = starR->rho;
            u[i] = starR->u;
            p[i] = starR->p;
          }
        }
      }
      else{
        /*------------------*/
        /* right shock      */
        /*------------------*/
        if (S>SR){
          /* we're outside the shock */
          rho[i] = right->rho;
          u[i] = right->u;
          p[i] = right->p;
        }
        else{
          /* we're in the star region */
          rho[i] = starR->rho;
          u[i] = starR->u;
          p[i] = starR->p;
        }
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
