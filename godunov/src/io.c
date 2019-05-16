#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "gas.h"
#include "godunov.h"
#include "params.h"

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

#define MAX_LINE_SIZE 200


extern double gamma;
extern params pars;

/* ======================================================== */
void read_cmdlineargs(int argc, char* argv[]){
/* ======================================================== */
  /* This function reads in the command line arguments and  */
  /* stores them in the params struct                       */
  /*--------------------------------------------------------*/

  if (argc < 3){
    printf("Too few arguments given. Run this program with ./riemann paramfile datafile\n");
    exit(600);
  }
  else {
    strcpy(pars.paramfilename, argv[1]);
    strcpy(pars.datafilename, argv[2]);
  };

  printf("Started computation for initial conditions file %s\n", pars.datafilename);

}



/* ======================================================== */
void read_ic(pstate* left, pstate* right){
/* ======================================================== */
  /* Read in initial conditions file, store read states.    */
  /*--------------------------------------------------------*/

  //open file
  FILE *par = fopen(pars.datafilename, "r");

  // check if file exists
  if (par == NULL) { 
    printf("Error: file '%s' not found.\n", pars.datafilename);
    exit(2);
  }

  char varname[80] ;
  char varvalue[80] ;
  char tempbuff[MAX_LINE_SIZE] ;

  


  while (fgets(tempbuff,MAX_LINE_SIZE,par))
  // fgets(str_buff, n,filepointer) :
  // gets n characters from file in filepointer and stores them
  // in str_buff.
  // returns 0 if EoF is reached.
  
  {
  
    // check whether tempbuff is empty line
    int isempty = 0;
    for (int i = 0; i<MAX_LINE_SIZE; i++){
      if (tempbuff[i] != ' '){
        if (tempbuff[i] == '\n'){
          isempty = 1;
        }
        break;
      }
    }

    if (isempty) continue;


    sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
    // reads formatted input from a string, writes it in
    // the variables given after the format string.
    // The format used is <string> separator <=> <string> ends with <;>
  
    if (strcmp(varname,"rho_L") == 0) {
      left->rho = atof(varvalue);
    } 
    else if (strcmp(varname, "u_L")==0){
      left->u = atof(varvalue);
    }
    else if (strcmp(varname, "p_L")==0){
      left->p = atof(varvalue);
    }
    else if (strcmp(varname, "rho_R") == 0) {
      right->rho = atof(varvalue);
    } 
    else if (strcmp(varname, "u_R")==0){
      right->u = atof(varvalue);
    }
    else if (strcmp(varname, "p_R")==0){
      right->p = atof(varvalue);
    }
    else if (strcmp(varname, "//")==0) {
      // ignore comments
      continue;
    }
    else if (strcmp(varname, "/*")==0) {
      // ignore comments
      continue;
    }
    else{
      printf("Unrecongized data : \"%s\" | \"%s\"\n", varname, varvalue);
    }
  }

  fclose(par);


}



/* ============================================================ */
void read_paramfile(){
/* ============================================================ */
  /* Read in parameter file, store read in global parameters.   */
  /*------------------------------------------------------------*/

  // open file
  FILE *par = fopen(pars.paramfilename, "r");

  // check if file exists
  if (par == NULL) { 
    printf("Error: file '%s' not found.\n", pars.paramfilename);
    exit(2);
  }

  char varname[80] ;
  char varvalue[80] ;
  char tempbuff[MAX_LINE_SIZE] ;

  


  while (fgets(tempbuff,MAX_LINE_SIZE,par))
  // fgets(str_buff, n,filepointer) :
  // gets n characters from file in filepointer and stores them
  // in str_buff.
  // returns 0 if EoF is reached.
  {
    // check whether tempbuff is empty line
    int isempty = 0;
    for (int i = 0; i<MAX_LINE_SIZE; i++){
      if (tempbuff[i] != ' '){
        if (tempbuff[i] == '\n'){
          isempty = 1;
        }
        break;
      }
    }

    if (isempty) continue;


    sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
    // reads formatted input from a string, writes it in
    // the variables given after the format string.
    // The format used is <string> separator <=> <string> ends with <;>
  

    if (strcmp(varname,"verbose") == 0) {
      pars.verbose = atoi(varvalue);
    } 
    else if (strcmp(varname, "nsteps")==0){
      pars.nsteps += atoi(varvalue);
    }
    else if (strcmp(varname, "tmax")==0){
      pars.tmax = atof(varvalue);
    }
    else if (strcmp(varname, "nx")==0){
      pars.nx = atoi(varvalue);
    }
    else if (strcmp(varname, "foutput")==0){
      pars.foutput = atoi(varvalue);
    }
    else if (strcmp(varname, "ccfl")==0){
      pars.ccfl = atof(varvalue);
    }
    else if (strcmp(varname, "gamma")==0){
      gamma = atof(varvalue);
    }
    else if (strcmp(varname, "//")==0) {
      // ignore comments
      continue;
    }
    else if (strcmp(varname, "/*")==0) {
      // ignore comments
      continue;
    }
    else{
      printf("Unrecongized parameter : \"%s\"\n", varname);
    }
  }

  fclose(par);

}



/* ========================================================= */
void write_output(int step, double t, double* x, pstate* s){
/* ========================================================= */
  /* Write output of step at time t.                         */
  /*---------------------------------------------------------*/

  char filename[MAX_FNAME_SIZE] = ""; 


  if (strlen(pars.outputfilename)==0) {
    /*---------------------------------------------------------*/
    /* Generate output directory filename based on ic filename */
    /*---------------------------------------------------------*/

    printf("Generating outputfilename step %d\n", step);

    int dot = 0;
    /* extract filename without suffix */
    for (int i = strlen(pars.datafilename); i > 0; i--){
      if (pars.datafilename[i] == '.'){
        dot = i;
        break;
      }
    }

    int slash = 0;
    /* remove possible directories paths from filename*/
    for (int i = 0; i < (int) strlen(pars.datafilename); i++){
      if (pars.datafilename[i] == '/'){
        slash = i;
      }
    }

    if (dot==0) dot = strlen(pars.datafilename);
    if (slash > 0) slash += 1;

    /* now copy the exact part that you want into filename string */
    char fname_extract[MAX_FNAME_SIZE] = "";
    strncpy(fname_extract, pars.datafilename+slash, dot-slash);


#ifdef RIEMANN_EXACT
    strcpy(filename, "EXACT");
#elif defined RIEMANN_TRRS
    strcpy(filename, "TRRS");
#elif defined RIEMANN_TSRS
    strcpy(filename, "TSRS");
#elif defined RIEMANN_HLL
    strcpy(filename, "HLL");
#elif defined RIEMANN_HLLC
    strcpy(filename, "HLLC");
#else
    strcpy(filename, "NONE");
#endif
    
    /* create output parent directory */
    mkdir(filename, 0777);

    strcat(filename, "/");
    strcat(filename, fname_extract);

    /* create output directory */
    mkdir(filename, 0777);

    strcat(filename, "/");
    strcat(filename, fname_extract);

    strcpy(pars.outputfilename,filename);
  }
  else{
    strcpy(filename, pars.outputfilename);
  }


  char sstep[3];
  if (step>99) printf("I'm not made to write outputs > 99; I'm probably gonna crash now.\n");
  sprintf(sstep, "%02d", step);
  strcat(filename, "-");
  strcat(filename, sstep);
  strcat(filename, ".out");


  
  /*-------------------*/
  /* Write output file */
  /*-------------------*/

  FILE *outfilep = fopen(filename, "w");
  fprintf(outfilep, "t = %10.4lf\n", t);
  fprintf(outfilep, "%12s %12s %12s %12s\n", "x", "rho", "u", "p");
  for (int i=NBC; i<pars.nx+NBC; i++){
    fprintf(outfilep, "%12.5lf %12.5lf %12.5lf %12.5lf\n", x[i], s[i].rho, s[i].u, s[i].p);
  }
  fclose(outfilep);

}
