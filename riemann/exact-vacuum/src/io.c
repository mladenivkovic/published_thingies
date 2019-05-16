#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "gas.h"
#include "params.h"
#include "riemann.h"

#define MAX_LINE_SIZE 200


extern double gamma;

/* ======================================================== */
void read_cmdlineargs(int argc, char* argv[], params* p){
/* ======================================================== */
  /* This function reads in the command line arguments and  */
  /* stores them in the globalparams struct                 */
  /*--------------------------------------------------------*/

  if (argc < 3){
    printf("Too few arguments given. Run this program with ./riemann paramfile datafile\n");
    exit(600);
  }
  else {
    strcpy(p->paramfilename, argv[1]);
    strcpy(p->datafilename, argv[2]);
  };

  printf("Started computation for initial conditions file %s\n", p->datafilename);

}



/* ======================================================== */
void read_ic(params* p, state* left, state* right){
/* ======================================================== */
  /* Read in initial conditions file, store read states.    */
  /*--------------------------------------------------------*/

  //open file
  FILE *par = fopen(p->datafilename, "r");

  // check if file exists
  if (par == NULL) { 
    printf("Error: file '%s' not found.\n", p->datafilename);
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
void read_paramfile(params* p){
/* ============================================================ */
  /* Read in parameter file, store read in global parameters.   */
  /*------------------------------------------------------------*/

  // open file
  FILE *par = fopen(p->paramfilename, "r");

  // check if file exists
  if (par == NULL) { 
    printf("Error: file '%s' not found.\n", p->paramfilename);
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
      p->verbose = atoi(varvalue);
    } 
    else if (strcmp(varname, "nsteps")==0){
      p->nsteps += atoi(varvalue);
    }
    else if (strcmp(varname, "tmax")==0){
      p->tmax = atof(varvalue);
    }
    else if (strcmp(varname, "nx")==0){
      p->nx = atoi(varvalue);
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



/* ================================================================================================= */
void write_output(int step, double t, double* x, double* rho, double* u, double* p, params* par){
/* ================================================================================================= */
  /* Write output of step at time t.                                                                 */
  /*-------------------------------------------------------------------------------------------------*/


  /*-------------------*/
  /* Generate filename */
  /*-------------------*/

  int dot = 0;
  for (int i = strlen(par->datafilename); i > 0; i--){
    if (par->datafilename[i] == '.'){
      dot = i;
      break;
    }
  }

  int slash = 0;
  for (int i = 0; i < (int) strlen(par->datafilename); i++){
    if (par->datafilename[i] == '/'){
      slash = i;
    }
  }

  if (dot==0) dot = strlen(par->datafilename);
  if (slash > 0) slash += 1;

  char filename[MAX_FNAME_SIZE] = "";
  strncpy(filename, par->datafilename+slash, dot-slash);

  /* create output directory */
  mkdir(filename, 0777);

  char fcpy[MAX_FNAME_SIZE] = "";
  strcpy(fcpy, filename);
  strcat(filename, "/");
  strcat(filename, fcpy);




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
  for (int i=0; i<par->nx; i++){
    fprintf(outfilep, "%12.5lf %12.5lf %12.5lf %12.5lf\n", x[i], rho[i], u[i], p[i]);
  }

  fclose(outfilep);

}
