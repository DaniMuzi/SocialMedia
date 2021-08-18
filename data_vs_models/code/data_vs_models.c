#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mt64.h"
#include "my_sort.h"
#include "models_algorithm.h"
#include "distrib_struct.h"
#include "stats_functions.h"
#include "basic_functions.h"
#include "fitting_protocol.h"



int main (int argc, char **argv)
{

  //sets initial seed for random number generator
  init_genrand64((unsigned) time(NULL) * getpid());

  if (argc != 5) 
  { 
    printf("Wrong number of inputs!\n"); 
    printf("Specify the name of the dataset, the value of the temporal resolution, the value of S_min and the file number!\n"); 
    return -1; 
  }

  double delta = atof(argv[2]);
  int Smin = atoi(argv[3]);
  int filenum  = atoi(argv[4]);


  int i;
  int M = 500;                 // Number of realization per configuration of each model

  double h = 0.1;              // Scale of the smoothing kernel
  double R, dR, z, dz;
  int Rvals = 108;
  int Zvals = 113;


  char out_file_name[1000];
  char in_file_name[1000];


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


  // Vectors for models and data


  int *TS = malloc(1*sizeof(int));
  int *TSsynt = malloc(1*sizeof(int));
  int *TSredux = malloc(1*sizeof(int));

  int *size_sequence = malloc(1*sizeof(int));
  double *fields = malloc(1*sizeof(double));

  double *Rs = malloc((Rvals+1)*sizeof(double));
  double *Zs = malloc((Zvals+1)*sizeof(double));


//////////////////////////////////////////////////////////////////////////////////


  // Variables for distributions


  struct distribution *Q1 = malloc((Rvals+1)*sizeof(struct distribution));

  for(i=0; i<=Rvals; i++)
  {
    Q1[i].x = malloc((1)*sizeof(int));
    Q1[i].pdf = malloc((1)*sizeof(double));
    Q1[i].cdf = malloc((1)*sizeof(double));
    Q1[i].pdf[0] = 0;

    Q1[i].num = 0;
    Q1[i].norm = 0.0;
  }


  struct distribution *Q2 = malloc((Zvals+1)*sizeof(struct distribution));

  for(i=0; i<=Zvals; i++)
  {
    Q2[i].x = malloc((1)*sizeof(int));
    Q2[i].pdf = malloc((1)*sizeof(double));
    Q2[i].cdf = malloc((1)*sizeof(double));
    Q2[i].pdf[0] = 0;

    Q2[i].num = 0;
    Q2[i].norm = 0.0;
  }



  struct distribution P;
  P.x = malloc((1)*sizeof(int));
  P.pdf = malloc((1)*sizeof(double));
  P.cdf = malloc((1)*sizeof(double));
  P.pdf[0] = 0;

  P.num = 0;
  P.norm = 0.0;



  struct distribution temp;
  temp.x = malloc((1)*sizeof(int));
  temp.pdf = malloc((1)*sizeof(double));
  temp.cdf = malloc((1)*sizeof(double));
  temp.pdf[0] = 0;

  temp.num = 0;
  temp.norm = 0.0;



  struct distribution submodel;
  submodel.x = malloc((1)*sizeof(int));
  submodel.pdf = malloc((1)*sizeof(double));
  submodel.cdf = malloc((1)*sizeof(double));
  submodel.pdf[0] = 0;

  submodel.num = 0;
  submodel.norm = 0.0;



  int *indices = malloc(1*sizeof(int));



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


  // Parameters


  R = 2.7;
  dR = 0.025;
  i = 1;
  while(R > 0.001)
  { 
    Rs[i] = R;
    R -= dR;
    i ++;
  }
  dR /= 2.0;
  

  z = 1.7;
  dz = 0.015;
  i = 1;
  while(z > 0.01)
  { 
    Zs[i] = z;
    z -= dz;
    i ++;
  }
  dz /= 2.0;


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


  // Input and Output files


  sprintf(out_file_name, "../res/%s_fit_Smin%d_delta%g_num%d.txt", argv[1], Smin, delta, filenum);
  FILE *fout = fopen(out_file_name, "w");
  if (!fout){
      printf("Cant open out file\n");
      return -1; 
  }



  sprintf(in_file_name, "../inp/%s_delta%g_model_selection_input_num%d.txt", argv[1], delta, filenum);
  FILE *fin = fopen(in_file_name, "r");
  if (!fin){
      printf("Cant open input file\n");
      return -1;
  }


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


  // Fitting protocol


  int success;

  int TSnum, N, Nold, Nav;
  int winParam1, winParam2;
  int fit1_success, fit2_success;


  double ks1, ks2; 
  double l1, l2;
  double p1, p2;
  double LogLikelRatio;


  int Smin_indx;

  Nold = -1;
  TSnum = 0;
  for(;;)
  {

    success = read_one_TS(fin, &TS, &TSredux, &TSsynt, &size_sequence, &indices, &fields, &N, &Nav);
    if (!success) break;



    PDF(&P, TS, indices);    
    q_sort_int_with_doublelst(P.x, P.pdf, 1, P.num);
    CDF(&P);
    reset_indices(&P, indices);



    if (N != Nold) build_theoretical_models(Q1, Q2, &temp, indices, size_sequence, fields, N, M, Rvals, Zvals, Rs, dR, Zs, dz, h);



    find_closest_index(P.x, Smin, 1, P.num, &Smin_indx);
    if (P.x[Smin_indx] < Smin) Smin_indx ++;


    if (P.num - Smin_indx > 1)
    {
      fit1_success = make_fit(Q1, &P, &l1, &ks1, Smin_indx, Smin, &winParam1, Rvals);
      fit2_success = make_fit(Q2, &P, &l2, &ks2, Smin_indx, Smin, &winParam2, Zvals);

      if (!fit1_success && !fit2_success) { fprintf(fout, "Both likelihoods are zero! %d %d\n", N, Nav); }
      else if (fit1_success && !fit2_success) { fprintf(fout, "RFIM fitted, BP not fitted: the selected model is RFIM. %d %g %.16f\n", N, Rs[winParam1], l1); }
      else if (!fit1_success && fit2_success) { fprintf(fout, "RFIM not fitted, BP fitted: the selected model is BP. %d %g %.16f\n", N, Zs[winParam2], l2); }
      else if (fit1_success && fit2_success)
      {
        p1 = single_fit_pval(Q1, &submodel, &temp, TS, TSredux, TSsynt, indices, ks1, Smin, winParam1, Rvals);
        p2 = single_fit_pval(Q2, &submodel, &temp, TS, TSredux, TSsynt, indices, ks2, Smin, winParam2, Zvals);

        LogLikelRatio = P.norm * (1.0 - P.cdf[Smin_indx-1]) * ( l2 - l1 );  
        fprintf(fout, "%d %g %g %.16f %.16f %g %g %.16f\n", N, Rs[winParam1], Zs[winParam2], l1, l2, p1, p2, LogLikelRatio);

      }
      
    } 


    fflush(fout);

    reset_norms(&P);
    TSnum ++;
    Nold = N;

  }

  fclose(fin);
  
  fflush(fout);
  fclose(fout);



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////



  for(i=0; i<=Rvals; i++)
  {
    free(Q1[i].x);
    free(Q1[i].pdf);
    free(Q1[i].cdf);
  }
  free(Q1);


  for(i=0; i<=Zvals; i++)
  {
    free(Q2[i].x);
    free(Q2[i].pdf);
    free(Q2[i].cdf);
  }
  free(Q2);

  free(P.x);
  free(P.pdf);
  free(P.cdf);

  free(temp.x);
  free(temp.pdf);
  free(temp.cdf);

  free(submodel.x);
  free(submodel.pdf);
  free(submodel.cdf);

  free(indices);

  free(size_sequence);
  free(fields);
  free(TS);
  free(TSsynt);
  free(TSredux);
  free(Rs);
  free(Zs);

  return 0;

}






