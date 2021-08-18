#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "my_sort.h"
#include "models_algorithm.h"
#include "distrib_struct.h"
#include "stats_functions.h"
#include "basic_functions.h"
#include "fitting_protocol.h"


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Construction of theoretical models


void build_theoretical_models(struct distribution *Q1, struct distribution *Q2, struct distribution *temp, int *indices, int *size_sequence, double *fields,\
                              int N, int M, int Rvals, int Zvals, double *Rs, double dR, double *Zs, double dz, double h)
{

  int i;

  for(i=1; i<=Rvals; i++) reset_norms(&Q1[i]);
  for(i=1; i<=Zvals; i++) reset_norms(&Q2[i]);

  build_RFIM(Q1, temp, indices, size_sequence, fields, N, M, Rvals, Rs, dR, h);
  build_BP(Q2, temp, indices, size_sequence, N, M, Zvals, Zs, dz, h);

}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Fit one time series to one dynamical model: maximize likelihood and compute KS distance


int make_fit(struct distribution *Q, struct distribution *P, double *l, double *ks, int Smin_indx, int Smin, int *winParam, int ParamVals)
{
  
  int indx;
  double Pnorm, Qnorm, Pshift, Qshift;

  *l = param_given_Smin(P, Q, winParam, Smin_indx, Smin, ParamVals);
  if(*winParam  == -1) return 0;


  Pnorm = 1.0 - P->cdf[Smin_indx-1];
  Pshift = P->cdf[Smin_indx-1];


  find_closest_index(Q[*winParam].x, Smin, 1, Q[*winParam].num, &indx);
  if (Q[*winParam].x[indx] < Smin ) indx ++;
  if (indx > Q[*winParam].num)
  { 
    *ks = 1.0 - ((P->cdf[Smin_indx] - Pshift) / Pnorm);
  }else{
    Qnorm = 1.0 - Q[*winParam].cdf[indx-1];  
    Qshift = Q[*winParam].cdf[indx-1];

    *ks = KS(P, &Q[*winParam], Pnorm, Qnorm, Pshift, Qshift, Smin_indx);
  }

  return 1;

}





double param_given_Smin(struct distribution *P, struct distribution *Q, int *winParam, int min_indxP, int Smin, int ParamVals)
{

  int i;
  int indxQ;

  double Pnorm, Qnorm;
  double l_val, l_min;


  // Smin = P->x[min_indxP];
  Pnorm = 1.0 - P->cdf[min_indxP - 1];

  l_min = 1e12;
  *winParam = -1;
  for(i=1; i<=ParamVals; i++)
  {

    find_closest_index(Q[i].x, Smin, 1, Q[i].num, &indxQ);
    if (Q[i].x[indxQ] != Smin && P->x[min_indxP] == Smin) 
    { 
      l_val = 1e30; 
      // printf("NO Smin match\n");
    }
    else{
      if(Q[i].x[indxQ] < Smin) indxQ ++;
      if (indxQ > Q[i].num) { l_val = 1e30; }
      else{
        Qnorm = 1.0 - Q[i].cdf[indxQ - 1];
        l_val = negative_log_likelihood(P, &Q[i], Pnorm, Qnorm, min_indxP);
      }
    }

    if (l_val < l_min)
    {
      l_min = l_val;
      *winParam = i;
    }
  }


  return l_min;

}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// p-value of one fit


double single_fit_pval(struct distribution *Q, struct distribution *submodel, struct distribution *temp, int *TS, int *TSredux, int *TSsynt, int *indices,\
 double ks_val, int Smin, int winParam, int ParamVals)
{
  int i, success;
  int SminIndx, submodel_len;
  int winParam_temp;

  double p, ntail;
  double kl_temp, ks_temp;


  fill_TSredux(TS, TSredux, Smin);


  find_closest_index(Q[winParam].x, Smin, 1, Q[winParam].num, &SminIndx);
  if (Q[winParam].x[SminIndx] < Smin) SminIndx ++;
  
  submodel_len = Q[winParam].num - SminIndx + 1;
  distrib_realloc(submodel, submodel_len+1);
  fill_submodel(Q[winParam], submodel, SminIndx);


  p = 0;
  ntail = 1.0 - (double)TSredux[0] / (double)TS[0];


  for(i=1; i<=100; i++)
  {

    generate_synt_sample(submodel, TSredux, TSsynt, ntail, TS[0]);

    PDF(temp, TSsynt, indices);
    q_sort_int_with_doublelst(temp->x, temp->pdf, 1, temp->num);
    CDF(temp);
    reset_indices(temp, indices);

    find_closest_index(temp->x, Smin, 1, temp->num, &SminIndx);
    if (temp->x[SminIndx] < Smin) SminIndx ++;

    if (SminIndx > temp->num) ks_temp = 2.0;
    else{
      success = make_fit(Q, temp, &kl_temp, &ks_temp, SminIndx, Smin, &winParam_temp, ParamVals);
      if (!success) ks_temp = 2.0;
    }

    if (ks_temp > ks_val) p += 1.0;

    reset_norms(temp);

  }

  reset_norms(submodel);

  p /= 100.0;

  return p;

}

