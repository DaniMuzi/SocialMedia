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



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Functions to build theoretical models


void build_RFIM(struct distribution *Q, struct distribution *temp, int *indices, int *size_sequence, double *fields, int N, int M, int Rvals, double *Rs, double dR, double h)
{

  int i, j;
  double R;

  for(i=1; i<=Rvals; i++)
  {

    for (j=0; j<M; j++)
    {
      size_sequence[0] = 0;
      R = (Rs[i] - dR) + genrand64_real3() *  2.0 * dR;
      RFIM_dynamics(N, R, size_sequence, fields);
      PDF(temp, size_sequence, indices);
    }

    q_sort_int_with_doublelst(temp->x, temp->pdf, 1, temp->num);
    CDF(temp);
    reset_indices(temp, indices);

    logbox_smoothing(&Q[i], temp, indices, h, N);
    q_sort_int_with_doublelst(Q[i].x, Q[i].pdf, 1, Q[i].num);
    CDF(&Q[i]);
    reset_indices(&Q[i], indices);

    reset_norms(temp);

  }

}




void build_BP(struct distribution *Q, struct distribution *temp, int *indices, int *size_sequence, int N, int M, int Zvals, double *Zs, double dz, double h)
{

  int i, j;
  double z;

  for(i=1; i<=Zvals; i++)
  {

    for (j=0; j<M; j++)
    {
      size_sequence[0] = 0;
      z = (Zs[i] - dz) + genrand64_real3() *  2.0 * dz;
      BP(N, z, size_sequence);
      PDF(temp, size_sequence, indices);
    }

    q_sort_int_with_doublelst(temp->x, temp->pdf, 1, temp->num);
    CDF(temp);
    reset_indices(temp, indices);

    logbox_smoothing(&Q[i], temp, indices, h, N);
    q_sort_int_with_doublelst(Q[i].x, Q[i].pdf, 1, Q[i].num);
    CDF(&Q[i]);
    reset_indices(&Q[i], indices);

    reset_norms(temp);

  }

}




//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Generic functions to read data and to allocate memory required for computations


int read_one_TS(FILE *fin, int **TS, int **TSredux, int **TSsynt, int **size_sequence, int **indices, double **fields, int *N, int *Nav)
{
  
  int i, q, success;
  int N_temp, Nav_temp;
  double val;

  q = fscanf(fin, "%d ", &N_temp);
  if (q<=0) return 0;

  q = fscanf(fin, "%d ", &Nav_temp);
  if (q<=0) { printf("Unexpected format in file! N has been read but Nav has not. Why?\n"); return -1; }

  *N = N_temp;
  *Nav = Nav_temp;
  success = make_allocations_for_data(TS, TSredux, TSsynt, size_sequence, indices, fields, *N, *Nav);
  if (!success) return 0;
  for(i=1; i<=*N; i++) (*indices)[i] = -1;


  (*TS)[0] = *Nav;
  for(i=1; i<=*Nav; i++)
  {
    q = fscanf(fin,"%lf ", &val);
    (*TS)[i] = val;
    if (q<=0) { printf("Unexpected format in file! N and Nav has been read but the TS has not been (fully). Why?\n"); return 0; };
  }

  return 1;

}





int make_allocations_for_data(int **TS, int **TSredux, int **TSsynt, int **size_seq, int **indices, double **fields, int N, int Nav)
{

  *TS = realloc(*TS, (Nav + 1)*sizeof(int));
  *TSsynt = realloc(*TSsynt, (Nav + 1)*sizeof(int));
  *TSredux = realloc(*TSredux, (Nav + 1)*sizeof(int));
  *size_seq = realloc(*size_seq, (N + 1)*sizeof(int));
  *indices = realloc(*indices, (N + 1)*sizeof(int));
  *fields = realloc(*fields, (N + 1)*sizeof(double));

  if (!*TS) {
    printf("TS realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
  if (!*TSsynt) {
    printf("TS realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
  if (!*TSredux) {
    printf("TS realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
  if (!*size_seq) {
    printf("Size_seq realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
  if (!*indices) {
    printf("Indices realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
  if (!*fields) {
    printf("Fields realloc: FAIL!\n");
    fflush(stdout);
    return 0;
  }
 
  return 1;

}



void distrib_realloc(struct distribution *distrib, int len)
{

  distrib->x = realloc(distrib->x, len*sizeof(int));
  distrib->pdf = realloc(distrib->pdf, len*sizeof(double));
  distrib->cdf = realloc(distrib->cdf, len*sizeof(double));

  if (!distrib->x) {
    printf("x realloc: FAIL!\n");
    fflush(stdout);
  }
  if (!distrib->pdf) {
    printf("pdf realloc: FAIL!\n");
    fflush(stdout);
  }
  if (!distrib->cdf) {
    printf("cdf realloc: FAIL!\n");
    fflush(stdout);
  }
  distrib->pdf[len-1] = 0.0;

}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Functions to construct synthetic samples


void fill_TSredux(int *TS, int *TSredux, int Smin)
{

  int i;

  TSredux[0] = 0;
  for(i=1; i<=TS[0]; i++)
  {
    if (TS[i] < Smin)
    {
      TSredux[0] ++;
      TSredux[TSredux[0]] = TS[i];
    }
  }

}



void fill_submodel(struct distribution Q, struct distribution *submodel, int SminIndx)
{

  int i;
  double Qnorm, Qshift, val;

  Qnorm = 1.0 - Q.cdf[SminIndx-1];
  Qshift = Q.cdf[SminIndx-1];

  submodel->num = 0;
  for(i=SminIndx; i<=Q.num; i++)
  {
    val = (Q.cdf[i] - Qshift) / Qnorm;
    submodel->num += 1;
    submodel->x[submodel->num] = Q.x[i];
    submodel->cdf[submodel->num] = val;
  }

}




void generate_synt_sample(struct distribution *submodel, int *TSredux, int *TSsynt, double ntail, int Nav)
{

  int i, indx, val;
  double p, r;

  TSsynt[0] = 0;
  for(i=1; i<=Nav; i++)
  {

    p = genrand64_real3();
    r = genrand64_real3();
    if (p <= ntail)
    {
      find_closest_index_double(submodel->cdf, r, 1, submodel->num, &indx);
      if (submodel->cdf[indx] < r) indx ++;
      val = submodel->x[indx];
    }else
    {
      indx = (int)(r*(double)TSredux[0])+1;
      if (indx > TSredux[0]) indx = 1;
      val = TSredux[indx];
    }

    TSsynt[0] ++;
    TSsynt[TSsynt[0]] = val;

  }

}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Binary search


void find_exact_index(int *list, int val, int low, int high, int *target)
{

  int mid;
  
  if ((int)list[low] == val)
  {
    target[0] = low;
    return;
  } 

  if ((int)list[high] == val)
  {
    target[0] = high;
    return;
  }
  
  mid = (int)((double)(low + high) / 2.0);
  if ((int)list[mid] == val || mid == low || mid == high)
  { 
    target[0] = mid;
    return;
  }

  if ((int)list[low] < val && val < list[mid])
  {
    high = mid-1;
    find_exact_index(list, val, low, high, target);
  }

  if ((int)list[mid] < val && val < list[high])
  {
    low = mid+1;
    find_exact_index(list, val, low, high, target);
  }


  return;

}



void find_closest_index(int *list, int val, int low, int high, int *target)
{

  int mid;

  if(val < list[low])
  {
    target[0] = low;
    return;
  }
  if(val > list[high])
  {
    target[0] = high;
    return;
  }

  while(low <= high)
  {
    mid = (int)((double)(low + high)/2.0);

    if(val < list[mid]){ high = mid - 1; }
    else if (val > list[mid]) { low  = mid + 1; }
    else 
    {
      target[0] = mid;
      return;
    }
  }

  if( (list[low] - val) < (val - list[high]) ) {target[0] = low;}
  else {target[0] = high;}
  
}





void find_closest_index_double(double *list, double val, int low, int high, int *target)
{

  int mid;

  if(val < list[low])
  {
    target[0] = low;
    return;
  }
  if(val > list[high])
  {
    target[0] = high;
    return;
  }

  while(low <= high)
  {
    mid = (int)((double)(low + high)/2.0);

    if(val < list[mid]){ high = mid - 1; }
    else if (val > list[mid]) { low  = mid + 1; }
    else 
    {
      target[0] = mid;
      return;
    }
  }

  if( (list[low] - val) < (val - list[high]) ) {target[0] = low;}
  else {target[0] = high;}
  
}


