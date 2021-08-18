#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "distrib_struct.h"
#include "stats_functions.h"
#include "basic_functions.h"




void PDF(struct distribution *distrib, int *X, int *indices)
{

  int i, j, val, len;


  distrib->x[0] = 0;
  for(i=1; i<=X[0]; i++)
  {

    val = X[i];
    distrib->norm += 1.0;

    if (indices[val] == -1)
    {
      distrib->num += 1;
      j = distrib->num;
      indices[val] = j;

      len = j+1;

      distrib_realloc(distrib, len);

    }else{ j = indices[val]; }

    distrib->x[j] = val;
    distrib->pdf[j] += 1.0;

  }

}


void logbox_smoothing(struct distribution *Q, struct distribution *temp, int *indices, double h, int N)
{

  int i, j, val, len;
  int left, right;

  double mid, w;
  double l1, l2, binlen, fact;


  Q->norm = 0.0;
  for (i=1; i<=temp->num; i++)
  {

    mid =  temp->x[i];

    l1 = log(mid) - h;
    l2 = log(mid) + h;

    left = (int)(exp(l1));
    right = (int)(exp(l2));

    if (left < 1) left = 1;
    if (right > N) right = N;

    binlen = right - left;
    if ((int)binlen == 0) binlen = 1.0;

    w = 1.0 / binlen;

    for(val=left; val<=right; val++)
    {

      fact = temp->pdf[i] * w;

      if (indices[val] == -1)
      {
        Q->num += 1;
        j = Q->num;
        indices[val] = j;

        len = j+1;

        distrib_realloc(Q, len);

      }else{ j = indices[val]; }


      Q->x[j] = val;
      Q->pdf[j] += fact;
      Q->norm += fact;


    }

  }


}


void CDF(struct distribution *distrib)
{

  int i;
  double n = 0;
  double norm = distrib->norm;

  for(i=1; i<=distrib->num; i++)
  {
    distrib->pdf[i] /= norm;
    n += distrib->pdf[i];
    distrib->cdf[i] = n;
  }
  distrib->x[0] = 0;
  distrib->pdf[0] = 0;
  distrib->cdf[0] = 0;

}



void reset_indices(struct distribution *distrib, int *indices)
{
  int i, j; 

  for(i=1; i<=distrib->num; i++)
  {
    j = distrib->x[i];
    indices[j] = -1;
  }

}


void reset_norms(struct distribution *distrib)
{

  distrib->num = 0;
  distrib->norm = 0;

}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// Kolmogorov-Smirnov distance (see SM M of XXX for a theoretical explanation of this approach to the KS calculation)


double KS(struct distribution *P, struct distribution *Q, double Pnorm, double Qnorm, double Pshift, double Qshift, int min_indx)
{

  int px;
  int i, indx;

  // int old_indx = min_indx;
  int old_indx = 1;
  
  double py, qy;
  double ks1, ks2, ksmax;

  ksmax = 0;
  for(i=min_indx; i<=P->num; i++)
  {

    px = P->x[i];
    find_exact_index(Q->x, px, old_indx, Q->num, &indx);
    old_indx = indx;

    py = (P->cdf[i] - Pshift) / Pnorm;
    qy = (Q->cdf[indx] - Qshift) / Qnorm;    
    ks1 = fabs(qy - py);



    py = (P->cdf[i-1] - Pshift) / Pnorm;
    qy = (Q->cdf[indx-1] - Qshift) / Qnorm;    
    ks2 = fabs(qy - py);



    if(ks1 > ksmax) ksmax = ks1; 
    if(ks2 > ksmax) ksmax = ks2;

  }

  return ksmax;

}






//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////





double negative_log_likelihood(struct distribution *P, struct distribution *Q, double Pnorm, double Qnorm, int min_indx)
{

  int i;
  int indx, old_indx;

  int px;
  double l, py, qy;


  l = 0.0;
  old_indx = 1;
  for(i=min_indx; i<=P->num; i++)
  {

    px = P->x[i];
    py = P->pdf[i] / Pnorm;


    find_closest_index(Q->x, px, old_indx, Q->num, &indx);
    if (Q->x[indx] != px) 
    {
      return 1e30;
    }

    old_indx = indx;

    qy = Q->pdf[indx] / Qnorm;
    l -= py * log(qy);

  }

  return l;

}






