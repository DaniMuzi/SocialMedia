#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "my_sort.h"
#include "percolation_algorithm.h"



void merge_left_to_right(int size_small, int small, int large, int indx, int *cluster_size, int *cluster_indx)
{
  int m = 1;
  while(m<=size_small)
  {
    cluster_indx[indx-m] = large;
    cluster_size[small] -= 1;
    cluster_size[large] += 1;
    m++;
  }

}


void merge_right_to_left(int size_small, int small, int large, int indx, int *cluster_size, int *cluster_indx)
{

  int m = 1;
  while(m<=size_small)
  {
    cluster_indx[indx+m-1] = large;
    cluster_size[small] -= 1;
    cluster_size[large] += 1;
    m++;
  }

}


void measure_order_param(double new_delta_value, double start, double scale, int new_S_value, int K, double **vectors, double *OP_values)
{

  int new_indx, old_indx, i;
  double unobserved_delta_value;

  double old_delta_value = OP_values[0];
  double old_S_value = OP_values[1];


  old_indx = (int)((old_delta_value - start)/scale);
  new_indx = (int)((new_delta_value - start)/scale);


  for (i=old_indx+1; i<new_indx; i++)
  {
    if (vectors[0][i] < 0) unobserved_delta_value = vectors[1][i]/vectors[0][i];
    else unobserved_delta_value = start + scale*(i+0.5);
    vectors[0][i] += 1.0;
    vectors[1][i] += unobserved_delta_value;
    vectors[2][i] += (double)old_S_value/(double)(K+1);
    vectors[3][i] += ((double)old_S_value/(double)(K+1)) * ((double)old_S_value/(double)(K+1));
  }

  i = new_indx;
  vectors[0][i] += 1.0;
  vectors[1][i] += new_delta_value;
  vectors[2][i] += (double)new_S_value/(double)(K+1);
  vectors[3][i] += ((double)new_S_value/(double)(K+1)) * ((double)new_S_value/(double)(K+1));


  OP_values[0] = new_delta_value;
  OP_values[1] = new_S_value;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void log_cluster_time_series (int K, double *time_series, int bins, double **order_param)
{
  double delta, old_delta, dt;
  int i, j, c, b, sc, sb, max_S;
  double log_start = log(1e-20);

  double *edge_weight = (double *)malloc((K+1)*sizeof(double));
  int *edge_indx = (int *)malloc((K+1)*sizeof(int));
  int *cluster_indx = (int *)malloc((K+1)*sizeof(int));
  int *cluster_size = (int *)malloc((K+1)*sizeof(int));
  double *OP_values = (double *)malloc(2*sizeof(double));

  cluster_size[0] = 1;
  cluster_indx[0] = 0;
  cluster_size[K] = 0;
  cluster_indx[K] = 0;

  for(i=1;i<K;i++)
  {
  	edge_weight[i] = log_start;
    if (time_series[i] > 0) edge_weight[i] = log(time_series[i]);

    edge_indx[i] = i;
    cluster_indx[i] = i;
    cluster_size[i] = 1;
  }

  q_sort_double_with_indx(edge_weight,edge_indx,1,K-1);

  dt = (log(1e30) - log_start) / (double)bins;

  max_S = 1;
  delta = old_delta = log_start; //log(1e-20);
  OP_values[0] = log_start;
  OP_values[1] = max_S;

  j=1;
  while(j<K)
  {

    i = edge_indx[j];
    delta = edge_weight[j];

    if(delta>old_delta)
    {
      if (old_delta >= log_start) measure_order_param(old_delta, log_start, dt, max_S, 0, order_param, OP_values);
    }

    c = cluster_indx[i-1];
    b = cluster_indx[i];
    sc = cluster_size[c];
    sb = cluster_size[b];

    if (sc<=sb)
    {
      merge_left_to_right(sc, c, b, i, cluster_size, cluster_indx);
      if (cluster_size[b] > max_S ) max_S = cluster_size[b];
    }
    else
    {
      merge_right_to_left(sb, b, c, i, cluster_size, cluster_indx);
      if (cluster_size[c] > max_S ) max_S = cluster_size[c];
    }

    old_delta = delta;
    j ++;

  }

  measure_order_param(delta, log_start, dt, max_S, 0, order_param, OP_values);         
  measure_order_param(log_start + bins*dt, log_start, dt, max_S, 0, order_param, OP_values);


  free(edge_weight);
  free(edge_indx);
  free(cluster_indx);
  free(cluster_size);
  free(OP_values);

}


