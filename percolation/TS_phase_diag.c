#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h> 
#include "percolation_algorithm.h"


int main (int argc, char **argv)
{

  int i, j, K, bins;
  FILE *fin, *fout;
  char input_file[400];

  sprintf(input_file,"../%s_ts.txt", argv[1]);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


  // int count = 0;
  K = 1;

  bins = (int)(5*1e03/1.);
  double *time_series = (double *)malloc((K+1)*sizeof(double));
  double **order_param = (double **)malloc(4*sizeof(double *));     //OP[2][i] = avg{ max size at temperature[i] }; OP[3][i] = avg{ (max size at temperature[i])^2 }
  double **OP = (double **)malloc(4*sizeof(double *));              //OP[2][i] = avg{ max size at temperature[i] }; OP[3][i] = avg{ (max size at temperature[i])^2 }

  for(i=0;i<4;i++)
  {
    order_param[i] = (double *)malloc((bins+1)*sizeof(double));
    OP[i] = (double *)malloc((bins+1)*sizeof(double));
    for(j=0;j<=bins;j++) {order_param[i][j] = 0.0; OP[i][j] = 0.0;}
  }



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



  int q;
  double val;
  fin = fopen(input_file, "r");
  if (!fin)
  {
    printf("Can not open input file.\n");
    exit(0);
  }

  for(;;)
  {

    q = fscanf(fin,"%d ",&K);
    if (q<=0) break;

    time_series = realloc(time_series, K*sizeof(double));
    for(j=0; j<K; j++)
    {
      q = fscanf(fin,"%lf ", &val);
      time_series[j] = val;
      if (q<=0) break;
    }

    if (K > 2)
    {
      log_cluster_time_series (K, time_series, bins, order_param);

      for(i=1;i<4;i++)
      {
        for(j=0;j<=bins;j++)
        {
          if (order_param[0][j] > 0.0)
          {
            OP[i][j] += order_param[i][j]/order_param[0][j];
          }
          order_param[i][j] = 0.0;
        }
      }
      for (j=0; j<=bins; j++)
      {
        if (order_param[0][j] > 0.0) OP[0][j] += 1;
        order_param[0][j] = 0.0;
      }
    }
  }

  fclose(fin);
  

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  // output_file_from_dataset_name(argv, filename);

  char filename[200];
  sprintf(filename,"./res/%s_PhaseDiag.txt",argv[1]);


  fout = fopen(filename,"w");
fprintf(fout, "# Order parameter of percolation in dimension 1.\n\
# The order parameter is measured as the size of largest cluster, normalized to its maximum value.\n\
# The maximum value correspond to the number of events (lattice sites).\n\
# The model is %s. Parameter setting is in the file name.\n\
# Columns are: -- measure counter (normalization) - threshold value - mean order parameter - mean squared order parameter --\n", argv[1]);
  for(j=0;j<=bins;j++)
  {
    if(OP[0][j]>0.0)
    {
    	fprintf(fout,"%.12f\t",exp(OP[1][j]/OP[0][j]));
    	fprintf(fout,"%.16f\t",OP[2][j]/OP[0][j]);
    	fprintf(fout,"%.16f\n",OP[3][j]/OP[0][j]);
    }
  }
  fflush(fout);
  fclose(fout);

  for(i=0;i<4;i++) {free(order_param[i]); free(OP[i]);}
  free(order_param);
  free(OP);
  free(time_series);

  return 0;
}
