#ifndef DISTRIBUTION_STRUCT_HEADER
#define DISTRIBUTION_STRUCT_HEADER

struct distribution
{
  int num;
  double norm;

  int *x;
  double *pdf;
  double *cdf;
};

#endif