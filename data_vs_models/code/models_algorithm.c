#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "my_sort.h"
#include "models_algorithm.h"


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// RFIM


void RFIM_dynamics(int N, double R, int *size_sequence, double *fields)
{

  double E;
  int S;

  initialization(fields, R, N);

  E = -1.0;
  while(fields[0] > 0)
  {

    // Dynamics. Our implementation proceeds sequentially: we sample avalanche-by-avalanche and each avalanche is sampled shell-by-shell
    next_avalanche(fields, &S, N, E);

    E += 2.0*(double)S/(double)N;

    size_sequence[0] += 1;
    size_sequence[size_sequence[0]] = S;

  }

}



void initialization(double *fields, double R, int N)
{
  int i;

  for(i=1; i<=N; i++) fields[i] = gaussian_number(0, R);
  q_sort_double(fields, 1, N); 
  fields[0] = N;  

}




double find_seed(double *fields, double E)
{
  double epsilon = 1e-06;
  double H;

  double E1 = fields[(int)fields[0]];
  fields[0] --;
  double E2 = fields[(int)fields[0]];

  double Dh = E1 - E2;
  if (Dh > epsilon || Dh < 0) { H = - E - E1 + epsilon; }
  else {H = - E - E1 + Dh/2; }

  return H;
}




void next_avalanche(double *fields, int *S, int N, double E)
{

  int s, shell_size, t;
  double H;

  s = 1;
  t = 0;

  H = find_seed(fields, E);
  E += 2.0/(double)N;

  shell_size = 1;
  // while(shell_size != 0)
  while(shell_size > 0)
  {
    shell_size = next_shell(fields, H, E);

    s += shell_size;
    t += 1;
    E += 2.0*shell_size/(double)N;
  }


  S[0] = s;

}





int next_shell(double *fields, double H, double E)
{
  int s = 0;
  double next_field = fields[(int)fields[0]];
  while(next_field + E + H > 0 && fields[0] >= 1)
  {
    s += 1;
    fields[0] --;
    next_field = fields[(int)fields[0]];
  }
  return s;
}





double gaussian_number(double mu, double sigma)
{

  double u1 = (double)genrand64_real3();
  double u2 = (double)genrand64_real3();

  double Z_0 = sqrt(-2*log(u1)) * cos(2*M_PI*u2);

  return Z_0*sigma + mu;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


// BP


void BP(int N, double z, int *size_sequence)
{

  int Stot, Stree, s;

  Stot = 0;
  while(Stot < N)
  {

    Stree = s = 1;
    Stot += s;
    while(s > 0) 
    {
      s = next_generation(s, z, Stot, N);
      Stree += s;
      Stot += s;
      if (Stot == N) break;
    }

    size_sequence[0] += 1;
    size_sequence[size_sequence[0]] = Stree;

  }

}



int next_generation(int s, double n, int Stot, int N)
{

  int i, num_offsprings, diff;


  int Sshell = 0;
  for(i=1; i<=s; i++)
  {

    num_offsprings = poisson_number(n);
    Sshell += num_offsprings;
    Stot += num_offsprings;
    if (Stot >= N)
    {
      diff = Stot - N;
      Sshell -= diff;
      break;
    }

  }

  return Sshell;

}



int poisson_number(double lambda)
{
  double x, p, s, u;

  x = 0.0;
  p = exp(-lambda);
  s = p;

  u = genrand64_real1();

  while(u > s)
  {
    x ++;
    p *= lambda/x;
    s += p;
  }

  return (int)x;

}
