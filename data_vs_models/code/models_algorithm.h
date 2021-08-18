void RFIM_dynamics(int N, double R, int *size_sequence, double *fields);

void initialization(double *fields, double R, int N);
double find_seed(double *fields, double E);
void next_avalanche(double *fields, int *avalanche, int N, double E);
int next_shell(double *fields, double H, double E);

double gaussian_number(double mu, double sigma);





void BP(int N, double z, int *size_sequence);
int next_generation(int s, double n, int Stot, int N);
int poisson_number(double lambda);