void build_RFIM(struct distribution *Q, struct distribution *temp, int *indices, int *size_sequence, double *fields, int N, int M, int Rvals, double *Rs, double dR, double h);
void build_BP(struct distribution *Q, struct distribution *temp, int *indices, int *size_sequence, int N, int M, int Zvals, double *Zs, double dz, double h);



int read_one_TS(FILE *fin, int **TS, int **TSredux, int **TSsynt, int **size_sequence, int **indices, double **fields, int *N, int *Nav);
int make_allocations_for_data(int **TS, int **TSredux, int **TSsynt, int **size_seq, int **indices, double **fields, int N, int Nav);
void distrib_realloc(struct distribution *distrib, int len);



void fill_TSredux(int *TS, int *TSredux, int Smin);
void fill_submodel(struct distribution Q, struct distribution *submodel, int SminIndx);
void generate_synt_sample(struct distribution *submodel, int *TSredux, int *TSsynt, double ntail, int Nav);




void find_exact_index(int *list, int val, int low, int high, int *target);
void find_closest_index(int *list, int val, int low, int high, int *target);
void find_closest_index_double(double *list, double val, int low, int high, int *target);
