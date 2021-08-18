void build_theoretical_models(struct distribution *Q1, struct distribution *Q2, struct distribution *temp, int *indices, int *size_sequence, double *fields,\
                              int N, int M, int Rvals, int Zvals, double *Rs, double dR, double *Zs, double dz, double h);

int make_fit(struct distribution *Q, struct distribution *P, double *l, double *ks, int Smin_indx, int Smin, int *winParam, int ParamVals);
double param_given_Smin(struct distribution *P, struct distribution *Q, int *winParam, int min_indxP, int Smin_val, int ParamVals);


double single_fit_pval(struct distribution *Q, struct distribution *submodel, struct distribution *temp, int *TS, int *TSredux, int *TSsynt, int *indices,\
 double ks_val, int Smin, int winParam, int ParamVals);


