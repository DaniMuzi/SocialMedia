void PDF(struct distribution *distrib, int *X, int *indices);
void logbox_smoothing(struct distribution *Q, struct distribution *test, int *indices, double h, int N);
void CDF(struct distribution *distrib);
void reset_indices(struct distribution *distrib, int *indices);
void reset_norms(struct distribution *distrib);


double KS(struct distribution *P, struct distribution *Q, double Pnorm, double Qnorm, double Pshift, double Qshift, int min_indx);



double negative_log_likelihood(struct distribution *P, struct distribution *Q, double Pnorm, double Qnorm, int min_indx);



