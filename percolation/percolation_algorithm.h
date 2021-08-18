void merge_left_to_right(int size_small, int small, int large, int indx, int *cluster_size, int *cluster_indx);
void merge_right_to_left(int size_small, int small, int large, int indx, int *cluster_size, int *cluster_indx);
void measure_order_param(double new_delta_value, double start, double scale, int new_S_value, int K, double **vectors, double *OP_values);


void log_cluster_time_series (int K, double *time_series, int bins, double **order_param);


