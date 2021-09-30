

#ifndef SPKMEANS_C_SPKMEANS_H
#define SPKMEANS_C_SPKMEANS_H

typedef struct {
    double ee;  /*  eigenvalue */
    int index;  /* the ee's (diagonal wise) index in final A' (used for stable sorting, in case when ee's are equal) */
    double *ve; /* eigenvector */
} ee_ve_tuple;
typedef struct {
    double **final_A_tag;
    double **V_matrix;
} FinalAtag_V_tuple;


double **w_matrix_generator(double **file_content_1, int n, int d);
double **Diagonal_Degree_Mat(double **W, int n );
void power_by_minus_half(double **D, int n);
double **Laplacian(double **D_halfed, double **W, int n);
FinalAtag_V_tuple convergence (double **A, int n);
ee_ve_tuple *sorting_ee_ve_tuples(double **final_A_tag, double **V, int n);
int determine_k (ee_ve_tuple *sorted_arr, int n);
void free_2D_double_array(double **arr, int len);





double **kmeanspp(double **file_content_2, double **centroids, int n, int d, int k, int max_iter);
int Python_flow_calculate_k (int external_k, int n, int d, double **file_content_1);
double** Python_flow_case_1 (int external_k, int n, int d, double **file_content_1);
void C_flow(int external_k, int goal, int n, int d, double **file_content_1);

#endif
