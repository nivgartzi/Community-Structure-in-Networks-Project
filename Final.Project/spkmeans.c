#define _CRT_SECURE_NO_WARNINGS

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define general_err "An Error Has Occured"
#define input_err "Invalid Input!"

/* checks if: -0.00005 < number < 0. if true -> returns 0, else -> returns val */
double fix_number (double number){
    if(number < 0 && number > -0.00005 ){
        return 0;
    }
    return number;
}

/* prints the items on the (main) diagonal */
void print_diagonal (double **matrix, int n)
{
    int i;
    double val=0;
    for(i=0;i<n;i++)
    {
        val = fix_number(matrix[i][i]);
        if(i!=n-1) {
            printf("%.4f,", val);
        }
        else {
            printf("%.4f", val);
        }
    }
}

/* prints matrix transpose */
void print_matrix_transpose(double **A, int n){
    int i,j;
    double val=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            val = fix_number(A[j][i]);
            printf("%.4f",val);
            if(j<n-1)
            {
                printf(",");
            }
        }
        if(i!=n-1){
            printf("\n");
        }
    }
}


/* prints nXm matrix */
void print_matrix(double **A, int n,int m){
    int i, j;
    double val = 0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            val = fix_number(A[i][j]);
            if(j<m-1)
            {
                printf("%.4f,", val);
            }

            else
            {
                printf("%.4f", val);
            }
        }
        if(i!=n-1){
            printf("\n");
        }
    }
}


/* generates and returns nXn matrix, initialized to c's */
double **init_nXn_continuous_matrix(int n, int c){
    double *p, **a;
    int i,j;
    a= (double**) calloc(n, sizeof (double *)); 
    assert(a!=NULL && general_err); 

    for(i=0; i<n; i++){

        p = (double*) calloc(n, sizeof (double)); 
        assert(p!=NULL && general_err); 
        a[i] = p;

        for(j=0; j<n; j++){
            a[i][j] = c;
        }
    }
    return a;
}

    
/* generates and returns nXm(!) matrix, initialized to 0's */
double **init_nXm_continuous_matrix(int n, int m){
    double *p, **a;
    int i,j;
    a=(double**) calloc(n, sizeof (double *)); 
    assert(a!=NULL && general_err); 

    for(i=0; i<n; i++){

        p =(double*) calloc(m, sizeof (double)); 
        assert(p!=NULL && general_err); 
        a[i] = p;

        for(j=0; j<m; j++){
            a[i][j] = 0;
        }
    }
    return a;
}


/* genrates and returns nXn Identity matrix */
double **init_eye_matrix(int n) {
    double **mat = init_nXn_continuous_matrix(n,0);
    int i=0;
    while(i<n)
    {
        mat[i][i]=1;
        i++;
    }
    return mat;
}

/* calculates the "connection's weight" between vi and vj, assuming i != j */
double weight_calculator(double* v1, double* v2, int d){
    int i;
    double sum = 0.0, curr, result, norm;
    for(i=0; i<d; i++){
        curr = pow(v1[i] - v2[i], 2);
        sum += curr;
    }
    norm = sqrt(sum);
    result = exp(-norm / 2);
    return result;
}


/* generates nXn "W" matrix */
double **w_matrix_generator(double **file_content_1, int n, int d){
    double **w_matrix;
    int i, j;
    double weight;
    w_matrix = init_nXn_continuous_matrix(n,0);

    assert(w_matrix!=NULL && general_err);

    for(i=0; i<n; i++)
    {
        for (j=0; j<i; j++)
        {
            weight = weight_calculator(file_content_1[i], file_content_1[j], d);
            w_matrix[i][j] = weight;
            w_matrix[j][i] = weight; /* from symetry */

        }
    }
    return w_matrix;
}

/* Faciliating function; sums up the items of the i_th row, at the 'W' matrix */
double Summary_Of_Row(double **W,int i, int n)
{
    double sum=0;
    int z;
    for(z=0;z<n;z++)
    {
        sum=sum+W[i][z];
    }
    return sum;
}

/* generates and returns nXn 'D' matrix */
double **Diagonal_Degree_Mat(double **W, int n ){

    int i=0;
    double d_ii;
    double **D= init_nXn_continuous_matrix(n,0);
    assert(D!=NULL && general_err);

    while(i<n)
    {
        d_ii = Summary_Of_Row(W,i,n);
        D[i][i]=d_ii;
        i++;
    }
    return D;
}

/* transforms 'D' matrix in-place to nXn 'D^-0.5' matrix */
void power_by_minus_half(double **D, int n) {
    int i=0;
    double base, result;
    while(i<n)
    {
        base = D[i][i];
        result= pow(base,(-0.5));
        D[i][i]=result;
        i++;
    }
}

/* transforms the given nXn matrix in-place to nXn I-A matrix */
void reduce_from_eye(double **A, int n){
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(i==j)
            {
                A[i][i]=1-A[i][i];
            }
            else
            {
                A[i][j]=0-A[i][j];
            }
        }
    }
    while(i<n)
    {

        i++;
    }
}

/* generates and returns nXn "L_norm" matrix */
double **Laplacian(double **D_halfed, double **W, int n){
    int i,j;
    double scalar_1 , scalar_2;
    double **L_norm = init_nXn_continuous_matrix(n,0);
    assert(L_norm!=NULL && general_err);

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            scalar_1=D_halfed[i][i];
            scalar_2=D_halfed[j][j];
            L_norm[i][j]=scalar_1*scalar_2*W[i][j];
        }
    }
    reduce_from_eye(L_norm,n);
    return L_norm;
}

/* indicates row and column of biggest absolute value, who's off diagnoal */
typedef struct {
    int row;
    int col;
} i_j_tuple;

/* returns row ('i') and column ('j') of biggest absolute value, who's off diagnoal */
i_j_tuple find_offdiagonal_indices(double** A, int n) {
    int row = 0;
    int col = 1;
    double pivot = fabs( A[0][1] );
    i_j_tuple indices;
    int i,j;
    double curr;

    for(j=2; j<n; j++)
    {
        for(i=0; i<j; i++) /* A and A' are symetric, thus we can go above the diagnoal only */
        {
            curr = fabs(A[i][j]);
            if(curr>pivot)
            {
                pivot = curr;
                row = i;
                col = j;
            }
        }
    }
    indices.row = row;
    indices.col = col;
    return indices;
}

/* calculates and returns theta */
double theta_calculator(double **A,int i, int j){
    double theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    return theta;
}

/* returns sign of a given double */
int sign(double num){
    if(num < 0.0)
    {
        return -1;
    }
    return 1;
}

/* calculates and returns t */
double t_calculator(double theta){
    double mark = sign(theta);
    double denominator = fabs(theta) + sqrt(pow(theta,2) + 1);
    double t = mark / denominator;
    return t;
}

/* calculates and returns c */
double c_calculator(double t){
    double denominator =  sqrt(pow(t,2) + 1);
    double c = 1 / denominator;
    return c;
}

/* calculates and returns s */
double s_calculator(double t, double c){
    return t*c;
}

/* returns sum of the squares of the items which are off the diagonal */
double off_diagonal_pow (double **mat, int n)
{
    int i,j;
    double sum=0;
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if(i != j)
            {
                sum += pow(mat[i][j],2);
            }
        }
    }
    return sum;
}


/* changing A to be A'(in-place) */
void modify_to_A_tag_matrix(double **A, double c, double s, int i, int j, int n){

    double a_ii = A[i][i];
    double a_jj = A[j][j];
    double a_ij = A[i][j];
    double a_ri, a_rj;
    int r;
    A[i][i] = pow(c,2) * a_ii + pow(s,2)*a_jj -2*s*c*a_ij;
    A[j][j] = pow(s,2) * a_ii + pow(c,2)*a_jj +2*s*c*a_ij;
    A[i][j]=0;
    A[j][i]= A[i][j]; /* A' symetric */
    for(r=0;r<n;r++)
    {
        if(r!=i && r!=j)
        {
            a_ri = A[r][i];
            a_rj = A[r][j];

            A[r][i]  = c*a_ri - s*a_rj;
            A[r][j]  = c*a_rj + s*a_ri;

            /* from symetry: */
            A[i][r] = A[r][i];
            A[j][r] = A[r][j];
        }
    }
}

typedef struct {
    double **final_A_tag;
    double **V_matrix;
} FinalAtag_V_tuple;


void free_2D_double_array(double **arr, int len){
    int i;
    for (i=0;i<len;i++){
        free(arr[i]);
    }
    free(arr);
} /* HW1 function */ 


void free_finalAtag_V_tuple(FinalAtag_V_tuple tuple, int len){
    free_2D_double_array(tuple.V_matrix,len);
}


FinalAtag_V_tuple convergence (double **A, int n){
    const double epsilon = 1.0e-15;
    double t,c,s;
    i_j_tuple current_indices;
    int i,j;
    int counter=1;
    int runner;
    double temp;
    FinalAtag_V_tuple output;
    double condition;
    output.V_matrix = init_eye_matrix(n);

    current_indices = find_offdiagonal_indices(A,n);
    i = current_indices.row;
    j= current_indices.col;
    t = t_calculator( theta_calculator(A,i,j) );
    c = c_calculator(t);
    s = s_calculator(t,c);

    for(runner=0; runner<n; runner++)
    {
        temp = output.V_matrix[runner][i];
        output.V_matrix[runner][i] = temp * c - output.V_matrix[runner][j] * s;
        output.V_matrix[runner][j] = temp * s + output.V_matrix[runner][j] * c;
    }

    /* condition = 2*fabs(pow(A[i][j],2)); */
    condition = 2*(pow(A[i][j],2));

    while(  condition > epsilon && counter<100   )
    {
        modify_to_A_tag_matrix(A,c,s,i,j,n);
        counter++;
        current_indices = find_offdiagonal_indices(A,n);
        i = current_indices.row;
        j= current_indices.col;
        t = t_calculator( theta_calculator(A,i,j) );
        c = c_calculator(t);
        s = s_calculator(t,c);

        for(runner=0; runner<n; runner++)
        {
        temp = output.V_matrix[runner][i];
        output.V_matrix[runner][i] = temp * c + output.V_matrix[runner][j] * (-s);
        output.V_matrix[runner][j] = temp * s + output.V_matrix[runner][j] * c;
        }

        condition  = fabs(2 * pow(A[i][j],2)) ;
    }
    output.final_A_tag=A;
    return output;
}


typedef struct {
    double ee;  /*  eigenvalue */
    int index;  /* the ee's (diagonal wise) index in final A' (used for stable sorting, in case when ee's are equal) */
    double *ve; /* eigenvector */
} ee_ve_tuple;

void free_ee_ve_tuple_array(ee_ve_tuple *tuple, int len){
    int i;
    for(i=0; i<len; i++){
        free(tuple[i].ve);
    }
    free(tuple);
}

void init_ee_ve_tuples(ee_ve_tuple *arr, int n) {
    int i;
    for(i=0;i<n;i++)
    {
        arr[i].ee=(-1);
        arr[i].index=(-1);
        arr[i].ve=(double*) calloc(n,sizeof(double));
        assert(arr[i].ve!=NULL && general_err); 
    }
}

void assign_ve_for_ee_ve_tuple(double *ve_base, int fixed_col, double **V, int n) {
    int i;
    for(i=0;i<n;i++)
    {
        ve_base[i] = V[i][fixed_col];
    }
}

/* (after the init) the func will set the self values and correnspondent self vectors (the i_th column of V will be the i_th item's field) */
void set_ee_ve_array(ee_ve_tuple *arr, double **final_A_tag, double **V, int n)
{
    int j, arr_index, fixed_col;
    double *ve_base;

    for(j=0;j<n;j++)
    {
        arr_index=j;
        arr[arr_index].ee = final_A_tag[j][j];
        arr[arr_index].index = j;
        ve_base = arr[arr_index].ve;
        fixed_col=j;
        assign_ve_for_ee_ve_tuple (ve_base, fixed_col, V, n);
    }
}

/* fills first n items starting from p - into vec (the so called array..) */
void swap_vectors_contents(double *vec1, double *vec2, int n) {
    int i;
    double temp;
    for(i=0;i<n;i++)
    {
        temp = vec1[i];
        vec1[i] = vec2[i];
        vec2[i] = temp;
    }
}

/* assuming: i1 < i2 */
void swap_ee_ve_tuples(ee_ve_tuple *arr, int i1, int i2, int n){
    double ee_temp;
    int index_temp;
    double *base_of_ve1, *base_of_ve2;

    ee_temp=arr[i1].ee;
    arr[i1].ee=arr[i2].ee;
    arr[i2].ee=ee_temp;

    index_temp=arr[i1].index;
    arr[i1].index=arr[i2].index;
    arr[i2].index=index_temp;

    base_of_ve1 =  arr[i1].ve;
    base_of_ve2 =  arr[i2].ve;
    swap_vectors_contents(base_of_ve1, base_of_ve2, n);
}

/* returns array of sorted ee ve tuples */
ee_ve_tuple *sorting_ee_ve_tuples(double **final_A_tag, double **V, int n) {

    double current_ee, next_ee;
    int current_ve_col_index, next_ve_col_index, i, j;
    ee_ve_tuple *arr;
    arr = (ee_ve_tuple*)calloc(n,sizeof(ee_ve_tuple));
    assert(arr!=NULL && general_err);
    init_ee_ve_tuples(arr,n);
    set_ee_ve_array(arr,final_A_tag,V,n);

    /* let's sort arr: */
    for(i=0;i<n-1;i++)
    {
        for(j=0;j<n-i-1;j++)
        {
            current_ee=arr[j].ee;
            next_ee=arr[j+1].ee;
            if(current_ee==next_ee) /* for stable sorting */
            {
                current_ve_col_index=arr[j].index;
                next_ve_col_index=arr[j+1].index;
                if(current_ve_col_index>next_ve_col_index)
                {
                    swap_ee_ve_tuples(arr,j,(j+1),n);
                }
            }
            if(current_ee>next_ee)
            {
                swap_ee_ve_tuples(arr,j,(j+1),n);
            }
        }
    }
    return arr;
}


/* The Eigengap Heuristic (determing the cluster's number) */
int determine_k (ee_ve_tuple *sorted_arr, int n) {
    int place=0;
    int i;
    int N; /* limit */
    double current_ee, next_ee, current_gap, max_gap, first_gap;

    if(n%2==0)
    {
        N=(n/2);
    }
    else
    {
        N=(n-1)/2;
    }

    first_gap = fabs(sorted_arr[0].ee - sorted_arr[1].ee);
    max_gap = first_gap;

    for(i=0;i<N;i++)
    {
        current_ee=sorted_arr[i].ee;
        next_ee=sorted_arr[i+1].ee;
        current_gap = fabs(current_ee - next_ee);

        if (current_gap>max_gap)   /*  "in case of equality - lowest index"; therefore: ' > ' ( and not ' >= ' )  */
        {
            max_gap = current_gap;
            place = i;
        }
    }
    return (place+1);
}

/* first k eigenvector's will be U's columns */
double **generate_U (ee_ve_tuple *sorted_arr, int n, int k) {
    double **U = init_nXm_continuous_matrix(n,k);
    int fixed_col, row;

    for(fixed_col=0;fixed_col<k;fixed_col++) /*   0 <=  fixed_col  < k ( going throught first k vectors )  */
    {
        for(row=0;row<n;row++)
        {
            U[row][fixed_col] = sorted_arr[fixed_col].ve[row];
        }
    }
    return U;
}

/* calculates sum of squared values in row i, assuming k == num of columns */
double Summary_Of_Row_Pow(double **matrix,int i, int k)
{
    double sum=0;
    int j;
    for(j=0;j<k;j++)
    {
        sum += pow(matrix[i][j],2);
    }
    return sum;
}

/* returns nXk matrix "T", based on matrix U */
double **form_T (double **U, int n, int k)
{
    double **T = init_nXm_continuous_matrix(n,k);
    int i,j;
    double denominator, mark;

    for(i=0; i<n; i++){
        denominator = sqrt(Summary_Of_Row_Pow(U,i,k)); /* one denominator for a whole row */
        for(j=0; j<k; j++){
            mark = U[i][j];
            T[i][j] = mark / denominator;
        }
    }
    return T;
}


/** HW1 **/
/*/////////////////////////////////// */
typedef struct
{
    int size;
    double* centroid;
    double* total_summary;
} cluster;

void init_clusters (cluster *clusters_arr, double **file_content_2, int k, int d)
{
    int i,j;
    for (i=0; i<k; i++)
    {
        cluster clt;

        clt.size=0;

        /* allocate memory */ 
        clt.centroid = (double*) calloc(d,sizeof(double));
        assert(clt.centroid!=NULL && general_err);

        clt.total_summary = (double*)calloc(d,sizeof(double)); /* vector of zeros, because  folder is empty at start*/
        assert(clt.total_summary!=NULL && general_err);

        for(j=0 ; j<d; j++)
        {
            clt.centroid[j] = (double)file_content_2[i][j];
            clt.total_summary[j]  = (double)0;

        }
        clusters_arr[i] = clt;
    }

}


int update_clusters_centroids_and_tell(cluster *clusters_arr, int d,int k)
{
    int i,c;
    double temp;
    int status = 0;
    for(i=0;i<k;i++)
    {
        for(c=0;c<d;c++)
        {
            temp = clusters_arr[i].centroid[c];
            clusters_arr[i].centroid[c] = (double)((double)clusters_arr[i].total_summary[c]/clusters_arr[i].size);
            if(temp!=clusters_arr[i].centroid[c]){
                status=1;
            }
        }
    }
    return status;
}

void reset_clusters (cluster *clusters_arr, int k, int d)
{
    int i=0, c=0;
    for(i=0; i<k;i++)
    {
        clusters_arr[i].size=0;
        for(c=0; c<d; c++)
        {
            double reset = 0.0;
            clusters_arr[i].total_summary[c]=reset;
        }
    }
}

double distance_calculator(double *v1, double *v2, int d)
{
    double distance = (double)0;
    int c;
    double coordinate;
    for(c=0;c<d;c++)
    {
        coordinate = (double) ((double)(v1[c] - v2[c]))*((double)(v1[c] - v2[c]));
        distance += coordinate;
    }
    return distance;
}


void print_centroids(cluster *clusters_arr, int k, int d){
    int i;
    int c;
    double val=0;
    for (i = 0; i < k; i++){
        for (c = 0; c < d; c++)
        {
            val = fix_number(clusters_arr[i].centroid[c]);
            printf("%.4f", val);

            if(c<d-1)
            {
                printf(",");
            }
        }

        if(i != k-1){
            printf("\n");
        }
    }
}


void free_clusters_array(cluster *clusters_arr, int k){
    int i;
    for(i=0;i<k;i++){
        free(clusters_arr[i].centroid);
        free(clusters_arr[i].total_summary);
    }
    free(clusters_arr); 
}

int is_int(char* num)
{
    int j, i;
    j = strlen(num);
    for (i = 0; i < j; i++){
        if (!isdigit(num[i]))
            return 0;
    }
    return 1;
}

/*/////////////////////////////////// */



/** HW2 **/
/*/////////////////////////////////// */
void init_clusters_assignment2 (cluster *clusters_arr, double **centroids, int k, int d)
{
    int i,j;
    for (i=0; i<k; i++)
    {
        cluster clt;

        clt.size=0;

        /* allocate memory */ 
        clt.centroid = (double*) calloc(d,sizeof(double));
        assert(clt.centroid!=NULL && general_err);

        clt.total_summary = (double*)calloc(d,sizeof(double));
        assert(clt.total_summary!=NULL && general_err);

        for(j=0 ; j<d; j++)
        {
            clt.centroid[j] = centroids[i][j];
            clt.total_summary[j]  = 0;

        }
        clusters_arr[i] = clt;
    }

}

double **kmeanspp(double **file_content_2, double **centroids, int n, int d, int k, int max_iter)
{

    int indicator=1;    /*1 means no convergence*/
    int counter=0;
    int folder_index;
    int i,j,c;
    double min_distance;
    double current_distance;
    cluster *clusters_arr;

    clusters_arr = (cluster*)calloc(k,sizeof(cluster));
    assert(clusters_arr!=NULL && general_err);

    init_clusters_assignment2(clusters_arr,centroids,k,d);

    while(counter<max_iter && indicator==1)  /* loop until convergence of all clusters */
    {
        indicator=0;
        for(i=0; i<n; i++)  /*for each xi*/
        {
            min_distance = distance_calculator(file_content_2[i],clusters_arr[0].centroid,d);
            folder_index=0;

            for(j=0; j<k; j++) /* for each centroid */
            {
                current_distance = distance_calculator(file_content_2[i],clusters_arr[j].centroid,d);
                if(current_distance<=min_distance)
                {
                    min_distance=current_distance;
                    folder_index=j;

                }
            }
            /* update fields of clusters_arr[folder_index]: */
            clusters_arr[folder_index].size += 1 ;
            for(c=0; c<d; c++)
            {
                clusters_arr[folder_index].total_summary[c] += file_content_2[i][c];
            }
        }
        counter++;
        indicator = update_clusters_centroids_and_tell(clusters_arr,d,k);

        reset_clusters(clusters_arr,k,d);
    }

    for(i=0;i<k;i++)
    {
        for(c=0;c<d;c++)
        {
            centroids[i][c]=clusters_arr[i].centroid[c];
        }
    }
    return centroids;
}
/*/////////////////////////////////// */


/** for goal=spk: returns T **/
double** Python_flow_case_1 (int external_k, int n, int d, double **file_content_1) {

    int k;
    double **W, **D, **L_norm, **U, **T;
    double **final_A_tag, **V;
    ee_ve_tuple *sorted_arr;
    FinalAtag_V_tuple important_tuple;
    W = w_matrix_generator(file_content_1,n,d);
    D= Diagonal_Degree_Mat(W,n);

    power_by_minus_half(D,n);
    /* D is now D^-(0.5) */

    L_norm= Laplacian(D,W,n);

    important_tuple = convergence(L_norm,n);
    final_A_tag = important_tuple.final_A_tag;
    V = important_tuple.V_matrix;

    sorted_arr = sorting_ee_ve_tuples(final_A_tag,V,n);

    if(external_k==0)
    {
        k = determine_k(sorted_arr,n);
    }
    else
    {
        k = external_k;
    }

    U = generate_U(sorted_arr,n,k);

    T = form_T(U,n,k);

    /* free all dynamic memory */ 
    free_2D_double_array(W, n);
    free_2D_double_array(D, n); 
    free_2D_double_array(L_norm, n); 
    free_ee_ve_tuple_array(sorted_arr,n); 
    free_2D_double_array(U,n);
    free_finalAtag_V_tuple(important_tuple,n);

    return T;
}

/** for goal=spk: returns k **/
int Python_flow_calculate_k (int external_k, int n, int d, double **file_content_1)
{
    int k;
    double **W, **D, **L_norm;
    double **final_A_tag, **V;
    ee_ve_tuple *sorted_arr;
    FinalAtag_V_tuple important_tuple;

    W = w_matrix_generator(file_content_1,n,d);  /* goal = wam */
    D= Diagonal_Degree_Mat(W,n);  /* goal = ddg */

    power_by_minus_half(D,n);
    /* D is now D^-(0.5) */

    L_norm= Laplacian(D,W,n);  /* goal = lnorm */

    important_tuple = convergence(L_norm,n);   /* contains: final_A_tag and V */
    final_A_tag = important_tuple.final_A_tag;
    V= important_tuple.V_matrix;

    sorted_arr = sorting_ee_ve_tuples(final_A_tag,V,n);


    if(external_k==0)
    {
        k = determine_k(sorted_arr,n);
    }
    else
    {
        k = external_k;
    }

    /* free all dynamic memory */ 

    free_2D_double_array(W, n); 
    free_2D_double_array(D, n); 
    free_2D_double_array(L_norm, n); 
    free_ee_ve_tuple_array(sorted_arr,n);  
    free_finalAtag_V_tuple(important_tuple,n);

    return k;
}


/** the end of the following function: implementation of k_means **/
void C_flow(int external_k, int goal, int n, int d, double **file_content_1) {

    int max_iter=300; /* given */
    double **file_content_2, *plained_arr_2;
    double **W, **D, **L_norm, **U, **T;
    double **final_A_tag, **V;
    ee_ve_tuple *sorted_arr;
    FinalAtag_V_tuple important_tuple;

    /* declreations from hw1: */

    int indicator=1;    /*1 means no convergence*/
    int counter=0;
    int folder_index;
    int i,j,c, k;
    double min_distance;
    double current_distance;
    cluster *clusters_arr;

    W = w_matrix_generator(file_content_1,n,d);  /* goal = wam */
    if(goal==2)
    {
        print_matrix(W,n,n); 

        free_2D_double_array(W, n);
        return;
    }

    D= Diagonal_Degree_Mat(W,n);  /* goal = ddg */
    if(goal==3)
    {
        print_matrix(D,n,n);

        free_2D_double_array(W, n);
        free_2D_double_array(D, n);
        return;
    }

    power_by_minus_half(D,n);
    /* D = D^-(0.5) */

    L_norm= Laplacian(D,W,n);  /* goal = lnorm */
    if(goal==4)
    {
        print_matrix(L_norm,n,n);

        free_2D_double_array(W, n);
        free_2D_double_array(D, n);
        free_2D_double_array(L_norm,n);
        return;
    }

    important_tuple = convergence(L_norm,n);   
    final_A_tag = important_tuple.final_A_tag;
    V= important_tuple.V_matrix;

    if(goal==5)
    {
        print_diagonal(final_A_tag,n); /*  goal=jacobi, part a | print's the self values */
        printf("\n");
        print_matrix_transpose(V,n);   /* goal=jacobi, part b | print's the self vectors (V columns) in lines */

        free_2D_double_array(W, n);
        free_2D_double_array(D, n);
        free_2D_double_array(L_norm,n);
        free_finalAtag_V_tuple(important_tuple,n);
        return;
    }

    /* if got this far, it means that goal = spk  */
    sorted_arr = sorting_ee_ve_tuples(final_A_tag,V,n);

    if(external_k==0)
    {
        k = determine_k(sorted_arr,n);
    }
    else
    {
        k = external_k;
    }

    U = generate_U(sorted_arr,n,k);
    T = form_T(U,n,k);

    /*
       treating each row of T as a vector (at R^k), and we cluster them into k clusters via the K-means algorithm.
       e.g: the argument d (the following 4_th argument, which is the dimension of the vectors) is k !
     */
    d=k;

    file_content_2 = (double**)malloc(n*sizeof(double*));
    assert(file_content_2!=NULL && general_err);

    plained_arr_2 = (double*)malloc(n*d*sizeof(double));
    assert(plained_arr_2!=NULL && general_err);

    for(i=0;i<n;i++)
    {
        for(j=0;j<d;j++)
        {
            plained_arr_2[ (d*i + j)  ] = T[i][j];
        }
    }

    for(i=0; i<n; i++)
    {
        file_content_2[i] = plained_arr_2 + i*d;
    }

    clusters_arr =  (cluster*)malloc(k*sizeof(cluster));
    assert(clusters_arr!=NULL && general_err);

    init_clusters(clusters_arr, file_content_2, k, d);  /* note: the dimension ('d') is k */

    while(counter<max_iter && indicator==1)  /* loop until convergence of all clusters */
    {
        indicator=0;
        for(i=0; i<n; i++)  /*for each xi*/
        {
            min_distance = distance_calculator(file_content_2[i],clusters_arr[0].centroid,d);
            folder_index=0;

            for(j=0; j<k; j++) /* for each centroid */
            {
                current_distance = distance_calculator(file_content_2[i],clusters_arr[j].centroid,d);
                if(current_distance<=min_distance)
                {
                    min_distance=current_distance;
                    folder_index=j;

                }
            }

            /*update the fields of clusters_arr[folder_index]: */
            clusters_arr[folder_index].size += 1 ;
            for(c=0; c<d; c++)
            {
                clusters_arr[folder_index].total_summary[c] += file_content_2[i][c];
            }
        }
        counter++;
        indicator = update_clusters_centroids_and_tell(clusters_arr,d,k);

        reset_clusters(clusters_arr,k,d);
    }

    print_centroids(clusters_arr,k,d); /* goal = spk printing */

    /* free all dynamic memory */ 
    free(plained_arr_2); 
    free(file_content_2);
    free_2D_double_array(W, n); 
    free_2D_double_array(D, n); 
    free_2D_double_array(L_norm,n);
    free_clusters_array(clusters_arr,k);
    free_finalAtag_V_tuple(important_tuple,n);
    free_ee_ve_tuple_array(sorted_arr,n); 
    free_2D_double_array(U,n); 
    free_2D_double_array(T,n); 

}      /*   1 < = goal  < = 5  (1-->spk,....., 5-->jacobi) */


int determine_goal(int enum_as_int)
{
    if (enum_as_int == 119) /* wam */
        return 2;
    if(enum_as_int==100) /* ddg */
        return 3;
    if(enum_as_int==108) /* lnorm */
        return 4;
    if(enum_as_int==106) /* jacobi */
        return 5;
    if(enum_as_int==115) /* spk */
        return 1;

    return 0;
}

int main (int argc, char *argv[])
{
    int goal, k;
    int enum_as_int = (int) argv[2][0];
    char c;
    int n = 0;
    int d = 0; /* the dimesion of the given external set of n points  */
    FILE *data;
    double **file_content_1;
    int commas = 0;  
    int n_minus_1 = 0;

    int counter, indicator, i=0, j=0;
    double number, first_number;

    assert(argc != 3 && input_err);

    /* setting k and goal */
    k = atoi(argv[1]);
    goal = determine_goal(enum_as_int);

    (void) kmeanspp; 

    /*  now we will determine n, and then put 'file_name' data into file_content_1 */
    data = fopen(argv[3], "r");
    assert(data);

    /* determining n and d */
    while ( (c = getc(data)) != '\n')
    {
        if (c == ',')
        {
            commas++;
        }
    }
    d = commas + 1;
    fseek(data, 0L, SEEK_SET);

    while  ( (c = getc(data)) != EOF)
    {
        if (c == '\n')
        {
            n_minus_1++;
        }
    }
    n = n_minus_1 + 1;
    fseek(data, 0L, SEEK_SET);

    /* put 'file_name' data into file_content_1 */
   file_content_1 = init_nXm_continuous_matrix(n,d);

    while ( ! feof(data) )
    {
        counter = 0, number = 0, indicator = 0;

        if (i == n)
        {
            break;
        }
        while (c != '.')
        {
            c = getc(data);

            if (c == '.')
            {
                break;
            }
            if (c == '-')
            {
                indicator = 1;
                continue;
            }
            first_number = c - '0';
            number = 10 * number + first_number;
        }
        c = getc(data);
        while ( c != ',' && c != '\n' && !feof(data) )
        {
            counter++;
            first_number = c - '0';
            number = number + pow(0.1, counter) * first_number;
            c = getc(data);
        }
        if (indicator == 1)
        {
            number = number * (-1);
        }

        file_content_1[i][j] = number;

        if (j == (d - 1))
        {
            if (c == '\n')
            {
                ++i;
                j = 0;
            }
        }
        else
        {
            if (c != '\n')
            {
                ++j;
            }
        }
    }
    fclose(data); 

    C_flow(k,goal,n,d,file_content_1);
 
    free_2D_double_array(file_content_1,n);

    return 0;
}









