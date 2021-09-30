#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
#include <assert.h>
#include <string.h>


#define general_err "An Error Has Occured"
#define input_err "Invalid Input!"

/* receives file_content_1, n, d_of_file_content_1, external_k and returns k*/
static PyObject* getk(PyObject *self, PyObject *args) 
{
    int external_k, n, d_of_file_content_1;
    int calculated_k;
    int i,j;
    double **file_content_1;
    int casting_failed = 0;
    PyObject *pythonVectors, *vector, *coordinate; 

    /* validate we got all parameters */
    if (!PyArg_ParseTuple(args, "Oiii", &pythonVectors, &n, &d_of_file_content_1, &external_k))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    /* allocate memory */
    file_content_1 = (double **)calloc(n, sizeof(double *));
    assert(file_content_1 != NULL && general_err);
    for (i = 0; i < n; i++){
        file_content_1[i] = (double *)calloc(d_of_file_content_1, sizeof(double));
        assert(file_content_1[i] != NULL && general_err);
    }

    /* file_content_1 will now be fed by what's inside &pythonVectors */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(pythonVectors, i);
        if (!PyList_Check(vector)){
            continue;
        }
        for (j = 0; j < d_of_file_content_1; j++) {
            coordinate = PyList_GetItem(vector, j);
            file_content_1[i][j] = PyFloat_AsDouble(coordinate); /** convert a Python float to C double */
            /* check casting from python float to C double **/
            if (PyErr_Occurred() && file_content_1[i][j]  == -1.0){
                casting_failed = 1;
                break;
            }
        }
        if (casting_failed == 1){
            break;
        }
    }
    if (casting_failed == 0)
    {
        calculated_k = Python_flow_calculate_k(external_k,n,d_of_file_content_1,file_content_1);
    }
    else
    {
        calculated_k = -1;
    }

    /* free allocated memory */ 
    free_2D_double_array(file_content_1,n);

    return Py_BuildValue("i", calculated_k); 
}

/*receive file_content_1, n, d_of_file_content_1, goal_id and, print matrix and return goal id */
static PyObject* other_cases(PyObject *self, PyObject *args) 
{
    int  n, d_of_file_content_1;
    int goal_id; /*between 2 to 5 */
    int i,j;
    double **file_content_1;
    int casting_failed = 0;
    PyObject *pythonVectors, *vector, *coordinate;

    /* validate we got all parameters */

    if (!PyArg_ParseTuple(args, "Oiii", &pythonVectors, &n, &d_of_file_content_1, &goal_id))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    /* allocate memory */
    file_content_1 = (double **)calloc(n, sizeof(double *));
    assert(file_content_1 != NULL && general_err);

    for (i = 0; i < n; i++){
        file_content_1[i] = (double *)calloc(d_of_file_content_1, sizeof(double));
        assert(file_content_1[i] != NULL && general_err);
    }

    /* file_content_1 will now be fed by what's inside &pythonVectors */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(pythonVectors, i);
        if (!PyList_Check(vector)){
            continue;
        }
        for (j = 0; j < d_of_file_content_1; j++) {
            coordinate = PyList_GetItem(vector, j);
            file_content_1[i][j] = PyFloat_AsDouble(coordinate); /** convert a Python float to C double */
            /* check casting from python float to C double **/
            if (PyErr_Occurred() && file_content_1[i][j]  == -1.0){
                casting_failed = 1;
                break;
            }
        }
        if (casting_failed == 1){
            break;
        }
    }
    if (casting_failed == 0)
    {
        C_flow(-1,goal_id,n,d_of_file_content_1,file_content_1);
    }

    /* free allocated memory */ 
    free_2D_double_array(file_content_1,n);

    return Py_BuildValue("i", goal_id);
}

/* receive file_content_1, n, d_of_file_content_1, external_k, calculated_k */
static PyObject* get_file_content_2(PyObject *self, PyObject *args) 
{
    int external_k, n, d_of_file_content_1;
    int calculated_k;
    int i,j;
    int d;
    double **file_content_1, **file_content_2;
    int casting_failed = 0;
    PyObject *python_file_content_2, *pythonVectors, *vector, *coordinate, *pythonVector;

    /* validate we got all parameters */
    if (!PyArg_ParseTuple(args, "Oiiii", &pythonVectors, &n, &d_of_file_content_1, &external_k,&calculated_k))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    /* allocate memory */
    file_content_1 = (double **)calloc(n, sizeof(double *));
    assert(file_content_1 != NULL && general_err);
    for (i = 0; i < n; i++){
        file_content_1[i] = (double *)calloc(d_of_file_content_1, sizeof(double));
        assert(file_content_1[i] != NULL && general_err);
    }

    /* file_content_1 will now be fed by what's inside &pythonVectors */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(pythonVectors, i);
        if (!PyList_Check(vector)){
            continue;
        }
        for (j = 0; j < d_of_file_content_1; j++) {
            coordinate = PyList_GetItem(vector, j);
            file_content_1[i][j] = PyFloat_AsDouble(coordinate); /** convert a Python float to C double **/
            /* check casting from python float to C double */
            if (PyErr_Occurred() && file_content_1[i][j]  == -1.0){
                casting_failed = 1;
                break;
            }
        }
        if (casting_failed == 1){
            break;
        }
    }
    file_content_2=NULL;
    if (casting_failed == 0)
    {
        file_content_2 = Python_flow_case_1(external_k,n,d_of_file_content_1,file_content_1);
    }
    python_file_content_2 = PyList_New(n);
    if (python_file_content_2 == NULL)
        return NULL;

    /* matrix T is nXk; i.e: "d" is k */
    d=calculated_k;

    for (i = 0; i < n; i++)
    {
        pythonVector = PyList_New(d);
        if (pythonVector == NULL)
            return NULL;
        for (j = 0; j < d; j++)
        {
            PyList_SET_ITEM(pythonVector, j, Py_BuildValue("d", file_content_2[i][j]));
        }
        PyList_SetItem(python_file_content_2, i, Py_BuildValue("O", pythonVector));
    }

    /* free allocated memory */ 
    free_2D_double_array(file_content_1,n);
    free_2D_double_array(file_content_2,n);   

    return python_file_content_2;
}


static PyObject* get_final_clusters(PyObject *self, PyObject *args)
{
    int k, n, d, max_iter, i, j;
    double **vectors, **centroids, **final_centroids;
    int casting_failed = 0;
    PyObject *pythonCentroids, *pythonVectors, *vector, *centroid, *coordinate, *pythonFinalCentroids, *pythonCentroid;

    /* validate we got all parameters */
    if (!PyArg_ParseTuple(args, "OOiiii", &pythonCentroids, &pythonVectors, &max_iter, &n, &d, &k))
        return NULL;
    if (!PyList_Check(pythonCentroids))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    /* initialize centoirds, vectors arrays */
    vectors = (double **)calloc(n, sizeof(double *));
    assert(vectors != NULL && general_err);

    centroids = (double **)calloc(k, sizeof(double *));
    assert(centroids != NULL && general_err);

    for (i = 0; i < n; i++){
        vectors[i] = (double *)calloc(d, sizeof(double));
        assert(vectors[i] != NULL && general_err);
    }
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(d, sizeof(double));
        assert(centroids[i] != NULL && general_err);
    }

    /* insert values from python objects to C arrays */
    for (i = 0; i < n; i++) {
        vector = PyList_GetItem(pythonVectors, i);
        if (!PyList_Check(vector)){
            continue;
        }
        for (j = 0; j < d; j++) {
            coordinate = PyList_GetItem(vector, j);
            vectors[i][j] = PyFloat_AsDouble(coordinate); /* convert a Python float to C double */
            /* check casting from python float to C double */
            if (PyErr_Occurred() && vectors[i][j]  == -1.0){
                casting_failed = 1;
                break;
            }
        }
        if (casting_failed == 1){
            break;
        }
    }

    for (i = 0; i < k; i++) {
        centroid = PyList_GetItem(pythonCentroids, i);
        if (!PyList_Check(centroid)){
            continue;
        }
        for (j = 0; j < d; j++) {
            coordinate = PyList_GetItem(centroid, j);
            centroids[i][j] = PyFloat_AsDouble(coordinate); /* convert a Python float to C double */
            /* check casting from python float to C double */
            if (PyErr_Occurred() && centroids[i][j]  == -1.0){
                casting_failed = 1;
                break;
            }
        }
        if (casting_failed == 1){
            break;
        }
    }
    final_centroids = NULL;
    if (casting_failed == 0){
        final_centroids = kmeanspp(vectors, centroids, n, d, k, max_iter);
    }
    pythonFinalCentroids = PyList_New(k);
    if (pythonFinalCentroids == NULL)
        return NULL;

    for (i = 0; i < k; i++)
    {
        pythonCentroid = PyList_New(d);
        if (pythonCentroid == NULL)
            return NULL;
        for (j = 0; j < d; j++)
        {
            PyList_SET_ITEM(pythonCentroid, j, Py_BuildValue("d", final_centroids[i][j]));
        }
        PyList_SetItem(pythonFinalCentroids, i, Py_BuildValue("O", pythonCentroid));
    }

    /* free allocated memory */
    free_2D_double_array(vectors,n);   
    free_2D_double_array(centroids,k);  

    return pythonFinalCentroids;
}


static PyMethodDef kmeansppMethods[] = {
        {"get_final_clusters", (PyCFunction)get_final_clusters, METH_VARARGS, PyDoc_STR("will give the final centroids resulted from k-means++ algo")},
        {"getk", (PyCFunction)getk, METH_VARARGS, PyDoc_STR("will give k")},
        {"other_cases", (PyCFunction)other_cases, METH_VARARGS, PyDoc_STR("will print other matrices, just like C")},
        {"get_file_content_2", (PyCFunction)get_file_content_2, METH_VARARGS, PyDoc_STR("will give file_content_2")},
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT, "mykmeanssp", NULL, -1, kmeansppMethods
};


PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}