/*
 * ============================================================================
 *
 *                C extensions for Bayesian post-processing codes
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>
#include "six.h"

#include <stdlib.h>

#define MAX_IN_FILES 128

double findMaxLogL(FILE * input, double maxLogL);

typedef struct tagBurnInInput{
		char** files;
		int nfiles;
		int spin;
		float deltaLogL;
		int thin;
		char* output_file_name;
}BurnInInput;

typedef struct tagBurnInOutput{
		double** pos;
		int nSamples;
		double bayesFactorFromHarmonicMean;
}BurnInOutput;

int BurnIn(BurnInInput*,BurnInOutput*);

static PyObject* _burnin(PyObject *self, PyObject *args) {

    /***DECLARATION***/

    //INPUT
    
    PyObject* py_inputfile_list=NULL ;// List of files containing posterior chains."
    Py_ssize_t nfiles;
    PyObject* py_spin=NULL;

    //RETURN
    PyObject* py_pos_array=NULL;
    PyObject* py_pos_obj=NULL;
    PyObject* py_bayes=NULL;
    
    //BurnIn input
    BurnInInput* input=(BurnInInput*)malloc(sizeof(BurnInInput));

    int ndims;

    /***PARSE/PROCESS INPUT***/

    if (!PyArg_ParseTuple(args,"O!O!ds",&PyList_Type,&py_inputfile_list,&PyBool_Type,&py_spin,&input->deltaLogL,&input->output_file_name))  return NULL;

    input->nfiles=PyList_Size(py_inputfile_list);

    input->files = (char**)calloc(MAX_IN_FILES, sizeof(char*));

    Py_ssize_t i;
    for(i=0;i<input->nfiles;i++){
    	char* ptemp=NULL;
#if PY_MAJOR_VERSION >= 3
        if((ptemp=PyUnicode_AsUTF8(PyList_GetItem(py_inputfile_list,i)))!=0){
#else
        if((ptemp=PyString_AsString(PyList_GetItem(py_inputfile_list,i)))!=0){
#endif
        	input->files[i]=ptemp;
        }
        else {
        	nfiles--;
        }
    }

    Py_INCREF(Py_True);
    if(py_spin==Py_True) {
        input->spin=1;
        ndims = 9;
    } else {
        input->spin=0;
        ndims=15;
    }
    
    Py_DECREF(Py_True);

    /* Create output struct*/
    BurnInOutput* output=(BurnInOutput*)malloc(sizeof(BurnInOutput));

    /** BEGIN ROUTINE**/
    int burnin_exit;
    burnin_exit=BurnIn(input,output);

    py_bayes=PyFloat_FromDouble(output->bayesFactorFromHarmonicMean);

    npy_intp pos_dims[2] = {output->nSamples,ndims};

    py_pos_array=PyArray_SimpleNewFromData(2,pos_dims,NPY_DOUBLE,(void*)output->pos);
    py_pos_obj=PyArray_Return((PyArrayObject*)py_pos_array);

    Py_INCREF(py_pos_obj);

    return Py_BuildValue("(OO)",py_pos_obj,py_bayes);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
    {"_burnin", _burnin, METH_VARARGS,
    "This function 'burns in' MCMC chains."
    },
    {NULL,}
};


static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_bayespputils",
    "This module provides C extensions for Bayesian analysis and post-processing codes.",
    -1, methods, NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC PyInit__bayespputils(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit__bayespputils(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
