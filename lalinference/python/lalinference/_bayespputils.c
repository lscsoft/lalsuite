/*
 * ============================================================================
 *
 *                C extensions for Bayesian post-processing codes
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdlib.h>

#define MAX_IN_FILES 128

#define MODULE_NAME "_bayespputils"

/* doc string */
const char BUDocstring[] =
"This module provides C extensions for Bayesian analysis and post-\n"
"processing codes.";

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
    PyObject* py_deltaLogL=NULL;
    PyObject* py_outputfilename=NULL ;

    //RETURN
    PyObject* py_pos_array=NULL;
    PyObject* py_pos_obj=NULL;
    PyObject* py_bayes=NULL;
    
    //BurnIn input
    BurnInInput* input=(BurnInInput*)malloc(sizeof(BurnInInput));

    int ndims;

    /***PARSE/PROCESS INPUT***/

    if (!PyArg_ParseTuple(args,"O!O!O!O!",&PyList_Type,&py_inputfile_list,&PyBool_Type,&py_spin,&PyFloat_Type,&py_deltaLogL,&PyString_Type,&py_outputfilename))  return NULL;

    input->deltaLogL=PyFloat_AsDouble(py_deltaLogL);

    input->nfiles=PyList_Size(py_inputfile_list);

    input->files = (char**)calloc(MAX_IN_FILES, sizeof(char*));

    Py_ssize_t i;
    for(i=0;i<input->nfiles;i++){
    	char* ptemp=NULL;
    	if((ptemp=PyString_AsString(PyList_GetItem(py_inputfile_list,i)))!=0){
        	input->files[i]=ptemp;
        }
        else {
        	nfiles--;
        }
    }
    
    input->output_file_name = PyString_AsString(py_outputfilename);

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


static struct PyMethodDef BUmethods[] = {
    {"_burnin", _burnin, METH_VARARGS,
    "This function 'burns in' MCMC chains."
    },
    {NULL,}
};


void init_bayespputils(void);
void init_bayespputils(void)
{
    (void)Py_InitModule3(MODULE_NAME,BUmethods,BUDocstring);
    import_array();
}
