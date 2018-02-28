/*
 * ============================================================================
 *
 *                C extensions for Bayesian post-processing codes
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>

#include "bayespputils.h"

#define MODULE_NAME "_bayespputils"

/* doc string */
const char BUDocstring[] =
"This module provides C extensions for Bayesian analysis and post-\n"
"processing codes.";

/* calculateConfidenceLevels
 * C extension function replacing pylal.bayespputils.calculateConfidenceLevels_slow .
 */
/* skyhist_cart
 * C extension function replacing pylal.bayespputils.skyhist_cart .
 */
static PyObject *_skyhist_cart(PyObject *self, PyObject *args)
{
    /***DECLARATION***/

    //INPUT
    PyArrayObject *SkyCartsArray=NULL,*SamplesArray=NULL;

    //OUTPUT
    PyArrayObject *BinsArray=NULL;

    //INTERNAL VARIABLES
    int skybin_idx,sample_no;//iterators
    int a,b;//SkyCarts dimensions
    int m,n;//Samples dimensions

    //containers for loop variables
    double xsample,ysample,zsample,xskybin,yskybin,zskybin;
    double longi,lat,maxbin_value,dot;
    int maxbin;

    /***PARSE/PROCESS INPUT***/

    if (!PyArg_ParseTuple(args, "O!O!",&PyArray_Type, &SkyCartsArray,&PyArray_Type, &SamplesArray))  return NULL;

    if(SkyCartsArray->nd != 2 || SkyCartsArray->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetObject(PyExc_ValueError, (PyObject *) SkyCartsArray);
        return NULL;
    }

    if( (SamplesArray->nd !=2 ) || SamplesArray->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetObject(PyExc_ValueError, (PyObject *) SamplesArray);
        return NULL;
    }

    //get dimensions
    a=SkyCartsArray->dimensions[0];
    b=SkyCartsArray->dimensions[1];

    m=SamplesArray->dimensions[0];
    n=SamplesArray->dimensions[1];


    npy_intp bin_dims[1] = {a};

    //create bins vector
    BinsArray=(PyArrayObject *) PyArray_SimpleNew(1,bin_dims,PyArray_INT);
    Py_INCREF(BinsArray);

    if(BinsArray==NULL) return NULL;

    //get iterators
    int SamplesArray_stride_i=SamplesArray->strides[0];
    int SamplesArray_stride_j=SamplesArray->strides[1];

    int SkyCartsArray_stride_i=SkyCartsArray->strides[0];
    int SkyCartsArray_stride_j=SkyCartsArray->strides[1];

    int BinsArray_stride_i=BinsArray->strides[0];

    int SkyCartsArray_iterx;

    //zero the bins
    for(skybin_idx=0;skybin_idx<a;skybin_idx++){
        *(int*)(BinsArray->data+skybin_idx*BinsArray_stride_i) = 0;
    }

    int SamplesArray_stride_j_RA_iter=0;
    int SamplesArray_stride_j_dec_iter=SamplesArray_stride_j;


    /***MAIN LOOP***/
    /* Now going to loop over samples. For each sapmle the max. dot product
     * corresponds to the sky bin that its going to be associated with.
     */

    for(sample_no=0;sample_no<m;sample_no++){

        longi = *(double*)(SamplesArray->data+sample_no*SamplesArray_stride_i+SamplesArray_stride_j_RA_iter);
        lat = *(double*)(SamplesArray->data+sample_no*SamplesArray_stride_i+SamplesArray_stride_j_dec_iter);

        //Cartesian co of sample for dot product
        xsample=cos(lat)*cos(longi);
        ysample=cos(lat)*sin(longi);
        zsample=sin(lat);
        maxbin=0;
        maxbin_value=-1.;
        for(skybin_idx=0;skybin_idx<a;skybin_idx++){

            SkyCartsArray_iterx=skybin_idx*SkyCartsArray_stride_i;

            //Cartesian co of bin for dot product
            xskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx);
            yskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx+SkyCartsArray_stride_j);
            zskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx+2*SkyCartsArray_stride_j);

            //dot (x,y,z) with each skybin
            dot=( xsample*xskybin + ysample*yskybin + zsample*zskybin );
            if(  dot > maxbin_value ){//if dot is the biggest so far update the maxbin info
                maxbin = skybin_idx;
                maxbin_value = dot;

            }

        }

        //add one to the max bin count for the highest scoring dot product
        *(int*)(BinsArray->data+maxbin*BinsArray_stride_i) += 1;
    }

    return PyArray_Return(BinsArray);

}

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

    py_pos_array=PyArray_SimpleNewFromData(2,pos_dims,PyArray_DOUBLE,(void*)output->pos);
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
    {"_skyhist_cart", _skyhist_cart, METH_VARARGS,
    "This function bins samples on a grid generated from pylal.skylocutils\n"
    "from the marginal posterior for the sky position by finding the maximum\n"
    "dot product from each sky bin for each sample."
    },
    {"_burnin", _burnin, METH_VARARGS,
    "This function 'burns in' MCMC chains."
    },
    {NULL,}
};


void init_bayespputils(void)
{
    (void)Py_InitModule3(MODULE_NAME,BUmethods,BUDocstring);
    import_array();
}
