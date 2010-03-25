#include "Python.h"
#include "structmember.h"

#include "daqc.h"
#include "daqc_internal.h"

#define MAX_CHANNEL_LIST 32768
#define NEW_VECT(type,dim) ((type*)malloc(dim*sizeof(type)))


/**********************************************
 * Wrapper around daq_t
 *********************************************/
typedef struct {
    PyObject_HEAD
    daq_t *daq;
} _nds2_DaqObject;

static void Daq_dealloc(_nds2_DaqObject *self)
{
    if (self->daq != NULL) {
        free(self->daq);
    }

    self->ob_type->tp_free((PyObject*) self);
}


static PyObject *Daq_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    _nds2_DaqObject *self;

    self = (_nds2_DaqObject *) type->tp_alloc(type, 0);

    self->daq = (daq_t *) malloc(sizeof(daq_t));

    return (PyObject *)self;
}


static int Daq_init(_nds2_DaqObject *self, PyObject *args, PyObject *kwds)
{
    self->daq = (daq_t *) malloc(sizeof(daq_t));

    return 0;
}

static PyTypeObject _nds2_DaqType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_nds2.daq_t",             /*tp_name*/
    sizeof(_nds2_DaqObject),   /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "daq_t objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Daq_init,        /* tp_init */
    0,                         /* tp_alloc */
    Daq_new,                   /* tp_new */
};



/**********************************************
 * Wrapper around daq_channel_t
 *********************************************/
typedef struct {
    PyObject_HEAD
    daq_channel_t *channel;
} _nds2_DaqChannelObject;


static void DaqChannel_dealloc(_nds2_DaqChannelObject *self)
{
    if (self->channel != NULL) {
        free(self->channel);
    }

    self->ob_type->tp_free((PyObject*) self);
}


static PyObject *DaqChannel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    _nds2_DaqChannelObject *self;

    self          = (_nds2_DaqChannelObject *) type->tp_alloc(type, 0);
    self->channel = (daq_channel_t *) malloc(sizeof(daq_channel_t));

    return (PyObject *)self;
}


static int DaqChannel_init(_nds2_DaqChannelObject *self, PyObject *args, PyObject *kwds)
{
    self->channel = (daq_channel_t *) malloc(sizeof(daq_channel_t));
    return 0;
}



static PyMemberDef DaqChannel_members[] = {
    {NULL}
};

static PyObject *DaqChannel_getname(_nds2_DaqChannelObject *self, void *closure) {
    PyObject *name;

    if (self->channel != NULL) {
        name = Py_BuildValue("s", (char *) self->channel);
    } else {
        name = Py_BuildValue("s","");
    }

    Py_INCREF(name);
    return name;
}


static int DaqChannel_setname(_nds2_DaqChannelObject *self, PyObject *value, void *closure) {
    char *name;

    name = PyString_AsString(value);
    
    strcpy((char *) self->channel, name);

    return 0;
}

static PyObject *DaqChannel_getrate(_nds2_DaqChannelObject *self, void *closure) {
    PyObject *rate;

    if (self->channel != NULL) {
        rate = Py_BuildValue("d", self->channel->rate);
    } else {
        rate = Py_BuildValue("d",0.0);
    }

    return rate;
}

static int DaqChannel_setrate(_nds2_DaqChannelObject *self, PyObject *value, void *closure) {
    self->channel->rate = PyFloat_AsDouble(value);

    return 0;
}


static PyObject *DaqChannel_gettpnum(_nds2_DaqChannelObject *self, void *closure) {
    PyObject *tpn;

    if (self->channel != NULL) {
        tpn = Py_BuildValue("i", self->channel->tpnum);
    } else {
        tpn = Py_BuildValue("i",0);
    }

    return tpn;
}

static int DaqChannel_settpnum(_nds2_DaqChannelObject *self, PyObject *value, void *closure) {
    self->channel->tpnum = PyInt_AsLong(value);

    return 0;
}


static PyObject *DaqChannel_getbps(_nds2_DaqChannelObject *self, void *closure) {
    PyObject *bps;

    if (self->channel != NULL) {
        bps = Py_BuildValue("i", self->channel->bps);
    } else {
        bps = Py_BuildValue("i",0);
    }

    return bps;
}

static int DaqChannel_setbps(_nds2_DaqChannelObject *self, PyObject *value, void *closure) {
    self->channel->bps = PyInt_AsLong(value);

    return 0;
}


static PyObject *DaqChannel_getchNum(_nds2_DaqChannelObject *self, void *closure) {
    PyObject *chNum;

    if (self->channel != NULL) {
        chNum = Py_BuildValue("i", self->channel->chNum);
    } else {
        chNum = Py_BuildValue("i",0);
    }

    return chNum;
}

static int DaqChannel_setchNum(_nds2_DaqChannelObject *self, PyObject *value, void *closure) {
    self->channel->chNum = PyInt_AsLong(value);

    return 0;
}

static PyGetSetDef DaqChannel_getseters[] = {
    {"name",  (getter) DaqChannel_getname, (setter) DaqChannel_setname, "name", NULL},
    {"rate",  (getter) DaqChannel_getrate, (setter) DaqChannel_setrate, "rate", NULL},
    {"tpnum", (getter) DaqChannel_gettpnum, (setter) DaqChannel_settpnum, "test point number", NULL},
    {"bps",   (getter) DaqChannel_getbps,   (setter) DaqChannel_setbps,   "bytes per sample", NULL},
    {"chNum", (getter) DaqChannel_getchNum, (setter) DaqChannel_getchNum, "channel number", NULL},
    {NULL}
};


static PyMethodDef daq_channel_t_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject _nds2_DaqChannelType = {
    PyObject_HEAD_INIT(NULL)
    0,                               /*ob_size*/
    "_nds2.daq_channel_t",           /*tp_name*/
    sizeof(_nds2_DaqChannelObject),  /*tp_basicsize*/
    0,                               /*tp_itemsize*/
    (destructor) DaqChannel_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "daq_channel_t objects",   /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0, /* DaqChannel_methods, */            /* tp_methods */
    0, /* DaqChannel_members, */            /* tp_members */
    DaqChannel_getseters,                    /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)DaqChannel_init,  /* tp_init */
    0,                         /* tp_alloc */
    DaqChannel_new,                 /* tp_new */
};




/**********************************************
 * Result object
 *********************************************/
typedef struct {
    PyObject_HEAD
    char *name;
    int chan_type;
    int rate;
    int data_type;
    int signal_gain;
    int signal_offset;
    int signal_slope;
    int signal_units;
    int exists;
    int start_gps_sec;
    int field_duration_sec;
    PyObject *data;
} _nds2_ReturnObject;


static void ReturnObject_dealloc(_nds2_ReturnObject *self)
{
    if (self->name) {
        free(self->name);
    }

    if (self->data) {
        /* Free all the objects in the data */
        /* Free the data array */
    }

    self->ob_type->tp_free((PyObject*) self);
}


static PyObject *ReturnObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    _nds2_ReturnObject *self;

    self = (_nds2_ReturnObject *) type->tp_alloc(type, 0);

    return (PyObject *)self;
}


static int ReturnObject_init(_nds2_ReturnObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

static PyObject *ReturnObject_get_name(_nds2_ReturnObject *self, void *closure) {
    PyObject *name;

    if (self->name != NULL) {
        name = Py_BuildValue("s", self->name);
    } else {
        name = Py_BuildValue("s","");
    }

    Py_INCREF(name);
    return name;
}

static int ReturnObject_set_name(_nds2_ReturnObject *self, PyObject *value, void *closure) {
    char *name = PyString_AsString(value);
    
    self->name = (char *) malloc(strlen(name)+1);

    strcpy((char *) self->name, name);

    return 0;
}

static PyObject *ReturnObject_get_rate(_nds2_ReturnObject *self, void *closure) {
    return Py_BuildValue("i", self->rate);
}

static int ReturnObject_set_rate(_nds2_ReturnObject *self, PyObject *value, void *closure) {
    self->rate = PyInt_AsLong(value);
    return 0;
}

static PyObject *ReturnObject_get_data(_nds2_ReturnObject *self, void *closure) {
    return self->data;
}

static int ReturnObject_set_data(_nds2_ReturnObject *self, PyObject *value, void *closure) {
    self->data = value;
    /* TODO: decrement count on data, inc on value */

    return 0;
}
static PyGetSetDef ReturnObject_getseters[] = {
    {"name",        (getter) ReturnObject_get_name, (setter) ReturnObject_set_name, "name", NULL},
    {"signal_rate", (getter) ReturnObject_get_rate, (setter) ReturnObject_set_rate, "rate", NULL},
    {"data",        (getter) ReturnObject_get_data, (setter) ReturnObject_set_data, "data", NULL},
    {NULL}
};


static PyMethodDef return_t_methods[] = {
    {NULL}  /* Sentinel */
};


static PyTypeObject _nds2_ReturnType = {
    PyObject_HEAD_INIT(NULL)
    0,                            /*ob_size*/
    "_nds2.return_t",             /*tp_name*/
    sizeof(_nds2_ReturnObject),   /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "return_t objects",        /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    ReturnObject_getseters,    /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ReturnObject_init,  /* tp_init */
    0,                            /* tp_alloc */
    ReturnObject_new,             /* tp_new */
};



/* Here's what we really want to do:
   (from: http://www.rmi.net/~lutz/newex.html)

    array_mod   = PyImport_AddModule("array");
    array_class = PyObject_GetAttrString(array_mod, "array");
    rargs       = Py_BuildValue("(s)","d");
    inst        = PyEval_CallObject(array_class, rargs);
    meth        = PyObject_GetAttrString(inst, "fromstring");
    rargs       = Py_BuildValue("(s)", "00000000");
    res         = PyEval_CallObject(meth, rargs);

    Except instead of 00000000 use
    ret.fromstring(PyString_FromStringAndSize(the data, len(the_data) * sizeof(double))

*/
PyObject *dump_data(daq_t *daq, chan_req_t *chan) {
    long N;
    PyObject *array_mod, *array_class, *data_string, *args, *array_inst, *method, *ret;

    /* Get a pointer to the array module.  This wont import the module, */
    /* so it has be loaded in init */
    array_mod = PyImport_AddModule("array");
    if (array_mod == NULL) {
        fprintf(stderr, "Unable to load array module\n");
        return NULL;
    }

    /* Get the array.array class */
    array_class = PyObject_GetAttrString(array_mod, "array");
    
    if(array_class == NULL) {
        fprintf(stderr, "Unable to get array.array class\n");
        return NULL;
    }

	N = chan->status;

	if (N < 0) {
	    printf("Channel: %s receive error (%li)\n", chan->name, N);	    
	    return NULL;
	}
	
	switch (chan->data_type) {
	    case _16bit_integer: 
            args = Py_BuildValue("(s)","h");
#ifdef DEBUG
            printf("_16bit_integer\n");
#endif
	        break;
	    case _32bit_integer: 
            args = Py_BuildValue("(s)","i");
#ifdef DEBUG
            printf("_32bit_integer\n");
#endif
	        break;
	    case _32bit_float: 
            args = Py_BuildValue("(s)","f");
#ifdef DEBUG
            printf("_32bit_float\n");
#endif
	        break;
	    case _64bit_double: 
            args = Py_BuildValue("(s)","d");
#ifdef DEBUG
            printf("_64bit_double\n");
#endif
	        break;
	    default:
	        fprintf(stderr, "Channel: %s data of unknown/unsupprted type (%i)\n", chan->name, chan->data_type);
	}

    /* Call the constructor */
    array_inst = PyEval_CallObject(array_class, args);
    Py_DECREF(array_class);
    Py_DECREF(args);

    if (array_inst == NULL) {
        fprintf(stderr,"Unable to call array.array()");
	    return NULL;
	}

    /* Find the fromstring method */
    method = PyObject_GetAttrString(array_inst, "fromstring");

    if (method == NULL) {
        fprintf(stderr,"Unable to find array.fromstring()");
        return NULL;
    }

    /* Make a string from the results */
    data_string = PyString_FromStringAndSize(daq->tb->data + chan->offset, N);

    /* Add it to the array */
    args = Py_BuildValue("(O)", data_string);


    PyEval_CallObject(method, args);

    /* Don't decriment here, the object needs to maintain a reference */
    /* Py_DECREF(data_string); */
    Py_DECREF(method);
    Py_DECREF(args);

#ifdef DEBUG
    printf("Count of array_class: %d\n", array_class->ob_refcnt);
    printf("Count of data_string: %d\n", data_string->ob_refcnt);
    printf("Count of method: %d\n", method->ob_refcnt);
    printf("Count of args: %d\n", args->ob_refcnt);
#endif

    return array_inst;
}

static PyObject *wrap_daq_startup(PyObject *self, PyObject *args)
{
    int ret;

    ret = daq_startup();

    return Py_BuildValue("i", ret);
}


static PyObject *wrap_daq_connect(PyObject *self, PyObject *args)
{
    enum nds_version vrsn = nds_v1;

    const char *node_id;
    int port_id;
    _nds2_DaqObject *daq;
    int rc;

    if (!PyArg_ParseTuple(args, "Osi", &daq, &node_id, &port_id))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }

    if (port_id != DAQD_PORT) vrsn = nds_v2;
    rc = daq_connect(daq->daq, node_id, port_id, vrsn);

    return Py_BuildValue("i", rc);
}

static PyObject *wrap_daq_recv_channels(PyObject *self, PyObject *args)
{
    int num;
    int i;
    int rc;
    PyObject *obj;
    _nds2_DaqObject *daq;
    _nds2_DaqChannelObject *inner;

    if (!PyArg_ParseTuple(args, "OO", &daq, &obj))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }
    
    num = PySequence_Length(obj);

    daq_channel_t* channel_list = NEW_VECT(daq_channel_t, num);
    int nChans = 0;

    rc = daq_recv_channels(daq->daq, channel_list, num, &nChans);

    if (num < nChans)
        nChans = num;

    if (rc == 0) {
        for (i=0; i < nChans; i++) {
            inner = (_nds2_DaqChannelObject *) PySequence_GetItem(obj, i);

            memcpy(inner->channel, channel_list + i, sizeof(daq_channel_t));
        }
    }

    return Py_BuildValue("i", rc);
}



static PyObject *wrap_daq_request_channel_from_chanlist(PyObject *self, PyObject *args)
{
    _nds2_DaqObject *daq;
    _nds2_DaqChannelObject *channel;
    int rc;

    if (!PyArg_ParseTuple(args, "OO", &daq, &channel))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }


    if (!PyObject_IsInstance((PyObject*) channel, (PyObject*) &_nds2_DaqChannelType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_channel_t required for first argument");
                return (PyObject *) -1;
    }

    rc = daq_request_channel_from_chanlist(daq->daq, channel->channel);

    return Py_BuildValue("i", rc);
}

static PyObject *wrap_daq_request_data(PyObject *self, PyObject *args)
{
    unsigned long gps_start, gps_end, delta;
    unsigned long gps;
    PyObject *ret;
    _nds2_DaqObject *daq;
    int rc;
    int i;

    if (!PyArg_ParseTuple(args, "Okkk", &daq, &gps_start, &gps_end, &delta))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }

    rc = daq_request_data(daq->daq, gps_start, gps_end, delta);

    for (i=0; i < daq->daq->num_chan_request; i++) {
	    chan_req_t* chan = daq->daq->chan_req_list + i;
        ret = dump_data(daq->daq, daq->daq->chan_req_list + i);
    
        for (gps=gps_start+delta; gps<gps_end; gps+=delta) {
            long count = daq_recv_next(daq->daq);

            if (count > 0) {
                ret = dump_data(daq->daq, daq->daq->chan_req_list + i);
            } else {
                printf("Error in daq_recv_block: %li\n", count);
                break;
            }
        }
    }

    return ret;
}

static PyObject *wrap_daq_disconnect(PyObject *self, PyObject *args)
{
    _nds2_DaqObject *daq;
    int rc;

    if (!PyArg_ParseTuple(args, "O", &daq))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }

    rc = daq_disconnect(daq->daq);
    return Py_BuildValue("i", rc);
}


static PyObject *wrap_daq_recv_shutdown(PyObject *self, PyObject *args)
{
    _nds2_DaqObject *daq;
    int rc;

    if (!PyArg_ParseTuple(args, "O", &daq))
        return NULL;

    if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
        PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
                return (PyObject *) -1;
    }

    rc = daq_recv_shutdown(daq->daq);

    return Py_BuildValue("i", rc);
}

static PyObject *get_data(PyObject *self, PyObject *args)
{
	_nds2_DaqObject *daq;
    PyObject        *requests;
    PyObject        *array_mod, *array_class, *rargs, *inst, *meth, *res, *pmeth;
	int             rc;
	PyObject *      ret;
	_nds2_ReturnObject *newRet;
	char            *name;
	chantype_t ctype = cUnknown;                                        
	double     rate  = 0.0; 

    int start_time;
    int end_time;

	int i;

	if (!PyArg_ParseTuple(args, "OOii", &daq, &requests, &start_time, &end_time))
		return NULL;

	if (!PyObject_IsInstance((PyObject*) daq, (PyObject*) &_nds2_DaqType)) {
		PyErr_SetString(PyExc_TypeError, "instance of daq_t required for first argument");
		return (PyObject *) -1;
	}

    for(i=0; i< PySequence_Length(requests); i++) {
        name = PyString_AsString(PySequence_GetItem(requests, i));
	    daq_request_channel(daq->daq, name, ctype, rate);
    }

	int err = daq_request_data(daq->daq, start_time, end_time, (end_time - start_time));

	if (err) {
        printf("daq_request_data failed\n"); 
        daq_clear_channel_list(daq->daq);
        return NULL;                
    }

	int nChanReq = daq->daq->num_chan_request;                                               

    ret = PyList_New(0);

	for(i=0; i < nChanReq; i++) {
        newRet          = (_nds2_ReturnObject *) _nds2_ReturnType.tp_alloc(&_nds2_ReturnType, 0);
        chan_req_t* req = daq->daq->chan_req_list + i;

        newRet->name          = (char *) malloc(strlen(req->name) + 1);
        strcpy(newRet->name, req->name);

		newRet->rate          = req->rate;
		newRet->signal_gain   = req->s.signal_gain;
		newRet->signal_offset = req->s.signal_offset;
		newRet->signal_slope  = req->s.signal_slope;

		if(req->status > 0) {
            newRet->data = dump_data(daq->daq, req);
		}

		PyList_Append(ret, (PyObject *) newRet);
	}

    daq_clear_channel_list(daq->daq);

	return ret;
}


static PyMethodDef NDS2Methods[] = {
    {"daq_startup",       wrap_daq_startup, METH_VARARGS, "Setup daq for use"},
    {"daq_connect",       wrap_daq_connect, METH_VARARGS, "Connect to server (you must kinit first)"},
    {"daq_recv_channels", wrap_daq_recv_channels, METH_VARARGS, "Receive a list of channels"},
    {"daq_request_channel_from_chanlist", wrap_daq_request_channel_from_chanlist, METH_VARARGS, "Request data from a listed channel"},
    {"daq_request_data",  wrap_daq_request_data, METH_VARARGS,  "Request data from a named channel"},
    {"daq_disconnect",    wrap_daq_disconnect, METH_VARARGS,    "Disconnect from server"},
    {"daq_recv_shutdown", wrap_daq_recv_shutdown, METH_VARARGS, "Shutdown and cleanup connection"},
    {"get_data",          get_data,               METH_VARARGS, "data request"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC init_nds2(void)
{
    PyObject *m;
    PyObject *array_mod;
    m = Py_InitModule("_nds2", NDS2Methods);

    _nds2_DaqChannelType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&_nds2_DaqChannelType) < 0)
        return;

    Py_INCREF(&_nds2_DaqChannelType);
    PyModule_AddObject(m, "daq_channel_t", (PyObject *)&_nds2_DaqChannelType);


    _nds2_DaqType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&_nds2_DaqType) < 0)
        return;

    Py_INCREF(&_nds2_DaqType);
    PyModule_AddObject(m, "daq_t", (PyObject *)&_nds2_DaqType);


    _nds2_ReturnType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&_nds2_ReturnType) < 0)
        return;

    Py_INCREF(&_nds2_ReturnType);
    PyModule_AddObject(m, "return_t", (PyObject *)&_nds2_ReturnType);

    array_mod = PyImport_ImportModule("array");

    if (! array_mod) 
    {
        fprintf(stderr,"Unable to load array module, data requests will not work.");
    }
}

