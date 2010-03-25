/**
 * nds_ext.cpp
 * Python wrapper for John Zweizig's NDS1/NDS2 library
 *
 * @author Leo Singer <leo.singer@ligo.org>
 *
 */

#include <boost/python.hpp>
#include <ostream>

extern "C" {
#include <daqc.h>
#include <daqc_internal.h>
#include <daqc_response.h>
}

#include <numpy/arrayobject.h>

using namespace boost::python;

// Can J. Zweizig add this to NDS library?
static const char* daq_strerror(int errornum)
{
    switch (errornum)
    {
        case DAQD_OK:
            return "OK";
        case DAQD_ERROR:
            return "unspecified error";
        case DAQD_NOT_CONFIGURED:
            return "not configured";
        case DAQD_INVALID_IP_ADDRESS:
            return "invalid IP address";
        case DAQD_INVALID_CHANNEL_NAME:
            return "invalid channel name";
        case DAQD_SOCKET:
            return "socket";
        case DAQD_SETSOCKOPT:
            return "setsockopt";
        case DAQD_CONNECT:
            return "connect";
        case DAQD_BUSY:
            return "busy";
        case DAQD_MALLOC:
            return "malloc";
        case DAQD_WRITE:
            return "write";
        case DAQD_VERSION_MISMATCH:
            return "version mismatch";
        case DAQD_NO_SUCH_NET_WRITER:
            return "no such net writer";
        case DAQD_NOT_FOUND:
            return "not found";
        case DAQD_GETPEERNAME:
            return "getpeername";
        case DAQD_DUP:
            return "dup";
        case DAQD_INVALID_CHANNEL_DATA_RATE:
            return "invalid channel data rate";
        case DAQD_SHUTDOWN:
            return "shutdown";
        case DAQD_NO_TRENDER:
            return "no trender";
        case DAQD_NO_MAIN:
            return "no main";
        case DAQD_NO_OFFLINE:
            return "no offline";
        case DAQD_THREAD_CREATE:
            return "thread create";
        case DAQD_TOO_MANY_CHANNELS:
            return "too many channels";
        case DAQD_COMMAND_SYNTAX:
            return "command syntax";
        case DAQD_SASL:
            return "sasl";
        case DAQD_NOT_SUPPORTED:
            return "not supported";
        default:
            return "unknown error";
    }
}

struct DaqError {
public:
    DaqError(int retval) : retval(retval) {}
    int retval;
};

struct NoMemoryError {
};

static void DaqErrorTranslator(const DaqError& exc)
{
    char* s;
    asprintf(&s, "error %d: %s", exc.retval, daq_strerror(exc.retval));
    PyErr_SetString(PyExc_RuntimeError, s);
    free(s);
}

static void NoMemoryErrorTranslator(const NoMemoryError& exc)
{
    PyErr_NoMemory();
}

static int numpy_typenum_for_daq_data_t(daq_data_t type)
{
    switch (type)
    {
        case _16bit_integer:
            return NPY_INT16;
        case _32bit_integer:
            return NPY_INT32;
        case _64bit_integer:
            return NPY_INT64;
        case _32bit_float:
            return NPY_FLOAT32;
        case _64bit_double:
            return NPY_FLOAT64;
        case _undefined:
        default:
            return NPY_VOID;
    }
}

struct _daq_t : daq_t {
public:
    
    _daq_t(const std::string& host, int port, nds_version version) : daq_t()
    {
        daq_startup();
        int retval = daq_connect(this, host.c_str(), port, version);
        if (retval) throw DaqError(retval);
    }
    
    _daq_t(const std::string& host, int port) : daq_t()
    {
        daq_startup();
        PyErr_Warn(PyExc_RuntimeWarning, "No protocol specified, attempting protocol nds_v2");
        int retval = daq_connect(this, host.c_str(), port, nds_v2);
        if (!retval)
            return;
        PyErr_Warn(PyExc_RuntimeWarning, "Protocol nds_v2 failed, falling back to nds_v1");
        retval = daq_connect(this, host.c_str(), port, nds_v1);
        if (retval)
            throw DaqError(retval);
    }
    
    ~_daq_t()
    {
        daq_disconnect(this);
    }
    
    void disconnect()
    {
        int retval = daq_disconnect(this);
        if (retval) throw DaqError(retval);
    }
    
    void clear_channel_list()
    {
        int retval = daq_clear_channel_list(this);
        if (retval) throw DaqError(retval);
    }
    
    void recv_next()
    {
        int retval = daq_recv_next(this);
        if (retval) throw DaqError(retval);
    }
    
    void request_data(time_t start, time_t end, time_t dt)
    {
        int retval = daq_request_data(this, start, end, dt);
        if (retval) throw DaqError(retval);
    }
    
    void request_channel(daq_channel_t* channel)
    {
        int retval = daq_request_channel_from_chanlist(this, channel);
        if (retval) throw DaqError(retval);
    }
    
    list unpack()
    {
        list l;
        for (chan_req_t* channel = this->chan_req_list; channel < &this->chan_req_list[this->num_chan_request]; channel++)
        {
            if (channel->status < 0)
                throw DaqError(- channel->status);
            int data_length = channel->status;
            npy_intp dim = data_length / data_type_size(channel->data_type);
            PyObject* array_obj = PyArray_SimpleNew(1, &dim, numpy_typenum_for_daq_data_t(channel->data_type));
            if (!array_obj) throw NoMemoryError();
            memcpy(PyArray_DATA(array_obj), this->tb->data + channel->offset, data_length);
            l.append(make_tuple(channel, handle<>(array_obj)));
        }
        return l;
    }
    
    list get_requested_channels()
    {
        list l;
        for (chan_req_t* channel = this->chan_req_list; channel < &this->chan_req_list[this->num_chan_request]; channel++)
            l.append(channel);
        return l;
    }
    
    list recv_channel_list(chantype channeltype = cUnknown)
    {
        daq_channel_t* channels = (daq_channel_t*) calloc(MAX_CHANNELS, sizeof(daq_channel_t));
        if (!channels) throw NoMemoryError();
        int nchannels_received;
        int retval = daq_recv_channel_list(this, channels, MAX_CHANNELS, &nchannels_received, 0, channeltype);
        if (retval)
        {
            free(channels);
            throw DaqError(retval);
        }
        list l;
        for (daq_channel_t* channel = channels ; channel < &channels[nchannels_received]; channel++)
            l.append(*channel);
        free(channels);
        return l;
    }
    
    list recv_channel_list_any()
    {
        return recv_channel_list(cUnknown);
    }
    
};

static str _signal_conv_get_units(const signal_conv_t* self)
{
    return str(self->signal_units);
}

static str _daq_channel_get_name(const daq_channel_t* self)
{
    return str(self->name);
}

static str _chan_req_get_name(const chan_req_t* self)
{
    return str((const char*)self->name);
}

template<class CHANNEL_T>
static str _channel_as_str(const CHANNEL_T* c)
{
    return str("%s (%dHz, %s)" % make_tuple(str((const char*)c->name), int(c->rate), c->type));
}

template<class CHANNEL_T>
static str _channel_as_repr(const CHANNEL_T* c)
{
    return str("<%s>" % _channel_as_str(c));
}

BOOST_PYTHON_MODULE(nds_ext)
{
    import_array();
    
    enum_<chantype>("channel_type")
        .value("unknown", cUnknown)
        .value("online", cOnline)
        .value("raw", cRaw)
        .value("reduced", cRDS)
        .value("second_trend", cSTrend)
        .value("minute_trend", cMTrend)
        .value("testpoint", cTestPoint);
    
    enum_<nds_version>("nds_version", "NDS protocol version")
        .value("v1", nds_v1)
        .value("v2", nds_v2);
    
    class_<signal_conv_t>("conversion", "unit conversion information")
        .def_readonly("gain", &signal_conv_t::signal_gain)
        .def_readonly("slope", &signal_conv_t::signal_slope)
        .def_readonly("offset", &signal_conv_t::signal_offset)
        .add_property("units", _signal_conv_get_units);
    
    class_<daq_channel_t>("channel")
        .add_property("name", _daq_channel_get_name)
        .def_readonly("rate", &daq_channel_t::rate)
        .def_readonly("type", &daq_channel_t::type)
        .def_readonly("conversion", &daq_channel_t::s)
        .def("__str__", &_channel_as_str<daq_channel_t>)
        .def("__repr__", &_channel_as_repr<daq_channel_t>);
    
    class_<chan_req_t>("channel_request")
        .add_property("name", _chan_req_get_name)
        .def_readonly("rate", &chan_req_t::rate)
        .def_readonly("type", &chan_req_t::type)
        .def_readonly("conversion", &chan_req_t::s)
        .def("__str__", &_channel_as_str<chan_req_t>)
        .def("__repr__", &_channel_as_repr<chan_req_t>);
    
    register_exception_translator<DaqError>(DaqErrorTranslator);
    
    register_exception_translator<NoMemoryError>(NoMemoryErrorTranslator);
    
    class_<_daq_t>("daq", init<const std::string&, int, nds_version>(args("host", "port", "nds_version")))
        .def(init<const std::string&, int>(args("host", "port")))
        .def("disconnect", &_daq_t::disconnect)
        .def("recv_channel_list", &_daq_t::recv_channel_list)
        .def("recv_channel_list", &_daq_t::recv_channel_list_any)
        .def("clear_channel_list", &_daq_t::clear_channel_list, "reset list of requested channels")
        .def("request_data", &_daq_t::request_data, args("gps_start", "gps_end", "stride"))
        .def("request_channel", &_daq_t::request_channel, "request a channel")
        .def("recv_next", &_daq_t::recv_next, "received next block")
        .def("unpack", &_daq_t::unpack, "unpack received data")
        .add_property("requested_channels", &_daq_t::get_requested_channels, "list of requested channels")
        .def_readonly("rev", &_daq_t::rev, "revision of the underlying NDS client library")
        .def_readonly("nds_version", &_daq_t::nds_versn, "version of NDS protocol used for this connection");
    
}
