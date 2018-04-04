"""
  MetaArray - a subclass of ndarray that holds metadata and preserves it across
            array operations.
  Metadata - a class for metadata stored in MetaArray
  MetaArrayList - a subclass of list that ollows element-wise MetaArray
                operations

  Spectrum - a subclass of MetaArray as demonstration
  SpectrumMetadata - a subclass of MetaData that strictly checks for metadata
                   compatibility between two Spectra
  SpectrumList - subclass of MetaArrayList that has some nice features specific
               to SpectrumMetadata
"""
from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"

import operator
import os
import sys
import copy

import numpy
from glue import segments

##############################################################################
# Method wrappers
##############################################################################

class _arraymethod(object):
    """
    Used to attach methods to MetaArrays.  Wrap ndarray methods that return
    a single array.  Merge metadata of all input Spectra.
    
    __init__ gets called when we define the MetaArray (sub-)class and attach
    methods.
    __get__ gets called on the object at method call time.
    __call__ is called immediately after __get__.
    """
    def __init__(self, methname):
        self._name = methname
        self.obj = None
        self.__doc__ = getattr(numpy.ndarray, self._name, None).__doc__
    
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    
    def __call__(self, *args, **params):
        data = self.obj.A
        cls = type(self.obj)
        result = getattr(data, self._name)(*args, **params).view(cls)
        result.metadata = self.obj.metadata.copy()
        for arg in args:
            result.metadata |= getattr(arg, 'metadata', None)
        return result

class _elementwise_method(object):
    """
    Used to attach methods to MetaArrayLists.  Wrap ndarray methods that
    operate upon arrays and apply them to each element of the list.
    
    __init__ gets called when we define the MetaArrayList class and attach
    methods.
    __get__ gets called on the list object at method call time.
    __call__ is called immediately after __get__.
    """
    def __init__(self, methname):
        self._name = methname
        self.list = None
        self.listtype = None
        self.itemtype = None
        # I don't think I can get the docstring, unfortunately
        # self.__doc__ = getattr(itemtype, self._name, None).__doc__

    def __get__(self, list, listtype=None):
        self.list = list
        self.listtype = type(list)
        self.itemtype = list._itemtype
        return self

    def __call__(self, *others, **params):
        result = self.listtype([])
        cls = self.itemtype
        # if we have two alike lists, just take pairwise action
        if len(others) == 1 and isinstance(others[0], self.listtype):
            result.extend([getattr(a, self._name)(b, **params).view(cls) \
                for a, b in zip(self.list, others[0])])
        # otherwise, just try passing the args along
        else:
            result.extend([getattr(a, self._name)(*others, **params).view(cls) for a in self.list])
        return result

##############################################################################
# Generic classes
##############################################################################

class Metadata(object):
    """
    Abstract class to hold metadata
    """
    # specify all allowed metadata here; type is required
    typemap = {}
    __slots__ = typemap.keys() + ['typemap']
    
    def __new__(cls, metadata=None):
        if isinstance(metadata, Metadata):
            return metadata
        return super(Metadata, cls).__new__(cls)
    
    def __init__(self, metadata):
        # Recycle if possible
        if isinstance(metadata, Metadata):
            return
        # user forgot to provide metadata
        elif metadata is None:
            return None
        # create Metadata from dictionary
        elif isinstance(metadata, dict):
            for key, val in metadata.items():
                try:
                    setattr(self, key, self.typemap[key](val))
                except KeyError, e:
                    raise KeyError, \
                        "This key is not in the typemap of %s: %s" \
                        % (type(self), str(e))
            
            # all fields must be filled
            for slot in self.__slots__:
                if not hasattr(self, slot):
                    raise ValueError, \
                        "Not enough metadata supplied; missing %s" % slot
        else:
            raise NotImplementedError
    
    def __str__(self):
        return str(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))
    
    def __repr__(self):
        return repr(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))
    
    def __ior__(self, other):
        """
        Merge metadata; this must be subclassed.
        """
        raise NotImplementedError
    
    def __or__(self, other):
        return self.copy().__ior__(other)
    
    def __ror__(self, other):
        return self.copy().__ior__(other)
    
    def __eq__(self, other):
        for slot in self.__slots__:
            if getattr(self, slot) != getattr(other, slot):
                return False
        return True
    
    def copy(self):
        return type(self)(dict([(slot, getattr(self, slot)) for slot \
            in self.__slots__ if slot != "typemap"]))
    
    def todict(self):
        return dict([(key, getattr(self, key)) for key in self.__slots__ \
            if key != "typemap"])

class MetaArray(numpy.ndarray):
    """
    An array containing a data and metadata.  Intended to be
    subclassed.
    
    On b = MetaArray(a) where a is a MetaArray, metadata is copied, but data
    is not.
    """
    __array_priority__ = 10.1  # ufuncs mixing array types return MetaArray
    _metadata_type = Metadata
    
    def __new__(subtype, data=None, metadata=None, dtype=None, copy=False,
        subok=True):
        # Case 1: data is an MetaArray, it has metadata already.
        if isinstance(data, MetaArray):
            if dtype is None:
                dtype = data.dtype
            else:
                dtype = numpy.dtype(dtype)
            if not copy and dtype==data.dtype and metadata is None:
                return data
            elif metadata is not None:
                # share memory, but not metadata
                new = numpy.array(data, dtype=dtype, copy=copy, subok=True)
                new.metadata = subtype._metadata_type(metadata)
                return new
            else:
                # copy, update metadata, then return
                new = data.astype(dtype)  # always copies
                new.metadata = data.metadata
                new._baseclass = data._baseclass
                return new
        
        # All other cases, create a new array and attach metadata
        # Unless you specified otherwise, we'll reuse memory from existing
        # arrays.
        new = numpy.array(data, dtype=dtype, copy=copy, subok=True)
        _baseclass = type(new)
        new = new.view(subtype)
        new.metadata = subtype._metadata_type(metadata)
        new._baseclass = _baseclass
        return new
    
    def __array_finalize__(self, obj):
        """
        Called anytime a MetaArray is returned; make sure that metadata
        is set to something.
        """
        self.metadata = getattr(obj, "metadata", None)
        self._baseclass = getattr(obj, "_baseclass", type(obj))
    
    def __array_wrap__(self, obj):
        """
        Called anytime a ufunc operates on a MetaArray and another object.
        The result of the ufunc is obj.  The MetaArray operand is self.
        """
        result = obj.view(type(self))
        result.metadata = self.metadata.copy()
        return result
    
    def __repr__(self):
        return "%s, %s)" % (repr(self.view(numpy.ndarray))[:-1], repr(self.metadata))
    
    def __str__(self):
        return "%s %s" % (str(self.view(numpy.ndarray)), str(self.metadata))
    
    # methods that return an array, wrapped to return a MetaArray
    __abs__ = _arraymethod('__abs__')
    __add__ = _arraymethod('__add__')
    __and__ = _arraymethod('__and__')
    __copy__ = _arraymethod('__copy__')
    __deepcopy__ = _arraymethod('__deepcopy__')
    __div__ = _arraymethod('__div__')
    __divmod__ = _arraymethod('__divmod__')
    __floordiv__ = _arraymethod('__floordiv__')
    __hex__ = _arraymethod('__hex__')
    __iadd__ = _arraymethod('__iadd__')
    __iand__ = _arraymethod('__iand__')
    __idiv__ = _arraymethod('__idiv__')
    __ifloordiv__ = _arraymethod('__ifloordiv__')
    __ilshift__ = _arraymethod('__ilshift__')
    __imod__ = _arraymethod('__imod__')
    __imul__ = _arraymethod('__imul__')
    __invert__ = _arraymethod('__invert__')
    __ior__ = _arraymethod('__ior__')
    __ipow__ = _arraymethod('__ipow__')
    __irshift__ = _arraymethod('__irshift__')
    __isub__ = _arraymethod('__isub__')
    __itruediv__ = _arraymethod('__itruediv__')
    __ixor__ = _arraymethod('__ixor__')
    __lshift__ = _arraymethod('__lshift__')
    __mul__ = _arraymethod('__mul__')
    __rmod__ = _arraymethod('__rmod__')
    __rmul__ = _arraymethod('__rmul__')
    __ror__ = _arraymethod('__ror__')
    __rpow__ = _arraymethod('__rpow__')
    __rrshift__ = _arraymethod('__rrshift__')
    __rshift__ = _arraymethod('__rshift__')
    __rsub__ = _arraymethod('__rsub__')
    __rtruediv__ = _arraymethod('__rtruediv__')
    __rxor__ = _arraymethod('__rxor__')
    __sub__ = _arraymethod('__sub__')
    __truediv__ = _arraymethod('__truediv__')
    __xor__ = _arraymethod('__xor__')
    astype = _arraymethod('astype')
    byteswap = _arraymethod('byteswap')
    choose = _arraymethod('choose')
    clip = _arraymethod('clip')
    compress = _arraymethod('compress')
    conj = _arraymethod('conj')
    conjugate = _arraymethod('conjugate')
    copy = _arraymethod('copy')
    cumprod = _arraymethod('cumprod')
    cumsum = _arraymethod('cumsum')
    diagonal = _arraymethod('diagonal')
    fill = _arraymethod('fill')
    flat = _arraymethod('flat')
    flatten = _arraymethod('flatten')
    repeat = _arraymethod('repeat')
    squeeze = _arraymethod('squeeze')
    transpose = _arraymethod('transpose')
    
    T = property(fget = lambda self: self.transpose())
    H = property(fget = lambda self: self.T.conj()) # Hermitian transpose
    
    def _get_data(self):
        return self.view(self._baseclass)
    A = property(fget=_get_data) # get at the underlying Array data
    
    # Pickling
    def __getstate__(self):
        """
        Returns the internal state of the object, for pickling purposes.
        """
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self.A.tostring(),
                 self.metadata.todict(),
                 )
        return state
    
    def __setstate__(self, state):
        """
        Restores the internal state of the masked array, for unpickling
        purposes. `state` is typically the output of the ``__getstate__``
        output, and is a 5-tuple:
        
          - class name
          - a tuple giving the shape of the data
          - a typecode for the data
          - a binary string for the data
          - a binary string for the mask.
        """
        (ver, shp, typ, isf, raw, meta) = state
        if ver != 1: raise NotImplementedError
        numpy.ndarray.__setstate__(self, (shp, typ, isf, raw))
        self.metadata = self._metadata_type(meta)
    
    def __reduce__(self):
        """
        Returns a 3-tuple for pickling a MetaArray
          - reconstruction function
          - tuple to pass reconstruction function
          - state, which will be passed to __setstate__
        """
        return (_mareconstruct,
                (self.__class__, self._baseclass, (0,), 'b', ),
                self.__getstate__())
    
def _mareconstruct(subtype, baseclass, baseshape, basetype):
    """
    Internal function that builds a new MaskedArray from the information
    stored in a pickle.
    """
    return subtype.__new__(subtype, [], dtype=basetype)

class MetaArrayList(list):
    _itemtype = MetaArray
    
    # these methods will act upon each element of this list
    __abs__ = _elementwise_method('__abs__')
    __add__ = _elementwise_method('__add__')
    __and__ = _elementwise_method('__and__')
    __copy__ = _elementwise_method('__copy__')
    __deepcopy__ = _elementwise_method('__deepcopy__')
    __div__ = _elementwise_method('__div__')
    __divmod__ = _elementwise_method('__divmod__')
    __floordiv__ = _elementwise_method('__floordiv__')
    __hex__ = _elementwise_method('__hex__')
    __iadd__ = _elementwise_method('__iadd__')
    __iand__ = _elementwise_method('__iand__')
    __idiv__ = _elementwise_method('__idiv__')
    __ifloordiv__ = _elementwise_method('__ifloordiv__')
    __ilshift__ = _elementwise_method('__ilshift__')
    __imod__ = _elementwise_method('__imod__')
    __imul__ = _elementwise_method('__imul__')
    __invert__ = _elementwise_method('__invert__')
    __ior__ = _elementwise_method('__ior__')
    __ipow__ = _elementwise_method('__ipow__')
    __irshift__ = _elementwise_method('__irshift__')
    __isub__ = _elementwise_method('__isub__')
    __itruediv__ = _elementwise_method('__itruediv__')
    __ixor__ = _elementwise_method('__ixor__')
    __lshift__ = _elementwise_method('__lshift__')
    __mul__ = _elementwise_method('__mul__')
    __rmod__ = _elementwise_method('__rmod__')
    __rmul__ = _elementwise_method('__rmul__')
    __ror__ = _elementwise_method('__ror__')
    __rpow__ = _elementwise_method('__rpow__')
    __rrshift__ = _elementwise_method('__rrshift__')
    __rshift__ = _elementwise_method('__rshift__')
    __rsub__ = _elementwise_method('__rsub__')
    __rtruediv__ = _elementwise_method('__rtruediv__')
    __rxor__ = _elementwise_method('__rxor__')
    __sub__ = _elementwise_method('__sub__')
    __truediv__ = _elementwise_method('__truediv__')
    __xor__ = _elementwise_method('__xor__')
    astype = _elementwise_method('astype')
    byteswap = _elementwise_method('byteswap')
    choose = _elementwise_method('choose')
    clip = _elementwise_method('clip')
    compress = _elementwise_method('compress')
    conj = _elementwise_method('conj')
    conjugate = _elementwise_method('conjugate')
    copy = _elementwise_method('copy')
    cumprod = _elementwise_method('cumprod')
    cumsum = _elementwise_method('cumsum')
    diagonal = _elementwise_method('diagonal')
    fill = _elementwise_method('fill')
    flat = _elementwise_method('flat')
    flatten = _elementwise_method('flatten')
    repeat = _elementwise_method('repeat')
    squeeze = _elementwise_method('squeeze')
    sum = _elementwise_method('sum')
    transpose = _elementwise_method('transpose')
    
    T = property(fget = lambda self: self.transpose()) # Transpose
    H = property(fget = lambda self: self.T.conj()) # Hermitian transpose
    
    # special case a few attribute accessors
    def _get_real(self):
        return type(self)([x.real for x in self])
    real = property(fget=_get_real)
    def _get_imag(self):
        return type(self)([x.imag for x in self])
    imag = property(fget=_get_imag)
    
    def _get_data(self):
        """
        Return a regular list of arrays.
        """
        return [x.view(x._baseclass) for x in self]
    A = property(fget=_get_data)


##############################################################################
# Define TimeSeries and associated structures
##############################################################################

class TimeSeriesMetadata(Metadata):
    """
    Hold the metadata associated with a spectrum, including frequency
    resolution, a segmentlist indicating what times were involved in taking
    the spectrum, a channel name, etc.
    """
    # specify all allowed metadata here; type is required
    typemap = {"name": str,
               "dt": float,
               "segments": segments.segmentlist,
               "comments": list}
    __slots__ = typemap.keys() + ['typemap']
    
    def __ior__(self, other):
        """
        Merge metadata.  No repeats.  Throw error on incompatible spectra.
        Let None act as identity.
        """
        if other is None:
            return self
        
        # check that metadata are compatible for merging
        assert numpy.alltrue([getattr(self, attr) == getattr(other, attr) \
            for attr in self.__slots__ \
            if attr not in ('segments', 'comments', 'typemap')])
        
        # add, but do not join, segments
        self.segments.extend([seg for seg in other.segments \
            if seg not in self.segments])
        self.segments.sort()
        
        # add only new comments
        self.comments.extend([comment for comment in other.comments \
            if comment not in self.comments])
        return self

class TimeSeries(MetaArray):
    """
    This is a MetaArray, but with the metadata typemap specified.
    """
    _metadata_type = TimeSeriesMetadata
    
    def ordinates(self):
        """
        Return an single-precision ndarray containing the times at
        which this TimeSeries is sampled.
        """
        m = self.metadata
        segs = m.segments
        segs.coalesce()
        
        # start with a uniformly spaced domain
        result = numpy.arange(len(self), dtype=numpy.float32)*m.dt
        
        # split domain by applying different offsets for each segment
        offsets = [seg[0] for seg in segs]
        seg_lengths = [int(abs(seg)/m.dt) for seg in segs]
        assert sum(seg_lengths) == len(self)
        boundaries = [0]+[sum(seg_lengths[:i]) for i in range(1, len(seg_lengths)+1)]
        assert boundaries[-1] == len(self)
        slices = [slice(a, b) for a,b in zip(boundaries[:-1], boundaries[1:])]
        
        for sl, offset in zip(slices, offsets):
            result[sl] += offset
        return result

class TimeSeriesList(MetaArrayList):
    """
    A list of TimeSeries.
    """
    _itemtype = TimeSeries

    def segments(self):
        """
        Return the (uncoalesced) list of segments represented by the TimeSeries
        in this TimeSeriesList.
        """
        segs = segments.segmentlist()
        for series in self:
            segs.extend(series.metadata.segments)
        return segs

    def merge_list(self):
        """
        Concatenate the list into one single TimeSeries.
        """
        meta = reduce(operator.or_, [ts.metadata for ts in self])
        return TimeSeries(numpy.concatenate(self.A), meta)

##############################################################################
# Define Spectrum and associated structures
##############################################################################

class SpectrumMetadata(Metadata):
    """
    Hold the metadata associated with a spectrum, including frequency
    resolution, a segmentlist indicating what times were involved in taking
    the spectrum, a channel name, etc.
    """
    # specify all allowed metadata here; type is required
    typemap = {"name": str,
               "df": float,
               "f_low": float,
               "segments": segments.segmentlist,
               "comments": list}
    __slots__ = typemap.keys() + ['typemap']
    
    def __ior__(self, other):
        """
        Merge metadata.  No repeats.  Throw error on incompatible spectra.
        Let None act as identity.
        """
        if other is None:
            return self
        if self is None:
            return other
        
        # check that metadata are compatible for merging
        assert numpy.alltrue([getattr(self, attr) == getattr(other, attr) \
            for attr in ("df", "f_low")])
        
        # add, but do not join segments
        self.segments.extend([seg for seg in other.segments \
            if seg not in self.segments])
        self.segments.sort()
        
        # add only new comments
        self.comments.extend([comment for comment in other.comments \
            if comment not in self.comments])
        return self

class Spectrum(MetaArray):
    """
    This is a MetaArray, but with the metadata typemap specified.
    """
    _metadata_type = SpectrumMetadata
    
    def ordinates(self):
        """
        Return an single-precision ndarray containing the frequencies at
        which this Spectrum is sampled.
        """
        m = self.metadata
        return numpy.arange(len(self), dtype=numpy.float32)*m.df + m.f_low

class SpectrumList(MetaArrayList):
    """
    A list of Spectra.
    """
    _itemtype = Spectrum
    
    def segments(self):
        """
        Return the (uncoalesced) list of segments represented by the Spectra
        in this SpectrumList.
        """
        segs = segments.segmentlist()
        for spectrum in self:
            segs.extend(spectrum.metadata.segments)
        return segs

    def ordinates(self):
        """
        Return an single-precision ndarray containing the frequencies at
        which these Spectrums are sampled (they are required to match).
        For an empty list, return a zero-length array.
        """
        if len(self) == 0:
            return numpy.array([], dtype=numpy.float32)
        m = self[0].metadata
        return numpy.arange(len(self[0]), dtype=numpy.float32)*m.df + m.f_low
    
    def sum_spectra(self):
        """
        Sum across Spectrums, returning a single Spectrum
        """
        if len(self) == 0:
            raise ValueError
        meta = reduce(operator.or_, [s.metadata for s in self])
        return Spectrum(numpy.sum(self.A, axis=0), meta)

class SpectrumDict(dict):
    """
    A dictionary allowing access to FFTs of different data streams and
    providing convenient batch operations.
    """
    pass
