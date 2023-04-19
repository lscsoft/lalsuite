import lal
from functools import wraps
from astropy import units as u


##############################################################################################################
# Mapping between LAL Errors and Python
##############################################################################################################

def mapexception(func):
    """
    Wrapper to map LAL erros to Python errors
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        exception = None
        try:
            return func(*args, **kwargs)
        except RuntimeError as e:
            lalerr = lal.GetBaseErrno()
            exception = exceptionmap[lalerr] if lalerr in exceptionmap else e
            raise exception from None
    return wrapper

# Map of errors raised by LAL to Python errors
exceptionmap = {
    lal.ENOENT:     FileNotFoundError("No such file or directory"),
    lal.EIO:        OSError("I/O Error"),
    lal.ENOMEM:     MemoryError("Memory Allocation Error"),
    lal.EFAULT:     ValueError("Invalid pointer"),
    lal.EINVAL:     ValueError("Invalid Argument"),
    lal.EDOM:       ValueError("Input domain error"),
    lal.ERANGE:     ArithmeticError("Output range error"),
    lal.ENOSYS:     NotImplementedError("Function not implemented"),
    lal.EFAILED:    RuntimeError("Generic Failure"),
    lal.EBADLEN:    ValueError("Inconsistent or invalid length"),
    lal.ESIZE:      ValueError("Wrong Size"),
    lal.EDIMS:      ValueError("Wrong dimensions"),
    lal.ETYPE:      TypeError("Wrong or unknown type"),
    lal.ETIME:      ValueError("Invalid time"),
    lal.EFREQ:      ValueError("Invalid frequency"),
    lal.EUNIT:      u.UnitsError("Invalid units"),
    lal.ENAME:      ValueError("Wrong name"),
    lal.EDATA:      ValueError("Invalid data"),
    lal.ESYS:       SystemError("System Error"),
    lal.EERR:        RuntimeError("Internal Error"),
    lal.EFPINVAL:   FloatingPointError("IEEE Invalid floating point operation, eg sqrt(-1), 0/0"),
    lal.EFPDIV0:    ZeroDivisionError("IEEE Division by zero floating point error"),
    lal.EFPOVRFLW:  OverflowError("IEEE Floating point overflow error"),
    lal.EFPUNDFLW:  FloatingPointError("IEEE Floating point underflow error"),
    lal.EFPINEXCT:  FloatingPointError("IEEE Floating point inexact error"),
    lal.EMAXITER:   ArithmeticError("Exceeded maximum number of iterations"),
    lal.EDIVERGE:   ArithmeticError("Series is diverging"),
    lal.ESING:      ArithmeticError("Apparent singularity detected"),
    lal.ETOL:       ArithmeticError("Failed to reach specified tolerance"),
    lal.ELOSS:      ArithmeticError("Loss of accuracy"),
    lal.EUSR0:      RuntimeError("User defined Error"),
    lal.EUSR1:      RuntimeError("User defined Error"),
    lal.EUSR2:      RuntimeError("User defined Error"),
    lal.EUSR3:      RuntimeError("User defined Error"),
    lal.EUSR4:      RuntimeError("User defined Error"),
    lal.EUSR5:      RuntimeError("User defined Error"),
    lal.EUSR6:      RuntimeError("User defined Error"),
    lal.EUSR7:      RuntimeError("User defined Error"),
    lal.EUSR8:      RuntimeError("User defined Error"),
    lal.EUSR9:      RuntimeError("User defined Error")
    }

##############################################################################################################
# Errors callable by GWSignal.
# Each error should be subclasses from existing python errors, for eg., ValueError, Exception etc.
##############################################################################################################

class WaveformGenerationError(Exception):
    """
    Raise Error when there is an issue with waveform generation.

    Parameters
    ----------
        message   : Output message for the error.
    Returns
    -------
        Raises waveform generation error exception

    """
    def __init__(self, message="Error during waveform generation"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}'


class DomainError(ValueError):
    """
    Raise Error when parameter passes is out of waveform domain

    Parameters
    ----------
        parameter : Parameter in question
        message   : Output message for the error.

    Returns
    -------
        Raises domain error exception
    """
    def __init__(self, parameter, message="Input parameter out of waveform domain"):
        self.parameter = parameter
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.parameter} -> {self.message}'


class PathError(FileNotFoundError):
    """
    Raise Exception when the needed file is not at the given path

    Parameters
    ----------
        file_name : Path to the file in question.
        message   : Output message for the error.

    Returns
    -------
        Raises File not found at given path error excetion.

    """
    def __init__(self, file_path, message="File not found at "):
        self.file_path = file_path
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} -> {self.file_path}'
