import lal
import lalinference as li
import collections.abc

class LIVariablesWrap(collections.abc.MutableMapping):
    def __init__(self,init=None):
        """
        Wrapper to present a LALInferenceVariable as a dict.

        Parameters
        ----------
        init : dict
          Initialise with the given dictionary.
          If init is itself a LALInferenceVariables C struct
          then the wrapper will wrap around it but not reallocate memory
        """
        self.owner=True # Python should manage the object's memory
        if isinstance(init,li.Variables):
            self.v=init
            self.owner=False
        else:
            self.v=li.Variables()
            if init:
                self.update(init)
    def __delitem__(self,key):
        if li.CheckVariable(self.v, key):
            li.RemoveVariable(self.v, key)
        else:
            raise KeyError(key)
    def __setitem__(self, key, value):
        if type(value)==float:
            li.AddREAL8Variable(self.v, key, value, li.LALINFERENCE_PARAM_LINEAR)
        elif type(value)==int:
            li.AddINT4Variable(self.v, key, value, li.LAINFERENCE_PARAM_LINEAR)
        else:
            raise TypeError('Unsupported type: ',key, self.type(key))
    def __getitem__(self, key):
        if li.CheckVariable(self.v, key):
            if self.type(key)==li.LALINFERENCE_REAL8_t:
                return li.GetREAL8Variable(self.v, key)
            elif self.type(key)==li.LALINFERENCE_INT4_t:
                return li.GetINT4Variable(self.v, key)
            elif self.type(key)==li.LALINFERENCE_UINT4_t:
                return li.GetUINT4Variable(self.v, key)
            else:
                raise(TypeError('Unsupported type: ',key,self.type(key)))
        else:
            raise KeyError(key)
    def __iter__(self):
        return _variterator(self.v)
    def __len__(self):
        return self.v.dimension
    def __del__(self):
        if self.owner:
            li.ClearVariables(self.v)
    def __repr__(self):
        return 'LIVariablesWrap('+repr(dict(self))+')'
    def __str__(self):
        return str(dict(self))
    def varyType(self, key):
        """
        Return the lalinference variable's varyType

        Parameters
        ----------
        key : str
            The name of the variable to look up

        Returns
        -------
        varytype : lalinference.varyType (e.g. lalinference.LALINFERENCE_PARAM_FIXED)
        """
        if not li.CheckVariable(self.v, key):
            raise KeyError(key)
        return li.GetVariableVaryType(self.v, key)
    def type(self, key):
        """
        Return the lalinference variable's varyType

        Parameters
        ----------
        key : str
            The name of the variable to look up

        Returns
        -------
        type : the LALInference type (e.g. lalinference.LALINFERENCE_REAL8_t)
        """
        if not li.CheckVariable(self.v, key):
            raise KeyError(key)
        return li.GetVariableType(self.v, key)

class _variterator(object):
      def __init__(self, var):
          self.varitem = var.head

      def __iter__(self):
          return self

      def __next__(self):
          if not self.varitem:
              raise StopIteration
          else:
              this = self.varitem
              self.varitem=self.varitem.next
              return(this.name)

      def next(self):
          return self.__next__()


class LALInferenceCBCWrapper(object):
    """
    Class to wrap a LALInference CBC analysis
    state, and expose the likelihood and prior
    methods to python programs
    """
    def __init__(self, argv):
        """
        Parameters
        ----------
        argv : list
            List of command line arguments that will be used to
            set up the lalinference state. (similar to argv)
        """
        strvec = lal.CreateStringVector(argv[0])
        for a in argv[1:]:
            strvec=lal.AppendString2Vector(strvec, a)
        procParams=li.ParseStringVector(strvec)
        self.state = li.InitRunState(procParams)
        self.state.commandLine=procParams
        li.InitCBCThreads(self.state,1)

        # This is what Nest does
        li.InjectInspiralSignal(self.state.data, self.state.commandLine)
        li.ApplyCalibrationErrors(self.state.data, procParams)
        if li.GetProcParamVal(procParams,'--roqtime_steps'):
            li.SetupROQdata(self.state.data, procParams)
        li.InitCBCPrior(self.state)
        li.InitLikelihood(self.state)
        li.InitCBCThreads(self.state,1)


    def log_likelihood(self,params):
        """
        Log-likelihood function from LALInference

        Parameters
        ----------
        params : dict
            Dict-like object of sampling parameters, will
            be automatically converted for lalinference

        Returns
        -------
        logL : float
            log-likelihood value
        """
        # Pick up non-sampled vars
        liv = LIVariablesWrap(self.state.threads.currentParams)
        # Update with proposed values
        liv.update(params)
        self.state.threads.model.currentParams=liv.v
        return li.MarginalisedPhaseLogLikelihood(liv.v, self.state.data, self.state.threads.model)

    def log_prior(self,params):
        """
        Log-prior function from LALInference

        Parameters
        ----------
        params : dict
            Dict-like object of sampling parameters, will
            be automatically converted for lalinference

        Returns
        -------
        logPr : float
            log-prior value
        """
        # Pick up non-sampled vars
        liv = LIVariablesWrap(self.state.threads.currentParams)
        # Update with proposed values
        liv.update(params)
        return li.InspiralPrior(self.state, liv.v, self.state.threads.model)

    def params(self):
        """
        Parameter names from the LALInference model. Includes
        those which are fixed

        Returns
        names : list
            A list of parameter names
        """
        LIV=LIVariablesWrap(self.state.threads.currentParams)
        return LIV.keys()

    def sampling_params(self):
        """
        Parameter names from the LALInference model. Includes
        only those which are varied in the sampling.

        Returns
        names : list
            A list of parameter names
        """
        pars = LIVariablesWrap(self.state.threads.currentParams)
        return [p for p in pars if pars.varyType(p)==li.LALINFERENCE_PARAM_LINEAR
                  or pars.varyType(p)==li.LALINFERENCE_PARAM_CIRCULAR
                  ]

    def prior_bounds(self):
        """
        Bounds of the sampling parameters.

        Returns
        bounds : dict
            A dict of (low,high) pairs, indexed by parameter name
            e.g. {'declination' : (0, 3.14159), ...}
        """
        bounds={}
        libounds = LIVariablesWrap(self.state.priorArgs)
        for p in self.sampling_params():
            try:
                low = libounds[p+'_min']
                high = libounds[p+'_max']
                bounds[p]=(low, high)
            except KeyError:
                pass
        return bounds
