import sys
import lalinference as li
from ctypes import c_char_p
import collections

class LIVariablesWrap(collections.MutableMapping):
    def __init__(self,init=None):
        """
        Wrapper for dictionaries to pass to LALInference functions
        init : dictionary type to set variables with
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
        if not li.CheckVariable(self.v, key):
            raise KeyError(key)
        return li.GetVariableVaryType(self.v, key)
    def type(self, key):
        if not li.CheckVariable(self.v, key):
            raise KeyError(key)
        return li.GetVariableType(self.v, key)

class _variterator(object):
      def __init__(self, var):
          self.varitem = var.head
      def __next__(self):
          if not self.varitem:
              raise StopIteration
          else:
              this = self.varitem
              self.varitem=self.varitem.next
              return(this.name)
      def next(self):
          return self.__next__()


class LALInferenceWrapper(object):
    def __init__(self, argv):
        procParams=li.ParseCommandLine(len(argv),argv)
        self.state = li.InitRunState(procParams)
        self.state.commandLine=procParams
        li.InitCBCThreads(self.state,1)

        li.InjectInspiralSignal(self.state.data, self.state.commandLine)
        li.ApplyCalibrationErrors(self.state.data, procParams)
        li.InitCBCPrior(self.state)
        li.InitLikelihood(self.state)
        li.InitCBCThreads(self.state,1)

    
    def log_likelihood(self,params):
        # Pick up non-sampled vars
        liv = LIVariablesWrap(self.state.threads.currentParams)
        # Update with proposed values
        liv.update(params)
        print(liv)
        self.state.threads.model.currentParams=liv.v
        print(self.state.data,self.state.threads.model)
        return li.MarginalisedPhaseLogLikelihood(liv.v, self.state.data, self.state.threads.model)
    
    def log_prior(self,params):
        # Pick up non-sampled vars
        liv = LIVariablesWrap(self.state.threads.currentParams)
        # Update with proposed values
        liv.update(params)
        return li.InspiralPrior(self.state, liv.v, self.state.threads.model)

    def params(self):
        LIV=LIVariablesWrap(self.state.threads.currentParams)
        return LIV.keys()
    
    def sampling_params(self):
        pars = LIVariablesWrap(self.state.threads.currentParams)
        return [p for p in pars if pars.varyType(p)==li.LALINFERENCE_PARAM_LINEAR
                  or pars.varyType(p)==li.LALINFERENCE_PARAM_CIRCULAR
                  ]
    
    def prior_bounds(self):
        bounds=[]
        libounds = LIVariablesWrap(self.state.priorArgs)
        for p in self.sampling_params():
            try:
                low = libounds[p+'_min']
                high = libounds[p+'_max']
                bounds.append((low, high))
            except KeyError:
                pass
        return bounds

import cpnest.model
import traceback

class LIModel(cpnest.model.Model):
    def __init__(self, *args, **kwargs):
        super(LIModel, self).__init__()
        self.limodel = LALInferenceWrapper(sys.argv)

        self.names = self.limodel.sampling_params()
        self.bounds = self.limodel.prior_bounds()
        print('Sampling in {0}'.format(self.names))
        print('Bounds: {0}'.format(self.bounds))
    def log_likelihood(self, x):
        print('Computing logL: '+str(x))
        traceback.print_stack()
        logl=self.limodel.log_likelihood(x)
        print(logl)
        return logl
    def log_prior(self, x):
        print('Computing prior: '+str(x))
        traceback.print_stack()
        logp=self.limodel.log_prior(x)
        print(logp)
        return logp


if __name__=='__main__':
    LIstate = LIModel(sys.argv)
    nest=cpnest.CPNest(LIstate, Nlive=100, Nthreads=1, verbose=1)
    print(li.MarginalisedPhaseLogLikelihood(nest.user.limodel.state.threads.model.params,nest.user.limodel.state.data, nest.user.limodel.state.threads.model))
    nest.run()
    
    
#if __name__=='__main__':
#    main()
    
