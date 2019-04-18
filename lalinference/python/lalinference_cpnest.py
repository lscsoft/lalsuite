import sys
import argparse
from lalinference.wrapper import LALInferenceCBCWrapper
import cpnest.model

class LIModel(cpnest.model.Model):
    def __init__(self, *args, **kwargs):
        super(LIModel, self).__init__()
        self.limodel = LALInferenceCBCWrapper(sys.argv)

        self.names = self.limodel.sampling_params()
        bounds_dict = self.limodel.prior_bounds()
        self.bounds = [bounds_dict[p] for p in self.names]
        print('Sampling in {0}'.format(self.names))
        print('Bounds: {0}'.format(self.bounds))
    def log_likelihood(self, x):
        logl=self.limodel.log_likelihood(x)
        return logl
    def log_prior(self, x):
        logp=self.limodel.log_prior(x)
        return logp

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Nested sampling for CBC analysis')
    parser.add_argument('--nlive',type=int,default=1000)
    parser.add_argument('--nthreads',type=int,default=1)
    parser.add_argument('--verbose',action='store_true',default=False)
    parser.add_argument('--outfile',required=True)
    parser.add_argument('--plot',default=False,const=True,nargs='?')
    parser.add_argument('--maxmcmc',default=5000,type=int)
    parser.add_argument('--poolsize',default=500,type=int)
    opts, args = parser.parse_known_args(sys.argv)
    print(args)
    LIstate = LIModel(sys.argv)
    nest=cpnest.CPNest(LIstate, Nlive=opts.nlive, Nthreads=opts.nthreads, verbose=opts.verbose, maxmcmc=opts.maxmcmc, Poolsize=opts.poolsize)
    nest.run()
    if opts.plot:
        nest.plot()

    
#if __name__=='__main__':
#    main()
    
