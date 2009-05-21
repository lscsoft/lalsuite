import numpy
import pylab

def posterior(VT, sigma, Lambda):
        '''
        This function implements the analytic marginalization in 
        Biswas, Creighton, Brady, Fairhurst (25)
        Cause hey, why not?
        This takes arrays of VT, sigma and Lambda to combine.
        '''

	length = 100000
        #FIXME does this need to be an integer???
        K = (VT / sigma)**2.0

        #FIXME, drew said this was cool?
        mu = numpy.arange(length) * 100.0 / VT.sum() / length

        #FIXME this is just ones with the right dtype
        post = numpy.ones(len(mu), dtype="float")

        for vt, k, lam in zip(VT, K, Lambda):
                post *= vt / (1.0 + lam) * ( (1.0 + mu * vt / k)**(-k-1) + (mu * vt * lam * (1.0 + 1.0/k) /(1.0 + mu * vt / k)**(k+2)) )

        return mu, post

def integrate_posterior(mu, post, conf):
	cumpost = post.cumsum()/post.sum()
	val = [idx for idx in range(len(cumpost)) if cumpost[idx] >= conf][0]
	return mu[val]

# test case

VT = numpy.array([10.0**8])
sigma = numpy.array([2.0*10**7])
Lambda = numpy.array([1.0])

mu, post = posterior(VT, sigma, Lambda)
#pylab.semilogy(mu, post.cumsum()/post.sum())
#pylab.show()
print integrate_posterior(mu, post, 0.90)
