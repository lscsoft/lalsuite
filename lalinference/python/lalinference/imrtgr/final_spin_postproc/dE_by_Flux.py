import numpy as np 
from numpy import log 
GAMMA = 0.577216
pi = np.pi

# Eqs.(3.2) of Ajith (2011) http://arxiv.org/pdf/1107.1267v2.pdf

""" compute the re-exapnded dEnergy/flux """ 
def denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, order):
    
	# different powers of v 
	v2 = v*v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v; v7 = v6*v; v9 = v7*v2 

	# initialize the cofficients
	dEbF0 = dEbF2 = dEbF3 = dEbF4 = dEbF5 = dEbF6 = dEbF6L = dEbF7 = 0.

	if order >= 0:
		dEbF0 = -5./(32.*eta)
	if order >= 2:
		dEbF2 = 2.2113095238095237 + (11*eta)/4. 
	if order >= 3:
		dEbF3 = (113*chiadL*delta)/12. + chisdL*(9.416666666666666 - (19.*eta)/3.) - 4.*pi 
	if order >= 4:
		dEbF4 = 3.010315295099521 + (233*chisdchia*delta)/48. - (719.*chiadL*chisdL*delta)/48. + \
			chiasqr*(2.4270833333333335 - 10.*eta) + pow(chisdL,2.)*(-7.489583333333333 - eta/24.) + \
			chissqr*(2.4270833333333335 + (7.*eta)/24.) + (5429.*eta)/1008. + (617*pow(eta,2.))/144. + \
			pow(chiadL,2.)*(-7.489583333333333 + 30.*eta) 
	if order >= 5:
			dEbF5 = chiadL*delta*(72.71676587301587 + (7*eta)/2.) + chisdL*(72.71676587301587 - \
						(1213*eta)/18. - (17*pow(eta,2))/2.) - (7729*pi)/672. + (13*eta*pi)/8. 
	if order >= 6:
		dEbF6 = -115.2253249962622 - (15211*pow(eta,2))/6912. + (25565*pow(eta,3))/5184. + \
					(32*pow(pi,2))/3. + eta*(258.1491854023631 - (451*pow(pi,2))/48.) + (1712*GAMMA)/105. 
		dEbF6L = 1712./105. 
	if order >= 7:
		dEbF7 = (-15419335.*pi)/1.016064e6 - (75703.*eta*pi)/6048. + (14809.*pow(eta,2)*pi)/3024. 

	return (dEbF0/v9)*(1. + dEbF2*v2 + dEbF3*v3 + dEbF4*v4 + dEbF5*v5 + (dEbF6+dEbF6L*log(4.*v))*v6 + dEbF7*v7) 
