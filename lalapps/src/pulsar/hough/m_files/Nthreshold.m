%  $Id$
% script to compute statistic values

 
 alpha =exp(-1.6);
 %alphaHough = 10^(-9);
 alphaHough = 10^(-10);
 
 N= 1886
 %N=671
 
 n_mean= N*alpha
 n_sigma = sqrt(N*alpha*(1-alpha) )
 
 n_th= n_mean +sqrt(2)*n_sigma* erfcinv(2*alphaHough)



 %  y=binopdf(HMhisto(:,1),M,a);
