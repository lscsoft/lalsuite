path = '/home/jsylvest/S2/S2v2.2_opti/'

Tau = 0.001;

t=0:1/16384:1-1/16384;
str = '';

for j=1:length(Tau)
tau = Tau(j);
str = [str 'G' num2str(tau) ','];

h0 = 1;
hp = h0*exp(-(t-0.5).^2/tau^2); 

ilwd_write([path 'G' num2str(tau) '_p.ilwd'],['g' num2str(tau)],['g' num2str(tau)],hp);

end


hc = zeros(size(t));
ilwd_write([path 'Zeros.ilwd'],'zeros','zeros',hc);
