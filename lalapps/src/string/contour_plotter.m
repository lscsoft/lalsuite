clear
load gamma.dat

T=365*24*3600;

epsi=gamma(:,3);
Gmu=gamma(:,4);
rate=gamma(:,5);

figure;loglog(10^-8,10^-8); hold on
p=logspace(-4,0, 40);

for i=1:size(gamma,1)
 
     if(gamma(i,5)*T > -log(1-.9))
        plot(gamma(i,4),gamma(i,3), 'ro')
     else
        plot(gamma(i,4),gamma(i,3), 'go')        
    end
 
end

Gmu2 = unique(gamma(:,4))';
eps2 = unique(gamma(:,3));
Gmumat2 = kron(Gmu2,ones(size(eps2)));
epsmat2 = kron(eps2,ones(size(Gmu2)));
rate = zeros(length(eps2),length(Gmu2));
for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = gamma(:,4) == Gmu2(kk) & gamma(:,3) == eps2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5);
    else
      disp('PROBLEM')
    end
  end
end

contour(Gmumat2,epsmat2,rate*T, [-log(1-0.997) -log(1-0.95)  -log(1-0.68)],'r*');



