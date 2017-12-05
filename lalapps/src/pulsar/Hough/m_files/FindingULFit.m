
MC=load('MultiMC_100_1000.25_coarse');
fmin=MC(:,1);
fmax=MC(:,2);
signi=MC(:,3);
hmin=MC(:,24);
hmax=MC(:,25);
UL=MC(:,26);

%plot( (fmax+fmin)/2, UL, 'b-o')

S=signi/sqrt(2) + erfcinv(0.1);
Nbands = length(S);

kk=load('MultiFit.txt');
xx=kk(1:Nbands ,2)*2^(0.25)/1800;
x=sqrt(S).*xx;

indices=find(UL<2.5e-23); % to avoind many non-converging MC-UL

f1=polyfit(x(indices),UL(indices),1)
m1=mean(UL(indices)./x(indices))
med2=median(UL(indices)./x(indices))
med1=median (UL./x)
std1=std(UL(indices)./x(indices))

plot((fmax(indices)+fmin(indices))/2, UL(indices)./x(indices), '.')

for k=1:9
  v=(k-1)*400+1:k*400;
  clear indices
  indices=find(UL(v)<2.5e-23);
  fv=polyfit(x(v(indices)),UL(v(indices)),1)
  mv=mean(UL(v(indices))./x(v(indices)))
  stdk=std(UL(v(indices))./x(v(indices)))
  med=median (UL(v)./x(v))
  med=median(UL(v(indices))./x(v(indices)))
  %figure
  %plot((fmax(v(indices))+fmin(v(indices)))/2, UL(v(indices))./x(v(indices)), '.')
end
 
 plot( (fmax+fmin)/2, hmax, 'c-')
 hold on
plot( (fmax+fmin)/2, UL, 'b-o')
hold on
plot( (fmax+fmin)/2, 11*x, 'r-x')
grid
xlabel('Frequency (Hz)')
ylabel('Hough coarse 95 % UL')
legend('not converged','MC-measured','predicted')
plot( (fmax+fmin)/2, hmax, 'c-')

 
 figure
 plot( (fmax+fmin)/2, hmax, 'c-')
 hold on
plot( (fmax+fmin)/2, UL, 'b-')
hold on
plot( (fmax+fmin)/2, 11*x, 'r-')
grid
xlabel('Frequency (Hz)')
ylabel('Hough coarse 95 % UL')
legend('not converged','MC-measured','predicted')

