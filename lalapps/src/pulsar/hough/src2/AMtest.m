kk=load('AMtest');

a=exp(-kk(:,8) );
N=kk(:,7);
r1=( kk(:,3) -N.*a)./sqrt(N.*a.*(1-a));
r2=(kk(:,4)-kk(:,9) )./kk(:,10);
plot(r1)
hold on
plot(r2,'r')
