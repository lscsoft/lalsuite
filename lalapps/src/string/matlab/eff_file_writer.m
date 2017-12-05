clear

a=-46.5023883108:1.841720510000044e-02:-39.1539234785;

b=zeros(size(a,2),1);
[i,j]=find(exp(a)>=1e-20);
b(j)=1.0;
c=zeros(size(a,2),1);

op(:,1)=a;
op(:,2)=b;
op(:,3)=c;

fid=fopen('efficiency2.txt','w');
save('efficiency2.txt','op','-ASCII','-TABS','-APPEND')
fclose(fid);