clear
dlna=0.01;
a=-1:dlna:15;

b=zeros(size(a,2),1);
[i,j]=find(a>=0.0);
b(j)=1.0;
c=zeros(size(a,2),1);

op(:,1)=a;
op(:,2)=b;
op(:,3)=c;

fid=fopen('efficiency2.txt','w');
save('efficiency2.txt','op','-ASCII','-TABS','-APPEND')
fclose(fid);