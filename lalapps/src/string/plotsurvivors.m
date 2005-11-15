clear
path(path,'/home/siemens/ligotools/matlab')

cd /scratch2/xavi/CosmicStrings0/slides/

H1ctrig_foreground=readMeta('H1-BURCA_H1L1_P_0_--1-0.xml','sngl_burst');

N=size(H1ctrig_foreground.amplitude,1);

cd /scratch2/xavi/CosmicStrings0/plotdata/

t=0:1/16384:15-1/16384;

N=1;

[i,j]=find(t>)

for k=1:N
   
    data_fileH1=['H1data',num2str(k),'.txt'];
    data_fileH2=['H2data',num2str(k),'.txt'];
    data_fileL1=['L1data',num2str(k),'.txt'];

    H1data=load (data_fileH1);
    H2data=load (data_fileH2);
    L1data=load (data_fileL1);

end

