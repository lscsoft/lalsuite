clear
path(path,'/home/siemens/ligotools/matlab')

cd /scratch2/xavi/CosmicStrings0/slides/

H1ctrig_foreground=readMeta('H1-BURCA_H1L1_P_0_--1-0.xml','sngl_burst');

N=size(H1ctrig_foreground.amplitude,1);

cd /scratch2/xavi/CosmicStrings0/plotdata/

t=0:1/16384:15-1/16384;


for k=1:N

    trigger_time=7+H1ctrig_foreground.peak_time_ns(k)*1e-9;
    [i1,j1]=find( t > (trigger_time-0.05) );
    [i2,j2]=find( t(j1) < (trigger_time+0.05) );
    
    data_fileH1=['H1data',num2str(k),'.txt'];
    data_fileH2=['H2data',num2str(k),'.txt'];
    data_fileL1=['L1data',num2str(k),'.txt'];

    H1data=load (data_fileH1);
    H2data=load (data_fileH2);
    L1data=load (data_fileL1);

    figure
    subplot (3,1,1)
    plot(t(j1(j2)),H1data(j1(j2)),'b')
    subplot (3,1,2)
    plot(t(j1(j2)),H2data(j1(j2)),'g')
    subplot (3,1,3)
    plot(t(j1(j2)),L1data(j1(j2)),'r')
end

