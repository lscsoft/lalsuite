chn = 'L1:LSC-AS_Q';
%bacfile = '/home/jsylvest/public_html/burstdso/Tuning/L1.130-400.S2_HPF_6_100.NoClust_gpar_first';
bacfile = '/home/jsylvest/public_html/burstdso/Tuning/L1.band2';

     Mat = {};

     for k=1:size(gpar,1)

     Mat{k} = [ chn ',' ];

     for j=1:size(gpar,2)-1
      Mat{k} = [Mat{k} num2str(gpar(k,j)) ','];
     end;

     Mat{k} = [Mat{k} num2str(gpar(k,end))];
     Mat{k} = TFClustersSTRFormat(Mat{k});

     end


ilwd_write_char([bacfile '.ilwd'],bacfile,bacfile,Mat);
