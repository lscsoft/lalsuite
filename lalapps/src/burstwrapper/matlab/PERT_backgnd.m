path = '/home/jsylvest/public_html/TFClusters/PERT/';

bpp = {'0.0210213','0.0215403','0.0220594','0.0225784','0.0230975','0.0236165','0.0241356','0.0246546','0.0251737','0.0256927'};

clf;

col = 'rgbkymc';
sym = '.O+*';

j = 1;

for k=1:length(bpp)

 [H1,H1Segs] = read_burstdso_files([path 'bpp_' bpp{k} '/jobH1*.0']);

 plot(H1Segs.time_s, H1Segs.Nevents,[col(1+mod(k,length(col))) sym(j)]);

 if(1+mod(k,length(col)) == length(col))
     j=j+1;
 end

 hold on

 drawnow

end

 xlabel('GPS time');
 ylabel('number of events per 300s segment');

 legend(bpp);
 
