function rplot(s)

clf

subplot(3,1,1)
     plot(s.H1Segs.time_s - min(s.H1Segs.time_s),s.H1Segs.Nevents,'.');
     ylabel('# H1 events');
     xlabel(['time - ' num2str(min(s.H1Segs.time_s)) ' s']);

subplot(3,1,2)
     plot(s.H2Segs.time_s - min(s.H2Segs.time_s),s.H2Segs.Nevents,'r.');
     ylabel('# H2 events');
     xlabel(['time - ' num2str(min(s.H2Segs.time_s)) ' s']);

subplot(3,1,3)
     plot(s.L1Segs.time_s - min(s.L1Segs.time_s),s.L1Segs.Nevents,'g.');
     ylabel('# L1 events');
     xlabel(['time - ' num2str(min(s.L1Segs.time_s)) ' s']);

