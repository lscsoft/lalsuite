
% $Id$
%
% script to produce the figures

load AllSkyMaxResults

figure
set(gcf,'PaperPosition', [0.25 0.5 8 10.7])

subplot(5,1,1);

  plot(FreqValues, MaxValues, 'b.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid;
  axis([200 210 400 800]);

subplot(5,1,2);

  plot(FreqValues, MaxValues, 'b.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid;
  axis([210 220 400 800]);


subplot(5,1,3);

  plot(FreqValues, MaxValues, 'b.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid;
  axis([220 230 400 800]);


subplot(5,1,4);

  plot(FreqValues, MaxValues, 'b.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid;
  axis([230 240 400 800]);

subplot(5,1,5);

  plot(FreqValues, MaxValues, 'b.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid;
  axis([240 250 400 800]);

print(gcf, '-djpeg', ['H1maxAll_200_250.jpeg'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(5,1,1);
  axis([250 260 400 800]);

subplot(5,1,2);
  axis([260 270 400 800]);

subplot(5,1,3);
  axis([270 280 400 800]);

subplot(5,1,4);
  axis([280 290 400 800]);

subplot(5,1,5);
  axis([290 300 400 800]);
  
print(gcf, '-djpeg', ['H1maxAll_250_300.jpeg'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,1);
  axis([300 310 400 800]);

subplot(5,1,2);
  axis([310 320 400 800]);

subplot(5,1,3);
  axis([320 330 400 800]);

subplot(5,1,4);
  axis([330 340 400 800]);

subplot(5,1,5);
  axis([340 350 400 800]);
    
print(gcf, '-djpeg', ['H1maxAll_300_350.jpeg'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,1);
  axis([350 360 400 800]);

subplot(5,1,2);
  axis([360 370 400 800]);

subplot(5,1,3);
  axis([370 380 400 800]);

subplot(5,1,4);
  axis([380 390 400 800]);

subplot(5,1,5);
  axis([390 400 400 800]);
    
print(gcf, '-djpeg', ['H1maxAll_350_400.jpeg'])

%print(gcf, '-djpeg', ['trend_01123_' channel '.jpeg'])
%print(gcf, '-dpsc', ['trend_' channel '.ps'])
%print(gcf, '-depsc', ['trend_' channel '.eps'])
%print(gcf, '-dtiff', ['trend_' channel '.tiff'])
