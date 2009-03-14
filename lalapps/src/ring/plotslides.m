% Plots of H1L1 background

bgH1inL1d=load('bgH1inL1_doub.mat');
bgL1inH1d=load('bgL1inH1_doub.mat');

for i=1:length(bgH1inL1d.id)
  tmp=bgH1inL1d.id(i);
  alsotmp=tmp{1};
  bgH1inL1d.slidenum(i)=str2num(alsotmp(33:36));
end

for i=1:length(bgH1inL1d.slidenum)
  if bgH1inL1d.slidenum(i)>5000 % negative time-slide
     H1L1event(i)=(bgH1inL1d.slidenum(i)-5000)*-1;
  else                              % positive timeslide
     H1L1event(i)=bgH1inL1d.slidenum(i);
  end
end

figure
hist(H1L1event,101)
hold on
h=findobj(gca,'Type','patch');
set(h,'FaceColor','r')  %set the time-slides to red
[A,B]=hist(H1L1event,101);
ave=sum(A)/100;
plot([-60,60],[ave,ave],'k--')
h_xlab=xlabel('SlideNumber');
h_ylab=ylabel('Number of Coincidences');
h_t=title('Number of background events in H1L1');
set(h_xlab,'FontSize',16);
set(h_ylab,'FontSize',16);
set(h_t,'FontSize',16);
grid on
axis([-60 60 0 320])
set(gca,'FontSize',16);
axis autoy          
saveas(gcf,'H1L1bghist.png')
ave_num_bg = ave


figure
hist(log10(bgH1inL1d.f),100)
h_xlab=xlabel('log_{10}(f)');
h_ylab=ylabel('Number of Coincidences');
h_t=title('Number of background events as a function of frequency');
set(gca,'FontSize',16);
grid on
set(h_xlab,'FontSize',16);
set(h_ylab,'FontSize',16);
set(h_t,'FontSize',16);
saveas(gcf,'H1L1bgfhist.png')
