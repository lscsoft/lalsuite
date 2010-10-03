clear

%files with rate data
%load /home/siemens/lscdocs/bursts/projects/cosmic_strings/s4cosmic_strings/gamma_S4sens.dat
load /archive/home/amergl/s5cosmicstring/gamma_S4UL.dat
%load /home/siemens/lscdocs/bursts/projects/cosmic_strings/s4cosmic_strings/gamma_LIGOsens.dat

%one year
%T=365*24*3600;

%S4 live-time 
%T=1363502.0;
T=342610.0;

%gammaS=gamma_S4sens;
gammaUL=gamma_S4UL;
%gammaL=gamma_LIGOsens;

%if size(gammaS) == size(gammaUL)
    %disp('Sizes of UL and sensitivity rate files agree. Good.')
%else
    %disp('Sizes of UL and sensitivity rate files do not agree. Exiting...')
    %exit
%end

numberofeventsUL=-log(1-0.90);
numberofeventsS=1.0000001;

epsi=gammaUL(:,3);
Gmu=gammaUL(:,4);

Gmu2 = unique(gammaUL(:,4))';
eps2 = unique(gammaUL(:,3));

rate = zeros(length(eps2),length(Gmu2));

for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = (gammaUL(:,4) == Gmu2(kk)) & (gammaUL(:,3) == eps2(jj));
    if sum(cut) == 1

      %rateS(jj,kk) = gammaS(cut,5);
      %rateSMin(jj,kk) = gammaS(cut,6);
      %rateSMax(jj,kk) = gammaS(cut,7);
      
      rateUL(jj,kk) = gammaUL(cut,5);
      rateULMin(jj,kk) = gammaUL(cut,6);
      rateULMax(jj,kk) = gammaUL(cut,7);
      
      %rateSLIGO(jj,kk)=gammaL(cut,5);
      
    else
      disp('PROBLEM')
    end
  end
end

%p=1 case
figure(1);
loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=1','FontSize',14)
contour(Gmu2,eps2,rateUL*T, [numberofeventsUL],'r*');
contour(Gmu2,eps2,rateULMin*T, [numberofeventsUL],'r:');
contour(Gmu2,eps2,rateULMax*T, [numberofeventsUL],'r:');
%contour(Gmu2,eps2,rateS*T, [numberofeventsS],'k*');
%contour(Gmu2,eps2,rateSMin*T, [numberofeventsS],'k:');
%contour(Gmu2,eps2,rateSMax*T, [numberofeventsS],'k:');

%contour(Gmu2,eps2,rateSLIGO*(365*24*3600), [numberofeventsS],'b*');

%p=1e-1 case
figure(2);
loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-1}','FontSize',14)
contour(Gmu2,eps2,rateUL*T/1e-1, [numberofeventsUL],'r*');
contour(Gmu2,eps2,rateULMin*T/1e-1, [numberofeventsUL],'r:');
contour(Gmu2,eps2,rateULMax*T/1e-1, [numberofeventsUL],'r:');
%contour(Gmu2,eps2,rateS*T/1e-1, [numberofeventsS],'k*');
%contour(Gmu2,eps2,rateSMin*T/1e-1, [numberofeventsS],'k:');
%contour(Gmu2,eps2,rateSMax*T/1e-1, [numberofeventsS],'k:');

%contour(Gmu2,eps2,rateSLIGO*(365*24*3600)/1e-1, [numberofeventsS],'b*');

%p=1e-2 case
figure(3);
loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-2}','FontSize',14)
contour(Gmu2,eps2,rateUL*T/1e-2, [numberofeventsUL],'r*');
contour(Gmu2,eps2,rateULMin*T/1e-2, [numberofeventsUL],'r:');
contour(Gmu2,eps2,rateULMax*T/1e-2, [numberofeventsUL],'r:');
%contour(Gmu2,eps2,rateS*T/1e-2, [numberofeventsS],'k*');
%contour(Gmu2,eps2,rateSMin*T/1e-2, [numberofeventsS],'k:');
%contour(Gmu2,eps2,rateSMax*T/1e-2, [numberofeventsS],'k:');

%contour(Gmu2,eps2,rateSLIGO*(365*24*3600)/1e-2, [numberofeventsS],'b*');

%p=1e-3 case
figure(4);
loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-3}','FontSize',14)
contour(Gmu2,eps2,rateUL*T/1e-3, [numberofeventsUL],'r*');
contour(Gmu2,eps2,rateULMin*T/1e-3, [numberofeventsUL],'r:');
contour(Gmu2,eps2,rateULMax*T/1e-3, [numberofeventsUL],'r:');
%contour(Gmu2,eps2,rateS*T/1e-3, [numberofeventsS],'k*');
%contour(Gmu2,eps2,rateSMin*T/1e-3, [numberofeventsS],'k:');
%contour(Gmu2,eps2,rateSMax*T/1e-3, [numberofeventsS],'k:');

%contour(Gmu2,eps2,rateSLIGO*(365*24*3600)/1e-3, [numberofeventsS],'b*');
