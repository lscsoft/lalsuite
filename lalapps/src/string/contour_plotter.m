clear
load gammaInitialLIGO.dat

gamma=gammaInitialLIGO;

T=365*24*3600;
%T=1363412.0;

numberofevents=-log(1-0.90);
%numberofevents=1.0000001

epsi=gamma(:,3);
Gmu=gamma(:,4);

Gmu2 = unique(gamma(:,4))';
eps2 = unique(gamma(:,3));

rate = zeros(length(eps2),length(Gmu2));

p=1;
for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = gamma(:,4) == Gmu2(kk) & gamma(:,3) == eps2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5)/p;
    else
      disp('PROBLEM')
    end
  end
end
figure(1);loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=1','FontSize',14)
contour(Gmu2,eps2,rate*T, [numberofevents],'r*');

p=1e-1;
for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = gamma(:,4) == Gmu2(kk) & gamma(:,3) == eps2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5)/p;
    else
      disp('PROBLEM')
    end
  end
end
figure(2);loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-1}','FontSize',14)
contour(Gmu2,eps2,rate*T, [numberofevents],'r*');

p=1e-2;
for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = gamma(:,4) == Gmu2(kk) & gamma(:,3) == eps2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5)/p;
    else
      disp('PROBLEM')
    end
  end
end
figure(3);loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-2}','FontSize',14)
contour(Gmu2,eps2,rate*T, [numberofevents],'r*');

p=1e-3;
for kk = 1:length(Gmu2)
  for jj = 1:length(eps2)
    cut = gamma(:,4) == Gmu2(kk) & gamma(:,3) == eps2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5)/p;
    else
      disp('PROBLEM')
    end
  end
end
figure(4);loglog(10^-12,1); hold on; loglog(10^-12,10^-12); loglog(10^-6,10^-12);loglog(10^-6,1);grid on;
xlabel('G\mu','FontSize',14)
ylabel('\epsilon','FontSize',14)
title('p=10^{-3}','FontSize',14)
contour(Gmu2,eps2,rate*T, [numberofevents],'r*');
