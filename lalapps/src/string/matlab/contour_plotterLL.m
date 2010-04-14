clear
load /home/siemens/LIGOCS-S4paper/gammaLL_S4UL.dat

gamma=gammaLL_S4UL;

%T=365*24*3600;
T=1363412.0;

%numberofevents=-log(1-0.90);
numberofevents=1.0000001

p=gamma(:,1);
Gmu=gamma(:,2);

Gmu2 = unique(gamma(:,2))';
p2 = unique(gamma(:,1));

rate = zeros(length(p2),length(Gmu2));


for kk = 1:length(Gmu2)
  for jj = 1:length(p2)
    cut = gamma(:,2) == Gmu2(kk) & gamma(:,1) == p2(jj);
    if sum(cut) == 1
      rate(jj,kk) = gamma(cut,5);
    else
      disp('PROBLEM')
    end
  end
end
figure(2);loglog(10^-12,10^-4); hold on; loglog(10^-12,1); loglog(10^-4,1);loglog(10^-4,10^-4);grid on;
xlabel('G\mu','FontSize',14)
ylabel('p','FontSize',14)
title('Large Loops','FontSize',14)
contour(Gmu2,p2,rate*T, [numberofevents],'r*');

