function bplot(dt,N,olap)

     clf

     subplot(1,2,1);

     plot(dt,N,'.-');
     
     xlabel('dt [s]');
     ylabel('number of events');

     yl = get(gca,'YLim');

     set(gca,'XLim',[min(dt) max(dt)]);

     xt = get(gca,'XTick');
     set(gca,'XTick',xt(1:end-1));



     g1 = gca;

     subplot(1,2,2);

     k = floor(min(yl)):1:ceil(max(yl));
     [n,y]=hist(N,k);

     barh(y,n,1);
     set(gca,'YLim',yl);

     xlabel('number of lags');

     g2 = gca;

     gp = get(g1,'Position');
     set(g1,'Position',[gp(1) gp(2) 0.55 gp(4)]);
     p2 = gp(1) + 0.55;

     gp = get(g2,'Position');
     set(g2,'Position',[p2 gp(2) 0.2 gp(4)]);

     set(gca,'YTickLabel',[]);


     hold on

     mr = mean(N);
     kk = min(k):max(k);
     pk = exp(log(mr)*kk - mr - log_factorial(kk));


if(k(2) > k(1)+1)
 pk = buffer(pk,k(2)-k(1));
 pk = sum(pk);
end

NSeg = length(N);
plot(NSeg*pk,k,'r');

if(nargin == 3)
 disp(['Mean rate: ' num2str(mean(N./olap)) ' Hz']);
 disp(['Std rate: ' num2str(std(N./olap)) ' Hz']);
end

 I = find(n>0);
disp(['chisq: ' num2str(sum((NSeg*pk(I)-n(I)).^2./n(I)))]);
disp(['DOF: ' num2str(length(I)-1)]);

xl = get(g2,'XLim');
yl = get(g2,'YLim');
text(xl(1) + 0.33*(xl(2)-xl(1)), yl(1)+0.95*(yl(2)-yl(1)),['chisq: ' num2str(sum((NSeg*pk(I)-n(I)).^2./n(I)))]);
text(xl(1) + 0.33*(xl(2)-xl(1)), yl(1)+0.9*(yl(2)-yl(1)),['DOF: ' num2str(length(I)-1)]);

text(xl(1) + 0.33*(xl(2)-xl(1)), yl(1)+0.85*(yl(2)-yl(1)),['mean: ' num2str(mean(N))]);
text(xl(1) + 0.33*(xl(2)-xl(1)), yl(1)+0.8*(yl(2)-yl(1)),['var: ' num2str(var(N))]);
