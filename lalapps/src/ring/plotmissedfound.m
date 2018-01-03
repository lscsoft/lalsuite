plottrigs

% Hanford effectve distance vs timefigure
ht=loglog(fsim.f(injH1t.ind),fsim.dh(injH1t.ind),'x');
hold on
hH1H2=loglog([fsimH1H2dt.f;fsimH1H2dd.f],[fsimH1H2dt.dh;fsimH1H2dd.dh],'g*');
hH1L1=loglog([fsimH1L1dt.f;fsimH1L1dd.f],[fsimH1L1dt.dh;fsimH1L1dd.dh],'c*');
hL1H2=loglog([fsimL1H2dt.f;fsimL1H2dd.f],[fsimL1H2dt.dh;fsimL1H2dd.dh],'m*');
hm=loglog([missed.f;missedd.f],[missed.dh;missedd.dh],'ro');
hl=plot([50,50],[1e-2,1e5],'k');
hu=plot([2000,2000],[1e-2,1e5],'k');
grid on
grid minor
h_xlab=xlabel('f / Hz');
h_ylab=ylabel('d_H / Mpc');
hleg=legend('triples','H1H2 doubles','H1L1 doubles','L1H2 doubles','missed');
set(h_xlab,'FontSize',16,'FontName','Times');
set(h_ylab,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16,'FontName','Times');
saveas(gcf,'mf_allvH.png')

