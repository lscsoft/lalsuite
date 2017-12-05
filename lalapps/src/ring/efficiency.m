clear
plottrigs

%%%%%%%%%%% Efficiency of triples in triple time %%%%%%%%%%%%%%%
N = 2.303;
rhoBL = 0.0198;

flow=[50,50,100,200,500,1000];
fhigh=[2000,100,200,500,1000,2000];
%flow=[45,45,200,700];
%fhigh=[2500,200,700,2500];
%nbins=[30,20,20,20,20,20];
nbins=[50,30,30,30,30,30];
%nbins=[50,30,30,30];
%flow =45;
%fhigh = 2500;
%nbins=30;
first=10^(-2);
last=10^3;

type = 'tintt';  % tintt, dintt, dintt

%for i=1:length(flow)
for i=1:1

  % make bins for histogramming
  binw(i) = (log10(last)-log10(first))/nbins(i);    % bin width
  binleft = (log10(first)+binw(i)*(0:(nbins(i))))'; % left edge of each bin. the last entry
                                                    %     is the far right edge of range
  bincent = 10.^(binleft+binw(i)/2);                % the center of the bins
  bincent = bincent(1:end-1);                       % get rid of last entry since this isn't
                                                    %     a bin (its just a boundary)

  if strcmp(type,'tintt')
    
    % get indices of injections that were found as triples
    tripind=intersect(injH1t.ind,fsim.ind(fsim.f>flow(i)&fsim.f<fhigh(i)));
    % bin the found injections
    Nfound=histc(log10(fsim.d(tripind)),binleft); 
    Nfound=Nfound(1:end-1);% get rid of the last bin as its not really a bin 
               % (this is necessary to do as histc doesnt count anything beyound the last entry 
               %    in binleft. Thus if the left edge of the last bi was the last entry in bin left
               %    none of these injections would be counted.

    % missed +  found = total number of injections made
    total=[fsim.d(fsim.f>flow(i)&fsim.f<fhigh(i));msim.d(msim.f>flow(i)&msim.f<fhigh(i))];
    Ntotal=histc(log10(total),binleft); % bin the injections
    Ntotal=Ntotal(1:end-1);
  
    % calculate the efficiency and binomial error bars
    eff=Nfound./Ntotal;
    err=sqrt(eff.*(1-eff)./Ntotal);

  elseif strcmp(type,'dintt')

    doubind=[fsimH1L1dt.d(fsimH1L1dt.f>flow(i)&fsimH1L1dt.f<fhigh(i));...
             fsimH1H2dt.d(fsimH1H2dt.f>flow(i)&fsimH1H2dt.f<fhigh(i));...
             fsimL1H2dt.d(fsimL1H2dt.f>flow(i)&fsimL1H2dt.f<fhigh(i));...
             fsimH1L1dd.d(fsimH1L1dd.f>flow(i)&fsimH1L1dd.f<fhigh(i));...
             fsimH1H2dd.d(fsimH1H2dd.f>flow(i)&fsimH1H2dd.f<fhigh(i));...
             fsimL1H2dd.d(fsimL1H2dd.f>flow(i)&fsimL1H2dd.f<fhigh(i))];

    Nfound=histc(log10(doubind),binleft);
    Nfound=Nfound(1:end-1);
    total=[fsim.d(fsim.f>flow(i)&fsim.f<fhigh(i));msim.d(msim.f>flow(i)&msim.f<fhigh(i))];
    Ntotal=histc(log10(total),binleft);
    Ntotal=Ntotal(1:end-1);
 
    eff=Nfound./Ntotal;
    err=sqrt(eff.*(1-eff)./Ntotal);

  elseif strcmp(type,'dindt')
   
    doubH1L1ind=[fsimH1L1dt.d(fsimH1L1dt.f>flow(i)&fsimH1L1dt.f<fhigh(i));...
                 fsimH1L1dd.d(fsimH1L1dd.f>flow(i)&fsimH1L1dd.f<fhigh(i))];

    doubH1H2ind=[fsimH1H2dt.d(fsimH1H2dt.f>flow(i)&fsimH1H2dt.f<fhigh(i));...
                 fsimH1H2dd.d(fsimH1H2dd.f>flow(i)&fsimH1H2dd.f<fhigh(i))];

    doubH2L1ind=[fsimL1H2dt.d(fsimL1H2dt.f>flow(i)&fsimL1H2dt.f<fhigh(i));...
                 fsimL1H2dd.d(fsimL1H2dd.f>flow(i)&fsimL1H2dd.f<fhigh(i))];

    NfoundH1L1=histc(log10(doubH1L1ind),binleft);
    NfoundH1L1=NfoundH1L1(1:end-1);

    NfoundH1H2=histc(log10(doubH1H2ind),binleft);
    NfoundH1H2=NfoundH1H2(1:end-1);

    NfoundH2L1=histc(log10(doubH2L1ind),binleft);
    NfoundH2L1=NfoundH2L1(1:end-1);

    total=[fsim.d(fsim.f>flow(i)&fsim.f<fhigh(i));msim.d(msim.f>flow(i)&msim.f<fhigh(i))];
    Ntotal=histc(log10(total),binleft);
    Ntotal=Ntotal(1:end-1);

    effH1L1=NfoundH1L1./Ntotal;
    errH1L1=sqrt(effH1L1.*(1-effH1L1)./Ntotal);

    effH1H2=NfoundH1H2./Ntotal;
    errH1H2=sqrt(effH1H2.*(1-effH1H2)./Ntotal);

    effH2L1=NfoundH2L1./Ntotal;
    errH2L1=sqrt(effH2L1.*(1-effH2L1)./Ntotal);

  end

  if strcmp(type,'tintt') || strcmp(type,'dintt')
    tintt.dV=4.*pi.*log(10).*binw(i).*bincent.^3 ;
    tintt.V(i)=sum(tintt.dV.*eff)
 
    % include errors 
    sigV=sqrt(sum(err.^2.*tintt.dV.^2));
    cal=tintt.V(i)*(1-1./(1+0.05)^3);
    toterr=sqrt(sigV^2+cal^2);
    err90pc=toterr*1.64;
    %tintt.V(i)=tintt.V(i)-err90pc;   % uncomment to include errors

    tintt.r(i)=(tintt.V(i)/pi*3/4)^(1/3);         % "effective reach"
   % tintt.T(i) = 0.041453;                        % live time
    tintt.T(i) = 0.0375;
    tintt.VT(i) = tintt.V(i).*tintt.T(i);
    tintt.R(i) = N / rhoBL ./ tintt.T(i) ./ tintt.V(i);

  elseif strcmp(type,'dindt')
    dindtH1H2.V(i)=sum(4.*pi.*log(10).*binw(i).*effH1H2.*bincent.^3);
    dindtH1H2.r(i)=(dindtH1H2.V(i)/pi*3/4)^(1/3);              
    dindtH1H2.T(i) = 0.014306;      

    dindtH1L1.V(i)=sum(4.*pi.*log(10).*binw(i).*effH1L1.*bincent.^3);
    dindtH1L1.r(i)=(dindtH1L1.V(i)/pi*3/4)^(1/3);
    dindtH1L1.T(i) = 0.005310;

    dindtH2L1.V(i)=sum(4.*pi.*log(10).*binw(i).*effH2L1.*bincent.^3);
    dindtH2L1.r(i)=(dindtH2L1.V(i)/pi*3/4)^(1/3);
    dindtH2L1.T(i) = 0.004405;

    dindtH1H2.VT(i) = dindtH1H2.V(i).*dindtH1H2.T(i);
    dindtH1L1.VT(i) = dindtH1L1.V(i).*dindtH1L1.T(i);
    dindtH2L1.VT(i) = dindtH2L1.V(i).*dindtH2L1.T(i);

    dindt.VT(i) = dindtH1H2.V(i).*dindtH1H2.T(i) + dindtH1L1.V(i).*dindtH1L1.T(i)...
                     + dindtH2L1.V(i).*dindtH2L1.T(i);

    dindtH1H2.R(i) = N / rhoBL ./ dindtH1H2.T(i) ./ dindtH1H2.V(i);
    dindtH1L1.R(i) = N / rhoBL ./ dindtH1L1.T(i) ./ dindtH1L1.V(i);
    dindtH2L1.R(i) = N / rhoBL ./ dindtH2L1.T(i) ./ dindtH2L1.V(i);

    dindt.R(i) = N / rhoBL ./ dindt.VT(i);
  end


  %plot efficiency curve
  if strcmp(type,'tintt') || strcmp(type,'dintt')
    figure
    h_p1=semilogx(bincent,eff,'bo');
    set(h_p1,'MarkerSize',7);
    hold on
    h_p2=semilogx(bincent,eff);
    set(h_p2,'LineWidth',2)
    for n=1:length(eff)    % plot the error bars
      h_pe=semilogx([bincent(n) bincent(n)],[eff(n)-err(n) eff(n)+err(n)],'r.-');
      set(h_pe,'MarkerSize',15,'LineWidth',2)
    end
    hold off
    h_xlab=xlabel('d / Mpc');
    set(h_xlab,'FontSize',16,'FontName','Times');
    h_ylab=ylabel('\epsilon');
    set(h_ylab,'FontSize',16,'FontName','Times');
    eval(['title(''Efficiency vs distance for triples in triple time in frequency band ' num2str(flow(i)) '-' num2str(fhigh(i)) 'Hz'')']);
    grid on
    set(gca,'FontSize',16,'FontName','Times');
axis([1e-3 1e3 0 1])
    eval(['saveas(gcf,''effvd_' type '_' num2str(flow(i)) '_' num2str(fhigh(i)) '.png'')'])
  end

%  clear Nfound Ntotal tripind  eff binw bincent binleft err total
end


