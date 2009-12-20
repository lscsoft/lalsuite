function injection_param_plots( injtype, file_list )

%
% NULL = injection_param_plots ( injtype, file_list )
%
% example: file_list = ['HL-INJECTIONS_1000_EOBNR_A-871147814-86400.xml';'HL-INJECTIONS_1002_EOBNR_A-871147814-86400.xml'];
%
% injtype = 'RINGDOWN', 'EOBNR', 'PHENOM'
%
%
% Sarah Caudill, Oct 28th 2009

  N_files=length(file_list);

% read in the injection files

  % create the ringdown parameters structure by reading in the first file
  eval(['rd_params=readMeta(''' file_list{1} ''',''sim_ringdown'',0,''quality,frequency,mass,spin'');'])

  % read in the rest of the injection files
  for i=2:N_files
    eval(['rd_paramsi=readMeta(''' file_list{i} ''',''sim_ringdown'',0,''quality,frequency,mass,spin'');'])
    rd_params.quality=[rd_params.quality;rd_paramsi.quality];
    rd_params.frequency=[rd_params.frequency;rd_paramsi.frequency];
    rd_params.mass=[rd_params.mass;rd_paramsi.mass];
    rd_params.spin=[rd_params.spin;rd_paramsi.spin];
  end

  %============== Injection Parameter Plots ==============%
 
  %%%%%%%%%%%%%%%%%% QUALITY VS FREQUENCY %%%%%%%%%%%%%%%%

  figure
  semilogx(rd_params.frequency,rd_params.quality,'k.')
  hold on
  grid on
  x_lab=xlabel('Freq_{inj} (Hz)');
  y_lab=ylabel('Q_{inj}');
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title(''' injtype ': Quality_{inj} versus Freq_{inj}'');'])
  set(plot_title,'FontSize',16,'FontWeight','b');
  eval(['saveas(gcf,''' injtype '_injectparams_qvsfreq.png'')'])


  %%%%%%%%%%%%%%%%%% SPIN VS MASS %%%%%%%%%%%%%%%%

  figure
  semilogx(rd_params.mass,rd_params.spin,'k.')
  hold on
  grid on
  x_lab=xlabel('Mass_{tot} (M_{sun})');
  y_lab=ylabel('a_{inj}');
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title(''' injtype ': Spin_{inj} versus Total Mass_{inj}'');'])
  set(plot_title,'FontSize',16,'FontWeight','b');
  eval(['saveas(gcf,''' injtype '_injectparams_spinvsmass.png'')'])


  if strcmp(injtype,'EOBNR')||strcmp(injtype,'PHENOM')

    % create the ringdown params structure and inspiral params structure by reading in the first file
    eval(['insp_params=readMeta(''' file_list{1} ''',''sim_inspiral'',0,''mass1,mass2'');'])

    % read in the rest of the injection files
    for i=2:N_files
      eval(['insp_paramsi=readMeta(''' file_list{i} ''',''sim_inspiral'',0,''mass1,mass2'');'])
      insp_params.mass1=[insp_params.mass1;insp_paramsi.mass1];
      insp_params.mass2=[insp_params.mass2;insp_paramsi.mass2];
    end

    %============== Injection Parameter Plots ==============%

    %%%%%%%%%%%%%%%%%% MASS1 VS MASS2 %%%%%%%%%%%%%%%%

    figure
    loglog(insp_params.mass2,insp_params.mass1,'k.')
    hold on
    grid on
    x_lab=xlabel('Mass2_{inj} (M_{sun})');
    y_lab=ylabel('Mass1_{inj} (M_{sun})');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' injtype ': Mass1_{inj} versus Mass2_{inj}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' injtype '_injectparams_mass1vsmass2.png'')'])

  end
