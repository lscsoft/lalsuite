% This is the outline of the over-arching ringdown plotting script.
% It hasn't been thoroughly test yet, but I wanted to get it into 
% the repository anyway.
% The idea is that the user will run this when the DAG has completed
% and will not need to copy or move any fields around. 
%
% usage example plot_ring(871147814,871234214,['CAT_2';'CAT_3'])
%
% author: Lisa Goggin 11/19/09


function plots_ring(tstart,tend,veto_type)

% Read in the injection details from a file having run the python script
[inj_type,dirs,seed]=textread('messedinjinfo.txt','%s %s %d');
INJ_type=upper(inj_type);

% group into the different types of injections
[types,unknown,type_ind]=unique(INJ_type);
for i=1:length(types)
  eval(['dir_type_' num2str(i) '=dirs(type_ind==i);'])
  eval(['seed_type_' num2str(i) '=seed(type_ind==i);'])
end
 
duration=tend-tstart;

% make a new directory and go there
if ~exist('plots')
  mkdir plots
end
cd plots

% Copy over the coinc and injection files from the injection directories  
% (just triples for the moment)
for j=1:length(inj_type)
  eval(['copyfile ../' dirs{j} '/H1H2L1*COINCRINGREAD*SUMMARY*.xml .'])
  eval(['copyfile ../' dirs{j} '/HL-INJECTIONS*.xml .'])
end
 
% copy over the playground, timeslide files
copyfile ../playground/H1H2L1*COINCRINGREAD_SUMMARY*.xml .
copyfile ../full_data/H1H2L1*COINCRINGREAD_SLIDE_SUMMARY*.xml .

% copy over veto files
copyfile ../segments/*CATEGORY*.txt .

ifo=['H1';'H2';'L1'];

% loop over each of the veto types
for i=1:length(veto_type)
  fprintf('veto_type=%s\n',veto_type(i,:))
  % read in timeslides
 eval(['separate(''bg'','''',veto_type(i,:),''H1H2L1-COINCRINGREAD_SLIDE_SUMMARY_FIRST_H1H2L1_FULL_DATA_' veto_type(i,:) '_VETO-' num2str(tstart) '-' num2str(duration) '.xml'')'])
 eval(['separate(''bg'','''',veto_type(i,:),''H1H2L1-COINCRINGREAD_SLIDE_SUMMARY_FIRST_H1H2L1_FULL_DATA_' veto_type(i,:) '_VETO-' num2str(tstart) '-' num2str(duration) '.xml'')'])
  
  % make plot of timeslide parameters  
  for k=1:length(ifo)-1
    eval(['background_param_plots(veto_type(i,:),''triple'', ifo(k,:), ifo(k+1,:),''' veto_type(i,:) '_bg_' ifo(k,:) 'trip.mat'','''  veto_type(i,:) '_bg_' ifo(k+1,:) 'trip.mat'' )'])
  end

  % read in the playground
  eval(['separate(''pg'','''',veto_type(i,:),''H1H2L1-COINCRINGREAD_SUMMARY_FIRST_H1H2L1_PLAYGROUND_' veto_type(i,:) '_VETO-' num2str(tstart) '-' num2str(duration) '.xml'')'])
  
  % plot playground/background snr/snr scatter plots
  for k=1:length(ifo)-1
    eval(['scatter_plots( veto_type(i,:) ,''triple'',ifo(k,:), ifo(k+1,:),''' veto_type(i,:) '_bg_' ifo(k,:) 'trip.mat'',''' veto_type(i,:) '_bg_' ifo(k+1,:) 'trip.mat'','''',''pg'',''' veto_type(i,:) '_pg_' ifo(k,:) 'trip.mat'',''' veto_type(i,:) '_pg_' ifo(k+1,:) 'trip.mat'')'])
   % eval(['scatter_plots( veto_type(i,:) ,''double'',ifo(k,:), ifo(k+1,:),''' veto_type(i,:) '_bg_' ifo(k,:) 'in' ifo(k+1,:)  'doub.mat'',''' veto_type(i,:) '_bg_' ifo(k+1,:) 'in' ifo(k,:) 'doub.mat'',types{j},''inj'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k,:) 'in' ifo(k+1,:) 'doub.mat'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k+1,:) 'in' ifo(k,:) 'doub.mat'')'])  
  end
 %{
  % plot playground/background histogram 
  eval(['trip_hist(veto_type(i,:),''' veto_type(i,:) '_bg_H1trip.mat'',''' veto_type(i,:) '_bg_H2trip.mat'',''' veto_type(i,:) '_bg_L1trip.mat'', ''pg'', ''' veto_type(i,:) '_pg_H1trip.mat'',''' veto_type(i,:) '_pg_H2trip.mat'',''' veto_type(i,:) '_pg_L1trip.mat'')'])
%}

  % loop over injection types
  for j=1:length(types)
    fprintf('veto_type=%s\n',types{j})
    % make plots of injected quantities
    eval(['InjFileDir=dir(''HL-INJ*' types{j} '*.xml'');'])
    inj_list={InjFileDir.name};
    injection_param_plots(types{j},inj_list)

    % separate the injections into triples and doubles
    eval(['FoundFileDir=dir(''H1H2L1*COINC*' types{j} '*FOUND*SUMM*' veto_type(i,:) '*.xml'');'])
    file_list={FoundFileDir.name};
    separate('inj',types{j},veto_type(i,:),file_list)

    % make the missed/found plots and write out the sim*.mat files
    eval(['MissedFileDir=dir(''H1H2L1*COINC*' types{j} '*MISSED*SUMM*' veto_type(i,:) '*.xml'');'])
    missed_list={MissedFileDir.name};
    if strcmp(veto_type(i,:),'CAT_2')
      eval(['veto_list=[''H1-CATEGORY_2_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt'';''H1-CATEGORY_2_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt'';''L1-CATEGORY_2_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt''];'])
    elseif strcmp(veto_type(i,:),'CAT_3')
      eval(['veto_list=[''H1-CATEGORY_3_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt'';''H1-CATEGORY_3_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt'';''L1-CATEGORY_3_VETO_SEGS-' num2str(tstart) '-' num2str(duration) '.txt''];'])
    end
    plottrigs(types{j},veto_type(i,:),file_list, missed_list,veto_list)

    % make plots of parameter accuracy
    ifo=['H1';'H2';'L1'];
    for k=1:length(ifo)
      eval(['injected_list = [''' veto_type(i,:) '_' types{j} '_simH1H2L1.mat''];'])
      eval(['detected_list = [''' veto_type(i,:) '_inj_' types{j} '_' ifo(k,:) 'trip.mat''];'])
      param_accuracy(veto_type(i,:), types{j}, ifo(k,:), injected_list, detected_list )    
    end
  
    % snr scatter plots
    for k=1:length(ifo)-1
      eval(['scatter_plots( veto_type(i,:) ,''triple'',ifo(k,:), ifo(k+1,:),''' veto_type(i,:) '_bg_' ifo(k,:) 'trip.mat'',''' veto_type(i,:) '_bg_' ifo(k+1,:) 'trip.mat'',types{j},''inj'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k,:) 'trip.mat'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k+1,:) 'trip.mat'')'])
   % eval(['scatter_plots( veto_type(i,:) ,''double'',ifo(k,:), ifo(k+1,:),''' veto_type(i,:) '_bg_' ifo(k,:) 'in' ifo(k+1,:)  'doub.mat'',''' veto_type(i,:) '_bg_' ifo(k+1,:) 'in' ifo(k,:) 'doub.mat'',types{j},''inj'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k,:) 'in' ifo(k+1,:) 'doub.mat'',''' veto_type(i,:) '_inj_' types{j} '_' ifo(k+1,:) 'in' ifo(k,:) 'doub.mat'')'])  
    end
  end
end
