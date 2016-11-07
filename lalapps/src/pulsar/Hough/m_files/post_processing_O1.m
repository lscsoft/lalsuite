function post_processing_O1(X_input,Y_input,grup_Band,BandNo,BandNf,Ntop,PixelFactor,p,sig_veto,rp_sig,r_W,r_W_cl,CHI_V,chisquare_STD_veto,population_cluster_veto,significance_veto_cluster,x_tref,x_tend,y_tref,y_tend,Extract_S,min_band)

format long

%% SETUP

% % Imput files from Toplist from 2 different data sets
% 
% X_input='/Users/wave/Collected/L1_Collected_%g_%g.dat';
% Y_input='/Users/wave/Collected/H1_Collected_%g_%g.dat';
% 
% grup_min=350; % initial group to be analysed and saved
% grup_Band=50;
% grup_max=350; % last group
% 
% BandNo=50;   %%% first band to analyse
% BandNf=66; %%% last band to analyse
% 
% Ntop=999; % Number of candidate per toplist
% 
% PixelFactor=2;
% p = 16; % degeneration
% 
% sig_veto=0; % threshold on significence
% 
% rp_sig=0.2;        % Pecentual distance in significance between candidates
% r_W=sqrt(14);    % Coincidence radious
% r_W_cl=sqrt(14); % Cluster radious
% 
% CHI_V=[[0,0.4902,1.414,0.3581,1.484];[100,0.2168,1.428,0.1902,1.499];[200,0.1187,1.470,0.06784,1.697]];
% %SETUP CHI2 VETO vector [Fin,A1,A2,B1,B2,Fend]
% chisquare_STD_veto = 5; % Number of STD we allow a candidate before vetoing
% 
% population_cluster_veto=0; % Minimum population of a cluster
% 
% significance_veto_cluster=0; % Minimum significance of a cluster
%
% min_band = 0.1 only one candidate per band
%
% Extract_S=[5,6,7]; select candidate per band per
% 5 = cluster mean significance of the mean 
% 6 = cluster mean significance of the hagmonic 
% 7 = cluster integrated significance of the hagmonic 
%
% %Detector == 'L':
% x_tref = 1126034147;
% x_tend = 1137248624;
% 
% %Detector == 'H':
% y_tref = 1125972653;
% y_tend = 1137250279;

 Tcoh=1800;
 
 Tobs=[x_tend,y_tend]-[x_tref,y_tref];
 
 df0=(1/Tcoh);
 df1=(df0./Tobs);
 
 f0_band=r_W*df0; %%% size of the Band to extract and analyse
 f0_rang_grup_band=10; %%% size of the Band to extract before analysis
 
 sigma_veto_toplist= @(x) sig_veto+0.*x;

    
    f0_grup=grup; % Starting frequency of the grup
    f0_grup_max=grup+grup_Band; % End frequency of the grup
    
    x_Toplist=[];
    y_Toplist=[];
    
    xmissing=[];
    ymissing=[];
    
    follow_up=[];
    Cluster=[];
    I_partial_cluster=[];
    
    %% SETUP CHI2 VETO
    [I,~]=find(CHI_V(:,1) < f0_grup);
    I=max(I);
    A1=CHI_V(I,2);
    A2=CHI_V(I,3);
    B1=CHI_V(I,4);
    B2=CHI_V(I,5);
    
    chi2_STD = @(x,y,p) (y-(p-1)-A1*x.^A2)./(sqrt(2*p-2)+B1*x.^B2);
    
    %% READING TOPLISTS
    for BandN=BandNo:BandNf;

        y_filename=sprintf(Y_input,grup,BandN);
        x_filename=sprintf(X_input,grup,BandN);
        
        x = [];
        y = [];
        
        x=load(x_filename);
        if ~isempty(x)
            x = [x,BandN*ones(length(x(:,1)),1)];
            [~,I]=sort(x(:,5),'descend'); x=x(I(1:Ntop),:); %Reduce toplist
            x_Toplist=cat(1,x_Toplist,x);
        else
            xmissing=cat(1,xmissing,BandN); %Safetycheck
        end
        
        y=load(y_filename);
        if ~isempty(y)
            y = [y,BandN*ones(length(y(:,1)),1)];
            [~,I]=sort(y(:,5),'descend'); y=y(I(1:Ntop),:); %Reduce toplist
            y_Toplist=cat(1,y_Toplist,y);
        else
            ymissing=cat(1,ymissing,BandN); %Safetycheck
        end
        
    end
    File_cn=sprintf('data_top_%g.mat',grup);
    save (File_cn,'x_Toplist','y_Toplist','xmissing','ymissing') % Save full group
    
    %% COINCIDENCES & GENERATED CENTRE LOOP
    g_Toplist=[];
    kmax=ceil((f0_grup_max-f0_grup)/f0_band-1);
    k0=0;
    if ~isempty(x_Toplist) && ~isempty(y_Toplist)
        
        x1_y = x_Toplist(:,1) - x_Toplist(:,4) * (x_tref - y_tref); % Translate frequency from one data set reftime to the other
        % Extract a portion f0_rang_grup_band of the candidates based on there frequency
        % to reduce the height of the toplist on the following step.
        % (We add a wings considering r_W/Tcoh the maximum possible distance in frequency)
        
        for k=1:kmax;
            if (k/(f0_rang_grup_band/f0_band)) == fix(k/(f0_rang_grup_band/f0_band)) || k==1;
                f0_rang_grup_band/f0_band
                I_rang_grup_band_x= find( f0_grup + k0*f0_rang_grup_band - r_W/Tcoh <= x1_y & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= x1_y);
                I_rang_grup_band_y= find( f0_grup + k0*f0_rang_grup_band <= y_Toplist(:,1) & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= y_Toplist(:,1));
                k0=k0+1;
            end
            
            % Extract a f0_band from the f0_rang_grup_band extraction to compute the coincidences
            % (We add wings considering r_W/Tcoh the maximum possible distance in frequency)
            
            L_band_x= ( f0_grup + k*f0_band - r_W/Tcoh <= x1_y(I_rang_grup_band_x,1)) &  (f0_grup + (k+1)*f0_band +  r_W/Tcoh >= x1_y(I_rang_grup_band_x,1));
            L_band_y= ( f0_grup + k*f0_band            <= y_Toplist(I_rang_grup_band_y,1)) &  (f0_grup + (k+1)*f0_band        >= y_Toplist(I_rang_grup_band_y,1));
            
            if sum(L_band_x) && sum(L_band_y)
                
                % We generate two toplist from the extracion
                x=[x1_y(I_rang_grup_band_x(L_band_x)),x_Toplist(I_rang_grup_band_x(L_band_x),2:end)];  y=y_Toplist(I_rang_grup_band_y(L_band_y),:);
                
                I_sigchi2_x=1:length(x(:,1));
                I_sigchi2_y=1:length(y(:,1));
                
                chi2_STD_x=chi2_STD(x(:,5),x(:,7),p);
                chi2_STD_y=chi2_STD(y(:,5),y(:,7),p);
                
                % We extract the veto and threshold candidates from the toplist
                
                I_sigchi2_x(chi2_STD_x>=chisquare_STD_veto | x(:,5)<sigma_veto_toplist(x(:,1)))=[];
                I_sigchi2_y(chi2_STD_y>=chisquare_STD_veto | y(:,5)<sigma_veto_toplist(y(:,1)))=[];
                
                if ~isempty(I_sigchi2_x) && ~isempty(I_sigchi2_y)
                    
                    % overwrite a toplist without the vetoed and thresholded
                    
                    x=x(I_sigchi2_x,:);  y=y(I_sigchi2_y,:);
                    
                    % Calculate the distace between candidates
                    
                    [m_D,m_D5,n0,n1]=Distance_candidates(x,y,df0,df1,Tcoh,PixelFactor);
                    [I_sigchi2_y_cn,I_sigchi2_x_cn]=find(m_D<r_W & m_D5<rp_sig); % candidates X,Y inside a radius r_W of each other
                    
                    if ~isempty(I_sigchi2_x_cn) && ~isempty(I_sigchi2_y_cn)
                        
                        % Calculate the harmonic mean for the significance, between all the coincidental candidates.
                        s_harm=2*(x(I_sigchi2_x_cn,5).*y(I_sigchi2_y_cn,5)./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5)));
                        
                        % Reduce the coincidences by sorting with the harmonic significance and thresholding by toplist max # of candidates
                        [~,I]=sort(s_harm,'descend');
                        if length(I)>Ntop; I=I(1:Ntop,:); end;
                        
                        I_sigchi2_x_cn=I_sigchi2_x_cn(I); I_sigchi2_x_cn=reshape(I_sigchi2_x_cn,1, numel(I_sigchi2_x_cn))';
                        I_sigchi2_y_cn=I_sigchi2_y_cn(I); I_sigchi2_y_cn=reshape(I_sigchi2_y_cn,1, numel(I_sigchi2_y_cn))';
                        
                        % Calculate the centre between each coincidetal candidates weighted by the significance
                        gf0=(x(I_sigchi2_x_cn,5).*x(I_sigchi2_x_cn,1)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,1))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                        
                        gn=(y(I_sigchi2_y_cn,5)*[1,1,1].*n0(I_sigchi2_y_cn,:)+x(I_sigchi2_x_cn,5)*[1,1,1].*n1(I_sigchi2_x_cn,:))./((x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5))*[1,1,1]);
                        gn=gn./(sqrt(sum(gn'.^2)')*[1,1,1]);
                        
                        galpha  = atan2(gn(:,2),gn(:,1));
                        gdelta  = real(asin(gn(:,3)));
                        
                        gf1=(x(I_sigchi2_x_cn,5).*x(I_sigchi2_x_cn,4)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,4))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                        s_mean=(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5))/2;
                        
                        % Array containing the results and the coincidental index, adding the results for each iteration
                        g_cn_data0 = [gf0,galpha,gdelta,gf1,s_mean,s_harm(I),I_sigchi2_x_cn,I_sigchi2_y_cn];
                        g_Toplist=cat(1,g_Toplist,g_cn_data0);
                        
                        g_cn_data0=[];
                    end
                end
            end
        end
    end
    
    
    %% PARTIAL CLUSTER IDENTIFICATION
    if ~isempty(g_Toplist)
        
        % Sort by frequency the results of the full group and save the generated Toplist
        [~,I]=sort(g_Toplist(:,1));  g_Toplist=g_Toplist(I,:);
        File_cn=sprintf('g_Toplist_%g.mat',grup);
        save (File_cn,'g_Toplist')
        
        k0=0;
        I_N_Cluster=0;
        I_partial_cluster=[];
        I_partial_cluster{1}=[];
        for k=1:kmax;
            
            % Extract a portion f0_rang_grup_band of the candidates based on there frequency
            % to reduce the height of the toplist on the following step.
            % (We add a wings considering r_W/Tcoh the maximum possible distance in frequency)
            
            if (k/(f0_rang_grup_band/f0_band)) == fix(k/(f0_rang_grup_band/f0_band)) || k==1;
                f0_grup + k0*f0_rang_grup_band;
                I_rang_grup_band_x= find( f0_grup + k0*f0_rang_grup_band - r_W/Tcoh <= g_Toplist(:,1) & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= g_Toplist(:,1));
                k0=k0+1;
            end
            
            gx_Toplist=g_Toplist(I_rang_grup_band_x,:);
            
            % Extract a f0_band from the f0_rang_grup_band extraction to compute the coincidences
            % (We add wings considering r_W/Tcoh the maximum possible distance in frequency)
            
            I_band =  find( f0_grup + k*f0_band - r_W_cl/Tcoh <= gx_Toplist(:,1) &  f0_grup + (k+1)*f0_band +  r_W_cl/Tcoh >= gx_Toplist(:,1));
            
            if ~isempty(I_band)
                
                x = gx_Toplist(I_band,:);
                min(x(:,1))
                % Calculate the distace between candidates
                [m_D,m_D5,~,~]=Distance_candidates(x,x,df0,df1,Tcoh,PixelFactor);
                
                % We get a matrix with a 1 for candidates that exist in a distance r_W_cl
                % of each other and 0 for the rest.
                L_m_D=(m_D<r_W_cl & m_D5<rp_sig);
                ones_x = ones(length(x(:,1)),1);
                
                % For matrix I_m each row correspond to candidate j and each column to candidate i.
                % If the condition is not satisfy we will set a 0, else the index i.
                I_m = ones_x * (1:1:size(ones_x,1)).*L_m_D;
                i0=[];
                for j=1:length(I_m(:,1))
                    I_v_ji_0=[]; I_v_ji=[];
                    % Jump*
                    if ~ismember(j,i0)
                        % We extract the neighbours for the candidates on the row j
                        I_v_ji=I_m(j,:);
                        I_v_ji(I_v_ji==0)=[];
                        % We will continue extracting until we find no new rows to extract from
                        while length(I_v_ji_0)~=length(I_v_ji);
                            I_v_ji_0=I_v_ji;
                            B=reshape(I_m(I_v_ji,:).' ,1,numel(I_m(I_v_ji,:)));
                            % Extract the neighbours for all the rows that
                            % were on j or on any subsequent extractions
                            I_v_ji=unique(B);
                            I_v_ji(I_v_ji==0)=[];
                        end
                        % We set a list for the rows we have extracted/been so
                        % we can (Jump*) avoid redoing the processes.
                        i0=cat(2,i0,I_v_ji);
                        % COUNT the # of partial clusters we have analyses per group
                        I_N_Cluster=I_N_Cluster+1;
                        % We translate the I_v_ji neighbours list to the generated toplist index.
                        I_partial_cluster{I_N_Cluster}=(I_band(I_v_ji))';
                    end
                end
            end
        end
        
        %% CLOUSING CLUSTERS
        N=1;
        x=g_Toplist;
        % First we define each centre in one cluster id 0
        N_Cluster=zeros(length(x(:,1)),1);
        I_Ci=[];
        % We extract the first partition*
        A=I_partial_cluster{1};
        
        if I_N_Cluster>1;
            for i=1:I_N_Cluster;
                if isempty(I_Ci);
                    % Set the cluster id N to partition A
                    N_Cluster(A)=N;
                    N=N+1;
                    N_Cluster(A)
                else
                    % Set a common id to partition A and any possible
                    % connection with different id assigned.
                    %(if any of the centres had been already assigned or that
                    % now are linking multiple clusters with different id
                    % assigned)
                    N_Cluster(A)=N_Cluster_i(I_Ci(1));
                    B = ismember(N_Cluster,N_Cluster_i(I_Ci));
                    N_Cluster(B)=N_Cluster_i(I_Ci(1));
                end
                if i<I_N_Cluster;
                    A=I_partial_cluster{i+1}; % * Subsecuent partition extractions
                    % Extract the assigned id
                    N_Cluster_i=unique(N_Cluster(A));
                    % Collect id diffect than 0
                    I_Ci = find(N_Cluster_i~=0);
                end
            end
        else
            % if we only have one partition
            N_Cluster(A)=1;
        end
        
        %% CLUSTER CENTRE CALCULATION
        j=1;
        for j=1:max(N_Cluster);
            
            I_Cluster_i=find(N_Cluster==j);
            s_mean=x(I_Cluster_i,6); % from the harmonic significance generated
            
            Cluster.sign_mean{j} =mean(x(I_Cluster_i,5));
            Cluster.ind{j}       =I_Cluster_i;
            Cluster.s_harm{j}    =mean(s_mean);
            Cluster.sign_sum{j}  =sum(s_mean);
            Cluster.parent_x{j}  =length(unique(x(I_Cluster_i,7)));
            Cluster.parent_y{j}  =length(unique(x(I_Cluster_i,8)));
            Cluster.f0{j}        =sum(x(I_Cluster_i,1).*s_mean)/Cluster.sign_sum{j};
            
            n0=sum(s_mean*[1,1,1].*[cos(x(I_Cluster_i,3)).*cos(x(I_Cluster_i,2)),cos(x(I_Cluster_i,3)).*sin(x(I_Cluster_i,2)),sin(x(I_Cluster_i,3))],1);
            n0=n0./sqrt(sum(n0.^2));
            
            Cluster.alpha{j}     =real(atan2(n0(:,2),n0(:,1)));
            Cluster.delta{j}     =real(asin(n0(:,3)));
            
            Cluster.f1{j}        =sum(x(I_Cluster_i,4).*s_mean)/Cluster.sign_sum{j};
            Cluster.length{j}    =length(I_Cluster_i);
            
            Cluster.sign_max{j}  =max(s_mean);
            
            if length(I_Cluster_i)<2; Cluster.noise{j}=1; else Cluster.noise{j}=0; end
        end
        follow_up=[[Cluster.f0{:}]',[Cluster.alpha{:}]',[Cluster.delta{:}]',[Cluster.f1{:}]',[Cluster.sign_mean{:}]',[Cluster.s_harm{:}]',...
            [Cluster.sign_sum{:}]',[Cluster.length{:}]',[Cluster.parent_x{:}]',[Cluster.parent_y{:}]'];
    end
    
    follow_up=[[Cluster.f0{:}]',[Cluster.alpha{:}]',[Cluster.delta{:}]',[Cluster.f1{:}]',[Cluster.sign_mean{:}]',[Cluster.s_harm{:}]',...
        [Cluster.sign_sum{:}]',[Cluster.length{:}]',[Cluster.parent_x{:}]',[Cluster.parent_y{:}]'];
    
    %% FOLLOW UP ASSIGMENT
    % veto population
    if ~isempty(follow_up)
        follow_up=follow_up(follow_up(:,7)>population_cluster_veto,:);
        follow_up=follow_up(follow_up(:,5)>significance_veto_cluster,:);
        
        f0=min(follow_up(:,1)):0.1:max(follow_up(:,1));
        follow_up0=[];
        for i=Extract_S;min_band
            Is=[];
            [~,Iz]=sort(follow_up(:,i),'descend'); follow_up=follow_up(Iz,:);
            for BandN=1:length(f0);
                I=find(abs(f0(BandN)-follow_up(:,1))<min_band);
                Is=[Is,min(I)];
            end
            follow_up0=[follow_up0;[follow_up(Is,:),i*ones(length(Is),1)]];
        end
        [~,I]=sort(follow_up0(:,1)); follow_up=follow_up0(I,:);
        
        % sigma veto cluster
        
        param=[sig_veto,chisquare_STD_veto,r_W,r_W_cl,population_cluster_veto,significance_veto_cluster,min_band];
    end
    follow_up
    File_cn=sprintf('follow_up_%g.mat',grup);
    save (File_cn,'follow_up','Cluster','param')
end

%% DISTACE FUCTION
function [m_D,m_D5,n0,n1] = Distance_candidates(x,y,df0,df1,Tcoh,PixelFactor)

m_W1 = df0; % Calculate distance the bins size freq
m_W4 = sum(df1)/2; % Calculate distance the bins size spin

[X,Y]=meshgrid(x(:,1),y(:,1));
m_D1 = (X-Y)/m_W1; % Distance in bins for the frequency between each candidate
m_W2 = (2./(Y + X))/(Tcoh*1e-4)/PixelFactor; % Calculate distance the bins sky

[X,Y]=meshgrid(x(:,4),y(:,4));
m_D4 = (X-Y)/m_W4; % Distance in bins for the spindown between each candidate

n0=[cos(y(:,3)).*cos(y(:,2)),cos(y(:,3)).*sin(y(:,2)),sin(y(:,3))];
n1=[cos(x(:,3)).*cos(x(:,2)),cos(x(:,3)).*sin(x(:,2)),sin(x(:,3))];

b=sqrt((n0(:,2)*n1(:,3)'-n0(:,3)*n1(:,2)').^2+(n0(:,3)*n1(:,1)'-n0(:,1)*n1(:,3)').^2+(n0(:,1)*n1(:,2)'-n0(:,2)*n1(:,1)').^2);
a=n0(:,1)*n1(:,1)'+n0(:,2)*n1(:,2)'+n0(:,3)*n1(:,3)';

m_D2 = atan2(b,a)./m_W2; % Distance in bins for the sky between each candidate
% Considering it the minimum distance over a great circle.

[X,Y]=meshgrid(x(:,5),y(:,5));
m_D5 = abs((X - Y)./(X + Y));

m_D = sqrt(m_D1.^2+m_D2.^2+m_D4.^2); % the total distance

end