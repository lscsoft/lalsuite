%{
 *  Copyright (C) 2018 Miquel Oliver, Rodrigo Tenorio, Pep Blai Covas
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA2
 %}

function post_processing_O2_ecl(X_input,Y_input,grup,grup_Band,BandNo,BandNf,Ntop,PixelFactor,p,sig_veto,rp_sig,r_W,r_W_cl,CHI_V,chisquare_STD_veto,population_cluster_veto,significance_veto_cluster,Extract_S,min_band,Clean_lines_I,Clean_lines_II,load_data,file1,file2,file3)

format long

if grup<300
    Tcoh = 3600;
    x_tref = 1167545503; %Detector == 'FirstSet':
    y_tref = 1179817563; %Detector == 'SecondSet':
elseif grup>=300 && grup<550
    Tcoh = 2700;
    x_tref = 1167545053; %Detector == 'FirstSet':
    y_tref = 1179817113; %Detector == 'SecondSet':
elseif grup>=550 && grup<1300
    Tcoh = 1800;
    x_tref = 1167545839; %Detector == 'FirstSet':
    y_tref = 1179816663; %Detector == 'SecondSet':
else
    Tcoh = 900;
    x_tref = 1167545053; %Detector == 'FirstSet':
    y_tref = 1179814469; %Detector == 'SecondSet':
end

%Tobs=[x_tend,y_tend]-[x_tref,y_tref];
Tobs=7919026;
df0=(1/Tcoh);
df1=(df0./Tobs);

f0_band=r_W*df0; %%% size of the Band to extract and analyse
f0_rang_grup_band=6000; %%% size of the Band to extract before analysis

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
[I,~]=find(CHI_V(:,1) <= f0_grup);
I=max(I);
A1=CHI_V(I,2);
A2=CHI_V(I,3);
B1=CHI_V(I,4);
B2=CHI_V(I,5);
chi2_STD = @(x,y,p) (y-(p-1)-A1*x.^A2)./(sqrt(2*p-2)+B1*x.^B2);

phi = 23.4 * pi / 180;

%% READING TOPLISTS
%if ~exist(file1) || load_data==0;
if load_data==0
    for BandN=BandNo:BandNf;
        if ((BandN/10) == fix(BandN/10)) || BandN==0;
            disp(['READING: ',num2str(BandN),'/',num2str(BandNf),'   ',datestr(now)])
        end
        
        y_filename=sprintf(Y_input,grup,BandN);
        x_filename=sprintf(X_input,grup,BandN);
        x = [];
        y = [];
        
        if exist(x_filename, 'file')
            x=load(x_filename);
            if ~isempty(x)
                x = [x,BandN*ones(length(x(:,1)),1)];
                %[~,I]=sort(x(:,5),'descend');
		%x=x(I(1:Ntop),:);
                x=x(1:Ntop,:);
                x_Toplist=cat(1,x_Toplist,x);
            else
                xmissing=cat(1,xmissing,BandN);
            end
        else
            xmissing=cat(1,xmissing,BandN);
            str = sprintf('FirstSet missing for %d',BandN)
        end
        
        if exist(y_filename, 'file')
            y=load(y_filename);
            if ~isempty(y)
                y = [y,BandN*ones(length(y(:,1)),1)];
                %[~,I]=sort(y(:,5),'descend');
		%y=y(I(1:Ntop),:);
                ll = length(y);
                if ll>=Ntop
                y=y(1:Ntop,:);
                end
                y_Toplist=cat(1,y_Toplist,y);
            else
                ymissing=cat(1,ymissing,BandN);
            end
        else
            ymissing=cat(1,ymissing,BandN);
            str	= sprintf('SecondSet missing for %d',BandN)
        end
        
    end

	y_filename=sprintf(Y_input,grup-50,499);
        x_filename=sprintf(X_input,grup-50,499);
       	x = [];
        y = [];

        if exist(x_filename, 'file')
            x=load(x_filename);
            if ~isempty(x)
                x = [x,-1*ones(length(x(:,1)),1)];
                %[~,I]=sort(x(:,5),'descend');
                %x=x(I(1:Ntop),:);
                x=x(1:Ntop,:);
                x_Toplist=cat(1,x_Toplist,x);
            else
                xmissing=cat(1,xmissing,BandN);
            end
        else
            xmissing=cat(1,xmissing,BandN);
        end

	if exist(y_filename, 'file')
            y=load(y_filename);
            if ~isempty(y)
                y = [y,-1*ones(length(y(:,1)),1)];
                %[~,I]=sort(y(:,5),'descend');
                %y=y(I(1:Ntop),:);
                y=y(1:Ntop,:);
                y_Toplist=cat(1,y_Toplist,y);
            else
                ymissing=cat(1,ymissing,BandN);
            end
        else
            ymissing=cat(1,ymissing,BandN);
        end

	y_filename=sprintf(Y_input,grup+50,0);
        x_filename=sprintf(X_input,grup+50,0);
       	x = [];
        y = [];

        if exist(x_filename, 'file')
            x=load(x_filename);
            if ~isempty(x)
                x = [x,500*ones(length(x(:,1)),1)];
                %[~,I]=sort(x(:,5),'descend');
                %x=x(I(1:Ntop),:);
                x=x(1:Ntop,:);
                x_Toplist=cat(1,x_Toplist,x);
            else
                xmissing=cat(1,xmissing,BandN);
            end
        else
            xmissing=cat(1,xmissing,BandN);
        end

	if exist(y_filename, 'file')
            y=load(y_filename);
            if ~isempty(y)
                y = [y,500*ones(length(y(:,1)),1)];
                %[~,I]=sort(y(:,5),'descend');
                %y=y(I(1:Ntop),:);
                y=y(1:Ntop,:);
                y_Toplist=cat(1,y_Toplist,y);
            else
                ymissing=cat(1,ymissing,BandN);
            end
        else
            ymissing=cat(1,ymissing,BandN);
        end

    disp(['Save: ',file1,'   ',datestr(now)])
    save (file1,'x_Toplist','y_Toplist','xmissing','ymissing') % Save full group
else
    load(file1)
    %x_Toplist = x_Toplist(1:Ntop,:);
    %y_Toplist=y_Toplist(1:Ntop,:);
end


%% COINCIDENCES & GENERATED CENTRE LOOP
g_Toplist=[];
kmax=ceil((f0_grup_max-f0_grup)/f0_band - 1);
k0=0;
if ~isempty(x_Toplist) && ~isempty(y_Toplist)
    
    %x1_y = x_Toplist(:,1) - x_Toplist(:,4) * (x_tref - y_tref); % Translate frequency from one data set reftime to the other
    y1_x = y_Toplist(:,1) - y_Toplist(:,4) * (y_tref - x_tref);

    % Extract a portion f0_rang_grup_band of the candidates based on there frequency to reduce the height of the toplist on the following step.
    % (We add a wings considering r_W/Tcoh the maximum possible distance in frequency)    
    for k=1:kmax+1;
        
        if (k/(f0_rang_grup_band)) == fix(k/(f0_rang_grup_band)) || k==1 || k==kmax
            disp(['COINCIDENCES: ',num2str(k),'/',num2str([kmax,grup+f0_band*(k-1)]),'(Hz)   ',datestr(now)])
            I_rang_grup_band_x= find( f0_grup + k0*floor(f0_rang_grup_band*f0_band) - f0_band <= x_Toplist(:,1) & f0_grup + (k0+1)*ceil(f0_rang_grup_band*f0_band) + f0_band >= x_Toplist(:,1));
            I_rang_grup_band_y= find( f0_grup + k0*floor(f0_rang_grup_band*f0_band)           <= y1_x           & f0_grup + (k0+1)*ceil(f0_rang_grup_band*f0_band) >= y1_x);
            k0=k0+1;
        end
        
        % Extract a f0_band from the f0_rang_grup_band extraction to compute the coincidences (We add wings considering r_W/Tcoh the maximum possible distance in frequency)
        L_band_x= ( f0_grup + (k-1)*f0_band - f0_band  <= x_Toplist(I_rang_grup_band_x,1))      &  (f0_grup + (k)*f0_band +  f0_band >= x_Toplist(I_rang_grup_band_x,1));
        L_band_y= ( f0_grup + (k-1)*f0_band            <= y1_x(I_rang_grup_band_y,1)) &  (f0_grup + (k)*f0_band            >= y1_x(I_rang_grup_band_y,1));
        
        if sum(L_band_x) && sum(L_band_y)
            
            % We generate two toplist from the extracion
            %x=[x1_y(I_rang_grup_band_x(L_band_x)),x_Toplist(I_rang_grup_band_x(L_band_x),2:end)];  y=y_Toplist(I_rang_grup_band_y(L_band_y),:);
            x=x_Toplist(I_rang_grup_band_x(L_band_x),:);  y=[y1_x(I_rang_grup_band_y(L_band_y)),y_Toplist(I_rang_grup_band_y(L_band_y),2:end)];
            
            I_sigchi2_x=1:length(x(:,1));
            I_sigchi2_y=1:length(y(:,1));
            
            chi2_STD_x=chi2_STD(x(:,5),x(:,7),p);
            chi2_STD_y=chi2_STD(y(:,5),y(:,7),p);
            
            % We extract the veto and threshold candidates from the toplist
            
            %I_sigchi2_x(chi2_STD_x>=chisquare_STD_veto | x(:,5)<sigma_veto_toplist(x(:,1)))=[];
            %I_sigchi2_y(chi2_STD_y>=chisquare_STD_veto | y(:,5)<sigma_veto_toplist(y(:,1)))=[];
            I_sigchi2_x(x(:,5)<sigma_veto_toplist(x(:,1)))=[];
            I_sigchi2_y(y(:,5)<sigma_veto_toplist(y(:,1)))=[];            

            if ~isempty(I_sigchi2_x) && ~isempty(I_sigchi2_y)
                
                % overwrite a toplist without the vetoed and thresholded
                
                x=x(I_sigchi2_x,:);  y=y(I_sigchi2_y,:);
                
                % Calculate the distace between candidates
                
                [m_D,m_D5]=Distance_candidates(x,y,df0,df1,Tcoh,PixelFactor,r_W);
                [I_sigchi2_y_cn,I_sigchi2_x_cn]=find(m_D<r_W); % & m_D5<rp_sig); % candidates X,Y inside a radius r_W of each other
                
                if ~isempty(I_sigchi2_x_cn) && ~isempty(I_sigchi2_y_cn)
                    
                    % Calculate the harmonic mean for the significance, between all the coincidental candidates.
                    s_harm=2*(x(I_sigchi2_x_cn,5).*y(I_sigchi2_y_cn,5)./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5)));
                    
                    % Reduce the coincidences by sorting with the harmonic significance and thresholding by toplist max # of candidates
                    [~,I]=sort(s_harm,'descend');
                    %if length(I)>Ntop; I=I(1:Ntop,:); end;
                    
                    I_sigchi2_x_cn=I_sigchi2_x_cn(I); I_sigchi2_x_cn=reshape(I_sigchi2_x_cn,1, numel(I_sigchi2_x_cn))';
                    I_sigchi2_y_cn=I_sigchi2_y_cn(I); I_sigchi2_y_cn=reshape(I_sigchi2_y_cn,1, numel(I_sigchi2_y_cn))';
                    
                    % Calculate the centre between each coincidetal candidates weighted by the significance
                    gf0=(x(I_sigchi2_x_cn,5).*x(I_sigchi2_x_cn,1)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,1))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                    
                    x_xecl = cos(x(I_sigchi2_x_cn,2)).*cos(x(I_sigchi2_x_cn,3));
                    y_xecl = cos(y(I_sigchi2_y_cn,2)).*cos(y(I_sigchi2_y_cn,3));                        
                    x_yecl = sin(x(I_sigchi2_x_cn,2)).*cos(x(I_sigchi2_x_cn,3)).*cos(phi) + sin(phi).*sin(x(I_sigchi2_x_cn,3));
                    y_yecl = sin(y(I_sigchi2_y_cn,2)).*cos(y(I_sigchi2_y_cn,3)).*cos(phi) + sin(phi).*sin(y(I_sigchi2_y_cn,3));                       
                    x_zecl = sin(x(I_sigchi2_x_cn,3)).*cos(phi) - sin(phi).*cos(x(I_sigchi2_x_cn,3)).*sin(x(I_sigchi2_x_cn,2));
                    y_zecl = sin(y(I_sigchi2_y_cn,3)).*cos(phi) - sin(phi).*cos(y(I_sigchi2_y_cn,3)).*sin(y(I_sigchi2_y_cn,2));
                        
                    gx = (x_xecl.*x(I_sigchi2_x_cn,5) + y_xecl.*y(I_sigchi2_y_cn,5))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                    gy = (x_yecl.*x(I_sigchi2_x_cn,5) + y_yecl.*y(I_sigchi2_y_cn,5))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                    gz = (x_zecl.*x(I_sigchi2_x_cn,5) + y_zecl.*y(I_sigchi2_y_cn,5))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
	       	    norm = sqrt(gx.^2 + gy.^2 + gz.^2);
		    gx = gx./norm;
                    gy = gy./norm;
                    gz = gz./norm;                    

                    gf1=(x(I_sigchi2_x_cn,5).*x(I_sigchi2_x_cn,4)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,4))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                    s_mean=(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5))/2;
                    
                    % Array containing the results and the coincidental index, adding the results for each iteration
%                    g_cn_data0 = [gf0,gx,gy,gz,gf1,s_mean,s_harm(I),I_sigchi2_x_cn,I_sigchi2_y_cn,x(I_sigchi2_x_cn,2),y(I_sigchi2_y_cn,2),x(I_sigchi2_x_cn,3),y(I_sigchi2_y_cn,3)];
		     g_cn_data0 = [gf0,gx,gy,gz,gf1,s_mean,s_harm(I),I_sigchi2_x_cn,I_sigchi2_y_cn,x(I_sigchi2_x_cn,1),y(I_sigchi2_y_cn,1),x(I_sigchi2_x_cn,2),y(I_sigchi2_y_cn,2),x(I_sigchi2_x_cn,3),y(I_sigchi2_y_cn,3),x(I_sigchi2_x_cn,4),y(I_sigchi2_y_cn,4),x(I_sigchi2_x_cn,10),y(I_sigchi2_y_cn,10)];
                    g_Toplist=cat(1,g_Toplist,g_cn_data0);
                    
                    g_cn_data0=[];
               end
           end
        end
    end
end

save (file2,'g_Toplist')

%% PARTIAL CLUSTER IDENTIFICATION
aaaab=1;
%if aaaab==2
if ~isempty(g_Toplist)
    
    % Sort by frequency the results of the full group and save the generated Toplist
    [~,I]=sort(g_Toplist(:,1));  g_Toplist=g_Toplist(I,:);
    disp(['Save: ',file2,'   ',datestr(now)])
    %save (file2,'g_Toplist')
    
    k0=0;
    I_N_Cluster=0;
    I_partial_cluster=[];
    I_partial_cluster{1}=[];
    for k=1:kmax+1;
        
        % Extract a portion f0_rang_grup_band of the candidates based on there frequency
        % to reduce the height of the toplist on the following step.
        % (We add a wings considering r_W/Tcoh the maximum possible distance in frequency)
        
        if (k/(f0_rang_grup_band)) == fix(k/(f0_rang_grup_band)) || k==1 || k==kmax
            disp(['CLUSTER: ',num2str(k),'/',num2str([kmax,grup+f0_band*k]),'(Hz)   ',datestr(now)])
            I_rang_grup_band_x= find( f0_grup + k0*floor(f0_rang_grup_band*f0_band) - r_W_cl/Tcoh <= g_Toplist(:,1) & f0_grup + (k0+1)*ceil(f0_rang_grup_band*f0_band) + r_W_cl/Tcoh >= g_Toplist(:,1));
            k0=k0+1;
        end
        
        gx_Toplist=g_Toplist(I_rang_grup_band_x,:);
        
        % Extract a f0_band from the f0_rang_grup_band extraction to compute the coincidences
        % (We add wings considering r_W/Tcoh the maximum possible distance in frequency)
        
        I_band =  find( f0_grup + (k-1)*f0_band - r_W_cl/Tcoh <= gx_Toplist(:,1) &  f0_grup + (k)*f0_band +  r_W_cl/Tcoh >= gx_Toplist(:,1));
        
        if ~isempty(I_band)
            x = gx_Toplist(I_band,:);
            % Calculate the distace between candidates
            [m_D,m_D5]=Distance_candidates_b(x,x,df0,df1,Tcoh,PixelFactor,r_W_cl);
            
            % We get a matrix with a 1 for candidates that exist in a distance r_W_cl of each other and 0 for the rest.
            L_m_D=(m_D<r_W_cl); % & m_D5<rp_sig);
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
                        % Extract the neighbours for all the rows that were on j or on any subsequent extractions
                        I_v_ji=unique(B);
                        I_v_ji(I_v_ji==0)=[];
                    end
                    % We set a list for the rows we have extracted/been so we can (Jump*) avoid redoing the processes.
                    i0=cat(2,i0,I_v_ji);
                    % COUNT the # of partial clusters we have analyses per group
                    I_N_Cluster=I_N_Cluster+1;
                    % We translate the I_v_ji neighbours list to the generated toplist index.
                    I_partial_cluster{I_N_Cluster}=(I_rang_grup_band_x(I_band(I_v_ji)))';
                end
            end
        end
    end
    
    %% CLOSING CLUSTERS
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
    id_clust=1;
    for j=1:max(N_Cluster);
        I_Cluster_i=find(N_Cluster==j);
        if ~isempty(I_Cluster_i)
            s_mean=x(I_Cluster_i,7); % from the harmonic significance generated
            Cluster.id{j} = id_clust;
            Cluster.sign_mean{j} =mean(x(I_Cluster_i,6));
            Cluster.ind{j}       =I_Cluster_i;
            Cluster.s_harm{j}    =mean(s_mean);
            Cluster.sign_sum{j}  =sum(s_mean);
            A_max=max(s_mean);
            if isempty(A_max);
                A_max=NaN;
            end
            Cluster.sign_max{j}  = A_max;
            Cluster.parent_x{j}  =length(unique(x(I_Cluster_i,8)));
            Cluster.parent_y{j}  =length(unique(x(I_Cluster_i,9)));
            Cluster.f0{j}        =sum(x(I_Cluster_i,1).*s_mean)/Cluster.sign_sum{j};

            gx_ecu = sum(s_mean.*x(I_Cluster_i,2))/Cluster.sign_sum{j};
            gy_ecu = sum(s_mean.*(cos(phi) .* x(I_Cluster_i,3) - sin(phi) .* x(I_Cluster_i,4)))/Cluster.sign_sum{j};
            gz_ecu = sum(s_mean.*(cos(phi) .* x(I_Cluster_i,4) + sin(phi) .* x(I_Cluster_i,3)))/Cluster.sign_sum{j};
            gr = sqrt(gx_ecu.^2 + gy_ecu.^2 + gz_ecu.^2);
	    gx_ecu = gx_ecu./gr;
            gy_ecu = gy_ecu./gr;
            gz_ecu = gz_ecu./gr;
            Cluster.alpha{j} = atan2(gy_ecu,gx_ecu);
            Cluster.delta{j} = asin(gz_ecu);

            Cluster.f1{j}        =sum(x(I_Cluster_i,5).*s_mean)/Cluster.sign_sum{j};
            Cluster.length{j}    =length(I_Cluster_i);
            if length(I_Cluster_i)<2;
                Cluster.noise{j}=1;
            else
                Cluster.noise{j}=0;
            end
        else
            id_clust=id_clust-1;
        end
        id_clust=id_clust+1;
    end
    
    follow_up=[[Cluster.f0{:}]',[Cluster.alpha{:}]',[Cluster.delta{:}]',[Cluster.f1{:}]',[Cluster.sign_mean{:}]',[Cluster.s_harm{:}]',...
        [Cluster.sign_sum{:}]',[Cluster.length{:}]',[Cluster.parent_x{:}]',[Cluster.parent_y{:}]',[Cluster.id{:}]'];
end

%% FOLLOW UP ASSIGMENT
disp(['Save: ',file3,'   ',datestr(now)])
param=[sig_veto,chisquare_STD_veto,r_W,r_W_cl,population_cluster_veto,significance_veto_cluster,min_band];

%if aaaab==2
if ~isempty(follow_up)
    %% WITH LINE CLEANUP
    follow_up_clean=[];
    try
        F_list=[];
        A=load(Clean_lines_II);
        for i=1:size(A,1);
            f0=(A(i,1)*[A(i,4):A(i,5)]'+A(i,3));
            F_list=[F_list;f0-A(i,6)-20*df0,f0+A(i,7)+20*df0];
        end
        A=load(Clean_lines_I);
        for i=1:size(A,1);
            f0=(A(i,1)*[A(i,4):A(i,5)]'+A(i,3));
            F_list=[F_list;f0-A(i,6)-20*df0,f0+A(i,7)+20*df0];
        end
        follow_up_clean=follow_up;
        for i=1:size(follow_up_clean,1)
            if sum(follow_up_clean(i,1) >= F_list(:,1) & follow_up_clean(i,1) <= F_list(:,2));
                follow_up_clean(i,:) = [];
            end
        end
        follow_up_clean = follow_up_clean(follow_up_clean(:,8)>population_cluster_veto(1) & follow_up_clean(:,9)>population_cluster_veto(2) & follow_up_clean(:,10)>population_cluster_veto(3),:);
        follow_up_clean = follow_up_clean(follow_up_clean(:,6)>significance_veto_cluster,:);
        follow_up_clean = follow_up_selection(follow_up_clean,Extract_S,min_band);
    catch
        %'no lines file'
    end
    
    %% WITHOUT CLEANUP
    follow_up = follow_up(follow_up(:,8)>=population_cluster_veto(1) & follow_up(:,9)>=population_cluster_veto(2) & follow_up(:,10)>=population_cluster_veto(3),:);
    follow_up = follow_up(follow_up(:,6)>significance_veto_cluster,:);
    follow_up = follow_up_selection(follow_up,Extract_S,min_band);
    follow_up

    save (file3,'follow_up','follow_up_clean','Cluster','param')
else
    %save (file3,'Cluster','param')
	aaaaa=1;
end

%exit
end


%% FOLLOW UP ASSIGMENT FUNCTION
function x = follow_up_selection(x,Extract_S,min_band)
f0_ext=(min(x(:,1))-0.1):0.1:(max(x(:,1))+0.1);
follow_up0=[];
for i=Extract_S;
    Is=[];
    [~,Iz]=sort(x(:,i),'descend');
    x=x(Iz,:);
    for BandN=1:length(f0_ext);
        I=find(abs(f0_ext(BandN)-x(:,1))<min_band);
        Is=[Is,min(I)];
    end
    Is=unique(Is);
    follow_up0=[follow_up0;[x(Is,:),i*ones(length(Is),1)]];
end
[~,I]=sort(follow_up0(:,1));
x=follow_up0(I,:);
end


%% DISTANCE FUCTION
function [m_D,m_D5] = Distance_candidates(x,y,df0,df1,Tcoh,PixelFactor,r_W)

phi = 23.4 * pi / 180;

m_W1 = df0; % Calculate distance the bins size freq
m_W4 = sum(df1)/2; % Calculate distance the bins size spin

[X,Y]=meshgrid(x(:,1),y(:,1));
m_D1 = (X-Y)/m_W1; % Distance in bins for the frequency between each candidate
m_W2 = (2./(Y + X))/(Tcoh*1e-4)/PixelFactor; % Calculate distance the bins sky

[X,Y]=meshgrid(x(:,4),y(:,4));
m_D4 = (X-Y)/m_W4; % Distance in bins for the spindown between each candidate

% Cartesianas Ecuatoriales (x, y)                   
x_ecux = cos(x(:,3)).*cos(x(:,2));
y_ecux = cos(x(:,3)).*sin(x(:,2));
z_ecux = sin(x(:,3));
                   
x_ecuy = cos(y(:,3)).*cos(y(:,2));
y_ecuy = cos(y(:,3)).*sin(y(:,2));
z_ecuy = sin(y(:,3));
                    
% Cartesianas Ecuatoriales -> Cartesianas Eclipticas
% Rotar un alngulo phi (ya definido arriba
x_eclx = x_ecux;
y_eclx = cos(phi) .* y_ecux + sin(phi) .* z_ecux;
z_eclx = cos(phi) .* z_ecux - sin(phi) .* y_ecux;
                             
x_ecly = x_ecuy;
y_ecly = cos(phi) .* y_ecuy + sin(phi) .* z_ecuy;
z_ecly = cos(phi) .* z_ecuy - sin(phi) .* y_ecuy;
                     
% Calculamos distancia en el plano
% SIN CONSIDERAR z_x, z_y
[X_x, X_y] = meshgrid(x_eclx, x_ecly);
[Y_x, Y_y] = meshgrid(y_eclx, y_ecly);

m_D2 = sqrt( (X_x - X_y).^2 + (Y_x - Y_y).^2 )./m_W2;

[X,Y]=meshgrid(x(:,5),y(:,5));
m_D5 = abs((X - Y)./(X + Y));

m_D = sqrt(m_D1.^2+m_D2.^2+m_D4.^2); % the total distance

% Calculamos la separacion en z en bines de distancia
% Si se pasa de r_W, fuera. 

[Z_x, Z_y] = meshgrid(z_eclx, z_ecly);
dz = abs(Z_x - Z_y) ./ m_W2;

m_D(dz > r_W) = r_W + 1.;

end


%% DISTANCE FUCTION B
function [m_D,m_D5] = Distance_candidates_b(x,y,df0,df1,Tcoh,PixelFactor,r_W_cl)

m_W1 = df0; % Calculate distance the bins size freq
m_W4 = sum(df1)/2; % Calculate distance the bins size spin

[X,Y]=meshgrid(x(:,1),y(:,1));
m_D1 = (X-Y)/m_W1; % Distance in bins for the frequency between each candidate
m_W2 = (2./(Y + X))/(Tcoh*1e-4)/PixelFactor; % Calculate distance the bins sky

[X,Y]=meshgrid(x(:,5),y(:,5));
m_D4 = (X-Y)/m_W4; % Distance in bins for the spindown between each candidate

[x_a, x_b] = meshgrid(x(:,2), y(:,2));
[y_a, y_b] = meshgrid(x(:,3), y(:,3));

m_D2 = sqrt( (x_a - x_b).^2 + (y_a-y_b).^2 )./m_W2;

m_D = sqrt(m_D1.^2+m_D2.^2+m_D4.^2); % the total distance

[z_a, z_b] = meshgrid(x(:,4), y(:,4));
dz = abs(z_a - z_b) ./ m_W2;

m_D(dz > r_W_cl ) = r_W_cl + 1.;

[X,Y]=meshgrid(x(:,6),y(:,6));
m_D5 = abs((X - Y)./(X + Y));

end
