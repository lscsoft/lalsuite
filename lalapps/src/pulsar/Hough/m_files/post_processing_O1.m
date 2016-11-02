clear
format long

grup_min=50; %f0_grup=grup;
grup_Band=50;
grup_max=450; %f0_grup_max=grup+grup_Band;

BandNo=1;
BandNf=500;

f0_rang_grup_band=10;
f0_band=0.1;

Ntop=999;

PixelFactor=2;
p = 16; %degeneration

sig=0;
sigma_veto_toplist= @(x) sig+0.*x;

r_W=sqrt(14);
r_W_cl=sqrt(14);

chisquare_STD_veto = 5;

population_cluster_veto=10;

sigma_veto_cluster=1.15;

%Detector == 'L':
x_tref = 1126034147;
x_tend = 1137248624;

%Detector == 'H':
y_tref = 1125972653;
y_tend = 1137250279;

Tcoh=1800;

Tobs=[x_tend,y_tend]-[x_tref,y_tref];

df0=(1/Tcoh);
df1=(df0./Tobs);

Y_input='/home/miquel.oliver/O1/Collected/H1/H1_Collected_%g_%g.dat';
X_input='/home/miquel.oliver/O1/Collected/L1/L1_Collected_%g_%g.dat';





for grup=grup_min:grup_Band:grup_max;
    
    f0_grup=grup;
    f0_grup_max=grup+grup_Band;
    
    x_Toplist=[];
    y_Toplist=[];
    
    xmissing=[];
    ymissing=[];
    
    follow_up=[];
	Cluster=[];
	I_partial_cluster=[];

if f0_grup>=50;
A1=0.4902;
A2=1.414;
B1=0.3581;
B2=1.484;
end

if f0_grup>=100;
A1=0.2168;
A2=1.428;
B1=0.1902;
B2=1.499;

end	

if f0_grup>=200;
A1 = 0.1187;
A2 = 1.470;
B1 = 0.06784;
B2 = 1.697;

end
	
	chi2_STD = @(x,y,p) (y-(p-1)-A1*x.^A2)./(sqrt(2*p-2)+B1*x.^B2);
 
    for BandN=BandNo:BandNf;
        
        y_filename=sprintf(Y_input,grup,BandN);
        x_filename=sprintf(X_input,grup,BandN);

        x = [];
        y = [];
        
        if exist(x_filename, 'file')
            x=load(x_filename);
            if ~isempty(x)
                x = [x,BandN*ones(length(x(:,1)),1)];
                [~,I]=sort(x(:,5),'descend'); x=x(I(1:Ntop),:);
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
                y = [y,BandN*ones(length(y(:,1)),1)];
                [~,I]=sort(y(:,5),'descend'); y=y(I(1:Ntop),:);
                y_Toplist=cat(1,y_Toplist,y);
            else
                ymissing=cat(1,ymissing,BandN);
            end
        else
            ymissing=cat(1,ymissing,BandN);
        end
    end
    File_cn=sprintf('data_top_%g.mat',grup);
    save (File_cn,'x_Toplist','y_Toplist','xmissing','ymissing')
    
    g_Toplist=[];
    kmax=ceil((f0_grup_max-f0_grup)/f0_band-1);
    k0=0;
    if ~isempty(x_Toplist) && ~isempty(y_Toplist)
        
        x1_y = x_Toplist(:,1) - x_Toplist(:,4) * (x_tref - y_tref);
       
        for k=1:kmax;
            if (k/(f0_rang_grup_band/f0_band)) == fix(k/(f0_rang_grup_band/f0_band)) || k==1;
                I_rang_grup_band_x= find( f0_grup + k0*f0_rang_grup_band - r_W/Tcoh <= x1_y & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= x1_y);
                I_rang_grup_band_y= find( f0_grup + k0*f0_rang_grup_band <= y_Toplist(:,1) & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= y_Toplist(:,1));
                k0=k0+1;
            end
            
            L_band_x= ( f0_grup + k*f0_band - r_W/Tcoh <= x1_y(I_rang_grup_band_x,1)) &  (f0_grup + (k+1)*f0_band +  r_W/Tcoh >= x1_y(I_rang_grup_band_x,1));
            L_band_y= ( f0_grup + k*f0_band            <= y_Toplist(I_rang_grup_band_y,1)) &  (f0_grup + (k+1)*f0_band        >= y_Toplist(I_rang_grup_band_y,1));
            
            if sum(L_band_x) && sum(L_band_y)
                
                x=[x1_y(I_rang_grup_band_x(L_band_x)),x_Toplist(I_rang_grup_band_x(L_band_x),2:end)];  y=y_Toplist(I_rang_grup_band_y(L_band_y),:);               
                
                I_sigchi2_x=1:length(x(:,1));
                I_sigchi2_y=1:length(y(:,1));
                
                chi2_STD_x=chi2_STD(x(:,5),x(:,7),p);
                chi2_STD_y=chi2_STD(y(:,5),y(:,7),p);
                
                I_sigchi2_x(chi2_STD_x>=chisquare_STD_veto | x(:,5)<sigma_veto_toplist(x(:,1)))=[];
                I_sigchi2_y(chi2_STD_y>=chisquare_STD_veto | y(:,5)<sigma_veto_toplist(y(:,1)))=[];
                
                if ~isempty(I_sigchi2_x) && ~isempty(I_sigchi2_y)
                    
                    x=x(I_sigchi2_x,:);  y=y(I_sigchi2_y,:);
                    

                    m_W1 = df0;          
                    m_W4 = sum(df1)/2;
                    
                    [X,Y]=meshgrid(x(:,1),y(:,1));
                    m_D1 = (X-Y)/m_W1;
                    m_W2 = (2./(Y + X))/(Tcoh*1e-4)/PixelFactor;                    
                    
                    [X,Y]=meshgrid(x(:,4),y(:,4));
                    m_D4 = (X-Y)/m_W4;
                    
                    n0=[cos(y(:,3)).*cos(y(:,2)),cos(y(:,3)).*sin(y(:,2)),sin(y(:,3))];
                    n1=[cos(x(:,3)).*cos(x(:,2)),cos(x(:,3)).*sin(x(:,2)),sin(x(:,3))];
                    
                    b=sqrt((n0(:,2)*n1(:,3)'-n0(:,3)*n1(:,2)').^2+(n0(:,3)*n1(:,1)'-n0(:,1)*n1(:,3)').^2+(n0(:,1)*n1(:,2)'-n0(:,2)*n1(:,1)').^2);
                    a=n0(:,1)*n1(:,1)'+n0(:,2)*n1(:,2)'+n0(:,3)*n1(:,3)';
                    
                    m_D2 = atan2(b,a)./m_W2;   
                                        
                    m_D = sqrt(m_D1.^2+m_D2.^2+m_D4.^2);
                    [I_sigchi2_y_cn,I_sigchi2_x_cn]=find(m_D<r_W);
                    
                    if ~isempty(I_sigchi2_x_cn) && ~isempty(I_sigchi2_y_cn)
                        
                        s_harm=2*(x(I_sigchi2_x_cn,5).*y(I_sigchi2_y_cn,5)./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5)));
                        s=(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5))/2;
                        
                        [~,I]=sort(s_harm,'descend');  if length(I)>Ntop; I=I(1:Ntop,:); end;
                        
                        I_sigchi2_x_cn=I_sigchi2_x_cn(I); I_sigchi2_x_cn=reshape(I_sigchi2_x_cn,1, numel(I_sigchi2_x_cn))';
                        I_sigchi2_y_cn=I_sigchi2_y_cn(I); I_sigchi2_y_cn=reshape(I_sigchi2_y_cn,1, numel(I_sigchi2_y_cn))';
                        
                        gf0=(x(I_sigchi2_x_cn,5).*x1_y(I_sigchi2_x_cn,1)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,1))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                        
                        gn=(y(I_sigchi2_y_cn,5)*[1,1,1].*n0(I_sigchi2_y_cn,:)+x(I_sigchi2_x_cn,5)*[1,1,1].*n1(I_sigchi2_x_cn,:))./((x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5))*[1,1,1]);
                        gn=gn./(sqrt(sum(gn'.^2)')*[1,1,1]);
                        
                        galpha  = atan2(gn(:,2),gn(:,1));
                        gdelta  = real(asin(gn(:,3)));
                        
                        gf1=(x(I_sigchi2_x_cn,5).*x(I_sigchi2_x_cn,4)+y(I_sigchi2_y_cn,5).*y(I_sigchi2_y_cn,4))./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5));
                        
                        s_harm=2*(x(I_sigchi2_x_cn,5).*y(I_sigchi2_y_cn,5)./(x(I_sigchi2_x_cn,5)+y(I_sigchi2_y_cn,5)));
                        g_cn_data0 = [gf0,galpha,gdelta,gf1,s(I),s_harm,I_sigchi2_x_cn,I_sigchi2_y_cn];
                        g_Toplist=cat(1,g_Toplist,g_cn_data0);
                        g_cn_data0=[];
                    end
                end
            end
        end       
        [~,I]=sort(g_Toplist(:,1));  g_Toplist=g_Toplist(I,:);        
        File_cn=sprintf('g_Toplist_%g.mat',grup);
        save (File_cn,'g_Toplist')
    end
    if ~isempty(g_Toplist)
        k0=0;
        I_N_Cluster=1;
        I_partial_cluster=[];
        I_partial_cluster{1}=[];
        for k=1:kmax;
            
            if (k/(f0_rang_grup_band/f0_band)) == fix(k/(f0_rang_grup_band/f0_band)) || k==1;
                f0_grup + k0*f0_rang_grup_band
                I_rang_grup_band_x= find( f0_grup + k0*f0_rang_grup_band - r_W/Tcoh <= g_Toplist(:,1) & f0_grup + (k0+1)*f0_rang_grup_band + f0_band >= g_Toplist(:,1));
                x=g_Toplist(I_rang_grup_band_x,:);
                k0=k0+1;
            end
            
            I_band =  find( f0_grup + k*f0_band - r_W_cl/Tcoh <= x(:,1) &  f0_grup + (k+1)*f0_band +  r_W_cl/Tcoh >= x(:,1));
            
            if ~isempty(I_band)
                
                ones_x = ones(1,length(x(I_band,1)));
                ones_y = ones(length(x(I_band,1)),1);
                
                m_W1 = df0;
                m_W2 = ((1./x(I_band))*ones_x + ones_y*(1./x(I_band))')/(Tcoh*1e-4)/PixelFactor/2;
                m_W4 = sum(df1)/2;
                
                m_D1 = ((x(I_band,1)*ones_x) - ones_y*x(I_band,1)')./m_W1;
                m_D4 = ((x(I_band,4)*ones_x) - ones_y*x(I_band,4)')./m_W4;
                
                n0=[cos(x(I_band,3)).*cos(x(I_band,2)),cos(x(I_band,3)).*sin(x(I_band,2)),sin(x(I_band,3))];
                n0n0=n0(:,1)*n0(:,1)'+n0(:,2)*n0(:,2)'+n0(:,3)*n0(:,3)';
                m_D2 = real(acos(n0n0)./m_W2);
                
                m_D= sqrt(m_D1.^2+m_D2.^2+m_D4.^2);
                                
                L_m_D=m_D<r_W_cl;
                I_m = ones_y * (1:1:size(ones_x,2)).*L_m_D;
                i0=[];
                for j=1:length(I_m(:,1))
                    I_v_ji_0=[]; I_v_ji=[];
                    if ~ismember(j,i0)
                        I_v_ji=I_m(j,:);
                        I_v_ji(I_v_ji==0)=[];
                        while length(I_v_ji_0)~=length(I_v_ji);
                            I_v_ji_0=I_v_ji;
                            B=reshape(I_m(I_v_ji,:).' ,1,numel(I_m(I_v_ji,:)));
                            I_v_ji=unique(B);
                            I_v_ji(I_v_ji==0)=[];
                        end
                        i0=cat(2,i0,I_v_ji);
                        I_partial_cluster{I_N_Cluster}=I_rang_grup_band_x(I_band(I_v_ji))';
                        I_N_Cluster=I_N_Cluster+1;
                    end
                end
            end
        end
        
        N=1;
        x=g_Toplist;
        N_Cluster=zeros(length(x(:,1)),1);
        I_Ci=[];
        A=I_partial_cluster{1};
        if I_N_Cluster>2;
            for i=1:(I_N_Cluster-1);
                if isempty(I_Ci);
                    N_Cluster(A)=N;
                    N=N+1;
                else
                    N_Cluster(A)=N_Cluster_i(I_Ci(1));
                    B = ismember(N_Cluster,N_Cluster_i(I_Ci));
                    N_Cluster(B)=N_Cluster_i(I_Ci(1));
                end
                if i<I_N_Cluster-1;
                    A=I_partial_cluster{i+1};
                    N_Cluster_i=unique(N_Cluster(A));
                    I_Ci = find(N_Cluster_i~=0);
                end
            end
        else
            N_Cluster(A)=1;
        end
        
        j=1;
        for j=1:max(N_Cluster);
            
            I_Cluster_i=find(N_Cluster==j);
            s=x(I_Cluster_i,6);
            
            Cluster.sign_mean{j} =mean(x(I_Cluster_i,5));
            Cluster.ind{j}       =I_Cluster_i;
            Cluster.s_harm{j}     =mean(s);
            Cluster.sign_sum{j}  =sum(s);
            Cluster.deg_x{j}     =length(unique(x(I_Cluster_i,7)));
            Cluster.deg_y{j}     =length(unique(x(I_Cluster_i,8)));
            Cluster.f0{j}        =sum(x(I_Cluster_i,1).*s)/Cluster.sign_sum{j};
            
            n0=sum(s*[1,1,1].*[cos(x(I_Cluster_i,3)).*cos(x(I_Cluster_i,2)),cos(x(I_Cluster_i,3)).*sin(x(I_Cluster_i,2)),sin(x(I_Cluster_i,3))],1);
            n0=n0./sqrt(sum(n0.^2));
            
            Cluster.alpha{j}     =real(atan2(n0(:,2),n0(:,1)));
            Cluster.delta{j}     =real(asin(n0(:,3)));
            
            Cluster.f1{j}        =sum(x(I_Cluster_i,4).*s)/Cluster.sign_sum{j};
            Cluster.length{j}    =length(I_Cluster_i);
            
            Cluster.sign_max{j}  =max(s);
            
            if length(I_Cluster_i)<2; Cluster.noise{j}=1; else Cluster.noise{j}=0; end
        end
        follow_up=[[Cluster.f0{:}]',[Cluster.alpha{:}]',[Cluster.delta{:}]',[Cluster.f1{:}]',[Cluster.sign_mean{:}]',[Cluster.s_harm{:}]',...
            [Cluster.sign_sum{:}]',[Cluster.length{:}]',[Cluster.deg_x{:}]',[Cluster.deg_y{:}]'];           
    end  
    
    if ~isempty(follow_up)
        follow_up=follow_up(follow_up(:,7)>0,:);
        
        [~,I]=sort(follow_up(:,1)); follow_up=follow_up(I,:);
        if ~isempty(follow_up)
            I_band_cumsum = cumsum([true; sum(abs(diff(follow_up(:,1),1,1)),2) > f0_band]);
            I_band=[];
            for i=1:max(I_band_cumsum);
                I_1=find(I_band_cumsum==i);
                [~,I]=max(follow_up(I_1,6)); %%%%%%%%%%%%
                I_band=cat(1,I_band,I_1(I));
            end
            follow_up=follow_up(I_band,:);
        end
        follow_up=follow_up(follow_up(:,5)>0,:);
        param=[sig,chisquare_STD_veto,r_W,r_W_cl,0,0];
    end
    
    File_cn=sprintf('follow_up_%g.mat',grup);
    save (File_cn,'follow_up','Cluster','param')
end



%                    [X,Y]=meshgrid(x(:,5),y(:,5));                    
%                     m_W5 = max(X,Y);
%                     m_D5 = abs((y(:,5)*ones_x) - ones_y*x(:,5)')./m_W5;                                         
%                     [I_sigchi2_y_cn,I_sigchi2_x_cn]=find(m_D<r_W & m_D5<0.4);