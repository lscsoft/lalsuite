F=[];
T=[];

grup_min=50;
grup_Band=50;
grup_max=800;

f0_band=0.1;

for grup=grup_min:grup_Band:grup_max;
    
	f0_grup=grup;

    follow_up_file=sprintf('../O1/follow_up_%g.mat',grup);
    data_top_file=sprintf('../O1/data_top_%g.mat',grup);
    g_top_file=sprintf('../O1/g_Toplist_%g.mat',grup);
    
    load(follow_up_file)
    load(data_top_file)
    load(g_top_file)
    
if f0_grup>=50;
A1=0.4902
A2=1.414
B1=0.3581
B2=1.484
chisquare_STD_veto = 5;
end

if f0_grup>=100;
A1=0.2168
A2=1.428
B1=0.1902
B2=1.499
chisquare_STD_veto = 6;
end

if f0_grup>=200;
A1 = 0.1187;
A2 = 1.470;
B1 = 0.06784;
B2 = 1.697;
chisquare_STD_veto = 5.3;
end

        chi2_STD = @(x,y,p) (y-(p-1)-A1*x.^A2)./(sqrt(2*p-2)+B1*x.^B2);
    
    kmax=ceil(grup_Band/f0_band-1);
    for k=1:500;
        data={'Name' 'f0' 'alpha' 'delta' 'f1' 'significance' 'chi2'};
        x=1;y=1;
        for iter=1:4;
            if ~isempty(x) && ~isempty(y)
                
		if iter==4;
                    I= ( f0_grup + k*f0_band <= follow_up(:,1)) & (f0_grup + (k+1)*f0_band >= follow_up(:,1));
                    a=follow_up(I,:);
                    [xmax,Ix]=max(x(:,5));
                    [ymax,Iy]=max(y(:,5));
                    name = 'F';
                    if ~isempty(a)
                        A=[a(1:5),0];
                        F=[F;a,percentaje_chi2x,percentaje_chi2y];
                        T=[T;x(Ix,1),xmax,x(Ix,2),x(Ix,3),min(x(:,5)),y(Iy,1),ymax,x(Iy,2),x(Iy,3),min(x(:,5)),percentaje_chi2x,percentaje_chi2y];
                    else
                        A=[0,0,0,0,0,0];
                        T=[T;x(Ix,1),min(x(:,5)),xmax,x(Ix,2),x(Ix,3),percentaje_chi2x,y(Iy,1),min(x(:,5)),ymax,x(Iy,2),x(Iy,3),percentaje_chi2y];
                    end
                end
                if iter==3;
                I= ( f0_grup + k*f0_band <= g_Toplist(:,1)) & (f0_grup + (k+1)*f0_band >= g_Toplist(:,1));
                    name = 'G';
                    if ~isempty(I)
                        A=[g_Toplist(I,1:5),zeros(length(g_Toplist(I,1)),1)];
                    else
                        A=[0,0,0,0,0,0];
                    end
                end
                if iter==2;
                    name = 'H1';
                    x=x_Toplist(x_Toplist(:,end)==k,:);
                    Ch1x=chi2_STD(x(:,5),x(:,7),16);
                    percentaje_chi2x=(sum(Ch1x>6))/(length(Ch1x))*100;
                    A=[x(:,1:5),x(:,7)];
                end
                if iter==1;
                    name = 'L1';
                    y=y_Toplist(y_Toplist(:,end)==k,:);
                    Ch1y=chi2_STD(y(:,5),y(:,7),16);
                    percentaje_chi2y=(sum(Ch1y>6))/(length(Ch1y))*100;
                    A=[y(:,1:5),y(:,7)];
                end                
                data=cat(1,data,[repmat({name},[length(A(:,1)) 1]),num2cell(A)]);
            end
	end	            
            name=sprintf('data/f0_%.1f.tsv',f0_grup + k*f0_band);
            fid = fopen(name,'wt');
            if fid>0
                for j=1:size(data,1)
                    if j==1;
                        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',data{j,:});
                    else
                        fprintf(fid,'%s\t%f\t%f\t%f\t%e\t%f\t%f\n',data{j,:});
                    end
                end
                fclose(fid);
            end
    end
    grup    
end

data={'id' 'f0' 'alpha' 'delta' 'f1' 'sig_mean' 's_old' 'sig_sum' 'length' 'xdeg' 'ydeg' 'xchi2' 'ychi2'};
if ~isempty(F)
data=cat(1,data,num2cell([[1:1:sum(F(:,1)>0)]',F(F(:,1)>0,:)]));

name=sprintf('data/follow_up.csv');
fig = fopen(name,'wt');
if fid>0
    
    for k=1:size(data,1)
        if k==1;
            fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',data{k,:});
        else
            fprintf(fid,'%.0f,%.8f,%.8f,%.8f,%.8e,%.2f,%.2f,%.2f,%.f,%.f,%.f,%.1f,%.1f\n',data{k,:});
        end
    end
    
    fclose(fid);
end

data={'xf0' 'xmin' 'xmax' 'xalpha' 'xdelta' 'xchi2' 'yf0' 'ymin' 'ymax' 'yalpha' 'ydelta' 'ychi2'};
data=cat(1,data,num2cell(T));

name=sprintf('data/Total.csv');
fid = fopen(name,'wt');
if fid>0
    
    for k=1:size(data,1)
        if k==1;
            fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',data{k,:});
        else
            fprintf(fid,'%.8f,%.3f,%.3f,%.8f,%.8f,%.3f,%.8f,%.3f,%.3f,%.8f,%.8f,%.3f\n',data{k,:});
        end
    end
    
    fclose(fid);
end
end
