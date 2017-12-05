clear
path(path,'/home/siemens/ligotools/matlab')

N=300;                %no of slides

cd /scratch2/xavi/CosmicStrings0/slides/

% This loop finds the minimum and maximum amplitudes
A1=[];
A2=[];
for i=1:N
   
    i
    file_neg=['H1_M',num2str(i),'.xml'];
    file_pos=['H1_P',num2str(i),'.xml'];

    H1ctrig_neg = readMeta(char(file_neg),'sngl_burst');
    H1ctrig_pos = readMeta(char(file_pos),'sngl_burst');
    
    A1=[A1;abs(H1ctrig_neg.amplitude)];
    A2=[A2;abs(H1ctrig_pos.amplitude)];
    
end

Amin=min(min(A1),min(A2))
Amax=max(max(A1),max(A2))

dA=0.5
A=Amin:dA:1.5*Amax;

% This loops finds the histograms
histA1=[];
histA2=[];
for i=1:N
       
    file_neg=['H1_M',num2str(i),'.xml'];
    file_pos=['H1_P',num2str(i),'.xml'];

    H1ctrig_neg = readMeta(char(file_neg),'sngl_burst');
    H1ctrig_pos = readMeta(char(file_pos),'sngl_burst');
    
    histA1=[histA1,histc(abs(H1ctrig_neg.amplitude),A)];
    histA2=[histA2,histc(abs(H1ctrig_pos.amplitude),A)];
    
end

histA=[histA1,histA2];

meanhistA=mean(histA,2);
stdhistA=std(histA,0,2); %the 0 is a flag to divide by N-1 

% Now make the plot
ymin=0.01;

figure
H1ctrig_foreground=readMeta('H1-BURCA_H1L1_P_0_--1-0.xml','sngl_burst');
histA_foreground=histc(abs(H1ctrig_foreground.amplitude),A);
%foreground
semilogy(A+0.5*dA,histA_foreground,'k.','MarkerSize',25)
%-- Add shaded bars for background spread
hold on
grid on

bkggray = 0.78;
cpatch = [bkggray bkggray bkggray];
for i=1:size(A,2)
  if meanhistA(i) > 0.0
    x1 = A(i);
    x2 = A(i)+dA;
    yu = meanhistA(i) + stdhistA(i);
    yl = max(meanhistA(i) - stdhistA(i),ymin);
    p = patch( [x1 x2 x2 x1], [yl yl yu yu], cpatch, 'LineStyle','none' );
  end
end

% -- add stairs
hold on
stairs(A,meanhistA,'k')

% -- errorbar for mean
for i=1:size(A,2)
 if meanhistA(i) > 0.0
    x = A(i)+0.5*dA;
    yu = meanhistA(i) + stdhistA(i)/sqrt(2*N);
    yl = max(meanhistA(i) - stdhistA(i)/sqrt(2*N),ymin);
    line([x x],[yl yu],'LineWidth',3.0,'Color',[0,0,0]);
 end
end

%foreground again
plot(A+0.5*dA,histA_foreground,'k.','MarkerSize',25)
ylabel('Events in S4','FontSize',14)
xlabel('Amplitude A/10^{-20} [s^{-1/3}]','FontSize',14)


%-- Add legend
legendfont = 12;
xmax = A(size(A,2));
ymax = 10;
legymax = 9;
legymin = 2.25;
bkggray = 0.78;
% There are also a few things hard-coded in the script below

%-- Histogram parameters
binwidth = 1;

xmark = 0.44*xmax;
xml = xmark - 0.5*binwidth;
xmu = xmark + 0.5*binwidth;
xtext = xmark + 0.1*xmax;
logspan = log(legymax)-log(legymin);

patch( [xmark-binwidth 0.98*xmax 0.98*xmax xmark-binwidth], ...
       [legymin legymin legymax legymax], [1 1 1] );

y = exp(log(legymin)+(5.0/6.0)*logspan);
z = line( xmark, y, 'Color',[0 0 0], 'LineStyle','none', ...
          'Marker','.', 'MarkerSize',25 );
text(xtext,y, 'Zero-lag event count in bin', 'FontSize',legendfont );

y = exp(log(legymin)+(3.0/6.0)*logspan);
z = line([xml,xmu],[y y],'Color',[0 0 0],'LineWidth',3.0);
ybox = [ exp(log(y)-logspan/12) exp(log(y)+logspan/12) ];
z = line([xmark,xmark],ybox,'Color',[0 0 0],'LineWidth',3.0);
text(xtext,y,'Estimated mean background','FontSize',legendfont);

y = exp(log(legymin)+(1.0/6.0)*logspan);
del = logspan/10.0;
ybox = [ exp(log(y)-del) exp(log(y)-del) exp(log(y)+del) exp(log(y)+del) ];
patch( [xml xmu xmu xml], ybox, cpatch, 'LineStyle','none' );
text(xtext,y,'Background spread (rms)','FontSize',legendfont);

line([0 8 8 0 0],[ymin ymin ymax ymax ymin], ...
     'LineWidth',1.5,'Color',[0 0 0]);

