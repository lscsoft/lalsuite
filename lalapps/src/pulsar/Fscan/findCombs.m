function findCombs(threshold,outputname,nround,varargin)
% findCombs(threshold,outputname,nround,varargin)
% If threshold < 1 then 1-threshold precent of the loudest powers are searched for combs; e.g., threshold=.999 using the top .0001 power values
% If threshold >=1 then po > threshold are serached for combs.
% The outputname get _combs.txt appended to it.
% Give nround equal to, e.g., 2, rounds to frequencies to 2 decimal place.
% Run findCombs(threshold,outputname,filename) if filename is an input file and [f,p] = textread(filename,'%f %f') is used
% Run findCombs(threshold,outputname,f,p) if frequency and power are given as inputs.

if nargin == 4,
    filename = varargin{1};
    [f,p] = textread(filename,'%f %f');
elseif nargin==5,
    f = varargin{1};
    p = varargin{2};
end

p_ = sort(p);
if threshold<1,
    ind999 = floor(threshold*numel(p));
    thresh = p_(ind999);
else
    thresh = threshold;
end
nround = round(nround);
nroundfactor = round(10^nround);
loudind = find(p>thresh);
loudf = f(loudind);
loudf = round(loudf(:)*nroundfactor)/nroundfactor;

% Open the output file:
outputCombTextFile = sprintf('%s_combs.txt',outputname);
fid = fopen(outputCombTextFile,'w');
fprintf(fid,'Combs given as delta f: teeth frequencies ... \n\n');

%Now we create a matrix of time-averaged differences. The row,col of each
%entry directly corresponds to the frequencies that created that
%difference. Almost directly from Greg's suggestion.

afreq = zeros(numel(loudf)); 
for i = 1:numel(loudf),
    for j = 1:numel(loudf),
        afreq(i,j) = round(1800*(loudf(j)-loudf(i))); 
    end
end

%Below, I replace the negative values with '0'. Didn't use abs so I could 
%avoid double counting. Additionally, since a comb is defined by something 
%that appears at least three times, we can get rid of most of the
%differences in the above matrix.  

afreq(afreq(:)<0) = 0; 
udiffs = nonzeros(unique(afreq)); 
udiff_max = max(udiffs)/3; 
udiffs(udiffs(:)>udiff_max) = [];
udiffs = sort(udiffs);

master = []; %Frequencies already due to one comb are added to this list.
freq_diffs = udiffs/1800; 
for q = 1:numel(udiffs),
    [row,col] = find(afreq == udiffs(q)); 
    freqs = [loudf(row) loudf(col)]; 
    [row_,col_] = find(freqs == 0);
    freqs(row_,:) = [];
    freqs = unique(freqs(:));
    %out1 = [freq_diffs(q)];
    deltaF = [freq_diffs(q)];
    out1 = [];
    %This next loop determines whether the frequency difference shows up 3
    %times in a row or not.  If not, we ignore it. If so, it's added to the
    %output list.
    for k = 1:numel(freqs),
        if and(ismember(freqs(k)+freq_diffs(q),freqs(:)),ismember(freqs(k)+2*freq_diffs(q),freqs(:))),
            out1 = [out1 freqs(k) freqs(k)+freq_diffs(q) freqs(k)+2*freq_diffs(q)];
        end
    end
    OUT = nonzeros(sort(unique(out1)));
    if numel(OUT)>=3,
        fprintf(fid,'%f: ', deltaF);
        fprintf(fid,'%f ', OUT);
        fprintf(fid,'\n\n\n');
        master = [master deltaF'];
        master = [master OUT'];
        loudf(ismember(loudf(:),master)) = 0; %Frequencies that have been 
        %printed are removed, ie cannot be a harmonic of some other frequency.
    end        
end

fclose(fid);
return;