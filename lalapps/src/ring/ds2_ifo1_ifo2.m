
function[ds2_out] = ds2_ifo1_ifo2(trig_type, ifo1, ifo2, file1, file2)

%
% ds2_ifo1_ifo2 calculates ds2 between 2 specified ifos, given the *.mat
% files output from separate.m
%
% It implements the time minimisation function time_min.m, based on the
% lalapps function written by Neil Cornish, and can be used to check the
% ds2 values for any triple or double coincidences for specified ifos.
%
% ds2_ifo1_ifo2 calls metric3d.m and time_min.m as external functions.
%
% trig_type = 'bg' or 'inj'
%              describes type of run: background (timeslides) or injections
%
% ifo1, ifo2 = 'H1', 'H2' or 'L1'
%
% file1, file2 = e.g. injH1_trip.mat, bgH2inL1_doub.mat etc.
%
% Usage Example:
%
% trig_type = 'bg';
%
% ifo1 = 'H1';
% ifo2 = 'L1';
%
% file1 = 'injH1inL1_doub.mat';
%
% file2 = 'injL1inH1_doub.mat';
%
% ds2H1L1_doub = ds2_ifo1_ifo(trig_type, ifo1, ifo2, file1, file2)
%
% *.mat files, time_min.m and metric3d.m should be in the working directory or
% accesssible path.
%
% Author: Fiona Speirits

% Load coincidence files

coincs1 = load(file1);
coincs2 = load(file2);

% Pre-allocate memory

ds2_out = zeros(length(coincs1.t));
dt_new = zeros(length(coincs1.t));
dt_min = zeros(length(coincs1.t));
dt_max = zeros(length(coincs1.t));
ds2_min = zeros(length(coincs1.t));
dt_best = zeros(length(coincs1.t));
ds2 = zeros(length(coincs1.t));

light_time = 1.001285e-2;

% Check files contain triggers, to add to unhelpful error messages from MatLab

if isempty(coincs1.t)
    
    disp('Error: Files do not contain any triggers')

end

% Check if doubles, in which case reslide the data

if strcmp(trig_type,'bg')
    dt = mod(coincs2.t - coincs1.t, 1 );
    dt = (dt - (dt > .5));
else
    dt = coincs2.t - coincs1.t;
end

% Check if H1H2, in which case time_min.m does not need to be called

if (strcmp(ifo1,'H1') && strcmp(ifo2,'H2'))  || (strcmp(ifo1,'H2') && strcmp(ifo2,'H1'))

    for i = 1:length(coincs1.t)

        ds2_out(i) = metric3d(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i),dt(i));
        dt_new(i) = dt(i);
    
    end

else

% Loop through each of the pairs, calling the functions metric3d and time_min.
for i = 1:length(coincs1.t)
    dt_min(i) = dt(i) - 1.00*light_time;
    dt_max(i) = dt(i) + 1.00*light_time;
    ds2_min(i) = metric3d(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i),dt(i));

    % Otherwise, solve for dt that minimises ds2

    dt_best(i) = time_min(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i));

    % Check that this is a valid time
    if dt_best(i) < dt_max(i) && dt_best(i) > dt_min(i)
        ds2_min(i) = metric3d(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i),dt_best(i));
        dt_new(i) = dt_best(i);
    else
        ds2_min(i) = metric3d(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i),dt_min(i));
        ds2(i) = metric3d(coincs1.f(i),coincs2.f(i),coincs1.q(i),coincs2.q(i),dt_max(i));
        if ds2(i) < ds2_min(i);
            ds2_min(i) = ds2(i);
            dt_new(i) = dt_max(i);
        else
            dt_new(i) = dt_min(i);
        end
    end

    ds2_out(i) = ds2_min(i);

end

end

end