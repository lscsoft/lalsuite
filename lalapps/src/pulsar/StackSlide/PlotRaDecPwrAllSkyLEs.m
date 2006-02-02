function PlotRaDecPwrAllSkyLEs(RA,DEC,PWR,decimation,parea,type)
% Usage: PlotRaDecPwrAllSkyLEs(RA,DEC,PWR,decimation,parea,type)
% RA,DEC,PWR: vectors with RA, DEC, and PWR from a search 
% decimation: factor to reduce amount of data
% (e.g., 1 = no reduction; 10 = reduce by factor of 10).
% parea: point area (in pts squared) of dots used in plot.
% type: type of plot (1 = scatter plot, 2 = use projected ra's.

vlength = length(RA);
ra = RA(1:decimation:vlength);
dec = DEC(1:decimation:vlength);
pwr = PWR(1:decimation:vlength);
if (type==1) 
   scatter(ra,dec,parea,pwr,'filled');
else
   projra = pi*(1 - cos(dec)) + ra.*cos(dec);
   scatter(projra,dec,parea,pwr,'filled');
end
xlabel('RA (radians)');
ylabel('DEC (radians)');
title('Power For Each Sky Position');
return;
