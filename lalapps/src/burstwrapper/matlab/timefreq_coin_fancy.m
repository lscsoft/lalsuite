function [Iref, Itest] = time_coin(ref_time, ref_dur, test_time, test_dur, ref_pfreq, test_pfreq, ref_bw, test_bw);

% [Iref, Itest] = time_coin(ref_time, ref_dur, test_time, test_dur)
% returns indices of reference triggers (Iref) which are overlaping in time with
% test triggers and |ref_pfreq - test_pfreq| <= 2*(ref_bw+test_bw)
% input must be time-ordered according to the start time of the intervals
% reference intervals must be non-overlapping

Iref=[];
Itest=[];


i1 = 1;

for j=1:length(ref_time)
   
   isin = 1;
   
   while(test_time(i1)+test_dur(i1) < ref_time(j))
      i1=i1+1;
      if(i1 >= length(test_time))
         i1 = length(test_time);
         break;
      end;
   end;
   
   for i2=i1:1:length(test_time)

     fmatch = 0;
     if(abs(ref_pfreq(j) - test_pfreq(i2)) <= 2.0*(ref_bw(j)+test_bw(i2)))
      fmatch = 1;
     elseif(ref_pfreq(j)==0 | test_pfreq(i2)==0)
      warning('Allowing frequency match because frequency = 0!');
      fmatch=1;
     end

     if(fmatch)
       if(ref_time(j) <= test_time(i2) + test_dur(i2) & ref_time(j)+ref_dur(j) >= test_time(i2))
         if(isin)
            Iref=[Iref;j];
            isin=0;
         end;
         
         Itest=[Itest;i2];
       end;
      
      end;

      if(test_time(i2) > ref_time(j)+ref_dur(j))
         break;
      end;
      
   end;
   
end;

Itest = unique(Itest);
         
      
      
