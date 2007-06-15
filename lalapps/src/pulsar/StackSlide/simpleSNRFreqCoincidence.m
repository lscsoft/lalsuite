function [f0H1out,FREQH1out,PWRH1out,f0L1out,FREQL1out,PWRL1out,iOutlier] = simpleSNRFreqCoincidence()

[f0H1,RAH1,DECH1,FREQH1,FDOTH1,SH1,SBinsH1,PWRH1,ULH1,ULUncH1,ConfH1,UncConfH1] = textread('resultsStackSlideS4H1WithVetoes50to1000Hz.txt','%f %f %f %f %f %f %f %f %f %f %f %f');
[f0L1,RAL1,DECL1,FREQL1,FDOTL1,SL1,SBinsL1,PWRL1,ULL1,ULUncL1,ConfL1,UncConfL1] = textread('resultsStackSlideS4L1WithVetoes50to1000Hz.txt','%f %f %f %f %f %f %f %f %f %f %f %f');

SNRH1 = (PWRH1 - 1.0)*sqrt(1004);
SNRL1 = (PWRL1 - 1.0)*sqrt(899);
          
fprintf('m,FREQH1(j),RAH1(j),DECH1(j),FDOTH1(j),SBinsH1(j),SNRH1(j),ULH1(j),FREQL1(j),RAL1(j),DECL1(j),FDOTL1(j),SBinsL1(j),SNRL1(j),ULL1(j)\n');

k = 0;
m = 0;
for j = 1:length(f0H1)
    if ( ( (SNRH1(j) <= 7.0) && (SNRL1(j) <= 7.0) ) ) 

         k = k + 1;
         f0H1out(k) = f0H1(j);
         FREQH1out(k) = FREQH1(j);
         PWRH1out(k) = PWRH1(j);
         f0L1out(k) = f0L1(j);
         FREQL1out(k) = FREQL1(j);
         PWRL1out(k) = PWRL1(j);

    else

       if ( ( (SNRH1(j) > 7.0) && (SNRL1(j) > 7.0) ) )
         %if ( (abs(FREQH1(j) - FREQL1(j)) <= 1.1e-4*f0H1(j)) )
         if ( (abs(FREQH1(j) - FREQL1(j)) <= 2.2e-4*FREQH1(j)) )
          k = k + 1;
          m = m + 1;
          iOutlier(m) = k;
          f0H1out(k) = f0H1(j);
          FREQH1out(k) = FREQH1(j);
          PWRH1out(k) = PWRH1(j);
          f0L1out(k) = f0L1(j);
          FREQL1out(k) = FREQL1(j);
          PWRL1out(k) = PWRL1(j);
          fprintf('%3i %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e\n',m,FREQH1(j),RAH1(j),DECH1(j),FDOTH1(j),SBinsH1(j),SNRH1(j),ULH1(j),FREQL1(j),RAL1(j),DECL1(j),FDOTL1(j),SBinsL1(j),SNRL1(j),ULL1(j));
         end
       end

       if (j > 1)        
          if ( ( (SNRH1(j) > 7.0) && (SNRL1(j-1) > 7.0) ) )
            %if ( (abs(FREQH1(j) - FREQL1(j-1)) <= 1.1e-4*f0H1(j)) )
            if ( (abs(FREQH1(j) - FREQL1(j-1)) <= 2.2e-4*FREQH1(j)) )
             k = k + 1;
             m = m + 1;
             iOutlier(m) = k;
             f0H1out(k) = f0H1(j);
             FREQH1out(k) = FREQH1(j);
             PWRH1out(k) = PWRH1(j);
             f0L1out(k) = f0L1(j-1);
             FREQL1out(k) = FREQL1(j-1);
             PWRL1out(k) = PWRL1(j-1);
             fprintf('%3i %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e\n',m,FREQH1(j),RAH1(j),DECH1(j),FDOTH1(j),SBinsH1(j),SNRH1(j),ULH1(j),FREQL1(j-1),RAL1(j-1),DECL1(j-1),FDOTL1(j-1),SBinsL1(j-1),SNRL1(j-1),ULL1(j-1));
            end
          end
       end

       if (j < length(f0H1))
          if ( ( (SNRH1(j) > 7.0) && (SNRL1(j+1) > 7.0) ) )
            %if ( (abs(FREQH1(j) - FREQL1(j+1)) <= 1.1e-4*f0H1(j)) )
            if ( (abs(FREQH1(j) - FREQL1(j+1)) <= 2.2e-4*FREQH1(j)) )
             k = k + 1;
             m = m + 1;
             iOutlier(m) = k;
             f0H1out(k) = f0H1(j);
             FREQH1out(k) = FREQH1(j);
             PWRH1out(k) = PWRH1(j);
             f0L1out(k) = f0L1(j+1);
             FREQL1out(k) = FREQL1(j+1);
             PWRL1out(k) = PWRL1(j+1);
             fprintf('%3i %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e %10.6f %8.6f %9.6f %10.3e %3i %7.2f %9.2e\n',m,FREQH1(j),RAH1(j),DECH1(j),FDOTH1(j),SBinsH1(j),SNRH1(j),ULH1(j),FREQL1(j+1),RAL1(j+1),DECL1(j+1),FDOTL1(j+1),SBinsL1(j+1),SNRL1(j+1),ULL1(j+1));
            end
          end
       end

    end
 end
 
 return
