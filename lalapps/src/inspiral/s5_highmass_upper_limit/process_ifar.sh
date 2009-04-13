LALAPPS_NEWCORSE=${1}
cat=${2}

if ${LALAPPS_NEWCORSE} -i H1,H2,L1 -b 0,50,85,inf -p thinca --veto-segments vetoes_${cat}.xml.gz --veto-segments-name=vetoes -v *${cat}.sqlite; then echo $(date) " : Done with corse on " ${cat} > successes.txt; else echo $(date) " : Failed corse on " ${cat} > errors.txt; exit; fi;

