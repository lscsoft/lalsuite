#!/bin/bash

echo dot ${1}.dot update > ${1}.tmp

# datafind 
echo "L datafind"
jid=1; for allid in `grep datafind.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep L | awk '{print $2}'` ; do sedcmd="s/$replid/datafind_L_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H datafind"
jid=1; for allid in `grep datafind.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H | awk '{print $2}'` ; do sedcmd="s/$replid/datafind_H_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G datafind"
jid=1; for allid in `grep datafind.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G | awk '{print $2}'` ; do sedcmd="s/$replid/datafind_H_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# tmpltbank

echo "L1 tmpltbank"
jid=1; for allid in `grep tmpltbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep L1 | awk '{print $2}'` ; do sedcmd="s/$replid/tmplt_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 tmpltbank"
jid=1; for allid in `grep tmpltbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/tmplt_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 tmpltbank"
jid=1; for allid in `grep tmpltbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/tmplt_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 tmpltbank"
jid=1; for allid in `grep tmpltbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/tmplt_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


# inspiral

echo "L1 inspiral"
jid=1; for allid in `grep inspiral_L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep L1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 inspiral"
jid=1; for allid in `grep inspiral_H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 inspiral"
jid=1; for allid in `grep inspiral_H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 inspiral"
jid=1; for allid in `grep inspiral_G1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# trigbank

echo "L1 trigbank"
jid=1; for allid in `grep trigbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="L1"'  | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


echo "H1 trigbank"
jid=1; for allid in `grep trigbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="H1"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 trigbank"
jid=1; for allid in `grep trigbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="H2"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 trigbank"
jid=1; for allid in `grep trigbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="G1"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


# inspiral veto 
echo "L1 inspiral veto"
jid=1; for allid in `grep inspiral_L1_veto.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep L1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_L1_V_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 inspiral veto"
jid=1; for allid in `grep inspiral_H1_veto.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_H1_V_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 inspiral veto"
jid=1; for allid in `grep inspiral_H2_veto.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_H2_V_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 inspiral veto"
jid=1; for allid in `grep inspiral_G1_veto.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/insp_G1_V_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# sinca 
echo "L1 sinca"
jid=1; for allid in `grep s_inca.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep L1 | awk '{print $2}'` ; do sedcmd="s/$replid/sinca_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 sinca"
jid=1; for allid in `grep s_inca.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/sinca_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 sinca"
jid=1; for allid in `grep s_inca.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/sinca_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 sinca"
jid=1; for allid in `grep s_inca.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/sinca_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca - 2 ifo
echo "G1H1 thinca"
jid=1; for allid in `grep thinca_G1H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 thinca"
jid=1; for allid in `grep thinca_G1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 thinca"
jid=1; for allid in `grep thinca_G1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 thinca"
jid=1; for allid in `grep thinca_H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 thinca"
jid=1; for allid in `grep thinca_H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 thinca"
jid=1; for allid in `grep thinca_H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca - 3 ifo
echo "G1H1H2 thinca"
jid=1; for allid in `grep thinca_G1H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 thinca"
jid=1; for allid in `grep thinca_G1H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 thinca"
jid=1; for allid in `grep thinca_G1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 thinca"
jid=1; for allid in `grep thinca_H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca - 4 ifo 
echo "G1H1H2L1 thinca"
jid=1; for allid in `grep thinca_G1H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done



# thinca_slides - 2 ifo
echo "G1H1 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 thinca_slides"
jid=1; for allid in `grep thinca_slides_H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca_slides - 3 ifo
echo "G1H1H2 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca_slides - 4 ifo 
echo "G1H1H2L1 thinca_slides"
jid=1; for allid in `grep thinca_slides_G1H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca_slides_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


# thinca2 - 2 ifo
echo "G1H1 thinca2"
jid=1; for allid in `grep thinca2_G1H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 thinca2"
jid=1; for allid in `grep thinca2_G1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 thinca2"
jid=1; for allid in `grep thinca2_G1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 thinca2"
jid=1; for allid in `grep thinca2_H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 thinca2"
jid=1; for allid in `grep thinca2_H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 thinca2"
jid=1; for allid in `grep thinca2_H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca2 - 3 ifo
echo "G1H1H2 thinca2"
jid=1; for allid in `grep thinca2_G1H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 thinca2"
jid=1; for allid in `grep thinca2_G1H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 thinca2"
jid=1; for allid in `grep thinca2_G1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 thinca2"
jid=1; for allid in `grep thinca2_H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca2 - 4 ifo 
echo "G1H1H2L1 thinca2"
jid=1; for allid in `grep thinca2_G1H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca2_slides - 2 ifo
echo "G1H1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H2 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca2_slides - 3 ifo
echo "G1H1H2 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H1H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H1L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep H1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# thinca2_slides - 4 ifo 
echo "G1H1H2L1 thinca2_slides"
jid=1; for allid in `grep thinca2_slides_G1H1H2L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep G1 | awk '{print $2}'` ; do sedcmd="s/$replid/thinca2_slides_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# sire
echo "sire"
jid=1; for allid in `grep sire.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/sire_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# sire_clust
echo "sire cluster"
jid=1; for allid in `grep sire_clust.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/sire_clust_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# cohbank 
echo "G1H1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H2-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1H2-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H2L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1H2 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1H2-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H2L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1H2L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1H2L1 cohbank"
jid=1; for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1H2L1-' | awk '{print $2}'` ; do sedcmd="s/$replid/cohbank_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} >> ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


# trigbank coherent

echo "L1 trigbank"
jid=1; for allid in `grep trigbank_coherent_L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="L1"'  | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_coherent_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 trigbank_coherent"
jid=1; for allid in `grep trigbank_coherent_H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="H1"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_coherent_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 trigbank_coherent"
jid=1; for allid in `grep trigbank_coherent_H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="H2"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_coherent_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 trigbank_coherent"
jid=1; for allid in `grep trigbank_coherent_G1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="G1"' | awk '{print $2}'` ; do sedcmd="s/$replid/trigbank_coherent_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# inspiral coherent

echo "L1 insp_coh"
jid=1; for allid in `grep insp_coh_L1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'L1'  | awk '{print $2}'` ; do sedcmd="s/$replid/insp_coh_L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1 insp_coh"
jid=1; for allid in `grep insp_coh_H1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'H1' | awk '{print $2}'` ; do sedcmd="s/$replid/insp_coh_H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2 insp_coh"
jid=1; for allid in `grep insp_coh_H2.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'H2' | awk '{print $2}'` ; do sedcmd="s/$replid/insp_coh_H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1 insp_coh"
jid=1; for allid in `grep insp_coh_G1.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep 'G1' | awk '{print $2}'` ; do sedcmd="s/$replid/insp_coh_G1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# coherent inspiral
echo "G1H1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H2-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1H2-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H2L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H2L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1H2 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1H2-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H1H2_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H1L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H2L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H2L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "H1H2L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_H1H2L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "G1H1H2L1 coherent inspiral"
jid=1; for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep '_G1H1H2L1-'  | awk '{print $2}'` ; do sedcmd="s/$replid/cohinspiral_G1H1H2L1_$jid/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

