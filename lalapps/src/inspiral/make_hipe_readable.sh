#!/bin/bash

echo dot ${1}.dot update > ${1}.tmp

# datafind 
for observatory in {G,H,L,V}; do 
  echo "${observatory} datafind"
  jid=1
  for allid in `grep datafind.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep \"${observatory}\" | awk '{print $2}'` ; do 
      sedcmd="s/$replid/datafind_${observatory}_${jid}/g"
      jid=$((jid+1))
      sed $sedcmd ${1} >> ${1}.tmp
      mv ${1}.tmp ${1} 
    done 
  done
done

for ifo in {G1,H1,H2,L1,V1}; do 
  # tmpltbank
  echo "${ifo} tmpltbank"
  jid=1; 
  for allid in `grep tmpltbank_${ifo}.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep ${ifo} | awk '{print $2}'` ; do
      sedcmd="s/$replid/tmplt_${ifo}_${jid}/g" 
      jid=$((jid+1))
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1}
    done
  done

  # inspiral
  echo "${ifo} inspiral"
  jid=1
  for allid in `grep inspiral_${ifo}.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep ${ifo} | awk '{print $2}'` ; do 
      sedcmd="s/$replid/insp_${ifo}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done 
  done

  # trigbank
  echo "${ifo} trigbank"
  jid=1
  for allid in `grep trigbank.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep macrooutputifo=\"${ifo}\"  | awk '{print $2}'` ; do 
      sedcmd="s/$replid/trigbank_ifo_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done 
  done

  # inspiral veto 
  echo "${ifo} inspiral veto"
  jid=1
  for allid in `grep inspiral_veto_${ifo}.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep ${ifo} | awk '{print $2}'` ; do 
      sedcmd="s/$replid/insp_${ifo}_Veto_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done
  done

  # sinca 
  echo "${ifo} sinca"
  jid=1
  for allid in `grep s_inca.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep ${ifo} | awk '{print $2}'` ; do 
      sedcmd="s/$replid/sinca_${ifo}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done 
  done

  # coherent trigbank
  echo "${ifo} trigbank_coherent"
  jid=1
  for allid in `grep trigbank_coherent_${ifo}.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep 'macrooutputifo="${ifo}"' | awk '{print $2}'` ; do 
      sedcmd="s/$replid/trigbank_coherent_${ifo}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp
      mv ${1}.tmp ${1} 
    done 
  done


  # inspiral coherent
  echo "${ifo} insp_coh"
  jid=1
  for allid in `grep insp_coh_${ifo}.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep ${ifo}  | awk '{print $2}'` 
    do sedcmd="s/$replid/insp_coh_${ifo}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done 
  done

done

for ifos in {G1H1,G1H2,G1L1,G1V1,H1H2,H1L1,H1V1,H2L1,H2V1,L1V1,G1H1H2,G1H1L1,G1H1V1,G1H2L1,G1H2V1,G1L1V1,H1H2L1,H1H2V1,H1L1V1,H2L1V1,G1H1H2L1,G1H1H2V1,G1H1L1V1,G1H2L1V1,H1H2L1V1,G1H1H2L1V1}; do
  # thinca
  echo "${ifos} thinca"
  jid=1
  for allid in `grep thinca_${ifos}.sub ${1} | awk '{print $2}'`; do 
    sedcmd="s/$allid/thinca_${ifos}_${jid}/g" 
    jid=$((jid+1)) 
    sed $sedcmd ${1} > ${1}.tmp 
    mv ${1}.tmp ${1} 
  done 

  # thinca slide
  echo "${ifos} thinca_slides"
  jid=1
  for allid in `grep thinca_slides_${ifos}.sub ${1} | awk '{print $2}'`; do 
    sedcmd="s/$allid/thinca_slides_${ifos}_${jid}/g" 
    jid=$((jid+1)) 
    sed $sedcmd ${1} > ${1}.tmp 
    mv ${1}.tmp ${1} 
  done

  # thinca2 
  echo "${ifos} thinca2"
  jid=1
  for allid in `grep thinca2_${ifos}.sub ${1} | awk '{print $2}'`; do
    sedcmd="s/$allid/thinca2_${ifos}_${jid}/g" 
    jid=$((jid+1)) 
    sed $sedcmd ${1} > ${1}.tmp 
    mv ${1}.tmp ${1} 
  done

  # thinca2_slides
  echo "${ifos} thinca2_slides"
  jid=1
  for allid in `grep thinca2_slides_${ifos}.sub ${1} | awk '{print $2}'`; do 
    sedcmd="s/$allid/thinca2_slides_${ifos}_${jid}/g" 
    jid=$((jid+1)) 
    sed $sedcmd ${1} > ${1}.tmp 
    mv ${1}.tmp ${1} 
  done

  # cohbank 
  echo "${ifos} cohbank"
  jid=1
  for allid in `grep cohbank.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep '_$cohbank-' | awk '{print $2}'` ; do 
      sedcmd="s/$replid/cohbank_${ifos}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} >> ${1}.tmp 
      mv ${1}.tmp ${1} 
    done
  done

  # coherent inspiral
  echo "${ifos} coherent inspiral"
  jid=1
  for allid in `grep cohinspiral.sub ${1} | awk '{print $2}'`; do 
    for replid in `grep "VARS $allid" ${1} | grep '_ifos-'  | awk '{print $2}'` ; do 
      sedcmd="s/$replid/cohinspiral_${ifos}_${jid}/g" 
      jid=$((jid+1)) 
      sed $sedcmd ${1} > ${1}.tmp 
      mv ${1}.tmp ${1} 
    done 
  done

done

# sire
echo "sire"
jid=1; for allid in `grep sire.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/sire_${jid}/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "sire_slide"
jid=1; for allid in `grep sire_coire.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/sire_${jid}/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

# coire
echo "coire"
jid=1; for allid in `grep coire.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/coire_${jid}/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done

echo "coire_slide"
jid=1; for allid in `grep coire_slide.sub ${1} | awk '{print $2}'`; do for replid in `grep "VARS $allid" ${1} | grep input | awk '{print $2}'` ; do sedcmd="s/$replid/coire_${jid}/g" ; jid=$((jid+1)) ; sed $sedcmd ${1} > ${1}.tmp ; mv ${1}.tmp ${1} ; done ; done


