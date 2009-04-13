LIGOLW_THINCA_TO_COINC=${1}
type=${2}
LIGOLW_SQLITE=${3}
SQLITE=${4}
cat=${5}

if ${LIGOLW_THINCA_TO_COINC} --ihope-cache=${type}${cat}.cache --veto-segments=vetoes_${cat}.xml.gz --veto-segments-name=vetoes --output-prefix=S5_HM --effective-snr-factor=50; then echo $(date) " : Done with thinca to coinc on FULL DATA" ${cat} > successes.txt; else echo $(date) ": Failed thinca to coinc failed for FULL DATA on " ${type}${cat}.cache > errors.txt; exit; fi;

if ${LIGOLW_SQLITE} -d ${type}${cat}.sqlite -t /tmp -v S5_HM_*${type}*${cat}*.xml.gz;  then echo $(date) " : Done inserting triggers into DB on FULL DATA" ${cat} > successes.txt; else echo $(date) ": Failed Inserting into DB failed for FULL DATA on " ${type}${cat}.cache > errors.txt; exit; fi;

if ${SQLITE} ${type}${cat}.sqlite < simplify.sql;  then echo $(date) " : Done simplifying FULL DATA" ${cat} > successes.txt; else echo $(date) ": Failed Simplifying failed for FULL DATA on " ${type}${cat}.cache > errors.txt; exit; fi;

if ${SQLITE} ${type}${cat}.sqlite < remove_h1h2.sql;  then echo $(date) " : Done removing h1h2 from FULL DATA" ${cat} > successes.txt; else echo $(date) ": Failed Removing H1H2 failed for FULL DATA on " ${type}${cat}.cache > errors.txt; exit; fi;

if ${SQLITE} ${type}${cat}.sqlite < cluster.sql;  then echo $(date) " : Done clustering FULL DATA" ${cat} > successes.txt; else echo $(date) ": Failed Clustering failed for FULL DATA on " ${type}${cat}.cache > errors.txt; exit; fi;

