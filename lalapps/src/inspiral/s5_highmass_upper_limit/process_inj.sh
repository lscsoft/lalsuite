LIGOLW_THINCA_TO_COINC=${1}
INJ_DESC=${2}
LIGOLW_SQLITE=${3}
INJ_FILE=${4}
SQLITE=${5}
LIGOLW_INSPINJFIND=${6}
cat=${7}


if ${LIGOLW_THINCA_TO_COINC} --ihope-cache=${INJ_DESC}${cat}.cache --veto-segments=vetoes_${cat}.xml.gz --veto-segments-name=vetoes --output-prefix=S5_HM_INJ --effective-snr-factor=50; then echo $(date) " : Done with thinca to coinc on " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed thinca to coinc for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${LIGOLW_SQLITE} -d ${INJ_DESC}${cat}.sqlite -t /tmp -v S5_HM_INJ*${INJ_DESC}*${cat}*.xml.gz; then echo $(date) " : Done adding triggers to DB for " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Adding triggers to DB failed for " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${LIGOLW_SQLITE} -d ${INJ_DESC}${cat}.sqlite -t /tmp -v ${INJ_FILE}; then echo $(date) " : Done adding sims to DB for " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Adding sims failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${SQLITE} ${INJ_DESC}${cat}.sqlite < simplify.sql; then echo $(date) " : Done simplifying " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Simplifying failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${SQLITE} ${INJ_DESC}${cat}.sqlite < remove_h1h2.sql; then echo $(date) " : Done removing h1h2 for " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed removing H1H2 triggers failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${SQLITE} ${INJ_DESC}${cat}.sqlite < cluster.sql; then echo $(date) " : Done clustering " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Clustering failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${LIGOLW_SQLITE} -d ${INJ_DESC}${cat}.sqlite -v -x ${INJ_DESC}${cat}.sqlite.xml.gz; then echo $(date) " : Done converting to XML for " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Converting injections to XML failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;
 
if ${LIGOLW_INSPINJFIND} -v ${INJ_DESC}${cat}.sqlite.xml.gz; then echo $(date) " : Done finding injections for " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Finding injections failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

if ${LIGOLW_SQLITE} -d ${INJ_DESC}${cat}.sqlite -v -r ${INJ_DESC}${cat}.sqlite.xml.gz; then echo $(date) " : Done converting injections back to XML " ${INJ_DESC}${cat}.cache >> successes.txt; else echo $(date) ": Failed Converting XML back to DB failed for injections on " ${INJ_DESC}${cat}.cache >> errors.txt; exit; fi;

