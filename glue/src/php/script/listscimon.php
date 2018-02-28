<?php
require "/usr1/ldbd/glue/var/php/script/glue_env.php";             //set up glue environment

//Set the content-type header to xml
$com = "ldbdc --query \"select distinct(process.username) from process where process_id in (select distinct(segment_summary.process_id) from segment_summary,segment_definer where segment_definer.segment_def_id=segment_summary.segment_def_id and segment_definer.creator_db=segment_summary.segment_def_cdb and segment_definer.name like 'SCI-%') \" 2>&1";


// run the command
$com = $gluepath . '/' . $com;
$com = $com . ' | ' . $gluepath . '/' . 'ligolw_print -t process -c username';
exec($com, $output, $returnval);

$count = count($output);
echo "<select name='scimon'>";
     echo "<option value='all' selected>All</option>";
     for ($i = 1; $i < $count; $i++) {
         echo "<option value='$output[$i]'>$output[$i]</option>";
     }
echo "</select>";
?>

