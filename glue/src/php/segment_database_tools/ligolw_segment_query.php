<?php

require '/usr1/ldbd/glue/var/php/script/time_conv_functions.php';  //include functions
require '/usr1/ldbd/glue/var/php/script/time_conv.php';            //execute actual convertion
require "/usr1/ldbd/glue/var/php/script/glue_env.php";             //set up glue environment
 
$com = "ligolw_segment_query ";
foreach($_GET as $key => $value) {
  //echo $key . " = " . $value;
  if ($key=="action") { #get action
     $com .= " --" . $value;
  } else if ($key=="segment-url") { #get datasource
     $com .= " --" . $key . " " . $value;
  } else if ($key=="dmt-files") {
     $com .= " --" . $key;
  } else if ($key=="gps-start-time") {
     $com .= " --" . $key . " " . $value;
  } else if ($key=="gps-end-time") {
     $com .= " --" . $key . " " . $value;
  } else if ($key=="include-segments") {
     $ins=explode(",",$value);   # split segments by comma 
     if ($ins){
        $instring  = "";
        $count = count($ins);
        for ($i=0; $i<$count; $i++) {
          if ($i==0) {
             $instring .= $ins[$i];
           }
           else {
             $instring .= ",$ins[$i]";
           }
        }
     #echo $instring . "<br/>";
     $com .= " --" . $key . " " . $instring;
     }
  } else if ($key=="exclude-segments") {
     $exs=explode(",",$value);  # split segments by comma
     if ($exs){
        $exstring  = "";
        $count = count($exs);
        for ($i=0; $i<$count; $i++) {
            if ($i==0) {
               $exstring .= $exs[$i];
            }
            else {
               $exstring .= ",$exs[$i]";
            }
        }
     #echo $exstring;
     $com .= " --" . $key . " " . $exstring;
    }
  }
}
#echo $com;

// comment showed in the source code
date_default_timezone_set('America/New_York');
$current_time = date("D M d, Y H:i:s", time());
$time_in_comment =  '<!-- Retrieved from ' . $segment_url . " on ". $current_time." EST -->";

// execute command
$com = $gluepath . '/' . $com . " 2>&1";
#echo $com;
exec($com, $output, $returnval);
$count = count($output);

// print xml results
header('Content-Type: text/xml');
for ($i=0; $i < $count; $i++) {
      echo $output[$i];
}
   echo $time_in_comment; // add retrieve time at the end


?>
