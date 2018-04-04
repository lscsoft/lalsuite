<?php

require '/usr1/ldbd/glue/var/php/script/time_conv_functions.php';  //include functions
require '/usr1/ldbd/glue/var/php/script/time_conv.php';            //execute actual convertion
require "/usr1/ldbd/glue/var/php/script/glue_env.php";             //set up glue environment
 
$com = "ligolw_dq_query ";
foreach($_GET as $key => $value) {
  $com .= " --" . $key . " " . $value;
}

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
if ($_GET["report"]) {
  sort($output);
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, start of hdr, tier3 navigation bar and setup glue e
  require "/usr1/ldbd/glue/var/php/script/segdb_tools_tier4.php";
  echo "<body id='dq_query'>";
  echo "<p><font color='blue'>There is a script runs periodically coalescing the segments. It is possible that DQ Report Page returns shorter segments than expected if the segments in your specified time range have not been coalesced. When in doubt, please run ligolw_segment_query to double check your results.</font></p>";
  echo "<hr/>";
  echo "<b>report = " . $_GET["report"] . "</b>";
  echo "<pre>";
  echo join("\n", $output);
  echo "</pre>";
  echo "</body></div></div></div></html>";
} else {
  header('Content-Type: text/xml');
  for ($i=0; $i < $count; $i++) {
      echo $output[$i];
  }  
  echo $time_in_comment; // add retrieve time at the end
}

?>
