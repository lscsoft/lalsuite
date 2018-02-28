<?php
require "/usr1/ldbd/glue/var/php/script/glue_env.php"; //setup glue environment

# get ifos
$ifos=$_GET['ifos'];
$segment_url=$_GET['segment_url'];

//Set the content-type header to xml
header("Content-type: text/xml");

$filter=trim($_GET["filter"]);

if ($ifos=="all") {
   $com = "ldbdc --segment-url $segment_url --query \"select rtrim(ifos) concat '|' concat rtrim(name) concat '|' concat cast(version as char) as ifo_name_version, comment from segment_definer ";
   if ($filter) {
      $com .= " WHERE name like '%$filter%' ";
   }
} else {
   $com = "ldbdc --segment-url $segment_url --query \"select rtrim(ifos) concat '|' concat rtrim(name) concat '|' concat cast(version as char) as ifo_name_version, comment from segment_definer where ifos='$ifos' ";
   if ($filter) {
      $com .= " AND name like '%$filter%' ";
   }
}


$com .= " order by ifos, name, version \" 2>&1";
#echo $com;

// run the command
$com = $gluepath . '/' . $com;
#echo $com;
exec($com, $output, $returnval);

// check the return code
if($returnval==0) {
   $count = count($output);
   for ($i = 0; $i < $count; $i++) {
     $j = str_replace('|',':',$output[$i]);
     #echo $output[$i] . "\n";
     echo $j . "\n";
   }
// comment showed in the source code
$servername = $_SERVER['SERVER_NAME'];
date_default_timezone_set('America/New_York');
$current_time = date("D M d, Y H:i:s", time());
echo "<!-- Retrieved from ".$servername. " on ". $current_time." EST -->";
}
else {
   echo chr(60).chr(63).'xml version="1.0" encoding="utf-8" '.chr(63).chr(62);
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
 <head>
 <title>Error querying scimon flags</title>
 </head>
 <body>
 <h1>Error querying scimon flags</h1>
 <tt>
<?php
  echo join('',$output);
?>
 </tt>
 </body>
 </html>
<?php
}
?>

