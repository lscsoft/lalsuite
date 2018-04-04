<?php
//Set the content-type header to xml
header("Content-type: text/xml");
$com = "ldbdc --query \"select segment.start_time,segment.end_time from segment,segment_definer where segment.segment_def_id=segment_definer.segment_def_id and segment.segment_def_cdb=segment_definer.creator_db and segment.start_time>=931035615 and segment_definer.version=2 and segment_definer.name='DMT-INJECTION_INSPIRAL' and segment_definer.ifos='L1' order by segment.start_time ASC\" 2>&1";

$gluepath = getenv('GLUEPATH');
$pythonpath = getenv('PYTHONPATH');
$ldlibpath = getenv('LD_LIBRARY_PATH');
$ldbdserver = getenv('LDBD_SERVER');
$x509_cert = getenv('X509_USER_CERT');
$x509_key = getenv('X509_USER_KEY');

putenv("PYTHONPATH=" . $pythonpath);
putenv("LD_LIBRARY_PATH=" . $ldlibpath);
putenv("LDBD_SERVER=" . $ldbdserver);
putenv("X509_USER_CERT=" . $x509_cert);
putenv("X509_USER_KEY=" . $x509_key);

// run the command
$com = $gluepath . '/' . $com;
exec($com, $output, $returnval);

// check the return code
if($returnval==0) {
   $count = count($output);
   for ($i = 0; $i < $count; $i++) {
     echo $output[$i] . "\n";
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

