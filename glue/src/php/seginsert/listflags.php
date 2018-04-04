<?php
//Set the content-type header to xml
header("Content-type: text/xml");
$com = "ldbdc --query \"select segment_definer.ifos,segment_definer.name,segment_summary.start_time,segment_summary.end_time,process.username as scimon,process.comment as scimon_comment, segment_summary.comment as elog_url from segment_summary,process,segment_definer where segment_definer.segment_def_id = segment_summary.segment_def_id and segment_summary.segment_def_cdb = segment_definer.creator_db and process.process_id = segment_summary.process_id and process.creator_db = segment_summary.creator_db and segment_definer.name like 'SCI-%ELOG' order by segment_definer.ifos,segment_definer.name,segment_summary.start_time\" 2>&1";

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
