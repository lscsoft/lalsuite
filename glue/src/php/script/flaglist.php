<?php
//Set the content-type header to xml
$com = "ldbdc --query \"select ifos,name,version from segment_definer order by ifos,name,version\" 2>&1";

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

$count = count($output);
$name_list = array();
$namevar = '';
for ($i = 1; $i < $count; $i++) {
   list($ifos,$flag,$version) = split(":",$output[$i]);
   $name = $ifos.":".$flag;
   if ($name not in $namevar) {
      $namevar.$name;
   #echo $output[$i] . "\n";
   echo $output[$i] . "<p/>";
   }
?>

