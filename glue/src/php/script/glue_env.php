<?php
$gluepath = getenv('GLUEPATH');
$pythonpath = getenv('PYTHONPATH');
$ldlibpath = getenv('LD_LIBRARY_PATH');
$ldbdserver = getenv('LDBD_SERVER');
$x509_cert = getenv('X509_USER_CERT');
$x509_key = getenv('X509_USER_KEY');
#$logname = getenv('LOGNAME');
$logname = "ldbd";
$onlinedq = "/online/DQ";

putenv("PYTHONPATH=" . $pythonpath);
putenv("LD_LIBRARY_PATH=" . $ldlibpath);
putenv("LDBD_SERVER=" . $ldbdserver);
putenv("X509_USER_CERT=" . $x509_cert);
putenv("X509_USER_KEY=" . $x509_key);
putenv("LOGNAME=" . $logname);
putenv("ONLINEDQ=" . $onlinedq);
?>
