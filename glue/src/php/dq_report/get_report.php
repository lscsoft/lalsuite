<?php
echo chr(60).chr(63).'xml version="1.0" encoding="utf-8" '.chr(63).chr(62);
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <?php $uTitle="DQ Report Page";
  require '../seginsert/scripts/styletitle.php'; ?>
</head>

<body>
<?php require './header.php';

  function getrepo($time,$isi){
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
      // Added -i flag, we should probably make this an option on the initial page.       
      if ($isi=="1") {
      $com = "ligolw_dq_query --segment-url https://segdb.ligo.caltech.edu -i --report ".$time ." 2>&1";
      $com = $gluepath . '/' . $com;
      } else {
      $com = "ligolw_dq_query --segment-url ldbd://segdb.ligo.caltech.edu:30020 --report ".$time ." 2>&1";
      $com = $gluepath . '/' . $com;
      }
?>
<p><pre>
<?php
echo "<br>time = ". $time."</br>";
?>
</pre><hr>
<?php
      exec($com, $output, $returnval);
      if($returnval==0) {
       sort ($output);
      } else {
       echo "<h3><font color='red'><center>Submit failed!</font></center></h3>";
      }     
      echo "<p><pre>";
      echo join("\n",$output);
      echo "</pre>";
  }
?>

<?php require '../seginsert/scripts/time_conv_functions.php' ?>

<?php

  $isi = $_POST['isi'];

  //convert start time
  if (!strcmp("time", $_POST['timevsgps'])) {
            if (!strcmp("pm", $_POST['type'])){
                $hour = (int)$_POST['hour'] + 12;
            } else {
                $hour = $_POST['hour'];
            }
            $month = $_POST['month'];
            $day = $_POST['day'];
            $year = $_POST['year'];
            $min = $_POST['min'];
            $sec = $_POST['sec'];
            $zone = $_POST['zone'];
            $time = $month.' '.$day.', ';
            $time .= $year.'&nbsp&nbsp '.$hour.':';
            $time .= $min.':'.$sec.'&nbsp&nbsp ';
            $time .= $zone;
            $UnixTime = time2unix($month, $day, $year, $hour, $min, $sec, $zone);
            $gps = unix2gps($UnixTime);
  } else {
            $gps = $_POST['gps'];
            $UnixTime = gps2unix($gps);
            $time = gmstrftime("%b %d, %Y&nbsp&nbsp %H:%M:%S", $UnixTime);
            $time .= '&nbsp&nbsp UTC';
        }

  // validate start and stop time
  if(strlen($gps)!=9)
      {
         echo "<p><font color='red'>GPS time must be 9 digits</font></p>";
      }
  else{
    getrepo($gps,$isi);
  }
?>


</body>
</html>
