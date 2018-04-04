<html>
<head>
<?php
$uTitle="Data Quality Flag Entry - Flag Info Check";
require './scripts/styletitle.php';
?>
</head>
<body>

<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
?>

<!-- Write info to file -->
<?php

   function addFlag($flagFile, $flagData){
      // Store file contents in an array
      // Open file for output
      if(($fh = fopen($flagFile, 'w+')) == FALSE){
         die('Failed to open ' . $flagFile . ' for writing!');
      }
      $flagData = "<?xml version='1.0'?>\n".
                  "<!DOCTYPE LIGO_LW SYSTEM 'http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt'>\n".
                  "<LIGO_LW>\n".
                  "$flagData".
                  "</LIGO_LW>";
              
      fwrite($fh, $flagData); 

      // insert xml file to the database
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

      $com = "dmtdq_seg_insert --file ".$flagFile . " 2>&1";
      $com = $gluepath . '/' . $com;
      exec($com, $output, $returnval);

      if($returnval==0) {
        fwrite($fh, $flagData);
        echo "<h3><center>Flag Submitted</center></h3>";
        echo '<center>You can check your flag in the <a href="listflags.php">list of all flags in the database.</a></center>';
      }
      else {
        echo "<h3><font color='blue'><center>Submit failed!</font></center></h3>";
        echo "<tt>\n";
        echo join('',$output);
        echo "</tt>\n";
      }
}

   //-------------get value to populate into the xml file----------
   // get flag name and its comment
   $flagarray = explode("(",$_POST[flag]);
   $name = rtrim($flagarray[0]);
   $seg_def_comment = $flagarray[1];

   // construct result filename
   $duration = (int)$_POST['stopgps'] - (int)$_POST['startgps'];
   $flagdatadir = getenv('SEGMENT_DATA_DIR');
   $filename = $_POST['site']."-"."SCIMON_DQ_".$name."-".$_POST['startgps']."-".$duration.".xml";
   if(eregi('/',$filename))
     {
       die("Error: attempt to insert a slash into DQ filename: request denied");
     }
   $filename = $flagdatadir."/".$filename;

   // construct username
   $split_username = explode("@",$_POST[user]);
   $username = $split_username[0]."@ligo.org";
   $pid = getmypid();
   $node = gethostbyaddr(getenv('REMOTE_ADDR'));

   $brief_desc = htmlspecialchars($_POST[comment]);

   $process ="  <Table Name='process:table'>\n".
"    <Column Name='process:program' Type='lstring'/>\n".
"    <Column Name='process:version' Type='lstring'/>\n".
"    <Column Name='process:cvs_repository' Type='lstring'/>\n".
"    <Column Name='process:comment' Type='lstring'/>\n".
"    <Column Name='process:node' Type='lstring'/>\n".
"    <Column Name='process:username' Type='lstring'/>\n".
"    <Column Name='process:unix_procid' Type='int_4s'/>\n".
"    <Column Name='process:start_time' Type='int_4s'/>\n".
"    <Column Name='process:end_time' Type='int_4s'/>\n".
"    <Column Name='process:process_id' Type='ilwd:char'/>\n".
"    <Column Name='process:ifos' Type='lstring'/>\n".
"    <Stream Name='process:table' Type='Local' Delimiter=','>\n".
'    "submitflag.php"'.','.
' "1.0"'.','.
' "/lalsuite/glue/php/submitflag.php"'.','.
'"'.$brief_desc.'"'.','.
'"'.$node.'"'.','.
'"'.$_POST[user].'"'.','.
" $pid".','.
" $_POST[startgps]".','.
" $_POST[stopgps]".','.
' "process:process_id:0"'.','.
' "'.$_POST[site].'"'."\n".
'    </Stream>'."\n".
'  </Table>'."\n";


  $segment_definer = " <Table Name='segment_definer:table'>\n".
"   <Column Name='segment_definer:process_id' Type='ilwd:char'/>\n".
"   <Column Name='segment_definer:segment_def_id' Type='ilwd:char'/>\n".
"   <Column Name='segment_definer:ifos' Type='lstring'/>\n".
"   <Column Name='segment_definer:name' Type='lstring'/>\n".
"   <Column Name='segment_definer:version' Type='int_4s'/>\n".
"   <Column Name='segment_definer:comment' Type='lstring'/>\n".
"   <Stream Name='segment_definer:table' Type='Local' Delimiter=','>\n".
' "process:process_id:0"'.','.
' "segment_definer:segment_def_id:0"'.','.
' "'.$_POST[site].'"'.','.
'"'.$name.'"'.",".
" 1".",".
'"'.$seg_def_comment.'"'.
"   </Stream>\n".
" </Table>\n";

$segment = " <Table Name='segment:table'>\n".
"   <Column Name='segment:segment_id' Type='ilwd:char'/>\n".
"   <Column Name='segment:start_time' Type='int_4s'/>\n".
"   <Column Name='segment:end_time' Type='int_4s'/>\n".
"   <Column Name='segment:segment_def_id' Type='ilwd:char'/>\n".
"   <Column Name='segment:process_id' Type='ilwd:char'/>\n".
"   <Stream Name='segment:table' Type='Local' Delimiter=','>\n".
' "segment:segment_id:0"'.','.
" $_POST[startgps]".",".
" $_POST[stopgps]".",".
' "segment_definer:segment_def_id:0"'.','.
' "process:process_id:0"'."\n".
"   </Stream>\n".
" </Table>\n";

$url_string = htmlspecialchars($_POST[url]);

$segment_summary = " <Table Name='segment_summary:table'>\n".
"   <Column Name='segment_summary:start_time' Type='int_4s'/>\n".
"   <Column Name='segment_summary:end_time' Type='int_4s'/>\n".
"   <Column Name='segment_summary:comment' Type='lstring'/>\n".
"   <Column Name='segment_summary:segment_def_id' Type='ilwd:char'/>\n".
"   <Column Name='segment_summary:process_id' Type='ilwd:char'/>\n".
"   <Column Name='segment_summary:segment_sum_id' Type='ilwd:char'/>\n".
"   <Stream Name='segment_summary:table' Type='Local' Delimiter=','>\n".
" $_POST[startgps]".",".
" $_POST[stopgps]".",".
'"'.$url_string.'"'.','.
' "segment_definer:segment_def_id:0"'.','.
' "process:process_id:0"'.','.
' "segment_summary:segment_sum_id:0"'."\n".
"   </Stream>\n".
" </Table>\n";

$flagData = $process.$segment_definer.$segment.$segment_summary;
addFlag($filename, $flagData);
?>
    <p>
     <center><input type="button" value="Enter Another Flag" onclick="history.go(-2)"></center>
    </p>

<?php require './scripts/footer.php'; ?>
</body>
</head>
