<html>
<head>
  <?php 
  # set up header and tier3 navigation bar
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  require '/usr1/ldbd/glue/var/php/script/header.php';
  ?>
</head>

<body id="seg_query">
<?php
require "/usr1/ldbd/glue/var/php/script/segdb_tools_tier4.php"; #include local level navigation bar
#################### 1. get argument value ####################
$error=0;

echo "<p><b>Please verify your request before submiting it:</b><p/>";
$array_to_run=array();

# get data source to run query against 
if ($_GET["datasource"] == "database") {
   $array_to_run["segment-url"]=$_GET["segment_url"];
   echo "<b>datasource</b> = segment database at " . $_GET["segment_url"] . "<p/>";
} else {
   $array_to_run["dmt-files"]="";
   echo " <b>datasource</b> = DMT XML files<p/>";
}

# get action
$array_to_run["action"]=$_GET["action"];
if ($_GET["action"] == "query-segments") {
   $hint=" (when was a given flag(s) active?)";
} else if ($_GET["action"] == "query-types") {
   $hint = "(when as a given flag(s) defined?)";
} else {
   $hint = "(what flags exist in the database in a given time range?)";
}
echo "<b>action</b> = " . $_GET["action"] . "<i> $hint</i><p/>";

# get --include-segments argument
$ins=$_GET['include_segments'];
if (($ins==NULL || $ins=="") && ($array_to_run["action"] != "show-types")) {
   $error +=1 ;
   echo "include-segments = <font color='red'> This field is required</font><p/>";
} else if (($ins) && ($array_to_run["action"] != "show-types")){
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
   $array_to_run["include-segments"]=$instring;
   echo "<b>include-segments</b> = " . $instring . "<p/>";
}

# get --exclude-segments argument
$exs=$_GET['exclude_segments'];
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
}
if ($exstring != "") {
   $array_to_run["exclude-segments"]=$exstring;
   echo "<b>exclude-segemnt</b> = " . $exstring . "<p/>";
}


# get --start-gps and --stop-gps
require "/usr1/ldbd/glue/var/php/script/time_conv_functions.php"; //include functions
require "/usr1/ldbd/glue/var/php/script/time_conv.php";           //execute actual convertion
if ($startgps != "") {
   if (strlen($startgps) != 9) {
      $error += 1;
      echo "<b>gps-start-time</b> = " .$startgps. " :: " .$starttime. "<font color='red'> gps time must be 9 digits</font><p/>";
   } else {
      $array_to_run["gps-start-time"] = $startgps;
      echo "<b>gps-start-time</b> = " .$startgps. " :: " .$starttime. "<p/>";   
   }
}
if ($stopgps != "") {
   if (strlen($stopgps) != 9) {
      $error += 1;
      $errmsg = "<b>gps-end-time</b> = " .$stopgps. " :: " .$stoptime. "<font color='red'> gps time must be 9 digits</font></p>";
   }
   if ($startgps > $stopgps) {
      $error += 1;
      $errmsg = "<b>gps-end-time</b> =" .$stopgps. " :: " .$stoptime. "<font color='red'> end time must be greater than start time</font><p/>";
   }
   if ($errmsg != "") {
      echo $errmsg;
   } else {
      $array_to_run["gps-end-time"]=$stopgps;
      echo "<b>gps-end-time</b> = " .$stopgps. " :: " .$stoptime ."<p/>";
   }
}
# get output format 
if ($_GET['xml']) {
  $array_to_run["output-format"]="xml";
} else if ($_GET['ascii']) {
  $array_to_run["output-format"]="ascii";
}
echo "<b>output-format</b> = " . $array_to_run["output-format"] . "<p/>";
########################## build url based on $array_to_run ##############################
$url=http_build_query($array_to_run);
$url="https://" . $_SERVER['SERVER_NAME'] . "/segment_database_tools/ligolw_segment_query.php?" . $url;
#echo $url. "<p/>";


########################### construct equivalent query on cml ###########################
# remove unnecessary arguments from the display array
unset($array_to_run["output-format"]);

# construct cml
$query_on_cml="ligolw_segment_query ";
foreach($array_to_run as $key => $value)
{
  if ($key=="action") {
     $query_on_cml .= " --" .$value;
  } else {
     $query_on_cml .= " --$key" . "  $value";
  }
}

if ($error == 0) {
   echo "<hr/>";
   echo "<b>Equivalent query on command line:</b><br/>";
   echo $query_on_cml;
   echo "<hr/>";
   echo "<center>";
   echo '<input type="button" value="Cancel" onClick="history.go(-1);return true;">';
   echo "&nbsp;&nbsp;";
   echo '<input type="button" value="Submit" onClick="location.href=' ."'" .$url.   "'" . '";>'; 
} else {
   echo "<hr/>";
   echo "<center>";
   echo '<input type="button" value="Back" onClick="history.go(-1);return true;">';
   echo "<p/>"; 
}
echo "</center>";
echo "<p/>";
?>

</div><!--/hdr-->
</div><!--container-->
</div><!--/wrapper-->

<? require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>
</body>
</html>
