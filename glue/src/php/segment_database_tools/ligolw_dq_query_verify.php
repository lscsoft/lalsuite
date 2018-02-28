<html>
<head>
  <?php 
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  require '/usr1/ldbd/glue/var/php/script/header.php';
  ?>
</head>

<body id="dq_query">
<?php
require "/usr1/ldbd/glue/var/php/script/segdb_tools_tier4.php";

echo "<p><b>Please verify your request before submitting it:</b></p>";

# get segment-url
$array_to_run["segment-url"]=$_GET["segment_url"];

# get action and its gps time
require '../script/time_conv_functions.php';  //include functions
require '../script/time_conv.php';            //execute actual convertion
#echo $startgps . "<p/>";
$action=explode('_', $_GET['action']);
$array_to_run[$action[0]]=$startgps;


# get --include-segments
if($_GET['include_segments']) {
  $ins=$_GET["include_segments"];
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
}

# get paddings
if ($_GET["start_pad"]) {
   if (!is_numeric($_GET["start_pad"])) {
      $error += 1;
      echo "start-pad = ". $_GET["start_pad"] . "<font color='red'> padding must be number(s)</font><p/>";
   } else {
      $array_to_run["start-pad"] = $_GET["start_pad"];
   }
}

if ($_GET["end_pad"]) {
   if (!is_numeric($_GET["end_pad"])) {
      $error += 1;
      echo "end-pad = " . $_GET["end_pad"] ." <font color='red'> padding must be number(s)</font><p/>";
   } else {
      $array_to_run["end-pad"] = $_GET["end_pad"];
   }
}



# get --in-segments-only
if ($_GET['action'] == "report") {
   if ($_GET["insegmentsonly"]) {
      $array_to_run["in-segments-only"]="";
   }
}


# construct url to run the actual command
$url=http_build_query($array_to_run);
$url="https://" . $_SERVER['SERVER_NAME'] . "/segment_database_tools/ligolw_dq_query.php?" . $url;

# construct equivalent query on command line
$cml="ligolw_dq_query ";
foreach($array_to_run as $key => $value) {
  if ($key == "in-segments-only") {
     echo $key . " = yes " . "<p/>";
  } elseif ($key=="active" || $key=="defined" || $key=="report") {
     echo $key . " = " . $value . " :: " . $starttime. "<p/>";
  } else {
     echo $key . " = " . $value . "<p/>";
  }
  $cml .= " --" . $key . " " . $value;
}


if ($error == 0) {
   echo "<hr/>";
   echo "<b>Equivalent query on command line:</b><br/>";
   echo $cml;
   echo "<hr/>";
   echo "<center>";
   echo '<input type="button" value="Cancel" onClick="history.go(-1);return true;">';
   echo "&nbsp;&nbsp;";
   echo '<input type="button" value="Submit" onClick="location.href=' ."'" .$url.   "'" . '";>'; 
} else {
   echo "<hr/>";
   echo "<center>";
   echo '<input type="button" value="Back" onClick="history.go(-1);return true;">';
}
echo "</center>";
echo "<p/>";
?>
</div><!--/hdr-->
</div><!--container-->
</div><!--/wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</body>
</html>
