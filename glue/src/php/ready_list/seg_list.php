<html>
<head>
  <?php
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  require "/usr1/ldbd/glue/var/php/script/glue_env.php"; //setup glue environment
  ?>
</head>


<body id="seg_list">
<!--setup header logo and navigation bar -->
<?php
 require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, start of hdr, tier3 navigation bar and setup glue environment
 require '/usr1/ldbd/glue/var/php/script/ready_list_tier4.php';
?>
</div> <!--/hdr-->
<div class="placeholder"></div>

<?php
$segment_url=$_GET['segment_url'];

$in=$_GET['include_segment'];
$pieces = explode(":", $in);
echo "<b>segment intervals of " .$in . "</b><p/>";

require "/usr1/ldbd/glue/var/php/script/glue_env.php";

# 2.2 construct command if xml output format is chosen
// comment showed in the source code
date_default_timezone_set('America/New_York');
$current_time = date("D M d, Y H:i:s", time());
$time_in_comment =  '<!-- Retrieved from ' . $segment_url . " on ". $current_time." EST -->";



$com = "ldbdc --segment-url $segment_url --query \"select segment.start_time,segment.end_time from segment,segment_definer where segment.segment_def_id=segment_definer.segment_def_id and segment.segment_def_cdb=segment_definer.creator_db and segment.start_time>=931035615 and segment_definer.version=$pieces[2] and segment_definer.name='$pieces[1]' and segment_definer.ifos='$pieces[0]' order by segment.start_time ASC\"";


// run the command
$com = $gluepath . '/' . $com;
$com = $com . " | " . $gluepath ."/ligolw_print -t segment -c start_time -c end_time -d ', ' | awk '{print $1\" \"$2\" \"$2-$1}' 2>&1";
#$com = $com . " | " . $gluepath ."/ligolw_print -t segment_summary -c start_time -c end_time -d ', ' | awk '{print $1\" \"$2\" \"$2-$1}' 2>&1";
exec($com, $output, $returnval);

// check the return code
if($returnval==0) {
   $count = count($output);
   echo "<p>start-time, &nbsp;&nbsp;&nbsp;end-time duration<br/>";
   echo "=========================<br/>";
   for ($i = 0; $i < $count; $i++) {
     echo $output[$i] . "<br/>";
   }
// comment showed in the source code
$servername = $_SERVER['SERVER_NAME'];
date_default_timezone_set('America/New_York');
$current_time = date("D M d, Y H:i:s", time());
//echo "<!-- Retrieved from ".$servername. " on ". $current_time." EST -->";
}
?>

<p/>
</div><!--/hdr-->
</div><!--/container-->
</div><!--/wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</body>
</html>

