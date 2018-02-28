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

/*Contains functions:
	getleaps()
	isleap()
	function countleaps($gpsTime, $dirFlag)
	function unix2gps($unixTime)
	function gps2unix($gpsTime)
	function time2unix($month, $day, $year, $hour, $min, $sec, $zone)
	function monthnum($month)
*/
require './scripts/time_conv_functions.php';

?>
    <h3><center>Please verify the information you entered</center></h3>

<?php
  //convert start time
        if (!strcmp("time", $_POST['starttimevsgps'])) {
            if (!strcmp("pm", $_POST['starttype'])){
                $starthour = (int)$_POST['starthour'] + 12;
            } else {
                $starthour = $_POST['starthour'];
            }
            $startmonth = $_POST['startmonth'];
            $startday = $_POST['startday'];
            $startyear = $_POST['startyear'];
            $startmin = $_POST['startmin'];
            $startsec = $_POST['startsec'];
            $startzone = $_POST['startzone'];
            $starttime = $startmonth.' '.$startday.', ';
            $starttime .= $startyear.'&nbsp&nbsp '.$starthour.':';
            $starttime .= $startmin.':'.$startsec.'&nbsp&nbsp ';
            $starttime .= $startzone;
            $startUnixTime = time2unix($startmonth, $startday, $startyear, $starthour, $startmin, $startsec, $startzone);
            $startgps = unix2gps($startUnixTime);
        } else {
            $startgps = $_POST['startgps'];
            $startUnixTime = gps2unix($startgps);
            $starttime = gmstrftime("%b %d, %Y&nbsp&nbsp %H:%M:%S", $startUnixTime);
            $starttime .= '&nbsp&nbsp UTC';
        }

  // convert stop time
        if (!strcmp("time", $_POST['stoptimevsgps'])){
            if (!strcmp("pm", $_POST['stoptype'])){
                $stophour = (int)$_POST['stophour'] + 12;
            } else {
                $stophour = $_POST['stophour'];
            }
            $stopmonth = $_POST['stopmonth'];
            $stopday = $_POST['stopday'];
            $stopyear = $_POST['stopyear'];
            $stopmin = $_POST['stopmin'];
            $stopsec = $_POST['stopsec'];
            $stopzone = $_POST['stopzone'];
            $stoptime = $stopmonth.' '.$stopday.', ';
            $stoptime .= $stopyear.'&nbsp&nbsp '.$stophour.':';
            $stoptime .= $stopmin.':'.$stopsec.'&nbsp&nbsp ';
            $stoptime .= $stopzone;
            $stopUnixTime = time2unix($stopmonth, $stopday, $stopyear, $stophour, $stopmin, $stopsec, $stopzone);
            $stopgps = unix2gps($stopUnixTime);
        } else {
            $stopgps = $_POST['stopgps'];
            $stopUnixTime = gps2unix($stopgps);
            $stoptime = gmstrftime("%b %d, %Y&nbsp&nbsp %H:%M:%S", $stopUnixTime);
            $stoptime .= '&nbsp&nbsp UTC';
        }
?>
<p> Affected Detector = <?php echo htmlspecialchars($_POST['site']); ?> </p>
    <?php $site = htmlspecialchars($_POST['site']); ?>
<p> Flag = <?php 
	if (!strcmp("OTHER_ELOG", $_POST['flag'])){
	    $flag = htmlspecialchars($_POST['otherflag']);
	} else {
	    $flag = $_POST['flag'];
	}
	echo $flag
 ?> </p>
<p>



<?php
$error = 0;
    // validate comment file
    $comment = $_POST['comment'];

    if(strlen($comment)==0 || strlen($comment)>255)
      {
         $error = $error + 1;
         echo "<p>Short description = <font color='red'>Comment is required, and length must not exceed 255 characters</font></p>";
      }
    else
      {
         $comment=htmlspecialchars($_POST['comment']);
         echo "<p>Short description = $comment</p>";
      }

   // validate start and stop time
    if($startgps>=$stopgps)
      {
         $error = $error + 1;
         echo "<p>Start time = <font color='red'>Start time must be earlier than stop time</font><p/>";
      }
    else
      {
         echo "<p>Start time = $starttime  ::  $startgps</p>";
         echo "<p>Stop time = $stoptime  ::  $stopgps</p>";
      }


   // validate elog url
    $url = $_POST['url'];
    // check url length
    if(strlen($url)>255)
      {
        $error = $error + 1;
        echo "<p>Elog url = <font color='red'>Length must not exceed 255 characters</font></p>";
      }

    // check url presence
    if(strpos($url,'http') !== 0)
      {
         echo "<p>Elog url = <font color='red'>Warning: no elog entry specified. Please give an elog URL, if possible.</font><p/>";
      }
    if(strpos($url,'http') === 0 && strlen($url)<=255)
      {
         $url=htmlspecialchars($_POST['url']);
         echo "<p>Please click on this <a href='$url' get='_blank'>link</a> to verify the elog url (opens in new window).</p>";
      }
   
    // validate username
    $username = $_POST['user'];
    if(strlen($username)==0 || strlen($username)>64) 
      {
         $error = $error + 1;
         echo "<p>User Name = <font color='red'>username must be specified, and length must not exceed 64 characters</font></p>";
      }
    else
      {
        $userlist = split( '@', $username );
        $user = $userlist[0];
        if(strlen($user)==0 || strpos($user,'.')==FALSE)
          {
             echo "<p>User Name = <font color='red'>username must be in the format of albert.einstein@LIGO.ORG or albert.einstein</font></p>";
             $error = $error + 1;
          }
        else
          {
            echo "User Name = $user";
          }
      }

?>
<p><center>
<form action="submitflag.php" method="post">
<input type="hidden" name="site" value="<?php echo $site; ?>">
<input type="hidden" name="flag" value="<?php echo $flag; ?>">
<input type="hidden" name="comment" value="<?php echo $comment; ?>">
<input type="hidden" name="starttime" value="<?php echo $starttime; ?>">
<input type="hidden" name="startgps" value="<?php echo $startgps; ?>">
<input type="hidden" name="stoptime" value="<?php echo $stoptime; ?>">
<input type="hidden" name="stopgps" value="<?php echo $stopgps; ?>">
<input type="hidden" name="url" value="<?php echo $url; ?>">
<input type="hidden" name="user" value="<?php echo $user; ?>">
<?php
if($error!=0)
      {
         echo "<p><center>";
         echo '<input type ="button" value="Back" onclick="history.back()">';
         echo '</form>';
         echo '</center></p>';
      }
else
    {
        echo "<p><center></>";
        echo '<input type ="button" value="Back" onclick="history.back()">';
        echo '<input type="submit" value="Write Flag to Database">';
        echo '</form>';
        echo '</center></p>';
    }
?>



<?php require './scripts/footer.php'; ?>
</body>
</html>
