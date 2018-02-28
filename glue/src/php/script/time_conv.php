<?php
//convert start time
        if (!strcmp("time", $_GET['starttimevsgps'])) {
            if (!strcmp("pm", $_GET['starttype'])){
                $starthour = (int)$_GET['starthour'] + 12;
            } else {
                $starthour = $_GET['starthour'];
            }
            $startmonth = $_GET['startmonth'];
            $startday = $_GET['startday'];
            $startyear = $_GET['startyear'];
            $startmin = $_GET['startmin'];
            $startsec = $_GET['startsec'];
            $startzone = $_GET['startzone'];
            $starttime = $startmonth.' '.$startday.', ';
            $starttime .= $startyear.'&nbsp&nbsp '.$starthour.':';
            $starttime .= $startmin.':'.$startsec.'&nbsp&nbsp ';
            $starttime .= $startzone;
            $startUnixTime = time2unix($startmonth, $startday, $startyear, $starthour, $startmin, $startsec, $startzone);
            $startgps = unix2gps($startUnixTime);
        } else {
            $startgps = $_GET['startgps'];
            $startUnixTime = gps2unix($startgps);
            $starttime = gmstrftime("%b %d, %Y&nbsp&nbsp %H:%M:%S", $startUnixTime);
            $starttime .= '&nbsp&nbsp UTC';
        }
#        echo $startgps . "<br/>";
  // convert stop time
        if (!strcmp("time", $_GET['stoptimevsgps'])){
            if (!strcmp("pm", $_GET['stoptype'])){
                $stophour = (int)$_GET['stophour'] + 12;
            } else {
                $stophour = $_GET['stophour'];
            }
            $stopmonth = $_GET['stopmonth'];
            $stopday = $_GET['stopday'];
            $stopyear = $_GET['stopyear'];
            $stopmin = $_GET['stopmin'];
            $stopsec = $_GET['stopsec'];
            $stopzone = $_GET['stopzone'];
            $stoptime = $stopmonth.' '.$stopday.', ';
            $stoptime .= $stopyear.'&nbsp&nbsp '.$stophour.':';
            $stoptime .= $stopmin.':'.$stopsec.'&nbsp&nbsp ';
            $stoptime .= $stopzone;
            $stopUnixTime = time2unix($stopmonth, $stopday, $stopyear, $stophour, $stopmin, $stopsec, $stopzone);
            $stopgps = unix2gps($stopUnixTime);
        } else {
            $stopgps = $_GET['stopgps'];
            $stopUnixTime = gps2unix($stopgps);
            $stoptime = gmstrftime("%b %d, %Y&nbsp&nbsp %H:%M:%S", $stopUnixTime);
            $stoptime .= '&nbsp&nbsp UTC';
        }
 #       echo $stopgps . "<br/>";
?>
