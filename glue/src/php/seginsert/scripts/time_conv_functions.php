<?php
// Define leap seconds
   function getleaps() {
      $leaps = array(46828800, 78364801, 109900802, 173059203, 252028804, 315187205, 346723206, 393984007, 425520008, 457056009, 504489610, 551750411, 599184012, 820108813, 914803214);
      return $leaps;
   }

// Test to see if a gps second is a leap second
   function isleap($gpsTime) {
      $isLeap = FALSE;
      $leaps = getleaps();
      $lenLeaps = count($leaps);
      for ($i = 0; $i < $lenLeaps; $i++) {
         if ($gpsTime == $leaps[$i]) {
            $isLeap = TRUE;
         }
      }
      return $isLeap;
   }

// Count number of leap seconds that have passed
   function countleaps($gpsTime, $dirFlag){
      $leaps = getleaps();
      $lenLeaps = count($leaps);
      $nleaps = 0;  // number of leap seconds prior to gpsTime
      for ($i = 0; $i < $lenLeaps; $i++) {
         if (!strcmp('unix2gps', $dirFlag)) {
            if ($gpsTime >= $leaps[$i] - $i) {
               $nleaps++;
            }
         } elseif (!strcmp('gps2unix', $dirFlag)) {
            if ($gpsTime >= $leaps[$i]) {
               $nleaps++;
            }
         } else {
            echo "ERROR Invalid Flag!";
         }
      }
      return $nleaps;
   }


// Convert Unix Time to GPS Time
   function unix2gps($unixTime){
      // Add offset in seconds
      if (fmod($unixTime, 1) != 0) {
         $unixTime = $unixTime - 0.5;
         $isLeap = 1;
      } else {
         $isLeap = 0;
      }
      $gpsTime = $unixTime - 315964800;
      $nleaps = countleaps($gpsTime, 'unix2gps');
      $gpsTime = $gpsTime + $nleaps + $isLeap;
      return $gpsTime;
   }

// Convert GPS Time to Unix Time
   function gps2unix($gpsTime){
     // Add offset in seconds
     $unixTime = $gpsTime + 315964800;
     $nleaps = countleaps($gpsTime, 'gps2unix');
     $unixTime = $unixTime - $nleaps;
     if (isleap($gpsTime)) {
        $unixTime = $unixTime + 0.5;
     }
     return $unixTime;
   }

// Convert a human-formatted time to unix time
   function time2unix($month, $day, $year, $hour, $min, $sec, $zone){
      $monthn = monthnum($month);
      if ($sec == 60) {
          $isLeap = TRUE;
          $sec = $sec - 1;
      } else {
          $isLeap = FALSE;
      }
      if (!strcmp("UTC", $zone)) {
          $unixTime = gmmktime($hour, $min, $sec, $monthn, $day, $year);
      } elseif (!strcmp("Central", $zone)) {
          putenv("TZ=US/Central");
          $unixTime = mktime($hour, $min, $sec, $monthn, $day, $year);
      } elseif (!strcmp("Pacific", $zone)) {
          putenv("TZ=US/Pacific");
          $unixTime = mktime($hour, $min, $sec, $monthn, $day, $year);
      } else {
          echo "ERROR: Invalid Time Zone! \n";
      }
      if ($isLeap) {
          $unixTime = $unixTime + .5;
      }
      return $unixTime;
   }

   function monthnum($month) {
      if (!strcmp("Jan", $month)) {
         $monthn = 1;
      } elseif (!strcmp("Feb", $month)) {
         $monthn = 2;
      } elseif (!strcmp("Mar", $month)) {
         $monthn = 3;
      } elseif (!strcmp("Apr", $month)) {
         $monthn = 4;
      } elseif (!strcmp("May", $month)) {
         $monthn = 5;
      } elseif (!strcmp("Jun", $month)) {
         $monthn = 6;
      } elseif (!strcmp("Jul", $month)) {
         $monthn = 7;
      } elseif (!strcmp("Aug", $month)) {
         $monthn = 8;
      } elseif (!strcmp("Sep", $month)) {
         $monthn = 9;
      } elseif (!strcmp("Oct", $month)) {
         $monthn = 10;
      } elseif (!strcmp("Nov", $month)) {
         $monthn = 11;
      } elseif (!strcmp("Dec", $month)) {
         $monthn = 12;
      } else {
         echo "ERROR: Invalid month! \n";
      }
      return $monthn;
}
?>
