<html>
<head>
  <!--include the stylesheet ./script/style.css-->
  <?php $uTitle="Scimon Data Quality Flag Entry Page";
  require './scripts/styletitle.php'; ?>
</head>

<body>
<?php
//Table displaying logos and $uTitle
require './scripts/header.php';
?>
<h3>Check current flags</h3>
<a href="listflags.php">List of scimon data quality flags entered into the database</a>
<hr>
<h3>Insert new flag</h3>
     <form action="flagcheck.php" method="post">
        <input type="reset"/ value="Clear insert form">
        <p>Affected Detector: <select name="site">
                    <option value="H1">H1</option>
	            <option value="L1">L1</option>
                 </select></p>
        <p>Flag:
	<select name="flag">
        <option value="SCI-AIRCRAFT_ELOG (Aircraft flyover described in the elog)">SCI-AIRCRAFT_ELOG (Aircraft flyover described in the elog)</option>
        <option value="SCI-CONLOG_TROUBLE_ELOG (Problem with the conlog described in the elog)">SCI-CONLOG_TROUBLE_ELOG (Problem with the conlog described in the elog)</option>
        <option value="SCI-EARTHQUAKE_ELOG (Earthquake described in the elog)">SCI-EARTHQUAKE_ELOG (Earthquake described in the elog)</option>
        <option value="SCI-EXPLOSION_ELOG (Nearby explosion described in the elog)">SCI-EXPLOSION_ELOG (Nearby explosion described in the elog)</option> 
        <option value="SCI-FLAKY_DAQ_ELOG (DAQ problems described in the elog)">SCI-FLAKY_DAQ_ELOG (DAQ problems described in the elog)</option>
        <option value="SCI-FLAKY_SERVO_ELOG (Servo problems described in the elog)">SCI-FLAKY_SERVO_ELOG (Servo problems described in the elog)</option>
        <option value="SCI-GLITCY_DATA_ELOG (Glitchy data cause by a known problem described in the elog)">SCI-GLITCY_DATA_ELOG (Glitchy data cause by a known problem described in the elog)</option>
        <option value="SCI-HEAVY_MACHINERY_ELOG (Heavy machinery operating on site described in elog)">SCI-HEAVY_MACHINERY_ELOG (Heavy machinery operating on site described in elog)</option>
        <option value="SCI-HIGH_MICROSEISMIC_ELOG (High microseismic levels described in the elog)">SCI-HIGH_MICROSEISMIC_ELOG (High microseismic levels described in the elog)</option>
        <option value="SCI-HIGH_WIND_ELOG (High winds described in the elog)">SCI-HIGH_WIND_ELOG (High winds described in the elog)</option>
        <option value="SCI-HUMAN_INTRUSION_ELOG (Human intrusion into the LVEA described in the elog)">SCI-HUMAN_INTRUSION_ELOG (Human intrusion into the LVEA described in the elog)</option>
        <option value="SCI-LOW_POWER_ELOG (low power laser operation described in the elog)">SCI-LOW_POWER_ELOG (low power laser operation described in the elog)</option>
        <option value="SCI-NONSTAND_CONFIG_ELOG (Nonstandard configuration described in the elog)">SCI-NONSTAND_CONFIG_ELOG (Nonstandard configuration described in the elog)</option>
        <option value="SCI-SEISMIC_ELOG (Seismic event of unknown cause described in the elog)">SCI-SEISMIC_ELOG (Seismic event described in the elog)</option>
        <option value="SCI-TRAIN_ELOG (Passing train causing high seismic levels described in the elog)">SCI-TRAIN_ELOG (Passing train causing high seismic levels described in the elog)</option>
        <option value="SCI-TUMBLEWEED_BAILING_ELOG (Tumbleweed bailer operating described in the elog)">SCI-TUMBLEWEED_BAILING_ELOG (Tumbleweed bailer operating described in the elog)</option>
        <option value="SCI-VEHICLE_ELOG (Vehicle driving on site described in the elog)">SCI-VEHICLE_ELOG (Vehicle driving on site described in the elog)</option>
        <option value="SCI-WATER_SKID_ELOG (Water skid described in the elog)">SCI-WATER_SKID_ELOG (Water skid described in the elog)</option>
        <option value="SCI-OTHER_ELOG (Other problem described in elog entry)">SCI-OTHER_ELOG (Other problem described in elog entry)</option>
    </select>
         <p>If you do not see an appropriate flag, please choose SCI-OTHER_ELOG and contact <a href="mailto:detchar@relativity.phys.lsu.edu">the detchar group</a> to suggest a new flag name.
         <p>Brief description (max 255 characters): <textarea rows="3" cols="60" name="comment" /></textarea></p>
 	 <p><input type="radio" checked="checked" name="starttimevsgps" value="time">
  	    Start time: <select name="startmonth"><?php require './scripts/form_month_list.php' ?></select>
                     <select name="startday"><?php require './scripts/form_day_list.php' ?></select>
		     <select name="startyear"><?php require './scripts/form_year_list.php' ?></select>
 		     &nbsp&nbsp&nbsp 		
            <input type="text" name="starthour" value=00 size=2 maxlength=2>
            :<input type="text" name="startmin" value=00 size=2 maxlength=2>
            :<input type="text" name="startsec" value=00 size=2 maxlength=2>
		    &nbsp&nbsp&nbsp
            <select name="starttype">
		    <option value="24Hour">24 Hour</option>
		    <option value="am">am</option>
		    <option value="pm">pm</option>
	    </select>
	    <select name="startzone">
		    <option value="UTC">UTC</option>
		    <option value="Central">Central</option>
		    <option value="Pacific">Pacific</option>
	    </select>
         </p>
	 <p><input type="radio" name="starttimevsgps" value="gps">
	   Start GPS:<input type="text" size=10 maxlength=10 name="startgps">
         </p>


        <p><input type="radio" checked="checked" name="stoptimevsgps" value="time">
	   Stop time: <select name="stopmonth"><?php require './scripts/form_month_list.php' ?></select>
		<select name="stopday"><?php require './scripts/form_day_list.php' ?></select>
		<select name="stopyear"><?php require './scripts/form_year_list.php' ?></select>
		    &nbsp&nbsp&nbsp
	<input type="text" name="stophour" value=00 size=2 maxlength=2>
	:<input type="text" name="stopmin"  value=00 size=2 maxlength=2>
	:<input type="text" name="stopsec"  value=00 size=2 maxlength=2>
		    &nbsp&nbsp&nbsp
	<select name="stoptype">
		    <option value="24Hour">24 Hour</option>
		    <option value="am">am</option>
		    <option value="pm">pm</option>
	</select>
	<select name="stopzone">
		    <option value="UTC">UTC</option>
		    <option value="Central">Central</option>
		    <option value="Pacific">Pacific</option>
	</select>
	 </p>
	<p><input type="radio" name="stoptimevsgps" value="gps">
	 Stop GPS:<input type="text" size=10 maxlength=10 name="stopgps">
       </p>

        <p>Elog URL with detailed description of problem:<input type="text" size=67 name="url" ></p>
	<p>Enter your @LIGO.ORG User Name (e.g albert.einstein):<input type="text" size=50 name="user" value="<?=$_SERVER['REMOTE_USER']?>"></p>
 	<p><input type="submit" value="Submit flag to database"/></p>
     </form>

<?php require './scripts/footer.php'; ?>
</body>
</html>
