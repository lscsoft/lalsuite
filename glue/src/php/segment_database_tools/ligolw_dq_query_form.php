<html>
<head>
  <?php require "/usr1/ldbd/glue/var/php/script/styletitle.php"; //linke to local style sheet ?>
  <script language="JavaScript" src="../script/optiontransfer.js"></script>
  <script language="JavaScript" src="../script/clear_form.js"></script>
  <script language="JavaScript">
	function toggle(value,id) {
		if ((value == "active_all") || (value == "defined_all") || (value == "report")) {
			document.getElementById(id).setAttribute("style","display:none");
                        document.getElementById('insecond').options.length=0;
		} else {
			document.getElementById(id).setAttribute("style","display:block");
		}
	}
  </script>
</head>

<body id="dq_query"  onLoad="clear_form('ligolw_dq_query_form')">
<!-- setup header logo and navigation bar -->
<?php
 require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, tier3 navigation bar and setup glue environment
 require '/usr1/ldbd/glue/var/php/script/segdb_tools_tier4.php';
?>
</div> <!--/hdr-->
<div class="placeholder"></div>

<!-- link to instruction wiki page -->
<div style="text-align:right;">
<b>Click to see <a href="https://www.lsc-group.phys.uwm.edu/daswg/wiki/S6OnlineGroup/ligolw_dq_query">usage instructions on command line and working examples</a></b>
</div>
<div class="placeholder"></div>



<form name="ligolw_dq_query_form" id="ligolw_dq_query_form" method="GET" action="ligolw_dq_query_verify.php">
<!-- select server -->
<div id="hdrTier4">
<b>Step 1. Choose server: </b><?php require "/usr1/ldbd/glue/var/php/script/dropdownserver.php"; ?>
</div>
<div class="placeholder2"></div>


<!-- select action -->
<div id="hdrTier4">
<b>Step 2. Choose action:</b><br/>
<table border=0>
  <tr>
    <td>ACTIVE:</td>
    <td><input type="radio" name="action" value="active_single" CHECKED onclick="toggle(this.value,'flags');"> is a given flag active at a given time?</input></td>
    <td><input type="radio" name="action" value="active_all" onclick="toggle(this.value,'flags');">what flags are active at a given time?</input><td>
  </tr>
  <tr>
    <td>DEFINED:</td>
    <td><input type="radio" name="action" value="defined_single" onclick="toggle(this.value,'flags');"> is a given flag defined at a given time?</input></td>
    <td><input type="radio" name="action" value="defined_all" onclick="toggle(this.value,'flags');">what flags are defined at a given time?</input></td>
  </tr>
  <tr>
    <td rowspan=2 valign="top">REPORT:</td>
    <td colspan=2><input type="radio" name="action" value="report" onclick="toggle(this.value,'flags');"> what's the status of all flags at a given time?</input></td>
  </tr>
  <tr>
    <td colspan=2>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type="checkbox" name="insegmentsonly" value="insegmentsonly" CHECKED>: If checked, report will only return segments that given times were within</input></td>
  </tr>
</table>
</div>
<div class="placeholder2"></div>


<!-- select time -->
<div id="hdrTier4">
<b>Step 2. Specify time range:</b><br/>
<?php require "/usr1/ldbd/glue/var/php/script/dropdowntime_single.php"; ?>
</div>
<div class="placeholder2"></div>

<!--specify padding-->
<div id="hdrTier4">
<b>Step 3.</b> <i> (optional) </i><b>Specify padding in seconds </b><br/>
Start pad: <input type="text" name="start_pad" size=9/>&nbsp;&nbsp;&nbsp;&nbsp;
End pad: <input type="text" name="end_pad" size=9/>
</div>
<div class="placeholder2"></div>


<!-- select include segments -->
<div id="hdrTier4">
<div id="flags">
   <b>Step 4. Choose include segment(s):</b> If no version number is specified, latest version will be returned for all time.<br/>
   <select name="include_first" id="infirst" size="5">
   <?php
     require "/usr1/ldbd/glue/var/php/script/dropdownflags_wv.php"; //list flags for selection
     echo $flags;
   ?>
   </select>
<input type="button" id="move" value=">" onclick="move_first_to_second('insecond',document.getElementById('infirst').value);">
<input type="button" id="remove" value="<" onclick="remove_elem('insecond',document.getElementById('insecond').selectedIndex);">
<select name="include_segments[]" id="insecond" size="5" multiple>
</select>
</div> <!--/flags-->
</div> <!--/hdrTier4-->

<p><center><input type="submit" value="Submit" /></center></p>
</form>
</div> <!--container-->
</div> <!--wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</body>
</html>
