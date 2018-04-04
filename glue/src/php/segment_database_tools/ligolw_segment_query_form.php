<html>
<head>
  <?php require '/usr1/ldbd/glue/var/php/script/styletitle.php'; //link to local style sheet ?>
  <script language="JavaScript" src="../script/optiontransfer.js"></script>
  <script language="JavaScript" src="../script/clear_form.js"></script>
  <script language="JavaScript">
        function toggle(value,id) {
                if (value == "show-types")  {
                        document.getElementById(id).setAttribute("style","display:none");                        document.getElementById('insecond').options.length=0;
                } else {
                        document.getElementById(id).setAttribute("style","display:block");
                }
        }
  </script>
</head>



<body id="seg_query"  onLoad="clear_form('ligolw_segment_query_form')">
<!-- setup header logo and navigation bar -->
<?php
 require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, start of hdr, tier3 navigation bar and setup glue environment
 require '/usr1/ldbd/glue/var/php/script/segdb_tools_tier4.php';
?>
</div> <!--/hdr-->
<div class="placeholder"></div>

<div style="text-align:right;">
<b>Click to see <a href="https://www.lsc-group.phys.uwm.edu/daswg/wiki/S6OnlineGroup/ligolw_segment_query">usage instructions on command line and working examples</a></b>
</div>
<div class="placeholder"></div>


<!-- form starts from here -->
<form name="ligolw_segment_query_form" id="ligolw_segment_query_form" method="GET" action="ligolw_segment_query_verify.php">
<div id="hdrTier4">
<b>Step 1: Choose data source accordingly:</b><br/>
<?php require "/usr1/ldbd/glue/var/php/script/datasource.php"; ?>
</div>
<div class="placeholder_bg"></div>


  <!-- select action -->
  <div id="hdrTier4">
  <b>Step 2: Choose action: </b><br/>
  <input type="radio" name="action" value="query-segments" CHECKED onclick="toggle(this.value, 'flags');">QUERY SEGMENTS: when was a given flag(s) active?</input><br/>
  <input type="radio" name="action" value="query-types" onclick="toggle(this.value, 'flags');">QUERY TYPES: when was a given flag(s) defined?</input><br/>
  <input type="radio" name="action" value="show-types" onclick="toggle(this.value, 'flags');">SHOW TYPES: what flags exist in the database in a given time range?</input>
  </div>
  <div class="placeholder2"></div>


  <!-- select time -->
  <div id="hdrTier4">
  <b>Step 3: Select time range: </b><br/>
  <?php require "/usr1/ldbd/glue/var/php/script/dropdowntime.php"; ?>
  </div>
  <div class="placeholder2"></div>

<!-- select include segments -->
<div id="hdrTier4">
<div id="flags">
   <b>Step 4: Choose include segment(s):</b> If no version number is specified, latest version will be returned for all time.<br/>
   <select name="include_first" id="infirst" size="5">
   <?php
     require "../script/dropdownflags_nov.php"; //list flags for selection
     echo $flags;
   ?>
   </select>
   <input type="button" id="move" value=">" onclick="move_first_to_second('insecond',document.getElementById('infirst').value);">
   <input type="button" id="remove" value="<" onclick="remove_elem('insecond',document.getElementById('insecond').selectedIndex);">
   <select name="include_segments[]" id="insecond" size="5" multiple>
   </select>
  <div class="placeholder_bg"></div>

   <!-- select exclude segments -->
   <b>Step 5: Choose exclude segment(s)</b><i><font size="-1"> optional</font></i>: If no version number is specified, latest version will be returned for all time.<br/>
   <select name="exclude_first" id="exfirst" size="5">
   <?php echo $flags; ?>
   </select>
   <input type="button" id="move" value=">" onclick="move_first_to_second('exsecond',document.getElementById('exfirst').value);">
   <input type="button" id="remove" value="<" onclick="remove_elem('exsecond',document.getElementById('exsecond').selectedIndex);">
   <select name="exclude_segments[]" id="exsecond" size="5" multiple>
   </select>
</div> <!--flags-->
</div><!--/hdrTier4-->
<div class="placeholder2"></div>

<p><center>
  <input type="submit" name="xml" id="xml" value="Get xml result" />
  <input type="submit" name="ascii" id="ascii" value="Get ascii result" />
</center></p>
</form>


</div><!--/hdr-->
</div><!--/container-->
</div><!--/wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</body>
</html>
