<html>
<head>
   <?php;
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  ?>

  <script language="JavaScript" src="../script/optiontransfer.js"></script>
  <script language="JavaScript" src="../script/clear_form.js"></script>
</head>


<body id="seg_list" onLoad="clear_form('seg_list_form')">
<!-- setup header logo and navigation bar -->
<?php
 require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, start of hdr, tier3 navigation bar and setup glue environment
 require '/usr1/ldbd/glue/var/php/script/ready_list_tier4.php';
?>
</div> <!--/hdr-->
<div class="placeholder"></div>


<form name="seg_list_form" id="seg_list_form" method="GET" action="seg_list.php">
This tool lists all the currently available segment intervals for the given flag type in the segment database.<p/>
<!-- select server -->
<b>Choose server: </b><?php require "/usr1/ldbd/glue/var/php/script/dropdownserver.php"; ?><p/>


<!--select include segments -->
<b>Select flag: </b><select name="include_segment">
<?php
  require "/usr1/ldbd/glue/var/php/script/dropdownflags_wv.php"; //list flags for selection
  echo $flags;
?>
</select>
</select>
<p/>

<center><input type="submit" value="submit" /></center>
</div><!--/hdr-->
</div><!--container-->
</div><!--wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</form>
</body>
</html>
