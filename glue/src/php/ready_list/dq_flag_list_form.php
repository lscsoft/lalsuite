<html>
<head>
   <?php;
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  ?>

  <script language="JavaScript" src="../script/optiontransfer.js"></script>
  <script language="JavaScript" src="../script/clear_form.js"></script>
</head>


<body id="dq_flag_list" onLoad="clear_form('dq_flag_list_form')">
<!-- setup header logo and navigation bar -->
<?php
 require "/usr1/ldbd/glue/var/php/script/header.php"; #includes header logo, wrapper, container, start of hdr, tier3 navigation bar and setup glue environment
 require '/usr1/ldbd/glue/var/php/script/ready_list_tier4.php';
?>
</div> <!--/hdr-->
<div class="placeholder"></div>


<form name="dq_flag_list_form" id="dq_flag_list_form" method="GET" action="dq_flag_list.php">
This tool lists all the available flag types and their brief explaination currently in the segment database.<p/>
<table>
<tr>
   <!-- select server -->
   <td><b>Choose server:</b></td>
   <td><?php require "/usr1/ldbd/glue/var/php/script/dropdownserver.php"; ?></td>
</tr>
<tr>
   <!--select ifo -->
   <td><b>Select ifo: </b></td>
   <td>
      <select name="ifos">
         <option value="H1" selected>H1</option>
         <option value="L1">L1</option>
         <option value="V1">V1</option>
         <option value="all">All</option>
      </select>
   </td>
</tr>

<tr>
   <td><b>Filter flag name with string:</b></td>
   <td><input type="text" name="filter"/></td>
</tr>
</table>


<center><input type="submit" value="submit" /></center>
</div><!--/hdr-->
</div><!--/container-->
</div><!--/wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</form>
</body>
</html>
