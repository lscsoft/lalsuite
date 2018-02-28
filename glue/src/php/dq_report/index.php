<html>
<head>
  <?php $uTitle="DQ Report Page";
  require '../seginsert/scripts/styletitle.php'; ?>
  <title>DQ Report Page</title>
</head>

<body>
<?php require './header.php'; ?>
<form action="get_report.php" method="post">
  <p><input type="radio" checked="checked" name="timevsgps" value="time">
     Time: 
       <select name="month"><?php require '../seginsert/scripts/form_month_list.php' ?></select>
       <select name="day"><?php require '../seginsert/scripts/form_day_list.php' ?></select>
       <select name="year"><?php require '../seginsert/scripts/form_year_list.php' ?></select>
       &nbsp&nbsp&nbsp
     <input type="text" name="hour" value=00 size=2 maxlength=2>
     :<input type="text" name="min" value=00 size=2 maxlength=2>
     :<input type="text" name="sec" value=00 size=2 maxlength=2>
       &nbsp&nbsp&nbsp
     <select name="type">
        <option value="24Hour">24 Hour</option>
        <option value="am">am</option>
        <option value="pm">pm</option>
     </select>
     <select name="zone">
        <option value="UTC">UTC</option>
        <option value="Central">Central</option>
        <option value="Pacific">Pacific</option>
     </select>
  </p>

  <p><input type="radio" name="timevsgps" value="gps">
      GPS:<input type="text" size=10 maxlength=10 name="gps">
  </p>

  <p>Setting for --in-segments-only flag: <input type="radio" checked="checked" name="isi" value="1">ON
     <input type="radio" name="isi" value="0">OFF
  </p>

<hr/>

<p><center><input type="submit" value="Submit"></center></p>
</form>
</body>
</html>
