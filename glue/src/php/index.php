<html>
<head>
  <?php
  $uTitle="HOME";
  require "/usr1/ldbd/glue/var/php/script/styletitle.php";
  ?>
</head>

<body class="home">
<?php require '/usr1/ldbd/glue/var/php/script/header.php'; ?>

<div id="hdrTier4">
<ul>
<li><a href="index_test.php#annoucement">Annoucement</a></li>
<li><a href="index_test.php#intro">Introduction</a></li>
<li><a href="index_test.php#design">Design of the Segment Database</a></li>
</ul>
</div>


<div style="padding:10px;">

<a name="annoucement"><h2>Annoucement</h2></a>
  <ul>
     <li><b>Mar.22, 2010 </b> Version 2 science, injection, calibrated, badgamma, light segments is automated via cronjob once per day at 5am Pacific time at LHO and LLO. These cronjobs allow a 2 hr latency.</li>
     <li><b>Apr. 2, 2010 </b> ldbd:// and ldbdi:// server will be permanently turned off. They will be replaced by https://segdb.ligo.caltech.edu. This change will not affect the S5 database hosted on ldbd://metaserver.phy.syr.edu</li>
  </ul>
<div class="placeholder"></div>
<hr/>

<a name="intro"><h2>Introduction</h2></a>
The LIGO S6 segment database is at https://segdb.ligo.caltech.edu. This database contains both the online DMT generated DQ flags and offline DQ flags and vetoes (such as SCI scimon flags). Currently there are no alternative S6 segment database available. The segment database is the default source for DQ segments and vetoes with more than a few hours of latency. For online analysis requiring low-latency please point the segment tools to the DMT XML files in $ONLINEDQ in command line tool.
<div class="placeholder"></div>


<a name="design"><h2>Design of the Segment Database</h2></a>
The S6 segment database schema (and hence the XML files) have been substantially simplified for S6. The five tables to be used for S6 are shown below. The process, process_params and segment_definer tables are unchanged from S5, with the exception that the domain column in the process table is used to store the distingished name of the user inserting the data. In S5, segments were either stored as on or off allowing users to distinguish between a segment being off or undefined, due to data not being analyzed. It was foung that the off/undefined query was used substantially less than the on query, and so in S6 this column is eliminated from the segment table. To ensure that off/undefined information is still available, the segment_summary is introduced in S6 to store time intervals when quality flags are defined. This further simplified the query "What versions were defined at what time?" allowing better use of version information in S6.
<div class="placeholder"></div>
<a href="./img/S6.jpg"><img src="./img/S6.jpg" width="60%" height="60%" alt="LIGO S6 Database Schema"></a><br/>
<a href="./img/S6.jpg">(Click to enlarge image)</a>


</div>
</div> <!--/hdr-->
</div> <!--/container-->
</div> <!--/wrapper-->
<?php require "/usr1/ldbd/glue/var/php/script/footer.php"; ?>

</body>
</html>
