<?php
$com = "ldbdc --query \"select rtrim(ifos) as ifos, rtrim(name) as name, version from segment_definer order by ifos,name,version\"";

// run the command
$flags = "<optgroup label='ifo:flag:version'>";
$com = $gluepath . '/' . $com;
$com = $com . ' | ' . $gluepath . '/' . 'ligolw_print -t segment_definer -c ifos -c name -c version --delimiter ":"';
exec($com, $output, $returnval);
if($returnval==0) {
  $count = count($output);
  $uniqueflags = array();
  for ($i=0; $i < $count; $i++) {
      $pieces = explode(":", $output[$i]);
      $flagname = $pieces[0] . ":" . $pieces[1];
      $key = array_search($flagname, $uniqueflags);
      if ($key) {
         $flags .= "<option value=\"$output[$i]\">$output[$i]</option>" . "\n";
      }
      else {
         $uniqueflags[] = $flagname;
         $flags .= "<option value=\"$flagname\">$flagname</option>" . "\n";  //add flag name w/o version number as a new option
         $flags .= "<option value=\"$output[$i]\">$output[$i]</option>" . "\n";
      }
  }
  $flags .= "</optgroup>";
}

?>
