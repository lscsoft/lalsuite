<?php
$com = "ldbdc --query \"select rtrim(ifos) concat '|' concat rtrim(name) concat '|' concat cast(version as char) as flag from segment_definer order by ifos, name, version\" 2>&1";


// run the command
$flags = "<optgroup label='ifo:flag:version'>";
$com = $gluepath . '/' . $com;
$com = $com . ' | ' . $gluepath . '/' . 'ligolw_print -t segment_definer -c flag';
exec($com, $output, $returnval);
for ($i=0; $i < count($output); $i++) {
    $j = str_replace('|', ':', $output[$i]);
    $flags .= "<option value=\"$j\">$j</option>" . "\n";
}
?>
