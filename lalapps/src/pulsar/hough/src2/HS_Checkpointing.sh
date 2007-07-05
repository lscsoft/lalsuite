cmd=`sed -n 's|.*<a class="code" href="\(.*\)">\(.*\)</a>.*|s%="\2"%="\1"%g;|p' | sort -u`
echo '</body></html>' | cat *.map - | sed "$cmd" | cat HS_Checkpointing.body - > HS_Checkpointing.html
