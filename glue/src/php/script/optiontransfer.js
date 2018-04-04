function move_first_to_second(id,value) {
  var len = document.getElementById(id).length;
  document.getElementById(id).options[len] = new Option(value,value,true,true);
}

function remove_elem(id,idx) {
  var len = document.getElementById(id).length;
  if (len <= 0) return;
  var combo = document.getElementById(id);
  var j=0;
  for (var i = 0; i < combo.options.length; i++) {
    if (i == idx) {
      combo.options[i].selected = false;
      combo.options[i] = null;
    } else {
      document.getElementById(id).options[j] = new Option(combo.options[i].value);
      j++;
    }
  }
}

function selectItemsInListbox(listbox)
{
  for (var i=0; i<listbox.options.length; i++)
      {
        listbox.options[i].selected = true;
      }
}
