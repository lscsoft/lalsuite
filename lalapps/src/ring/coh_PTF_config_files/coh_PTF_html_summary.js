function toggleVisible(division) {
  if (document.getElementById("div_" + division).style.display == "none") {
    document.getElementById("div_" + division).style.display = "block";
    document.getElementById("input_" + division).checked = true;
  } else {
    document.getElementById("div_" + division).style.display = "none";
    document.getElementById("input_" + division).checked = false;
  } 
}
function gotoSection(section) {
  document.getElementById("div_" + section).style.display = "block";
  document.getElementById("input_" + section).checked = true;
  window.location.hash = section;
}

function extractPageName(hrefString) {
  var arr = hrefString.split('/');
  return  (arr.length < 2) ? hrefString : arr[arr.length-2].toLowerCase() + arr[
arr.length-1].toLowerCase();
}

function setActiveMenu(arr, crtPage) {
  for (var i=0; i < arr.length; i++) {    if(extractPageName(arr[i].href) == crtPage) {
      arr[i].parentNode.className = "selected";
    }  }
}

function setPage() {
  hrefString = document.location.href ? document.location.href : document.location;
  if (document.getElementById("menubar") !=null )
    setActiveMenu(document.getElementById("menubar").getElementsByTagName("a"), extractPageName(hrefString));
}

function setSubPage() {  hrefString  = document.location.href ? document.location.href : document.location;
  if (document.getElementById("menubar") !=null )
    setActiveMenu(document.getElementById("menubar").getElementsByTagName("a"), extractPageName(hrefString));
    var path = hrefString.split('/')
    index = path[path.length-2]
    setActiveMenu(document.getElementById("menubar").getElementsByTagName("a"), index)
}


