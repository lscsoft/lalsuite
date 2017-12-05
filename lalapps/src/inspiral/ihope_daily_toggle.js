function findButtonValue(theForm) {
    var theButtons = theForm.getElementsByTagName('input')

    for (var i = 0; i < theButtons.length; i++) {
        if (theButtons[i].checked) {
            return(theButtons[i].value);
        }
    }
}

function switchView() {
    var theIFO       = findButtonValue(document.getElementById('whichIFO'));
    var theVetoLevel = findButtonValue(document.getElementById('whichVetoLevel'));
    var theContent   = findButtonValue(document.getElementById('whichContents'));
    var theCluster   = findButtonValue(document.getElementById('whichCluster'));

    var name = theIFO + '_CAT' + theVetoLevel + '_' + theCluster + '_' + theContent;
    var divs = document.getElementsByTagName('div'); 

    for (i = 0; i < divs.length; i++) {
        if(divs[i].id[2] == '_') {
            if (divs[i].id == name) {
                divs[i].style.visibility = 'visible';
            } else {
                divs[i].style.visibility = 'hidden'
            }
        }
    }
}


function switchViewIframe() {
    var theIFO       = findButtonValue(document.getElementById('whichIFO'));
    var theVetoLevel = findButtonValue(document.getElementById('whichVetoLevel'));
    var theContent   = findButtonValue(document.getElementById('whichContents'));
    var theCluster   = findButtonValue(document.getElementById('whichCluster'));

    var name = theIFO + '_CAT' + theVetoLevel + '_' + theCluster + '_' + theContent + '.html';

    var theFrame     = document.getElementsByName('theframe')[0]

    theFrame.src     = name
}


