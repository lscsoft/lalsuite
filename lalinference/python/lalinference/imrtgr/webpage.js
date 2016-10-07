$(document).ready(function(){
	//large imge show on click the thumbnail image
	$('.graph-image').click(function() {
		$(this).parents('td').find('.graph-image-big').show();
	})

	//large image hide
	$('.graph-image-big').click(function() {
		$(this).hide();
	})

	//Scrolling to particular section
	$('a[href^="#"]').on('click', function(event) {
		var target = $( $(this).attr('href') );
		if( target.length ) {
	        event.preventDefault();
	        $('html, body').animate({
	            scrollTop: target.offset().top
	        }, 1000);
	    }

	});
})