function setup() {
    parseTables();

    /* simplify table titles (disabled)
    $('h2').each(function() {
       var s = $(this).html().replace(/^(.*):table/,'$1');
       $(this).html(s);
    }); */

    // simplify column titles
    $('th').each(function() {
       var s = $(this).html().replace(/^.*:(.*)\s(.*)/,'$1 $2');
       $(this).html(s);
    });

    // make table titles toggle their table
    $('h2.title').click(function() {
        $(this).parents('div').find('table').toggle(100);
    });

    // make toc titles toggle and go their table
    $('div.toc li').click(function() {
        // go to the table
        window.location = '#' + $(this).html();
        
        var table = $("div.table a[name='" + $(this).html() +"']").parents('div.table').find('table.content');
        table.toggle(100);
        
        // display the table if hidden
        // if(table.is(':hidden')) { table.show(100); }
    });
    
    // color every other table row
    $('table tr:odd').addClass('blue');

    // add thumbnails for images
    $('table.content a').each(function() {
        var href = $(this).attr('href');

        if(href.match(/.*\.png/i)) {
            var td = $(this).parent();            
            td.prepend('<a href="' + href + '" target="_blank"><img src="' + href + '" width="' + td.width() + '"/></a>');
        }
    });
    
    // hide all tables
    $('.content').hide();
}

function parseToken(ret) {
    // strip string
    ret = ret.replace(/^\s*([^\s].+)/,'$1');
    ret = ret.replace(/(.+[^\s])\s*$/,'$1');

    // de-escape commas
    ret = ret.replace('\\,',',');

    // add hyperlinks for selected filetypes
    if(ret.match(/^.*\.(xml|txt|png|gz)$/i)) {
        ret = '<a href="' + ret + '" target="_blank">' + ret + '</a>';
    }
    
    // normalize zeros
    if(ret.match(/^0\.0*$/) || ret.match(/^0\.0*e\+?0+$/i)) {
        ret = '0';
    }
    
    return ret;
}

function parseTables() {
    $('pre').each(function() {
        var table = $(this).parent().find('table:first'); // couldn't make siblings() work
        var cols = $('tr:first th',table).size();
        
        var res = '';
        var count = 0;
        
        var tokenizer = new $.tokenizer([
            // may not handle some special cases with escaped commas or quotes
            /\s*"(.*?)",/,  /([^\s].*?),/,  /\s*"(.*?)"/,   /([^\s].*)/
            ],function(text,isToken,regex) {
                if(isToken) {
                    var token = text.replace(regex,'$1');

                    if(count % cols == 0) res += '<tr>';                    
                    res += '<td>' + parseToken(token) + '</td>';
                    count++;
                    if(count % cols == 0) res += '</tr>';   
                }
            });

        var tokens = tokenizer.parse(this.innerHTML);
        table.append(res);
        
        $(this).remove();
    });
}

/**
 * Tokenizer/jQuery.Tokenizer
 * Copyright (c) 2007-2008 Ariel Flesler - aflesler(at)gmail(dot)com | http://flesler.blogspot.com
 * Dual licensed under MIT and GPL.
 * Date: 2/29/2008
 *
 * @projectDescription JS Class to generate tokens from strings.
 * http://flesler.blogspot.com/2008/03/string-tokenizer-for-javascript.html
 *
 * @author Ariel Flesler
 * @version 1.0.1
 */
;(function(){
	
	var Tokenizer = function( tokenizers, doBuild ){
		if( !(this instanceof Tokenizer ) )
			return new Tokenizer( tokenizers, onEnd, onFound );
			
		this.tokenizers = tokenizers.splice ? tokenizers : [tokenizers];
		if( doBuild )
			this.doBuild = doBuild;
	};
	
	Tokenizer.prototype = {
		parse:function( src ){
			this.src = src;
			this.ended = false;
			this.tokens = [ ];
			do this.next(); while( !this.ended );
			return this.tokens;
		},
		build:function( src, real ){
			if( src )
				this.tokens.push(
					!this.doBuild ? src :
					this.doBuild(src,real,this.tkn)
				);	
		},
		next:function(){
			var self = this,
				plain;
				
			self.findMin();
			plain = self.src.slice(0, self.min);
			
			self.build( plain, false );
				
			self.src = self.src.slice(self.min).replace(self.tkn,function( all ){
				self.build(all, true);
				return '';
			});
			
			if( !self.src )
				self.ended = true;
		},
		findMin:function(){
			var self = this, i=0, tkn, idx;
			self.min = -1;
			self.tkn = '';
			
			while(( tkn = self.tokenizers[i++]) !== undefined ){
				idx = self.src[tkn.test?'search':'indexOf'](tkn);
				if( idx != -1 && (self.min == -1 || idx < self.min )){
					self.tkn = tkn;
					self.min = idx;
				}
			}
			if( self.min == -1 )
				self.min = self.src.length;
		}
	};
	
	if( window.jQuery ){
		jQuery.tokenizer = Tokenizer;//export as jquery plugin
		Tokenizer.fn = Tokenizer.prototype;
	}else
		window.Tokenizer = Tokenizer;//export as standalone class
	
})();