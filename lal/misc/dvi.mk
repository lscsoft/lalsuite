TEXSOURCES = $(SOURCES) $(HEADERS) $(EXTRA_TEXSOURCES)
.dvi-dep: $(TEXSOURCES)
	@ test -d ../doc/.adoc || mkdir ../doc/.adoc
	-rm -f *.tex
	@ for file in $(TEXSOURCES) ; do \
	    echo "processing $$file with laldoc" ; \
	    $(top_builddir)/doc/laldoc/laldoc $$file laldoc.log ../doc/.adoc $(srcdir) ; \
	  done
	@ for file in *.tex ; do \
	    filepath="`pwd`/$$file" ; \
	    if test ! -f $(top_builddir)/doc/autodoc/$$file ; then \
	      cmd="$(LN_S) $$filepath $(top_builddir)/doc/autodoc" ; \
	      $$cmd ; \
	    fi ; \
	    if test ! -f ../doc/.adoc/$$file ; then \
	      cmd="$(LN_S) $$filepath ../doc/.adoc" ; \
	      $$cmd ; \
	    fi ; \
	  done
	@ touch .dvi-dep
