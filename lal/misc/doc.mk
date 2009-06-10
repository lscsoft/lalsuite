.dvi-dep: ../include/.dvi-dep ../src/.dvi-dep ../test/.dvi-dep
	@ test -d .adoc || mkdir .adoc
	@ if test ! -f .adoc/main.tex ; then \
	    pkgdir=`cd .. && pwd` ; \
	    package=`basename $$pkgdir` ; \
	    sed "s/@package@/$$package/" $(top_srcdir)/misc/main.tex.in > .adoc/main.tex ; \
	  fi
	@ for file in $(DOCSOURCEFILES) ; do \
	    filepath=`cd $(srcdir) && pwd`/$$file ; \
	    if test ! -f $(top_builddir)/doc/autodoc/$$file ; then \
	      cmd="$(LN_S) $$filepath $(top_builddir)/doc/autodoc" ; \
	      $$cmd ; \
	    fi ; \
	    if test ! -f .adoc/$$file ; then \
	      cmd="$(LN_S) $$filepath .adoc" ; \
	      $$cmd ; \
	    fi ; \
	  done
	@ cd .adoc && ${LATEX} main ${TEXFLAGS} && ${MKIND} main \
	  && ${LATEX} main ${TEXFLAGS} && ${LATEX} main ${TEXFLAGS}
	@ for file in main.dvi main.pdf main.ps ; do \
	    if test -f .adoc/$$file ; then \
	      test -f $$file || $(LN_S) .adoc/$$file . ; \
	    fi; \
	  done
	@ touch .dvi-dep

../include/.dvi-dep:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	(cd ../include && $(MAKE) $(AM_MAKEFLAGS) dvi ) \
         || case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
	test -z "$$fail"
../src/.dvi-dep:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	(cd ../src && $(MAKE) $(AM_MAKEFLAGS) dvi ) \
         || case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
	test -z "$$fail"
../test/.dvi-dep:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	(cd ../test && $(MAKE) $(AM_MAKEFLAGS) dvi ) \
         || case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
	test -z "$$fail"
