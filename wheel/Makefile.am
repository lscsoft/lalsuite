# FIXME: this command overrides ${prefix}. If you have run `make install`
# beforehand, you may need to run `make clean` for the wheel to build
# successfully.

wheel: setup.py
	$(MAKE) prefix=$(abs_builddir)/build/inst pkgdatadir=$(abs_builddir)/build/inst/share -C .. install-strip && \
	$(PYTHON) setup.py bdist_wheel -d .

clean-local:
	rm -rf build destdir *.whl

EXTRA_DIST = setup.py.in
