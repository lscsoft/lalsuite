# FIXME: this command overrides ${prefix}. If you have run `make install`
# beforehand, you may need to run `make clean` for the wheel to build
# successfully.

# NOTE: setting LALSUITE_PYTHON_SYS_PREFIX is required to install Python scripts
# directly into ${prefix}/bin without wrapper scripts; see the
# pybin-install-exec-local rule in gnuscripts/lalsuite_python.am.

wheel: setup.py
	prefix="$(abs_builddir)/build/inst"; \
	$(MAKE) \
		prefix="$${prefix}" \
		pkgdatadir="$${prefix}/share" \
		PYTHON_PREFIX="$${prefix}" \
		PYTHON_EXEC_PREFIX="$${prefix}" \
		LALSUITE_PYTHON_SYS_PREFIX="$${prefix}" \
		-C .. install-strip \
	&& \
	$(PYTHON) setup.py bdist_wheel -d .

clean-local:
	rm -rf build destdir *.whl

EXTRA_DIST = setup.py.in
