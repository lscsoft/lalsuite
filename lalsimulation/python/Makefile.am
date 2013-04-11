pkgpython_PYTHON = \
	__init__.py \
	git_version.py

EXTRA_DIST = \
	git_version.py

BUILT_SOURCES =
MOSTLYCLEANFILES =

if HAVE_GIT_REPO
BUILT_SOURCES += git_version.py
MOSTLYCLEANFILES += git_version.py
git_version.py .PHONY : gen-git_version.py
gen-git_version.py :
	$(AM_V_at)$(GENERATE_VCS_INFO) git_version.py $(srcdir)/git_version.py.git
endif
