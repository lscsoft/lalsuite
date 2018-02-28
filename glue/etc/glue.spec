%define _glue_prefix /usr
%define _sysconfdir %{_glue_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?glue_python_sitearch: %define glue_python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_glue_prefix}')")}
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1)")}



Name: 		glue
Summary:	The Grid LSC User Environment
Version:	1.51.0
Release:	1%{?dist}
License:	None
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	python-cjson m2crypto glue-common glue-segments python >= 2.6
BuildRequires:  python-devel
Prefix:         %{_glue_prefix}
%description
Glue (Grid LSC User Environment) is a suite of python modules and programs to
allow users to run LSC codes on the grid.

%package common
Summary:	The common files needed for all sub-packages
Group: 		Development/Libraries
Requires: 	python 
%description common
This is for the files that are common across the glue subpackages, namely git_version, iterutils and __init__.py

%package segments
Summary:        The segments subpackage
Group:          Development/Libraries
Requires:       python glue-common
%description segments
This is for the segments subpackage, written by Kipp.

%prep
%setup 

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install -O1 \
        --skip-build \
        --root=%{buildroot} \
        --prefix=%{_glue_prefix}
rm -rf $RPM_BUILD_ROOT/usr/lib64/python2.?/site-packages/glue-1.51.0-py2.?.egg-info

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{glue_python_sitearch}/glue/
%{_glue_prefix}/bin/*
%exclude %{_glue_prefix}/etc/
%exclude %{_glue_prefix}/var/
%exclude %{_glue_prefix}/share/nmi/lalsuite-build*
%exclude %{glue_python_sitearch}/glue/__init__.py
%exclude %{glue_python_sitearch}/glue/__init__.pyc
%exclude %{glue_python_sitearch}/glue/segments.py
%exclude %{glue_python_sitearch}/glue/iterutils.py
%exclude %{glue_python_sitearch}/glue/git_version.py
%exclude %{glue_python_sitearch}/glue/segments.pyc
%exclude %{glue_python_sitearch}/glue/iterutils.pyc
%exclude %{glue_python_sitearch}/glue/git_version.pyc
%exclude %{glue_python_sitearch}/glue/__segments.so
#%exclude %{_glue_prefix}/src/segments/
#%exclude %{_glue_prefix}/test/segment_verify.py
#%exclude %{_glue_prefix}/test/segmentsUtils_verify.py
#%exclude %{_glue_prefix}/test/verifyutils.py

%files segments
%{glue_python_sitearch}/glue/segments.py
%{glue_python_sitearch}/glue/segments.pyc
%{glue_python_sitearch}/glue/__segments.so
#%{glue_python_sitearch}/src/segments/
#%{glue_python_sitearch}/test/segment_verify.py
#%{glue_python_sitearch}/test/segmentsUtils_verify.py
#%{glue_python_sitearch}/test/verifyutils.py

%files common
%{glue_python_sitearch}/glue/__init__.py
%{glue_python_sitearch}/glue/__init__.pyc
%{glue_python_sitearch}/glue/iterutils.pyc
%{glue_python_sitearch}/glue/iterutils.py
%{glue_python_sitearch}/glue/git_version.py
%{glue_python_sitearch}/glue/git_version.pyc

%changelog
* Thu Jul 23 2015 Ryan Fisher <rpfisher@syr.edu>
- Pre-ER8 release, attempt 2.

* Wed Jul 22 2015 Ryan Fisher <rpfisher@syr.edu>
- Pre-ER8 release.

* Fri May 22 2015 Ryan Fisher <rpfisher@syr.edu>
- ER7 release.

* Wed Nov 19 2014 Ryan Fisher <rpfisher@syr.edu>
- ER6 pre-release bug fix for dmt files method of ligolw_segment_query.

* Thu Nov 13 2014 Ryan Fisher <rpfisher@syr.edu>
- ER6 pre-release.

* Tue May 6 2014 Ryan Fisher <rpfisher@syr.edu>
- Version update to match git tag.

* Tue May 6 2014 Ryan Fisher <rpfisher@syr.edu>
- Sub-version release to add datafind to debian package.

* Tue Dec 3  2013 Ryan Fisher <rpfisher@syr.edu>
- ER5 release.

* Tue Jul 2 2013 Ryan Fisher <rpfisher@syr.edu>
- ER4 release, matching spec file.

* Tue Jul 2 2013 Ryan Fisher <rpfisher@syr.edu>
- ER4 release.

* Thu Mar 7 2013 Ryan Fisher <rpfisher@syr.edu>
- Post ER3 release of glue for pegasus 4.2 transition, added new RFC Proxy
    support.

* Fri Mar 1 2013 Ryan Fisher <rpfisher@syr.edu>
- Post ER3 release of glue for pegasus 4.2 transition.

* Mon Nov 19 2012 Ryan Fisher <rpfisher@syr.edu>
- New Release of glue for ER3 with updates to ligolw and lal codes.

* Tue Sep 4 2012 Ryan Fisher <rpfisher@syr.edu>
- New Release of glue with upgrades and bugfixes to segment database infrastructure. 

* Fri May 18 2012 Ryan Fisher <rpfisher@syr.edu>
- Bugfix release of 1.39 labelled 1.39.2.  Includes fix to ligolw for URL reading, and packaging fixes. 

* Fri May 11 2012 Ryan Fisher <rpfisher@syr.edu>
- Bugfix release of 1.39 labelled 1.39.1 

* Thu May 10 2012 Ryan Fisher <rpfisher@syr.edu>
- New release of glue to replace Apr 12 near-release.  This includes ligolw changes and updates for job submission over remote pools.

* Thu Apr 12 2012 Ryan Fisher <rpfisher@syr.edu>
- New release of glue with updates to ligolw library, including some bug fixes for ligowl_sqlite and ligolw_print.  

* Wed Nov 16 2011 Ryan Fisher <rpfisher@syr.edu>
- New release of glue with glue-segments and glue-common split from glue, lvalerts, lars and gracedb removed.

* Mon Oct 10 2011 Ryan Fisher <rpfisher@syr.edu>
- New release of glue to fix build issues called 1.36. 

* Fri Oct 7 2011 Ryan Fisher <rpfisher@syr.edu>
- New release of glue with Kipp's fixes to ligolw_sqlite bugs, Kipp's checksums added, and Peter and my change to the coalescing script for segment databases.

* Thu Sep 29 2011 Ryan Fisher <rpfisher@syr.edu>
- New release of glue with speedup to string to xml conversion and 10 digit gps fixes.

* Wed Sep 15 2010 Peter Couvares <pfcouvar@syr.edu>
- New release of glue with GEO publishing

* Thu Apr 22 2010 Duncan Brown <dabrown@physics.syr.edu>
- Third S6 release of glue

* Wed Nov 11 2009 Duncan Brown <dabrown@physics.syr.edu>
- Second S6 release of glue

* Mon Jul 27 2009 Duncan Brown <dabrown@physics.syr.edu>
- First S6 release of glue

* Wed Jul 01 2009 Duncan Brown <dabrown@physics.syr.edu>
- Pre S6 release of glue

* Wed Jun 24 2009 Duncan Brown <dabrown@physics.syr.edu>
- Post E14 release of glue

* Tue Jun 11 2009 Duncan Brown <dabrown@physics.syr.edu>
- Allow segment tools to see multiple ifos

* Tue Jun 10 2009 Duncan Brown <dabrown@physics.syr.edu>
- Restored LSCdataFindcheck and fixed debian control files

* Tue Jun 09 2009 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.19-1

* Tue Jun 24 2008 Ping Wei <piwei@syr.edu>
- Build for glue 1.18-1

* Wed Jun 19 2008 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.17

* Fri Nov 04 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Build for glue 1.6

* Tue Aug 23 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Initial build for glue 1.0
