%define _prefix /opt/lscsoft/glue
%define _sysconfdir %{_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_prefix}')")}


Name: 		glue
Summary:	The Grid LSC User Environment
Version:	1.36
Release:	1.lscsoft
License:	None
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	python python-cjson m2crypto pyxmpp
BuildRequires:  python-devel
Prefix:         %{_prefix}

%description
Glue (Grid LSC User Environment) is a suite of python modules and programs to
allow users to run LSC codes on the grid.

%prep
%setup 

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install -O1 \
        --skip-build \
        --root=%{buildroot} \
        --prefix=%{_prefix}


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{python_sitearch}/glue/
%{_prefix}/bin/
%{_prefix}/etc/
%{_prefix}/var/
%{_prefix}/share/nmi/lalsuite-build*

%changelog
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
