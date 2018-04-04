Name: 		python-pylal
Summary:	Python LSC Algorithm Library
Version:	0.10.0
Release:	1%{?dist}
License:	See file LICENSE
Group:		Development/Libraries
Source:		pylal-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/pylal.html
BuildRoot:	%{_tmppath}/pylal-%{version}-root
Requires:	python glue glue-segments lal lal-python lalframe lalframe-python lalsimulation lalsimulation-python lalinspiral lalinspiral-python lalburst lalburst-python numpy scipy python-matplotlib
BuildRequires:  python-devel lal-devel lalmetaio-devel lalframe-devel lalsimulation-devel lalinspiral-devel lalburst-devel numpy pkgconfig
%description
The PyLAL package is best described as the Python LIGO Algorithm Library. It was originally a Python wrapping of parts of the LAL library, and although it continues to provide that function it has acquired a large collection of independent code of its own so that it is no longer exclusively a Python interface to LAL. In this package you will find convenience code to assist with manipulating XML documents using the glue.ligolw I/O library, you will find a wrapping to libframe to enable GWF frame-file reading, you will find binning and smoothing code, and you will find (partial) wrappings of LAL's burstsearch, date, inject, tools, and window packages. Additionally, you will find most, if not all, of the inspiral pipeline's follow-up and summary tools, and several burst-related trigger post-production tools.

%prep
%setup -n pylal-%{version}

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install -O1 \
        --skip-build \
        --root=%{buildroot} \
        --prefix=/usr \
        --record=INSTALLED_FILES

%clean
[ ${RPM_BUILD_ROOT} != "/" ] && rm -Rf ${RPM_BUILD_ROOT}
rm -Rf ${RPM_BUILD_DIR}/pylal-%{version}

%files -f INSTALLED_FILES
%defattr(-,root,root)
