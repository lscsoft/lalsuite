%define name              ligo-gracedb
%define version           1.7
%define unmangled_version 1.7
%define release           1
 
Summary:   Gravity Wave Candidate Event Database
Name:      %{name}
Version:   %{version}
Release:   %{release}
Source0:   %{name}-%{unmangled_version}.tar.gz
License:   GPL
Group:     Development/Libraries
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Prefix:    %{_prefix}
BuildArch: noarch
Vendor:    Brian Moe <brian.moe@ligo.org>
Requires:  ligo-common m2crypto
Url:       http://www.lsc-group.phys.uwm.edu/daswg/gracedb.html
 
%description
The gravitational-wave candidate event database (GraCEDb) is a prototype 
system to organize candidate events from gravitational-wave searches and 
to provide an environment to record information about follow-ups. A simple 
client tool is provided to submit a candidate event to the database.
 
%prep
%setup -n %{name}-%{unmangled_version}
 
%build
python setup.py build
 
%install
python setup.py install --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES
 
%clean
rm -rf $RPM_BUILD_ROOT
 
%files -f INSTALLED_FILES
%defattr(-,root,root)
%exclude %{python_sitelib}/ligo/gracedb/*pyo
%exclude %{python_sitelib}/ligo/gracedb/test/*pyo
