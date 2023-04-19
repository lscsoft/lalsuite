Parameters inputs & units
==========================

This page should provide a user with an overview on how :code:`GWSignal` manages the various parameters entering the waveform. 

Parameters dictionary
---------------------
Every parameter should always have :code:`astropy` units. They can be added manually as in the following::

	import astropy.units as u

	m1 = 60.*u.solMass
	m2 = 30.*u.solMass
	s1x = 0.*u.dimensionless_unscaled
	s1y = 0.*u.dimensionless_unscaled
	s1z = 0.*u.dimensionless_unscaled
	s2x = 0.*u.dimensionless_unscaled
	s2y = 0.*u.dimensionless_unscaled
	s2z = 0.*u.dimensionless_unscaled

	deltaT = 1./1024.*u.s
	f_min = 20.*u.Hz
	f_ref = 20.*u.Hz
	distance = 1000.*u.Mpc
	inclination = 0.*u.rad

	phiRef = 0.*u.rad
	eccentricity = 0.*u.dimensionless_unscaled
	longAscNodes = 0.*u.rad
	meanPerAno = 0.*u.rad

	approximant = 'IMRPhenomPv2_NRTidalv2'

	lambda1 = 500*u.dimensionless_unscaled
	lambda2 = 100*u.dimensionless_unscaled

	python_dict = {'mass1' : m1,
	               'mass2' : m2,
	               'spin1x' : s1x,
	               'spin1y' : s1y,
	               'spin1z' : s1z,
	               'spin2x' : s2x,
	               'spin2y' : s2y,
	               'spin2z' : s2z,
	               'deltaT' : deltaT,
	               'f22_start' : f_min,
	               'f22_ref': f_ref,
	               'phi_ref' : phiRef,
	               'distance' : distance,
	               'inclination' : inclination,
	               'eccentricity' : eccentricity,
	               'longAscNodes' : longAscNodes,
	               'meanPerAno' : meanPerAno,
	               'lambda1':lambda1,
	               'lambda2':lambda2,
	               'condition':0
	               }



The units to a python dictionary can be set automatically while using the function :code:`lalsimulation.gwsignal.core.utils.add_params_units`. The :code:`units_sys` can be set 
either to 'S.I.' or 'Cosmo'. For example::

	import lalsimulation 
	m1 = 60.
	m2 = 30.
	s1x = 0.
	s1y = 0.
	s1z = 0.
	s2x = 0.
	s2y = 0.
	s2z = 0.

	deltaT = 1./1024.
	f_min = 20.
	f_ref = 20.
	distance = 1000.
	inclination = 0.

	phiRef = 0.*u.rad
	eccentricity = 0.
	longAscNodes = 0.*u.rad
	meanPerAno = 0.*u.rad

	approximant = 'IMRPhenomPv2_NRTidalv2'

	lambda1 = 500
	lambda2 = 100

	python_dict = {'mass1' : m1,
	               'mass2' : m2,
	               'spin1x' : s1x,
	               'spin1y' : s1y,
	               'spin1z' : s1z,
	               'spin2x' : s2x,
	               'spin2y' : s2y,
	               'spin2z' : s2z,
	               'deltaT' : deltaT,
	               'f22_start' : f_min,
	               'f22_ref': f_ref,
	               'phi_ref' : phiRef,
	               'distance' : distance,
	               'inclination' : inclination,
	               'eccentricity' : eccentricity,
	               'longAscNodes' : longAscNodes,
	               'meanPerAno' : meanPerAno,
	               'lambda1':lambda1,
	               'lambda2':lambda2,
	               'condition':0
	               }
	python_dict = lalsimulation.gwsignal.core.utils.add_params_units(python_dict, units_sys='Cosmo')


Extra parameters
------------------

In the case where your waveform require the usage of extra non-standard parameters in gravitational waves (refer to :ref:`Reference to standard GW parameters` for a the standard gw parameters as implemented in :code:`GWSignal`), these should be added as part of the metadata of your waveform generator as a dictionary in the following way::

	gen.metadata["extra_parameters"] = {'param_1':u.dimensionless_unscaled,'param_2':u.solMass}
	
Then the parameters enter the waveform dictionary as the standard ones::

	python_dict = {'param_1':2*u.dimensionless_unscaled,
	               'param_2':30*u.solMass}
