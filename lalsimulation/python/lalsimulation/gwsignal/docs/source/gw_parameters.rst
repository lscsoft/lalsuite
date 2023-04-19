Reference to standard GW parameters
====================================

In this section we walk through the parameters that are considered to be standard and their units as implemented in :code:`GWSignal`. 

Masses & spins
------------------

Masses should be specified using any parameters of the following list::

		mass_parameters = ["mass1", "mass2","total_mass","chirp_mass","mass_difference", 
		                   "reduced_mass","mass_ratio","sym_mass_ratio"]

They need to be such that the tuple :math:`(m_1,m_2)` is unambiguously specified. Those parameters that are dimensionless should be given the :code:`u.dimensionless_unscaled` unit to be accepted by :code:`GWSignal`. Spins can also be used in cartesian or spherical coordinates with names::

		spin_parameters = ["spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z",
		                   "spin1_norm","spin1_tilt","spin1_phi","spin2_norm","spin2_tilt","spin2_phi"]

All of them being dimensionless parameters but the spin angles which are :code:`u.rad`.


Extrinsic parameters
--------------------------------

Those parameters extrinsic to the binary::

	extrinsic_params = ["distance", "inclination",  "longAscNodes", "meanPerAno"]


 
Waveform generation parameters
--------------------------------

Parameters that completely specify how the waveform is generated::

	gen_params = ["deltaT", "deltaF", "f22_start", "f_max", "phi_ref", "f22_ref"]

where any frequency should be given :code:`u.Hz` and any time :code:`u.s` or :code:`astropy` units compatible with those.


Other parameters
--------------------------------

Other parameters include :code:`eccentricity` but also array-like parameters::

	arr_params = ["ModeArray","ModeArrayJframe"]

and the parameter condition, that switches he conditioning routines on/off, and which can take the values :math:`1` or :math:`0`. 

Tidal parameters
--------------------------------

Tidal parameters accepted by waveforms that take into account matter effects::

		tidal_params = ["lambda1","lambda2","TidalOctupolarLambda1","TidalOctupolarLambda2",
		                "TidalHexadecapolarLambda1","TidalHexadecapolarLambda2",
		                "TidalQuadrupolarFMode1","TidalQuadrupolarFMode2",
		                "TidalOctupolarFMode1","TidalOctupolarFMode2"]

All of them are dimensionless.


Non GR parameters
--------------------------------

Standard parameters used in testing gr and related fields::

		nongr_params = ["phi1","phi2","phi3","phi4","dchi0","dchi1","dchi2","dchi3","dchi4",
		                "dchi5","dchi5l","dchi6","dchi6l","dchi7","dxi1","dxi2","dxi3","dxi4",
		                "dxi5","dxi6","dsigma1","dsigma2","dsigma3","dsigma4","dalpha1",
		                "dalpha2","dalpha3","dalpha4","dalpha5","dbeta1","dbeta2","dbeta3",
		                "alphaPPE","betaPPE","alphaPPE0","betaPPE0","alphaPPE1","betaPPE1",
		                "alphaPPE2","betaPPE2","alphaPPE3","betaPPE3","alphaPPE4","betaPPE4",
		                "alphaPPE5","betaPPE5","alphaPPE6","betaPPE6","alphaPPE7","betaPPE7",
		                "liv","log10lambda_eff","LIV_A_sign","nonGR_alpha"]

You can check the definitions `here <https://arxiv.org/abs/1703.01076>`__ 