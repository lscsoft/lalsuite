The code `output_evolved_spins.py` evolves the posterior samples of spins resulting from parameter estimation runs. It reads in the posterior samples of masses (m1_source, m2_source), spin magnitudes (a1, a2) and orientation w.r.t. to orbital angular momentum LN (tilt1, tilt2 and phi12) at some initial frequency `v0`. Here, tilt1 and tilt2 are zenith angles of spins S1 and S2 while phi12 is the difference of azimuthal angles of S1 and S2 on orbital plane. The code outputs the posterior samples of tilt1, tilt2 and phi12 at some final frequency `v_final`. It also prints the initial and final samples to the stdout.

It imports a module `pneqs.py` for the evolution of posterior samples. The function `find_tilts_and_phi12_at_freq`, which takes v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final and dt as inputs, returns the values of tilt1, tilt2 and phi12 at v_final. The ACTUAL solving of the set of differential equations is done by the function `evolve_spins_dt` which takes same inputs as `find_tilts_and_phi12_at_freq`. The spin and orbital angular momentum vectors are evolved from a starting frequency `v0` to a final frequency `v_final` in the time steps of `dt`. 

The set of differential equations are provided by the function `precession_eqns` which contains precessional equations for S1, S2 and LN vectors as well as a differential equation for `v`. The precessional equations for S1, S2 and LN are given by Eqs.(3.6)-(3.8) in Ajith 2011 (Phys. Rev. D 84, 084037). The differential equation for v is given by Eqs.(3.2) and (3.5) of the same paper. The expression for (dE/dv)/Flux in Eq.(3.5) is coded up in file `dE_by_Flux.py`.

----------------------------------------------------
Description of how to evolve spin posterior samples
----------------------------------------------------
output_evolved_spins.py: main code

pneqns.py: contains functions that have been used in the main code.

dE_by_Flux.py: called by pneqns.py and contains a function `denergy_by_flux` to compute re-exapnded dEnergy/flux 


----------------------------------
Usage of output_evolved_spins.py
----------------------------------
python output_evolved_spins.py <posterior_samples.dat file> <tag for output> <location to write output>

Here there are three required command-line inputs:
1. The posterior_samples.dat file to run on
2. The tag to assign to the output
3. The location to write the output


