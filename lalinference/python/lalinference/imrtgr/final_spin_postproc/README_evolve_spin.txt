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

