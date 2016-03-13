Description of NKJ-M's postprocessing script and associated files:

calc_Mf_af_post.py: Main script

Mfaf_posterior.py: Calculates the posteriors for the final mass and spin using various fits

nr_fits_Mf_af.py: Contains the numerical relativity fits used

aux.py: Contains some of Archisman's functions that are used to compute confidence intervals

-----

Usage of calc_Mf_af_post.py:

python calc_Mf_af_post.py <posterior_samples.dat file> <tag for output> <location to write output (default: current directory)> [--outputsamples]

Here there are three required command-line inputs:

1. The posterior_samples.dat file to run on
2. The tag to assign to the output
3. The location to write the output (the current directory, if not specified)

There is also one optional flag:

--outputsamples: Output samples for the final mass and spin using the most accurate fit we have found to date, the HLZ augmented with in-plane spins using M (HLZinplaneIMRPhenomP).

-----

Notes: The "BR" and "BR aligned" tag results for the final mass and spin use the HLZ fit for the final spin--this is just done to simplify the coding, since the HLZ final mass fit
is quite accurate even for precessing spins, so we have not coded up the BR final mass fit.