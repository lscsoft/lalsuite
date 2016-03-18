As it takes days to evolve large number of spin posterior samples, we can split the big posterior_sample.dat file into small files 
and evolve them parallely using condor. First, `split_posterior_samples.py` splits the big file into small files. Then `run_condor_output_evolved_spins.py` takes each of these sud-posterior-sample files and evolve them upto some v_final. After the output file corresponding to each input sub file is created, those are combined using the code `combine_evolved_posterior_samples.py`.

-------------------------------------
Usage of `split_posterior_samples.py`

python split_posterior_samples.py <Big_posterior_sample_file> <tag name for the small files> <directory path where the small files will be created> <No of files>
-------------------------------------

---------------------------------------------
Usage of `run_condor_output_evolved_spins.py`

python run_condor_output_evolved_spins.py <tag for the sub file> <output directory> <Node number>
----------------------------------------------


------------------------------------------------
Usage of `combine_evolved_posterior_samples.py`

python combine_evolved_posterior_samples.py <tag name of the small files> <directory path where the big file will be created> <output file name> <number of files to be combined>
----------------------------------------------
