Description of the NR final mass and spin comparison script. This script compares the predictions of final mass and spin fits to the values from the metadata of the publicly
available SXS runs. It plots histograms of the fractional and absolute errors and also plots these errors versus the mass ratio, total spins, in-plane spin components, and total in-plane spin, normalized to the total mass

---------------

Main code: calc_final_mass_spin_error_NKJ-M_combined_hist.py

Data file: SXS_catalog_params.txt

---------------

Usage:

Edit the code and select the appropriate fit_tag for the fit you want to test. Set evolve_spins to 0 or 1 depending on whether or not you want to evolve the spins (1 for evolution,
0 for no evolution). If you are evolving the spins and do not want to evolve to the ISCO frequency, change the value of v_final appropriately. Set v0L to 0 or 1 depending on
whether you want to calculate the initial velocity for the spin evolution from M omega0 or from the magnitude of the orbital angular momentum (using the consistent 1PN expression).

Then run

python calc_final_mass_spin_error_NKJ-M_combined_hist.py
