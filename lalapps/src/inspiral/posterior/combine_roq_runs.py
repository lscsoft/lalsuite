import numpy as np
import h5py as h5
import argparse

##############################################################################################################################################################

parser = argparse.ArgumentParser(description='posterior samples from individual ROQ runs')
parser.add_argument("-p", "--posteriors", default=[], type=str, nargs='+',
                    help="hdf5 file containing posterior samples")


def add_h5_link(file_name, group_name, comb_file):


	comb_file.create_group(group_name)
	comb_file[group_name+'/posterior'] = h5.ExternalLink(file_name, "/")


	nest_meta_data = comb_file[group_name+'/posterior']['lalinference']['lalinference_nest'].attrs
	log_evidence = nest_meta_data['log_evidence']
	log_prior_volume = nest_meta_data['log_prior_volume']
	
	comb_file[group_name].create_dataset("evidence_scale_factor", (1,))
	comb_file[group_name].attrs['evidence_scale_factor'] = np.exp(log_prior_volume)/np.exp(comb_file.attrs["total_log_prior_volume"])

	comb_file[group_name].attrs['evidence_volume_prior'] = log_prior_volume + log_evidence - comb_file.attrs["total_log_prior_volume"]

	comb_file[group_name].attrs['PDF_weight'] = np.exp(comb_file[group_name].attrs['evidence_volume_prior'] - comb_file.attrs["evidence_volume_prior_max"])
	print comb_file[group_name].attrs['PDF_weight']	

def combine_posterior_files(posterior_files):

	comb_file = h5.File('combined_roq_posteriors.hdf5','w')

	compute_total_prior_vol(posterior_files, comb_file)

	for (i,post_file) in enumerate(posterior_files):

		group_name = str(i)
		
		add_h5_link(post_file, group_name, comb_file) # this is where the combination actually gets done
	
	comb_file.close()	

def compute_total_prior_vol(posterior_files, comb_file):

	evp_max = -np.inf
	total_prior_volume = 0.0
	for post_file in posterior_files:

		file = h5.File(post_file, 'r')

		log_prior_volume = file['lalinference']['lalinference_nest'].attrs['log_prior_volume']
		log_evidence = file['lalinference']['lalinference_nest'].attrs['log_evidence']

		total_prior_volume += np.exp(log_prior_volume)
		evp = log_prior_volume + log_evidence
		if evp > evp_max:
			evp_max = evp
		 
		file.close()

	comb_file.attrs["total_log_prior_volume"] = np.log(total_prior_volume)
	comb_file.attrs["evidence_volume_prior_max"] = evp_max - comb_file.attrs["total_log_prior_volume"]

##############################################################################################################################################################

args = parser.parse_args()

combine_posterior_files(args.posteriors)



