from numpy import loadtxt, logaddexp, log, mean
from optparse import OptionParser
import gzip
import os
import os.path
import sys
from functools import reduce

from lalinference import LALInferenceHDF5PosteriorSamplesDatasetName as posterior_dset_name
from lalinference import LALInferenceHDF5NestedSamplesDatasetName as nested_dset_name
from lalinference.nest2pos import draw_posterior_many, draw_N_posterior_many, compute_weights

usage = '''%prog [-N Nlive] [-p posterior.hdf5] [-H header.txt] [--npos Npos] [--non-strict-versions] datafile1.hdf5 [datafile2.hdf5 ...]

%prog takes at least one nested sampling output file and outputs posterior
\tsamples. If more than one input file is specified, each file is converted,
\tthen posterior samples drawn according to the evidence of each.
\tIf the --npos option is used the algorithm
\twill draw approximately that number of samples from the posterior. This may
\tgive repeated samples in the output file. By default, the non-repeating
\talgorithm is used, but that may not produce enough samples.
\tThe input and output files may be in either HDF5 or ASCII format, with
\tASCII tables being deprecated. The type will be chosen based on the file extensions.
'''


def read_nested_from_hdf5(nested_path_list, strict_versions=True):
    headers = None
    input_arrays = []
    metadata = {}
    log_noise_evidences = []
    log_max_likelihoods = []
    nlive = []
    from lalinference.io import read_samples, extract_metadata

    for path in nested_path_list:
        if not os.path.isfile(path):
                print('Unable to open %s, skipping file'%(path))
                continue
        try:
                tab = read_samples(path,tablename=nested_dset_name)
                input_arrays.append(tab)
        except KeyError:
                print('Unable to read table from %s, skipping'%(path))
                continue

        # N.B.: This appends to metadata, log_noise_evidences, log_max_likelihoods, nlive, in addition to outputting run_identifier
        run_identifier = extract_metadata(path, metadata, log_noise_evidences, log_max_likelihoods, nlive, nested_dset_name, True, strict_versions)

    if len(input_arrays)==0:
        print('No nested samples could be read from %s'.format(str(nested_path_list)))
        raise IOError

    # for metadata which is in a list, take the average.
    for level in metadata:
        for key in metadata[level]:
            if isinstance(metadata[level][key], list) and all(isinstance(x, (int,float)) for x in metadata[level][key]):
                metadata[level][key] = mean(metadata[level][key])
            elif isinstance(metadata[level][key], list) and all(isinstance(x, (str)) for x in metadata[level][key]):
                print("Warning: only printing the first of the %d entries found for metadata %s/%s. You can find the whole list in the headers of individual hdf5 output files\n"%(len(metadata[level][key]),level,key))
                metadata[level][key] = metadata[level][key][0]
    log_noise_evidence = reduce(logaddexp, log_noise_evidences) - log(len(log_noise_evidences))
    log_max_likelihood = max(log_max_likelihoods)

    return input_arrays, log_noise_evidence, log_max_likelihood, metadata, nlive, run_identifier


def read_nested_from_ascii(nested_path_list):
    print('Warning: ASCII files are deprecated in favor of HDF5.')
    default_header_path = nested_path_list[0]+'_params.txt'
    if opts.headers is not None:
        header_path = opts.headers
        if not os.access(header_path, os.R_OK):
            raise OSError('Unable to open header file: %r' % header_path)
    elif os.access(default_header_path, os.R_OK):
        header_path = default_header_path
    else:
        header_path = None
    if header_path is None:
        if opts.verbose:
            print('No header file found, assuming inspnest default')
        header = 'mchirp eta time phi0 dist ra dec psi iota logL'
    else:
        if opts.verbose:
            print('Reading headers from %r' % header_path)
        with open(header_path, 'r') as f:
            header = f.readline()
    headers = header.split()
    input_arrays = map(loadtxt, nested_path_list)

    log_noise_evidences = []
    log_max_likelihoods = []
    for path in nested_path_list:
        path = path.replace('.gz',  '')+'_B.txt'
        if os.access(path, os.R_OK):
            content = loadtxt(path)
            log_noise_evidences.append(content[2])
            log_max_likelihoods.append(content[3])
    if log_noise_evidences:
        log_noise_evidence = reduce(logaddexp, log_noise_evidences)
        log_max_likelihood = max(log_max_likelihoods)
    else:
        log_noise_evidence = 0
        log_max_likelihood = 0

    return headers, input_arrays, log_noise_evidence, log_max_likelihood


def write_posterior_to_hdf(posterior_path, headers, posterior, metadata,
                           run_identifier):
  from lalinference.io import write_samples
  write_samples(posterior, posterior_path, path='/'.join(['','lalinference',run_identifier,posterior_dset_name]), metadata=metadata, overwrite=True)

def write_posterior_to_ascii(posterior_path, headers, posterior,
                             log_bayes_factor, log_evidence,
                             log_noise_evidence, log_max_likelihood):
    print('Warning: HDF5 output is inactive since %r was passed as '
          'output file. ASCII output is deprecated.' % opts.pos)

    # set the output number format
    if 1 <= opts.prec <= 20:
        prec = opts.prec
    else:
        prec = 14
        print('Using default precision instead of %r.' % opts.prec)

    if posterior_path is None:
        posterior_file = sys.stdout
    else:
        if opts.gz:
            posterior_file = gzip.open(posterior_path, 'wb')
        else:
            posterior_file = open(posterior_path, 'w')
    for field in headers:
        posterior_file.write('%s\t' % field)
    posterior_file.write('\n')

    for row in posterior:
        for value in row:
            posterior_file.write('{1:.{0}e}\t'.format(prec, value))
        posterior_file.write('\n')
    if opts.pos is not None:
        posterior_file.close()

    # Write out an evidence file
    if opts.pos is not None:
        evidence_file_path = opts.pos + '_B.txt'
        with open(evidence_file_path, 'w') as out:
            out.write('%f %f %f %f\n' % (log_bayes_factor, log_evidence,
                                         log_noise_evidence, log_max_likelihood))
    # Write out an header (git info and command line) file
    if opts.pos is not None:
        strout = ''
        for dfile in datafiles:
            fin = '%s_header.txt' % dfile
            if os.path.isfile(fin):
                with open(fin, 'r') as f:
                    for l in f.readlines():
                        strout += l
                strout += '\n\n'
        if strout != '':
            with open(opts.pos + '_header.txt', 'w') as fout:
                fout.write(strout)


def is_hdf5(path):
    extension = os.path.splitext(path)[1]
    if '.h5' in extension or '.hdf' in extension:
        return True
    elif extension.strip().lower() == '.dat':
        return False
    else:
        raise ValueError("Extension %r not recognized as HDF5 or '.dat': %r"
                         % (extension, path))


if __name__ == '__main__':
    parser = OptionParser(usage)
    parser.add_option(
        '-N', '--Nlive', action='store', type='int', dest='Nlive', default=None,
        help='Number of live points in each chain loaded', metavar='NUM')
    parser.add_option(
        '-p', '--pos', action='store', type='string', default=None,
        help='Output file for posterior samples', metavar='posterior.dat')
    parser.add_option(
        '--npos', action='store', type='int', default=None,
        help='Draw a specific number of posteriors samples. May give '
             'repetitions. Disabled by default.', metavar='NUM')
    parser.add_option(
        '-H', '--headers', action='store', type='string', default=None,
        help='Header file explaining columns in data file',
        metavar='file.dat_params.txt')
    parser.add_option(
        '-d', '--prec', action='store', type='int', dest='prec', default=14,
        help='Number of decimal place required for output posterior samples. '
             'Default is 14.', metavar='NUM')
    parser.add_option(
        '-z', '--gzip', action="store_true", dest='gz', default=False,
        help='Gzip compress the output posterior samples (this will append '
             '.gz onto the posterior file). Default is no compression.')
    parser.add_option(
        '-v', '--verbose', action='store_true', default=False,
        help='Print some additional information')
    parser.add_option(
        '--non-strict-versions', action='store_true', dest='nonstrict',
        default=False, help='Set this flag to allow combining of nested '
                            'sample HDF5 files that do no necessarily have '
                            'the same verion information')
    opts, args = parser.parse_args()

    # Argument checking
    datafiles = args
    if len(datafiles) < 1:
        print('No input file specified, exiting')
        sys.exit(1)

    if all(map(is_hdf5, datafiles)):
        hdf_input = True
    elif not any(map(is_hdf5, datafiles)):
        hdf_input = False
    else:
        raise ValueError('Input files appear to be mixed between HDF5 and ASCII.')

    if opts.pos is None:
        hdf_output = False
        print('No output file given, writing to stdout.')
    else:
        hdf_output = is_hdf5(opts.pos)

    if opts.Nlive is None and not hdf_input:
        print('Must specify number of live points using the --Nlive option if '
              'ASCII tables are used.')
        sys.exit(1)
    elif hdf_input and opts.Nlive is not None:
        print('Nlive is ignored in favor of the value in the HDF metadata.')

    # Read the input file
    return_values = read_nested_from_hdf5(datafiles,
                                          strict_versions=(not opts.nonstrict))
    (input_arrays, log_noise_evidence, log_max_likelihood,
       metadata, nlive, run_identifier) = return_values
    if len(input_arrays)==0:
      print('Error: No input file were read')
      sys.exit(1)
    headers = input_arrays[0].dtype.names
    nlive = list(map(int, nlive))

    if opts.npos is not None:
        def sampler(datas, Nlives, **kwargs):
            return draw_N_posterior_many(datas, Nlives, opts.npos, **kwargs)
    else:
        sampler = draw_posterior_many

    # Create the posterior from nested samples
    # inarrays has shape (nfiles, nsamples. nfields)
    posterior = sampler(input_arrays,
                        nlive,
                        verbose=opts.verbose)
    # posterior is a list of lists/array with overall shape (nsamples, nfields)
    # posterior = array(posterior)

    log_evs, log_wts = zip(*[compute_weights(data['logL'], n)
                             for data, n in list(zip(input_arrays, nlive))])
    if opts.verbose:
        print('Log evidences from input files: %s' % str(log_evs))

    # Compute the evidence
    log_evidence = reduce(logaddexp, log_evs) - log(len(log_evs))
    log_bayes_factor = log_evidence - log_noise_evidence

    # Update the metadata with the new evidences
    run_level = '/lalinference/'+run_identifier
    if run_level not in metadata:
        metadata[run_level] = {}
    metadata[run_level]['log_bayes_factor'] = log_bayes_factor
    metadata[run_level]['log_evidence'] = log_evidence
    metadata[run_level]['log_noise_evidence'] = log_noise_evidence
    metadata[run_level]['log_max_likelihood'] = log_max_likelihood
    if hdf_output:
        write_posterior_to_hdf(opts.pos, headers, posterior, metadata,
                               run_identifier)
    else:
        write_posterior_to_ascii(opts.pos, headers, posterior,
                                 log_bayes_factor, log_evidence,
                                 log_noise_evidence, log_max_likelihood)
