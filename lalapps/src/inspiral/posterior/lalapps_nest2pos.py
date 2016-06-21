from numpy import loadtxt, logaddexp, log, array
from optparse import OptionParser
import sys
import os
import gzip
import h5py

from lalinference import LALInferenceHDF5PosteriorSamplesGroupName as posterior_grp_name
from lalinference import LALInferenceHDF5NestedSamplesGroupName as nested_grp_name
from lalapps.nest2pos import draw_posterior_many, draw_N_posterior_many, compute_weights

usage = '''%prog -N Nlive [-p posterior.dat] [-H header.txt] [--npos Npos] datafile1.dat [datafile2.dat ...]

%prog takes at least one nested sampling output file and outputs posterior samples.
\tIf more than one input file is specified, each file is converted, then posterior samples drawn
\taccording to the evidence of each. Will output to stdout if no -p option given.
\tIf the --npos option is used the algorithm will draw approximately that number of samples
\tfrom the posterior. This may give repeated samples in the output file. By default, the
\tnon-repeating algorithm is used, but that may not produce enough samples.
'''


def getBfile(Bfile):
    Bfile = Bfile
    if os.access(Bfile, os.R_OK):
        outstat = loadtxt(Bfile)
        return outstat
    else:
        return None


def read_hdf5(nested_samples_path):
    attr_dict = {}
    with h5py.File(nested_samples_path, 'r') as hdf:
        # We step through the file level by level and check uniqueness and
        # availability. We also store the metadata, i.e. the attr dict with
        # the corresponding path for each level and dataset.
        group = 'lalinference'
        assert group in hdf
        attr_dict[group] = {key: value for key, value in hdf[group].attrs.items()}

        if len(hdf[group].keys()) != 1:
            raise KeyError('Multiple run-identifiers found: %r' % list(hdf[group].keys()))
        # we ensured above that there is only one identifier in the group.
        run_identifier = list(hdf[group].keys())[0]
        group += '/' + run_identifier
        attr_dict[group] = {key: value for key, value in hdf[group].attrs.items()}

        assert 'nested_samples' in hdf[group]
        group += '/nested_samples'
        nested_samples_group = hdf[group]
        # The following line is commented out since nest_specific metadata
        # doesn't make much sense for the posterior samples.
        attr_dict[group] = {key: value for key, value in hdf[group].attrs.items()}

        input_data_dict = {}
        for key in nested_samples_group:
            input_data_dict[key] = nested_samples_group[key][...]
            attr_dict[group+'/'+key] = {key: value for key, value in hdf[group+'/'+key].attrs.items()}
    return input_data_dict, attr_dict, run_identifier


if __name__ == '__main__':
    headerfilename = None
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
    opts, args = parser.parse_args()

    # Argument checking
    datafiles = args
    if len(datafiles) < 1:
        print('No input file specified, exiting')
        sys.exit(1)

    if opts.Nlive is None:
        print('Must specify number of live points using the --Nlive option')
        sys.exit(1)

    # Determine the file format: HDF5 as the current standard or ASCII table
    # and header file as legacy support.
    # Load the data and generate a list of fields in the appropriate way.
    if all('.hdf' in dfile or '.h5' in dfile for dfile in datafiles):
        headers = None
        inarrays = []
        for dfile in datafiles:
            input_data_dict, attr_dict, run_identifier = read_hdf5(dfile)
            current_headers = sorted(input_data_dict.keys())
            if headers is None:
                headers = current_headers
            else:
                # We have the fields already, but check to avoid inconsistencies
                if headers != current_headers:
                    # fields in headers but not current_headers
                    diff = set(headers) - set(current_headers)
                    # add fields in current_headers but not headers
                    diff = diff | (set(current_headers) - set(headers))
                    raise ValueError('Data fields do not match between '
                                     'inputfiles. Unique fields: %r'
                                     % sorted(diff))
            # check that dimensions are equal for each field
            shape_set = set(dset.shape for dset in input_data_dict.values())
            assert len(shape_set) == 1, repr(shape_set)
            # convert to numpy array and store in the 'inarrays' list
            # loadtxt was of shape (nsamples, nfields), so we transpose
            inarrays.append(array([input_data_dict[key] for key in headers]).T)

    elif all(dfile.endswith('.dat') for dfile in datafiles):
        print('Warning: ASCII files are deprecated in favor of HDF5.')

        if opts.headers is None:
            defaultparamfilename = args[0]+'_params.txt'
            if os.access(defaultparamfilename, os.R_OK):
                headerfilename = defaultparamfilename
        else:
            headerfilename = opts.headers
        if headerfilename is not None:
            if not os.access(headerfilename, os.R_OK):
                print('Unable to open header file %s!' % headerfilename)
                sys.exit(1)
            elif opts.verbose:
                print('Reading headers from %s' % headerfilename)
            headerfile = open(headerfilename, 'r')
            headerstr = headerfile.readline()
            headerfile.close()
        else:  # old inspnest defaults
            if opts.verbose:
                print('No header file found, assuming inspnest default')
            headerstr = 'mchirp\t eta\t time\t phi0\t dist\t ra\t dec\t psi\t iota\t logL'
        headers = headerstr.lower().split()

        # Create posterior samples for each input file
        inarrays = map(loadtxt, datafiles)
    else:
        raise ValueError('Input files appear to be mixed between hdf5 and ASCII.')

    if opts.npos is not None:
        def sampler(datas, Nlives, logLcols, **kwargs):
            return draw_N_posterior_many(datas, Nlives, opts.npos,
                                         logLcols=logLcols, **kwargs)
    else:
        sampler = draw_posterior_many

    #   Preperations
    if 'logl' not in [field.lower() for field in headers]:
        raise KeyError("Key 'logl' not found. Available: %r" % headers)
    logLcol = [field.lower() for field in headers].index('logl')

    # set the output number format
    if 1 <= opts.prec <= 20:
        prec = opts.prec
        # precision is outside a reasonable range, so just set to default of 14
    else:
        prec = 14
    outformat = '%%.%de' % prec

    # Create the posterior from nested samples
    # inarrays has shape (nfiles, nsamples. nfields)
    posterior = sampler(inarrays,
                        [int(opts.Nlive) for d in datafiles],
                        logLcols=[logLcol for d in datafiles],
                        verbose=opts.verbose)
    # posterior is a list of lists/array with overall shape (nsamples, nfields)
    posterior = array(posterior)

    log_evs, log_wts = zip(*[compute_weights(data[:, logLcol], opts.Nlive) for data in inarrays])
    if opts.verbose:
        print('Log evidences from input files: %s' % str(log_evs))

    # Compute the evidence
    if all('.hdf' in dfile or '.h5' in dfile for dfile in datafiles):
        raise NotImplementedError
    elif all(dfile.endswith('.dat') for dfile in datafiles):
        # read evidence from _B files
        Bs = map(getBfile, [d.replace('.gz',  '')+'_B.txt' for d in datafiles])
        Bs = [b for b in Bs if b is not None]
        meanZ = reduce(logaddexp, log_evs)-log(len(log_evs))
        if len(Bs) > 0:  # If there were files to load, use them
            noiseZ = reduce(logaddexp, [b[2] for b in Bs]) - log(len(Bs))
            maxL = max([b[3] for b in Bs])
        else:  # Otherwise fill in some placeholder data
            noiseZ = 0
            maxL = 0
        meanB = meanZ-noiseZ



    #   Write the posterior file
    if opts.pos is None or not ('.h5' in opts.pos or '.hdf' in opts.pos):
        print('Warning: HDF5 output is inactive since %r was passed as output '
              'file.' % opts.pos)
        if opts.pos is not None:
            if opts.gz:
                posfile = gzip.open(opts.pos+'.gz', 'wb')
            else:
                posfile = open(opts.pos, 'w')
        else:
            posfile = sys.stdout
        for h in headers:
            posfile.write('%s\t' % h)
        posfile.write('\n')
        for row in posterior:
            for i in row:
                outval = outformat % i
                posfile.write('%s\t' % outval)
            posfile.write('\n')
        if opts.pos is not None:
            posfile.close()

        # Write out an evidence file
        if opts.pos is not None:
            Bfile = opts.pos+'_B.txt'
            outB = open(Bfile, 'w')
            outB.write('%f %f %f %f\n' % (meanB, meanZ, noiseZ, maxL))
            outB.close()
        # Write out an header (git info and command line) file
        if opts.pos is not None:
            strout = ""
            for i in datafiles:
                fin = "%s_header.txt" % i
                if os.path.isfile(fin):
                    f = open(fin)
                    for l in f.readlines():
                        strout += l
                    strout += '\n\n'
            if strout != "":
                fout = open(opts.pos+"_header.txt", 'w')
                fout.write(strout)
                fout.close
    else:  # HDF5
        try:
            run_identifier
        except:
            print("Input was not read from HDF5, so no run_identifier was "
                  "read. Will fall back to 'default_run_identifier'.")
            run_identifier = 'default_run_identifier'
        try:
            attr_dict
        except:
            print("Input was not read from HDF5, so no metadata was read.")
            attr_dict = {}
        if opts.pos is None:
            raise ValueError('Posterior output path required.')

        with h5py.File(opts.pos, 'w') as hdf:
            group = hdf.create_group('lalinference')
            group = group.create_group(run_identifier)
            group.attrs.update({'log_bayes_factor': meanB,
                                'log_evidence': meanZ,
                                'log_noise_evidence': noiseZ,
                                'log_max_likelihood': maxL})
            group = group.create_group('posterior_samples')
            for i, key in enumerate(headers):
                group.create_dataset(key, data=posterior[:, i], shuffle=True, compression='gzip')
            for internal_path, attributes in attr_dict.items():
                internal_path = internal_path.replace('nested_samples', 'posterior_samples')
                for key, value in attributes.items():
                    hdf[internal_path].attrs[key] = value
        # metadata, compression
