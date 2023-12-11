# LALSuite

This is the main LALSuite development repository.
If you would like to just use a released version,
see [here](https://computing.docs.ligo.org/conda/)
for IGWN conda environments including LALSuite,
or the project pages
on [conda-forge](https://anaconda.org/conda-forge/lalsuite)
and on [PyPI](https://pypi.org/project/lalsuite/).

## Acknowledgment

We request that any academic report, publication, or other academic
disclosure of results derived from the use of this software acknowledge
the use of the software by an appropriate acknowledgment or citation.

The whole LALSuite software suite can be cited with the DOI
[10.7935/GT1W-FZ16][doi]. An example BibTeX entry could look like this:

     @misc{lalsuite,
           author         = "{LIGO Scientific Collaboration} and {Virgo Collaboration} and {KAGRA Collaboration}",
           title          = "{LVK} {A}lgorithm {L}ibrary - {LALS}uite",
           howpublished   = "Free software (GPL)",
           doi            = "10.7935/GT1W-FZ16",
           year           = "2018"
     }

The [Python/Octave interfaces to LALSuite][swiglal] are described in this paper:

     @article{swiglal,
              title     = "{SWIGLAL: Python and Octave interfaces to the LALSuite gravitational-wave data analysis libraries}",
              author    = "Karl Wette",
              journal   = "SoftwareX",
              volume    = "12",
              pages     = "100634",
              year      = "2020",
              doi       = "10.1016/j.softx.2020.100634"
     }

In addition, some codes contained in this package may be directly based
on one or several scientific papers, which should be cited when using
those specific codes; some of these can be discovered through the
documentation.

## Cloning the Repository

We now utilize [Git LFS][gitlfs] for the managament of large files and
as such `git-lfs` needs to be installed and configured to correctly
clone this repository. After installing `git-lfs` it can be configured
using:

     $ git lfs install

This only needs to be done once for each machine you access the
repository. It can then be cloned using:

     $ git clone git@git.ligo.org:lscsoft/lalsuite.git

## Building from Source

The recommended way to build LALSuite from source is in a `conda` environment.
[A recipe file](common/conda/environment.yml) is available with all main dependencies.
This can serve as the base for custom recipes,
or be used directly via:

     $ conda env create -f common/conda/environment.yml

Pulling in dependencies may take a while depending on your internet connection.
After the environment setup succeeded, you can activate it with:

     $ conda activate lalsuite-dev

The first time building LALSuite from source, you need to set up the build
system by running (from the top level of the Git repository):

     $ ./00boot

The next step is to configure the LALSuite build system, which determines what
components of LALSuite are built, and with what features. This is done by
running (from the top level of the Git repository):

     $ ./configure [options...]

A full list of options is available by running `./configure --help=recursive`.
Some commonly-used options are:

- `--prefix=<path>`: Sets which `<path>` LALSuite will be installed into. By
  default, LALSuite is installed into a directory `./_inst/` at the top level of
  the Git repository.

- `--enable-all-lal`: Build all component libraries of LALSuite, including
  `lalapps`; the default.

- `--disable-all-lal`: Build only the `lal` library component of LALSuite.

- `--enable-lal<name>`: Build the `lal<name>` component of LALSuite. This option
  can be combined with the `--disable-all-lal` option to selectively build
  LALSuite components, e.g.  `--disable-all-lal --enable-lalframe
  --enable-lalinference` will build only the `lal`, `lalframe`, and
  `lalinference` components.

- `--enable-swig-python`, `--enable-swig-octave`: Build Python and/or Octave
  interfaces to the LALSuite libraries (using the SWIG tool), which allow
  LALSuite library functions to be called directly from Python/Octave. By
  default, only the Python interface is built. Use `--enable-swig` to build both
  interfaces.

- `--enable-fftw`, `--enable-intelfft`: Perform Fast Fourier Transforms using
  either the FFTW library (the default) or the Intel MKL library. To enable
  aligned memory optimisations with FFTW, use the `--enable-fftw3-memalign`
  option.

- `--enable-framel`, `--enable-framec`: Read/write gravitational wave frame
  files using either the FrameL or FrameCPPC libraries. By default, the LALSuite
  build system will use whichever library is available on your system; if
  neither are available, frame file reading/writing will be disabled.

- `--with-hdf5`: Build support for reading/writing to HDF5 files.

- `--with-cfitsio`: Build support for reading/writing to FITS files.

- `--with-cuda=<path>`: Build CUDA code which runs on graphics processing
  units. `<path>` should specify the top-level directory of a local install of
  the CUDA Toolkit.

- `PYTHON=<pythonX.Y>`: set the Python interpreter used by LALSuite to run
  Python scripts and compile Python modules. By default the LALSuite build
  system will find the default system Python interpreter, so this is only needed
  if you want to build against a specific Python version `X.Y`.

Once LALSuite has been configured, you can then build and install LALSuite by
running (from the top level of the Git repository):

     $ make [-j<processes>] && make install

where `<processes>` optionally specifies the number of processes used for
building in parallel. If the build is successful, `make install` will install
LALSuite into the directory given by `--prefix`.

After pulling updates or making your own changes, you will usually only need to
call `make && make install` again, as reconfiguration and re-running `00boot`
should be handled automatically if needed.

If you prefer managing dependencies yourself without conda,
see [here](https://wiki.ligo.org/Computing/LALSuite#Dependencies).

## Contributing to LALSuite

The [guide to contributing to LALSuite][contributing] explains how to
report issues and contribute fixes or new features using the fork and
pull workflow. Please read and follow these directions.

## Nightly Documentation

The [LALSuite Doxygen documentation][nightlydocs] is built under
GitLab-CI every night.

## Notes on Ancient History

LALSuite was transferred to `git.ligo.org` in December 2017. Older
history has been imported, though commit hashes were rewritten during
the [Git LFS][gitlfs] conversion. Please note:

1. The `Original:` commit IDs quoted in each commit message can be used
   to compare with the [archived reference repo][oldlalsuite], old issue
   discussions on the [Redmine tracker][oldredmine], review wiki pages
   etc.

1. Commits before December 2017 may also include references to issues
   (`#number`). These refer to the corresponding [Redmine
   issue][oldredmine] (LVC-authorized access only), and any clickable
   link the internal GitLab web interface produces for those old commits
   will therefore be spurious.

## For More Information

Please visit the [LALSuite project page][project].

[doi]:          https://doi.org/10.7935/GT1W-FZ16
[swiglal]:      https://lscsoft.docs.ligo.org/lalsuite/lalsuite/swiglal_tutorial.html
[gitlfs]:       https://wiki.ligo.org/Computing/GitLFS#Install_the_git_LFS_client
[contributing]: https://lscsoft.docs.ligo.org/lalsuite/lalsuite/lalsuite_contributing.html
[nightlydocs]:  https://lscsoft.docs.ligo.org/lalsuite
[oldlalsuite]:  https://git.ligo.org/lscsoft/lalsuite-archive
[oldredmine]:   https://bugs.ligo.org/redmine/projects/lalsuite
[project]:      https://wiki.ligo.org/Computing/LALSuite
