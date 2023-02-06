# Conda build for LALSuite subpackages

The `conda/` subdirectory for each lalsuite subpackage defines the
Conda package builds for that package, at least insofar as they
are executed in the Continuous Integration pipelines.

## Build files

The following files are used in the build:

| File | Purpose |
| ---- | ------- |
| `meta.yaml.in.in` | Defines the conda build, including all build/host/run dependencies |
| `build.sh` | Builds the basic C library and C-only executables, but doesn't install anything |
| `install-data.sh` | Installs the data files only (only for some libraries) |
| `install-lib.sh` | Installs the C libraries and development files only |
| `build-python.sh` | Builds the Python extensions, and installs the Python libraries |
| `install-bin.sh` | Installs the C/Python executables, and other user environment files |

Note, each sub-package may have a different layout than the general one, so
don't be too surprised if there are subtle differences in file naming or
purpose.

## Migrations

The `conda/migrations/` directory enables mimicing migrations that
have been introduce by Conda-forge in the LALSuite CI builds.
When obsolete they will be automatically ignored by `conda-build`, but
should be periodically purged from the LALSuite git project.
