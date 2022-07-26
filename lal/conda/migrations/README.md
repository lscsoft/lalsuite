# Conda migrations

This directory contains migrations that are to be applied to the
conda feedstock generaed during the CI builds of LALSuite.
These migration files should be copied from the conda-forge version
of the relevant LALSuite subpackage feedstock to enable that migration
in the CI builds while the wider conda-forge migration is ongoing.

Introducing migrations in this way should bring the CI conda builds
of LALSuite subpackages closer to the production conda-forge builds
of the same subpackages.
