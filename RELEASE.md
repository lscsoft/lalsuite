# How to create a LALSuite release

The document attempts to describe the workflow for creating a release of
a subpackage of LALSuite.
Note that only authorised persons will be able to push directly to the
relevant release branch.

The workflow is (roughly) as follows:

1. cherry-pick/merge everything required onto the release branch
   (e.g. `o3b-release`) via a Merge Request; please use the `-x` and `--gpg-sign`
   option to `git-cherry-pick`
1. for each subpackage (e.g. `lal`):
   1. update the following for the release:
      - `/lal/configure.ac` (package and library version numbers, and
        upstream LALSuite dependencies -- carefully read the instructions on
        how to change the library version number)
      - `/lal/debian/control.in` (changelog entry)
      - `/lal/lalsimulation.spec.in` (changelog entry)
   1. commit the resulting changes onto the release branch
      (one file at a time, if you like)
   1. tag `lal-vX.Y.Z` on the release branch
      ```console
      git tag --annotate --sign lal-vX.Y.Z
      ```
   1. generate the release tarball
      ```console
      git checkout lal-vX.Y.Z
      ./00boot
      ./configure --enable-swig
      make dist -C lal
      git checkout <release-branch>
      ```
1. update the following for the top-level `lalsuite` release:
   - `/configure.ac` (lalsuite version number)
1. commit the resulting changes onto the release branch
1. tag `lalsuite-vX.Y` on the release branch
   ```
   git tag --annotate --sign lalsuite-vX.Y
   ```
1. push the new references up to the repository
   ```console
   git push upstream lal-vX.Y.Z  # (repeat as required for each subpackage)
   git push upstream lalsuite-vX.Y
   ```
1. upload the release tarballs to <http://software.ligo.org/lscsoft/source/>
1. [open an SCCB ticket](https://sccb.docs.ligo.org/requests/)
1. merge the release branch onto `master` via a Merge Request
   - use a local branch on your fork to merge the two branches
   - update the package versions in `/configure.ac` and all
     `/<subpackage>/configure.ac` files for released subpackages to
     `X.Y.Z.1` to indicate that they are back in development mode
