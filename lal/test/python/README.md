# Testing python code

## Adding new tests

To add a new test script:

1. Create a new test script in this directory. The script should be directly executable via `python <script>` and should exit with a non-zero exit code to declare a failure.
2. Add the script to the `test_scripts` variable in the `Makefile.am` file

If your new test script requires extra packages to run, make sure and declare the new dependencies in the following places

- `/lal/debian/control.in`: `Build-Depends` (`python-` and `python3-`)
- `/lal/lal.spec.in`: `BuildRequires`

## Running the test as part of `make check`

To run the new tests you will need to rerun `00boot` to rebuild the makefiles:

```shell
./00boot
./configure --enable-swig-python --enable-python <args>
make
make check
```
