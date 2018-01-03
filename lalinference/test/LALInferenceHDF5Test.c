#include <lal/XLALError.h>
#include <lal/LALInferenceHDF5.h>
#include <gsl/gsl_test.h>

int main(int argc, char **argv)
{
  /* Not used */
  (void)argc;
  (void)argv;
  XLALSetErrorHandler(XLALExitErrorHandler);

  /* Populate example variables array. */
  UINT4 N = 64;
  LALInferenceVariables **vars_array = XLALCalloc(
    N, sizeof(LALInferenceVariables *));
  for (UINT4 i = 0; i < N; i ++)
  {
    LALInferenceVariables *vars = XLALCalloc(1, sizeof(LALInferenceVariables));
    vars_array[i] = vars;
    LALInferenceAddREAL8Variable(vars, "abc", i, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddREAL8Variable(vars, "def", i, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddREAL8Variable(vars, "ghi", 5, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddREAL8Variable(vars, "ijk", i, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable (vars, "lmn", i, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddINT4Variable (vars, "opq", i, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddINT4Variable (vars, "rst", 5, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddINT4Variable (vars, "uvw", i, LALINFERENCE_PARAM_OUTPUT);
  }

  /* Open HDF5 file for writing and create group structure. */
  LALH5File *file = XLALH5FileOpen("test.hdf5", "w");
  LALH5File *group = LALInferenceH5CreateGroupStructure(
    file, "lalinference", "lalinference_mcmc");

  /* Write variables array to dataset. */
  LALInferenceH5VariablesArrayToDataset(
    group, vars_array, N, LALInferenceHDF5PosteriorSamplesDatasetName);

  /* Free variables array. */
  for (UINT4 i = 0; i < N; i ++)
  {
    LALInferenceClearVariables(vars_array[i]);
    XLALFree(vars_array[i]);
  }
  XLALFree(vars_array);

  /* Close file. */
  XLALH5FileClose(group);
  XLALH5FileClose(file);

  /* Open file for reading and find dataset. */
  file = XLALH5FileOpen("test.hdf5", "r");
  LALH5Dataset *dataset = XLALH5DatasetRead(
    file, "lalinference/lalinference_mcmc/posterior_samples");

  /* Read dataset back to variables array. */
  N = 0;
  vars_array = NULL;
  LALInferenceH5DatasetToVariablesArray(dataset, &vars_array, &N);

  gsl_test_int(N, 64, "number of rows read back");
  for (UINT4 i = 0; i < N; i ++)
  {
    LALInferenceVariables *vars = vars_array[i];
    gsl_test_int(LALInferenceGetVariableDimension(vars), 8,
      "number of columns read back");
    gsl_test_abs(LALInferenceGetREAL8Variable(vars, "abc"), i, 0,
      "value of column abc");
    gsl_test_abs(LALInferenceGetREAL8Variable(vars, "def"), i, 0,
      "value of column def");
    gsl_test_abs(LALInferenceGetREAL8Variable(vars, "ghi"), 5, 0,
      "value of column ghi");
    gsl_test_abs(LALInferenceGetREAL8Variable(vars, "ijk"), i, 0,
      "value of column ijk");
    gsl_test_int(LALInferenceGetINT4Variable (vars, "lmn"), i,
      "value of column lmn");
    gsl_test_int(LALInferenceGetINT4Variable (vars, "opq"), i,
      "value of column opq");
    gsl_test_int(LALInferenceGetINT4Variable (vars, "rst"), 5,
      "value of column rst");
    gsl_test_int(LALInferenceGetINT4Variable (vars, "uvw"), i,
      "value of column uvw");

    gsl_test_int(LALInferenceGetVariableVaryType(vars, "abc"),
      LALINFERENCE_PARAM_LINEAR, "vary type of column abc");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "def"),
      LALINFERENCE_PARAM_CIRCULAR, "vary type of column def");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "ghi"),
      LALINFERENCE_PARAM_FIXED, "vary type of column ghi");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "ijk"),
      LALINFERENCE_PARAM_OUTPUT, "vary type of column ijk");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "lmn"),
      LALINFERENCE_PARAM_LINEAR, "vary type of column lmn");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "opq"),
      LALINFERENCE_PARAM_CIRCULAR, "vary type of column opq");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "rst"),
      LALINFERENCE_PARAM_FIXED, "vary type of column rst");
    gsl_test_int(LALInferenceGetVariableVaryType(vars, "uvw"),
      LALINFERENCE_PARAM_OUTPUT, "vary type of column uvw");
  }

  /* Free variables array. */
  for (UINT4 i = 0; i < N; i ++)
  {
    LALInferenceClearVariables(vars_array[i]);
    XLALFree(vars_array[i]);
  }
  XLALFree(vars_array);

  /* Close dataset. */
  XLALH5DatasetFree(dataset);

  /* Close file. */
  XLALH5FileClose(file);

  /* Check for memory leaks. */
  LALCheckMemoryLeaks();

  /* Done! */
  return gsl_test_summary();
}
