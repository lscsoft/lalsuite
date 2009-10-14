/* Scans OpenCL platforms & devices */
/* Bernd Machenschalk 2009 */
/* gcc -framework opencl QueryOpenCLDevices.c -o QueryOpenCLDevices */

#include <stdio.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif 

main() {
  cl_int err = 0;

  /* platforms */
  const cl_uint max_num_platforms = 8;
  cl_platform_id platform_ids[max_num_platforms];
  cl_uint num_platforms;
  unsigned int platform_id;

  /* info */
  cl_uint infoUInt;
  size_t infoSize;
  cl_bool infoBool;
  cl_device_type infoType;

  char infoStr[1024];
  size_t infoStrLen;

  /* devices */
  const cl_uint max_num_devices = 8;
  cl_device_id device_ids[max_num_devices];
  cl_uint num_devices;
  unsigned int device_id;

  /* get available platforms */
  err = clGetPlatformIDs(max_num_platforms, platform_ids, &num_platforms);
  if (err != CL_SUCCESS) {
    fprintf(stderr,"ERROR %d calling clGetPlatformIDs().\n", err);
    return -1;
  }
  printf("Found %d platforms:\n", num_platforms);
  
  /* loop over available platforms */
  for(platform_id = 0; platform_id < num_platforms; platform_id++) {
    printf("  Platform %d:\n", platform_id);

    err = clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_PROFILE, sizeof(infoStr), infoStr, &infoStrLen );
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetPlatformInfo().\n", err);
      return -1;
    }
    printf("    Profile: %s\n", infoStr);

    err = clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_VERSION, sizeof(infoStr), infoStr, &infoStrLen );
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetPlatformInfo().\n", err);
      return -1;
    }
    printf("    Version: %s\n", infoStr);

    err = clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_NAME, sizeof(infoStr), infoStr, &infoStrLen );
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetPlatformInfo().\n", err);
      return -1;
    }
    printf("    Name: %s\n", infoStr);

    err = clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_VENDOR, sizeof(infoStr), infoStr, &infoStrLen );
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetPlatformInfo().\n", err);
      return -1;
    }
    printf("    Vendor: %s\n", infoStr);

    err = clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_EXTENSIONS, sizeof(infoStr), infoStr, &infoStrLen );
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetPlatformInfo().\n", err);
      return -1;
    }
    printf("    Extensions: %s\n", infoStr);

    /* get the list of available devices per platform */
    err = clGetDeviceIDs( platform_ids[platform_id], CL_DEVICE_TYPE_ALL, max_num_devices, device_ids, &num_devices);
    if (err != CL_SUCCESS) {
      fprintf(stderr,"ERROR %d calling clGetDeviceIDs().\n", err);
      return -1;
    }
    printf("    Found %d devices:\n", num_devices);

    /* loop over available devices */
    for(device_id = 0; device_id < num_devices; device_id++) {
      printf("      Device %d:\n", device_id);

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_TYPE, sizeof(infoType), &infoType, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Type:");
      if (infoType & CL_DEVICE_TYPE_CPU)
	printf(" CPU");
      if (infoType & CL_DEVICE_TYPE_GPU)
	printf(" GPU");
      if (infoType & CL_DEVICE_TYPE_ACCELERATOR)
	printf(" ACCELERATOR");
      if (infoType & CL_DEVICE_TYPE_DEFAULT)
	printf(" DEFAULT");
      printf("\n");

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_PROFILE, sizeof(infoStr), infoStr, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Profile: %s\n", infoStr);

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_VERSION, sizeof(infoStr), infoStr, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Version: %s\n", infoStr);

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_NAME, sizeof(infoStr), infoStr, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Name: %s\n", infoStr);

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_VENDOR, sizeof(infoStr), infoStr, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Vendor: %s\n", infoStr);

      err = clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_EXTENSIONS, sizeof(infoStr), infoStr, &infoStrLen );
      if (err != CL_SUCCESS) {
	fprintf(stderr,"ERROR %d calling clGetDeviceInfo().\n", err);
	return -1;
      }
      printf("        Extensions: %s\n", infoStr);

    } /* loop over devices */
  } /* loop over platforms */
}
