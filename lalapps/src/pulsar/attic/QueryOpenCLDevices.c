/* Scans OpenCL platforms & devices */
/* Bernd Machenschalk 2009 */
/* gcc -framework opencl QueryOpenCLDevices.c -o QueryOpenCLDevices */

#include <stdio.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif 

#define OPENCL_TRY(name,fun)					\
  {								\
    cl_int err = fun;						\
    if (err != CL_SUCCESS) {					\
      fprintf(stderr,"ERROR %d calling %s().\n", err,name);	\
      return -1;						\
    }								\
  }

main() {
  /* platforms */
  const cl_uint max_num_platforms = 8;
  cl_platform_id platform_ids[max_num_platforms];
  cl_uint num_platforms;
  unsigned int platform_id;

  /* info */
  cl_uint infoUInt;
  cl_ulong infoULong;
  size_t infoSize;
  cl_bool infoBool;
  cl_device_type infoType;
  char infoStr[1024];
  size_t infoLen;

  /* devices */
  const cl_uint max_num_devices = 8;
  cl_device_id device_ids[max_num_devices];
  cl_uint num_devices;
  unsigned int device_id;

  /* get available platforms */
  OPENCL_TRY("clGetPlatformIDs", clGetPlatformIDs(max_num_platforms, platform_ids, &num_platforms));
  printf("Found %d platforms:\n", num_platforms);
  
  /* loop over available platforms */
  for(platform_id = 0; platform_id < num_platforms; platform_id++) {
    printf("  Platform %d:\n", platform_id);

    OPENCL_TRY("clGetPlatformInfo", clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_PROFILE, sizeof(infoStr), infoStr, &infoLen ));
    printf("    Profile: %s\n", infoStr);

    OPENCL_TRY("clGetPlatformInfo", clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_VERSION, sizeof(infoStr), infoStr, &infoLen ));
    printf("    Version: %s\n", infoStr);

    OPENCL_TRY("clGetPlatformInfo", clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_NAME, sizeof(infoStr), infoStr, &infoLen ));
    printf("    Name: %s\n", infoStr);

    OPENCL_TRY("clGetPlatformInfo", clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_VENDOR, sizeof(infoStr), infoStr, &infoLen ));
    printf("    Vendor: %s\n", infoStr);

    OPENCL_TRY("clGetPlatformInfo", clGetPlatformInfo ( platform_ids[platform_id], CL_PLATFORM_EXTENSIONS, sizeof(infoStr), infoStr, &infoLen ));
    printf("    Extensions: %s\n", infoStr);

    /* get the list of available devices per platform */
    OPENCL_TRY("clGetDeviceIDs", clGetDeviceIDs( platform_ids[platform_id], CL_DEVICE_TYPE_ALL, max_num_devices, device_ids, &num_devices));
    printf("    Found %d devices:\n", num_devices);

    /* loop over available devices */
    for(device_id = 0; device_id < num_devices; device_id++) {
      printf("      Device %d:\n", device_id);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_TYPE, sizeof(infoType), &infoType, &infoLen ));
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

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_PROFILE, sizeof(infoStr), infoStr, &infoLen ));
      printf("        PROFILE: %s\n", infoStr);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_VERSION, sizeof(infoStr), infoStr, &infoLen ));
      printf("        VERSION: %s\n", infoStr);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_NAME, sizeof(infoStr), infoStr, &infoLen ));
      printf("        NAME: %s\n", infoStr);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_VENDOR, sizeof(infoStr), infoStr, &infoLen ));
      printf("        VENDOR: %s\n", infoStr);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_EXTENSIONS, sizeof(infoStr), infoStr, &infoLen ));
      printf("        EXTENSIONS: %s\n", infoStr);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_VENDOR_ID, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        VENDOR_ID: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_COMPUTE_UNITS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_WORK_ITEM_DIMENSIONS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_CLOCK_FREQUENCY: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_ADDRESS_BITS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        ADDRESS_BITS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_MEM_ALLOC_SIZE , sizeof(infoULong), &infoULong, &infoLen ));
      printf("        MAX_MEM_ALLOC_SIZE: %lu\n", infoULong);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_IMAGE_SUPPORT , sizeof(infoBool), &infoBool, &infoLen ));
      printf("        IMAGE_SUPPORT: %s\n", infoBool?"TRUE":"FALSE");

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_READ_IMAGE_ARGS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_READ_IMAGE_ARGS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_WRITE_IMAGE_ARGS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_WRITE_IMAGE_ARGS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_SAMPLERS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_SAMPLERS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MEM_BASE_ADDR_ALIGN: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MIN_DATA_TYPE_ALIGN_SIZE: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        GLOBAL_MEM_CACHELINE_SIZE: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_GLOBAL_MEM_CACHE_SIZE , sizeof(infoULong), &infoULong, &infoLen ));
      printf("        GLOBAL_MEM_CACHE_SIZE: %lu\n", infoULong);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_GLOBAL_MEM_SIZE , sizeof(infoULong), &infoULong, &infoLen ));
      printf("        GLOBAL_MEM_SIZE: %lu\n", infoULong);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE , sizeof(infoULong), &infoULong, &infoLen ));
      printf("        MAX_CONSTANT_BUFFER_SIZE: %lu\n", infoULong);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_MAX_CONSTANT_ARGS, sizeof(infoUInt), &infoUInt, &infoLen ));
      printf("        MAX_CONSTANT_ARGS: %u\n", infoUInt);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_LOCAL_MEM_SIZE , sizeof(infoULong), &infoULong, &infoLen ));
      printf("        LOCAL_MEM_SIZE: %lu\n", infoULong);

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_ERROR_CORRECTION_SUPPORT , sizeof(infoBool), &infoBool, &infoLen ));
      printf("        ERROR_CORRECTION_SUPPORT: %s\n", infoBool?"TRUE":"FALSE");

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_ENDIAN_LITTLE , sizeof(infoBool), &infoBool, &infoLen ));
      printf("        ENDIAN_LITTLE: %s\n", infoBool?"TRUE":"FALSE");

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_AVAILABLE , sizeof(infoBool), &infoBool, &infoLen ));
      printf("        AVAILABLE: %s\n", infoBool?"TRUE":"FALSE");

      OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_COMPILER_AVAILABLE , sizeof(infoBool), &infoBool, &infoLen ));
      printf("        COMPILER_AVAILABLE: %s\n", infoBool?"TRUE":"FALSE");

      {
	cl_device_fp_config infoBitField;
	OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_SINGLE_FP_CONFIG, sizeof(infoBitField), &infoBitField, &infoLen ));
	printf("        CL_DEVICE_SINGLE_FP_CONFIG:");
	if(infoBitField & CL_FP_DENORM)
	  printf(" CL_FP_DENORM");
	if(infoBitField & CL_FP_INF_NAN)
	  printf(" CL_FP_INF_NAN");
	if(infoBitField & CL_FP_ROUND_TO_NEAREST)
	  printf(" CL_FP_ROUND_TO_NEAREST");
	if(infoBitField & CL_FP_ROUND_TO_ZERO)
	  printf(" CL_FP_ROUND_TO_ZERO");
	if(infoBitField & CL_FP_ROUND_TO_INF)
	  printf(" CL_FP_ROUND_TO_INF");
	if(infoBitField & CL_FP_FMA)
	  printf(" CL_FP_FMA");
	printf("\n");
      }

      {
	cl_device_mem_cache_type infoBitField;
	OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_GLOBAL_MEM_CACHE_TYPE, sizeof(infoBitField), &infoBitField, &infoLen ));
	printf("        CL_DEVICE_GLOBAL_MEM_CACHE_TYPE:");
	if(infoBitField & CL_NONE)
	  printf(" CL_NONE");
	if(infoBitField & CL_READ_ONLY_CACHE)
	  printf(" CL_READ_ONLY_CACHE");
	if(infoBitField & CL_READ_WRITE_CACHE)
	  printf(" CL_READ_WRITE_CACHE");
	printf("\n");
      }

      {
	cl_device_exec_capabilities infoBitField;
	OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_EXECUTION_CAPABILITIES, sizeof(infoBitField), &infoBitField, &infoLen ));
	printf("        CL_DEVICE_EXECUTION_CAPABILITIES:");
	if(infoBitField & CL_EXEC_KERNEL )
	  printf(" CL_EXEC_KERNEL");
	if(infoBitField & CL_EXEC_NATIVE_KERNEL)
	  printf(" CL_EXEC_NATIVE_KERNEL");
	printf("\n");
      }

      {
	cl_command_queue_properties infoBitField;
	OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_QUEUE_PROPERTIES, sizeof(infoBitField), &infoBitField, &infoLen ));
	printf("        CL_DEVICE_QUEUE_PROPERTIES:");
	if(infoBitField & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)
	  printf(" CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE");
	if(infoBitField & CL_QUEUE_PROFILING_ENABLE)
	  printf(" CL_QUEUE_PROFILING_ENABLE");
	printf("\n");
      }

      {
	cl_device_local_mem_type infoBitField;
	OPENCL_TRY("clGetDeviceInfo", clGetDeviceInfo ( device_ids[device_id], CL_DEVICE_LOCAL_MEM_TYPE, sizeof(infoBitField), &infoBitField, &infoLen ));
	printf("        CL_DEVICE_LOCAL_MEM_TYPE:");
	if(infoBitField == CL_LOCAL)
	  printf(" CL_LOCAL");
	if(infoBitField == CL_GLOBAL)
	  printf(" CL_GLOBAL");
	printf("\n");
      }

    } /* loop over devices */
  } /* loop over platforms */
}
