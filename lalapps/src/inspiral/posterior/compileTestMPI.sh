openmpicc -Wall -ggdb -std=gnu99 -o InferenceTest LALInference.c InferenceTest.c LALInferenceTemplate.c LALInferenceReadData.c LALInferenceMCMCMPI.c `pkg-config --cflags --libs lal libframe lalsupport libmetaio lalmetaio lalframe fftw3 lalinspiral lalpulsar`

