gcc -Wall -ggdb -std=gnu99 -o InferenceTest LALInference.c InferenceTest.c LALInferenceTemplate.c LALInferenceReadData.c LALInferenceMCMC.c `pkg-config --cflags --libs lal libframe lalsupport libmetaio lalmetaio lalframe fftw3 lalinspiral lalpulsar`

