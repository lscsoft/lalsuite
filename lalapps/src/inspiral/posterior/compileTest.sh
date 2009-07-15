gcc -Wall -g -o InferenceTest LALInference.c InferenceTest.c LALInferenceTemplate.c LALInferenceReadData.c `pkg-config --cflags --libs lal libframe lalsupport libmetaio lalmetaio lalframe fftw3`

