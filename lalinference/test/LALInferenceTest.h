// Some utilities for quick and dirty debugging.
#define ECHO_F(x) printf("%f\n", x)
#define ECHO_S(x) printf("%s\n", x)
#define ECHO_I(x) printf("%i\n", x)

// Define test utility macros.
#define TEST_RUN(name, count) \
{ \
	count += name(); \
	printf("\n"); \
}
#define TEST_HEADER() \
	printf("Testing: %s\n", __func__); \
	int test_failure_count = 0
#define TEST_FOOTER() \
{ \
	if (test_failure_count > 0) \
		printf("%i test(s) failed.\n", test_failure_count); \
	else \
		printf("All tests passed.\n"); \
	return test_failure_count; \
}
#define TEST_FAIL(...) \
{ \
	fprintf(stderr, "FAIL - %s (%s:%i): ", __func__, __FILE__, __LINE__); \
	fprintf(stderr, __VA_ARGS__); \
	fprintf(stderr, "\n"); \
	test_failure_count += 1; \
}

int compareFloats(REAL8 x, REAL8 y, REAL8 epsilon);

int compareFloats(REAL8 x, REAL8 y, REAL8 epsilon)
{
	return fabs(x - y) < epsilon;
}
