#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Computes the PN Coefficients for using in the TaylorT2 phasing equation.
 *
 * Terms given in equation 3.8a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_0PNCoeff(
	REAL8 eta)
{
	return -1./(32.*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_2PNCoeff(
	REAL8 eta)
{
	return 3.715/1.008 + 5.5/1.2 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -10. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_4PNCoeff(
	REAL8 eta)
{
	return 15.293365/1.016064 + 27.145/1.008 * eta + 30.85/1.44 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_5PNCoeff(
	REAL8 eta)
{
	return (386.45/6.72 - 65./8. * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_6PNCoeff(
	REAL8 eta)
{
	return 1234.8611926451/1.8776862720 - 160./3. * LAL_PI*LAL_PI - 171.2/2.1 * LAL_GAMMA
		+ (225.5/4.8 * LAL_PI*LAL_PI - 1573.7765635/1.2192768) * eta
		+ 76.055/6.912 * eta*eta - 127.825/5.184 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -85.6/2.1;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_7PNCoeff(
	REAL8 eta)
{
	return (77.096675/2.032128 + 37.8515/1.2096 * eta - 74.045/6.048 * eta*eta) * LAL_PI;
}


/**
 * Computes the PN Coefficients for using in the TaylorT2 timing equation.
 *
 * Terms given in equation 3.8b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_0PNCoeff(
	REAL8 totalmass,
	REAL8 eta)
{
	totalmass *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */
	return -5.*totalmass/(256.*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_2PNCoeff(
	REAL8 eta)
{
	return 7.43/2.52 + 11./3. * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -32./5. * LAL_PI;;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_4PNCoeff(
	REAL8 eta)
{
	return 30.58673/5.08032 + 54.29/5.04*eta + 61.7/7.2*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_5PNCoeff(
	REAL8 eta)
{
	return -(77.29/2.52 -13./3.*eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_6PNCoeff(
	REAL8 eta)
{
	return -1005.2469856691/2.3471078400 + 128./3. * LAL_PI*LAL_PI + 68.48/1.05 * LAL_GAMMA
		+ (3147.553127/3.048192 - 45.1/1.2 * LAL_PI*LAL_PI) * eta
		- 15.211/1.728 * eta*eta + 25.565/1.296 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return 34.24/1.05;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_7PNCoeff(
	REAL8 eta)
{
	return (-154.19335/1.27008 - 757.03/7.56 * eta + 148.09/3.78 * eta*eta) * LAL_PI;
}


/**
 * Computes the PN Coefficients for using in the TaylorT3 phasing equation.
 *
 * Terms given in equation 3.10a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */


static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_0PNCoeff(
	REAL8 eta)
{
	return -1./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_2PNCoeff(
	REAL8 eta)
{
	return 3.715/8.064 + 5.5/9.6 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -3./4. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_4PNCoeff(
	REAL8 eta)
{
	return 9.275495/14.450688 + 2.84875/2.58048 * eta + 1.855/2.048 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_5PNCoeff(
	REAL8 eta)
{
	return (3.8645/2.1504 - 6.5/25.6 * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_6PNCoeff(
	REAL8 eta)
{
	return 83.1032450749357/5.7682522275840 - 5.3/4.0 * LAL_PI*LAL_PI - 10.7/5.6 * LAL_GAMMA
		+ (-126.510089885/4.161798144 + 2.255/2.048 * LAL_PI*LAL_PI) * eta
		+ 1.54565/18.35008 * eta*eta - 1.179625/1.769472 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -10.7/5.6;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_7PNCoeff(
	REAL8 eta)
{
	return (1.88516689/1.73408256 + 4.88825/5.16096 * eta - 1.41769/5.16096 * eta*eta) * LAL_PI;
}


/**
 * Computes the PN Coefficients for using in the TaylorT3 frequency equation.
 *
 * Terms given in equation 3.10b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_0PNCoeff(
	REAL8 totalmass)
{
	totalmass *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */
	return 1. / (8. * LAL_PI * totalmass);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_2PNCoeff(
	REAL8 eta)
{
	return 7.43/26.88 + 1.1/3.2 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -3./10. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_4PNCoeff(
	REAL8 eta)
{
	return 1.855099/14.450688 + 5.6975/25.8048 * eta + 3.71/20.48 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_5PNCoeff(
	REAL8 eta)
{
	return (-7.729/21.504 + 1.3/25.6 * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_6PNCoeff(
	REAL8 eta)
{
	return -7.20817631400877/2.88412611379200 + 5.3/20.0 * LAL_PI*LAL_PI + 1.07/2.80 * LAL_GAMMA
		+ (25.302017977/4.161798144 - 4.51/20.48 * LAL_PI*LAL_PI) * eta
		- 3.0913/183.5008 * eta*eta + 2.35925/17.69472 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return 1.07/2.80;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_7PNCoeff(
	REAL8 eta)
{
	return (-1.88516689/4.33520640 - 9.7765/25.8048 * eta + 1.41769/12.90240 * eta*eta) * LAL_PI;
}

