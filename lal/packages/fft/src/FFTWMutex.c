/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <lal/FFTWMutex.h>

#if defined(LAL_PTHREAD_LOCK) && defined(LAL_FFTW3_ENABLED)
static pthread_mutex_t lalFFTWMutex = PTHREAD_MUTEX_INITIALIZER;
#endif


/**
 * Aquire LAL's FFTW wisdom lock.  This lock must be held when creating or
 * destroying FFTW plans.  This function is a no-op if LAL has been
 * compiled without pthread support or with an FFT backend other than FFTW.
 *
 * See also:  XLALFFTWWisdomUnlock()
 */

void XLALFFTWWisdomLock(void)
{
#if defined(LAL_PTHREAD_LOCK) && defined(LAL_FFTW3_ENABLED)
    pthread_mutex_lock( &lalFFTWMutex );
#endif
}


/**
 * Release LAL's FFTW wisdom lock.  This function is a no-op if LAL has
 * been compiled without pthread support or with an FFT backend other than
 * FFTW.
 *
 * See also:  XLALFFTWWisdomLock()
 */

void XLALFFTWWisdomUnlock(void)
{
#if defined(LAL_PTHREAD_LOCK) && defined(LAL_FFTW3_ENABLED)
    pthread_mutex_unlock( &lalFFTWMutex );
#endif
}
