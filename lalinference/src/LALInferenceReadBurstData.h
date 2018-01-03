/*
 *
 *  LALInference:             LAL Inference library
 *  LALInferenceReadData.h    Utility functions for handling IFO data
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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

#ifndef LALInferenceReadNonCBCData_h
#define LALInferenceReadNonCBCData_h

#include <lal/LALInference.h>

/**
 * \defgroup LALInferenceReadData_h Header LALInferenceReadData.h
 * \ingroup pkg_LALInference
 * \brief Utility functions for handling IFO data
 */
/*@{*/

/** \brief Read IFO data according to command line arguments.
 * This function reads command line arguments and returns a \c LALInferenceIFOData linked
 * list.
 * \param IFOdata [out] A new \c LALInferenceIFOData linked list containing the data, or NULL upon error
 * \param commandLine [in] Pointer to a ProcessParamsTable containing command line arguments
 * \author John Veitch
 */
void LALInferenceInjectBurstSignal(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine);
void LALInferenceBurstInjectionToVariables(SimBurst *theEventTable, LALInferenceVariables *vars);
/*@}*/
#endif
