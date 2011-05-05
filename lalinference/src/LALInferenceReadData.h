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

/**
 * \file LALInferenceReadData.h
 * \brief Utility functions for handling IFO data
 */

#ifndef LALInferenceReadData_h
#define LALInferenceReadData_h


#include <lal/LALInference.h>

/** Read IFO data according to command line arguments */
struct tagLALInferenceIFOData * LALInferenceReadData (ProcessParamsTable * commandLine);

/** Convenience function to inject a signal into the data, using a SimInspiralTable */
void LALInferenceInjectInspiralSignal(struct tagLALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine);

#endif
