#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "datacondAPI/DatacondCaller.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "lal/LALDatatypes.h"

  int TranslateREAL4TimeSeries( translation_direction Direction,
				void** UserData,
				void** DatacondData,
				void* AuxData );

  int TranslateREAL8TimeSeries( translation_direction Direction,
				void** UserData,
				void** DatacondData,
				void* AuxData );

  int TranslateCOMPLEX8TimeSeries( translation_direction Direction,
				   void** UserData,
				   void** DatacondData,
				   void* AuxData );

  int TranslateCOMPLEX16TimeSeries( translation_direction Direction,
				    void** UserData,
				    void** DatacondData,
				    void* AuxData );

  int TranslateREAL4Sequence( translation_direction Direction,
			      void** UserData,
			      void** DatacondData,
			      void* AuxData );

  int TranslateCHARVectorSequence( translation_direction Direction,
				   void** UserData,
				   void** DatacondData,
				   void* AuxData );

  int TranslateCOMPLEX8FrequencySeries( translation_direction Direction,
					void** UserData,
					void** DatacondData,
					void* AuxData );

  int TranslateCOMPLEX16FrequencySeries( translation_direction Direction,
					 void** UserData,
					 void** DatacondData,
					 void* AuxData );

  int TranslateREAL4FrequencySeries( translation_direction Direction,
				     void** UserData,
				     void** DatacondData,
				     void* AuxData );

  int TranslateREAL8FrequencySeries( translation_direction Direction,
				     void** UserData,
				     void** DatacondData,
				     void* AuxData );

#ifdef __cplusplus
}
#endif

#endif /* REAL4TimeSeries_H */
