// include string for makeTclPointer kludge
#include <string>

#include "Translation.h"

#include "datacondAPI/TimeSeries.hh"
#include "datacondAPI/TimeBoundedFreqSequenceUDT.hh"

#include "general/types.hh"

enum {
  TRANS_OK=1,
  TRANS_FAIL=0
};


using General::GPSTime;
using namespace datacondAPI;

// begin makeTclPointer kludge
std::string
makeTclPointer(void const*, const std::string &)
{
  return "";
}
// end makeTclPointer kludge


int
TranslateREAL4TimeSeries( translation_direction Direction,
			  void** UserData,
			  void** DatacondData,
			  void* AuxData )
{
  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
  {
    // One of the pointers cannot be dereferenced
    return TRANS_FAIL;
  }
  switch ( Direction )
  {
  case DATACOND_SYMBOL_INPUT:
    {
      // Importing data into the algorithm section
      REAL4TimeSeries*	in =
	reinterpret_cast< REAL4TimeSeries* >(*UserData);
      //-----------------------------------------------------------------
      // Basic sanity checks
      //-----------------------------------------------------------------
      if ( ( in->data == (REAL4Sequence*)NULL ) ||
	   ( in->data->data == (REAL4*)NULL ) )
      {
	return TRANS_FAIL;
      }
      //-----------------------------------------------------------------
      // Create a Time Series UDT.
      //-----------------------------------------------------------------
      try
      {
	std::auto_ptr< TimeSeries< REAL4 > >
	  dc( new TimeSeries< REAL4 >( Sequence< REAL4 >( in->data->data,
							  in->data->length ),
				       in->f0,
				       GPSTime( in->epoch.gpsSeconds,
						in->epoch.gpsNanoSeconds ) ) );
	dc->SetName( in->name );
	*DatacondData = dc.release( );
      }
      catch( ... )
      {
	//---------------------------------------------------------------
	// Some error occurred. Do cleanup and leave
	//---------------------------------------------------------------
	return TRANS_FAIL;
      }
    }
    break;
  case DATACOND_SYMBOL_OUTPUT:
    {
      // Exporting data from the algorithm section
      udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
      if ( ! udt::IsA< TimeSeries< REAL4 > >(*udt_dc) )
      {
	// Wrong input type
	return TRANS_FAIL;
      }
      TimeSeries< REAL4 >&
	dc( *( dynamic_cast< TimeSeries< REAL4 >* >( udt_dc ) ) );
      //-----------------------------------------------------------------
      // Fill in the header section of the REAL4TimeSeries
      //-----------------------------------------------------------------
      REAL4TimeSeries*	user =
	(REAL4TimeSeries*)malloc( sizeof( REAL4TimeSeries ) );
      strcpy( user->name, dc.name( ).c_str( ) );
      user->epoch.gpsSeconds = dc.GetStartTime( ).GetSeconds( );
      user->epoch.gpsNanoSeconds = dc.GetStartTime( ).GetNanoseconds( );
      user->deltaT = 1.0 / dc.GetSampleRate( );
      user->f0 = dc.GetSampleRate( );
      // :TODO: user->sampleUnits = ;
      user->data = (REAL4Sequence*)malloc( sizeof( REAL4Sequence ) );
      //-----------------------------------------------------------------
      // Fill in the REAL4Sequence section
      //-----------------------------------------------------------------
      user->data->length = dc.size( );
      user->data->data = (REAL4*)calloc( sizeof(REAL4), dc.size( ) );
      std::copy( &( dc[ 0 ] ), &( dc[ dc.size( ) ] ), user->data->data );
      *UserData = user;
    }
    break;
  }
  return TRANS_OK;
 }



int TranslateREAL4Sequence( translation_direction Direction,
			    void** UserData,
			    void** DatacondData,
			    void* AuxData ) {
  // NOTE: AuxData contains name info.
  // Must be pre-allocated for OUTPUT.

  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
  {
    // One of the pointers cannot be dereferenced
    return TRANS_FAIL;
  }
  switch ( Direction )
  {
  case DATACOND_SYMBOL_INPUT:
    {
      // Importing data into the algorithm section
      REAL4Sequence*	in =
	reinterpret_cast< REAL4Sequence* >(*UserData);
      //-----------------------------------------------------------------
      // Basic sanity checks
      //-----------------------------------------------------------------
      if ( ( in->data == NULL ) )
      {
	return TRANS_FAIL;
      }
      //-----------------------------------------------------------------
      // Create a Time Series UDT.
      //-----------------------------------------------------------------
      try
      {
	std::auto_ptr< Sequence< REAL4 > >
	  dc( new Sequence< REAL4 >( in->data,
				     in->length ) );
	dc->SetName( (char *)AuxData );
	*DatacondData = dc.release( );
      }
      catch( ... )
      {
	//---------------------------------------------------------------
	// Some error occurred. Do cleanup and leave
	//---------------------------------------------------------------
	return TRANS_FAIL;
      }
    }
    break;
  case DATACOND_SYMBOL_OUTPUT:
    {
      // Exporting data from the algorithm section
      udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
      if ( ! udt::IsA< Sequence< REAL4 > >(*udt_dc) )
      {
	// Wrong input type
	return TRANS_FAIL;
      }
      Sequence< REAL4 >&
	dc( *( dynamic_cast< Sequence< REAL4 >* >( udt_dc ) ) );
      //-----------------------------------------------------------------
      // Fill in the header section of the REAL4Sequence
      //-----------------------------------------------------------------
      REAL4Sequence*	user =
	(REAL4Sequence*)malloc( sizeof( REAL4Sequence ) );

      //      strcpy( (char *)AuxData, dc.name( ).c_str( ) );
      user->length = dc.size( );
      user->data = (REAL4*)calloc( sizeof(REAL4), dc.size( ) );
      std::copy( &( dc[ 0 ] ), &( dc[ dc.size( ) ] ), user->data );
      *UserData = user;
    }
    break;
  }
  return TRANS_OK;
}

int TranslateCHARVectorSequence( translation_direction Direction,
				 void** UserData,
				 void** DatacondData,
				 void* AuxData ) {
  return TRANS_FAIL;
}

int TranslateREAL8FrequencySeries( translation_direction Direction,
				   void** UserData,
				   void** DatacondData,
				   void* AuxData ) {
  return TRANS_FAIL;

}

int TranslateREAL4FrequencySeries( translation_direction Direction,
				   void** UserData,
				   void** DatacondData,
				   void* AuxData ) {
  return TRANS_FAIL;

}

int TranslateCOMPLEX8FrequencySeries( translation_direction Direction,
				      void** UserData,
				      void** DatacondData,
				      void* AuxData ) {

  // NOTE: on input, AuxData is LIGOTimeGPS for end of time series

  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
    {
      // One of the pointers cannot be dereferenced
      return TRANS_FAIL;
    }
  switch ( Direction )
    {
    case DATACOND_SYMBOL_INPUT:
      {
	// Importing data into the algorithm section
	COMPLEX8FrequencySeries*	in =
	  reinterpret_cast< COMPLEX8FrequencySeries* >(*UserData);
	//-----------------------------------------------------------------
	// Basic sanity checks
	//-----------------------------------------------------------------
	if ( ( in->data == (COMPLEX8Sequence*)NULL ) ||
	     ( in->data->data == (COMPLEX8*)NULL ) )
	  {
	    return TRANS_FAIL;
	  }
	//-----------------------------------------------------------------
	// Create a Frequency Series UDT.
	//-----------------------------------------------------------------
	try
	  {
	    LIGOTimeGPS *endEpoch = (LIGOTimeGPS *)AuxData;
	    
	    std::auto_ptr< TimeBoundedFreqSequence< COMPLEX_8 > >
	      dc( new TimeBoundedFreqSequence< COMPLEX_8 >( 
			      Sequence< COMPLEX_8 >( in->data->length ),
			      in->f0,
			      in->deltaF,
			      GPSTime( in->epoch.gpsSeconds,
				       in->epoch.gpsNanoSeconds ),
			      GPSTime( endEpoch->gpsSeconds, 
				       endEpoch->gpsNanoSeconds) )
		  );
	    COMPLEX8* data = in->data->data;
	    for( unsigned int i=0, end = in->data->length;
		 i < end;
		 ++i, ++data) {
	      (*dc)[i] = std::complex<REAL4>(data->re, data->im);
	    }

	    dc->SetName( in->name );
	    *DatacondData = dc.release( );
	  }
	catch( ... )
	  {
	    //---------------------------------------------------------------
	    // Some error occurred. Do cleanup and leave
	    //---------------------------------------------------------------
	    return TRANS_FAIL;
	  }
      }
      break;
    case DATACOND_SYMBOL_OUTPUT:
      {
	// Exporting data from the algorithm section
	udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
	if ( ! udt::IsA< TimeBoundedFreqSequence< COMPLEX_8 > >(*udt_dc) )
	  {
	    // Wrong input type
	    return TRANS_FAIL;
	  }
	TimeBoundedFreqSequence< COMPLEX_8 >&
	  dc( *( dynamic_cast< TimeBoundedFreqSequence< COMPLEX_8 >* >( udt_dc ) ) );
	//-----------------------------------------------------------------
	// Fill in the header section of the COMPLEX8FrequencySeries
	//-----------------------------------------------------------------
	COMPLEX8FrequencySeries*	user =
	  (COMPLEX8FrequencySeries*)malloc( sizeof( COMPLEX8FrequencySeries ) );
	strcpy( user->name, dc.name( ).c_str( ) );
	user->epoch.gpsSeconds = dc.GetStartTime( ).GetSeconds( );
	user->epoch.gpsNanoSeconds = dc.GetStartTime( ).GetNanoseconds( );
	user->deltaF = dc.GetFrequencyDelta( );
	user->f0 = 0.0; // CAREFUL!!
	// :TODO: user->sampleUnits = ;
	user->data = (COMPLEX8Sequence*)malloc( sizeof( COMPLEX8Sequence ) );
	//-----------------------------------------------------------------
	// Fill in the COMPLEX8Sequence section
	//-----------------------------------------------------------------
	user->data->length = dc.size( );
	user->data->data = (COMPLEX8*)calloc( sizeof(COMPLEX8), dc.size( ) );

	
	COMPLEX8* data = user->data->data;
	COMPLEX_8* dc_data = &(dc[0]);
	for( unsigned int i=0, end = user->data->length;
	     i < end;
	     ++i, ++data, ++dc_data) {
	  data->re = dc_data->real();
	  data->im = dc_data->imag();
	}

	*UserData = user;
    }
    break;
  }
  return TRANS_OK;

}

int TranslateCOMPLEX16FrequencySeries( translation_direction Direction,
				   void** UserData,
				   void** DatacondData,
				   void* AuxData ) {
  return TRANS_FAIL;

}


int
TranslateREAL8TimeSeries( translation_direction Direction,
			  void** UserData,
			  void** DatacondData,
			  void* AuxData )
{
  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
  {
    // One of the pointers cannot be dereferenced
    return TRANS_FAIL;
  }
  switch ( Direction )
  {
  case DATACOND_SYMBOL_INPUT:
    {
      // Importing data into the algorithm section
      REAL8TimeSeries*	in =
	reinterpret_cast< REAL8TimeSeries* >(*UserData);
      //-----------------------------------------------------------------
      // Basic sanity checks
      //-----------------------------------------------------------------
      if ( ( in->data == (REAL8Sequence*)NULL ) ||
	   ( in->data->data == (REAL8*)NULL ) )
      {
	return TRANS_FAIL;
      }
      //-----------------------------------------------------------------
      // Create a Time Series UDT.
      //-----------------------------------------------------------------
      try
      {
	std::auto_ptr< TimeSeries< REAL8 > >
	  dc( new TimeSeries< REAL8 >( Sequence< REAL8 >( in->data->data,
							  in->data->length ),
				       in->f0,
				       GPSTime( in->epoch.gpsSeconds,
						in->epoch.gpsNanoSeconds ) ) );
	dc->SetName( in->name );
	*DatacondData = dc.release( );
      }
      catch( ... )
      {
	//---------------------------------------------------------------
	// Some error occurred. Do cleanup and leave
	//---------------------------------------------------------------
	return TRANS_FAIL;
      }
    }
    break;
  case DATACOND_SYMBOL_OUTPUT:
    {
      // Exporting data from the algorithm section
      udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
      if ( ! udt::IsA< TimeSeries< REAL8 > >(*udt_dc) )
      {
	// Wrong input type
	return TRANS_FAIL;
      }
      TimeSeries< REAL8 >&
	dc( *( dynamic_cast< TimeSeries< REAL8 >* >( udt_dc ) ) );
      //-----------------------------------------------------------------
      // Fill in the header section of the REAL8TimeSeries
      //-----------------------------------------------------------------
      REAL8TimeSeries*	user =
	(REAL8TimeSeries*)malloc( sizeof( REAL8TimeSeries ) );
      strcpy( user->name, dc.name( ).c_str( ) );
      user->epoch.gpsSeconds = dc.GetStartTime( ).GetSeconds( );
      user->epoch.gpsNanoSeconds = dc.GetStartTime( ).GetNanoseconds( );
      user->deltaT = 1.0 / dc.GetSampleRate( );
      user->f0 = dc.GetSampleRate( );
      // :TODO: user->sampleUnits = ;
      user->data = (REAL8Sequence*)malloc( sizeof( REAL8Sequence ) );
      //-----------------------------------------------------------------
      // Fill in the REAL8Sequence section
      //-----------------------------------------------------------------
      user->data->length = dc.size( );
      user->data->data = (REAL8*)calloc( sizeof(REAL8), dc.size( ) );
      std::copy( &( dc[ 0 ] ), &( dc[ dc.size( ) ] ), user->data->data );
      *UserData = user;
    }
    break;
  }
  return TRANS_OK;
 }




int TranslateCOMPLEX8TimeSeries( translation_direction Direction,
				 void** UserData,
				 void** DatacondData,
				 void* AuxData ) {
  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
    {
      // One of the pointers cannot be dereferenced
      return TRANS_FAIL;
    }
  switch ( Direction )
    {
    case DATACOND_SYMBOL_INPUT:
      {
	// Importing data into the algorithm section
	COMPLEX8TimeSeries*	in =
	  reinterpret_cast< COMPLEX8TimeSeries* >(*UserData);
	//-----------------------------------------------------------------
	// Basic sanity checks
	//-----------------------------------------------------------------
	if ( ( in->data == (COMPLEX8Sequence*)NULL ) ||
	     ( in->data->data == (COMPLEX8*)NULL ) )
	  {
	    return TRANS_FAIL;
	  }
	//-----------------------------------------------------------------
	// Create a Time Series UDT.
	//-----------------------------------------------------------------
	try
	  {
	    std::auto_ptr< TimeSeries< COMPLEX_8 > >
	      dc( new TimeSeries< COMPLEX_8 >( Sequence< COMPLEX_8 >( in->data->length ),

					       in->f0,
					       GPSTime( in->epoch.gpsSeconds,
							in->epoch.gpsNanoSeconds ) ) );
	    COMPLEX8* data = in->data->data;
	    for( unsigned int i=0, end = in->data->length;
		 i < end;
		 ++i, ++data) {
	      (*dc)[i] = std::complex<REAL4>(data->re, data->im);
	    }

	    dc->SetName( in->name );
	    *DatacondData = dc.release( );
	  }
	catch( ... )
	  {
	    //---------------------------------------------------------------
	    // Some error occurred. Do cleanup and leave
	    //---------------------------------------------------------------
	    return TRANS_FAIL;
	  }
      }
      break;
    case DATACOND_SYMBOL_OUTPUT:
      {
	// Exporting data from the algorithm section
	udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
	if ( ! udt::IsA< TimeSeries< COMPLEX_8 > >(*udt_dc) )
	  {
	    // Wrong input type
	    return TRANS_FAIL;
	  }
	TimeSeries< COMPLEX_8 >&
	  dc( *( dynamic_cast< TimeSeries< COMPLEX_8 >* >( udt_dc ) ) );
	//-----------------------------------------------------------------
	// Fill in the header section of the COMPLEX8TimeSeries
	//-----------------------------------------------------------------
	COMPLEX8TimeSeries*	user =
	  (COMPLEX8TimeSeries*)malloc( sizeof( COMPLEX8TimeSeries ) );
	strcpy( user->name, dc.name( ).c_str( ) );
	user->epoch.gpsSeconds = dc.GetStartTime( ).GetSeconds( );
	user->epoch.gpsNanoSeconds = dc.GetStartTime( ).GetNanoseconds( );
	user->deltaT = 1.0 / dc.GetSampleRate( );
	user->f0 = dc.GetSampleRate( );
	// :TODO: user->sampleUnits = ;
	user->data = (COMPLEX8Sequence*)malloc( sizeof( COMPLEX8Sequence ) );
	//-----------------------------------------------------------------
	// Fill in the COMPLEX8Sequence section
	//-----------------------------------------------------------------
	user->data->length = dc.size( );
	user->data->data = (COMPLEX8*)calloc( sizeof(COMPLEX8), dc.size( ) );

	
	COMPLEX8* data = user->data->data;
	COMPLEX_8* dc_data = &(dc[0]);
	for( unsigned int i=0, end = user->data->length;
	     i < end;
	     ++i, ++data, ++dc_data) {
	  data->re = dc_data->real();
	  data->im = dc_data->imag();
	}

	*UserData = user;
    }
    break;
  }
  return TRANS_OK;
}






int TranslateCOMPLEX16TimeSeries( translation_direction Direction,
				 void** UserData,
				 void** DatacondData,
				 void* AuxData ) {
  //---------------------------------------------------------------------
  // Get the basic information
  //---------------------------------------------------------------------
  if ( ( UserData == (void**)NULL ) || ( DatacondData == (void**)NULL ) )
    {
      // One of the pointers cannot be dereferenced
      return TRANS_FAIL;
    }
  switch ( Direction )
    {
    case DATACOND_SYMBOL_INPUT:
      {
	// Importing data into the algorithm section
	COMPLEX16TimeSeries*	in =
	  reinterpret_cast< COMPLEX16TimeSeries* >(*UserData);
	//-----------------------------------------------------------------
	// Basic sanity checks
	//-----------------------------------------------------------------
	if ( ( in->data == (COMPLEX16Sequence*)NULL ) ||
	     ( in->data->data == (COMPLEX16*)NULL ) )
	  {
	    return TRANS_FAIL;
	  }
	//-----------------------------------------------------------------
	// Create a Time Series UDT.
	//-----------------------------------------------------------------
	try
	  {
	    std::auto_ptr< TimeSeries< COMPLEX_16 > >
	      dc( new TimeSeries< COMPLEX_16 >( Sequence< COMPLEX_16 >( in->data->length ),

					       in->f0,
					       GPSTime( in->epoch.gpsSeconds,
							in->epoch.gpsNanoSeconds ) ) );
	    COMPLEX16* data = in->data->data;
	    for( unsigned int i=0, end = in->data->length;
		 i < end;
		 ++i, ++data) {
	      (*dc)[i] = std::complex<REAL8>(data->re, data->im);
	    }

	    dc->SetName( in->name );
	    *DatacondData = dc.release( );
	  }
	catch( ... )
	  {
	    //---------------------------------------------------------------
	    // Some error occurred. Do cleanup and leave
	    //---------------------------------------------------------------
	    return TRANS_FAIL;
	  }
      }
      break;
    case DATACOND_SYMBOL_OUTPUT:
      {
	// Exporting data from the algorithm section
	udt* udt_dc = reinterpret_cast< udt* >( *DatacondData );
	if ( ! udt::IsA< TimeSeries< COMPLEX_16 > >(*udt_dc) )
	  {
	    // Wrong input type
	    return TRANS_FAIL;
	  }
	TimeSeries< COMPLEX_16 >&
	  dc( *( dynamic_cast< TimeSeries< COMPLEX_16 >* >( udt_dc ) ) );
	//-----------------------------------------------------------------
	// Fill in the header section of the COMPLEX16TimeSeries
	//-----------------------------------------------------------------
	COMPLEX16TimeSeries*	user =
	  (COMPLEX16TimeSeries*)malloc( sizeof( COMPLEX16TimeSeries ) );
	strcpy( user->name, dc.name( ).c_str( ) );
	user->epoch.gpsSeconds = dc.GetStartTime( ).GetSeconds( );
	user->epoch.gpsNanoSeconds = dc.GetStartTime( ).GetNanoseconds( );
	user->deltaT = 1.0 / dc.GetSampleRate( );
	user->f0 = dc.GetSampleRate( );
	// :TODO: user->sampleUnits = ;
	user->data = (COMPLEX16Sequence*)malloc( sizeof( COMPLEX16Sequence ) );
	//-----------------------------------------------------------------
	// Fill in the COMPLEX16Sequence section
	//-----------------------------------------------------------------
	user->data->length = dc.size( );
	user->data->data = (COMPLEX16*)calloc( sizeof(COMPLEX16), dc.size( ) );

	
	COMPLEX16* data = user->data->data;
	COMPLEX_16* dc_data = &(dc[0]);
	for( unsigned int i=0, end = user->data->length;
	     i < end;
	     ++i, ++data, ++dc_data) {
	  data->re = dc_data->real();
	  data->im = dc_data->imag();
	}

	*UserData = user;
    }
    break;
  }
  return TRANS_OK;
}

