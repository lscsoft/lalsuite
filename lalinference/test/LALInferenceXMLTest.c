#include <lal/LALInference.h>
#include <lal/LALInferenceXML.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>

#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>


int testLALInferenceVariables(void);

int testLALInferenceVariables(void){
  LALInferenceVariables var;
  char *xmlString = NULL;
  xmlDocPtr xmlDocument = NULL;
  var.dimension=0;
  var.head=NULL;
  REAL8 r8test=42.0;
  xmlNodePtr xmlFragment;
  
  LALInferenceAddVariable(&var, "real8 test", (void *)&r8test, 
	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  
  printf("Initial LALInferenceVariables:\n");
  LALInferencePrintVariables(&var);
  printf( "--> Serializing into XML string ... ");
  xmlFragment=XLALInferenceVariables2VOTNode(&var, "test");
  /* convert VOTable tree into XML string */
   if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
      XLALPrintError ("%s: XLALCreateVOTStringFromTree() failed.\n", __func__);
      return 1;
    }
    printf ("ok.\n");
     /* display serialized structure */
    printf( "----------------------------------------------------------------------\n");
    printf( "%s", xmlString );
    printf( "----------------------------------------------------------------------\n");
   /* ---------- parse XML string back and validate ---------- */
    /* convert VOTable string back into VOTable document */
    printf ("--> Parsing XML string into xmlDoc ... ");
    if ( (xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
      XLALPrintError( "%s: XLALXMLString2Doc() failed.\n", __func__);
      return 1;
    }
    printf ("ok.\n");

    
    xmlFreeDoc(xmlDocument);
    xmlFreeNode ( xmlFragment );
    XLALFree ( xmlString );

    return 0;
}

int main(void)
{
  return(testLALInferenceVariables());
}