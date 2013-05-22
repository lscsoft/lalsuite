#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceXML.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>

int testLALInferenceVariables(void);

int testLALInferenceVariables(void){
  LALInferenceVariables var;
  LALInferenceVariables *vars;
  char *xmlString = NULL;
  xmlDocPtr xmlDocument = NULL;
  var.dimension=0;
  var.head=NULL;
  REAL8 r8test=42.8,r8test2=101.0;
  INT4 i4test=12;
	xmlNodePtr xmlFragment;
  xmlNodePtr xmlTable;
  
  LALInferenceAddVariable(&var, "real8 test", (void *)&r8test, 
	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(&var, "param test", (void *)&r8test,
	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(&var, "field test 2",(void *)&r8test2,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(&var, "int test",(void *)&i4test,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(&var, "int field test",(void *)&i4test,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);

	
  printf("Initial LALInferenceVariables:\n");
  LALInferencePrintVariables(&var);
  printf( "--> Serializing into XML string ... ");
  xmlFragment=XLALInferenceVariables2VOTParamNode(&var);
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
    
    /* Convert array of variables into table */
    printf( "--> Serializing array of variables into XML Table ... ");

    vars=XLALCalloc(3,sizeof(LALInferenceVariables));
    int i;
    for(i=0;i<3;i++)
    {
      printf("Copying %i\n",i);
      LALInferenceCopyVariables(&var,&(vars[i]));
    }
    printf("Creating XML Table...\n");
      xmlTable=XLALInferenceVariablesArray2VOTTable(vars, 3, "Test table");
      printf("Created XML Table, tree = %lx ...\n",(long unsigned int)xmlTable);
      xmlString = XLALCreateVOTStringFromTree ( xmlTable );
      printf("Created XML String %s\n",xmlString);
      if(xmlString == NULL ) {
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
    
		FILE *outDoc=fopen("test_vot.xml","w");
	fprintf(outDoc,"%s",xmlString);
	fclose(outDoc);
    
    xmlFreeDoc(xmlDocument);
    xmlFreeNode ( xmlTable );
    XLALFree ( xmlString );

    return 0;
}

int main(void)
{
  return(testLALInferenceVariables());
}
