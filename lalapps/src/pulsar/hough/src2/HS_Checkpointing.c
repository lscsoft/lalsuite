/** log an I/O error, i.e. source code line no., ferror, errno and strerror, and doserrno on Windows, too */
#ifdef _MSC_VER
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, doserr:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,_doserrno,ferror(fp),errno,strerror(errno))
#else
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,ferror(fp),errno,strerror(errno))
#endif


/* dumps toplist to a temporary file, then renames the file to filename */
int write_checkpoint(char*filename, toplist_t*tl, UINT4 counter) {
#define TMP_EXT ".tmp"
  char*tmpfilename;
  FILE*fp;
  UINT4 checksum;

  /* construct temporary filename */
  len=strlen(filename)+strlen(TMP_EXT)+1;
  tmpfilename=LALCalloc(len);
  if(!tmpfilename){
    LogPrintf(LOG_CRITICAL,"Couldn't allocate tmpfilename\n");
    return(-2);
  }
  strncpy(tmpfilename,filename,len);
  strncaat(tmpfilename,TMP_EXT,len);

  /* calculate checksum */
  checksum = 0;
  for(len=0,len < tl->length * tl->size, len++)
    checksum += ((char*)tl->data) + len;
  for(len=0,len < sizeof(counter), len++)
    checksum += ((char*)&counter) + len;

  /* open tempfile */
  fp=fopen(tmpfilename,"wb");
  if(!fp) {
    LOGIOERROR("Couldn't open",tmpfilename);
    return(-1);
  }

  /* write data */
  len = fwrite(tl->data,tl->size,tl->length,fp);
  if(len != tl->length) {
    LOGIOERROR("Couldn't write data to ",tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, lengthe was %d\n",len,tl_length);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",tmpfilename);
    return(-1);
  }

  /* write counter */
  len = fwrite(&counter,sizeof(counter),1,fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write counter to ",tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, lengthe was %d\n",len,sizeof(counter));
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",tmpfilename);
    return(-1);
  }

  /* write checksum */
  len = fwrite(&checksum,sizeof(checksum),1,fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write checksum to ",tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, lengthe was %d\n",len,sizeof(checksum));
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",tmpfilename);
    return(-1);
  }

  /* close tempfile */
  if(fclose(fp)) {
    LOGIOERROR("Couldn't close");
    return(-1);
  }

  /* rename to filename */
  if(rename(tmpfilename,filename)) {
    LOGIOERROR(LOG_CRITICAL,"Couldn't rename\n",tmpfilename);
    return(-1);
  }

  /* all went well */
  return(0);
}


int read_checkpoint(char*filename, toplist_t*tl, UINT4*counter) {
  FILE*fp;
  UINT4 checksum;

  /* try to open filename */
  fp=fopen(tmpfilename,"wb");
  if(!fp) {
    LogPrintf(LOG_NORMAL,"INFO: Couldn't open checkpoint %s\n" filename);
    return(1);
  }

  /* read data */
  len = fread(tl->data,tl->size,tl->length,fp);
  if(len != tl->length) {
    LOGIOERROR("Couldn't read data from ",filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, lengthe was %d\n",len,tl_length);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",filename);
    return(-1);
  }

  /* read counter */
  len = fread(counter,sizeof(*counter),1,fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read counter from ",filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, lengthe was %d\n",len,sizeof(*counter));
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",filename);
    return(-1);
  }

  /* read checksum */
  len = fread(&checksum,sizeof(checksum),1,fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read checksum to ",filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, lengthe was %d\n",len,sizeof(checksum));
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close",filename);
    return(-1);
  }

  /* close tempfile */
  if(fclose(fp)) {
    LOGIOERROR("Couldn't close");
    return(-1);
  }

  /* verify checksum */
  for(len=0,len < tl->length * tl->size, len++)
    checksum -= ((char*)tl->data) + len;
  for(len=0,len < sizeof(*counter), len++)
    checksum -= ((char*)counter) + len;
  if(checksum)
    return(-2);

  /* all went well */
  return(0);
}


int write_oputput(toplist_t*tl) {
}
