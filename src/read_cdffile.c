/**
 *  These C routines read data held CDF files, one of the data formats
 *  used by a well known manufacturer of oligonucleotides arrays.
 *
 *  Permission is granted to use this for non-commercial purposes and within the affy
 *  package.
 *  This was done in the hope it will be usefull, not more. There is no warranty of any kind.
 *
 *
 * In the begining, Laurent (laurent@cbs.dtu.dk 2001)
 *
 * ...fixes were made since that time, some of them by biocore members...
 *
 *  Aug 7, 2003 - Increased SIZE_LINE to 1024 (bmb)
 *
 */

/** --- includes --- */
#include <R.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef Macintosh
#include <unistd.h>
#endif

/** --- defines --- */
#ifdef Macintosh
#define index strchr
#endif

#if (defined(WIN32) || defined(Macintosh))
#define index strchr
#define rindex strrchr
#endif

/** --- contants --- */
#define SIZE_CHAR sizeof(char)
#define SIZE_LINE  1024        /* was 160 */ /* line buffer length...  too short ? YES!!! for dat header*/

/** Here for debugging */
/* Uncomment the line below for a (very) verbose output during parsing */
/* #define DEBUG */

/** Here for debugging */
/* #define HAVE_ZLIB */ 

#if defined(HAVE_ZLIB)
# include <zlib.h>
#endif


typedef struct {
  char *filepath;
  int lineno;     /* current position in the file */
  int compress;   /* compressed file (or not ?) */
  void *stream;   /* can be a FILE or a gzFile*/
} affy_file;
/* structure as handle on an Affymetrix file (currently CDF) */



/******************EXPORT*****************/
SEXP getInfo(SEXP filename, SEXP filetype, SEXP unitR, SEXP propertyR, SEXP compressR);
/*
 * get misc. info. from a CDF file
 */

SEXP readQC(SEXP filename, SEXP startunitR, SEXP indexvalR, SEXP compressR);
/* this one is experimental */
/*****************************************/


int static openCDFfile(affy_file *affyFile, char *buffy);
/**
 * Open a CDF file (making sure it looks like a CDF file)
 */

int static goToUnit(const char *unit, affy_file *affyFile, char *retour);
/*
 A unit is a stuff between '[' and ']'
 PARAMETERS:
 -- unit: the name of the unit to go to
 -- input_stream: a FILE
 -- buffy: a character buffer
 RETURNS:
 -- 1 found
 -- 0 not found
 */


int static goToUnitStartingWith(const char *unit, affy_file *affyFile, 
				char *retour);
/**
 * Same as 'goToUnit' but with partial match of the unit name
 */

char static *getProperty(const char *unit, affy_file *affyFile, char *retour);
/**
 * Return the property or NULL if not found
 */

void static close_affy_file(affy_file *affyFile);

char static *readline_affy_file(char *retour, int sizeLine, affy_file *affyFile);

int static fillCharFromLine(char *retour, int nCols, int index, SEXP cdffile,
			    char *buffer);
/**
 * Fills the data structures
 * PARAMETER(s):
 * -- retour: a char pointer (a line in the input file)
 * -- nCols: number of columns in the matrix (cf cdffile)
 * -- index: column index number in the the input file
 * -- cdffile: a matrix (R object)
 * -- buffer: a char* buffer
 */

/*
 * The function called from R 
 * - filename is a character with the path to the file
 * - indexR is an integer to indicate what should be returned
 * - compressR is an integer to indicate whether the file is compressed or not
 */

SEXP readCDFfile(SEXP filename, SEXP indexR, SEXP compressR) {
  SEXP cdffile,dim;
  affy_file affyFile;
  
  const char *param_unit = "Chip";
  const char *start_unit = "Unit";
  char *retour;
  char *buffer;
  int numUnits;
  int numBlocks;
  int numCells;
  int nRows;
  int nCols;
  int tmp;
  int ii; /* iterator across units */
  int jj; /* iterator across blocks (in a given 'unit') */
  int kk; /* iterator across cells (in a given 'block') */
  int indexVal; /* */
  
  affyFile.compress = INTEGER(compressR)[0];
  affyFile.filepath = CHAR(STRING_ELT(filename,0));
  indexVal = INTEGER(indexR)[0];
  
  retour = R_alloc(SIZE_LINE, SIZE_CHAR);
  buffer = R_alloc(SIZE_LINE, SIZE_CHAR);
  
  tmp = openCDFfile(&affyFile, retour);
  
  if (tmp == 0) {
    
    close_affy_file(&affyFile);
    error("The file %s does not appear to be a CDF file.", affyFile.filepath);
  }
  if (tmp == -1) {
    
    error("Cannot open the file %s.", affyFile.filepath);
  }
  tmp = goToUnit(param_unit, &affyFile, retour);
  /*DEBUG: should test tmp*/
  nCols = atoi(getProperty("Rows", &affyFile, retour));
  nRows = atoi(getProperty("Cols", &affyFile, retour));
  numUnits = atoi(getProperty("NumberOfUnits", &affyFile, retour));
  
  
  PROTECT(cdffile = NEW_CHARACTER(nRows*nCols));
  PROTECT(dim = NEW_INTEGER(2));
  INTEGER_POINTER(dim)[0] = nRows;
  
  INTEGER_POINTER(dim)[1] = nCols;
  SET_DIM(cdffile,dim);

  for (ii=0; ii<numUnits; ii++) {
    goToUnitStartingWith(start_unit, &affyFile, retour);
    /* Weirdness of the CDF file format... <sigh> 
     */
    /* NumCells in the Unit */
    retour = getProperty("NumCells", &affyFile, retour);
    
    /* number of blocks in the unit */
    retour = getProperty("NumberBlocks", &affyFile, retour);

    if (retour == NULL) {
	close_affy_file(&affyFile);
	
	UNPROTECT(2);
	error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
	break;
	
    }
    
    numBlocks = atoi(retour);

    for (jj=0; jj<numBlocks; jj++) {
      
      /* NumCells in the Block */
      retour = getProperty("NumCells", &affyFile, retour);
      
      if (retour == NULL) {
	close_affy_file(&affyFile);
	
	UNPROTECT(2);
	error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
	break;
	
      }
      
      numCells = atoi(retour);
      
      retour = getProperty("CellHeader", &affyFile, retour);
      
      /** not needed any longer ? */
      if (retour == NULL) {
	
	close_affy_file(&affyFile);
	UNPROTECT(2);
	
	error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
	break;
      }
    
      for (kk=0; kk<numCells; kk++) {
	readline_affy_file(retour, SIZE_LINE, &affyFile);
	
	/* CDF format madness  */
	if (retour == NULL) {
	  close_affy_file(&affyFile);
	  
	  UNPROTECT(2);
	  error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
	  break;
	}
	if (strlen(retour) <= 1) {
	  readline_affy_file(retour, SIZE_LINE, &affyFile);
	}
	if (fillCharFromLine(retour, nCols, indexVal, cdffile, buffer) == 0) {
	  
	  close_affy_file(&affyFile);
	  error("File %s corrupted (may be around line %i).", affyFile.filepath, affyFile.lineno);
	}
      }
    }
  }
  
  close_affy_file(&affyFile); 
  UNPROTECT(2);
  return cdffile;
}


SEXP getInfo(SEXP filename, SEXP filetype, SEXP unitR, SEXP propertyR, SEXP compressR) {
  SEXP value;
  affy_file affyFile;
  int tmp;
  char *retour;
  char *unit;
  char *property;

  unit = CHAR(STRING_ELT(unitR,0));
  property = CHAR(STRING_ELT(propertyR,0));
  
  retour = R_alloc(SIZE_LINE, SIZE_CHAR);
  
  affyFile.compress = INTEGER(compressR)[0];
  affyFile.filepath = CHAR(STRING_ELT(filename,0));

  tmp = 0;
  
  if (strncmp(CHAR(STRING_ELT(filetype,0)), "CDF", 2) == 0) {
    tmp = openCDFfile(&affyFile, retour);
    if (tmp == 0) {
      close_affy_file(&affyFile);
      error("The file %s does not appear to be a CDF file.", affyFile.filepath);
    }
  }
  
  if (tmp == 0) {
    error("Unknown filetype !");
  }
  
  if (tmp == -1) {
    error("Cannot open the file %s.", affyFile.filepath);
  }
  
  tmp = goToUnit(unit, &affyFile, retour);
  
  if (tmp == 0) {
    error("Unit %s not found !", unit);
  }
  
  retour = getProperty(property, &affyFile, retour);
  
  PROTECT(value = NEW_CHARACTER(1));
  SET_STRING_ELT(value, 0, mkChar(retour));
  UNPROTECT(1);
  return(value);
}

int static openCDFfile(affy_file *affyFile, char *buffy) {
  int ok;
  char *cdffileHEADER;
  
  
  cdffileHEADER = "[CDF]";
  ok = 0;

  if (affyFile->compress == 1) {
#if defined HAVE_ZLIB
    affyFile->stream = gzopen(affyFile->filepath, "rb");
    
    if (affyFile->stream == NULL) {
      ok = -1;
    } else {
      
      gzgets(affyFile->stream, buffy, SIZE_LINE);
      if (strncmp(cdffileHEADER, buffy, 4) == 0) {
	ok = 1;
	
	gzrewind(affyFile->stream);
      }
    }
#else
    
error("Compression not supported on your system.");
#endif
  } else {    
    affyFile->stream = fopen(affyFile->filepath, "r");
    
    if (affyFile->stream == NULL) {
      ok = -1;
    } else {
      
      fgets(buffy, SIZE_LINE, affyFile->stream);
      if (strncmp(cdffileHEADER, buffy, 4) == 0) {
	ok = 1;
	
	rewind(affyFile->stream);
      }
    }
  }
  affyFile->lineno = 0;
  return ok;
}


void static close_affy_file(affy_file *affyFile) {
  if (affyFile->compress == 1) {
#if defined HAVE_ZLIB
    gzclose(affyFile->stream);
#endif
  } else {
    
    fclose(affyFile->stream);
  }
}
			      
/* Wrapper for fgets: read the next line from 'affyfile' and write into the buffer 
   provided by 'retour'. At most 'sizeLine' characters are written.
   On EOF, this function returns NULL; otherwise, the value of 'retour'.  */
char static *readline_affy_file(char *retour, int sizeLine, affy_file *affyFile) {
  char *res;
  
  if (affyFile->compress == 1) {
#if defined HAVE_ZLIB
    res = gzgets(affyFile->stream, retour, sizeLine);    
#else
    error("No compression library !");
#endif
  } else {
    
    res = fgets(retour, SIZE_LINE, affyFile->stream);
  }
  /* DEBUG: test EOF */
  affyFile->lineno++;

#if defined DEBUG
  printf("%i : %s", affyFile->lineno, retour);
#endif
  return(res);
}

int static fillCharFromLine(char *retour, int nCols, int indexVal, SEXP cdffile, char *buffer) {
  int x;
  int y;
  int ii;
  int stringLen;

  /* Les voies du format sont impenetrables (sic).... */
  retour = (char *)index(retour,'=') + 1;
  x = atoi(retour);
 
  retour = (char *)index(retour,'\t') + 1;
  y = atoi(retour);
  
  for (ii=0; ii<indexVal; ii++) {
    retour = (char *)index(retour,'\t') + 1;
  }
  
  stringLen = strcspn(retour, "\t");

  buffer = strncpy(buffer, retour, stringLen);
  /* to end the line...*/
  buffer[stringLen] = '\0';
  
  
  /*SET_ELEMENT(celfile, x + nCols * y, val); */
  
  /*NUMERIC_POINTER(celfile)[x + nCols * y] = val;*/
  
  SET_STRING_ELT(cdffile, x + nCols * y, mkChar(buffer));
  return 1;
}


int static goToUnit(const char *unit, affy_file *affyFile, char *buffy) {
  int size;
  int found;
  
  char *search_string;
  search_string = R_alloc( strlen(unit)+2, SIZE_CHAR);
  search_string[0]= '\0';
  strcat(search_string,"[");
  strcat(search_string,unit);
  strcat(search_string,"]");
  size = strlen(search_string);

  found = 0;

# if defined DEBUG
  printf("Going to: [%s] \n", unit);
# endif

  buffy = readline_affy_file(buffy, SIZE_LINE, affyFile);

  while(buffy != NULL) {
    if (strncmp(search_string, buffy, size) == 0) {
      found = 1;
      break;
    }
    buffy = readline_affy_file(buffy, SIZE_LINE, affyFile);
  }

/*   /\* for now the error is fired from here *\/ */
/*   if (buffy == NULL) { */
/*     close_affy_file(affyFile);     */
/*     error("File %s is corrupted\n(Cannot find '%s')", affyFile->filepath, search_string); */
/*   } */

# if defined DEBUG
  printf("Returning: [%i \n", found);
# endif
  
  return found;
} 

int static goToUnitStartingWith(const char *unit, affy_file *affyFile, char *buffy) {
  int size;
  int found;
 
  char *search_string;
  search_string = R_alloc( strlen(unit)+1, SIZE_CHAR);
  search_string[0]= '\0';
  strcat(search_string,"[");
  strcat(search_string,unit);
  /* DEBUG: check this really a UNIT !!!!!!!!*/
  /* strcat(search_string,"]"); */
  
  size = strlen(search_string);

# if defined DEBUG
  printf("Going to something starting like: [%s \n", unit);
# endif

  
  found = 0;

  buffy = readline_affy_file(buffy, SIZE_LINE, affyFile);

  while(buffy != NULL) {
    /*      if ((strlen(buffy)-2 == size) */
    /*  &&	(strncmp(search_string, buffy, size) == 0)) { */
    
    if (strncmp(search_string, buffy, size) == 0) {
      found = 1;
      break;
    }
    buffy = readline_affy_file(buffy, SIZE_LINE, affyFile);
  }  

  /* for now the error is fired from here */
  if (buffy == NULL) {
    close_affy_file(affyFile);    
    error("File %s is corrupted\n(Cannot find '%s')", affyFile->filepath, search_string);
  }
  
# if defined DEBUG
  printf("Returning: [%i \n", found);
# endif
  
  return found;
} 


char static *getProperty(const char *unit, affy_file *affyFile, char *buffy) 
{
  int size;
  char *result_string, *rval;
  char *search_string;
  
  search_string = R_alloc( strlen(unit)+1, SIZE_CHAR);
  
  result_string = R_alloc( SIZE_LINE, SIZE_CHAR);
 
  
  search_string[0]= '\0';
  strcat(search_string, unit);
  
  strcat(search_string, "=");
  
  size = strlen(search_string);
  
  /* result_string = NULL; */

# if defined DEBUG
  printf("Going to property: %s \n", unit);
# endif

  rval = readline_affy_file(buffy, SIZE_LINE, affyFile);
  
  while(rval != NULL) {
    if (strncmp(search_string, buffy, size) == 0) {
      result_string = R_alloc( SIZE_LINE-size, SIZE_CHAR );
      strcpy(result_string, buffy+size);
      break;
    }
    rval = readline_affy_file(buffy, SIZE_LINE, affyFile);
    
  }
  
  /* for now the error is fired from here */
  if (rval == NULL) {
    close_affy_file(affyFile);    
    error("File %s is corrupted\n(Cannot find '%s')", 
	  affyFile->filepath, search_string);
  }
  result_string[strlen(result_string)-1] = '\0';
  return result_string;
}


/* DEBUG */

/* 
 * If more than one indexvalR it segfaults... no idea why... it seems to be exactly when
 * buffer2=buffer
 */ 
SEXP readQC(SEXP filename, SEXP startunitR, SEXP indexvalR, SEXP compressR) {
  const char *param_unit = "Chip";
   
  SEXP qcindex,dim;
  affy_file affyFile;
  
  char *start_unit; /* a unit to find in the CDF */
  char *retour;     /* */
  char *buffer;
  char *buffer2;
  char *field;
  
  int numVals;
  int numCells;
  int typeQC;
  int tmp;
  int ii, jj, kk;
  int stringLen;

  int x, y;
  int indexval;

  affyFile.filepath = CHAR(STRING_ELT(filename,0));  
  start_unit = CHAR(STRING_ELT(startunitR,0));
  numVals = GET_LENGTH(indexvalR);
  affyFile.compress = INTEGER(compressR)[0];
  
  retour = R_alloc(SIZE_LINE, SIZE_CHAR);
  buffer = R_alloc(SIZE_LINE, SIZE_CHAR);
  buffer2 = R_alloc(SIZE_LINE, SIZE_CHAR);
  field = R_alloc(SIZE_LINE, SIZE_CHAR);

  retour[0] = '\0';
  buffer[0] = '\0';
  buffer2[0] = '\0';
  
  /* --- open the CDF + sanity checks --- */
  tmp = openCDFfile(&affyFile, retour);   
  if (tmp == 0) {
    error("The file %s does not appear to be a CDF file.", affyFile.filepath);
  }
  if (tmp == -1) {
    error("Cannot open file %s", affyFile.filepath);
  }

  /* --- got the wished unit --- */
  tmp = goToUnit(start_unit, &affyFile, retour);
  if (tmp == 0) {
    close_affy_file(&affyFile);
    error("File %s corrupted.", affyFile.filepath);
  }
  /* --- get the characteristics of the unit --- */
  typeQC = atoi(getProperty("type", &affyFile, retour));
  numCells = atoi(getProperty("NumberCells", &affyFile, retour));
  retour = getProperty("CellHeader", &affyFile, retour);
  
  
  /* --- qcindex will contain the indexes of the QC cells --- */
  PROTECT(qcindex = NEW_NUMERIC(numCells * (numVals+2)));
  for (ii=0; ii<(numCells * (numVals+2)); ii++) {
    NUMERIC_POINTER(qcindex)[ii] = (float)0.0;
  }
  
  PROTECT(dim = NEW_INTEGER(2));
  
  INTEGER_POINTER(dim)[0] = numCells;
  INTEGER_POINTER(dim)[1] = numVals+2;
  SET_DIM(qcindex,dim);
 
  /* --- loop over the number of cells in the unit --- */
  for (ii=0; ii<numCells; ii++) {
    
    readline_affy_file(retour, SIZE_LINE, &affyFile);
    
    if (retour == NULL) {
      close_affy_file(&affyFile);
      UNPROTECT(2);
      error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
      break;
    }
    
    while (strlen(retour) <= 1) {
      /* just to make sure it was not an empty line */
      readline_affy_file(retour, SIZE_LINE, &affyFile);
    }
    
    /* DEBUG */
    /* forced to use 'buffer'... can't figure out why...*/
    buffer = (char *)(index(retour,'=') + SIZE_CHAR);
    
    x = atof(buffer);
    
    buffer = (char *)(index(buffer,'\t') + SIZE_CHAR);
    
    y = atof(buffer);
    
    NUMERIC_POINTER(qcindex)[ii + (numCells) * 0] = x;
    NUMERIC_POINTER(qcindex)[ii + (numCells) * 1] = y;
    
    /* --- loop over the number of values to get --- */
    for (jj=0; jj<numVals; jj++) {
      
      buffer2 = buffer;
      indexval = INTEGER(indexvalR)[jj];
      
      /* printf("---%d-->%s",indexval,buffer); */
      for (kk=0; kk<indexval; kk++) {
	buffer2 = (char *)(index(buffer2,'\t') + SIZE_CHAR);
	
	/* printf("-->%s",buffer2); */
	if (buffer2 == NULL) {
	 
	  close_affy_file(&affyFile);
	  UNPROTECT(2);
	  
	  error("Invalid index number, Sir...");
	}
      }
      
      stringLen = strcspn(buffer2, "\t");
      
      field = strncpy(field, buffer2, stringLen);
      /* to end the line...*/
   
      buffer2[stringLen] = '\0';
      /* printf("%s\n",field); */
      
   
      NUMERIC_POINTER(qcindex)[ii + (numCells) * (jj+2)] = atof(field);
      
    }
  }
  
  close_affy_file(&affyFile);
  
  UNPROTECT(2);
  return qcindex;  
}

