/**
 *  These C routines read data held CDF and CEL files, two of the data formats
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
/* structure as handle on an Affymetrix file (currently CEL or CDF) */



/******************EXPORT*****************/
SEXP readCELfile(SEXP filename, SEXP indexR, SEXP compressR);
/*
 * - filename is a character with the path to the file
 * - indexR is an integer to indicate what should be returned
 * - compressR is an integer to indicate whether the file is compressed or not
 */

SEXP getIndexExtraFromCEL (SEXP filename, SEXP unitR, SEXP compressR);
/*
 * get the masks and outliers 
 */

SEXP readCDFfile(SEXP filename, SEXP indexR, SEXP compressR);

SEXP getInfo(SEXP filename, SEXP filetype, SEXP unitR, SEXP propertyR, SEXP compressR);
/*
 * get misc. info. from a CEL or CDF file
 */

SEXP readQC(SEXP filename, SEXP startunitR, SEXP indexvalR, SEXP compressR);
/* this one is experimental */
/*****************************************/


int static openCDFfile(affy_file *affyFile, const char *mode, char *buffy);
/**
 * Open a CDF file (making sure it looks like a CDF file)
 */

int static openCELfile(affy_file *affyFile, const char *mode, char *buxffy);
/**
 * Open a CEL file (making sure it is a CEL file)
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

int static fillFromLine(char *retour, int nCols, int index, SEXP celfile);
/**
 * Fills the data structures
 * PARAMETER(s):
 * -- retour: a char pointer (a line read from input file)
 * -- nCols: number of columns in the matrix (cf celfile)
 * -- index: index for the column to use in the the input file
 * -- celfile: a matrix of mode 'numeric' (R object)
 */

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

SEXP readCELfile(SEXP filename, SEXP indexR, SEXP compressR) {
  SEXP celfile,dim;

  affy_file affyFile; /* handle structure */
  
  const char *mode = "r";
  const char *param_unit = "HEADER";
  const char *start_unit = "INTENSITY";
  const char *axis_invert_X_prop = "Axis-invertX";
  const char *axis_invert_Y_prop = "AxisInvertY";
  const char *swap_XY_prop = "swapXY";
  const char *numCells_prop = "NumberCells";
  
  char *retour; /* buffer for each line read */
  int numCells; /* number of cells in the CEL */
  int nRows;
  int nCols;
  int tmp;
  int ii;       /* iterator through cells */
  int indexVal; /* in which column are the values to return */
  
  affyFile.compress = INTEGER(compressR)[0];
  affyFile.filepath = CHAR(STRING_ELT(filename,0));
  
  indexVal = INTEGER(indexR)[0];
  
  retour = R_alloc(SIZE_LINE, SIZE_CHAR);
  
  tmp = openCELfile(&affyFile, mode, retour);
  
  if (tmp == 0) {
    close_affy_file(&affyFile);
    error("The file %s does not appear to be a CEL file.", affyFile.filepath);
  }
  
  if (tmp == -1) {
    error("Cannot open the file %s", affyFile.filepath);
  }
  
  tmp = goToUnit(param_unit, &affyFile, retour);
  if (tmp == 0) {
    close_affy_file(&affyFile);
    error("File %s corrupted.", affyFile.filepath);
  }
  /* useless for the moment... cross-check with numCells
   * under investigation... */
  nCols = atoi(getProperty("Cols", &affyFile, retour));
  nRows = atoi(getProperty("Rows", &affyFile, retour));
  
  /* scary properties in CEL files... */
  if (strncmp("0", getProperty(axis_invert_X_prop, &affyFile, retour), 1) != 0) {
    error("inverted X axis (not handled )!");
  }
  if (strncmp("0", getProperty(axis_invert_Y_prop, &affyFile, retour), 1) != 0) {
    error("inverted Y axis (not handled )!");
  }
  if (strncmp("0", getProperty(swap_XY_prop, &affyFile, retour), 1) != 0) {
    error("swap XY (not handled )!");
  }

  tmp = goToUnit(start_unit, &affyFile, retour);
  
  if (tmp == 0) {
    close_affy_file(&affyFile);
    error("File %s corrupted.", affyFile.filepath);
  }
  numCells = atoi(getProperty(numCells_prop, &affyFile, retour));
  
  PROTECT(celfile = NEW_NUMERIC(nRows * nCols));
  PROTECT(dim = NEW_INTEGER(2));
  INTEGER_POINTER(dim)[0] = nRows;
  
  INTEGER_POINTER(dim)[1] = nCols;
  SET_DIM(celfile, dim);
  
  /*NOTE: headers are here...*/
  readline_affy_file(retour, SIZE_LINE, &affyFile);
  
  for (ii=0; ii<numCells; ii++) {
    
    readline_affy_file(retour, SIZE_LINE, &affyFile);
    if (retour == NULL) {
      close_affy_file(&affyFile);
      UNPROTECT(2);
      
      error("Unexpected and premature end of the file %s\n(truncated CEL file ?).", 
	    affyFile.filepath);
      break;
    }
    if (fillFromLine(retour, nCols, indexVal, celfile) == 0) {
      
      UNPROTECT(2);
      close_affy_file(&affyFile);
      error("File %s corrupted (may be around line %i).", affyFile.filepath, affyFile.lineno);
      
    }
  }
  
  close_affy_file(&affyFile);
  UNPROTECT(2);
  
  return celfile;
}


/** get the 'extra' stuff held in a CEL file 
(i.e. MASKED, OUTLIERS ...) */
SEXP getIndexExtraFromCEL (SEXP filename, SEXP unitR, SEXP compressR) {
  SEXP extraindex, dim;
  
  affy_file affyFile;
  
  const char *mode = "r";

  char *start_unit;
  char *retour;
  char *buffer;

  int tmp;
  int ii,x,y;
  int numCells;

  affyFile.filepath = CHAR(STRING_ELT(filename,0));

  start_unit = CHAR(STRING_ELT(unitR,0));

  affyFile.compress = INTEGER(compressR)[0];
  
  retour = R_alloc(SIZE_LINE, SIZE_CHAR);
  buffer = R_alloc(SIZE_LINE, SIZE_CHAR);
  
  retour[0] = '\0';
  
  buffer[0] = '\0';
  
  tmp = openCELfile(&affyFile, mode, retour);
  
  if (tmp == 0) {
    close_affy_file(&affyFile);
 
   error("The file %s does not appear to be a CEL file.", 
	 affyFile.filepath);
  }
  
  if (tmp == -1) {
    
    error("Cannot open the file %s.", affyFile.filepath);
  }
  
  tmp = goToUnit(start_unit, &affyFile, retour);
  if (tmp == 0) {
    /* close_affy_file(&affyFile); */
    warning("File %s has no unit '%s'.",  affyFile.filepath, start_unit);
    numCells = 0;
  } else {
    numCells = atoi(getProperty("NumberCells", &affyFile, retour));
  }
  
  
  PROTECT(dim = NEW_INTEGER(2));
  if (numCells == 0) {
    /* DEBUG: I could not figure out how to make a matrix of dim (0,0)
     * neither how to set NULL to the slot in Cel... the hack belows
     * could be very dependant on the R version... */

    PROTECT(extraindex = NEW_LOGICAL(1));
    LOGICAL_POINTER(extraindex)[0] = NA_LOGICAL;
    INTEGER_POINTER(dim)[0] = 1;
    INTEGER_POINTER(dim)[1] = 1;
  } else {
    /* what's the use ? */
    retour = getProperty("CellHeader", &affyFile, retour);
    
    PROTECT(extraindex = NEW_NUMERIC(numCells * (2)));
    /* init needed ?*/
    for (ii=0; ii<(numCells * (2)); ii++) {
      NUMERIC_POINTER(extraindex)[ii] = (float)0.0;
    }
    INTEGER_POINTER(dim)[0] = numCells;
    INTEGER_POINTER(dim)[1] = 2; /* only two columns to hold x and y indexes */
  } 

  SET_DIM(extraindex, dim);  

  for (ii=0; ii<numCells; ii++) {
    
    readline_affy_file(retour, SIZE_LINE, &affyFile);
    
    if (retour == NULL) {
      
      close_affy_file(&affyFile);
      UNPROTECT(2);
      
      error("Unexpected and premature end of the file %s\n (truncated CDF file ?).", affyFile.filepath);
      break;
    }
    
    if (strlen(retour) <= 1) {
      readline_affy_file(retour, SIZE_LINE, &affyFile);
    }
    
    /* DEBUG */
    /* forced to use 'buffer'... can't figure out why...*/
      
    buffer = retour;
    x = atof(buffer);
    
    buffer = (char *)index(buffer,'\t') + 1;
    
    y = atof(buffer);
    
    
    NUMERIC_POINTER(extraindex)[ii + (numCells) * 0] = x + 1; /* +1 one needed as indexing */
    NUMERIC_POINTER(extraindex)[ii + (numCells) * 1] = y + 1; /* start at 1 in R */
  }
  
  
  close_affy_file(&affyFile);

  UNPROTECT(2);

  return extraindex;  
}

/*
 * The function called from R 
 * - filename is a character with the path to the file
 * - indexR is an integer to indicate what should be returned
 * - compressR is an integer to indicate whether the file is compressed or not
 */

SEXP readCDFfile(SEXP filename, SEXP indexR, SEXP compressR) {
  SEXP cdffile,dim;
  affy_file affyFile;
  
  const char *mode = "r";
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
  
  tmp = openCDFfile(&affyFile, mode, retour);
  
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

  const char *mode = "r";

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
  
  /** "filetype" can only be "CDF" or "CEL" */
  if (strncmp(CHAR(STRING_ELT(filetype,0)), "CDF", 2) == 0) {
    tmp = openCDFfile(&affyFile, mode, retour);
    if (tmp == 0) {
      close_affy_file(&affyFile);
      error("The file %s does not appear to be a CDF file.", affyFile.filepath);
    }
  }
  
  if (strncmp(CHAR(STRING_ELT(filetype, 0)), "CEL", 2) == 0) {
    tmp = openCELfile(&affyFile, mode, retour);
    if (tmp == 0) {
      close_affy_file(&affyFile);
      error("The file %s does not appear to be a CEL file.", affyFile.filepath);
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

int static openCELfile(affy_file *affyFile, const char *mode, char *buffy) {
  int ok;
  char *celfileHEADER;
  
  celfileHEADER = "[CEL]";
  ok = 0;
  
  if (affyFile->compress == 1) {
#if defined HAVE_ZLIB
    affyFile->stream = gzopen(affyFile->filepath, mode);
    
    if (affyFile->stream == NULL) {
      ok = -1;
    } else {
      gzgets(affyFile->stream, buffy, SIZE_LINE);
      if (strncmp(celfileHEADER, buffy, 4) == 0) {
	ok = 1;
	gzrewind(affyFile->stream);
      }
    }
#else    
error("Compression not supported on your system.");
#endif
  } else { 
    affyFile->stream = fopen(affyFile->filepath, mode);
 
   if (affyFile->stream == NULL) {
      ok = -1;
    } else {
      fgets(buffy, SIZE_LINE, affyFile->stream);
    }
    if (strncmp(celfileHEADER, buffy, 4) == 0) {
      ok = 1;      
      rewind(affyFile->stream);
    }
  }
  affyFile->lineno = 0;
  
  return ok;
}

int static openCDFfile(affy_file *affyFile, const char *mode, char *buffy) {
  int ok;
  char *cdffileHEADER;
  
  
  cdffileHEADER = "[CDF]";
  ok = 0;

  if (affyFile->compress == 1) {
#if defined HAVE_ZLIB
    affyFile->stream = gzopen(affyFile->filepath, mode);
    
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
    affyFile->stream = fopen(affyFile->filepath, mode);
    
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

int static fillFromLine(char *retour, int nCols, int indexVal, SEXP celfile) {
  
  int x;
  int y;
  double val;
  
  x = atoi(retour);
  
  retour = (char *)index(retour,'\t');
  
  if (retour == NULL) 
    { return 0; }
  
  y = atoi(++retour);
  
  retour = (char *)index(retour,'\t');

  if (retour == NULL) { return 0; }
  
  if (indexVal == 2) {
    retour = (char *)index(++retour, '\t');
  }
  
  if (retour == NULL) { return 0; }
  
  val = atof(++retour);
  /*SET_ELEMENT(celfile, x + nCols * y, val); */
  
  NUMERIC_POINTER(celfile)[x + nCols * y] = val;
  return 1;
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
 
  const char *mode = "r";
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
  tmp = openCDFfile(&affyFile, mode, retour);   
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

