/****************************************************************
 **
 ** File: read_cdffile2.c
 **
 ** Implementation by: B. M. Bolstad
 **
 ** A parser designed to read text CDF files into an R List structure
 **
 ** Note this version only parses GC3.0 version text files (which should
 ** be almost all text CDF files currently used)
 **
 ** Note that the original text CDF parser (from which this file is not in
 ** anyway based) was written by Laurent Gautier. That file was named
 ** read_cdffile.c
 **
 ** Implemented based on documentation available from Affymetrix
 **
 ** Implementation begun 2005.
 **
 ** Modification Dates
 ** Jul 24 - Initial version
 **
 **
 **
 **
 *******************************************************************/

#include <R.h>
#include <Rdefines.h>

#include "stdlib.h"
#include "stdio.h"


#define BUFFER_SIZE 1024


/***

A structure for holding information in the 
"CDF" and "Chip" sections (basically header information)

***/



typedef struct {

  char *version;
  char *name;
  int rows,cols;
  int numberofunits;
  int maxunit;
  int NumQCUnits;
  char *chipreference;
} cdf_text_header;



/*** A structure for holding QC probe information
 *** Note the "CYCLES" item is ignored and never parsed
 ***/

typedef struct {
  int x;
  int y;
  char *probe;
  int plen;
  int atom;
  int index;
  int match;
  int bg;
} cdf_text_qc_probe;







/** A structure for holding QC information */



typedef struct{
  unsigned short type;
  unsigned int n_probes;

  cdf_text_qc_probe *qc_probes;
  
} cdf_text_qc_unit;


/*** A structure for holding probe information for unit_blocks **/



typedef struct{
  int x;
  int y;
  char *probe;
  char *feat;
  char *qual;
  int expos;
  int pos;
  char cbase;
  char pbase;
  char tbase;
  int atom;
  int index;
  int codonid;
  int codon;
  int regiontype;
  char* region;
} cdf_text_unit_block_probe;




/***

A structure holding Unit_blocks

***/

typedef struct{
  char *name;
  int blocknumber;
  int num_atoms;
  int num_cells;
  int start_position;
  int stop_position;
  int direction;
  cdf_text_unit_block_probe *probes;

} cdf_text_unit_block;






/***

A structure for holding "Units" AKA known as probesets

***/


typedef struct{
  char *name;
  int direction;
  int num_atoms;
  int num_cells;
  int unit_number;
  int unit_type;
  int numberblocks;
  int MutationType;
  cdf_text_unit_block *blocks;
} cdf_text_unit;



/***

 A structure for holding a text CDF file

***/

typedef struct{
  cdf_text_header header;
  cdf_text_qc_unit *qc_units;
  cdf_text_unit *units;
} cdf_text;


/**************************************************************
 **
 ** The following code is for tokenizing strings 
 ** originally included in read_abatch.c from the affy package.
 **
 *************************************************************/

/***************************************************************
 **
 ** tokenset
 ** 
 ** char **tokens  - a array of token strings
 ** int n - number of tokens in this set.
 **
 ** a structure to hold a set of tokens. Typically a tokenset is
 ** created by breaking a character string based upon a set of 
 ** delimiters.
 **
 **
 **************************************************************/

typedef struct{
  char **tokens;
  int n;
} tokenset;



/******************************************************************
 **
 ** tokenset *tokenize(char *str, char *delimiters)
 **
 ** char *str - a string to break into tokens
 ** char *delimiters - delimiters to use in breaking up the line
 **
 **
 ** RETURNS a new tokenset
 **
 ** Given a string, split into tokens based on a set of delimitors
 **
 *****************************************************************/

static tokenset *tokenize(char *str, char *delimiters){

  int i=0;

  char *current_token;
  tokenset *my_tokenset = Calloc(1,tokenset);
  my_tokenset->n=0;
  
  my_tokenset->tokens = NULL;

  current_token = strtok(str,delimiters);
  while (current_token != NULL){
    my_tokenset->n++;
    my_tokenset->tokens = Realloc(my_tokenset->tokens,my_tokenset->n,char*);
    my_tokenset->tokens[i] = Calloc(strlen(current_token)+1,char);
    strcpy(my_tokenset->tokens[i],current_token);
    i++;
    current_token = strtok(NULL,delimiters);
  }

  return my_tokenset; 
}


/******************************************************************
 **
 ** int tokenset_size(tokenset *x)
 **
 ** tokenset *x - a tokenset
 ** 
 ** RETURNS the number of tokens in the tokenset 
 **
 ******************************************************************/

static int tokenset_size(tokenset *x){
  return x->n;
}


/******************************************************************
 **
 ** char *get_token(tokenset *x, int i)
 **
 ** tokenset *x - a tokenset
 ** int i - index of the token to return
 ** 
 ** RETURNS pointer to the i'th token
 **
 ******************************************************************/

static char *get_token(tokenset *x,int i){
  return x->tokens[i];
}

/******************************************************************
 **
 ** void delete_tokens(tokenset *x)
 **
 ** tokenset *x - a tokenset
 ** 
 ** Deallocates all the space allocated for a tokenset 
 **
 ******************************************************************/

static void delete_tokens(tokenset *x){
  
  int i;

  for (i=0; i < x->n; i++){
    Free(x->tokens[i]);
  }
  Free(x->tokens);
  Free(x);
}

/*******************************************************************
 **
 ** int token_ends_with(char *token, char *ends)
 ** 
 ** char *token  -  a string to check
 ** char *ends_in   - we are looking for this string at the end of token
 **
 **
 ** returns  0 if no match, otherwise it returns the index of the first character
 ** which matchs the start of *ends.
 **
 ** Note that there must be one additional character in "token" beyond 
 ** the characters in "ends". So
 **
 **  *token = "TestStr"
 **  *ends = "TestStr"   
 **  
 ** would return 0 but if 
 **
 ** ends = "estStr"
 **
 ** we would return 1.
 **
 ** and if 
 ** 
 ** ends= "stStr"
 ** we would return 2 .....etc
 **
 **
 ******************************************************************/

static int token_ends_with(char *token, char *ends_in){
  
  int tokenlength = strlen(token);
  int ends_length = strlen(ends_in);
  int start_pos;
  char *tmp_ptr;
  
  if (tokenlength <= ends_length){
    /* token string is too short so can't possibly end with ends */
    return 0;
  }
  
  start_pos = tokenlength - ends_length;
  
  tmp_ptr = &token[start_pos];

  if (strcmp(tmp_ptr,ends_in)==0){
    return start_pos;
  } else {
    return 0;
  }
}


/******************************************************************
 **
 ** The following code, also from read_abatch.c is more about locating
 ** sections in the file and reading it in.
 **
 ******************************************************************/


/**
 ** This reads a line from the specified file stream
 **
 **
 **/


static void ReadFileLine(char *buffer, int buffersize, FILE *currentFile){
  if (fgets(buffer, buffersize, currentFile) == NULL){
    error("End of file reached unexpectedly. Perhaps this file is truncated.\n");
  }
}



/******************************************************************
 **
 ** void findStartsWith(FILE *my_file,char *starts, char *buffer)
 **
 ** FILE *my_file - an open file to read from
 ** char *starts - the string to search for at the start of each line
 ** char *buffer - where to place the line that has been read.
 **
 **
 ** Find a line that starts with the specified character string.
 ** At exit buffer should contain that line
 **
 *****************************************************************/


static void  findStartsWith(FILE *my_file,char *starts, char *buffer){

  int starts_len = strlen(starts);
  int match = 1;

  do {
    ReadFileLine(buffer, BUFFER_SIZE, my_file);
    match = strncmp(starts, buffer, starts_len);
  } while (match != 0);
}


/******************************************************************
 **
 ** void AdvanceToSection(FILE *my_file,char *sectiontitle, char *buffer)
 **
 ** FILE *my_file - an open file
 ** char *sectiontitle - string we are searching for
 ** char *buffer - return's with line starting with sectiontitle
 **
 **
 *****************************************************************/

static void AdvanceToSection(FILE *my_file,char *sectiontitle, char *buffer){
  findStartsWith(my_file,sectiontitle,buffer);
}










/*
typedef struct {

  char *version;
  char *name;
  int rows,cols;
  int numberofunits;
  int maxunit;
  int NumQCUnits;
  char *chipreference;
} cdf_text_header;


*/


void read_cdf_header(FILE *infile,  cdf_text *mycdf, char* linebuffer){

  tokenset *cur_tokenset;

  /* move to the Chip section */
  AdvanceToSection(infile,"[Chip]",linebuffer);
 
  findStartsWith(infile,"Name",linebuffer);
  
  /* Read the Name */
  cur_tokenset = tokenize(linebuffer,"=\r\n");
  mycdf->header.name = Calloc(sizeof(get_token(cur_tokenset,1))+1,char);
  strcpy(mycdf->header.name,get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  /* Read the Rows and Cols, Number of units etc  */

  findStartsWith(infile,"Rows",linebuffer);  
  cur_tokenset = tokenize(linebuffer,"=");
  mycdf->header.rows = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(infile,"Cols",linebuffer);
  cur_tokenset = tokenize(linebuffer,"=");
  mycdf->header.cols = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(infile,"NumberOfUnits",linebuffer);
  cur_tokenset = tokenize(linebuffer,"=");
  mycdf->header.numberofunits = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(infile,"MaxUnit",linebuffer);
  cur_tokenset = tokenize(linebuffer,"=");
  mycdf->header.maxunit = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);

  findStartsWith(infile,"NumQCUnits",linebuffer);
  cur_tokenset = tokenize(linebuffer,"=");
  mycdf->header.NumQCUnits = atoi(get_token(cur_tokenset,1));
  delete_tokens(cur_tokenset);
  
  findStartsWith(infile,"ChipReference",linebuffer);
  cur_tokenset = tokenize(linebuffer,"=\r\n");
  if (cur_tokenset->n > 1){
    mycdf->header.chipreference = Calloc(sizeof(get_token(cur_tokenset,1))+1,char);
    strcpy(mycdf->header.chipreference,get_token(cur_tokenset,1));
  } else {
    mycdf->header.chipreference = NULL;
  }


  delete_tokens(cur_tokenset);



}




void read_cdf_QCUnits(FILE *infile,  cdf_text *mycdf, char* linebuffer){
  
  tokenset *cur_tokenset;
  int i;

  mycdf->qc_units = Calloc(mycdf->header.NumQCUnits,cdf_text_qc_unit);


  for (i =0; i < mycdf->header.NumQCUnits; i++){
    /* move to the next QC section */
    AdvanceToSection(infile,"[QC",linebuffer);
    

  }







}




int read_cdf_text(char *filename, cdf_text *mycdf){

  FILE *infile;
  int i;
  char linebuffer[BUFFER_SIZE];  // a character buffer
  tokenset *cur_tokenset;
  
  if ((infile = fopen(filename, "r")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }
  


  /* Check that is is a text CDF file */
  ReadFileLine(linebuffer, BUFFER_SIZE, infile);
  if (strncmp("[CDF]", linebuffer, 5) != 0){
      error("The file %s does not look like a text CDF file",filename);
  }
  
  /* Read the version number */
  ReadFileLine(linebuffer, BUFFER_SIZE, infile);

  cur_tokenset = tokenize(linebuffer,"=\r\n");
  if (strncmp("GC3.0", get_token(cur_tokenset,1), 5) != 0){
    error("The file %s does not look like a version GC3.0 CDF file",filename);
  } else {
    mycdf->header.version = Calloc(sizeof(get_token(cur_tokenset,1))+1,char);
    strcpy(mycdf->header.version,get_token(cur_tokenset,1));
  }
  delete_tokens(cur_tokenset);


  read_cdf_header(infile,mycdf,linebuffer);
  read_cdf_QCUnits(infile,mycdf,linebuffer);



  return 1;
}





/****
     
 this function should be called from R. When supplied the name
 of a text cdf file it first parses it into a C data structure.

 An R list structure is then constructed from the C data structure

 The R list is then returned.

 Note no special effort is made to reduce down the information in
 the text CDF file. Instead almost everything is returned, even somewhat redundant information.

****/

SEXP ReadtextCDFFileIntoRList(SEXP filename){

  SEXP CDFInfo;  /* this is the object that will be returned */
  SEXP CDFInfoNames;
  SEXP HEADER;  /* The file header */
  SEXP HEADERNames;
  SEXP TEMPSXP;
						

  cdf_text my_cdf;

  char *cur_file_name;
  cur_file_name = CHAR(VECTOR_ELT(filename,0));

  if(!read_cdf_text(cur_file_name, &my_cdf)){
    error("Problem reading text cdf file %s. Possibly corrupted or truncated?\n",cur_file_name);
  }


  /* Now build the R list structure */


   /* return the full structure */
  PROTECT(CDFInfo = allocVector(VECSXP,1));
  PROTECT(CDFInfoNames = allocVector(STRSXP,1));
  SET_VECTOR_ELT(HEADERNames,0,mkChar("Chip"));

  setAttrib(CDFInfo,R_NamesSymbol,CDFInfoNames);
  UNPROTECT(1);

  /* Deal with the HEADER */
  PROTECT(HEADER = allocVector(VECSXP,8));
  PROTECT(HEADERNames = allocVector(STRSXP,8));
  SET_VECTOR_ELT(HEADERNames,0,mkChar("Version"));
  SET_VECTOR_ELT(HEADERNames,1,mkChar("Name"));
  SET_VECTOR_ELT(HEADERNames,2,mkChar("Rows"));
  SET_VECTOR_ELT(HEADERNames,3,mkChar("Cols"));
  SET_VECTOR_ELT(HEADERNames,4,mkChar("NumberOfUnits"));
  SET_VECTOR_ELT(HEADERNames,5,mkChar("MaxUnit"));
  SET_VECTOR_ELT(HEADERNames,6,mkChar("NumQCUnits"));
  SET_VECTOR_ELT(HEADERNames,7,mkChar("ChipReference"));
  setAttrib(HEADER,R_NamesSymbol,HEADERNames);
  UNPROTECT(1);
  
  PROTECT(TEMPSXP = allocVector(STRSXP,1));
  SET_VECTOR_ELT(TEMPSXP,0,mkChar(my_cdf.header.version));
  SET_VECTOR_ELT(HEADER,0,TEMPSXP); 
  UNPROTECT(1);
  
  PROTECT(TEMPSXP = allocVector(STRSXP,1));
  SET_VECTOR_ELT(TEMPSXP,0,mkChar(my_cdf.header.name));
  SET_VECTOR_ELT(HEADER,1,TEMPSXP); 
  UNPROTECT(1);

  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  NUMERIC_POINTER(TEMPSXP)[0] = (double)my_cdf.header.rows;
  SET_VECTOR_ELT(HEADER,2,TEMPSXP);
  UNPROTECT(1);
  
  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  NUMERIC_POINTER(TEMPSXP)[0] = (double)my_cdf.header.cols;
  SET_VECTOR_ELT(HEADER,3,TEMPSXP);
  UNPROTECT(1);
 
  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  NUMERIC_POINTER(TEMPSXP)[0] = (double)my_cdf.header.numberofunits;
  SET_VECTOR_ELT(HEADER,4,TEMPSXP);
  UNPROTECT(1);

  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  NUMERIC_POINTER(TEMPSXP)[0] = (double)my_cdf.header.maxunit;
  SET_VECTOR_ELT(HEADER,5,TEMPSXP);
  UNPROTECT(1);

  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  NUMERIC_POINTER(TEMPSXP)[0] = (double)my_cdf.header.NumQCUnits;
  SET_VECTOR_ELT(HEADER,6,TEMPSXP);
  UNPROTECT(1);
  
  PROTECT(TEMPSXP = allocVector(REALSXP,1));
  if (my_cdf.header.chipreference !=NULL){
    SET_VECTOR_ELT(TEMPSXP,0,mkChar(my_cdf.header.chipreference));
    SET_VECTOR_ELT(HEADER,7,TEMPSXP); 
  }
  UNPROTECT(1);
  
  SET_VECTOR_ELT(CDFInfo,0,HEADER);

  

  UNPROTECT(2);
  return CDFInfo;
}





