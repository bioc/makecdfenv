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
 ** Sep 20 - Continued Implementation
 ** Sep 21 - Continued Implementation and debugging
 ** Sep 22 - Continued Implementation and testing
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
  int type;
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
  mycdf->header.name = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
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
    mycdf->header.chipreference = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
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
    findStartsWith(infile,"Type",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->qc_units[i].type = (unsigned short)atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    findStartsWith(infile,"NumberCells",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->qc_units[i].n_probes = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    mycdf->qc_units[i].qc_probes = Calloc(mycdf->qc_units[i].n_probes,cdf_text_qc_probe);
  }







}



void read_cdf_unit_block(FILE *infile,  cdf_text *mycdf, char* linebuffer, int unit){
  tokenset *cur_tokenset;
  int i;
  

  
  for (i=0; i < mycdf->units[unit].numberblocks; i++){ 

    mycdf->units[unit].blocks[i].blocknumber;
    findStartsWith(infile,"Name",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=\r\n");
    mycdf->units[unit].blocks[i].name = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
    strcpy(mycdf->units[unit].blocks[i].name,get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    // Rprintf("%s\n",mycdf->units[unit].blocks[i].name);
    


    findStartsWith(infile,"BlockNumber",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[unit].blocks[i].blocknumber = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    //  Rprintf("%d %d %d\n",unit,i,mycdf->header.numberofunits);

    findStartsWith(infile,"NumAtoms",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[unit].blocks[i].num_atoms = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    findStartsWith(infile,"NumCells",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[unit].blocks[i].num_cells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);


    findStartsWith(infile,"StartPosition",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[unit].blocks[i].start_position = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    findStartsWith(infile,"StopPosition",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[unit].blocks[i].stop_position = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    if (mycdf->units[unit].unit_type == 2){
      findStartsWith(infile,"Direction",linebuffer);
      cur_tokenset = tokenize(linebuffer,"=");
      mycdf->units[unit].blocks[i].direction = atoi(get_token(cur_tokenset,1));
      delete_tokens(cur_tokenset); 
    } else {
      mycdf->units[unit].blocks[i].direction = mycdf->units[unit].direction;
    }
    
    mycdf->units[unit].blocks[i].probes = Calloc(mycdf->units[unit].blocks[i].num_cells,cdf_text_unit_block_probe);

  }
  



}








void read_cdf_Units(FILE *infile,  cdf_text *mycdf, char* linebuffer){
  tokenset *cur_tokenset;
  int i;

  mycdf->units = Calloc(mycdf->header.numberofunits,cdf_text_unit);

  for (i =0; i < mycdf->header.numberofunits; i++){
    /* move to the next Unit section */
    AdvanceToSection(infile,"[Unit",linebuffer);
    findStartsWith(infile,"Name",linebuffer); 
    cur_tokenset = tokenize(linebuffer,"=\r\n");
    mycdf->units[i].name = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
    strcpy(mycdf->units[i].name,get_token(cur_tokenset,1));
 
    delete_tokens(cur_tokenset);
  
    
    
    findStartsWith(infile,"Direction",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].direction = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    
    findStartsWith(infile,"NumAtoms",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].num_atoms = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);
    
    findStartsWith(infile,"NumCells",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].num_cells = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    findStartsWith(infile,"UnitNumber",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].unit_number = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    findStartsWith(infile,"UnitType",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].unit_type = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset);

    findStartsWith(infile,"NumberBlocks",linebuffer);
    cur_tokenset = tokenize(linebuffer,"=");
    mycdf->units[i].numberblocks = atoi(get_token(cur_tokenset,1));
    delete_tokens(cur_tokenset); 

    /*Skip MutationType since only appears on one type of array */
    
    mycdf->units[i].blocks = Calloc(mycdf->units[i].numberblocks,cdf_text_unit_block);

  
    read_cdf_unit_block(infile,mycdf,linebuffer,i); 
    //    AdvanceToSection(infile,"[Unit",linebuffer);
    Rprintf("%d\n",i);
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
    mycdf->header.version = Calloc(strlen(get_token(cur_tokenset,1))+1,char);
    strcpy(mycdf->header.version,get_token(cur_tokenset,1));
  }
  delete_tokens(cur_tokenset);


  read_cdf_header(infile,mycdf,linebuffer);
  read_cdf_QCUnits(infile,mycdf,linebuffer);
  read_cdf_Units(infile,mycdf,linebuffer);


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
  SEXP TEMPSXP2;
  SEXP TEMPSXP3;
  SEXP TEMPSXP4;

  SEXP QCUNITS;
  SEXP UNITS;

  int i,j;
						

  cdf_text my_cdf;

  char *cur_file_name;
  cur_file_name = CHAR(VECTOR_ELT(filename,0));

  if(!read_cdf_text(cur_file_name, &my_cdf)){
    error("Problem reading text cdf file %s. Possibly corrupted or truncated?\n",cur_file_name);
  }


  /* Now build the R list structure */


   /* return the full structure */
  PROTECT(CDFInfo = allocVector(VECSXP,3));
  PROTECT(CDFInfoNames = allocVector(STRSXP,3));
  SET_VECTOR_ELT(CDFInfoNames,0,mkChar("Chip"));
  SET_VECTOR_ELT(CDFInfoNames,1,mkChar("QC"));
  SET_VECTOR_ELT(CDFInfoNames,2,mkChar("Unit"));

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

  PROTECT(QCUNITS = allocVector(VECSXP,my_cdf.header.NumQCUnits));
  for (i=0; i <my_cdf.header.NumQCUnits; i++){
    PROTECT(TEMPSXP=allocVector(VECSXP,3));
    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.qc_units[i].type;
    SET_VECTOR_ELT(TEMPSXP,0,TEMPSXP2);
    UNPROTECT(1);
    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.qc_units[i].n_probes;
    SET_VECTOR_ELT(TEMPSXP,1,TEMPSXP2);
    UNPROTECT(1);
    PROTECT(TEMPSXP2=allocVector(STRSXP,3));
    SET_VECTOR_ELT(TEMPSXP2,0,mkChar("Type"));
    SET_VECTOR_ELT(TEMPSXP2,1,mkChar("NumberCells"));
    SET_VECTOR_ELT(TEMPSXP2,2,mkChar("QCCells"));
    setAttrib(TEMPSXP,R_NamesSymbol,TEMPSXP2);
    UNPROTECT(1);
    SET_VECTOR_ELT(QCUNITS,i,TEMPSXP);
   
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(CDFInfo,1,QCUNITS);
  UNPROTECT(1);


  PROTECT(UNITS = allocVector(VECSXP,my_cdf.header.numberofunits));
  for (i=0; i < my_cdf.header.numberofunits; i++){
    PROTECT(TEMPSXP=allocVector(VECSXP,8));
    PROTECT(TEMPSXP2=allocVector(STRSXP,1));
      
    SET_VECTOR_ELT(TEMPSXP2,0,mkChar(my_cdf.units[i].name));
    SET_VECTOR_ELT(TEMPSXP,0,TEMPSXP2);
    UNPROTECT(1);
    
    
    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].direction;
    SET_VECTOR_ELT(TEMPSXP,1,TEMPSXP2);
    UNPROTECT(1);

    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].num_atoms;
    SET_VECTOR_ELT(TEMPSXP,2,TEMPSXP2);
    UNPROTECT(1);

    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].num_cells;
    SET_VECTOR_ELT(TEMPSXP,3,TEMPSXP2);
    UNPROTECT(1);

   
    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].unit_number;
    SET_VECTOR_ELT(TEMPSXP,4,TEMPSXP2);
    UNPROTECT(1); 

    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].unit_type;
    SET_VECTOR_ELT(TEMPSXP,5,TEMPSXP2);
    UNPROTECT(1); 

    PROTECT(TEMPSXP2 = allocVector(REALSXP,1));
    NUMERIC_POINTER(TEMPSXP2)[0] = (double)my_cdf.units[i].numberblocks;
    SET_VECTOR_ELT(TEMPSXP,6,TEMPSXP2);
    UNPROTECT(1); 

    PROTECT(TEMPSXP2 = allocVector(VECSXP,my_cdf.units[i].numberblocks));

    for (j=0; j <my_cdf.units[i].numberblocks; j++){
      PROTECT(TEMPSXP3 = allocVector(VECSXP,8));

      
      PROTECT(TEMPSXP4=allocVector(STRSXP,1));
      
      SET_VECTOR_ELT(TEMPSXP4,0,mkChar(my_cdf.units[i].blocks[j].name));
      SET_VECTOR_ELT(TEMPSXP3,0,TEMPSXP4);
      UNPROTECT(1);

      
      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].blocknumber;
      SET_VECTOR_ELT(TEMPSXP3,1,TEMPSXP4);
      UNPROTECT(1);

      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].num_atoms;
      SET_VECTOR_ELT(TEMPSXP3,2,TEMPSXP4);
      UNPROTECT(1);
      
      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].num_cells;
      SET_VECTOR_ELT(TEMPSXP3,3,TEMPSXP4);
      UNPROTECT(1);
      

      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].start_position;
      SET_VECTOR_ELT(TEMPSXP3,4,TEMPSXP4);
      UNPROTECT(1);
      
      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].stop_position;
      SET_VECTOR_ELT(TEMPSXP3,5,TEMPSXP4);
      UNPROTECT(1);


      PROTECT(TEMPSXP4=allocVector(REALSXP,1));
      NUMERIC_POINTER(TEMPSXP4)[0] = (double)my_cdf.units[i].blocks[j].direction;
      SET_VECTOR_ELT(TEMPSXP3,6,TEMPSXP4);
      UNPROTECT(1);



      PROTECT(TEMPSXP4=allocVector(STRSXP,8));
      SET_VECTOR_ELT(TEMPSXP4,0,mkChar("Name"));
      SET_VECTOR_ELT(TEMPSXP4,1,mkChar("BlockNumber"));
      SET_VECTOR_ELT(TEMPSXP4,2,mkChar("NumAtoms"));
      SET_VECTOR_ELT(TEMPSXP4,3,mkChar("NumCells"));
      SET_VECTOR_ELT(TEMPSXP4,4,mkChar("StartPosition"));
      SET_VECTOR_ELT(TEMPSXP4,5,mkChar("StopPosition"));
      SET_VECTOR_ELT(TEMPSXP4,6,mkChar("Direction"));
      SET_VECTOR_ELT(TEMPSXP4,7,mkChar("Unit_Block_Cells"));
      setAttrib(TEMPSXP3,R_NamesSymbol,TEMPSXP4);
      UNPROTECT(1);

      SET_VECTOR_ELT(TEMPSXP2,j,TEMPSXP3);
      UNPROTECT(1);
    }






    SET_VECTOR_ELT(TEMPSXP,7,TEMPSXP2);
    UNPROTECT(1); 

    


    PROTECT(TEMPSXP2 = allocVector(STRSXP,8));
    SET_VECTOR_ELT(TEMPSXP2,0,mkChar("Name"));
    SET_VECTOR_ELT(TEMPSXP2,1,mkChar("Direction"));
    SET_VECTOR_ELT(TEMPSXP2,2,mkChar("NumAtoms"));
    SET_VECTOR_ELT(TEMPSXP2,3,mkChar("NumCells"));
    SET_VECTOR_ELT(TEMPSXP2,4,mkChar("UnitNumber"));
    SET_VECTOR_ELT(TEMPSXP2,5,mkChar("UnitType"));
    SET_VECTOR_ELT(TEMPSXP2,6,mkChar("NumberBlocks"));
    SET_VECTOR_ELT(TEMPSXP2,7,mkChar("Unit_Block"));
    setAttrib(TEMPSXP,R_NamesSymbol,TEMPSXP2);
    UNPROTECT(1);





    SET_VECTOR_ELT(UNITS,i,TEMPSXP);
    UNPROTECT(1);
    


  }
  SET_VECTOR_ELT(CDFInfo,2,UNITS);
  UNPROTECT(1);


  

  UNPROTECT(2);
  return CDFInfo;
}





