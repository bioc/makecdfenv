/****************************************************************
 **
 ** File: read_cdf_xda.c
 **
 ** A parser designed to read the binary format cdf files
 **
 ** Implemented based on documentation available from Affymetrix
 **
 ** Modification Dates
 ** Feb 4 - Initial version
 ** Feb 5 - A bunch of hacks for SNP chips.
 **
 ****************************************************************/

/** --- includes --- */
#include <R.h>
#include <Rdefines.h>

#include "stdlib.h"
#include "stdio.h"


//#define READ_CDF_DEBUG
//#define READ_CDF_DEBUG_SNP
#define READ_CDF_NOSNP



/************************************************************************
 **
 ** Structures for holding the CDF file information
 **
 ************************************************************************/

typedef struct {
  int magicnumber;
  int version_number;
  short rows,cols;
  int n_units,n_qc_units;
  int len_ref_seq;
  int i;
  char *ref_seq;
} cdf_xda_header;


/*   QC information, repeated for each QC unit:

Type - unsigned short
Number of probes - integer

Probe information, repeated for each probe in the QC unit:
X coordinate - unsigned short
Y coordinate - unsigned short
Probe length - unsigned char
Perfect match flag - unsigned char
Background probe flag - unsigned char
*/

typedef struct{
  unsigned short x;
  unsigned short y;
  unsigned char probelength;
  unsigned char pmflag;
  unsigned char bgprobeflag;

} cdf_qc_probe;

typedef struct{
  unsigned short type;
  unsigned int n_probes;

  cdf_qc_probe *qc_probes;
  
} cdf_qc_unit;


/* Unit information, repeated for each unit:

UnitType - unsigned short (1 - expression, 2 - genotyping, 3 - CustomSeq, 3 - tag)
Direction - unsigned char
Number of atoms - integer
Number of blocks - integer (always 1 for expression units)
Number of cells - integer
Unit number (probe set number) - integer
Number of cells per atom - unsigned char

Block information, repeated for each block in the unit:

Number of atoms - integer
Number of cells - integer
Number of cells per atom - unsigned char
Direction - unsigned char
The position of the first atom - integer
<unused integer value> - integer
The block name - char[64]

Cell information, repeated for each cell in the block:

Atom number - integer
X coordinate - unsigned short
Y coordinate - unsigned short
Index position (relative to sequence for resequencing units, for expression and mapping units this value is just the atom number) - integer
Base of probe at substitution position - char
Base of target at interrogation position - char

*/


typedef struct{
  int atomnumber;
  unsigned short x;
  unsigned short y;
  int indexpos;
  char pbase;
  char tbase;
} cdf_unit_cell;


typedef struct{
  int natoms;
  int ncells;
  unsigned char ncellperatom;
  unsigned char direction;
  int firstatom;
  int unused;         /* in the docs this is called "unused" but by the looks of it it is actually the lastatom */
  char blockname[64];

  cdf_unit_cell  *unit_cells;

} cdf_unit_block;


typedef struct{
  unsigned short unittype;
  unsigned char direction;
  int natoms;
  int nblocks;
  int ncells;
  int unitnumber;
  unsigned char ncellperatom;

  cdf_unit_block *unit_block;
  
} cdf_unit;







typedef struct {
  
  cdf_xda_header header;  /* Header information */
  char **probesetnames;   /* Names of probesets */
  
  int *qc_start;          /* These are used for random access */
  int *units_start;
  
  cdf_qc_unit *qc_units;
  cdf_unit *units;


} cdf_xda;









/*************************************************************************
 **
 ** Code for reading from the binary files, doing bit flipping if
 ** necessary (on big-endian machines)
 **
 **
 ************************************************************************/


static size_t fread_int32(int *destination, int n, FILE *instream){

  size_t result;

  result = fread(destination,sizeof(int),n,instream);

#ifdef WORDS_BIGENDIAN
  /* bit flip since all Affymetrix binary files are little endian */

  *destination=(((*destination>>24)&0xff) | ((*destination&0xff)<<24) |
                ((*destination>>8)&0xff00) | ((*destination&0xff00)<<8));
#endif
  return result;
}



static size_t fread_uint32(unsigned int *destination, int n, FILE *instream){


  size_t result;

  result = fread(destination,sizeof(unsigned int),n,instream);


#ifdef WORDS_BIGENDIAN
  /* bit flip since all Affymetrix binary files are little endian */
  *destination=(((*destination>>24)&0xff) | ((*destination&0xff)<<24) |
            ((*destination>>8)&0xff00) | ((*destination&0xff00)<<8));

#endif
  return result;
}



static size_t fread_int16(short *destination, int n, FILE *instream){
   size_t result;

   result = fread(destination,sizeof(short),n,instream);

#ifdef WORDS_BIGENDIAN
  /* bit flip since all Affymetrix binary files are little endian */
   *destination=(((*destination>>8)&0xff) | ((*destination&0xff)<<8));

#endif
   return result;

}




static size_t fread_uint16(unsigned short *destination, int n, FILE *instream){
   size_t result;

   result = fread(destination,sizeof(unsigned short),n,instream);

#ifdef WORDS_BIGENDIAN
  /* bit flip since all Affymetrix binary files are little endian */
   *destination=(((*destination>>8)&0xff) | ((*destination&0xff)<<8));

#endif
   return result;

}



static void swap_float_4(float *tnf4)              /* 4 byte floating point numbers */
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}



static size_t fread_float32(float *destination, int n, FILE *instream){

  size_t result;



  result = fread(destination,sizeof(float),n,instream);

#ifdef WORDS_BIGENDIAN
  swap_float_4(destination);
#endif

  return result;
}

static size_t fread_char(char *destination, int n, FILE *instream){

  int i=0;
  size_t result;

  result = fread(destination,sizeof(char),n,instream);

#ifdef WORDS_BIGENDIAN
  /* Probably don't need to do anything for characters */

#endif

  return result;

}

static size_t fread_uchar(unsigned char *destination, int n, FILE *instream){

  int i=0;
  size_t result;

  result = fread(destination,sizeof(unsigned char),n,instream);

#ifdef WORDS_BIGENDIAN
  /* Probably don't need to do anything for characters */
  destination = ~destination;
#endif

  return result;

}




int read_cdf_qcunit(cdf_qc_unit *my_unit,int filelocation,FILE *instream){
  
  int i,j;


  fseek(instream,filelocation,SEEK_SET);

  fread_uint16(&(my_unit->type),1,instream);
  fread_int32(&(my_unit->n_probes),1,instream);


  my_unit->qc_probes = Calloc(my_unit->n_probes,cdf_qc_probe);

  for (i=0; i < my_unit->n_probes; i++){
    fread_uint16(&(my_unit->qc_probes[i].x),1,instream);
    fread_uint16(&(my_unit->qc_probes[i].y),1,instream);
    fread_uchar(&(my_unit->qc_probes[i].probelength),1,instream);
    fread_uchar(&(my_unit->qc_probes[i].pmflag),1,instream);
    fread_uchar(&(my_unit->qc_probes[i].bgprobeflag),1,instream);
    
     }
  return 1;
}


int read_cdf_unit(cdf_unit *my_unit,int filelocation,FILE *instream){

  int i,j;

  fseek(instream,filelocation,SEEK_SET);

  fread_uint16(&(my_unit->unittype),1,instream);
  fread_uchar(&(my_unit->direction),1,instream);
  

  fread_int32(&(my_unit->natoms),1,instream);
  fread_int32(&(my_unit->nblocks),1,instream);
  fread_int32(&(my_unit->ncells),1,instream);
  fread_int32(&(my_unit->unitnumber),1,instream);
  fread_uchar(&(my_unit->ncellperatom),1,instream);

  my_unit->unit_block = Calloc(my_unit->nblocks,cdf_unit_block);

  for (i=0; i < my_unit->nblocks; i++){
    fread_int32(&(my_unit->unit_block[i].natoms),1,instream);
    fread_int32(&(my_unit->unit_block[i].ncells),1,instream);
    fread_uchar(&(my_unit->unit_block[i].ncellperatom),1,instream);
    fread_uchar(&(my_unit->unit_block[i].direction),1,instream);
    fread_int32(&(my_unit->unit_block[i].firstatom),1,instream);
    fread_int32(&(my_unit->unit_block[i].unused),1,instream);
    fread_char(my_unit->unit_block[i].blockname,64,instream); 

    my_unit->unit_block[i].unit_cells = Calloc(my_unit->unit_block[i].ncells,cdf_unit_cell);

    for (j=0; j < my_unit->unit_block[i].ncells; j++){
      fread_int32(&(my_unit->unit_block[i].unit_cells[j].atomnumber),1,instream);
      fread_uint16(&(my_unit->unit_block[i].unit_cells[j].x),1,instream);
      fread_uint16(&(my_unit->unit_block[i].unit_cells[j].y),1,instream);
      fread_int32(&(my_unit->unit_block[i].unit_cells[j].indexpos),1,instream);
      fread_char(&(my_unit->unit_block[i].unit_cells[j].pbase),1,instream);
      fread_char(&(my_unit->unit_block[i].unit_cells[j].tbase),1,instream);
    }


  }


  return 1;

}


static void dealloc_cdf_xda(cdf_xda *my_cdf){

  int i;

  for (i=0; i < my_cdf->header.n_units; i++){
    Free(my_cdf->probesetnames[i]);
  }              
  Free(my_cdf->probesetnames);

  Free(my_cdf->qc_start);
  Free(my_cdf->units_start);

  for (i=0; i < my_cdf->header.n_qc_units; i++){
    Free(my_cdf->qc_units[i].qc_probes);
  }

  Free(my_cdf->qc_units);


  for (i=0; i < my_cdf->header.n_units; i++){
    Free(my_cdf->units[i].unit_block);
  }
  Free(my_cdf->units);
  Free(my_cdf->header.ref_seq);

} 



/*************************************************************
 **
 ** int read_cdf_xda(char *filename)
 **
 ** filename - Name of the prospective binary cel file
 **
 ** Returns 1 if the file was completely successfully parsed
 ** otherwise 0 (and possible prints a message to screen)
 ** 
 **
 **
 **
 *************************************************************/

static int read_cdf_xda(char *filename,cdf_xda *my_cdf){

  FILE *infile;

  int i;

  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }

  if (!fread_int32(&my_cdf->header.magicnumber,1,infile)){
    return 0;
  }

  if (!fread_int32(&my_cdf->header.version_number,1,infile)){
    return 0;
  }


  if (my_cdf->header.magicnumber != 67){
    Rprintf("Magic number is not 67. This is probably not a binary cdf file.\n");
    return 0;
  }

  if (my_cdf->header.version_number != 1){
    Rprintf("Don't know if version %d binary cdf files can be handled.\n",my_cdf->header.version_number);
    return 0;
  } 
  if (!fread_uint16(&my_cdf->header.cols,1,infile)){
    return 0;
  }
  if (!fread_uint16(&my_cdf->header.rows,1,infile)){
    return 0;
  }
 
  if (!fread_int32(&my_cdf->header.n_units,1,infile)){
    return 0;
  }

  if (!fread_int32(&my_cdf->header.n_qc_units,1,infile)){
    return 0;
  }

  
  if (!fread_int32(&my_cdf->header.len_ref_seq,1,infile)){
    return 0;
  }
  
  my_cdf->header.ref_seq = Calloc(my_cdf->header.len_ref_seq,char);

  fread_char(my_cdf->header.ref_seq, my_cdf->header.len_ref_seq, infile);
  my_cdf->probesetnames = Calloc(my_cdf->header.n_units,char *);


  for (i =0; i < my_cdf->header.n_units;i++){
    my_cdf->probesetnames[i] = Calloc(64,char);
    if (!fread_char(my_cdf->probesetnames[i], 64, infile)){
      return 0;
    }
  }



  my_cdf->qc_start = Calloc(my_cdf->header.n_qc_units,int);
  my_cdf->units_start = Calloc(my_cdf->header.n_units,int);

  if (!fread_int32(my_cdf->qc_start,my_cdf->header.n_qc_units,infile) 
      || !fread_int32(my_cdf->units_start,my_cdf->header.n_units,infile)){
    return 0;
  }
  

  /* We will read in all the QC and Standard Units, rather than  
     random accessing what we need */
  my_cdf->qc_units = Calloc(my_cdf->header.n_qc_units,cdf_qc_unit);
  
  
  for (i =0; i < my_cdf->header.n_qc_units; i++){
    if (!read_cdf_qcunit(&my_cdf->qc_units[i],my_cdf->qc_start[i],infile)){
      return 0;
    }
  }
    
  my_cdf->units = Calloc(my_cdf->header.n_units,cdf_unit);


  for (i=0; i < my_cdf->header.n_units; i++){
    if (!read_cdf_unit(&my_cdf->units[i],my_cdf->units_start[i],infile)){
      return 0;
    }
  }
  

#ifdef READ_CDF_DEBUG
   Rprintf("%d %d %d %d  %d\n",my_cdf->header.cols,my_cdf->header.rows,my_cdf->header.n_units,my_cdf->header.n_qc_units,my_cdf->header.len_ref_seq);
  for (i =0; i < my_cdf->header.n_units;i++){
    Rprintf("%s\n",my_cdf->probesetnames[i]);
  }

  for (i =0; i < my_cdf->header.n_qc_units;i++){
    Rprintf("%d\n",my_cdf->qc_start[i]);
  }
  
  for (i =0; i < my_cdf->header.n_qc_units;i++){
    Rprintf("%d\n",my_cdf->units_start[i]);
  }

  Rprintf("%d %d\n",my_cdf->qc_units[0].type,my_cdf->qc_units[0].n_probes);
  
  for (i=0; i < my_cdf->qc_units[0].n_probes; i++){
    Rprintf("%d %d %d %u %d\n",my_cdf->qc_units[0].qc_probes[i].x,my_cdf->qc_units[0].qc_probes[i].y,
	    my_cdf->qc_units[0].qc_probes[i].probelength,
	    my_cdf->qc_units[0].qc_probes[i].pmflag,
	    my_cdf->qc_units[0].qc_probes[i].bgprobeflag);

  }
  

  Rprintf("%u %u %d %d %d %d %u\n",my_cdf->units[0].unittype,my_cdf->units[0].direction,
	  my_cdf->units[0].natoms,
	  my_cdf->units[0].nblocks,
	  my_cdf->units[0].ncells,
	  my_cdf->units[0].unitnumber,
	  my_cdf->units[0].ncellperatom);

    Rprintf("%d %d %u %u %d %d %s\n",my_cdf->units[0].unit_block[0].natoms,my_cdf->units[0].unit_block[0].ncells,
	    my_cdf->units[0].unit_block[0].ncellperatom,
	    my_cdf->units[0].unit_block[0].direction,
	    my_cdf->units[0].unit_block[0].firstatom,
	    my_cdf->units[0].unit_block[0].unused,
	    my_cdf->units[0].unit_block[0].blockname);
    
    for (i=0; i <my_cdf->units[0].unit_block[0].ncells  ; i++){
      Rprintf("%d %u %u %d %c %c\n",
	      my_cdf->units[0].unit_block[0].unit_cells[i].atomnumber,
	      my_cdf->units[0].unit_block[0].unit_cells[i].x,
	      my_cdf->units[0].unit_block[0].unit_cells[i].y,
	      my_cdf->units[0].unit_block[0].unit_cells[i].indexpos,
	      my_cdf->units[0].unit_block[0].unit_cells[i].pbase,
	      my_cdf->units[0].unit_block[0].unit_cells[i].tbase);
    }
#endif
    
  fclose(infile);
  return 1;

  // fseek()
}






static int check_cdf_xda(char *filename){

  FILE *infile;

  int i;
  int magicnumber,version_number;

  if ((infile = fopen(filename, "rb")) == NULL)
    {
      error("Unable to open the file %s",filename);
      return 0;
    }

  if (!fread_int32(&magicnumber,1,infile)){
    error("File corrupt or truncated?");
    return 0;
  }

  if (!fread_int32(&version_number,1,infile)){ 
    error("File corrupt or truncated?");
    return 0;
  }


  if (magicnumber != 67){
    // error("Magic number is not 67. This is probably not a binary cdf file.\n");
    return 0;
  }

  if (version_number != 1){
    // error("Don't know if version %d binary cdf files can be handled.\n",my_cdf->header.version_number);
    return 0;
  } 

  return 1;

}


static int isPM(char pbase,char tbase){
  /*
  if (Pbase.Cmp(Tbase) == 0){
    *isPM = false;
  } else if (((Pbase.Cmp("A")== 0) && (Tbase.Cmp("T") != 0)) || ((Pbase.Cmp("T")
== 0) && (Tbase.Cmp("A") != 0))){
    *isPM = false;
  } else if (((Pbase.Cmp("C")== 0) && (Tbase.Cmp("G") != 0)) || ((Pbase.Cmp("G")
== 0) && (Tbase.Cmp("C") != 0))){
    *isPM = false;
  } else {
    *isPM = true;
  }
  */

  pbase = toupper(pbase);
  tbase = toupper(tbase);

  if (pbase == tbase){
    return 0;
  } else if ((( pbase == 'A') && (tbase != 'T')) || (( pbase == 'T') && (tbase != 'A'))){
    return 0;
  } else if ((( pbase == 'C') && (tbase != 'G')) || (( pbase == 'G') && (tbase != 'C'))){
    return 0;
  }
  return 1;


}



SEXP CheckCDFXDA(SEXP filename){
  SEXP tmp;
  int good;
  char *cur_file_name;
  
  cur_file_name = CHAR(VECTOR_ELT(filename,0));
  
  good = check_cdf_xda(cur_file_name);
  
  PROTECT(tmp= allocVector(INTSXP,1));

  INTEGER(tmp)[0] = good;

  UNPROTECT(1);
  return tmp;
}






SEXP ReadCDFFile(SEXP filename){
  
  SEXP CDFInfo;
  SEXP Dimensions;
  SEXP LocMap,tempLocMap;
  SEXP CurLocs;
  SEXP PSnames,tempPSnames;
  SEXP ColNames;
  SEXP dimnames;

  cdf_xda my_cdf;
  char *cur_file_name;
  char *tmp_name;

  int i,j,k;
  int cur_blocks,cur_cells, cur_atoms;
  int which_probetype;
  int which_psname=0;

  cdf_unit_cell *current_cell;

  double *curlocs;
  
  int nrows, ncols;

 
  cur_file_name = CHAR(VECTOR_ELT(filename,0));

  if (!read_cdf_xda(cur_file_name,&my_cdf)){
    error("Problem reading binary cdf file %s. Possibly corrupted or truncated?\n",cur_file_name);
  }
  

  /* We output:
     nrows, ncols in an integer vector, plus a list of probesets PM MM locations (in the BioC style) */
  PROTECT(CDFInfo = allocVector(VECSXP,2));
  PROTECT(Dimensions = allocVector(REALSXP,2));

  if (my_cdf.units[0].unittype ==1){ 
    PROTECT(LocMap = allocVector(VECSXP,my_cdf.header.n_units));
    PROTECT(PSnames = allocVector(STRSXP,my_cdf.header.n_units));
  } else {
    PROTECT(tempLocMap = allocVector(VECSXP,2*my_cdf.header.n_units));
    PROTECT(tempPSnames = allocVector(STRSXP,2*my_cdf.header.n_units));
  }

  NUMERIC_POINTER(Dimensions)[0] = (double)my_cdf.header.rows;
  NUMERIC_POINTER(Dimensions)[1] = (double)my_cdf.header.cols;
  

  for (i=0; i < my_cdf.header.n_units; i++){
#ifdef READ_CDF_DEBUG
    printf("%d\n",i);
#endif
    cur_blocks = my_cdf.units[i].nblocks;

#ifdef READ_CDF_DEBUG
    Rprintf("New Block: ");
#endif
    if (my_cdf.units[i].unittype ==1){
      /* Expression analysis */
      for (j=0; j < cur_blocks; j++){
	
#ifdef READ_CDF_DEBUG
	Rprintf("%s ",my_cdf.units[i].unit_block[j].blockname);
#endif

	cur_cells = my_cdf.units[i].unit_block[j].ncells;
	cur_atoms = my_cdf.units[i].unit_block[j].natoms; 
	
	SET_VECTOR_ELT(PSnames,i,mkChar(my_cdf.units[i].unit_block[j].blockname));
	
	PROTECT(CurLocs = allocMatrix(REALSXP,cur_atoms,2));
	PROTECT(ColNames = allocVector(STRSXP,2));
	PROTECT(dimnames = allocVector(VECSXP,2));
	SET_VECTOR_ELT(ColNames,0,mkChar("pm"));
	SET_VECTOR_ELT(ColNames,1,mkChar("mm"));
	
	curlocs = NUMERIC_POINTER(AS_NUMERIC(CurLocs));
	
	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[j].unit_cells[k]);
	  
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber+ cur_atoms] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	}
	

      
	SET_VECTOR_ELT(dimnames,1,ColNames);
	setAttrib(CurLocs, R_DimNamesSymbol, dimnames);
	SET_VECTOR_ELT(LocMap,i,CurLocs);
	UNPROTECT(3);
      }
    } else if (my_cdf.units[i].unittype == 2){
      /* Genotyping array */

#ifndef READ_CDF_NOSNP
      if (cur_blocks == 1){
	
	cur_cells = my_cdf.units[i].unit_block[0].ncells;
	cur_atoms = my_cdf.units[i].unit_block[0].natoms; 
	
	SET_VECTOR_ELT(tempPSnames,which_psname,mkChar(my_cdf.units[i].unit_block[0].blockname));
	
	PROTECT(CurLocs = allocMatrix(REALSXP,cur_atoms,2));
	PROTECT(ColNames = allocVector(STRSXP,2));
	PROTECT(dimnames = allocVector(VECSXP,2));
	SET_VECTOR_ELT(ColNames,0,mkChar("pm"));
	SET_VECTOR_ELT(ColNames,1,mkChar("mm"));
	
	curlocs = NUMERIC_POINTER(AS_NUMERIC(CurLocs));
	
	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[0].unit_cells[k]);
	  
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber+ cur_atoms] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	}
	

      
	SET_VECTOR_ELT(dimnames,1,ColNames);
	setAttrib(CurLocs, R_DimNamesSymbol, dimnames);
	SET_VECTOR_ELT(tempLocMap,which_psname,CurLocs);
	UNPROTECT(3);
	which_psname++;

      } else if (cur_blocks == 4){
	for (j=0; j < cur_blocks; j++){
#ifdef READ_CDF_DEBUG_SNP
	  Rprintf("%s %s\n",my_cdf.probesetnames[i],my_cdf.units[i].unit_block[j].blockname);
#endif
	}
	
	j = 0;
	cur_cells = my_cdf.units[i].unit_block[0].ncells;
	cur_atoms = my_cdf.units[i].unit_block[0].natoms; 
	if (strlen(my_cdf.units[i].unit_block[j].blockname) == 1){
	  tmp_name = Calloc(strlen(my_cdf.probesetnames[i])+2,char);
	  tmp_name = strcpy(tmp_name,my_cdf.probesetnames[i]);
	  tmp_name = strcat(tmp_name,my_cdf.units[i].unit_block[j].blockname);
	  SET_VECTOR_ELT(tempPSnames,which_psname,mkChar(tmp_name));
	  Free(tmp_name);
	} else {
	  SET_VECTOR_ELT(tempPSnames,which_psname,mkChar(my_cdf.units[i].unit_block[0].blockname));
	}

	PROTECT(CurLocs = allocMatrix(REALSXP,2*cur_atoms,2));
	PROTECT(ColNames = allocVector(STRSXP,2));
	PROTECT(dimnames = allocVector(VECSXP,2));
	SET_VECTOR_ELT(ColNames,0,mkChar("pm"));
	SET_VECTOR_ELT(ColNames,1,mkChar("mm"));
	
	curlocs = NUMERIC_POINTER(AS_NUMERIC(CurLocs));


	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[0].unit_cells[k]);
	  //	  Rprintf("%d %d  %u %u \n",cur_cells, current_cell->atomnumber,current_cell->x,current_cell->y);
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber+ 2*cur_atoms] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	  if (current_cell->x + current_cell->y*(my_cdf.header.rows) + 1 == 370737){
	    Rprintf("%d %c %c",isPM(current_cell->pbase,current_cell->tbase),current_cell->pbase,current_cell->tbase);
	  }
	}

	j=2;
	cur_cells = my_cdf.units[i].unit_block[2].ncells;
	cur_atoms = my_cdf.units[i].unit_block[2].natoms; 
	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[2].unit_cells[k]);
	  //Rprintf("half : %d %d  %u %u \n",cur_cells, current_cell->atomnumber,current_cell->x,current_cell->y);
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber - (cur_atoms)] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber - (cur_atoms)+ 2*cur_atoms] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	}
	
	SET_VECTOR_ELT(dimnames,1,ColNames);
	setAttrib(CurLocs, R_DimNamesSymbol, dimnames);
	SET_VECTOR_ELT(tempLocMap,which_psname,CurLocs);	
	UNPROTECT(3);
	which_psname++;





	j = 1;	
	cur_cells = my_cdf.units[i].unit_block[1].ncells;
	cur_atoms = my_cdf.units[i].unit_block[1].natoms; 
	if (strlen(my_cdf.units[i].unit_block[j].blockname) == 1){
	  tmp_name = Calloc(strlen(my_cdf.probesetnames[i])+2,char);
	  tmp_name = strcpy(tmp_name,my_cdf.probesetnames[i]);
	  tmp_name = strcat(tmp_name,my_cdf.units[i].unit_block[j].blockname);
	  SET_VECTOR_ELT(tempPSnames,which_psname,mkChar(tmp_name));
	  Free(tmp_name);
	} else {
	  SET_VECTOR_ELT(tempPSnames,which_psname,mkChar(my_cdf.units[i].unit_block[1].blockname));
	}
	PROTECT(CurLocs = allocMatrix(REALSXP,2*cur_atoms,2));
	PROTECT(ColNames = allocVector(STRSXP,2));
	PROTECT(dimnames = allocVector(VECSXP,2));
	SET_VECTOR_ELT(ColNames,0,mkChar("pm"));
	SET_VECTOR_ELT(ColNames,1,mkChar("mm"));
	curlocs = NUMERIC_POINTER(AS_NUMERIC(CurLocs));
	
	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[1].unit_cells[k]);
	  //Rprintf("Dual : %d %d  %u %u \n",cur_cells, current_cell->atomnumber,current_cell->x,current_cell->y);
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber - (cur_atoms)] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber - (cur_atoms)+ 2*cur_atoms] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	}
		
	j=3;
	cur_cells = my_cdf.units[i].unit_block[3].ncells;
	cur_atoms = my_cdf.units[i].unit_block[3].natoms; 
	for (k=0; k < cur_cells; k++){
	  current_cell = &(my_cdf.units[i].unit_block[3].unit_cells[k]);
	  //Rprintf("half deux : %d %d  %d %u %u \n",cur_cells, current_cell->atomnumber, cur_atoms,current_cell->x,current_cell->y);
	  if(isPM(current_cell->pbase,current_cell->tbase)){
	    curlocs[current_cell->atomnumber - (2*cur_atoms)] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1;           //"y*", sizex, "+x+1";
	  } else {
	    curlocs[current_cell->atomnumber] = current_cell->x + current_cell->y*(my_cdf.header.rows) + 1; 
	  }
	}
	
	SET_VECTOR_ELT(dimnames,1,ColNames);
	setAttrib(CurLocs, R_DimNamesSymbol, dimnames);
	SET_VECTOR_ELT(tempLocMap,which_psname,CurLocs);	
	UNPROTECT(3);
	which_psname++;















      } else {
	error("makecdfenv does not currently know how to handle cdf files of this type (genotyping with blocks != 1 or 4.)"); 	
      }
#else
      error("makecdfenv does not currently know how to handle cdf files of this type (genotyping).");
#endif




    } else { 
      error("makecdfenv does not currently know how to handle cdf files of this type (ie not expression or genotyping)"); 
    }

    
#ifdef READ_CDF_DEBUG
    Rprintf("\n");    
#endif
  }

  if (my_cdf.units[0].unittype ==2){
    PROTECT(PSnames = allocVector(STRSXP,which_psname));
    PROTECT(LocMap = allocVector(VECSXP,which_psname));
    for (i =0; i < which_psname; i++){
      SET_VECTOR_ELT(PSnames,i,mkChar(CHAR(VECTOR_ELT(tempPSnames,i))));
      SET_VECTOR_ELT(LocMap,i,VECTOR_ELT(tempLocMap,i));
    }
    
  }
#ifdef READ_CDF_DEBUG
  Rprintf("%d \n",which_psname);
#endif
  setAttrib(LocMap,R_NamesSymbol,PSnames);
  SET_VECTOR_ELT(CDFInfo,0,Dimensions);
  SET_VECTOR_ELT(CDFInfo,1,LocMap);
  if (my_cdf.units[0].unittype ==2){
    UNPROTECT(6);
  } else {
    UNPROTECT(4);
  }

  dealloc_cdf_xda(&my_cdf);
  return CDFInfo;

}








