/* cddio.c:  Basic Input and Output Procedures for cdd.c
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.55a, December 18, 1994 
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. Jan.23 ,1994 or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


void SetInputFile(FILE **f, boolean *success)
{
  boolean opened=FALSE,stop, quit=FALSE;
  int i,dotpos=0, semipos=0;
  char ch;
  char *tempname;
  
  *success=FALSE;
  while (!opened && !quit) {
    if (FileInputMode!=Auto){
      printf("\n>> Input file (*.ine) : ");
      scanf("%s",inputfile);
      ch=getchar();
    }
    stop=FALSE;
    for (i=0; i<filenamelen && !stop; i++){
      ch=inputfile[i];
      switch (ch) {
        case '.': 
          dotpos=i+1;
          break;
        case ';':  case ' ':  case '\0':  case '\n':  case '\t':     
          if (ch==';'){
            semipos=i+1;
            FileInputMode=SemiAuto;   
            /* semicolon at the end of the filename
            -> output file names will be creeated with defaults. */
          }
          stop=TRUE;
          tempname=(char*)calloc(filenamelen,sizeof ch);
          strncpy(tempname, inputfile, i);
          strcpy(inputfile,tempname);
          break;
      }
    }
    if (dotpos>0){
      strncpy(ifilehead, inputfile, dotpos-1);
    }else{
      strcpy(ifilehead, inputfile);
    }
    if (debug){
      printf("inputfile name: %s\n", inputfile);  
      printf("inputfile name head: %s\n", ifilehead);  
      printf("semicolon pos: %d\n", semipos);
    }  
    if ( ( *f = fopen(inputfile, "r") )!= NULL) {
      if (DynamicWriteOn) printf("input file %s is open\n", inputfile);
      opened=TRUE;
      *success=TRUE;
    }
    else{
      printf("The file %s not found\n",inputfile);
      if (FileInputMode==Auto) {
        quit=TRUE;
      }
    }
  }
}

void SetWriteFile(FILE **f, DataFileType fname, char cflag, char *fscript)
{
  boolean quit=FALSE;
  char *extension;
  
  switch (cflag) {
    case 'o':
      extension=".ext";break;
    case 'a':
      extension=".adj";break;
    case 'i':
      extension=".icd";break;
    case 'l':
      extension=".ddl";break;
    case 'd':
      extension=".dex";break;   /* decomposition output */
    case 'p':
      extension="sub.ine";break;
    case 'v':
      extension=".solved";break;
    default:
      extension=".xxx";break;
  }
  if (FileInputMode==Manual){
    while (!quit) {
      printf("\n>> %s file name (*%s)   : ",fscript, extension);
      scanf("%s",fname);
      if (fname[0]==';'|| fname[0]<'0'){
        quit=TRUE;  /* default file name */
      } 
      else if (strcmp(inputfile, fname)!=0){
        *f = fopen(fname, "w");
        if (DynamicWriteOn) printf("write file %s is open\n",fname);
        goto _L99;
      }
      else {
        printf("%s file %s must have a name different from inputfile.\n",fscript,fname);
      }
    }
  }
  /* Auto or SemiAuto FileInput */
  strcpy(fname,ifilehead); 
  strcat(fname,extension); 
  if (strcmp(inputfile, fname)==0) {
    strcpy(fname,inputfile); 
    strcat(fname,extension); 
  }
  *f = fopen(fname, "w");
  if (DynamicWriteOn) printf("%s file %s is open\n",fscript,fname);
_L99:;
}


void SetNumberType(char *line)
{
  if (strncmp(line, "integer", 7)==0) {
    Number = Integer;
    return;
  }
  else if (strncmp(line, "rational", 8)==0) {
    Number = Rational;
    Error=ImproperInputFormat;  /* Rational Input not supported */
    return;
  }
  else if (strncmp(line, "real", 4)==0) {
    Number = Integer;
    return;
  }
  else { 
    Number=Unknown;
    Error=ImproperInputFormat;
  }
}

void ProcessCommandLine(char *line)
{
  colrange j;
  long var,msize;
  double cost;

  if (strncmp(line, "dynout_off", 10)==0) {
    DynamicRayWriteOn = FALSE;
    return;
  }
  if (strncmp(line, "stdout_off", 10)==0) {
    DynamicRayWriteOn = FALSE;
    DynamicWriteOn = FALSE;
    return;
  }
  if (strncmp(line, "logfile_on", 10)==0) {
    LogWriteOn = TRUE;
    return;
  }
  if (strncmp(line, "logfile_off", 11)==0) {
    LogWriteOn = FALSE;
    return;
  }
  if (strncmp(line, "hull", 4)==0) {
    Conversion = ExtToIne;
    return;
  }
  if (strncmp(line, "incidence", 9)==0) {
    IncidenceOutput = IncSet;
    return;
  }
  if (strncmp(line, "#incidence", 10)==0) {
    IncidenceOutput = IncCardinality;
    return;
  }
  if (strncmp(line, "adjacency", 9)==0) {
    AdjacencyOutput = OutputAdjacency;
    return;
  }
/*  algebraic option is not efficient in most cases and deleted from Version 051 */
/*
  if (strncmp(line, "algebraic", 9)==0) {
    AdjacencyTest = Algebraic;
    return;
  }
*/
  if (strncmp(line, "nondegenerate", 14)==0) {
    NondegAssumed = TRUE;
    return;
  }
  if (strncmp(line, "minindex", 8)==0) {
    HyperplaneOrder = MinIndex;
    return;
  }
  if (strncmp(line, "maxindex", 8)==0) {
    HyperplaneOrder = MaxIndex;
    return;
  }
  if (strncmp(line, "mincutoff", 9)==0) {
    HyperplaneOrder = MinCutoff;
    return;
  }
  if (strncmp(line, "maxcutoff", 9)==0) {
    HyperplaneOrder = MaxCutoff;
    return;
  }
  if (strncmp(line, "mixcutoff", 9)==0) {
    HyperplaneOrder = MixCutoff;
    return;
  }
  if (strncmp(line, "lexmin", 6)==0) {
    HyperplaneOrder = LexMin;
    return;
  }
  if (strncmp(line, "lexmax", 6)==0) {
    HyperplaneOrder = LexMax;
    return;
  }
  if (strncmp(line, "random", 6)==0) {
    HyperplaneOrder = RandomRow;
    fscanf(reading,"%d", &rseed);
    if (rseed <= 0) rseed = 1;
    return;
  }
  if (strncmp(line, "lineshell", 9)==0) {
    HyperplaneOrder = LineShelling;
    return;
  }
  if (strncmp(line, "initbasis_at_bottom", 19)==0) {
    InitBasisAtBottom = TRUE;
    return;
  }
  if (strncmp(line, "verify_input", 12)==0) {
    VerifyInput = TRUE;
    return;
  }
  if (strncmp(line, "debug", 5)==0) {
    debug = TRUE;
    return;
  }
  if ((strncmp(line, "partial_enum", 12)==0 || strncmp(line, "equality", 8)==0) 
    && RestrictedEnumeration==FALSE) {
    fscanf(reading,"%ld", &msize);
    for (j=1;j<=msize;j++) {
      fscanf(reading,"%ld",&var);
      EqualityIndex[var]=1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Partial Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    RestrictedEnumeration=TRUE;
    return;
  }
  if (strncmp(line, "strict_ineq", 11)==0 && RelaxedEnumeration==FALSE) {
    fscanf(reading,"%ld", &msize);
    for (j=1;j<=msize;j++) {
      fscanf(reading,"%ld",&var);
      EqualityIndex[var]=-1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Partial Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    RelaxedEnumeration=TRUE;
    return;
  }
  if (strncmp(line, "preprojection", 13)==0 && Conversion != Projection) {
    set_initialize(&projvars,nn);
    fscanf(reading,"%ld", &projdim);
    if (debug) printf("dimension of projection = %ld  in variables:\n",projdim);
    for (j=1;j<=projdim;j++) {
      fscanf(reading,"%ld",&var);
      if (debug) printf(" %ld",var);
      if (Inequality==NonzeroRHS)
        set_addelem(projvars,var+1);
      else
        set_addelem(projvars,var);
    }
    Conversion=Projection;
    return;
  }
  if (strncmp(line, "maximize", 8)==0 && Conversion != LPmax) {
    if (debug) printf("linear maximization is chosen.\n");
    for (j=0;j<nn;j++) {
      fscanf(reading,"%lf",&cost);
      LPcost[j]=cost;
      if (debug) printf(" cost[%ld] = %.9E\n",j,LPcost[j]);
    }
    Conversion=LPmax;
    if (debug) {
      printf("\n");
    }
    return;
  }
  if (strncmp(line, "minimize", 8)==0 && Conversion != LPmin) {
    if (debug) printf("linear minimization is chosen.\n");
    for (j=0;j<nn;j++) {
      fscanf(reading,"%lf",&cost);
      LPcost[j]=cost;
      if (debug) printf(" cost[%ld] = %.9E\n",j,LPcost[j]);
    }
    Conversion=LPmin;
    if (debug) {
      printf("\n");
    }
    return;
  }
  if (strncmp(line, "find_interior", 13)==0 && Conversion != InteriorFind) {
    printf("Interior finding option is chosen.\n");
    Conversion=InteriorFind;
    return;
  }
  if (strncmp(line, "row_decomp", 10)==0 && !RowDecomposition) {
    printf("Row decomposition is chosen.\n");
    RowDecomposition=TRUE;
    return;
  }
}

void AmatrixInput(boolean *successful)
{
  long i,j;
  double value;
  long value1,value2;
  boolean found=FALSE,decided=FALSE, fileopened;
  char command[wordlenmax], numbtype[wordlenmax], stemp[wordlenmax], line[linelenmax];

  *successful = FALSE;

  SetInputFile(&reading, &fileopened);
  if (!fileopened){
    Error=FileNotFound;
    goto _L99;
  }

  while (!found)
  {
    if (fscanf(reading,"%s",command)==EOF) {
      Error=ImproperInputFormat;
      goto _L99;
    }
    else if (strncmp(command, "begin", 5)==0) {
      found=TRUE;
    }
  }
  fscanf(reading, "%ld %ld %s", &minput, &ninput, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", minput, ninput, numbtype);
  SetNumberType(numbtype);
  if (Number==Unknown || Number == Rational) {
      goto _L99;
    } 
  Inequality=ZeroRHS;
  for (i=1; i<=minput && !decided; i++) {
    fscanf(reading,"%lf", &value);
    if (fabs(value) > zero) {
      Inequality = NonzeroRHS;
      decided=TRUE;
    }
    for (j=2; j<=ninput; j++) {
      fscanf(reading,"%lf", &value);
    }
    fgets(line,linelenmax,reading);
    if (debug) printf("comments to be skipped: %s\n",line);
    if (debug) putchar('\n');
  }
  if (Inequality==NonzeroRHS) {
  	printf("Nonhomogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
  	nn = ninput;
  }
  else {
    printf("Homogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
  	nn = ninput-1;
  }
  if (nn > NMAX || minput > MMAX) {
    Error = DimensionTooLarge;
    goto _L99;
  }
  
  RestrictedEnumeration = FALSE;
  RelaxedEnumeration = FALSE;
  EqualityIndex=(long *)calloc(minput+2, sizeof *EqualityIndex);
  for (i = 0; i <= minput+1; i++) EqualityIndex[i]=0;

  while (!feof(reading)) {
    fscanf(reading,"%s", command);
    ProcessCommandLine(command);
  } 
  fclose(reading);

  reading=fopen(inputfile, "r");
  found=FALSE;
  while (!found)
  {
  	if (!feof(reading)) {
  	  fscanf(reading,"%s",command);
  	  if (strncmp(command, "begin", 5)==0) {
  		  found=TRUE;
  	  }
  	}
  	else {
  	  Error=ImproperInputFormat;
  	  goto _L99;
  	}
  }
  fscanf(reading, "%ld %ld %s", &value1, &value2, stemp);
  for (i = 1; i <= minput; i++) {
    AA[i-1]=(double *) calloc(ninput, sizeof value);
    for (j = 1; j <= ninput; j++) {
      fscanf(reading, "%lf", &value);
      if (Inequality==NonzeroRHS) 
      	AA[i-1][j - 1] = value;
      else if (j>=2) {
        AA[i-1][j - 2] = value;
	  }
	  if (debug) printf("a(%3ld,%5ld) = %10.4f\n",i,j,value);
    }  /*of j*/
    fgets(line,linelenmax,reading);
    if (debug) printf("comments to be skipped: %s\n",line);
    if (debug) putchar('\n');
  }  /*of i*/
  if (fscanf(reading,"%s",command)==EOF) {
   	 Error=ImproperInputFormat;
  	 goto _L99;
  }
  else if (strncmp(command, "end", 3)!=0) {
     if (debug) printf("'end' missing or illegal extra data: %s\n",command);
   	 Error=ImproperInputFormat;
  	 goto _L99;
  }
  
  switch (Conversion) {
  case IneToExt:
    if (Inequality==NonzeroRHS){
      mm = minput + 1;
      AA[mm-1]=(double *) calloc(ninput, sizeof value);
      for (j = 1; j <= ninput; j++) {   /*artificial row for x_1 >= 0*/
        if (j == 1)
          AA[mm - 1][j - 1] = 1.0;
        else
          AA[mm - 1][j - 1] = 0.0;
      }
    } else {
      mm = minput;
    }
    break;

  case ExtToIne:
    mm = minput;
    break;

  case LPmax:  case LPmin:
    mm = minput + 1;
    OBJrow=mm;
    RHScol=1L;
    AA[mm-1]=(double *) calloc(ninput, sizeof value);
    for (j = 1; j <= ninput; j++) {   /*objective row */
 	   AA[mm - 1][j - 1] = LPcost[j-1];
 	}
	break;

  default:
    mm = minput;
  }
  SetInequalitySets(EqualityIndex);
  *successful = TRUE;
_L99: ;
  if (reading!=NULL) fclose(reading);
}

void WriteRayRecord(FILE *f, RayRecord *RR)
{
  long j;
  double scaler;

  if (Inequality == ZeroRHS) {
    fprintf(f, " %2d", 0);
    for (j = 0; j < nn; j++)
      WriteReal(f, RR->Ray[j]);
  } 
  else {
    scaler = fabs(RR->Ray[0]);
    if (scaler > zero) {
      if (RR->Ray[0] > 0) {
        if (Conversion == IneToExt) {
	      fprintf(f, " %2d", 1);
	      for (j = 1; j < nn; j++)
	      WriteReal(f, RR->Ray[j] / scaler);
        } 
        else {
	      /* hull computation is chosen */
          for (j = 0; j < nn; j++)
	        WriteReal(f, RR->Ray[j]);
	    }
      }
      else {
        /* hull computation must have been chosen, since RR->Ray[0] is negative */
	    for (j = 0; j < nn; j++)
	      WriteReal(f, RR->Ray[j]);
      }
    } 
    else {
      fprintf(f, " %2d", 0);
      for (j = 1; j < nn; j++)
        WriteReal(f, RR->Ray[j]);
    }
  }
  if (IncidenceOutput==IncCardinality) {
    fprintf(f," %4ld", set_card(RR->ZeroSet));
  }
  putc('\n', f);
}


void WriteRayRecord2(FILE *f, RayRecord *RR)
{
  long j;

  fprintf(f, " Ray = ");
  for (j = 0; j < nn; j++)
    fprintf(f, "%6.2f", RR->Ray[j]);
  putchar('\n');
  fprintf(f, " ZeroSet =");
  WriteSetElements(f, RR->ZeroSet);
  putc('\n', f);
}

void WriteSubMatrixOfAA(FILE *f,long *ChosenRow, long *ChosenCol,
      InequalityType ineq)
{
  long i,j;

  fprintf(f, "begin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",mm, set_card(ChosenCol)+1);
  else
    fprintf(f, "  %ld   %ld    real\n",mm, set_card(ChosenCol));
  for (i=1; i <= mm; i++) {
    if (set_member(i,ChosenRow)) {
      if (ineq==ZeroRHS){  /* If RHS==0, output zero first */
        WriteReal(f, 0);
        for (j=1; j <= nn; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteReal(f, AA[i-1][j-1]);
          }
        }
      }
      else {
        for (j=1; j <= nn; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteReal(f, AA[i-1][j-1]);
          }
        }
      }
      fprintf(f,"\n");
    }
  }
  fprintf(f,"end\n"); 
}

void WriteAmatrix(FILE *f, Amatrix A, long rowmax, long colmax,
      InequalityType ineq)
{
  long i,j;

  fprintf(f, "begin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",rowmax, colmax+1);
  else
    fprintf(f, "  %ld   %ld    real\n",rowmax, colmax);
  for (i=1; i <= rowmax; i++) {
    if (ineq==ZeroRHS)
      WriteReal(f, 0);  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j <= colmax; j++) {
      WriteReal(f, A[i-1][j-1]);
    }
    fprintf(f,"\n");
  }
  fprintf(f, "end\n");
}

void WriteBmatrix(FILE *f, Bmatrix T)
{
  colrange j1, j2;

  for (j1 = 0; j1 < nn; j1++) {
    for (j2 = 0; j2 < nn; j2++) {
      fprintf(f, "%15.7f ", T[j1][j2]);
    }  /*of j2*/
    putc('\n', f);
  }  /*of j1*/
  putc('\n', f);
}

void WriteReal(FILE *f, double x)
{
  long ix1,ix2,ix;

  ix1= fabs(x) * 10000. + 0.5;
  ix2= (fabs(x) + 0.5);
  ix2= ix2*10000;
  if ( ix1 == ix2) {
    if (x>0) {
      ix = x + 0.5;
    } else {
      ix = -x + 0.5;
      ix = -ix;
    }
    fprintf(f, " %2ld", ix);
  } else
    fprintf(f, " % .9E", x);
}


void WriteSetElements(FILE *f, long *S)
{
  rowrange i;

  for (i = 1; i <= mm; i++) {
    if (set_member(i, S))
      fprintf(f, " %4ld", i);
  }
}

void WriteIncidence(FILE *f, RayRecord *RR)
{
  rowset cset;
  long zcar;

  set_initialize(&cset,mm);
  zcar = set_card(RR->ZeroSet);
  switch (IncidenceOutput) {

  case IncCardinality:
    fprintf(f, "%8ld", zcar);
    break;

  case IncSet:
    if (mm - zcar >= zcar) {
      fprintf(f, "%8ld", zcar);
      WriteSetElements(f, RR->ZeroSet);
    } else {
      set_diff(cset, GroundSet, RR->ZeroSet);
      fprintf(f, "%8ld", zcar - mm);
      WriteSetElements(f, cset);
    }
    break;

  case IncOff:
    break;
  }
  putc('\n', f);
  set_free(cset);
}

void WriteHeading(void)
{
  WriteProgramDescription(stdout);
  printf("----------------------------------------------\n");
  printf(" Enumeration of all vertices and extreme rays\n");
  printf(" of a convex polyhedron P={ x : b - A x >= 0}\n");
  printf(" Use hull option for convex hull computation!\n");
  printf("----------------------------------------------\n");
}

void WriteProgramDescription(FILE *f)
{
  fprintf(f, "* cdd: Double Description Method C-Code:%s\n", DDVERSION);
  fprintf(f,"* %s\n",COPYRIGHT);
}

void WriteRunningMode(FILE *f)
{
  if (Conversion==IneToExt || Conversion==ExtToIne){ 
    switch (HyperplaneOrder) {

    case MinIndex:
      fprintf(f, "*HyperplaneOrder: MinIndex\n");
      break;

    case MaxIndex:
      fprintf(f, "*HyperplaneOrder: MaxIndex\n");
      break;

    case MinCutoff:
      fprintf(f, "*HyperplaneOrder: MinCutoff\n");
      break;

    case MaxCutoff:
      fprintf(f, "*HyperplaneOrder: MaxCutoff\n");
    break;

    case MixCutoff:
      fprintf(f, "*HyperplaneOrder: MixCutoff\n");
      break;

    case LexMin:
      fprintf(f, "*HyperplaneOrder: LexMin\n");
      break;

    case LexMax:
      fprintf(f, "*HyperplaneOrder: LexMax\n");
      break;

    case RandomRow:
      fprintf(f, "*HyperplaneOrder: Random,  Seed = %d\n",rseed);
      break;

    case LineShelling:
      fprintf(f, "*HyperplaneOrder: LineShelling\n");
      break;
    }
    if (NondegAssumed) {
      fprintf(f, "*Degeneracy preknowledge for computation: NondegenerateAssumed\n");
    }
    else {
      fprintf(f, "*Degeneracy preknowledge for computation: None (possible degeneracy)\n");
    }
  }
  switch (Conversion) {
    case ExtToIne:
      fprintf(f, "*Hull computation is chosen.\n");
      break;
    
    case IneToExt:
      fprintf(f, "*Vertex/Ray enumeration is chosen.\n");
      break;
    
    case LPmax:
      fprintf(f, "*LP (maximization) is chosen.\n");
      break;

    case LPmin:
      fprintf(f, "*LP (minimization) is chosen.\n");
      break;

    case Projection:
      fprintf(f, "*Preprojection is chosen.\n");
      break;
  
    default: break;
  }
  if (RestrictedEnumeration) {
    fprintf(f, "*The equality option is chosen.\n* => Permanently active rows are:");
    WriteSetElements(f,EqualitySet);
    fprintf(f,"\n");
  }
  if (RelaxedEnumeration) {
    fprintf(f, "*The strict_inequality option is chosen.\n* => Permanently nonactive rows are:");
    WriteSetElements(f,NonequalitySet);
    fprintf(f,"\n");
  }
}

void WriteRunningMode2(FILE *f)
{
  long j;

  if (Conversion==IneToExt || Conversion==ExtToIne){ 
    switch (HyperplaneOrder) {

    case MinIndex:
      fprintf(f, "minindex\n");
      break;

    case MaxIndex:
      fprintf(f, "maxindex\n");
      break;

    case MinCutoff:
      fprintf(f, "mincutoff\n");
      break;

    case MaxCutoff:
      fprintf(f, "maxcutoff\n");
    break;

    case MixCutoff:
      fprintf(f, "mixcutoff\n");
      break;

    case LexMin:
      fprintf(f, "lexmin\n");
      break;

    case LexMax:
      fprintf(f, "lexmax\n");
      break;

    case RandomRow:
      fprintf(f, "random  %d\n",rseed);
      break;

    case LineShelling:
      fprintf(f, "lineshelling\n");
      break;
    }
    if (NondegAssumed) {
      fprintf(f, "nondegenerate\n");
    }
  }
  switch (Conversion) {
    case ExtToIne:
      fprintf(f, "hull\n");
      break;
    
    case LPmax:
      fprintf(f, "maximize\n");
      for (j=0; j<nn; j++) WriteReal(f,LPcost[j]);
      fprintf(f,"\n");
      break;

    case LPmin:
      fprintf(f, "minimize\n");
      for (j=0; j<nn; j++) WriteReal(f,LPcost[j]);
      fprintf(f,"\n");
      break;

    case Projection:
      fprintf(f, "*Preprojection is chosen.\n");
      break;
  
    default: break;
  }
  if (RestrictedEnumeration) {
    fprintf(f, "equality ");
    WriteSetElements(f,EqualitySet);
    fprintf(f,"\n");
  }
  if (RelaxedEnumeration) {
    fprintf(f, "strict_inequality ");
    WriteSetElements(f,NonequalitySet);
    fprintf(f,"\n");
  }
}

void WriteCompletionStatus(FILE *f)
{
  if (Iteration<mm && CompStatus==AllFound) {
    fprintf(f,"*Computation completed at Iteration %4ld.\n", Iteration);
  } 
  if (CompStatus == RegionEmpty) {
    fprintf(f,"*Computation completed at Iteration %4ld because the region found empty.\n", Iteration);
  }   
}


void WriteTimes(FILE *f)
{ 
  long ptime,ptime_sec,ptime_minu, ptime_hour;
  
  /* ptime=difftime(endtime,starttime); */   /* This function is ANSI standard, but not available sometime */
  ptime=endtime-starttime;      /* This is to replace the line above, but it may not give correct time in seconds */ 
  ptime_hour=ptime/3600;
  ptime_minu=(ptime-ptime_hour*3600)/60;
  ptime_sec=ptime%60;
  fprintf(f, "*Computation starts     at %s", asctime(localtime(&starttime)));
  fprintf(f, "*            terminates at %s", asctime(localtime(&endtime)));
  fprintf(f, "*Total processor time = %ld seconds\n", ptime);
  fprintf(f, "*                     = %ld hour %ld min %ld sec\n", ptime_hour,ptime_minu,ptime_sec);
}


void WriteAdjacency(FILE *f)
{
  RayRecord *RayPtr1, *RayPtr2;
  long pos1, pos2, degree;
  boolean adj;
  node *headnode, *tailnode, *newnode, *prevnode;

  headnode=NULL; tailnode=NULL;
  switch (Conversion) {
  case IneToExt:
    if (AdjacencyOutput==OutputAdjacency)
      fprintf(f,"*Adjacency List of output (=vertices/rays)\n");
    break;
  case ExtToIne:
    if (AdjacencyOutput==OutputAdjacency)
      fprintf(f,"*Adjacency List of output (=inequalities=facets)\n");
      break;

  default:
    break;
  }
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, minput, ninput);
  fprintf(f, "*cdd output file: %s\n", outputfile);
  fprintf(f,"begin\n");
  fprintf(f,"  %ld\n",RayCount);
  if (RayCount==0){
    goto _L99;
  }
  LastRay->Next=NULL;
  for (RayPtr1=FirstRay, pos1=1;RayPtr1 != NULL; RayPtr1 = RayPtr1->Next, pos1++){
    for (RayPtr2=FirstRay, pos2=1,degree=0; RayPtr2 != NULL; RayPtr2 = RayPtr2->Next, pos2++){
      if (RayPtr1!=RayPtr2){
        CheckAdjacency2(&RayPtr1, &RayPtr2, &adj);
        if (adj) {
          degree++;
          if (degree==1){
            newnode=(node *)malloc(sizeof *newnode);
            newnode->key=pos2;
            newnode->next=NULL;
            headnode=newnode;
            tailnode=newnode;
          }
          else{
            newnode=(node *)malloc(sizeof *newnode);
            newnode->key=pos2;
            newnode->next=NULL;
            tailnode->next=newnode;
            tailnode=newnode;
          }
        }
      }
    }
    fprintf(f," %ld  %ld :", pos1, degree);
    for (newnode=headnode; newnode!=NULL; newnode=newnode->next, free(prevnode)){
      prevnode=newnode;
      fprintf(f," %ld", newnode->key);
    }
    fprintf(f,"\n");
  }
_L99:;
  fprintf(f,"end\n");
}

void WriteDDResult(void)
{
  RayRecord *TempPtr;

  if (!debug) writing=freopen(outputfile,"w",writing);
  time(&endtime);
  WriteProgramDescription(writing);
  fprintf(writing, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, minput, ninput);
  WriteRunningMode(writing);
  WriteCompletionStatus(writing);
  WriteTimes(writing);
  if (Conversion == ExtToIne)
    fprintf(writing,
      "*Since hull computation is chosen, the output is a minimal inequality system\n");
  fprintf(writing, "*FINAL RESULT:\n");
  if (DynamicWriteOn)
    printf("*Computation complete.\n");
  if (Conversion == IneToExt) {
	if (Inequality==ZeroRHS && Conversion == IneToExt){
      fprintf(writing, "*Number of Rays =%8ld\n", FeasibleRayCount);
	  fprintf(writing, "*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n");
      if (DynamicWriteOn){
        printf("*Number of Rays =%8ld\n", FeasibleRayCount);
	    printf("*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n");
	  }
    }else{
      if (DynamicWriteOn)
        printf("*Number of Vertices =%8ld,   Rays =%8ld\n",
	       VertexCount, FeasibleRayCount - VertexCount);
      fprintf(writing, "*Number of Vertices =%8ld,   Rays =%8ld\n",
	      VertexCount, FeasibleRayCount - VertexCount);
	}
  } else {
    if (DynamicWriteOn)
      printf("*Number of Facets =%8ld\n", FeasibleRayCount);
    fprintf(writing, "*Number of Facets =%8ld\n", FeasibleRayCount);
  }
  fprintf(writing, "begin\n");
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, nn + 1);
    break;
  case NonzeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, nn);
    break;
  }
  if (IncidenceOutput == IncSet) {
    writing_icd=freopen(icdfile,"w",writing_icd);
    switch (Conversion) {
    case IneToExt:
      fprintf(writing_icd,
	    "*Incidences of output(=vertices/rays) and input (=hyperplanes)\n");
      fprintf(writing_icd,
        "*   for each output, #incidence and the set of hyperplanes containing it\n");
      fprintf(writing_icd,
	    "*   or its complement with its cardinality with minus sign\n");
      break;
    case ExtToIne:
      fprintf(writing_icd,
	    "*Incidences of output(=facets) and input (=points)\n");
      fprintf(writing_icd,
        "*   for each output, #incidence and the set of points lying on it\n");
      fprintf(writing_icd,
	    "*   or its complement with its cardinality with minus sign\n");
      break;

    default:
      break;
    }
    fprintf(writing_icd, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, minput, ninput);
    fprintf(writing_icd, "*cdd output file: %s\n", outputfile);
    fprintf(writing_icd, "begin\n");
    fprintf(writing_icd, "%8ld%5ld%5ld\n", FeasibleRayCount, minput, mm);
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(writing, TempPtr);
      if (IncidenceOutput == IncSet)
        WriteIncidence(writing_icd, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  fprintf(writing, "end\n");
  if (Conversion == IneToExt) fprintf(writing, "hull\n");
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout);
    WriteTimes(stdout);
  }
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteRunningMode(writing_log);
    WriteCompletionStatus(writing_log);
    if (PreOrderedRun) {
      fprintf(writing_log, "*set_intersection total#, effective#, loss# = %ld   %ld   %ld\n", 
      count_int, count_int_good, count_int_bad);
    }
    WriteTimes(writing_log);
  }
  if (IncidenceOutput == IncSet)
    fprintf(writing_icd, "end\n");
  if (AdjacencyOutput != AdjOff){
    if (DynamicWriteOn) printf("Writing the adjacency file...\n");
    WriteAdjacency(writing_adj);
  }
}

void WriteDecompResult(void)
{
  RayRecord *TempPtr;
  static int callnumber=0;

  callnumber++;
  if (callnumber==1) {
    WriteProgramDescription(writing_dex);
    fprintf(writing_dex, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, minput, ninput);
    if (Conversion == ExtToIne)
      fprintf(writing_dex,
        "*Since hull computation is chosen, the output is a minimal inequality system\n");
  }
  fprintf(writing_dex, "begin\n");
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing_dex, " %8ld  %5ld    real\n", FeasibleRayCount, nn + 1);
    break;
  case NonzeroRHS:
    fprintf(writing_dex, " %8ld  %5ld    real\n", FeasibleRayCount, nn);
    break;
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(writing_dex, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  fprintf(writing_dex, "end\n");
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteCompletionStatus(writing_log);
    if (PreOrderedRun) {
      fprintf(writing_log, "*set_intersection total#, effective#, loss# = %ld   %ld   %ld\n", 
      count_int, count_int_good, count_int_bad);
    }
    WriteTimes(writing_log);
  }
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout);
    WriteTimes(stdout);
  }
}

void WriteProjRayRecord(FILE *f, RayRecord *RR, long *dbrow)
{
  long i,j,k;
  static double *vec;
  static rowset dbset;
  static long mprev=0;
  
  if (dbset==NULL || mprev<mm){
    set_initialize(&dbset,mm);  
    vec=(double *)calloc(mm, sizeof *vec);
    /* initialize only for the first time or when a larger space is needed */
    mprev=mm;
    if (debug) printf("mprev is replaced with  = %ld\n", mprev);
  }
  for (j = 1; j <= mm-nn; j++){
    i=dbrow[j];
    set_addelem(dbset,i);
    if (debug) printf("index %ld is added to dbset\n",i);
    vec[i-1]=0;
    for (k=1; k<=nn; k++) {
      vec[i-1]+= (RR->Ray[k-1])*AA[j-1][k-1];
    }
    if (debug) printf("vec[ %ld]= %.5E \n",i-1, vec[i-1]);
  }
  i=1;
  for (j = 1; j <= mm; j++){
    if (!set_member(j,dbset)){
      vec[j-1]=RR->Ray[i-1];
      i++;
    }
  }
  fprintf(f, " %2d", 0);
  for (j = 0; j < mm; j++)
    WriteReal(f, vec[j]);
  putc('\n', f);
}

void WriteProjResult(long *dbrow)
{
  RayRecord *TempPtr;

  if (!debug) writing=freopen(outputfile,"w",writing);
  time(&endtime);
  WriteProgramDescription(writing);
  fprintf(writing, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, minput, ninput);
  WriteRunningMode(writing);
  WriteCompletionStatus(writing);
  WriteTimes(writing);
  fprintf(writing, "*FINAL RESULT:\n");
  if (DynamicWriteOn)
    printf("*Computation complete.\n");
  if (DynamicWriteOn)
     printf("*Number of Vertices =%8ld,   Rays =%8ld\n",
	   VertexCount, FeasibleRayCount - VertexCount);
  fprintf(writing, "*Number of Vertices =%8ld,   Rays =%8ld\n",
	 VertexCount, FeasibleRayCount - VertexCount);
  fprintf(writing, "begin\n");
  fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, mm + 1);
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    WriteProjRayRecord(writing, TempPtr, dbrow);
    TempPtr = TempPtr->Next;
  }
  fprintf(writing, "end\n");
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout);
    WriteTimes(stdout);
  }
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteRunningMode(writing_log);
    WriteCompletionStatus(writing_log);
    WriteTimes(writing_log);
  }
  if (AdjacencyOutput != AdjOff){
    if (DynamicWriteOn) printf("Writing the adjacency file...\n");
    WriteAdjacency(writing_adj);
  }
}

void InitialWriting(void)
{
  if (LogWriteOn) {
    fprintf(writing_log, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, minput, ninput);
	fprintf(writing_log,"*Initial set of hyperplanes: ");
    WriteSetElements(writing_log, AddedHyperplanes);
    fprintf(writing_log,"\n");
    fprintf(writing_log, "begin\n");
    fprintf(writing_log, "%5ld %3d\n", mm - nn, 5);
  }
}

void WriteErrorMessages(FILE *f)
{
  switch (Error) {

  case LowColumnRank:
    if (Conversion==IneToExt) {
      fprintf(f,"*Input Error: Input matrix (b, -A) is not column full rank => no vertices.\n");
      break;
    } else {
      fprintf(f,"*Input Error: The polytope (convex hull) is not full dimensional.\n");
      break;
    }
 
  case DimensionTooLarge:
    fprintf(f, "*Input Error: Input matrix is too large:\n");
    fprintf(f, "*Please increase MMAX and/or NMAX in the source code and recompile.\n");
    break;

  case DependentMarkedSet:
    fprintf(f, "*Input Error: Active (equality) rows are linearly dependent.\n");
    fprintf(f, "*Please select independent rows for equality restriction.\n");
    break;

  case FileNotFound:
    fprintf(f, "*Input Error: Specified input file does not exist.\n");
    break;

  case ImproperInputFormat:
    if (Number == Rational) {
     }
    else {
      fprintf(f,"*Input Error: Input format is not correct.\n");
      fprintf(f,"*Format:\n");
      fprintf(f," begin\n");
      fprintf(f,"   m   n  NumberType(real, rational or integer)\n");
      fprintf(f,"   b  -A\n");
      fprintf(f," end\n");
    }
    break;

  case None:
    fprintf(f,"*No Error found.\n");
    break;
  }
}

void WriteLPResult(FILE *f, LPStatusType LPS, double optval,
   Arow sol, Arow dsol, colindex NBIndex, rowrange re, colrange se,
   long iter)
{
  long j;

  time(&endtime);
  fprintf(f,"*cdd LP Solver (Criss-Cross Method) Result\n");  
  fprintf(f,"*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, minput, ninput);
  if (Conversion==LPmax)
    fprintf(f,"*maximization is chosen.\n");  
  else if (Conversion==LPmin)
    fprintf(f,"*minimization is chosen.\n");
  else if (Conversion==InteriorFind){
    fprintf(f,"*inerior point computation is chosen.\n");
    fprintf(f,"*the following is the result of solving the LP:\n");
    fprintf(f,"*   maximize      x_{d+1}\n");
    fprintf(f,"*   s.t.    A x + x_{d+1}  <=  b.\n");
    fprintf(f,"*Thus, the optimum value is zero     if the polyhedron has no interior point.\n");
    fprintf(f,"*      the optimum value is negative if the polyhedron is empty.\n");
    fprintf(f,"*      the LP is dual inconsistent   if the polyhedron admits unbounded inscribing balls.\n");
  }
  if (Conversion==LPmax||Conversion==LPmin){
    fprintf(f,"*Objective function is\n");  
    for (j=0; j<nn; j++){
      if (j>0 && LPcost[j]>=0 ) fprintf(f," +",j);
      if (j>0 && (j % 5) == 0) fprintf(f,"\n");
      WriteReal(f, LPcost[j]);
      if (j>0) fprintf(f," X[%3ld]",j);
    }
    fprintf(f,"\n");
  }

  switch (LPS){
  case Optimal:
    fprintf(f,"*LP status: a dual pair (x, y) of optimal solutions found.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_solution\n");
    for (j=1; j<nn; j++) {
      fprintf(f,"  %3ld : ",j);
      WriteReal(f,sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"  dual_solution\n");
    for (j=1; j<nn; j++){
      fprintf(f,"  %3ld : ",NBIndex[j+1]);
      WriteReal(f,dsol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"  optimal_value : % .9E\n", optval);
    fprintf(f,"end\n");
    break;

  case Inconsistent:
    fprintf(f,"*LP status: LP is inconsistent.\n");
    fprintf(f,"*The positive combination of original inequalities with\n");
    fprintf(f,"*the following coefficients will prove the inconsistency.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  dual_direction\n");
    fprintf(f,"  %3ld : \n",re);
    WriteReal(f,1.0);
    for (j=1; j<nn; j++){
      fprintf(f,"  %3ld : ",NBIndex[j+1]);
      WriteReal(f,dsol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  case DualInconsistent:
    fprintf(f,"*LP status: LP is dual inconsistent.\n");
    fprintf(f,"*The linear combination of columns with\n");
    fprintf(f,"*the following coefficients will prove the dual inconsistency.\n");
    fprintf(f,"*(It is also an unbounded direction for the primal LP.)\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_direction\n");
    for (j=1; j<nn; j++) {
      fprintf(f,"  %3ld : ",j);
      WriteReal(f,sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  default:
    break;
  }
  fprintf(f,"*number of pivot operations = %ld\n", iter);
  WriteTimes(f);
}


void WriteSolvedProblem(FILE *f)
{
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, minput, ninput);
  fprintf(f, "*The input data has been interpreted as the following.\n");
  WriteAmatrix(writing_ver,AA,mm,nn,Inequality);
  WriteRunningMode2(writing_ver);
}

/* end of cddio.c */

