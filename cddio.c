/* cddio.c:  Basic Input and Output Procedures for cdd.c
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.61, December 1, 1997
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
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
      printf("\n>> Input file (*.ine or *.ext) : ");
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
  DataFileType newname;

  switch (cflag) {
    case 'o':
      switch (Conversion) {
        case ExtToIne:
          extension=".ine"; break;     /* output file for ine data */
        case IneToExt:   case Projection:
          extension=".ext"; break;     /* output file for ext data */
        case LPmax:  case LPmin:  case InteriorFind:
          extension=".lps"; break;     /* output file for LPmax, LPmin, interior_find */
        default:
        extension=".out";break;
      }
      break;
    case 'a':
      if (Conversion==IneToExt)
        extension=".ead";       /* adjacency file for ext data */
      else
        extension=".iad";       /* adjacency file for ine data */
      break;
    case 'i':
      if (Conversion==IneToExt)
        extension=".ecd";       /* ext incidence file */
      else
        extension=".icd";       /* ine incidence file */
      break;
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
    strcpy(newname,fname);  
    strcat(newname,".old");
    rename(fname,newname);
    if (DynamicWriteOn){
      printf("Default output file %s exists.\n",fname);
      printf("Caution: The old file %s is renamed to %s!\n",fname,newname);
      printf("Create a new file %s.\n",fname);
    }
/*    strcpy(fname,inputfile); 
      strcat(fname,extension);
*/
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

void ProcessCommandLine(rowrange m_input, colrange n_input, char *line)
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
    AdjacencyOutput = AdjacencyList;
    return;
  }
  if (strncmp(line, "#adjacency", 10)==0) {
    AdjacencyOutput = AdjacencyDegree;
    fscanf(reading,"%ld", &msize);
    set_initialize(&CheckPoints, m_input+1);
    for (j=1;j<=msize;j++) {
      fscanf(reading,"%ld",&var);
      if (var <= m_input+1) set_addelem(CheckPoints,var);
    }
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
  if (strncmp(line, "output_reordered", 16)==0 && !OutputReordered){
    printf("* The reordered problem will be generated.\n");
    OutputReordered=TRUE;
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
    set_initialize(&projvars,n_input);
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
    for (j=0;j<n_input;j++) {
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
    for (j=0;j<n_input;j++) {
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
  if (strncmp(line, "dual_simplex", 12)==0 && LPsolver != DualSimplex) {
    if (DynamicWriteOn) printf("Use the dual simplex method.\n");
    LPsolver=DualSimplex;
    return;
  }
  if (strncmp(line, "criss-cross", 9)==0 && LPsolver != CrissCross) {
    if (DynamicWriteOn) printf("Use the crss-cross method.\n");
    LPsolver=CrissCross;
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

void AmatrixInput(rowrange *m_input, colrange *n_input,
  rowrange *m_size, colrange *n_size, Amatrix A, boolean *successful)
{
  long i,j;
  double value;
  long value1,value2;
  boolean found=FALSE,decided=FALSE, fileopened, newformat=FALSE;
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
    else {
      if (strncmp(command, "V-representation", 16)==0) {
        Conversion = ExtToIne; newformat=TRUE;
      } else if (strncmp(command, "H-representation", 16)==0){
          Conversion =IneToExt; newformat=TRUE;
        } else if (strncmp(command, "begin", 5)==0) found=TRUE;
    }
  }
  fscanf(reading, "%ld %ld %s", m_input, n_input, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", *m_input, *n_input, numbtype);
  SetNumberType(numbtype);
  if (Number==Unknown || Number == Rational) {
      goto _L99;
    } 
  Inequality=ZeroRHS;
  for (i=1; i<= *m_input && !decided; i++) {
    fscanf(reading,"%lf", &value);
    if (fabs(value) > zero) {
      Inequality = NonzeroRHS;
      decided=TRUE;
    }
    for (j=2; j<= *n_input; j++) {
      fscanf(reading,"%lf", &value);
    }
    fgets(line,linelenmax,reading);
    if (debug) printf("comments to be skipped: %s\n",line);
    if (debug) putchar('\n');
  }
  if (Inequality==NonzeroRHS) {
  	printf("Nonhomogeneous system with  m = %5ld  n = %5ld\n", *m_input, *n_input);
  	*n_size = *n_input;
  }
  else {
    printf("Homogeneous system with  m = %5ld  n = %5ld\n", *m_input, *n_input);
  	*n_size = *n_input-1;
  }
  if (*n_size > NMAX || *m_input > MMAX) {
    Error = DimensionTooLarge;
    goto _L99;
  }
  
  RestrictedEnumeration = FALSE;
  RelaxedEnumeration = FALSE;
  EqualityIndex=(long *)calloc(*m_input+2, sizeof *EqualityIndex);
  for (i = 0; i <= *m_input+1; i++) EqualityIndex[i]=0;

  while (!feof(reading)) {
    fscanf(reading,"%s", command);
    ProcessCommandLine(*m_input, *n_input, command);
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
  for (i = 1; i <= *m_input; i++) {
    A[i-1]=(double *) calloc(*n_input, sizeof value);
    for (j = 1; j <= *n_input; j++) {
      fscanf(reading, "%lf", &value);
      if (Inequality==NonzeroRHS) 
      	A[i-1][j - 1] = value;
      else if (j>=2) {
        A[i-1][j - 2] = value;
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
      *m_size = *m_input + 1;
      A[*m_size-1]=(double *) calloc(*n_input, sizeof value);
      for (j = 1; j <= *n_input; j++) {   /*artificial row for x_1 >= 0*/
        if (j == 1)
          A[*m_size - 1][j - 1] = 1.0;
        else
          A[*m_size - 1][j - 1] = 0.0;
      }
    } else {
      *m_size = *m_input;
    }
    break;

  case ExtToIne:
    *m_size = *m_input;
    break;

  case LPmax:  case LPmin:
    *m_size = *m_input + 1;
    OBJrow=*m_size;
    RHScol=1L;
    A[*m_size-1]=(double *) calloc(*n_input, sizeof value);
    for (j = 1; j <= *n_input; j++) {   /*objective row */
 	   A[*m_size - 1][j - 1] = LPcost[j-1];
 	}
	break;

  default:
    *m_size = *m_input;
  }
  SetInequalitySets(*m_size, EqualityIndex);
  *successful = TRUE;
_L99: ;
  if (reading!=NULL) fclose(reading);
}

void WriteRayRecord(FILE *f, colrange n_size, RayRecord *RR)
{
  long j;
  double scaler;

  if (Inequality == ZeroRHS) {
    fprintf(f, " %2d", 0);
    for (j = 0; j < n_size; j++)
      WriteReal(f, RR->Ray[j]);
  } 
  else {
    scaler = fabs(RR->Ray[0]);
    if (scaler > zero) {
      if (RR->Ray[0] > 0) {
        if (Conversion == IneToExt) {
	      fprintf(f, " %2d", 1);
	      for (j = 1; j < n_size; j++)
	      WriteReal(f, RR->Ray[j] / scaler);
        } 
        else {
	      /* hull computation is chosen */
          for (j = 0; j < n_size; j++)
	        WriteReal(f, RR->Ray[j]);
	    }
      }
      else {
        /* hull computation must have been chosen, since RR->Ray[0] is negative */
	    for (j = 0; j < n_size; j++)
	      WriteReal(f, RR->Ray[j]);
      }
    } 
    else {
      fprintf(f, " %2d", 0);
      for (j = 1; j < n_size; j++)
        WriteReal(f, RR->Ray[j]);
    }
  }
  if (IncidenceOutput==IncCardinality) {
    fprintf(f," : %1ld", set_card(RR->ZeroSet));
  }
  putc('\n', f);
}


void WriteRayRecord2(FILE *f, colrange n_size, RayRecord *RR)
{
  long j;

  fprintf(f, " Ray = ");
  for (j = 0; j < n_size; j++)
    fprintf(f, "%6.2f", RR->Ray[j]);
  putchar('\n');
  fprintf(f, " ZeroSet =");
  set_fwrite(f, RR->ZeroSet);
  putc('\n', f);
}

void WriteSubMatrixOfA(FILE *f, rowrange m_size, colrange n_size, Amatrix A,
    rowset ChosenRow, colset ChosenCol, InequalityType ineq)
{
  long i,j;

  fprintf(f, "H-representation\nbegin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",m_size, set_card(ChosenCol)+1);
  else
    fprintf(f, "  %ld   %ld    real\n",m_size, set_card(ChosenCol));
  for (i=1; i <= m_size; i++) {
    if (set_member(i,ChosenRow)) {
      if (ineq==ZeroRHS){  /* If RHS==0, output zero first */
        WriteReal(f, 0);
        for (j=1; j <= n_size; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteReal(f, A[i-1][j-1]);
          }
        }
      }
      else {
        for (j=1; j <= n_size; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteReal(f, A[i-1][j-1]);
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

void WriteAmatrix2(FILE *f, Amatrix A, long rowmax, long colmax,
      InequalityType ineq, rowindex OV, rowrange iteration)
{ /* Writing an A matrix with a given ordering and size */
  long i,j,imax;

  imax=iteration;
  if (imax>rowmax) imax=rowmax;
  fprintf(f, "begin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",imax, colmax+1);
  else
    fprintf(f, "  %ld   %ld    real\n",imax, colmax);
  for (i=1; i <= imax; i++) {
    if (ineq==ZeroRHS)
      WriteReal(f, 0);  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j <= colmax; j++) {
      WriteReal(f, A[OV[i]-1][j-1]);
    }
    fprintf(f,"  #%1ld\n", OV[i]);
  }
  fprintf(f, "end\n");
}

void WriteBmatrix(FILE *f, colrange n_size, Bmatrix T)
{
  colrange j1, j2;

  for (j1 = 0; j1 < n_size; j1++) {
    for (j2 = 0; j2 < n_size; j2++) {
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


void WriteIncidence(FILE *f, rowrange m_size, colrange n_size, RayRecord *RR)
{
  rowset cset;
  long zcar;

  set_initialize(&cset,m_size);
  zcar = set_card(RR->ZeroSet);
  switch (IncidenceOutput) {

  case IncCardinality:
    fprintf(f, " : %1ld", zcar);
    break;

  case IncSet:
    if (m_size - zcar >= zcar) {
      fprintf(f, " %1ld : ", zcar);
      set_fwrite(f, RR->ZeroSet);
    } else {
      set_diff(cset, GroundSet, RR->ZeroSet);
      fprintf(f, " %1ld : ", zcar - m_size);
      set_fwrite(f, cset);
    }
    break;

  case IncOff:
    break;
  }
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
  if (Conversion==IneToExt || Conversion==ExtToIne || Conversion==Projection){ 
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
    set_fwrite(f,EqualitySet);
    fprintf(f,"\n");
  }
  if (RelaxedEnumeration) {
    fprintf(f, "*The strict_inequality option is chosen.\n* => Permanently nonactive rows are:");
    set_fwrite(f,NonequalitySet);
    fprintf(f,"\n");
  }
}

void WriteRunningMode2(FILE *f, rowrange m_size, colrange n_size)
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
      for (j=0; j<n_size; j++) WriteReal(f,LPcost[j]);
      fprintf(f,"\n");
      break;

    case LPmin:
      fprintf(f, "minimize\n");
      for (j=0; j<n_size; j++) WriteReal(f,LPcost[j]);
      fprintf(f,"\n");
      break;

    case Projection:
      fprintf(f, "*Preprojection is chosen.\n");
      break;
  
    default: break;
  }
  if (RestrictedEnumeration) {
    fprintf(f, "equality ");
    set_fwrite(f,EqualitySet);
    fprintf(f,"\n");
  }
  if (RelaxedEnumeration) {
    fprintf(f, "strict_inequality ");
    set_fwrite(f,NonequalitySet);
    fprintf(f,"\n");
  }
}

void WriteRunningMode0(FILE *f)
{
  switch (Conversion) {
    case ExtToIne:
      fprintf(f, "V-representation\n");
      break;
    
    case IneToExt: case LPmax: case LPmin: case Projection:
      fprintf(f, "H-representation\n");
      break;

    default: break;
  }
}



void WriteCompletionStatus(FILE *f, rowrange m_size, colrange n_size, rowrange iter)
{
  if (iter<m_size && CompStatus==AllFound) {
    fprintf(f,"*Computation completed at Iteration %4ld.\n", iter);
  } 
  if (CompStatus == RegionEmpty) {
    fprintf(f,"*Computation completed at Iteration %4ld because the region found empty.\n", iter);
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


void WriteAdjacencyDegree(FILE *f, rowrange m_input, colrange n_input,
    rowrange m_size, colrange n_size, Amatrix A, rowrange iter)
{
  RayRecord *RayPtr1, *RayPtr2;
  long pos1, pos2, deg, feasdeg, icd, feasicd;
  double totaldeg=0.0, totalfeasdeg=0.0, totalfeasicd=0.0, averfeasicd, averdeg, averfeasdeg;
  boolean adj;

  switch (Conversion) {
  case IneToExt:
    fprintf(f,"*Adjacency Degree of output (=vertices/rays)\n");
    fprintf(f,"*vertex/ray#, icd#, curr icd#, adj#, feas adj#\n");
    break;
  case ExtToIne:
    fprintf(f,"*Adjacency Degree of output (=inequalities=facets)\n");
    fprintf(f,"*facet#, icd#, curr icd#, adj#, feas adj#\n");
    break;

  default:
    break;
  }
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, m_input, n_input);
  fprintf(f, "*At iteration =  %ld\n", iter);  
  fprintf(f,"begin\n");
  fprintf(f,"  %ld\n",RayCount);
  if (RayCount==0){
    goto _L99;
  }
  LastRay->Next=NULL;
  for (RayPtr1=FirstRay, pos1=1;RayPtr1 != NULL; RayPtr1 = RayPtr1->Next, pos1++){
    for (RayPtr2=FirstRay, pos2=1,deg=0,feasdeg=0;RayPtr2 != NULL; RayPtr2 = RayPtr2->Next, pos2++){
      if (RayPtr1!=RayPtr2){
        CheckAdjacency2(m_size, n_size, A, &RayPtr1, &RayPtr2, &adj);
        if (adj) {
          deg++;
          if (RayPtr2->feasible) feasdeg++;
        }
      }
    }
    if (RayPtr1->feasible) fprintf(f," f"); else fprintf(f," i");
    icd = set_card(RayPtr1->ZeroSet);
    set_int(Face, RayPtr1->ZeroSet, AddedHyperplanes);
    feasicd = set_card(Face);
    totalfeasicd = totalfeasicd + feasicd;
    totalfeasdeg = totalfeasdeg + feasdeg;
    totaldeg = totaldeg + deg;
    fprintf(f," %4ld %4ld %4ld %4ld %4ld", pos1,icd, feasicd,deg, feasdeg);
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
  averfeasicd = totalfeasicd / (double) RayCount;
  averdeg = totaldeg / (double) RayCount;
  averfeasdeg = totalfeasdeg / (double) RayCount;
  fprintf(f,"*average (curr icd#, adj#, feasible adj#) = %5.2f %5.2f %5.2f\n\n", averfeasicd, averdeg, averfeasdeg);
_L99:;
}

void WriteAdjacency(FILE *f, rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A)
{
  RayRecord *RayPtr1, *RayPtr2;
  long pos1, pos2, degree;
  boolean adj;
  node *headnode, *tailnode, *newnode, *prevnode;

  headnode=NULL; tailnode=NULL;
  switch (Conversion) {
  case IneToExt:
    if (AdjacencyOutput==AdjacencyList)
      fprintf(f,"*Adjacency List of output (=vertices/rays)\n");
    break;
  case ExtToIne:
    if (AdjacencyOutput==AdjacencyList)
      fprintf(f,"*Adjacency List of output (=inequalities=facets)\n");
      break;

  default:
    break;
  }
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, m_input, n_input);
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
        CheckAdjacency2(m_size, n_size, A, &RayPtr1, &RayPtr2, &adj);
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

void WriteDDResult(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A, rowrange iter)
{
  RayRecord *TempPtr;

  if (!debug) writing=freopen(outputfile,"w",writing);
  time(&endtime);
  WriteProgramDescription(writing);
  fprintf(writing, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, m_input, n_input);
  WriteRunningMode(writing);
  WriteCompletionStatus(writing, m_size, n_size, iter);
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
  if (Conversion == ExtToIne) fprintf(writing, "H-representation\n");
  if (Conversion == IneToExt) fprintf(writing, "V-representation\n");
  fprintf(writing, "begin\n");
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, n_size + 1);
    break;
  case NonzeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, n_size);
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
	  inputfile, m_input, n_input);
    fprintf(writing_icd, "*cdd output file: %s\n", outputfile);
    fprintf(writing_icd, "begin\n");
    fprintf(writing_icd, "%8ld%5ld%5ld\n", FeasibleRayCount, m_input, m_size);
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(writing, n_size, TempPtr);
      if (IncidenceOutput == IncSet)
        WriteIncidence(writing_icd, m_size, n_size, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  fprintf(writing, "end\n");
  if (Conversion == IneToExt) fprintf(writing, "hull\n");
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout, m_size, n_size, iter);
    WriteTimes(stdout);
  }
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteRunningMode(writing_log);
    WriteCompletionStatus(writing_log, m_size, n_size, iter);
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
    WriteAdjacency(writing_adj, m_input, n_input, m_size, n_size, A);
  }
}

void WriteDecompResult(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, rowrange iter)
{
  RayRecord *TempPtr;
  static int callnumber=0;

  callnumber++;
  if (callnumber==1) {
    WriteProgramDescription(writing_dex);
    fprintf(writing_dex, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, m_input, n_input);
    if (Conversion == ExtToIne)
      fprintf(writing_dex,
        "*Since hull computation is chosen, the output is a minimal inequality system\n");
  }
  fprintf(writing_dex, "begin\n");
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing_dex, " %8ld  %5ld    real\n", FeasibleRayCount, n_size + 1);
    break;
  case NonzeroRHS:
    fprintf(writing_dex, " %8ld  %5ld    real\n", FeasibleRayCount, n_size);
    break;
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(writing_dex, n_size, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  fprintf(writing_dex, "end\n");
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteCompletionStatus(writing_log, m_size, n_size, iter);
    if (PreOrderedRun) {
      fprintf(writing_log, "*set_intersection total#, effective#, loss# = %ld   %ld   %ld\n", 
      count_int, count_int_good, count_int_bad);
    }
    WriteTimes(writing_log);
  }
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout, m_size, n_size, iter);
    WriteTimes(stdout);
  }
}

void WriteProjRayRecord(FILE *f, rowrange m_size, colrange n_size, 
    Amatrix A, RayRecord *RR, long *dbrow)
{
  long i,j,k;
  static double *vec;
  static rowset dbset;
  static long mprev=0;
  
  if (dbset==NULL || mprev<m_size){
    set_initialize(&dbset,m_size);  
    vec=(double *)calloc(m_size, sizeof *vec);
    /* initialize only for the first time or when a larger space is needed */
    mprev=m_size;
    if (debug) printf("mprev is replaced with  = %ld\n", mprev);
  }
  for (j = 1; j <= m_size-n_size; j++){
    i=dbrow[j];
    set_addelem(dbset,i);
    if (debug) printf("index %ld is added to dbset\n",i);
    vec[i-1]=0;
    for (k=1; k<=n_size; k++) {
      vec[i-1]+= (RR->Ray[k-1])*A[j-1][k-1];
    }
    if (debug) printf("vec[ %ld]= %.5E \n",i-1, vec[i-1]);
  }
  i=1;
  for (j = 1; j <= m_size; j++){
    if (!set_member(j,dbset)){
      vec[j-1]=RR->Ray[i-1];
      i++;
    }
  }
  fprintf(f, " %2d", 0);
  for (j = 0; j < m_size; j++)
    WriteReal(f, vec[j]);
  putc('\n', f);
}

void WriteProjResult(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A, long *dbrow, rowrange iter)
{
  RayRecord *TempPtr;

  if (!debug) writing=freopen(outputfile,"w",writing);
  time(&endtime);
  WriteProgramDescription(writing);
  fprintf(writing, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, m_input, n_input);
  WriteRunningMode(writing);
  WriteCompletionStatus(writing, m_size, n_size, iter);
  WriteTimes(writing);
  fprintf(writing, "*FINAL RESULT:\n");
  if (DynamicWriteOn)
    printf("*Computation complete.\n");
  if (DynamicWriteOn)
     printf("*Number of Vertices =%8ld,   Rays =%8ld\n",
	   VertexCount, FeasibleRayCount - VertexCount);
  fprintf(writing, "*Number of Vertices =%8ld,   Rays =%8ld\n",
	 VertexCount, FeasibleRayCount - VertexCount);
  fprintf(writing, "V-representation\nbegin\n");
  fprintf(writing, " %8ld  %5ld    real\n", FeasibleRayCount, m_size + 1);
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    WriteProjRayRecord(writing, m_size, n_size, A, TempPtr, dbrow);
    TempPtr = TempPtr->Next;
  }
  fprintf(writing, "end\n");
  if (DynamicWriteOn) {
    WriteCompletionStatus(stdout, m_size, n_size, iter);
    WriteTimes(stdout);
  }
  if (LogWriteOn) {
    fprintf(writing_log, "end\n");
    WriteRunningMode(writing_log);
    WriteCompletionStatus(writing_log, m_size, n_size, iter);
    WriteTimes(writing_log);
  }
  if (AdjacencyOutput != AdjOff){
    if (DynamicWriteOn) printf("Writing the adjacency file...\n");
    WriteAdjacency(writing_adj, m_input, n_input, m_size, n_size, A);
  }
}

void InitialWriting(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size)
{
  if (LogWriteOn) {
    fprintf(writing_log, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, m_input, n_input);
	fprintf(writing_log,"*Initial set of hyperplanes: ");
    set_fwrite(writing_log, AddedHyperplanes);
    fprintf(writing_log,"\n");
    fprintf(writing_log, "begin\n");
    fprintf(writing_log, "%5ld %3d\n", m_size - n_size, 5);
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


void WriteSolvedProblem(FILE *f, rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A)
{
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, m_input, n_input);
  fprintf(f, "*The input data has been interpreted as the following.\n");
  WriteRunningMode0(writing_ver);
  WriteAmatrix(writing_ver,A,m_size,n_size,Inequality);
  WriteRunningMode2(writing_ver, m_size, n_size);
}

void WriteSolvedSubProblem(FILE *f, rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A, rowindex OV, 
    rowrange iteration)
{ rowrange imax=iteration;

  if (imax<1 || imax>m_size) {
    fprintf(f, "*The subproblem %ld does not exist, thus the whole problem will be ouput.\n",imax);
    imax=m_size;
  }
  fprintf(f, "*cdd input file : %s   (%4ld  x %4ld)\n",
	  inputfile, m_input, n_input);
  fprintf(f, "*The %ld th subproblem with the options:\n",imax);
  WriteRunningMode2(writing_ver, m_size, n_size);
  fprintf(f, "*is the following problem.\n");
  WriteAmatrix2(writing_ver,A,m_size,n_size,Inequality, OV, imax);
  if (Conversion==ExtToIne) fprintf(f, "hull\n");
  fprintf(f, "minindex\n");
}


/* end of cddio.c */

