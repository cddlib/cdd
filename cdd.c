#define COPYRIGHT   "Copyright (C) 1993, Komei Fukuda, fukuda@dma.epfl.ch"
#define DDVERSION   "Version C0.27 (December 8, 1993)"

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


/* cdd: C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
   P= {x :  b - A x >= 0}.  Please read the manual cddman.tex for detail.
*/

/* The first version C0.21 was created , November 10, 1993 
   with Dave Gillespie's p2c translator 
   from the Pascal program pdd.p written by Komei Fukuda. 
*/

#include "setoper.h"      /* set operation library header (Ver. Dec.8,1993 or later) */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define mmax            601   /* USER'S CHOICE: max row size of A plus one */
#define nmax            24    /* USER'S CHOICE: max column size of A plus one */

#define SETBITS 32        /* Important Constant: Number of bits in a long integer    */
#define rowsetsize mmax   /* The size of the column index set */
#define colsetsize nmax   /* The size of the row index set */
#define rowsetblocks (rowsetsize-1)/SETBITS+2   /* Number of (long) blocks for row set */
#define colsetblocks (colsetsize-1)/SETBITS+2   /* Number of (long) blocks for column set */

#define datawidth       10
#define filenamelen     255 
#define wordlenmax      128 
#define linelenmax      255

#define FALSE 0
#define TRUE 1

#define zero            1.0e-5   /*real zero*/

typedef char boolean;
typedef long rowrange;
typedef long colrange;
typedef long rowset[rowsetblocks];
typedef long colset[colsetblocks];
typedef double Amatrix[mmax][nmax];
typedef double Arow[nmax];
typedef double Bmatrix[nmax][nmax];
typedef char DataFileType[filenamelen];

typedef struct RayRecord {
  Arow Ray;
  rowset ZeroSet;
  double ARay;   /*temporary area to store some row of A*Ray*/
  struct RayRecord *Next;
} RayRecord;

typedef struct AdjacencyRecord {
  RayRecord *Ray1, *Ray2;
  struct AdjacencyRecord *next;
} AdjacencyRecord;

typedef enum {
  Combinatorial, Algebraic
} AdjacencyTestType;
typedef enum {
  LargestIndex, LeastIndex, MinCutoff, MaxCutoff, MixCutoff, LexMin, LexMax
} HyperplaneOrderType;
typedef enum {
  Real, Rational, Integer, Unknown
} NumberType;
typedef enum {
  ZeroRHS, NonzeroRHS
} InequalityType;
typedef enum {
  IneToExt, ExtToIne
} ConversionType;
typedef enum {
  IncOff, IncCardinality, IncSet
} IncidenceOutputType;
typedef enum {
  DimensionTooLarge, LowColumnRank, ImproperInputFormat, None
} ErrorType;
typedef enum {
  InProgress, AllFound, RegionEmpty
} CompStatusType;

typedef char LineType[linelenmax];
typedef char WordType[wordlenmax];

long minput, ninput;   /*size of input data [b -A] */
long mm, nn;   /*size of the homogenous system to be solved by dd*/
rowset MarkedSet, GroundSet;
rowrange Iteration, hh;
rowset AddedHyperplanes, InitialHyperplanes;
long RayCount, FeasibleRayCount, TotalRayCount, VertexCount;
boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
Amatrix AA;
Bmatrix InitialRays;
RayRecord *ArtificialRay, *FirstRay, *LastRay;
boolean found, inputsuccessful;
HyperplaneOrderType HyperplaneOrder;
AdjacencyTestType AdjacencyTest;
NumberType Number;
InequalityType Inequality;
boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
CompStatusType CompStatus;     /* Computation Status */
ConversionType Conversion;
IncidenceOutputType IncidenceOutput;
ErrorType Error;
DataFileType inputfile,outputfile,icdfile,logfile;
FILE *reading, *writing, *writing_icd, *writing_log;
time_t starttime, endtime;

void SelectNextHyperplane(HyperplaneOrderType, long *, rowrange *);

void SetInputFile(FILE **f)
{
  boolean opened=FALSE;
  
  while (!opened) {
    printf("\n>> Input file (*.ine) : ");
    scanf("%s",inputfile);
    if ( ( *f = fopen(inputfile, "r") )!= NULL) {
      printf("input file %s is open\n", inputfile);
      opened=TRUE;
    }
    else printf("The file %s not found\n",inputfile);
  }
}

void SetWriteFile(FILE **f)
{
  boolean opened=FALSE;
  char ch;
  long i;

  while (!opened) {
    printf("\n>> Output file (*.ext)   : ");
    scanf("%s",outputfile);
    if (strcmp(inputfile, outputfile)!=0) {
      *f = fopen(outputfile, "w");
      printf("write file %s is open\n",outputfile);
      opened=TRUE;
    }
    else {
      printf("write file %s must have a name different from inputfile.\n",outputfile);
    }
  }
}


void SetIncidenceFile(FILE **f)
{
  boolean opened=FALSE;
  char ch;
  long i;

  while (!opened) {
    printf("\n>> Incidence file (*.icd): ");
    scanf("%s",icdfile);
    if (strcmp(inputfile, icdfile)!=0) {
      *f = fopen(icdfile, "w");
      printf("Incidence file %s is open\n",icdfile);
      opened=TRUE;
    }
    else {
      printf("The name of incidence file %s must be different from inputfile.\n",icdfile);
    }
  }
}

void SetLogFile(FILE **f)
{
  boolean opened=FALSE;
  char ch;
  long i;

  while (!opened) {
    printf("\n>> Log file (*.ddl)      : ");
    scanf("%s",logfile);
    if (strcmp(inputfile, logfile)!=0) {
      *f = fopen(logfile, "w");
      printf("Log file %s is open\n",logfile);
      opened=TRUE;
    }
    else {
      printf("The name of log file %s must be different from inputfile.\n",logfile);
    }
  }
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
  if (strncmp(line, "algebraic", 9)==0) {
    AdjacencyTest = Algebraic;
    return;
  }
  if (strncmp(line, "nondegenerate", 14)==0) {
    NondegAssumed = TRUE;
    return;
  }
  if (strncmp(line, "minindex", 8)==0) {
    HyperplaneOrder = LeastIndex;
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
  if (strncmp(line, "initbasis_at_bottom", 19)==0) {
    InitBasisAtBottom = TRUE;
    return;
  }
  if (strncmp(line, "debug", 5)==0) {
    debug = TRUE;
    return;
  }
}

void AmatrixInput(boolean *successful)
{
  long i,j;
  double value,value1,value2;
  boolean found=FALSE,decided=FALSE;
  char command[wordlenmax], numbtype[wordlenmax], stemp[wordlenmax];

  *successful = FALSE;
  SetInputFile(&reading);
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
  		fscanf(reading,"%ld", &value);
  	}
  }
  if (Inequality==NonzeroRHS) {
  	printf("Nonhomogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
  	nn = ninput;
  }
  else {
    printf("Homogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
  	nn = ninput-1;
  }
  if (nn > nmax || minput > mmax) {
    Error = DimensionTooLarge;
    goto _L99;
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
    for (j = 1; j <= ninput; j++) {
      fscanf(reading, "%lf", &value);
      if (Inequality==NonzeroRHS) 
      	AA[i - 1][j - 1] = value;
      else if (j>=2) {
        AA[i - 1][j - 2] = value;
	  }
	  if (debug)
	    printf("a(%3ld,%5ld) = %10.4lf\n",i,j,value);
    }  /*of j*/
    if (debug)
      putchar('\n');
  }  /*of i*/
  while (!feof(reading)) {
    fscanf(reading,"%s", command);
    ProcessCommandLine(command);
   } 
  if (Conversion == IneToExt) {
    mm = minput + 1;
    for (j = 1; j <= ninput; j++) {   /*artificial row for x_1 >= 0*/
      if (j == 1)
	    AA[mm - 1][j - 1] = 1.0;
      else
	    AA[mm - 1][j - 1] = 0.0;
    }
  } else {
    mm = minput;
  }
  set_initialize(MarkedSet, rowsetsize);
  *successful = TRUE;
_L99: ;
  fclose(reading);
}


long Cardinality(long *S)
{
  rowrange i;
  long car;

  car = 0;
  for (i = 1; i <= mm; i++) {
    if (set_member(i, S))
      car++;
  }
  return car;
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
    fprintf(f," %4ld", Cardinality(RR->ZeroSet));
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


double AValue(double *p, rowrange i)
{
  /*return the ith component of the vector  A x p */
  colrange j;
  double temp;

  temp = 0.0;
  for (j = 0; j < nn; j++)
    temp += AA[i - 1][j] * p[j];
  return temp;
}

void WriteIncidence(FILE *f, RayRecord *RR)
{
  rowset cset;
  long zcar;

  set_initialize(cset,rowsetsize);
  zcar = Cardinality(RR->ZeroSet);
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
  }
  putc('\n', f);
}



void StoreRay(double *p, RayRecord *RR, boolean *feasible)
{
  rowrange i;
  colrange j;
  double temp;
  rowset SET;

  *feasible = TRUE;
  set_initialize(RR->ZeroSet,rowsetsize);
  RR->ARay = 0.0;
  for (j = 0; j < nn; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= mm; i++) {
    temp = AValue(p, i);
    if (fabs(temp) < zero)
      set_addelem(RR->ZeroSet, i);
    if (temp <= -zero)
      *feasible = FALSE;
  }
}


void AddRay(double *p)
{  
  boolean feasible;

  if (FirstRay == NULL) {
    FirstRay = (struct RayRecord *) malloc(sizeof *FirstRay);
    if (debug)
      printf("Create the first ray pointer\n");
    LastRay = FirstRay;
    ArtificialRay->Next = FirstRay;
  } else {
    LastRay->Next = (struct RayRecord *) malloc(sizeof *FirstRay);
    if (debug)
      printf("Create a new ray pointer\n");
    LastRay = LastRay->Next;
  }
  LastRay->Next = NULL;
  RayCount++;
  TotalRayCount++;
  if (DynamicWriteOn) {
    if (TotalRayCount % 100 == 0) {
      fprintf(writing,
	      "*Rays (Total, Currently Active, Feasible) =%8ld%8ld%8ld\n",
	      TotalRayCount, RayCount, FeasibleRayCount);
      printf("*Rays (Total, Currently Active, Feasible) =%8ld%8ld%8ld\n",
	     TotalRayCount, RayCount, FeasibleRayCount);
    }
  }
  StoreRay(p, LastRay, &feasible);
  if (!feasible)
    return;
  FeasibleRayCount++;
  if (fabs(LastRay->Ray[0]) > zero && Inequality == NonzeroRHS)
    VertexCount++;
  if (DynamicRayWriteOn) {
    WriteRayRecord(writing, LastRay);
    WriteRayRecord(stdout, LastRay);
    if (IncidenceOutput==IncSet && writing_icd != NULL) {
    	WriteIncidence(writing_icd, LastRay);
    }
  }
}


void AddArtificialRay(void)
{  
  Arow zerovector;
  long j;
  boolean feasible;

  if (ArtificialRay != NULL) {
    printf("Warning !!!  FirstRay in not nil.  Illegal Call\n");
    return;
  }
  ArtificialRay = (struct RayRecord *) malloc(sizeof *ArtificialRay);
  if (debug)
    printf("Create the artificial ray pointer\n");
  for (j = 0; j < nn; j++)
    zerovector[j] = 0.0;
  StoreRay(zerovector, ArtificialRay, &feasible);
  ArtificialRay->Next = NULL;
}


void Normalize(double *V)
{
  long j;
  double min, temp;

  min = 1.0e+20;
  for (j = 0; j < nn; j++) {
    temp = fabs(V[j]);
    if (temp > zero && temp < min)
      min = temp;
  }
  for (j = 0; j < nn; j++)
    V[j] /= min;
}


void ZeroIndexSet(double *x, long *ZS)
{
  rowrange i;
  rowset tempset;
  double temp;
  rowset SET;

  set_initialize(tempset,rowsetsize);
  for (i = 1; i <= mm; i++) {
    temp = AValue(x, i);
    if (fabs(temp) < zero)
      set_addelem(tempset, i);
  }
  set_copy(ZS, tempset);
}

void CopyBmatrix(double (*T)[nmax], double (*TCOPY)[nmax])
{
  colrange j1, j2;

  for (j1 = 1; j1 <= nn; j1++) {
    for (j2 = 1; j2 <= nn; j2++) {
      TCOPY[j1 - 1][j2 - 1] = T[j1 - 1][j2 - 1];
    }
  }
}


void SelectPivot1(double (*X)[nmax], HyperplaneOrderType roworder,
            rowrange rowmax, long *NopivotRow,
			long *NopivotCol, rowrange *r, colrange *s,
			boolean *selected)
/* Select a position (*r,*s) in the matrix X such that X[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  boolean stop;
  rowrange rtemp;
  colrange stemp;
  rowset rowexcluded;

  stop = FALSE;
  set_initialize(rowexcluded,rowsetsize);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    SelectNextHyperplane(roworder, rowexcluded, &rtemp);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= nn && !*selected) {
        if (!set_member(*s,NopivotCol) && fabs(X[*r - 1][*s - 1]) > zero) {
          *selected = TRUE;
          stop = TRUE;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = TRUE;
    }
  } while (!stop);
}

double TableauEntry(double (*X)[nmax], double (*T)[nmax],
				rowrange r, colrange s)
/* Compute the (r,s) entry of X.T   */
{
  colrange j;
  double temp;
  
  temp=0;
  for (j=0; j< nn; j++) {
    temp = temp + X[r-1][j] * T[j][s-1];
  }
  return temp;
}

void SelectPivot2(double (*X)[nmax], double (*T)[nmax],
            HyperplaneOrderType roworder,
            rowrange rowmax, long *NopivotRow,
			long *NopivotCol, rowrange *r, colrange *s,
			boolean *selected)
/* Select a position (*r,*s) in the matrix X.T such that (X.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  boolean stop;
  rowrange rtemp;
  colrange stemp;
  rowset rowexcluded;
  double Xtemp;

  stop = FALSE;
  set_initialize(rowexcluded,rowsetsize);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    SelectNextHyperplane(roworder, rowexcluded, &rtemp);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= nn && !*selected) {
        Xtemp=TableauEntry(X,T,*r,*s);
        if (!set_member(*s,NopivotCol) && fabs(Xtemp) > zero) {
          *selected = TRUE;
          stop = TRUE;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = TRUE;
    }
  } while (!stop);
}

void GausianColumnPivot1(double (*X)[nmax], rowrange r, colrange s,
				rowrange rowmax)
/* Make a column pivot operation on X on (r,s)  */
{
  long i, j;
  double Xtemp0, Xtemp;

  Xtemp0 = X[r - 1][s - 1];
  for (j = 0; j < nn; j++) {
    if (j + 1 != s) {
      Xtemp = X[r - 1][j];
      for (i = 0; i < rowmax; i++) {
        if (i + 1 != r)
        X[i][j] -= X[i][s - 1] * Xtemp / Xtemp0;
      }
      X[r - 1][j] = 0.0;
    }
  }
  for (i = 0; i < rowmax; i++) {
    if (i + 1 != r)
      X[i][s - 1] /= Xtemp0;
  }
  X[r - 1][s - 1] = 1.0;
}


void GausianColumnPivot2(double (*X)[nmax], double (*T)[nmax],
				rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure does make any pivot on the matrix X
 */
{
  long i, j, j1;
  Arow Rtemp;
  double Xtemp0, Xtemp;

  for (j=1; j<=nn; j++) Rtemp[j-1]=TableauEntry(X, T, r,j);
  Xtemp0 = Rtemp[s-1];
  for (j = 1; j <= nn; j++) {
    if (j != s) {
      Xtemp = Rtemp[j-1];
      for (j1 = 1; j1 <= nn; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= nn; j++)
    T[j-1][s - 1] /= Xtemp0;
}


void InitializeBmatrix(double (*T)[nmax])
{
  colrange j1, j2;

  for (j1 = 1; j1 <= nn; j1++) {
    for (j2 = 1; j2 <= nn; j2++) {
      if (j1 == j2)
        T[j1 - 1][j2 - 1] = 1.0;
      else
        T[j1 - 1][j2 - 1] = 0.0;
    }
  }
}


void WriteBmatrix(FILE *f, double (*T)[nmax])
{
  colrange j1, j2;

  for (j1 = 0; j1 < nn; j1++) {
    for (j2 = 0; j2 < nn; j2++) {
      fprintf(f, "%5.2f ", T[j1][j2]);
    }  /*of j2*/
    putc('\n', f);
  }  /*of j1*/
  putc('\n', f);
}

void ComputeRank(double (*A1)[nmax], long *TargetRows, long *rank)
/* Compute the rank of the submatrix of a Amatrix A1 indexed by RowSelected.
   This procedure does not change the matrix A1.
 */
{
  boolean stop, chosen;
  rowrange i,r;
  colrange s;
  rowset NoPivotRow;
  colset ColSelected;
  Bmatrix Btemp;   /* dual basis inverse */
  
  *rank = 0;
  stop = FALSE;
  set_initialize(NoPivotRow, rowsetsize);
  set_compl(NoPivotRow,TargetRows);
  set_initialize(ColSelected, colsetsize);
  InitializeBmatrix(Btemp);
  if (debug) WriteBmatrix(stdout,Btemp);
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, Btemp, LeastIndex, mm, NoPivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(NoPivotRow, r);
        set_addelem(ColSelected, s);
        (*rank)++;
        GausianColumnPivot2(A1,Btemp, r, s);
        if (debug) {
          WriteBmatrix(stdout,Btemp);
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
      } else {
        stop=TRUE;
      }
  } while (!stop);
}

void ComputeBInverse(double (*A1)[nmax], long lastrow,
       double (*InvA1)[nmax], long *rank)
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset RowSelected;
  colset ColSelected;

  *rank = 0;
  stop = FALSE;
  InitializeBmatrix(InvA1);
  set_initialize(RowSelected, rowsetsize);
  set_initialize(ColSelected, colsetsize);
  do {
    SelectPivot2(A1, InvA1, LeastIndex, lastrow, RowSelected, ColSelected, &r, &s, &chosen);
    if (chosen) {
      (*rank)++;
      if (debug)
        printf("%3ldth pivot on%3d, %3d\n", *rank, r, s);
      GausianColumnPivot2(A1, InvA1, r, s);
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
    } else
      stop = TRUE;
  } while (!stop);
}

void FindBasis(double (*A1)[nmax], 
              HyperplaneOrderType roword, long *lastrow, long *RowSelected,
              double (*BasisInverse)[nmax], long *rank)
{
  boolean stop, chosen;
  long rowsize;
  rowrange r;
  colrange s;
  colset ColSelected;

  *rank = 0;
  stop = FALSE;
  set_initialize(RowSelected, rowsetsize);
  set_initialize(ColSelected, colsetsize);
  InitializeBmatrix(BasisInverse);
  if (debug) WriteBmatrix(stdout,BasisInverse);
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, BasisInverse, roword, mm, RowSelected, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        (*rank)++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
        if (debug) {
          WriteBmatrix(stdout,BasisInverse);
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
      } else {
        stop=TRUE;
      }
      if (*rank==nn) stop = TRUE;
  } while (!stop);
}


void FindInitialRays(long *InitialHyperplanes,
			    double (*InitRays)[nmax], boolean *found)
{
  Bmatrix BInverse;
  rowset CandidateRows;
  long i, lastrow, rank;
  HyperplaneOrderType roworder;

  *found = FALSE;
  set_initialize(CandidateRows, rowsetsize);
  switch (HyperplaneOrder) {
  case LargestIndex:
    roworder = LargestIndex;
    break;

  case LeastIndex:
    roworder = LeastIndex;
    break;

  case MinCutoff:
    roworder = LexMin;
    break;

  case MaxCutoff:
    roworder = LexMin;
    break;

  case MixCutoff:
    roworder = LexMin;
    break;

  case LexMin:
    roworder = LexMin;
    break;

  case LexMax:
    roworder = LexMax;
    break;
  }
  if (InitBasisAtBottom==TRUE) roworder=LargestIndex;
  for (i = 1; i <= mm; i++)
    set_addelem(CandidateRows, i);      /*all rows are candidates for initial cone*/
  if (DynamicWriteOn)
    printf("*Computing an initial set of rays\n");
  lastrow=mm;
  FindBasis(AA, roworder, &lastrow, InitialHyperplanes, BInverse, &rank);
  if (debug) printf("nn = %ld, rank = %ld\n",nn,rank);
  if (rank < nn) {
    Error = LowColumnRank;
    return;
  }
  *found = TRUE;
  if (debug) {
    WriteBmatrix(stdout, BInverse);
  }
  CopyBmatrix(BInverse,InitRays);
  if (debug) 
    WriteBmatrix(stdout, InitRays);
}


void CheckAdjacency1(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  rowset Face, SET;
  long lastrow, rank;

  *adjacent = TRUE;
  set_initialize(SET,rowsetsize);
  set_initialize(Face,rowsetsize);
  set_int(SET, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, SET, AddedHyperplanes);
  if (debug)
    printf("Check adjacency\n");
  if (Cardinality(Face) < nn - 2) {
    *adjacent = FALSE;
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  ComputeRank(AA,Face,&rank);
  if (rank < nn - 2)
    *adjacent = FALSE;
}


void CheckAdjacency2(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  rowset Face;
  RayRecord *TempRay;
  rowset SET;

  *adjacent = TRUE;
  set_initialize(SET,rowsetsize);
  set_initialize(Face,rowsetsize);
  set_int(SET, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, SET, AddedHyperplanes);
  if (debug)
    printf("Check adjacency\n");
  if (Cardinality(Face) < nn - 2) {
    *adjacent = FALSE;
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  TempRay = FirstRay;
  while (TempRay != NULL && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
    	set_int(SET, TempRay->ZeroSet, AddedHyperplanes);
      	if (set_subset(Face, SET))
			*adjacent = FALSE;
    }
    TempRay = TempRay->Next;
  }
}


void CheckEquality(RayRecord **RP1, RayRecord **RP2, boolean *equal)
{
  long j;

  if (debug)
    printf("Check equality of two rays\n");
  *equal = TRUE;
  j = 1;
  while (j <= nn && *equal) {
    if (fabs((*RP1)->Ray[j - 1] - (*RP2)->Ray[j - 1]) > 2.0 * zero)
      *equal = FALSE;
    j++;
  }
  if (*equal)
    printf("Equal records found !!!!\n");
}


void Eliminate(RayRecord **Ptr)
{
  /*eliminate the record pointed by Ptr^.Next*/
  RayRecord *TempPtr;

  if (debug) {
    printf("            Delete:");
    WriteRayRecord(stdout, (*Ptr)->Next);
  }
  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == FirstRay)   /*Update the first pointer*/
    FirstRay = (*Ptr)->Next;
  if (TempPtr == LastRay)   /*Update the last pointer*/
    LastRay = *Ptr;
  free(TempPtr);
  RayCount--;
}


void CreateNewRay(RayRecord *Ptr1, RayRecord *Ptr2, rowrange ii)
{
  /*Create a new ray by taking a linear combination of two rays*/
  colrange j;
  Arow NewRay;
  double v1, v2;

  v1 = fabs(AValue(Ptr1->Ray, ii));
  v2 = fabs(AValue(Ptr2->Ray, ii));
  for (j = 0; j < nn; j++)
    NewRay[j] = Ptr1->Ray[j] * v2 + Ptr2->Ray[j] * v1;
  Normalize(NewRay);
  AddRay(NewRay);
}


void EvaluateARay(rowrange i)
/* Evaluate the ith component of the vector  A x RD.Ray 
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First 
 */
{
  colrange j;
  double temp;
  RayRecord *Ptr, *PrevPtr, *TempPtr;

  Ptr = FirstRay;
  PrevPtr = ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    printf("Error.  Artificial Ray does not point to FirstRay!!!\n");
  }
  while (Ptr != NULL) {
    temp = 0.0;
    for (j = 0; j < nn; j++)
      temp += AA[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
    if ( temp <= -zero && Ptr != FirstRay) {
      /* printf("Moving an infeasible record w.r.t. %ld to FirstRay\n",i); */
      if (Ptr==LastRay) LastRay=PrevPtr;
      TempPtr=Ptr;
      Ptr = Ptr->Next;
      PrevPtr->Next = Ptr;
      ArtificialRay->Next = TempPtr;
      TempPtr->Next = FirstRay;
      FirstRay = TempPtr;
    }
    else {
      PrevPtr = Ptr;
      Ptr = Ptr->Next;
    }
  }
}


void FeasibilityIndices(long *fnum, long *infnum, rowrange i)
{
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  colrange j;
  double temp;
  RayRecord *Ptr;

  *fnum = 0;
  *infnum = 0;
  Ptr = FirstRay;
  while (Ptr != NULL) {
    temp = 0.0;
    for (j = 0; j < nn; j++)
      temp += AA[i - 1][j] * Ptr->Ray[j];
    if (temp >= 0)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}

void SelectNextHyperplane0(long *excluded, rowrange *hnext)
{
  /*A natural way to choose the next hyperplane.  Simply the largest index*/
  long i;
  boolean determined;

  i = mm;
  determined = FALSE;
  do {
    if (set_member(i, excluded))
      i--;
    else
      determined = TRUE;
  } while (!determined && i>=1);
  *hnext = i;
}


void SelectNextHyperplane1(long *excluded, rowrange *hnext)
{
  /*Natural way to choose the next hyperplane.  Simply the least index*/
  long i;
  boolean determined;

  i = 1;
  determined = FALSE;
  do {
    if (set_member(i, excluded))
      i++;
    else
      determined = TRUE;
  } while (!determined);
  if (i<=mm) *hnext = i;
  else *hnext=0;
}


void SelectNextHyperplane2(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmin, fi;   /*feasibility and infeasibility numbers*/

  infmin = RayCount + 1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (inf < infmin) {
	infmin = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (DynamicWriteOn) {
    printf("*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infmin, fi);
    fprintf(writing, "*infeasible rays (min) =%5ld, #feas rays =%5ld\n",
	    infmin, fi);
  }
}


void SelectNextHyperplane3(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax, fi;   /*feasibility and infeasibility numbers*/

  infmax = -1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (inf > infmax) {
	infmax = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (DynamicWriteOn) {
    printf("*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infmax, fi);
    fprintf(writing, "*infeasible rays (max) =%5ld, #feas rays =%5ld\n",
	    infmax, fi);
  }
}


void SelectNextHyperplane4(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with the most unbalanced cut*/
  long i, fea, inf, max, tmax, fi, infi;
      /*feasibility and infeasibility numbers*/

  max = -1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (fea <= inf)
	tmax = inf;
      else
	tmax = fea;
      if (tmax > max) {
	max = tmax;
	fi = fea;
	infi = inf;
	*hnext = i;
      }
    }
  }
  if (!DynamicWriteOn)
    return;
  if (max == fi) {
    printf("*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infi, fi);
    fprintf(writing, "*infeasible rays (min) =%5ld, #feas rays =%5ld\n",
	    infi, fi);
  } else {
    printf("*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infi, fi);
    fprintf(writing, "*infeasible rays (max) =%5ld, #feas rays =%5ld\n",
	    infi, fi);
  }
}


boolean LexSmaller(double *v1, double *v2)
{
  boolean determined, smaller;
  colrange j;

  smaller = FALSE;
  determined = FALSE;
  j = 1;
  do {
    if (fabs(v1[j - 1]-v2[j - 1])>zero) {
      if (v1[j - 1] < v2[j - 1]) {
	    smaller = TRUE;
	  }
      determined = TRUE;
    } else
      j++;
  } while (!(determined) && (j < nn));
  return smaller;
}


boolean LexLarger(double *v1, double *v2)
{
  Arow u1, u2;
  colrange j;

  for (j = 1; j <= nn; j++) {
    u1[j-1] = -v1[j-1];
    u2[j-1] = -v2[j-1];
  }
  return (LexSmaller(u1, u2));
}

void CopyArow(double *vcopy, double *v)
{
 colrange j;

  for (j = 1; j <= nn; j++) {
    vcopy[j-1] = v[j-1];
  }
}


void SelectNextHyperplane5(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-min*/
  long i, minindex;
  colrange j;
  Arow v1, v2;

  minindex = 0;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      for (j = 1; j <= nn; j++)
	    v2[j-1] = AA[i - 1][j - 1];
      if (minindex == 0) {
	    minindex = i;
	    CopyArow(v1,v2);
      } else if (LexSmaller(v2,v1)) {
        minindex = i;
	    CopyArow(v1,v2);
      }
    }
  }
  *hnext = minindex;
}


void SelectNextHyperplane6(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  colrange j;
  Arow v1, v2;

  maxindex = 0;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
       for (j = 1; j <= nn; j++)
        v2[j - 1] = AA[i - 1][j - 1];
      if (maxindex == 0) {
        maxindex = i;
        CopyArow(v1,v2);
      } else if (LexLarger(v2, v1)) {
        maxindex = i;
        CopyArow(v1,v2);
     }
    }
  }
  *hnext = maxindex;
}


void SelectNextHyperplane(HyperplaneOrderType ho, 
         long *excluded, rowrange *hh)
{
  switch (ho) {

  case LargestIndex:
    SelectNextHyperplane0(excluded, hh);
    break;

  case LeastIndex:
    SelectNextHyperplane1(excluded, hh);
    break;

  case MinCutoff:
    SelectNextHyperplane2(excluded,hh);
    break;

  case MaxCutoff:
    SelectNextHyperplane3(excluded, hh);
    break;

  case MixCutoff:
    SelectNextHyperplane4(excluded, hh);
    break;

  case LexMin:
    SelectNextHyperplane5(excluded, hh);
    break;

  case LexMax:
    SelectNextHyperplane6(excluded, hh);
    break;
  }
}


void AddNewHyperplane(rowrange hnew)
{
  RayRecord *RayPtr0, *RayPtr1, *RayPtr2, *RayPtr2s, *RayPtr3;
  long pos1, pos2;
  double prevprogress, progress, value1, value2;
  boolean adj, equal, completed;

  EvaluateARay(hnew);        /*Check feasibility of rays w.r.t. hnew 
                               and put all infeasible ones consecutively */
  RayPtr0 = ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  value1 = FirstRay->ARay;
  if (value1 > -zero) {
    if (RayCount==FeasibleRayCount) CompStatus=AllFound;
    goto _L99;               /* Sicne there is no hnew-infeasible ray and nothing to do */
  }
  else {
    RayPtr2s = RayPtr1->Next;/* RayPtr2s must point the first feasible ray */
    pos2=1;
    while (RayPtr2s!=NULL && RayPtr2s->ARay <= -zero) {
      RayPtr2s = RayPtr2s->Next;
      pos2++;
    }
  }
  if (RayPtr2s==NULL) {
    FirstRay=NULL;
    ArtificialRay->Next=FirstRay;
    CompStatus=RegionEmpty;
    RayCount=0;
    VertexCount=0;
     goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;   /*2nd i-feasible ray to scan and compare with 1st*/
  RayPtr3 = LastRay;    /*Last feasible for scanning*/
  prevprogress=0;
  pos1 = 1;
  completed=FALSE;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    value1 = RayPtr1->ARay;
    value2 = RayPtr2->ARay;
    if (debug) {
      WriteRayRecord2(stdout, RayPtr1);
      WriteRayRecord2(stdout, RayPtr2);
      printf("check feasibility%8.3f%8.3f\n", value1, value2);
    }
    CheckEquality(&RayPtr1, &RayPtr2, &equal);
    if ((value1 >= zero && value2 <= -zero) || (value2 >= zero && value1 <= -zero)) {
      switch (AdjacencyTest) {

      case Algebraic:
		CheckAdjacency1(&RayPtr1, &RayPtr2, &adj);
		break;

      case Combinatorial:
		CheckAdjacency2(&RayPtr1, &RayPtr2, &adj);
		break;
      }
      if (adj)
		CreateNewRay(RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
    if (value1 <= -zero || equal) {
      Eliminate(&RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed=TRUE;
    }
    pos1++;
    progress = 100 * ((double)pos1 / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (DynamicWriteOn && pos1%10==0 && progress-prevprogress>=10 ) {
      printf("*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	     Iteration, mm, pos1, pos2, progress);
      fprintf(writing,
	  "*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	  Iteration, mm, pos1, pos2, progress);
	  prevprogress=progress;
	  fflush(writing);
	  if (writing_icd != NULL) fflush(writing_icd);
    }
  }
  if (RayCount==FeasibleRayCount) CompStatus=AllFound;
  _L99:;
}


void WriteProgramDescription(FILE *f)
{
  fprintf(f, "* cdd: Double Description Method C-Code:%s\n", DDVERSION);
  fprintf(f,"* %s\n",COPYRIGHT);
}

void WriteRunningMode(FILE *f)
{
  switch (HyperplaneOrder) {

  case LeastIndex:
    fprintf(f, "*HyperplaneOrder: LeastIndex\n");
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
  }
  switch (AdjacencyTest) {

  case Combinatorial:
    fprintf(f, "*AdjacencyTest: Combinatorial\n");
    break;

  case Algebraic:
    fprintf(f, "*AdjacencyTest: Algebraic\n");
    break;
  }
  if (NondegAssumed) {
    fprintf(f, "*Degeneracy preknowledge for computation: NondegenerateAssumed\n");
   }
  else {
    fprintf(f, "*Degeneracy preknowledge for computation: None (possible degeneracy)\n");
  }
  if (Conversion == ExtToIne) {
    fprintf(f, "*Hull computation is chosen.\n");
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


void WriteFinalResult(void)
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
    if (DynamicWriteOn)
      printf("*Number of Vertices =%8ld,   Rays =%8ld\n",
	     VertexCount, RayCount - VertexCount);
    fprintf(writing, "*Number of Vertices =%8ld,   Rays =%8ld\n",
	    VertexCount, RayCount - VertexCount);
  } else {
    if (DynamicWriteOn)
      printf("*Number of Facets =%8ld\n", RayCount);
    fprintf(writing, "*Number of Facets =%8ld\n", RayCount);
  }
  fprintf(writing, "begin\n");
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", RayCount, nn + 1);
    break;
  case NonzeroRHS:
    fprintf(writing, " %8ld  %5ld    real\n", RayCount, nn);
    break;
  }
  if (IncidenceOutput == IncSet) {
    writing_icd=freopen(icdfile,"w",writing_icd);
    fprintf(writing_icd, "*Input File:%.*s   (%4ld  x %4ld)\n",
	  filenamelen, inputfile, minput, ninput);
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
    }
    fprintf(writing_icd, "begin\n");
    fprintf(writing_icd, "%8ld%5ld%5ld\n", RayCount, minput, mm);
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    WriteRayRecord(writing, TempPtr);
    if (IncidenceOutput == IncSet)
      WriteIncidence(writing_icd, TempPtr);
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
  if (IncidenceOutput == IncSet)
    fprintf(writing_icd, "end\n");
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
  if (DynamicWriteOn) {
    printf("*Initial set of hyperplanes: ");
    WriteSetElements(stdout, AddedHyperplanes);
    putchar('\n');
    fprintf(writing,"*Initial set of hyperplanes: ");
    WriteSetElements(writing, AddedHyperplanes);
    fprintf(writing,"\n");
  }
}


void Main(void)
{
  long SET[mmax / 32 + 2];

  Iteration = nn + 1;
  while (Iteration <= mm) {
    SelectNextHyperplane(HyperplaneOrder, AddedHyperplanes, &hh);
    if (DynamicWriteOn) {
      fprintf(writing,
	      "*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	      Iteration, hh);
      printf("*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	     Iteration, hh);
    }
    AddNewHyperplane(hh);
    if (CompStatus==AllFound||CompStatus==RegionEmpty) goto _L99;
    if (LogWriteOn)
      fprintf(writing_log, "%3ld %5ld %6ld %6ld %6ld\n",
	      Iteration, hh, TotalRayCount, RayCount, FeasibleRayCount);
    set_addelem(AddedHyperplanes, hh);
    Iteration++;
  }
  _L99:;
  WriteFinalResult();
}


void WriteErrorMessages(FILE *f)
{
  switch (Error) {

  case LowColumnRank:
    if (Conversion==IneToExt) {
      fprintf(f,"*Input Error: Input matrix (b, -A) is not column full rank => no vertices and rays.\n");
      break;
    } else {
      fprintf(f,"*Input Error: Input matrix A is not column full rank.=> no vertices and rays.\n");
      break;
    }
      

  case DimensionTooLarge:
    fprintf(f, "*Input Error: Input matrix is too large:\n");
    fprintf(f, "*Please increase mmax and/or nmax in the source code and recompile.\n");
    break;

  case ImproperInputFormat:
    if (Number == Rational) {
      fprintf(f,"*Sorry, rational input is not supported by this version of cdd.\n");
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
  }
}

void Initialization(void)
/* Initialization of global set variables */
{
  set_initialize(InitialHyperplanes,rowsetsize);
  set_initialize(AddedHyperplanes,rowsetsize);
  set_initialize(GroundSet, rowsetsize);
  Error=None;
  CompStatus=InProgress;
}

void InitialDataSetup(void)
{
  long j, r;
  Arow Vector;
  rowset ZSet;

  time(&starttime);
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  VertexCount = 0;
  ArtificialRay = NULL;
  FirstRay = NULL;
  LastRay = NULL;
  set_initialize(ZSet,rowsetsize);
  for (j = 1; j <= mm; j++)
    set_addelem(GroundSet, j);
  AddArtificialRay();
  Iteration = nn;   /*Initially,we have already  nn  hyperplanes */
  set_copy(AddedHyperplanes, InitialHyperplanes);
  for (r = 0; r < nn; r++) {
    for (j = 0; j < nn; j++)
      Vector[j] = InitialRays[j][r];
    Normalize(Vector);
    ZeroIndexSet(Vector, ZSet);
    if (set_subset(MarkedSet, ZSet))
      AddRay(Vector);
  }
}


void DefaultOptionSetup(void)
{
  debug = FALSE;
  DynamicWriteOn = TRUE;
  DynamicRayWriteOn = TRUE;
  LogWriteOn = FALSE;
  HyperplaneOrder = LexMin;
  AdjacencyTest = Combinatorial;
  NondegAssumed = FALSE;
  Conversion = IneToExt;
  IncidenceOutput = IncOff;
  InitBasisAtBottom = FALSE;
}


void WriteHeading(void)
{
  WriteProgramDescription(stdout);
  printf("------------------------------------------\n");
  printf(" Vertex & Extremal Ray Enumeration for\n");
  printf(" the polyhedron P = { x :  b - A x >= 0 }\n");
  printf("------------------------------------------\n");
}


main(int argc, char *argv[])
{
  writing_log = NULL;
  writing_icd = NULL;
  writing = NULL;
  reading = NULL;
  WriteHeading();
  DefaultOptionSetup();
  Initialization();
  AmatrixInput(&inputsuccessful);
  if (inputsuccessful) {
    SetWriteFile(&writing);
    if (IncidenceOutput == IncSet)
      SetIncidenceFile(&writing_icd);
    if (LogWriteOn)
      SetLogFile(&writing_log);
    if (DynamicWriteOn) {
       WriteRunningMode(stdout);
       WriteRunningMode(writing);
    }
    FindInitialRays(InitialHyperplanes, InitialRays, &found);
    if (found) {
      InitialDataSetup();
      InitialWriting();
      Main();
    } else {
      WriteErrorMessages(stdout);
      WriteErrorMessages(writing);
    }
  } else {
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
  if (writing != NULL)
    fclose(writing);
  if (writing_icd != NULL)
    fclose(writing_icd);
  if (writing_log != NULL)
    fclose(writing_log);
}

/* End. */
