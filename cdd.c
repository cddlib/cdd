#define COPYRIGHT   "Copyright (C) 1994, Komei Fukuda, fukuda@dma.epfl.ch"
#define DDVERSION   "Version C0.33 (January 16, 1994)"

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

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for details.
*/

/* The first version C0.21 was created , November 10, 1993 
   with Dave Gillespie's p2c translator 
   from the Pascal program pdd.p written by Komei Fukuda. 
*/

#include "setoper.h"     /* set operation library header (Dec.8,1993 version or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

long minput, ninput;   /*size of input data [b -A] */
long mm, nn;   /*size of the homogenous system to be solved by dd*/
long projdim;  /*dimension of orthogonal preprojection */
colrange RHScol;
rowrange OBJrow;
colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
rowset MarkedSet, GroundSet;
rowrange Iteration, hh;
rowset AddedHyperplanes, InitialHyperplanes;
long RayCount, FeasibleRayCount, TotalRayCount, VertexCount;
boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
Amatrix AA;
Bmatrix InitialRays;
Arow LPcost;  /* LP cost vector to be maximized  */
RayRecord *ArtificialRay, *FirstRay, *LastRay;
boolean found, inputsuccessful;
HyperplaneOrderType HyperplaneOrder;
AdjacencyTestType AdjacencyTest;
NumberType Number;
InequalityType Inequality;
boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
boolean PartialEnumeration; /* Partial enumeration Switch (TRUE if it is restricted on the intersection of MarkedSet hyperplanes) */
CompStatusType CompStatus;  /* Computation Status */
ConversionType Conversion;
IncidenceOutputType IncidenceOutput;
ErrorType Error;
DataFileType inputfile,outputfile,projfile,icdfile,logfile;
FILE *reading, *writing, *writing_proj, *writing_icd, *writing_log;
time_t starttime, endtime;


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

void SetProjFile(FILE **f)
{
  boolean opened=FALSE;
  char ch;
  long i;

  while (!opened) {
    printf("\n>> Inequality output file of projection variables (*.ine)   : ");
    scanf("%s",projfile);
    if (strcmp(inputfile, projfile)!=0) {
      *f = fopen(projfile, "w");
      printf("write file %s is open\n",projfile);
      opened=TRUE;
    }
    else {
      printf("write file %s must have a name different from inputfile.\n",projfile);
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
  switch (Conversion) {
    case ExtToIne:
      fprintf(f, "*Hull computation is chosen.\n");
      break;
    
    case IneToExt:
      fprintf(f, "*Vertex/Ray enumeration is chosen.\n");
      break;
    
    case LPmax:  case LPmin:
      fprintf(f, "*Linear optimization is chosen.\n");
      break;

    case Projection:
      fprintf(f, "*Preprojection is chosen.\n");
      break;
  
    default: break;
  }
  if (PartialEnumeration) {
    fprintf(f, "*Partial enumeration is chosen.  The permanently active rows are:");
    WriteSetElements(f,MarkedSet);
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

void WriteProjRayRecord(FILE *f, RayRecord *RR, long *dbrow)
{
  long i,j,k;
  double vec[MMAX];
  rowset dbset;

  set_initialize(dbset,rowsetsize);
  for (j = 1; j <= mm-nn; j++){
    i=dbrow[j];
    set_addelem(dbset,i);
    if (debug) printf("index %ld is added to dbset\n",i);
    vec[i-1]=0;
    for (k=1; k<=nn; k++) {
      vec[i-1]+= (RR->Ray[k-1])*AA[j-1][k-1];
    }
    if (debug) printf("vec[ %ld]= %lg \n",i-1, vec[i-1]);
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
	   VertexCount, RayCount - VertexCount);
  fprintf(writing, "*Number of Vertices =%8ld,   Rays =%8ld\n",
	 VertexCount, RayCount - VertexCount);
  fprintf(writing, "begin\n");
  fprintf(writing, " %8ld  %5ld    real\n", RayCount, mm + 1);
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


void DDMain(void)
{
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
    fprintf(f, "*Please increase MMAX and/or NMAX in the source code and recompile.\n");
    break;

  case DependentMarkedSet:
    fprintf(f, "*Input Error: Marked rows are linearly dependent.\n");
    fprintf(f, "*Please select independent rows for partial enumeration.\n");
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
    if (set_subset(MarkedSet, ZSet)){
      if (debug) {
        printf("add an initial ray with zero set:");
        set_write(ZSet);
      }
      AddRay(Vector);
    }
  }
}

void DDEnumerate(void)
{
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
    DDMain();
    WriteDDResult();
  } else {
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}

void PreProjection(void)
{
  rowset subrows1,subrows2,DBrows;
  colset subcols1,subcols2;  /* subcols1:projvars,  subcols2:rest */
  rowrange i;
  colrange j,k;
  colindex pivrow;
  Bmatrix DBinv;  /* dual basis matrix inverse */
  long DBrank;
  
  if (IncidenceOutput == IncSet)
    SetIncidenceFile(&writing_icd);
  if (LogWriteOn)
    SetLogFile(&writing_log);
  if (DynamicWriteOn) {
    WriteRunningMode(stdout);
    WriteRunningMode(writing);
  }
  set_initialize(subrows1,rowsetsize);
  set_initialize(subrows2,rowsetsize);
  set_initialize(subcols1,colsetsize);  /* subcol1 : projvar & RHS columns */
  set_initialize(subcols2,colsetsize);  /* subcol2 : remaining columns */
  SetProjFile(&writing_proj);  
  for (j=1;j<=nn;j++){
    if (set_member(j,projvars) || (j==1 && Inequality==NonzeroRHS))
      set_addelem(subcols1,j);
    else
      set_addelem(subcols2,j);
  }
  for (i=1; i<=mm; i++) set_addelem(subrows1,i);
  if (DynamicWriteOn){
    WriteSubMatrixOfAA(stdout,subrows1,subcols1,Inequality);
  }
  WriteSubMatrixOfAA(writing_proj,subrows1,subcols1,Inequality);
  Inequality=ZeroRHS;
  ReduceAA(subrows1,subcols2);
    /* Extract the submatrix of AA index by subcols2. 
       subcols2 is changed to a consecutive sequence starting from 1 */
  if (debug) {
    WriteAmatrix(stdout,AA,mm,nn,NonzeroRHS);
    WriteAmatrix(writing,AA,mm,nn,NonzeroRHS);
  }
  FindBasis(AA,LeastIndex,DBrows,pivrow,DBinv,&DBrank);
    /* DBrows stores the rows associated with a dual basis */
  if (debug){
    printf("rank of the new (deletion col) matrix is %ld\n", DBrank);
    printf("dual basis rows ="); set_write(DBrows);
  }
  for (j=1;j<=nn;j++) fprintf(writing,"pivot row at col %ld = %ld\n",j, pivrow[j]);
  set_diff(subrows2,subrows1,DBrows); 
    /* subrows2 stores the rows not in DBrows */
  for (j=1; j<=nn;j++){
    if (pivrow[j]==0) {
      set_delelem(subcols2,j);
      fprintf(writing,"Warning: col %ld is a linear combination of the other colums. The column linear dependency must be deleted for ray computation\n",j);
      for (k=j; k<=nn-1; k++){ /* shifting all pivrow information */
        pivrow[j]=pivrow[j+1];
      }
      pivrow[nn]=0;
      nn--;
    }
  }
  if (debug)  {
    printf("rows for ray enumeration:");set_write(subrows2);
    printf("cols for ray enumeration:");set_write(subcols2);
  }
  ReduceAA(subrows2,subcols2); 
    /* subrows2 is changed to a consecutive sequence starting from 1 */
  DualizeAA(DBinv);
  if (Error==DimensionTooLarge) goto _L99;
  if (debug) {
    WriteAmatrix(stdout,AA,mm,nn,ZeroRHS);
    WriteAmatrix(writing,AA,mm,nn,ZeroRHS);
  }
  if (DynamicWriteOn) {
    WriteRunningMode(stdout);
    WriteRunningMode(writing);
  }
  FindInitialRays(InitialHyperplanes, InitialRays, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting();
    DDMain();
    WriteProjResult(pivrow);
  } else {
    _L99:;
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}


void WriteHeading(void)
{
  WriteProgramDescription(stdout);
  printf("------------------------------------------\n");
  printf(" Vertex & Extremal Ray Enumeration for\n");
  printf(" the polyhedron P = { x :  b - A x >= 0 }\n");
  printf("------------------------------------------\n");
}


void main(int argc, char *argv[])
{
  Arow LPsol;
  
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
    switch (Conversion) {
    case ExtToIne: case IneToExt: /* vertex/facets enumeration is chosen */
      DDEnumerate();
      break;
    
    case LPmax:                   /* LP maximization is chosen */
      if (Inequality==NonzeroRHS)
        CrissCrossSolve(AA, InitialRays, OBJrow, RHScol, LPsol);
      else
        printf("Sorry, LP maximization is not implemented for RHS==0.\n");
      break;

    case Projection:              /* preprojection is chosen */
      PreProjection();
      break;
  
    default: break;
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
