/* cdd.c: Main program of the sofware cdd
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.60, August 21, 1996
   Standard ftp site: ifor13.ethz.ch(129.132.154.13), Directory: pub/fukuda/cdd
*/

/* cdd : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

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

/* The first version C0.21 was created on November 10,1993 
   with Dave Gillespie's p2c translator 
   from the Pascal program pdd.p written by Komei Fukuda. 
*/

#include "setoper.h" 
  /* set operation library header (March 16, 1995 version or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* #include <profile.h>    THINK C PROFILER */
/* #include <console.h>    THINK C PROFILER */

long projdim;  /*dimension of orthogonal preprojection */
colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1, CheckPoints;
rowindex EqualityIndex;  
  /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
  TotalRayCount, VertexCount, ZeroRayCount;
long EdgeCount, TotalEdgeCount;
long count_int=0,count_int_good=0,count_int_bad=0;
boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
colrange RHScol;   /* LP RHS column */
rowrange OBJrow;   /* LP OBJ row */
Arow LPcost;  /* LP cost vector to be maximized  */
RayRecord *ArtificialRay, *FirstRay, *LastRay;
RayRecord *PosHead, *ZeroHead, *NegHead, *PosLast, *ZeroLast, *NegLast;
AdjacencyRecord *Edges[MMAX];  /* adjacency relation storage for iteration k */
boolean RecomputeRowOrder, found, inputsuccessful;
HyperplaneOrderType HyperplaneOrder;
AdjacencyTestType AdjacencyTest;
NumberType Number;
InequalityType Inequality;
boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
boolean RestrictedEnumeration; /* Restricted enumeration switch (TRUE if it is restricted on the intersection of EqualitySet hyperplanes) */
boolean RelaxedEnumeration; /* Relaxed enumeration switch (TRUE if NonequalitySet inequalities must be satisfied with strict inequality) */
boolean RowDecomposition; /* Row decomposition enumeration switch */
boolean VerifyInput; /* Verification switch for the input data */
boolean PreOrderedRun; 
  /* TRUE if the rows are ordered before execution & all necessary adjacencies are stored */
CompStatusType CompStatus;  /* Computation Status */
ConversionType Conversion;
LPsolverType LPsolver;
IncidenceOutputType IncidenceOutput;
AdjacencyOutputType AdjacencyOutput;
ErrorType Error;
FileInputModeType FileInputMode;
DataFileType inputfile,ifilehead,ifiletail,
  outputfile,projfile,icdfile,adjfile,logfile,dexfile,verfile;
FILE *reading, *writing, *writing_proj, 
  *writing_icd, *writing_adj,*writing_log,*writing_dex,*writing_ver,*reading_dex;
time_t starttime, endtime;
unsigned int rseed=1;  /* random seed for random row permutation */

void DefaultOptionSetup(void)
{
  debug = FALSE;
  DynamicWriteOn = TRUE;
  DynamicRayWriteOn = TRUE;
  LogWriteOn = FALSE;
  HyperplaneOrder = LexMin;
  AdjacencyTest = Combinatorial;
  NondegAssumed = FALSE;
  RecomputeRowOrder=TRUE;
  PreOrderedRun=TRUE;
  VerifyInput=FALSE;
  Conversion = IneToExt;
  LPsolver = DualSimplex;
    
  IncidenceOutput = IncOff;
  AdjacencyOutput = AdjOff;
  InitBasisAtBottom = FALSE;
}

void CheckAdjacency1(rowrange m_size, colrange n_size, Amatrix A, rowindex ordervec,
    RayRecord **RP1, RayRecord **RP2, boolean *adjacent)
{
  long rank;

  *adjacent = TRUE;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (debug)
    printf("Check adjacency\n");
  if (set_card(Face)< n_size - 2) {
    *adjacent = FALSE;
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  ComputeRank(m_size, n_size, A, Face, ordervec, &rank);
  if (rank < n_size - 2){
    *adjacent = FALSE;
  }
}

void CheckAdjacency2(rowrange m_size, colrange n_size, Amatrix A,
    RayRecord **RP1, RayRecord **RP2, boolean *adjacent)
{
  RayRecord *TempRay;
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  *adjacent = TRUE;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (localdebug){
    printf("Check adjacency of\n");
    WriteRayRecord(stdout, n_size, *RP1);
    WriteRayRecord(stdout, n_size, *RP2);    
  }
  if (set_card(Face)< n_size - 2) {
    *adjacent = FALSE;
    if (localdebug) {
      printf("non adjacent: set_card(face) %ld < %ld = n_size.\n",
        set_card(Face),n_size);
    }
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  TempRay = FirstRay;
  while (TempRay != NULL && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
    	set_int(Face1, TempRay->ZeroSet, AddedHyperplanes);
      	if (set_subset(Face, Face1)) *adjacent = FALSE;
    }
    TempRay = TempRay->Next;
  }
}

void Eliminate(colrange n_size, RayRecord **Ptr)
{
  /*eliminate the record pointed by Ptr^.Next*/
  RayRecord *TempPtr;

  if (debug) {
    printf("            Delete:");
    WriteRayRecord(stdout, n_size, (*Ptr)->Next);
  }
  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == FirstRay)   /*Update the first pointer*/
    FirstRay = (*Ptr)->Next;
  if (TempPtr == LastRay)   /*Update the last pointer*/
    LastRay = *Ptr;
  free(TempPtr->Ray);          /* free the ray vector memory */
  set_free(TempPtr->ZeroSet);  /* free the ZeroSet memory */
  free(TempPtr);   /* free the RayRecord structure memory */
  RayCount--; 
}

void CompileDecompResult(rowrange m_size, colrange n_size, Amatrix A, rowindex ordervec)
{
  long i,j,k;
  double value;
  long mray,nray;
  char numbtype[wordlenmax],command[wordlenmax];
  static double *vec;
  static long mprev=0;
  boolean localdebug=FALSE;
  
  if (mprev<m_size){
    vec=(double *)calloc(m_size, sizeof *vec);
    /* initialize only for the first time or when a larger space is needed */
    mprev=m_size;
    if (localdebug) printf("mprev is replaced with  = %ld\n", mprev);
  }
  AddArtificialRay(m_size, n_size, A, ordervec);
  if (writing_dex != NULL){
    fclose(writing_dex);
    if (DynamicWriteOn) printf("closing the file %s\n",dexfile);
  }
  reading_dex = fopen(dexfile, "r");
  for (i=1; i<=m_size-n_size+2;i++){
    found=FALSE;
    while (!found)
    {
      if (fscanf(reading_dex,"%s",command)==EOF) {
       Error=ImproperInputFormat;
       goto _L99;
      }
      else if (strncmp(command, "begin", 5)==0) {
        found=TRUE;
      }
    }
    fscanf(reading_dex, "%ld %ld %s", &mray, &nray, numbtype);
    if (localdebug) printf("decomp size = %ld x %ld\nNumber Type = %s\n", mray, nray, numbtype);
    for (k=1; k<=mray;k++){
      for (j=1; j<=nray; j++){
        fscanf(reading_dex, "%lf", &value);
        if (Inequality==NonzeroRHS) {
          vec[j - 1] = value;
        } else if (j>=2) {
          vec[j - 2] = value;
        }
        if (localdebug) WriteReal(stdout, value);
      }
      if (localdebug) printf("\n");
      AddRay(m_size, n_size, A, vec, ordervec);
    }
  }
_L99:;
}


void DDInit(rowrange m_size, colrange n_size, Amatrix A, Bmatrix InitialRays,
    rowindex ordervec)
{
  Error=None;
  CompStatus=InProgress;
  SetInequalitySets(m_size, EqualityIndex);
  set_initialize(&InitialHyperplanes,m_size);
  set_initialize(&AddedHyperplanes,m_size);
  set_initialize(&WeaklyAddedHyperplanes,m_size);
  set_initialize(&Face, m_size);   /* used in CheckAdjacency  */
  set_initialize(&Face1, m_size);  /* used in CheckAdjacency  */
  ComputeRowOrderVector(m_size, n_size, A, ordervec, HyperplaneOrder);
  RecomputeRowOrder=FALSE;
  InitializeBmatrix(n_size, InitialRays);
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  WeaklyFeasibleRayCount = 0;
  VertexCount = 0;
  EdgeCount=0; /* active edge count */
  TotalEdgeCount=0; /* active edge count */
}

void DDMain(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A, rowrange *Iteration, rowindex ordervec)
{
  rowrange hh;

  *Iteration = n_size + 1;
  while (*Iteration <= m_size) {
    SelectNextHyperplane(m_size, n_size, A, HyperplaneOrder, 
       WeaklyAddedHyperplanes, &hh, &RecomputeRowOrder, ordervec);
    if (DynamicWriteOn) {
      fprintf(writing,
	     "*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	      *Iteration, hh);
      printf("*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	     *Iteration, hh);
    }
    if (set_member(hh,NonequalitySet)){  /* Skip the row hh */
      if (DynamicWriteOn) {
        fprintf(writing,"*The row # %3ld should be inactive and thus skipped.\n", hh);
        printf("*The row # %3ld should be inactive and thus skipped.\n", hh);
      }
      set_addelem(WeaklyAddedHyperplanes, hh);
    } else {
      if (PreOrderedRun)
        AddNewHyperplane2(m_size, n_size, A, hh, *Iteration, ordervec);
      else
        AddNewHyperplane1(m_size, n_size, A, hh, *Iteration, ordervec);
      set_addelem(AddedHyperplanes, hh);
      set_addelem(WeaklyAddedHyperplanes, hh);
    }
    if (LogWriteOn)
      fprintf(writing_log, "%3ld %5ld %6ld %6ld %6ld\n",
        *Iteration, hh, TotalRayCount, RayCount, FeasibleRayCount);
    if (AdjacencyOutput==AdjacencyDegree && set_member(*Iteration,CheckPoints))
      WriteAdjacencyDegree(writing_adj, m_input, n_input, m_size, n_size, A, *Iteration);
    if (CompStatus==AllFound||CompStatus==RegionEmpty) {
      set_addelem(AddedHyperplanes, hh);
      goto _L99;
    }
    (*Iteration)++;
  }
  _L99:;
}



void Initialization(int ARGC, char *ARGV[])
/* Initialization of global variables */
{
  writing_log = NULL;
  writing_icd = NULL;
  writing_adj = NULL;
  writing_proj = NULL;
  writing = NULL;
  reading = NULL;

  Error=None;
  CompStatus=InProgress;
  if (ARGC>1){
    FileInputMode=Auto;
    strcpy(inputfile,ARGV[1]);
  }
  else{
    FileInputMode=Manual;
  }
}

void InitialDataSetup(rowrange m_size, colrange n_size, 
    Amatrix A, Bmatrix InitialRays, colindex InitialRayIndex, rowindex ordervec)
{
  long j, r;
  Arow Vector1,Vector2;
  rowset ZSet;

  RecomputeRowOrder=FALSE;
  ArtificialRay = NULL;
  FirstRay = NULL;
  LastRay = NULL;
  set_initialize(&ZSet,m_size);
  AddArtificialRay(m_size, n_size, A, ordervec);
  set_copy(AddedHyperplanes, InitialHyperplanes);
  set_copy(WeaklyAddedHyperplanes, InitialHyperplanes);
  UpdateRowOrderVector(m_size, n_size, InitialHyperplanes, ordervec);
  for (r = 1; r <= n_size; r++) {
    for (j = 0; j < n_size; j++){
      Vector1[j] = InitialRays[j][r-1];
      Vector2[j] = -InitialRays[j][r-1];
    }
    Normalize(n_size, Vector1);
    Normalize(n_size, Vector2);
    ZeroIndexSet(m_size, n_size, A, Vector1, ZSet);
    if (set_subset(EqualitySet, ZSet)){
      if (debug) {
        printf("add an initial ray with zero set:");
        set_write(ZSet);
      }
      AddRay(m_size, n_size, A, Vector1, ordervec);
      if (InitialRayIndex[r]==0) {
        AddRay(m_size, n_size, A, Vector2, ordervec);
        if (debug) {
          printf("and add its negative also.\n");
        }
      }
    }
  }
  CreateInitialEdges(m_size, n_size, ordervec);
  set_free(ZSet);
}

void DDEnumerate(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A)
{
  Bmatrix InitRays;
  colindex InitRayIndex; /* 0 if the corr. ray is for generator of an extreme line */ 
  rowrange Iteration=0;
  rowindex OrderVector;
  
  if (IncidenceOutput == IncSet && writing_icd == NULL)
    SetWriteFile(&writing_icd, icdfile, 'i', "incidence");
  if (AdjacencyOutput != AdjOff && writing_adj == NULL)
    SetWriteFile(&writing_adj, adjfile, 'a', "adjacency");
  if (LogWriteOn && writing_log == NULL)
    SetWriteFile(&writing_log, logfile, 'l', "log");
  OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);
  DDInit(m_size, n_size, A, InitRays, OrderVector);
  time(&starttime);
  FindInitialRays(m_size, n_size, A, OrderVector, 
    InitialHyperplanes, InitRays, InitRayIndex, &found);
  if (found) {
    InitialDataSetup(m_size, n_size, A, InitRays, InitRayIndex, OrderVector);
    InitialWriting(m_input, n_input, m_size, n_size);
    DDMain(m_input, n_input, m_size, n_size, A, &Iteration, OrderVector);
    WriteDDResult(m_input, n_input, m_size, n_size, A, Iteration);
    FreeDDMemory(OrderVector);
  } else {
    WriteDDResult(m_input, n_input, m_size, n_size, A, Iteration);
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}

void DecompositionCore(rowrange m_input, colrange n_input, 
  rowrange m_size, colrange n_size, 
  Amatrix A, Bmatrix InitRays, colindex InitRayIndex, 
  rowrange *Iteration)
{
  rowindex OrderVector;

  OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);  
  DDInit(m_size, n_size, A, InitRays, OrderVector);
  time(&starttime);
  FindInitialRays(m_size, n_size, A, OrderVector, InitialHyperplanes, 
    InitRays, InitRayIndex, &found);
  if (found) {
    InitialDataSetup(m_size, n_size, A, InitRays, InitRayIndex, OrderVector);
    InitialWriting(m_input, n_input, m_size, n_size);
    DDMain(m_input, n_input, m_size, n_size, A, Iteration, OrderVector);
    time(&endtime);
    WriteDecompResult(m_input, n_input, m_size, n_size, *Iteration);
    FreeDDMemory(OrderVector);
  } else {
    time(&endtime);
    WriteDecompResult(m_input, n_input, m_size, n_size, *Iteration);
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}

void DDRowDecomposition(rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix A)
{
  rowrange i,k;
  long FeasibleRaySum=0;
  time_t starttime_save;
  Bmatrix InitRays;
  colindex InitRayIndex; /* 0 if the corr. ray is for generator of an extreme line */ 
  rowrange Iteration;
  rowindex OrderVector;
  
  time(&starttime_save);
  if (RowDecomposition && writing_dex == NULL)
    SetWriteFile(&writing_dex, dexfile, 'd', "decomposition");
  if (IncidenceOutput == IncSet && writing_icd == NULL)
    SetWriteFile(&writing_icd, icdfile, 'i', "incidence");
  if (AdjacencyOutput != AdjOff && writing_adj == NULL)
    SetWriteFile(&writing_adj, adjfile, 'a', "adjacency");
  if (LogWriteOn && writing_log == NULL)
    SetWriteFile(&writing_log, logfile, 'l', "log");
  RestrictedEnumeration=TRUE;
  RelaxedEnumeration=TRUE;
  for (i = 0; i <= m_size; i++) EqualityIndex[i]=0;
  for (k = 1; k <= m_size-n_size+2; k++){
    EqualityIndex[k]=1;   /* Equality for k-th inequality */
    if (k>=2) EqualityIndex[k-1]=-1;  /* Strict inequality for 1,2,...,(k-1)st inequalities */
    if (DynamicWriteOn) {
      fprintf(writing, "* Decomposition problem number =%3ld(/%3ld)\n", k, m_size-n_size+2);
      fprintf(stdout, "* Decomposition problem number =%3ld(/%3ld)\n", k, m_size-n_size+2);
    }
    DecompositionCore(m_input, n_input, m_size, n_size, A, InitRays, InitRayIndex, &Iteration);
    FeasibleRaySum=FeasibleRaySum+FeasibleRayCount;
  }
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing_dex, "*Total outputs = %8ld  %5ld    real\n",FeasibleRaySum, n_size + 1);
    fprintf(stdout, "*Total outputs = %8ld  %5ld    real\n",FeasibleRaySum, n_size + 1);
    break;
  case NonzeroRHS:
    fprintf(writing_dex, "*Total outputs = %8ld  %5ld    real\n", FeasibleRaySum, n_size);
    fprintf(stdout, "*Total outputs = %8ld  %5ld    real\n", FeasibleRaySum, n_size);
    break;
  }
  OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);
  DDInit(m_size, n_size, A, InitRays, OrderVector);
  for (i = 0; i <= m_size; i++) {
    EqualityIndex[i]=0;
    set_addelem(AddedHyperplanes,i);
  }
  CompileDecompResult(m_size, n_size, A, OrderVector);
  starttime=starttime_save;
  RestrictedEnumeration=FALSE;
  RelaxedEnumeration=FALSE;
  WriteDDResult(m_input, n_input, m_size, n_size, A, Iteration);
  FreeDDMemory(OrderVector);
}

void PreProjection(rowrange m_input, colrange n_input,
    rowrange m_size, colrange n_size, Amatrix A)
{
  rowset subrows1,subrows2,DBrows;
  colset subcols1,subcols2;  /* subcols1:projvars,  subcols2:rest */
  rowrange i;
  colrange j,k;
  colindex pivrow;
  Bmatrix DBinv;  /* dual basis matrix inverse */
  long DBrank;
  Bmatrix InitRays;
  colindex InitRayIndex; /* 0 if the corr. ray is for generator of an extreme line */ 
  rowrange Iteration;
  rowindex OrderVector;

  time(&starttime);
  if (IncidenceOutput == IncSet)
    SetWriteFile(&writing_icd, icdfile, 'i', "incidence");
  if (AdjacencyOutput != AdjOff)
    SetWriteFile(&writing_adj, adjfile, 'a', "adjacency");
  if (LogWriteOn)
    SetWriteFile(&writing_log, logfile, 'l', "log");
  set_initialize(&subrows1,m_size);
  set_initialize(&subrows2,m_size);
  set_initialize(&DBrows,m_size);
  set_initialize(&subcols1,n_size);  /* subcol1 : projvar & RHS columns */
  set_initialize(&subcols2,n_size);  /* subcol2 : remaining columns */
  SetWriteFile(&writing_proj, projfile, 'p', "preprojection variable subsystem");
  for (j=1;j<=n_size;j++){
    if (set_member(j,projvars) || (j==1 && Inequality==NonzeroRHS))
      set_addelem(subcols1,j);
    else
      set_addelem(subcols2,j);
  }
  for (i=1; i<=m_size; i++) set_addelem(subrows1,i);
  if (DynamicWriteOn){
    WriteSubMatrixOfA(stdout, m_size, n_size, A, subrows1,subcols1,Inequality);
  }
  WriteSubMatrixOfA(writing_proj, m_size, n_size, A, subrows1,subcols1,Inequality);
  Inequality=ZeroRHS;
  ReduceA(&m_size, &n_size, A, subrows1,subcols2);
    /* Extract the submatrix of A index by subcols2. 
       subcols2 is changed to a consecutive sequence starting from 1 */
  if (debug) {
    WriteAmatrix(stdout, A, m_size, n_size, NonzeroRHS);
    WriteAmatrix(writing, A, m_size, n_size, NonzeroRHS);
  }
  PreOrderedRun=FALSE;
  InitializeBmatrix(n_size, DBinv);
  FindBasis(m_size, n_size, A, MinIndex, OrderVector, DBrows,pivrow,DBinv,&DBrank);
    /* DBrows stores the rows associated with a dual basis */
  if (debug){
    printf("rank of the new (deletion col) matrix is %ld\n", DBrank);
    printf("dual basis rows ="); set_write(DBrows);
  }
  for (j=1;j<=n_size;j++) fprintf(writing,"pivot row at col %ld = %ld\n",j, pivrow[j]);
  set_diff(subrows2,subrows1,DBrows); 
    /* subrows2 stores the rows not in DBrows */
  for (j=1; j<=n_size;j++){
    if (pivrow[j]==0) {
      set_delelem(subcols2,j);
      fprintf(writing,"Warning: col %ld is a linear combination of the other colums. The column linear dependency must be deleted for ray computation\n",j);
      for (k=j; k<=n_size-1; k++){ /* shifting all pivrow information */
        pivrow[j]=pivrow[j+1];
      }
      pivrow[n_size]=0;
      n_size--;
    }
  }
  if (debug)  {
    printf("rows for ray enumeration:");set_write(subrows2);
    printf("cols for ray enumeration:");set_write(subcols2);
  }
  ReduceA(&m_size, &n_size, A, subrows2,subcols2); 
    /* subrows2 is changed to a consecutive sequence starting from 1 */
  DualizeA(&m_size, &n_size, A, DBinv);
  if (Error==DimensionTooLarge) goto _L99;
  if (debug) {
    WriteAmatrix(stdout,A,m_size,n_size,ZeroRHS);
    WriteAmatrix(writing,A,m_size,n_size,ZeroRHS);
  }
  if (DynamicWriteOn) {
    WriteRunningMode(stdout);
    WriteRunningMode(writing);
  }
  OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);
  DDInit(m_size, n_size, A, InitRays,OrderVector);
  FindInitialRays(m_size, n_size, A, OrderVector,
    InitialHyperplanes, InitRays, InitRayIndex, &found);
  if (found) {
    InitialDataSetup(m_size, n_size, A, InitRays, InitRayIndex, OrderVector);
    InitialWriting(m_input, n_input, m_size, n_size);
    DDMain(m_input, n_input, m_size, n_size, A,  &Iteration, OrderVector);
    WriteProjResult(m_input, n_input, m_size, n_size, A, pivrow, Iteration);
  } else {
    _L99:;
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
  set_free(subrows1);
  set_free(subrows2);
  set_free(DBrows);
  set_free(subcols1);
  set_free(subcols2);
}

void LPMain(rowrange m_input, colrange n_input,
    rowrange m_size, colrange n_size, Amatrix A)
{
  colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  double ov;  /* LP optimum value */
  long LPiter;
  boolean UsePrevBasis=FALSE;
  Bmatrix BasisInverse;
  dp_LPSolverType solver;
  dp_LPConversionType lpconv;
  dp_ErrorType error;
  dp_LPStatusType LPStatus;

  time(&starttime);
  dp_InitializeBmatrix(n_size, BasisInverse);
  OBJrow=m_size; RHScol=1L;

  if (Inequality==ZeroRHS){
    printf("Sorry, LP optimization is not implemented for RHS==0.\n");
    goto _L99;
  }
  switch (Conversion) {
    case LPmax: lpconv=dp_LPmax;break;
    case LPmin: lpconv=dp_LPmin;break;
    default: lpconv=dp_LPmax;
  }  
  switch (LPsolver) {
    case DualSimplex: solver=dp_DualSimplex;break;
    case CrissCross:  solver=dp_CrissCross;break;
    default:          solver=dp_DualSimplex;
  }
  dp_LPSolve(lpconv, solver, m_size, n_size, A, BasisInverse, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter, &error);
  dp_WriteLPResult(writing, lpconv, solver, m_size, n_size, A, OBJrow, RHScol,
    LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
  if (DynamicWriteOn)
    dp_WriteLPResult(stdout, lpconv, solver, m_size, n_size, A, OBJrow, RHScol,
      LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
_L99:;
}

void InteriorFindMain(rowrange m_input, colrange n_input,
    rowrange m_size, colrange n_size, Amatrix A)
{
  colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  double ov;  /* LP optimum value */
  long LPiter;
  Bmatrix BasisInverse;
  boolean UsePrevBasis;
  dp_LPSolverType solver;
  dp_LPConversionType lpconv;
  dp_ErrorType error;
  dp_LPStatusType LPStatus;

  lpconv=dp_LPmax;
  switch (LPsolver) {
    case DualSimplex: solver=dp_DualSimplex;break;
    case CrissCross:  solver=dp_CrissCross;break;
    default:          solver=dp_DualSimplex;
  }

  if (Inequality==ZeroRHS){
    printf("Sorry, find_interior is not implemented for RHS==0.\n");
    goto _L99;
  }
  EnlargeAforInteriorFinding(&m_input, &n_input, A);
  dp_InitializeBmatrix(n_size, BasisInverse);

  time(&starttime);
  OBJrow=m_size; RHScol=1;
  UsePrevBasis=FALSE;

  fprintf(writing,"*inerior point computation is chosen.\n");
  fprintf(writing,"*the following is the result of solving the LP:\n");
  fprintf(writing,"*   maximize      x_{d+1}\n");
  fprintf(writing,"*   s.t.    A x + x_{d+1}  <=  b.\n");
  fprintf(writing,"*Thus, the optimum value is zero     if the polyhedron has no interior point.\n");
  fprintf(writing,"*      the optimum value is negative if the polyhedron is empty.\n");
  fprintf(writing,"*      the LP is dual inconsistent   if the polyhedron admits unbounded inscribing balls.\n");
 
  dp_LPSolve(lpconv, solver, m_size, n_size, A, BasisInverse, OBJrow, RHScol, UsePrevBasis, 
      &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter, &error);
  dp_WriteLPResult(writing, lpconv, solver, m_size, n_size, A, OBJrow, RHScol,
    LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
  if (DynamicWriteOn)
    dp_WriteLPResult(stdout, lpconv, solver, m_size, n_size, A, OBJrow, RHScol,
      LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
_L99:;
}


void main(int argc, char *argv[])
{
  rowrange minput;
  colrange ninput;   /*size of input data [b -A] */
  rowrange mm;
  colrange nn;   /*size of the homogenous system to be solved by cdd */
  Amatrix AA;

  WriteHeading();
  DefaultOptionSetup();
  Initialization(argc, argv);
  AmatrixInput(&minput, &ninput, &mm, &nn, AA, &inputsuccessful);

  /* InitProfile(200,200);                THINK C PROFILER */
  /* cecho2file("cdd profile", 0,stdout); THINK C PROFILER */
  /* _trace=0;                            THINK C PROFILER */

  if (inputsuccessful) {
    SetWriteFile(&writing,outputfile,'o',"output");
    if (VerifyInput){
      SetWriteFile(&writing_ver,verfile,'v',"input verification");
      WriteSolvedProblem(writing_ver, minput, ninput, mm, nn, AA);
      fclose(writing_ver);
      if (DynamicWriteOn) printf("closing the file %s\n",verfile);
     }
    if (DynamicWriteOn) {
      WriteRunningMode(stdout);
      WriteRunningMode(writing);
    }
    switch (Conversion) {
    case ExtToIne: case IneToExt: /* vertex/facets enumeration is chosen */
      if (RowDecomposition) DDRowDecomposition(minput, ninput, mm, nn, AA);
      else DDEnumerate(minput, ninput, mm, nn, AA);
      break;
    
    case LPmax:  case LPmin:      /* LP is chosen */
      LPMain(minput, ninput, mm, nn, AA);
      
      break;

    case Projection:              /* preprojection is chosen */
      PreProjection(minput, ninput, mm, nn, AA);
      break;

    case InteriorFind:      /* Interior point search is chosen */
      InteriorFindMain(minput, ninput, mm, nn, AA);
      break;
  
    default: break;
    }
  } else {
    WriteErrorMessages(stdout);
    if (writing!=NULL) WriteErrorMessages(writing);
  }
  if (writing != NULL){
    fclose(writing);
    if (DynamicWriteOn) printf("closing the file %s\n",outputfile);
  }
  if (writing_dex != NULL){
    fclose(writing_dex);
    if (DynamicWriteOn) printf("closing the file %s\n",dexfile);
  }
  if (writing_icd != NULL){
    fclose(writing_icd);
    if (DynamicWriteOn) printf("closing the file %s\n",icdfile);
  }
  if (writing_adj != NULL){
    fclose(writing_adj);
    if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
  }
  if (writing_log != NULL){
    fclose(writing_log);
    if (DynamicWriteOn) printf("closing the file %s\n",logfile);
  }
  /* DumpProfile();    THINK C PROFILER */
  /* exit(0);          THINK C PROFILER */
}

/* end of cdd.c */
