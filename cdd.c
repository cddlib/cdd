/* cdd.c: Main program of the sofware cdd
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.55a, December 18, 1994
   Standard ftp site: ftp.epfl.ch,  Directory: incoming/dma
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
  /* set operation library header (Dec 5, 1994 version or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* #include <profile.h>    THINK C PROFILER */
/* #include <console.h>    THINK C PROFILER */

long minput, ninput;   /*size of input data [b -A] */
long mm, nn;   /*size of the homogenous system to be solved by dd*/
long projdim;  /*dimension of orthogonal preprojection */
colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1;
rowrange Iteration, hh;
rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
rowindex EqualityIndex;  
  /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
  TotalRayCount, VertexCount, ZeroRayCount;
long EdgeCount, TotalEdgeCount;
long count_int=0,count_int_good=0,count_int_bad=0;
boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
Amatrix AA;
Bmatrix InitialRays;
colindex InitialRayIndex; /* 0 if the corr. ray is for generator of an extreme line */ 
colrange RHScol;   /* LP RHS column */
rowrange OBJrow;   /* LP OBJ row */
LPStatusType LPStatus;
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
boolean QPivotOn; /* QPivot Switch (TRUE if Q-Pivot scheme of Jack Edmonds is chosen) */
CompStatusType CompStatus;  /* Computation Status */
ConversionType Conversion;
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
    
  IncidenceOutput = IncOff;
  AdjacencyOutput = AdjOff;
  InitBasisAtBottom = FALSE;
  QPivotOn=FALSE;
}

void CheckAdjacency1(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  long rank;

  *adjacent = TRUE;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (debug)
    printf("Check adjacency\n");
  if (set_card(Face)< nn - 2) {
    *adjacent = FALSE;
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = TRUE;
  	return;
  }
  ComputeRank(AA,Face,&rank);
  if (rank < nn - 2){
    *adjacent = FALSE;
  }
}

void CheckAdjacency2(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  RayRecord *TempRay;
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  *adjacent = TRUE;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (localdebug){
    printf("Check adjacency of\n");
    WriteRayRecord(stdout, *RP1);
    WriteRayRecord(stdout, *RP2);    
  }
  if (set_card(Face)< nn - 2) {
    *adjacent = FALSE;
    if (localdebug) {
      printf("non adjacent: set_card(face) %ld < %ld = nn.\n",
        set_card(Face),nn);
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
  free(TempPtr->Ray);          /* free the ray vector memory */
  set_free(TempPtr->ZeroSet);  /* free the ZeroSet memory */
  free(TempPtr);   /* free the RayRecord structure memory */
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
  if (determined) 
    *hnext = i;
  else
    *hnext = 0;
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
  } while (!determined && i<=mm);
  if (determined) 
    *hnext = i;
  else 
    *hnext=0;
}

void SelectNextHyperplane2(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmin, fi=0;   /*feasibility and infeasibility numbers*/

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
    fprintf(writing, 
      "*infeasible rays (min) =%5ld, #feas rays =%5ld\n", infmin, fi);
  }
}

void SelectNextHyperplane3(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax, fi=0;   /*feasibility and infeasibility numbers*/

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
    fprintf(writing,
      "*infeasible rays (max) =%5ld, #feas rays =%5ld\n", infmax, fi);
  }
}

void SelectNextHyperplane4(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with the most unbalanced cut*/
  long i, fea, inf, max, tmax, fi=0, infi=0;
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
  double *v1, *v2;

  minindex = 0;
  v1 = NULL;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
	  v2 = AA[i - 1];
      if (minindex == 0) {
	    minindex = i;
	    v1=v2;
      } else if (LexSmaller(v2,v1,nn)) {
        minindex = i;
	    v1=v2;
      }
    }
  }
  *hnext = minindex;
}


void SelectNextHyperplane6(long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  double *v1, *v2;

  maxindex = 0;
  v1 = NULL;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      v2= AA[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1=v2;
      } else if (LexLarger(v2, v1, nn)) {
        maxindex = i;
        v1=v2;
     }
    }
  }
  *hnext = maxindex;
}

long Partition(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  double *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  i=p-1;
  j=r+1;
  while (TRUE){
    do{
      j--;
    } while (LexLarger(A[OV[j]-1],x,nmax));
    do{
      i++;
    } while (LexSmaller(A[OV[i]-1],x,nmax));
    if (i<j){
      ovi=OV[i];
      OV[i]=OV[j];
      OV[j]=ovi;
    }
    else{
      return j;
    }
  }
}

void QuickSort(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  long q;
  
  if (p < r){
    q = Partition(OV, p, r, A, nmax);
    QuickSort(OV, p, q, A, nmax);
    QuickSort(OV, q+1, r, A, nmax);
  }
}

void LineShellingOrder(rowindex OV, double *z, double *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  double temp1,temp2,infinity=10.0e+20;
  static double *beta;
  static long mlast=0;
  boolean localdebug=FALSE;
  
  if ( mlast<mm ){
    if (beta!=NULL) free(beta);
    beta=(double *)calloc(mm, sizeof *beta);
    /* initialize only for the first time or when last mm is smaller */
    if (localdebug) printf("Initialize the component beta[%ld],\n", i-1);
    mlast=mm;
  }
  for (i=1; i<= mm; i++) beta[i-1]=AA[i-1][0]; /* store the first column in beta */
  for (i=1; i<= mm; i++){
    temp1 = 0.0;
    temp2 = 0.0;
    for (j = 1; j <= nn; j++){
      temp1 += AA[i - 1][j-1] * z[j-1];
      temp2 += AA[i - 1][j-1] * d[j-1];
    }
    if (abs(temp1)>zero) AA[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > 0) AA[i-1][0]= infinity;
    else AA[i-1][0]= -infinity;
     /* use the first column of AA tentatively */
  }
  if (localdebug) 
    for (i=1; i<= mm; i++){
      printf("set AA[%ld] = %lg\n", i, AA[i-1][0]);
    }
  QuickSort(OV, 1, mm, AA, 1);
  for (i=1; i<= mm; i++) {
    AA[i-1][0]=beta[i-1]; 
     /* restore the first column of AA */ 
    if (localdebug) printf("restore AA[%ld] with %lg\n", i, AA[i-1][0]);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void RandomPermutation(rowindex OV, long t, unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  boolean localdebug=FALSE;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=j*u +1;
    k=xk;
    if (localdebug) printf("u=%lg, k=%ld, r=%lg, randmax= %lg\n",u,k,r,rand_max);
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
    if (localdebug) printf("row %ld is exchanged with %ld\n",j,k); 
  }
}

void ComputeRowOrderVector(rowindex OV, HyperplaneOrderType ho)
{
  long i,itemp,j;
  Arow zvec, dvec;
  
  OV[0]=0;
  switch (ho){
  case MaxIndex:
    for(i=1; i<=mm; i++) OV[i]=mm-i+1;
    break;

  case MinIndex: 
    for(i=1; i<=mm; i++) OV[i]=i;
    break;

  case LexMin: case MinCutoff: case MixCutoff: case MaxCutoff:
    for(i=1; i<=mm; i++) OV[i]=i;
    QuickSort(OV, 1, mm, AA, nn);
    break;

  case LexMax:
    for(i=1; i<=mm; i++) OV[i]=i;
    QuickSort(OV, 1, mm, AA, nn);
    for(i=1; i<=mm/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[mm-i+1];
      OV[mm-i+1]=itemp;
    }
    break;

  case RandomRow:
    for(i=1; i<=mm; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    RandomPermutation(OV, mm, rseed);
    break;

  case LineShelling:
    for(i=1; i<=mm; i++) OV[i]=i;
    zvec[0]=1;
    dvec[0]=0;
    if (rseed<=0) rseed=1;
    srand(rseed);
    for(j=2; j<=nn; j++){
      zvec[j-1]=0;
      dvec[j-1]=nn-j+1;
      /* dvec[j-1]=rand(); */
    }
    LineShellingOrder(OV, zvec, dvec);
    break;
  }
}

void UpdateRowOrderVector(long *PriorityRows)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  rowrange i,j,k,i1,j1,oj;
  long rr;
  boolean found, localdebug=FALSE;
  
  if (debug) localdebug=TRUE;
  found=TRUE;
  rr=set_card(PriorityRows);
  if (localdebug) set_write(PriorityRows);
  for (i=1; i<=rr; i++){
    found=FALSE;
    for (j=i; j<=mm && !found; j++){
      oj=OrderVector[j];
      if (set_member(oj, PriorityRows)){
        found=TRUE;
        if (localdebug) printf("%ldth in sorted list (row %ld) is in PriorityRows\n", j, oj);
        j1=j;
      }
    }
    if (found){
      if (j1>i) {
        /* shift everything lower: OV[i]->OV[i+1]..OV[j1-1]->OV[j1] */
        for (k=j1; k>=i; k--) OrderVector[k]=OrderVector[k-1];
        OrderVector[i]=oj;
        if (localdebug){
          printf("OrderVector updated to:\n");
          for (j = 1; j <= mm; j++) printf(" %2ld", OrderVector[j]);
          printf("\n");
        }
      }
    } else {
      printf("UpdateRowOrder: Error.\n");
      goto _L99;
    }
  }
_L99:;
}

void SelectPreorderedNext(long *excluded, rowindex OV, rowrange *hnext)
{
  rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=mm && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k, excluded)) *hnext=k ;
  }
}

void SelectNextHyperplane(HyperplaneOrderType ho, 
         long *excluded, rowrange *hh, boolean *RefreshOrderVector)
{
  if (PreOrderedRun){
    if (debug) {
      printf("debug SelectNextHyperplane: Use PreorderNext\n");
    }
    SelectPreorderedNext(excluded, OrderVector, hh);
  }
  else {
    if (debug) {
      printf("debug SelectNextHyperplane: Use DynamicOrderedNext\n");
    }

    switch (ho) {

    case MaxIndex:
      SelectNextHyperplane0(excluded,hh);
      break;

    case MinIndex:
      SelectNextHyperplane1(excluded,hh);
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

    default:
      SelectNextHyperplane0(excluded,hh);
      break;
    }
  }
}

void CompileDecompResult(void)
{
  long i,j,k;
  double value;
  long mray,nray;
  char numbtype[wordlenmax],command[wordlenmax],line[linelenmax];
  static double *vec;
  static long mprev=0;
  boolean localdebug=FALSE;
  
  if (mprev<mm){
    vec=(double *)calloc(mm, sizeof *vec);
    /* initialize only for the first time or when a larger space is needed */
    mprev=mm;
    if (localdebug) printf("mprev is replaced with  = %ld\n", mprev);
  }
  AddArtificialRay();
  if (writing_dex != NULL){
    fclose(writing_dex);
    if (DynamicWriteOn) printf("closing the file %s\n",dexfile);
  }
  reading_dex = fopen(dexfile, "r");
  for (i=1; i<=mm-nn+2;i++){
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
      AddRay(vec);
    }
  }
_L99:;
}



void DDInit(void)
{
  long i;

  Error=None;
  CompStatus=InProgress;
  SetInequalitySets(EqualityIndex);
  set_initialize(&InitialHyperplanes,mm);
  set_initialize(&AddedHyperplanes,mm);
  set_initialize(&WeaklyAddedHyperplanes,mm);
  set_initialize(&Face, mm);   /* used in CheckAdjacency  */
  set_initialize(&Face1, mm);  /* used in CheckAdjacency  */
  OrderVector=(long *)calloc(mm+1, sizeof *OrderVector);
  ComputeRowOrderVector(OrderVector, HyperplaneOrder);
  RecomputeRowOrder=FALSE;
  InitializeBmatrix(InitialRays);
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  WeaklyFeasibleRayCount = 0;
  VertexCount = 0;
  EdgeCount=0; /* active edge count */
  TotalEdgeCount=0; /* active edge count */
}

void DDMain(void)
{
  rowrange i;

  Iteration = nn + 1;
  while (Iteration <= mm) {
    SelectNextHyperplane(HyperplaneOrder, WeaklyAddedHyperplanes, &hh, &RecomputeRowOrder);
    if (DynamicWriteOn) {
      fprintf(writing,
	     "*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	      Iteration, hh);
      printf("*----------  Iteration =%3ld :   add  row # %3ld ----------\n",
	     Iteration, hh);
    }
    if (set_member(hh,NonequalitySet)){  /* Skip the row hh */
      if (DynamicWriteOn) {
        fprintf(writing,"*The row # %3ld should be inactive and thus skipped.\n", hh);
        printf("*The row # %3ld should be inactive and thus skipped.\n", hh);
      }
      set_addelem(WeaklyAddedHyperplanes, hh);
    } else {
      if (PreOrderedRun)
        AddNewHyperplane2(hh);
      else
        AddNewHyperplane1(hh);
      set_addelem(AddedHyperplanes, hh);
      set_addelem(WeaklyAddedHyperplanes, hh);
    }
    if (LogWriteOn)
      fprintf(writing_log, "%3ld %5ld %6ld %6ld %6ld\n",
	      Iteration, hh, TotalRayCount, RayCount, FeasibleRayCount);
    if (CompStatus==AllFound||CompStatus==RegionEmpty) {
      set_addelem(AddedHyperplanes, hh);
      goto _L99;
    }
    Iteration++;
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

void InitialDataSetup(void)
{
  long j, r;
  Arow Vector1,Vector2;
  rowset ZSet;

  RecomputeRowOrder=FALSE;
  ArtificialRay = NULL;
  FirstRay = NULL;
  LastRay = NULL;
  set_initialize(&ZSet,mm);
  AddArtificialRay();
  Iteration = nn;   /*Initially,we have already  nn  hyperplanes */
  set_copy(AddedHyperplanes, InitialHyperplanes);
  set_copy(WeaklyAddedHyperplanes, InitialHyperplanes);
  UpdateRowOrderVector(InitialHyperplanes);
  for (r = 1; r <= nn; r++) {
    for (j = 0; j < nn; j++){
      Vector1[j] = InitialRays[j][r-1];
      Vector2[j] = -InitialRays[j][r-1];
    }
    Normalize(Vector1);
    Normalize(Vector2);
    ZeroIndexSet(Vector1, ZSet);
    if (set_subset(EqualitySet, ZSet)){
      if (debug) {
        printf("add an initial ray with zero set:");
        set_write(ZSet);
      }
      AddRay(Vector1);
      if (InitialRayIndex[r]==0) {
        AddRay(Vector2);
        if (debug) {
          printf("and add its negative also.\n");
        }
      }
    }
  }
  CreateInitialEdges();
  set_free(ZSet);
}

void DDEnumerate(void)
{
  if (IncidenceOutput == IncSet && writing_icd == NULL)
    SetWriteFile(&writing_icd, icdfile, 'i', "incidence");
  if (AdjacencyOutput != AdjOff && writing_adj == NULL)
    SetWriteFile(&writing_adj, adjfile, 'a', "adjacency");
  if (LogWriteOn && writing_log == NULL)
    SetWriteFile(&writing_log, logfile, 'l', "log");
  DDInit();
  time(&starttime);
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting();
    DDMain();
    WriteDDResult();
    FreeDDMemory();
  } else {
    WriteDDResult();
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}

void DecompositionCore(void)
{
  DDInit();
  time(&starttime);
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting();
    DDMain();
    time(&endtime);
    WriteDecompResult();
    FreeDDMemory();
  } else {
    time(&endtime);
    WriteDecompResult();
    WriteErrorMessages(stdout);
    WriteErrorMessages(writing);
  }
}

void DDRowDecomposition(void)
{
  rowrange i,k;
  long FeasibleRaySum=0;
  time_t starttime_save;
  
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
  for (i = 0; i <= mm; i++) EqualityIndex[i]=0;
  for (k = 1; k <= mm-nn+2; k++){
    EqualityIndex[k]=1;   /* Equality for k-th inequality */
    if (k>=2) EqualityIndex[k-1]=-1;  /* Strict inequality for 1,2,...,(k-1)st inequalities */
    if (DynamicWriteOn) {
      fprintf(writing, "* Decomposition problem number =%3ld(/%3ld)\n", k, mm-nn+2);
      fprintf(stdout, "* Decomposition problem number =%3ld(/%3ld)\n", k, mm-nn+2);
    }
    DecompositionCore();
    FeasibleRaySum=FeasibleRaySum+FeasibleRayCount;
  }
  switch (Inequality) {
  case ZeroRHS:
    fprintf(writing_dex, "*Total outputs = %8ld  %5ld    real\n",FeasibleRaySum, nn + 1);
    fprintf(stdout, "*Total outputs = %8ld  %5ld    real\n",FeasibleRaySum, nn + 1);
    break;
  case NonzeroRHS:
    fprintf(writing_dex, "*Total outputs = %8ld  %5ld    real\n", FeasibleRaySum, nn);
    fprintf(stdout, "*Total outputs = %8ld  %5ld    real\n", FeasibleRaySum, nn);
    break;
  }
  DDInit();
  for (i = 0; i <= mm; i++) {
    EqualityIndex[i]=0;
    set_addelem(AddedHyperplanes,i);
  }
  CompileDecompResult();
  starttime=starttime_save;
  RestrictedEnumeration=FALSE;
  RelaxedEnumeration=FALSE;
  WriteDDResult();
  FreeDDMemory();
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
    SetWriteFile(&writing_icd, icdfile, 'i', "incidence");
  if (AdjacencyOutput != AdjOff)
    SetWriteFile(&writing_adj, adjfile, 'a', "adjacency");
  if (LogWriteOn)
    SetWriteFile(&writing_log, logfile, 'l', "log");
  time(&starttime);
  set_initialize(&subrows1,mm);
  set_initialize(&subrows2,mm);
  set_initialize(&DBrows,mm);
  set_initialize(&subcols1,nn);  /* subcol1 : projvar & RHS columns */
  set_initialize(&subcols2,nn);  /* subcol2 : remaining columns */
  SetWriteFile(&writing_proj, projfile, 'p', "preprojection variable subsystem");
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
  PreOrderedRun=FALSE;
  InitializeBmatrix(DBinv);
  FindBasis(AA,MinIndex,DBrows,pivrow,DBinv,&DBrank);
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
  DDInit();
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
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
  set_free(subrows1);
  set_free(subrows2);
  set_free(DBrows);
  set_free(subcols1);
  set_free(subcols2);
}

void LPMain(void)
{
  colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  double ov;  /* LP optimum value */
  long LPiter;

  time(&starttime);
  if (Inequality==ZeroRHS){
    printf("Sorry, LP optimization is not implemented for RHS==0.\n");
    goto _L99;
  }
  if (Conversion==LPmax){
    CrissCrossMaximize(AA, InitialRays, OBJrow, RHScol, 
      &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  }
  else if (Conversion==LPmin){
    CrissCrossMinimize(AA, InitialRays, OBJrow, RHScol, 
      &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  }
  WriteLPResult(writing, LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
  if (DynamicWriteOn)
    WriteLPResult(stdout,LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
_L99:;
}

void InteriorFindMain(void)
{
  colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  double ov;  /* LP optimum value */
  long LPiter;

  if (Inequality==ZeroRHS){
    printf("Sorry, find_interior is not implemented for RHS==0.\n");
    goto _L99;
  }
  EnlargeAAforInteriorFinding();
  time(&starttime);
  OBJrow=mm; RHScol=1;
  CrissCrossMaximize(AA, InitialRays, OBJrow, RHScol, 
    &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  WriteLPResult(writing, LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
  if (DynamicWriteOn)
    WriteLPResult(stdout,LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
_L99:;
}


void main(int argc, char *argv[])
{
  WriteHeading();
  DefaultOptionSetup();
  Initialization(argc, argv);
  AmatrixInput(&inputsuccessful);

  /* InitProfile(200,200);                THINK C PROFILER */
  /* cecho2file("cdd profile", 0,stdout); THINK C PROFILER */
  /* _trace=0;                            THINK C PROFILER */

  if (inputsuccessful) {
    SetWriteFile(&writing,outputfile,'o',"output");
    if (VerifyInput){
      SetWriteFile(&writing_ver,verfile,'v',"input verification");
      WriteSolvedProblem(writing_ver);
      fclose(writing_ver);
      if (DynamicWriteOn) printf("closing the file %s\n",verfile);
     }
    if (DynamicWriteOn) {
      WriteRunningMode(stdout);
      WriteRunningMode(writing);
    }
    switch (Conversion) {
    case ExtToIne: case IneToExt: /* vertex/facets enumeration is chosen */
      if (RowDecomposition) DDRowDecomposition();
      else DDEnumerate();
      break;
    
    case LPmax:  case LPmin:      /* LP is chosen */
      LPMain();
      
      break;

    case Projection:              /* preprojection is chosen */
      PreProjection();
      break;

    case InteriorFind:      /* Interior point search is chosen */
      InteriorFindMain();
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
