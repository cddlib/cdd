/* cddarith.c:  Floating Point Arithmetic Procedures for cdd.c
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.51c, March 15, 1994 
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
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
/*  algebraic option is not efficint in most cases and deleted from Version 051 */
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
  if (strncmp(line, "initbasis_at_bottom", 19)==0) {
    InitBasisAtBottom = TRUE;
    return;
  }
  if (strncmp(line, "debug", 5)==0) {
    debug = TRUE;
    return;
  }
  if (strncmp(line, "partial_enum", 12)==0 && !PartialEnumeration) {
    fscanf(reading,"%ld", &msize);
    for (j=1;j<=msize;j++) {
      fscanf(reading,"%ld",&var);
       set_addelem(MarkedSet,var);
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Partial Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    PartialEnumeration=TRUE;
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
}

void AmatrixInput(boolean *successful)
{
  long i,j;
  double value;
  long value1,value2;
  boolean found=FALSE,decided=FALSE, fileopened;
  char command[wordlenmax], numbtype[wordlenmax], stemp[wordlenmax];

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
  
  PartialEnumeration = FALSE;
  set_initialize(&MarkedSet, minput+1);

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
    if (debug) putchar('\n');
  }  /*of i*/
  
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
    } else{
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


void StoreRay(double *p, RayRecord *RR, boolean *feasible)
{
  rowrange i,k,fii=mm+1;
  colrange j;
  double temp;

  *feasible = TRUE;
  set_initialize(&(RR->ZeroSet),mm);
  RR->ARay = 0.0;
  for (j = 0; j < nn; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= mm; i++) {
    k=OrderVector[i];
    temp = AValue(p, k);
    if (fabs(temp) < zero)
      set_addelem(RR->ZeroSet, k);
    if (temp < -zero){
      *feasible = FALSE;
      if (fii>mm) fii=i;  /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex=fii;
  if (debug) printf("store ray with fii= %ld\n", fii);
}


void AddRay(double *p)
{  
  boolean feasible;
  double x;

  if (FirstRay == NULL) {
    FirstRay = (struct RayRecord *) malloc(sizeof *FirstRay);
    FirstRay->Ray = (double *) calloc(nn, sizeof x);
    if (debug)
      printf("Create the first ray pointer\n");
    LastRay = FirstRay;
    ArtificialRay->Next = FirstRay;
  } else {
    LastRay->Next = (struct RayRecord *) malloc(sizeof *FirstRay);
    LastRay->Next->Ray = (double *) calloc(nn, sizeof x);
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
  double x;

  if (ArtificialRay != NULL) {
    printf("Warning !!!  FirstRay in not nil.  Illegal Call\n");
    return;
  }
  ArtificialRay = (struct RayRecord *) malloc(sizeof *ArtificialRay);
  ArtificialRay->Ray = (double *) calloc(nn, sizeof x);
  if (debug)
    printf("Create the artificial ray pointer\n");
  for (j = 0; j < nn; j++)
    zerovector[j] = 0.0;
  StoreRay(zerovector, ArtificialRay, &feasible);
  ArtificialRay->Next = NULL;
}

void ConditionalAddEdge(RayRecord *Ray1, RayRecord *Ray2, RayRecord *ValidFirstRay)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  boolean adjacent,lastchance;
  RayRecord *TempRay,*Rmin,*Rmax;
  AdjacencyRecord *NewEdge;
  boolean localdebug=FALSE;
  long *ZSmin, *ZSmax;
  
  fii1=Ray1->FirstInfeasIndex;
  fii2=Ray2->FirstInfeasIndex;
  if (fii1<fii2){
    fmin=fii1; fmax=fii2;
    Rmin=Ray1;
    Rmax=Ray2;
  }
  else{
    fmin=fii2; fmax=fii1;
    Rmin=Ray2;
    Rmax=Ray1;
  }
  ZSmin = Rmin->ZeroSet;
  ZSmax = Rmax->ZeroSet;
  if (localdebug) printf("ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
    fmin,OrderVector[fmin], fmax);
  if (fmin==fmax){
    if (localdebug) printf("ConditionalAddEdge: equal FII value-> No edge added\n");
  }
  else if (set_member(OrderVector[fmin],ZSmax)){
    if (localdebug) printf("ConditionalAddEdge: No strong separation -> No edge added\n");
  }
  else {  /* the pair will be separated at the iteration fmin */
    lastchance=TRUE;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    count_int++;
    if (localdebug){
      printf("Face: ");
      for (it=1; it<=mm; it++) {
        it_row=OrderVector[it];
        if (set_member(it_row, Face1)) printf("%ld ",it_row);
      }
      printf("\n");
    }
    for (it=Iteration+1; it < fmin && lastchance; it++){
      it_row=OrderVector[it];
      if (set_member(it_row, Face1)){
        lastchance=FALSE;
        count_int_bad++;
        if (localdebug){
          printf("There will be another chance iteration %ld (row %ld) to store the pair\n", it, it_row);
        }
      }
    }
    if (lastchance){
      adjacent = TRUE;
      count_int_good++;
      /* adjacent checking */
      set_int(Face, Face1, AddedHyperplanes);
      if (localdebug){
        printf("Check adjacency\n");
        printf("AddedHyperplanes: "); set_write(AddedHyperplanes);
        printf("Face: ");
        for (it=1; it<=mm; it++) {
          it_row=OrderVector[it];
          if (set_member(it_row, Face)) printf("%ld ",it_row);
        }
        printf("\n");
      }
      if (Cardinality(Face)< nn - 2) {
        adjacent = FALSE;
      }
      else if (NondegAssumed) {
    	adjacent = TRUE;
      }
      else{
        TempRay = ValidFirstRay;  /* the first ray for adjacency checking */
        while (TempRay != NULL && adjacent) {
          if (TempRay != Ray1 && TempRay != Ray2) {
            set_int(Face1, TempRay->ZeroSet, AddedHyperplanes);
            if (set_subset(Face, Face1)) {
              if (localdebug) set_write(Face1);
              adjacent = FALSE;
            }
          }
          TempRay = TempRay->Next;
        }
      }
      if (adjacent){
        if (localdebug) printf("The pair is adjacent and the pair must be stored for iteration %ld (row%ld)\n",
          fmin, OrderVector[fmin]);
        NewEdge=(struct AdjacencyRecord *) malloc(sizeof *NewEdge);
        NewEdge->Ray1=Rmax;  /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2=Rmin;  /* save the one deleted in iteration fmin in the second */
        NewEdge->Next=NULL;
        EdgeCount++;
        if (Edges[fmin]==NULL){
          Edges[fmin]=NewEdge;
          if (localdebug) printf("Create a new edge list of %ld\n", fmin);
        }else{
          NewEdge->Next=Edges[fmin];
          Edges[fmin]=NewEdge;
        }
      }
    }
  }
}

void CreateInitialEdges(void)
{
  RayRecord *Ptr1, *Ptr2;
  rowrange fii1,fii2;
  long count=0;
  boolean localdebug=FALSE;
  
  if (FirstRay ==NULL || LastRay==NULL){
    printf("Error found: CreateInitialEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=FirstRay;
  while(Ptr1!=LastRay && Ptr1!=NULL){
    fii1=Ptr1->FirstInfeasIndex;
    Ptr2=Ptr1->Next;
    while(Ptr2!=NULL){
      fii2=Ptr2->FirstInfeasIndex;
      count++;
      if (localdebug) printf("CreateInitialEdges: edge %ld \n",count);
      if (fii1!=fii2) ConditionalAddEdge(Ptr1,Ptr2,FirstRay);
      Ptr2=Ptr2->Next;
    }
    Ptr1=Ptr1->Next;
  }
_L99:;  
}


void UpdateEdges(RayRecord *RRbegin, RayRecord *RRend)
/* This procedure must be called after the ray list is sorted
   by EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  RayRecord *Ptr1, *Ptr2begin, *Ptr2;
  rowrange fii1;
  boolean ptr2found,quit,localdebug=FALSE;
  long count=0,pos1,pos2;
  float workleft,prevworkleft=100,totalpairs;

  totalpairs=(ZeroRayCount-1.0)*(ZeroRayCount-2.0)+1.0;
  Ptr2begin = NULL; 
  if (RRbegin ==NULL || RRend==NULL){
    if (1) printf("Warning: UpdateEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=RRbegin;
  pos1=1;
  do{
    ptr2found=FALSE;
    quit=FALSE;
    fii1=Ptr1->FirstInfeasIndex;
    pos2=2;
    for (Ptr2=Ptr1->Next; !ptr2found && !quit; Ptr2=Ptr2->Next, pos2++){
      if  (Ptr2->FirstInfeasIndex > fii1){
        Ptr2begin=Ptr2;
        ptr2found=TRUE;
      }
      else if (Ptr2==RRend) quit=TRUE;
    }
    if (ptr2found){
      quit=FALSE;
      for (Ptr2=Ptr2begin; !quit ; Ptr2=Ptr2->Next){
        count++;
        if (localdebug) printf("UpdateEdges: edge %ld \n",count);
        ConditionalAddEdge(Ptr1,Ptr2,RRbegin);
        if (Ptr2==RRend || Ptr2->Next==NULL) quit=TRUE;
      }
    }
    Ptr1=Ptr1->Next;
    pos1++;
    workleft = 100 * (ZeroRayCount-pos1) * (ZeroRayCount - pos1-1) / totalpairs;
    if (ZeroRayCount>=200 && DynamicWriteOn && pos1%10==0 && prevworkleft-workleft>=10 ) {
      printf("*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
	     Iteration, mm, pos1, ZeroRayCount, workleft);
      fprintf(writing,
	  "*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
	  Iteration, mm, pos1, ZeroRayCount, workleft);
	  prevworkleft=workleft;
	  fflush(writing);
	  if (writing_icd != NULL) fflush(writing_icd);
    }    
  }while(Ptr1!=RRend && Ptr1!=NULL);
_L99:;  
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


void ZeroIndexSet(double *x, rowset ZS)
{
  rowrange i;
  double temp;

  set_emptyset(ZS);
  for (i = 1; i <= mm; i++) {
    temp = AValue(x, i);
    if (fabs(temp) < zero)
      set_addelem(ZS, i);
  }
}

void CopyBmatrix(Bmatrix T, Bmatrix TCOPY)
{
  colrange j;

  for (j=0; j < nn; j++) {
    TCOPY[j] = T[j];
  }
}


void SelectPivot1(Amatrix X, HyperplaneOrderType roworder,
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
  rowset rowexcluded;

  stop = FALSE;
  set_initialize(&rowexcluded,mm);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    SelectNextHyperplane(roworder, rowexcluded, &rtemp, &RecomputeRowOrder);
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
  set_free(rowexcluded);
}

double TableauEntry(Amatrix X, Bmatrix T,
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

void WriteTableau(FILE *f, Amatrix X, Bmatrix T,
  InequalityType ineq)
/* Write the tableau  X.T   */
{
  colrange j;
  rowrange i;
  
  fprintf(f, "begin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",mm, nn+1);
  else
    fprintf(f, "  %ld   %ld    real\n",mm, nn);
  for (i=1; i<= mm; i++) {
    if (ineq==ZeroRHS)
      WriteReal(f, 0);  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j<= nn; j++) {
      fprintf(f," %5.2f",TableauEntry(X,T,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
}

void SelectPivot2(Amatrix X, Bmatrix T,
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
  rowrange i,rtemp;
  rowset rowexcluded;
  double Xtemp;

  stop = FALSE;
  set_initialize(&rowexcluded,mm);
  set_copy(rowexcluded,NopivotRow);
  if (debug) {
    printf("select pivot2: rowexcluded=");
    set_fwrite(stdout,rowexcluded);
  }
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    i=1;rtemp=0;
    while (i<=mm && rtemp==0) {  /* MarkedSet vars have highest priorities */
      if (set_member(i,MarkedSet) && !set_member(i,rowexcluded)){
        if (debug) printf("marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) SelectNextHyperplane(roworder, rowexcluded, &rtemp, &RecomputeRowOrder);
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
  set_free(rowexcluded);
}

void GausianColumnPivot1(Amatrix X, rowrange r, colrange s,
				rowrange rowmax)
/* Make a column pivot operation in Amatrix X on position (r,s)  */
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


void GausianColumnPivot2(Amatrix X, Bmatrix T,
				rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  long j, j1;
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


void InitializeBmatrix(Bmatrix T)
{
  colrange j;
  double x;

  for (j = 0; j < nn; j++) {
    T[j]=(double *)calloc(nn, sizeof x);
  }
}

void free_Bmatrix(Bmatrix T)
{
  colrange j;

  for (j = 0; j < nn; j++) {
    free(T[j]);
  }
}

void SetToIdentity(Bmatrix T)
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

void ReduceAA(long *ChosenRow, long *ChosenCol)
/* Set the matrix AA to be the submatrix of AA with chosen 
   rows & columns, and change mm and nn accordingly 
*/
{
  long i,j,inew,jnew,mnew=0,nnew=0;
  Amatrix Acopy;
  double x;

  mnew=set_card(ChosenRow);
  nnew=set_card(ChosenCol);
  for (i=0; i<mnew; i++){
    Acopy[i]=(double *)calloc(nnew,sizeof x);
  }
  inew=0;
  for (i=1; i <= mm; i++) {
    if (set_member(i,ChosenRow)) {
      inew++;
      jnew=0;
      for (j=1; j <= nn; j++) {
        if (set_member(j, ChosenCol)){
          jnew++;
          Acopy[inew-1][jnew-1]=AA[i-1][j-1];
          if (debug) WriteReal(stdout, AA[i-1][j-1]);
        }
      }
      if (debug) fprintf(stdout,"\n");
    }
  }
  for (i=1;i<=mm;i++) {
    if (i<=mnew) 
      set_addelem(ChosenRow,i);
    else
      set_delelem(ChosenRow,i);
  }
  for (j=1;j<=nn;j++) {
    if (j<=nnew) 
      set_addelem(ChosenCol,j);
    else
      set_delelem(ChosenCol,j);
  }
  if (debug) {
    fprintf(stdout, "new row indices:");set_write(ChosenRow);
    fprintf(stdout, "new col indices:");set_write(ChosenCol);
  }
  for (i=0;i<mm;i++){
    free(AA[i]);
  }
  for (i=0;i<mnew;i++){
    AA[i]=Acopy[i];
  }
  mm=mnew;  nn=nnew;
}

void DualizeAA(Bmatrix T)
/* Set the matrix AA to be the transpose of the matrix [-AA.T | I],
   and change mm and nn acordingly 
*/
{
  long i,j,mnew,nnew;
  Amatrix Acopy;
  double x;

  mnew=nn+mm;
  nnew=mm;
  for (i=0; i<mm; i++){
    Acopy[i]=(double *)calloc(nn,sizeof x);
  }
  if (mnew > MMAX) {
    printf("MMAX  is too small for ray computation. MMAX must be >= %ld.\n",mnew);
    Error = DimensionTooLarge;
    goto _L99;
  }
  if (nnew > NMAX) {
    printf("NMAX  is too small for ray computation. NMAX must be >= %ld.\n",nnew);
    Error = DimensionTooLarge;
    goto _L99;
  }
  for (i=1;i<=mm;i++){
    for (j=1;j<=nn;j++){
      Acopy[i-1][j-1]=TableauEntry(AA,T,i,j);
    }
  }
  for (i=0;i<mm;i++){
    free(AA[i]);
  }
  for (i=0; i<mnew; i++){
    AA[i]=(double *)calloc(nnew,sizeof x);
  }
  for (i=1;i<=nn;i++){
    for (j=1;j<=mm;j++){
      AA[i-1][j-1]=-Acopy[j-1][i-1];
    }
  }
  for (i=1;i<=mm;i++){
    for (j=1;j<=mm;j++){
      if (i==j) AA[nn+i-1][j-1]=1;
      else AA[nn+i-1][j-1]=0;
    }
  }
  for (i=0; i<mm; i++){
    free(Acopy[i]);
  }
  mm=mnew;  nn=nnew;
  _L99:;
}

void EnlargeAAforInteriorFinding(void)
/* Add an extra column with all minus ones to the matrix AA, 
   add an objective row with (0,...,0,1), and 
   rows & columns, and change mm and nn accordingly 
*/
{
  long i,j,mnew=0,nnew=0;
  Amatrix Acopy;
  double x;

  mnew=mm+1;
  nnew=nn+1;
  for (i=0; i<mnew; i++){
    Acopy[i]=(double *)calloc(nnew,sizeof x);
  }
  for (i=1; i <= mm; i++) {
    for (j=1; j <= nn; j++) {
      Acopy[i-1][j-1]=AA[i-1][j-1];
      if (debug) WriteReal(stdout, AA[i-1][j-1]);
    }
    if (debug) fprintf(stdout,"\n");
  }
  for (i=1;i<=mm;i++) {
    Acopy[i-1][nn]=-1.0;  /* new column with all minus one's */
  }
  for (j=1;j<=nn;j++) {
    Acopy[mm][j-1]=0.0;  /* new row with (0,...,0,1) */
  }
  Acopy[mm][nn]=1.0;  /* new row with (0,...,0,1) */
  for (i=0;i<mm;i++){
    free(AA[i]);
  }
  for (i=0;i<mnew;i++){
    AA[i]=Acopy[i];
  }
  mm=mnew;  nn=nnew;
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
      fprintf(f, "%5.2f ", T[j1][j2]);
    }  /*of j2*/
    putc('\n', f);
  }  /*of j1*/
  putc('\n', f);
}

void ComputeRank(Amatrix A1, long *TargetRows, long *rank)
/* Compute the rank of the submatrix of a Amatrix A1 indexed by TargetRows.
   This procedure does not change the matrix A1.
 */
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset NoPivotRow;
  colset ColSelected;
  Bmatrix Btemp;   /* dual basis inverse */
  
  *rank = 0;
  stop = FALSE;
  set_initialize(&NoPivotRow, mm);
  set_compl(NoPivotRow,TargetRows);
  set_initialize(&ColSelected, nn);
  InitializeBmatrix(Btemp);
  SetToIdentity(Btemp);
  if (debug) WriteBmatrix(stdout,Btemp);
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, Btemp, MinIndex, mm, NoPivotRow, ColSelected, &r, &s, &chosen);
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
  set_free(NoPivotRow);
  set_free(ColSelected);
  free_Bmatrix(Btemp);
}

void ComputeBInverse(Amatrix A1, long lastrow,
       Bmatrix InvA1, long *rank)
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset RowSelected;
  colset ColSelected;

  *rank = 0;
  stop = FALSE;
  SetToIdentity(InvA1);
  set_initialize(&RowSelected, mm);
  set_initialize(&ColSelected, nn);
  do {
    SelectPivot2(A1, InvA1, MinIndex, lastrow, RowSelected, ColSelected, &r, &s, &chosen);
    if (chosen) {
      (*rank)++;
      if (debug)
        printf("%3ldth pivot on%3ld, %3ld\n", *rank, r, s);
      GausianColumnPivot2(A1, InvA1, r, s);
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
    } else
      stop = TRUE;
  } while (!stop);
  set_free(RowSelected);
  set_free(ColSelected);
}


void FindBasis(Amatrix A1, 
              HyperplaneOrderType roword, 
              rowset RowSelected,colset ColInd,
              Bmatrix BasisInverse, long *rank)
{
  boolean stop, chosen;
  colset ColSelected;
  rowrange r;
  colrange j,s;

  *rank = 0;
  stop = FALSE;
  for (j=0;j<=nn;j++) ColInd[j]=0;
  set_emptyset(RowSelected);
  set_initialize(&ColSelected, nn);
  SetToIdentity(BasisInverse);
  if (debug) WriteBmatrix(stdout,BasisInverse);
  if (DynamicWriteOn && !debug){
    printf("*Initial set of rows:");
    fprintf(writing,"*Initial set of rows:");
  }
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, BasisInverse, roword, mm, RowSelected, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        ColInd[s]=r;    /* ColInd[s] stores the corr. row index */
        (*rank)++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
        if (debug) {
          WriteBmatrix(stdout,BasisInverse);
          WriteTableau(stdout,A1,BasisInverse,NonzeroRHS),
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
	    if (DynamicWriteOn && !debug){
	      printf(" %ld",  r);
	      fprintf(writing, " %ld",  r);
	    }
      } else {
        stop=TRUE;
      }
      if (*rank==nn) stop = TRUE;
  } while (!stop);
  if (DynamicWriteOn && !debug){
    printf("\n");
    fprintf(writing,"\n");
  }
  set_free(ColSelected);
}


void SelectCrissCrossPivot(Amatrix X, Bmatrix T,
            long bflag[], rowrange objrow, colrange rhscol,
            rowrange *r, colrange *s,
			boolean *selected, LPStatusType *lps)
{
  boolean colselected=FALSE, rowselected=FALSE, stop=FALSE;
  rowrange i;
  double val;
  
  *selected=FALSE;
  *lps=LPSundecided;
  while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=mm; i++) {
      if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=TableauEntry(X,T,i,rhscol);
        if (val < -zero) {
          rowselected=TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=TableauEntry(X,T,objrow,bflag[i]);
        if (val > zero) {
          colselected=TRUE;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=mm; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          val=TableauEntry(X,T,*r,bflag[i]);
          if (val > zero) {
            colselected=TRUE;
            *s=bflag[i];
            *selected=TRUE;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=mm; i++) {
        if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=TableauEntry(X,T,i,*s);
          if (val < -zero) {
            rowselected=TRUE;
            *r=i;
            *selected=TRUE;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=DualInconsistent;
    }
    else if (!colselected) {
      *lps=Inconsistent;
    }
  }
}

void CrissCrossMinimize(Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
{
   colrange j;
   
   for (j=1; j<=nn; j++)
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   CrissCrossMaximize(A1,BasisInverse, OBJrow, RHScol, 
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter);
   *optvalue=-*optvalue;
   for (j=1; j<=nn; j++)
     dsol[j-1]=-dsol[j-1];
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
}

void CrissCrossMaximize(Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  boolean stop, chosen;
  long rank;
  rowrange i,r,entering,leaving;
  colrange j,s;
  colset ColSelected;
  rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0;

  if (BasisFlag==NULL || mlast!=mm){
     if (mlast!=mm) free(BasisFlag);   /* called previously with different mm */
     BasisFlag=(long *) calloc(mm+1, sizeof *BasisFlag);  
     /* initialize only for the first time or when a larger space is needed */
     mlast=mm;
  }
  *re=0; *se=0; *iter=0;
  rank = 0;
  stop = FALSE;
  InitializeBmatrix(InitialRays);
  set_initialize(&Cobasis,mm);
  set_initialize(&Basis,mm);
  set_initialize(&RowSelected, mm);
  set_initialize(&ColSelected, nn);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  for (i=0; i<=mm; i++) BasisFlag[i]=0;
  for (j=0; j<=NMAX; j++) NBIndex[j]=0;
  for (i=1; i<=mm; i++) {
    set_addelem(Basis,i);
    BasisFlag[i]=-1;    /* basic variable has index -1 */
  }
  BasisFlag[OBJrow]= 0;
  SetToIdentity(BasisInverse);
  if (debug) WriteBmatrix(stdout,BasisInverse);
  do {   /* Find a LP basis */
      SelectPivot2(A1, BasisInverse, MinIndex
      , mm, RowSelected, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        set_addelem(Cobasis, r);
        set_delelem(Basis,r);
        if (debug) {
          fprintf(stdout, "CC: find initial basis: nonbasis = ");   
          WriteSetElements(stdout,Cobasis); fprintf(stdout,"\n");
          fprintf(stdout, "CC: find initial basis: basis = ");
          WriteSetElements(stdout,Basis);fprintf(stdout,"\n");
        }
        BasisFlag[r]=s;   /* the nonbasic variable r corresponds to column s */
        NBIndex[s]=r;     /* the nonbasic variable on s column is r */
        if (debug) fprintf(stdout,"nonbasic variable %ld has index %ld\n", r, BasisFlag[r]);
        rank++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
        if (debug) {
          WriteBmatrix(stdout,BasisInverse);
          WriteTableau(stdout,A1,BasisInverse,Inequality),
	      fprintf(stdout, "%3ldth row added to the initial set (%ldth elem)\n",  r, rank);
	    }
      } else {
        stop=TRUE;
      }
      if (rank==nn-1) stop = TRUE;
  } while (!stop);
  stop=FALSE;
  do {   /* Criss-Cross Method */
    SelectCrissCrossPivot(A1, BasisInverse, BasisFlag,
       OBJrow, RHScol, &r, &s, &chosen, LPS);
    if (debug && chosen) fprintf(stdout, "Procedure Criss-Cross: pivot on (r,s) =(%ld, %ld).\n", r, s);
    if (chosen) {
      entering=NBIndex[s];
      leaving=r;
      set_addelem(Cobasis, leaving);
      set_delelem(Cobasis, entering);
      set_delelem(Basis,leaving);
      set_addelem(Basis,entering);
      BasisFlag[leaving]=s;
      BasisFlag[entering]=-1;
      NBIndex[s]=leaving;
      if (debug) {
        fprintf(stdout, "nonbasis = "); set_write(Cobasis);
        fprintf(stdout, "basis = "); set_write(Basis);
        fprintf(stdout, "new nonbasic variable %ld has index %ld\n", leaving, BasisFlag[leaving]);
      }
      GausianColumnPivot2(A1,BasisInverse, r, s);
      (*iter)++;
      if (debug) {
        WriteBmatrix(stdout,BasisInverse);
        WriteTableau(stdout,A1,BasisInverse,Inequality),
	    fprintf(stdout, "%3ldth row added to the initial set (%ldth elem)\n",  r, rank);
	  }
    } else {
      switch (*LPS){
        case Inconsistent: *re=r;
        case DualInconsistent: *se=s;
        default: break;
      }
      stop=TRUE;
    }
  } while(!stop);
  switch (*LPS){
  case Optimal:
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,OBJrow,j);
      *optvalue=TableauEntry(A1,BasisInverse,OBJrow,RHScol);
      if (debug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    break;
  case Inconsistent:
    if (debug) printf("CrissCrossSolve: LP is inconsistent.\n");
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,*re,j);
      if (debug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    break;
  case DualInconsistent:
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][*se-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,OBJrow,j);
      if (debug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    if (debug) printf("CrissCrossSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }
  set_free(ColSelected);
  set_free(RowSelected);
  set_free(Basis);
  set_free(Cobasis);
}

void WriteLPResult(FILE *f, LPStatusType LPS, double optval,
   Arow sol, Arow dsol, colindex NBIndex, rowrange re, colrange se,
   long iter)
{
  long j;

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
}


void FindInitialRays(rowset InitHyperplanes,
			    Bmatrix InitRays, boolean *found)
{
  Bmatrix BInverse;
  rowset CandidateRows;
  colindex PivRow;
  long i, rank;
  HyperplaneOrderType roworder;

  *found = FALSE;
  
  set_initialize(&CandidateRows, mm);
  switch (HyperplaneOrder) {
  case MaxIndex:
    roworder = MaxIndex;
    break;

  case MinIndex:
    roworder = MinIndex;
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
  if (InitBasisAtBottom==TRUE) roworder=MaxIndex;
  for (i = 1; i <= mm; i++)
    set_addelem(CandidateRows, i);      /*all rows are candidates for initial cone*/
  if (DynamicWriteOn)
    printf("*Computing an initial set of rays\n");
  InitializeBmatrix(BInverse);
  FindBasis(AA, roworder, InitHyperplanes, PivRow, BInverse, &rank);
  if (debug) {
    printf("FindInitialBasis: InitHyperplanes=");
    set_fwrite(stdout,InitHyperplanes);
  }
  if (debug) printf("nn = %ld, rank = %ld\n",nn,rank);
  if (rank < nn) {
    Error = LowColumnRank;
    return;
  }
  if (!set_subset(MarkedSet,InitHyperplanes)) {
    Error = DependentMarkedSet;
    return;
  }
  *found = TRUE;
  if (debug) {
    WriteBmatrix(stdout, BInverse);
  }
  /* free_Bmatrix(InitRays);  */
  CopyBmatrix(BInverse,InitRays);
  if (debug) WriteBmatrix(stdout, InitRays);
  set_free(CandidateRows);
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
  if (debug){
    WriteRayRecord(stdout, Ptr1);
    WriteRayRecord(stdout, Ptr2);  
    printf("create a new ray by eliminating %ld:\n",ii);
    WriteRayRecord(stdout, LastRay);
  }
}

void EvaluateARay1(rowrange i)
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

void EvaluateARay2(rowrange i)
/* Evaluate the ith component of the vector  A x RD.Ray 
   and rearrange the linked list so that
   the infeasible rays with respect to  i  will be
   placed consecutively from First. Also for all feasible rays,
   "positive" rays and "zero" rays will be placed consecutively.
 */
{
  colrange j;
  double temp;
  RayRecord *Ptr, *NextPtr;
  boolean zerofound=FALSE,negfound=FALSE,posfound=FALSE;
  boolean localdebug=FALSE;

  PosHead=NULL;ZeroHead=NULL;NegHead=NULL;
  PosLast=NULL;ZeroLast=NULL;NegLast=NULL;
  Ptr = FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=NULL;     /* then clear the Next pointer */
    temp = 0.0;
    for (j = 0; j < nn; j++)
      temp += AA[i - 1][j] * Ptr->Ray[j];
    Ptr->ARay = temp;
    if ( temp < -zero) {
      if (!negfound){
        negfound=TRUE;
        NegHead=Ptr;
        NegLast=Ptr;
      }
      else{
        Ptr->Next=NegHead;
        NegHead=Ptr;
      }
    }
    else if (temp > zero){
      if (!posfound){
        posfound=TRUE;
        PosHead=Ptr;
        PosLast=Ptr;
      }
      else{  
        Ptr->Next=PosHead;
        PosHead=Ptr;
       }
    }
    else {
      if (!zerofound){
        zerofound=TRUE;
        ZeroHead=Ptr;
        ZeroLast=Ptr;
      }
      else{
        Ptr->Next=ZeroHead;
        ZeroHead=Ptr;
      }
    }
    Ptr=NextPtr;
  }
  /* joining three neg, pos and zero lists */
  if (negfound){                 /* -list nonempty */
    FirstRay=NegHead;
    if (posfound){               /* -list & +list nonempty */
      NegLast->Next=PosHead;
      if (zerofound){            /* -list, +list, 0list all nonempty */
        PosLast->Next=ZeroHead;
        LastRay=ZeroLast;
      } 
      else{                      /* -list, +list nonempty but  0list empty */
        LastRay=PosLast;      
      }
    }
    else{                        /* -list nonempty & +list empty */
      if (zerofound){            /* -list,0list nonempty & +list empty */
        NegLast->Next=ZeroHead;
        LastRay=ZeroLast;
      } 
      else {                      /* -list nonempty & +list,0list empty */
        LastRay=NegLast;
      }
    }
  }
  else if (posfound){            /* -list empty & +list nonempty */
    FirstRay=PosHead;
    if (zerofound){              /* -list empty & +list,0list nonempty */
      PosLast->Next=ZeroHead;
      LastRay=ZeroLast;
    } 
    else{                        /* -list,0list empty & +list nonempty */
      LastRay=PosLast;
    }
  }
  else{                          /* -list,+list empty & 0list nonempty */
    FirstRay=ZeroHead;
    LastRay=ZeroLast;
  }
  ArtificialRay->Next=FirstRay;
  LastRay->Next=NULL;
}

void DeleteNegativeRays(void)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the RayRecord list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  rowrange fii,fiitest;
  double temp;
  RayRecord *Ptr, *PrevPtr, *NextPtr, *ZeroPtr1, *ZeroPtr0;
  boolean found, completed, zerofound=FALSE,negfound=FALSE,posfound=FALSE;
  boolean localdebug=FALSE;
  
  PosHead=NULL;ZeroHead=NULL;NegHead=NULL;
  PosLast=NULL;ZeroLast=NULL;NegLast=NULL;

  /* Delete the infeasible rays  */
  PrevPtr= ArtificialRay;
  Ptr = FirstRay;
  ZeroRayCount=0;
  if (PrevPtr->Next != Ptr) 
    printf("Error at DeleteNegativeRays: ArtificialRay does not point the FirstRay.\n");
  completed=FALSE;
  while (Ptr != NULL && !completed){
    if ( (Ptr->ARay) < -zero ){
      Eliminate(&PrevPtr);
      Ptr=PrevPtr->Next;
    }
    else{
      completed=TRUE;
    }
  }
  
  /* Sort the zero rays */
  Ptr = FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    temp = Ptr->ARay;
    if (localdebug) printf("Ptr->ARay : %5.3f \n", temp);
    if ( temp < -zero) {
      if (!negfound){
        printf("Error: An infeasible ray found after their removal\n");
        negfound=TRUE;
      }
    }
    else if (temp > zero){
      if (!posfound){
        posfound=TRUE;
        PosHead=Ptr;
        PosLast=Ptr;
      }
      else{  
        PosLast=Ptr;
       }
    }
    else {
      ZeroRayCount++;
      if (!zerofound){
        zerofound=TRUE;
        ZeroHead=Ptr;
        ZeroLast=Ptr;
        ZeroLast->Next=NULL;
      }
      else{/* Find a right position to store the record sorted w.r.t. FirstInfeasIndex */
        fii=Ptr->FirstInfeasIndex; 
        found=FALSE;
        ZeroPtr1=NULL;
        for (ZeroPtr0=ZeroHead; !found && ZeroPtr0!=NULL ; ZeroPtr0=ZeroPtr0->Next){
          fiitest=ZeroPtr0->FirstInfeasIndex;
          if (fiitest >= fii){
            found=TRUE;
          }
          else ZeroPtr1=ZeroPtr0;
        }
        /* printf("insert position found \n %d  index %ld\n",found, fiitest); */
        if (!found){           /* the new record must be stored at the end of list */
          ZeroLast->Next=Ptr;
          ZeroLast=Ptr;
          ZeroLast->Next=NULL;
        }
        else{
          if (ZeroPtr1==NULL){ /* store the new one at the head, and update the head ptr */
            /* printf("Insert at the head\n"); */
            Ptr->Next=ZeroHead;
            ZeroHead=Ptr;
          }
          else{                /* store the new one inbetween ZeroPtr1 and 0 */
            /* printf("Insert inbetween\n");  */
            Ptr->Next=ZeroPtr1->Next;
            ZeroPtr1->Next=Ptr;
          }
        }
        /*
        Ptr->Next=ZeroHead;
        ZeroHead=Ptr;
        */
      }
    }
    Ptr=NextPtr;
  }
  /* joining the pos and zero lists */
  if (posfound){            /* -list empty & +list nonempty */
    FirstRay=PosHead;
    if (zerofound){              /* +list,0list nonempty */
      PosLast->Next=ZeroHead;
      LastRay=ZeroLast;
    } 
    else{                        /* 0list empty & +list nonempty */
      LastRay=PosLast;
    }
  }
  else{                          /* +list empty & 0list nonempty */
    FirstRay=ZeroHead;
    LastRay=ZeroLast;
  }
  ArtificialRay->Next=FirstRay;
  LastRay->Next=NULL;
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
  } while (!(determined) && (j <= nn));
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

void AddNewHyperplane1(rowrange hnew)
/* This procedure 1 must be used with PreorderedRun=FALSE 
   This procedure is the most elementary implementation of
   DD and can be used with any type of ordering, including
   dynamic ordering of rows, e.g. MaxCutoff, MinCutoff.
   The memory requirement is minimum because it does not
   store any adjacency among the rays.
*/
{
  RayRecord *RayPtr0, *RayPtr1, *RayPtr2, *RayPtr2s, *RayPtr3;
  long pos1, pos2;
  double prevprogress, progress, value1, value2;
  boolean adj, equal, completed;

  EvaluateARay1(hnew);  
   /*Check feasibility of rays w.r.t. hnew 
     and put all infeasible ones consecutively */

  RayPtr0 = ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  value1 = FirstRay->ARay;
  if (value1 > -zero) {
    if (RayCount==FeasibleRayCount) CompStatus=AllFound;
    goto _L99;        /* Sicne there is no hnew-infeasible ray and nothing to do */
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
  RayPtr2 = RayPtr2s;   /*2nd feasible ray to scan and compare with 1st*/
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

void AddNewHyperplane2(rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  RayRecord *RayPtr0,*RayPtr1, *RayPtr2;
  AdjacencyRecord *EdgePtr, *EdgePtr0;
  long pos1;
  rowrange fii1, fii2;
  boolean localdebug=FALSE;

  EvaluateARay2(hnew);
   /* Check feasibility of rays w.r.t. hnew 
      and sort them. ( -rays, +rays, 0rays)*/

  if (PosHead==NULL && ZeroHead==NULL) {
    FirstRay=NULL;
    ArtificialRay->Next=FirstRay;
    CompStatus=RegionEmpty;
    RayCount=0;
    VertexCount=0;
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }

  if (localdebug){
    pos1=0;
    printf("(pos, FirstInfeasIndex, A Ray)=\n");
    for (RayPtr0=FirstRay; RayPtr0!=NULL; RayPtr0=RayPtr0->Next){
      pos1++;
      printf("(%ld,%ld,",pos1,RayPtr0->FirstInfeasIndex);
      WriteReal(stdout,RayPtr0->ARay); 
      printf(") ");
   }
    printf("\n");
  }
  
  if (ZeroHead==NULL) ZeroHead=LastRay;

  EdgePtr=Edges[Iteration];
  while (EdgePtr!=NULL){
    RayPtr1=EdgePtr->Ray1;
    RayPtr2=EdgePtr->Ray2;
    fii1=RayPtr1->FirstInfeasIndex;   
    CreateNewRay(RayPtr1, RayPtr2, hnew);
    fii2=LastRay->FirstInfeasIndex;
    if (fii1 != fii2) ConditionalAddEdge(RayPtr1,LastRay,PosHead);
    EdgePtr0=EdgePtr;
    EdgePtr=EdgePtr->Next;
    free(EdgePtr0);
    EdgeCount--;
  }
  
  DeleteNegativeRays();
    
  set_addelem(AddedHyperplanes, hnew);

  if (Iteration<mm){
    if (ZeroHead!=NULL && ZeroHead!=LastRay){
      if (ZeroRayCount>=200 && DynamicWriteOn) printf("*New edges being scanned...\n");
      UpdateEdges(ZeroHead,LastRay);
    }
    if (localdebug) printf("*Edges currently stored = %ld\n", EdgeCount);
  }

  if (RayCount==FeasibleRayCount) CompStatus=AllFound;
_L99:;
}

/* end of cddarith.c */

