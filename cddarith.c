/* cddarith.c:  Floating Point Arithmetic Procedures for cdd.c
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.60, August 21, 1996
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


void SetInequalitySets(rowrange m_size, rowindex id)
{
  long i;
  static long mprev=0;
  boolean localdebug=FALSE;
  
  if (mprev!=m_size){
    if (localdebug) printf("SetInequalitySets: initializing inequality sets.\n");
    if (GroundSet!=NULL) set_free(GroundSet);
    if (EqualitySet!=NULL) set_free(EqualitySet);
    if (NonequalitySet!=NULL) set_free(NonequalitySet);
    set_initialize(&EqualitySet, m_size);
    set_initialize(&NonequalitySet, m_size);
    set_initialize(&GroundSet, m_size);
  }
  if (localdebug && mprev==m_size) printf("SetInequalitySets: Resetting inequality sets.\n");
  set_emptyset(GroundSet);
  set_emptyset(EqualitySet);
  set_emptyset(NonequalitySet);  
  for (i = 1; i <= m_size; i++){
    set_addelem(GroundSet, i);
    if (id[i]==1) set_addelem(EqualitySet,i);
    if (id[i]==-1) set_addelem(NonequalitySet,i);
  }
  mprev=m_size;
}


double AValue(colrange n_size, Amatrix A, double *p, rowrange i)
{
  /*return the ith component of the vector  A x p */
  colrange j;
  double temp;

  temp = 0.0;
  for (j = 0; j < n_size; j++)
    temp += A[i - 1][j] * p[j];
  return temp;
}

void StoreRay1(rowrange m_size, colrange n_size, Amatrix A, 
    double *p, RayRecord *RR, rowindex ordervec, boolean *feasible)
{  /* Original ray storing routine when RelaxedEnumeration is FALSE */
  rowrange i,k,fii=m_size+1;
  colrange j;
  double temp;
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  *feasible = TRUE;
  set_initialize(&(RR->ZeroSet),m_size);
  RR->ARay = 0.0;
  for (j = 0; j < n_size; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= m_size; i++) {
    k=ordervec[i];
    temp = AValue(n_size, A, p, k);
    if (fabs(temp) < zero)
      set_addelem(RR->ZeroSet, k);
    if (temp < -zero){
      *feasible = FALSE;
      if (fii>m_size) fii=i;  /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
  if (localdebug) 
    if (fii<= m_size)
      printf("StoreRay1:store ray with fii= %ld (row =%ld)\n", fii,ordervec[fii]);
    else
      printf("StoreRay1:store ray with fii= %ld (feasible=%d)\n", fii,*feasible);    
}

void StoreRay2(rowrange m_size, colrange n_size, Amatrix A, 
    double *p, RayRecord *RR, rowindex ordervec, 
    boolean *feasible, boolean *weaklyfeasible)
   /* Ray storing routine when RelaxedEnumeration is TRUE.
       weaklyfeasible is true iff it is feasible with
       the strict_inequality conditions deleted. */
{
  rowrange i,k,fii=m_size+1;
  colrange j;
  double temp;
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  *feasible = TRUE;
  *weaklyfeasible = TRUE;
  set_initialize(&(RR->ZeroSet),m_size);
  RR->ARay = 0.0;
  for (j = 0; j < n_size; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= m_size; i++) {
    k=ordervec[i];
    temp = AValue(n_size, A, p, k);
    if (fabs(temp) < zero){
      set_addelem(RR->ZeroSet, k);
      if (EqualityIndex[k]==-1) 
        *feasible=FALSE;  /* strict inequality required */
    }
    if (temp < -zero){
      *feasible = FALSE;
      if (fii>m_size && EqualityIndex[k]>=0) {
        fii=i;  /* the first violating inequality index */
        *weaklyfeasible=FALSE;
      }
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
  if (localdebug) {
    if (fii<= m_size)
      printf("StoreRay2:store ray with fii= %ld (row =%ld)\n", fii,ordervec[fii]);
    else
      printf("StoreRay2:store ray with fii= %ld (weaklyfeasible=%d)\n", fii,*weaklyfeasible);    
    if (*weaklyfeasible) WriteRayRecord(stdout, n_size, RR);
  }
}


void AddRay(rowrange m_size, colrange n_size, Amatrix A, double *p, rowindex ordervec)
{  
  boolean feasible, weaklyfeasible;
  double x;

  if (FirstRay == NULL) {
    FirstRay = (struct RayRecord *) malloc(sizeof *FirstRay);
    FirstRay->Ray = (double *) calloc(n_size, sizeof x);
    if (debug)
      printf("Create the first ray pointer\n");
    LastRay = FirstRay;
    ArtificialRay->Next = FirstRay;
  } else {
    LastRay->Next = (struct RayRecord *) malloc(sizeof *FirstRay);
    LastRay->Next->Ray = (double *) calloc(n_size, sizeof x);
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
  if (RelaxedEnumeration){
    StoreRay2(m_size, n_size, A, p, LastRay, ordervec, &feasible, &weaklyfeasible);
    if (weaklyfeasible) WeaklyFeasibleRayCount++;
  } else {
    StoreRay1(m_size, n_size, A, p, LastRay, ordervec, &feasible);
    if (feasible) WeaklyFeasibleRayCount++;
    /* weaklyfeasible is equiv. to feasible in this case. */
  }
  if (!feasible) return;
  else {
    FeasibleRayCount++;
    if (fabs(LastRay->Ray[0]) > zero && Inequality == NonzeroRHS)
      VertexCount++;
    if (DynamicRayWriteOn) {
      WriteRayRecord(writing, n_size, LastRay);
      WriteRayRecord(stdout, n_size, LastRay);
      if (IncidenceOutput==IncSet && writing_icd != NULL) {
     	WriteIncidence(writing_icd, m_size, n_size, LastRay);
      }
    }
  }
}

void AddArtificialRay(rowrange m_size, colrange n_size, Amatrix A, rowindex ordervec)
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
  ArtificialRay->Ray = (double *) calloc(n_size, sizeof x);
  if (debug)
    printf("Create the artificial ray pointer\n");
  for (j = 0; j < n_size; j++)
    zerovector[j] = 0.0;
  StoreRay1(m_size, n_size, A, zerovector, ArtificialRay, ordervec, &feasible);
  ArtificialRay->Next = NULL;
}

void ConditionalAddEdge(rowrange m_size, colrange n_size, 
    RayRecord *Ray1, RayRecord *Ray2, RayRecord *ValidFirstRay, 
    rowrange iter, rowindex ordervec)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  boolean adjacent,lastchance;
  RayRecord *TempRay,*Rmin,*Rmax;
  AdjacencyRecord *NewEdge;
  boolean localdebug=FALSE;
  rowset ZSmin, ZSmax;
  
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
  if (localdebug) {
    printf("ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
      fmin,ordervec[fmin], fmax);
  }
  if (fmin==fmax){
    if (localdebug) printf("ConditionalAddEdge: equal FII value-> No edge added\n");
  }
  else if (set_member(ordervec[fmin],ZSmax)){
    if (localdebug) printf("ConditionalAddEdge: No strong separation -> No edge added\n");
  }
  else {  /* the pair will be separated at the iteration fmin */
    lastchance=TRUE;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    count_int++;
    if (localdebug){
      printf("Face: ");
      for (it=1; it<=m_size; it++) {
        it_row=ordervec[it];
        if (set_member(it_row, Face1)) printf("%ld ",it_row);
      }
      printf("\n");
    }
    for (it=iter+1; it < fmin && lastchance; it++){
      it_row=ordervec[it];
      if (EqualityIndex[it_row]>=0 && set_member(it_row, Face1)){
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
        for (it=1; it<=m_size; it++) {
          it_row=ordervec[it];
          if (set_member(it_row, Face)) printf("%ld ",it_row);
        }
        printf("\n");
      }
      if (set_card(Face)< n_size - 2) {
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
          fmin, ordervec[fmin]);
        NewEdge=(struct AdjacencyRecord *) malloc(sizeof *NewEdge);
        NewEdge->Ray1=Rmax;  /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2=Rmin;  /* save the one deleted in iteration fmin in the second */
        NewEdge->Next=NULL;
        EdgeCount++; 
        TotalEdgeCount++;
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

void CreateInitialEdges(rowrange m_size, colrange n_size, rowindex ordervec)
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
      if (fii1!=fii2) 
        ConditionalAddEdge(m_size, n_size, Ptr1, Ptr2, FirstRay, n_size, ordervec);
      Ptr2=Ptr2->Next;
    }
    Ptr1=Ptr1->Next;
  }
_L99:;  
}


void UpdateEdges(rowrange m_size, colrange n_size, 
    RayRecord *RRbegin, RayRecord *RRend, rowrange iter, rowindex ordervec)
/* This procedure must be called after the ray list is sorted
   by EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  RayRecord *Ptr1, *Ptr2begin, *Ptr2;
  rowrange fii1;
  boolean ptr2found,quit,localdebug=FALSE;
  long count=0,pos1, pos2;
  float workleft,prevworkleft=110.0,totalpairs;

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
    for (Ptr2=Ptr1->Next; !ptr2found && !quit; Ptr2=Ptr2->Next,pos2++){
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
        ConditionalAddEdge(m_size, n_size, Ptr1,Ptr2,RRbegin, iter, ordervec);
        if (Ptr2==RRend || Ptr2->Next==NULL) quit=TRUE;
      }
    }
    Ptr1=Ptr1->Next;
    pos1++;
    workleft = 100.0 * (ZeroRayCount-pos1) * (ZeroRayCount - pos1-1.0) / totalpairs;
    if (ZeroRayCount>=500 && DynamicWriteOn && pos1%10 ==0 && prevworkleft-workleft>=10 ) {
      printf("*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
	     iter, m_size, pos1, ZeroRayCount, workleft);
      fprintf(writing,
	 "*Work of iteration %5ld(/%ld): %4ld/%4ld => %4.1f%% left\n",
	 iter, m_size, pos1, ZeroRayCount, workleft);
      prevworkleft=workleft;
      fflush(writing);
      if (writing_icd != NULL) fflush(writing_icd);
    }    
  }while(Ptr1!=RRend && Ptr1!=NULL);
_L99:;  
}

void FreeDDMemory(rowindex ordervec)
{
  RayRecord *Ptr, *PrevPtr;
  long count;
  boolean localdebug=FALSE;
  
  PrevPtr=ArtificialRay;
  count=0;
  for (Ptr=ArtificialRay->Next; Ptr!=NULL; Ptr=Ptr->Next){
    free(PrevPtr->Ray);
    free(PrevPtr->ZeroSet);
    free(PrevPtr);
    count++;
    PrevPtr=Ptr;
  };
  LastRay=NULL;
  FirstRay=NULL;
  ArtificialRay=NULL;
  if (localdebug) printf("%ld ray storage spaces freed\n",count);
  
  set_free(InitialHyperplanes);
  set_free(AddedHyperplanes);
  set_free(GroundSet);
  set_free(Face);
  set_free(Face1);
  free(ordervec);
  
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  WeaklyFeasibleRayCount = 0;
  VertexCount = 0;
}

void Normalize(colrange n_size, double *V)
{
  long j;
  double min, temp;

  min = 1.0e+20;
  for (j = 0; j < n_size; j++) {
    temp = fabs(V[j]);
    if (temp > zero && temp < min)
      min = temp;
  }
  for (j = 0; j < n_size; j++)
    V[j] /= min;
}


void ZeroIndexSet(rowrange m_size, colrange n_size, Amatrix A, double *x, rowset ZS)
{
  rowrange i;
  double temp;

  set_emptyset(ZS);
  for (i = 1; i <= m_size; i++) {
    temp = AValue(n_size, A, x, i);
    if (fabs(temp) < zero)
      set_addelem(ZS, i);
  }
}

void CopyBmatrix(colrange n_size, Bmatrix T, Bmatrix TCOPY)
{
  colrange j;

  for (j=0; j < n_size; j++) {
    TCOPY[j] = T[j];
  }
}


void SelectPivot1(rowrange m_size, colrange n_size, 
    Amatrix X, HyperplaneOrderType roworder, rowindex ordervec,
    rowrange rowmax, rowset NopivotRow,
    colset NopivotCol, rowrange *r, colrange *s,
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
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=m_size;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    SelectNextHyperplane(m_size, n_size, X, roworder, rowexcluded, 
      &rtemp, &RecomputeRowOrder, ordervec);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= n_size && !*selected) {
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

double TableauEntry(rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
				rowrange r, colrange s)
/* Compute the (r,s) entry of X.T   */
{
  colrange j;
  double temp;
  
  temp=0;
  for (j=0; j< n_size; j++) {
    temp = temp + X[r-1][j] * T[j][s-1];
  }
  return temp;
}

void WriteTableau(FILE *f, rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
  InequalityType ineq)
/* Write the tableau  X.T   */
{
  colrange j;
  rowrange i;
  
  fprintf(f, "begin\n");
  if (ineq==ZeroRHS)
    fprintf(f, "  %ld   %ld    real\n",m_size, n_size+1);
  else
    fprintf(f, "  %ld   %ld    real\n",m_size, n_size);
  for (i=1; i<= m_size; i++) {
    if (ineq==ZeroRHS)
      fprintf(f," %5d", 0);  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j<= n_size; j++) {
      fprintf(f," %12.6f",TableauEntry(m_size, n_size, X,T,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
}


void SelectPivot2(rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
            HyperplaneOrderType roworder, rowindex ordervec,
            rowrange rowmax, rowset NopivotRow,
            colset NopivotCol, rowrange *r, colrange *s,
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
  boolean localdebug=FALSE;

  if (debug) localdebug=TRUE;
  stop = FALSE;
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  if (localdebug) {
    switch (roworder) {

    case MinIndex:
      fprintf(stdout, "*SelectPivot2: MinIndex\n");
      break;

    case MaxIndex:
      fprintf(stdout, "*SelectPivot2: MaxIndex\n");
      break;

    case MinCutoff:
      fprintf(stdout, "*SelectPivot2: MinCutoff\n");
      break;

    case MaxCutoff:
      fprintf(stdout, "*SelectPivot2: MaxCutoff\n");
    break;

    case MixCutoff:
      fprintf(stdout, "*SelectPivot2: MixCutoff\n");
      break;

    case LexMin:
      fprintf(stdout, "*SelectPivot2: LexMin\n");
      break;

    case LexMax:
      fprintf(stdout, "*SelectPivot2: LexMax\n");
      break;

    case RandomRow:
      fprintf(stdout, "*SelectPivot2: Random,  Seed = %d\n",rseed);
      break;

    case LineShelling:
      fprintf(stdout, "*SelectPivot2: LineShelling\n");
      break;
    }
    printf("select pivot2: rowexcluded=");
    set_fwrite(stdout,rowexcluded);
  }
  for (rtemp=rowmax+1;rtemp<=m_size;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = FALSE;
  do {
    i=1;rtemp=0;
    while (i<=m_size && rtemp==0) {  /* EqualitySet vars have highest priorities */
      if (set_member(i,EqualitySet) && !set_member(i,rowexcluded)){
        if (localdebug) printf("marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) SelectNextHyperplane(m_size, n_size, X,
       roworder, rowexcluded, &rtemp, &RecomputeRowOrder, ordervec);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= n_size && !*selected) {
        Xtemp=TableauEntry(m_size, n_size,X,T,*r,*s);
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

void GausianColumnPivot1(rowrange m_size, colrange n_size, 
    Amatrix X, rowrange r, colrange s,
    rowrange rowmax)
/* Make a column pivot operation in Amatrix X on position (r,s)  */
{
  long i, j;
  double Xtemp0, Xtemp;

  Xtemp0 = X[r - 1][s - 1];
  for (j = 0; j < n_size; j++) {
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


void GausianColumnPivot2(rowrange m_size, colrange n_size, 
    Amatrix X, Bmatrix T, rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  long j, j1;
  Arow Rtemp;
  double Xtemp0, Xtemp;

  for (j=1; j<=n_size; j++) Rtemp[j-1]=TableauEntry(m_size, n_size, X, T, r,j);
  Xtemp0 = Rtemp[s-1];
  for (j = 1; j <= n_size; j++) {
    if (j != s) {
      Xtemp = Rtemp[j-1];
      for (j1 = 1; j1 <= n_size; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= n_size; j++)
    T[j-1][s - 1] /= Xtemp0;
}


void InitializeBmatrix(colrange n_size, Bmatrix T)
{
  colrange j;
  double x;

  for (j = 0; j < n_size; j++) {
    T[j]=(double *)calloc(n_size, sizeof x);
  }
}

void free_Bmatrix(colrange n_size, Bmatrix T)
{
  colrange j;

  for (j = 0; j < n_size; j++) {
    free(T[j]);
  }
}

void SetToIdentity(colrange n_size, Bmatrix T)
{
  colrange j1, j2;

  for (j1 = 1; j1 <= n_size; j1++) {
    for (j2 = 1; j2 <= n_size; j2++) {
      if (j1 == j2)
        T[j1 - 1][j2 - 1] = 1.0;
      else
        T[j1 - 1][j2 - 1] = 0.0;
    }
  }
}

void ReduceA(rowrange *m_size, colrange *n_size, Amatrix A, rowset ChosenRow, colset ChosenCol)
/* Set the matrix AA to be the submatrix of AA with chosen 
   rows & columns, and change m_size and n_size accordingly 
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
  for (i=1; i <= *m_size; i++) {
    if (set_member(i,ChosenRow)) {
      inew++;
      jnew=0;
      for (j=1; j <= *n_size; j++) {
        if (set_member(j, ChosenCol)){
          jnew++;
          Acopy[inew-1][jnew-1]=A[i-1][j-1];
          if (debug) WriteReal(stdout, A[i-1][j-1]);
        }
      }
      if (debug) fprintf(stdout,"\n");
    }
  }
  for (i=1;i<= *m_size;i++) {
    if (i<=mnew) 
      set_addelem(ChosenRow,i);
    else
      set_delelem(ChosenRow,i);
  }
  for (j=1;j<= *n_size;j++) {
    if (j<=nnew) 
      set_addelem(ChosenCol,j);
    else
      set_delelem(ChosenCol,j);
  }
  if (debug) {
    fprintf(stdout, "new row indices:");set_write(ChosenRow);
    fprintf(stdout, "new col indices:");set_write(ChosenCol);
  }
  for (i=0;i< *m_size;i++){
    free(A[i]);
  }
  for (i=0;i< mnew;i++){
    A[i]=Acopy[i];
  }
  *m_size=mnew;  *n_size=nnew;
}

void DualizeA(rowrange *m_size, colrange *n_size, Amatrix A, Bmatrix T)
/* Set the matrix A to be the transpose of the matrix [-A.T | I],
   and change m_size and n_size accordingly 
*/
{
  long i,j,mnew,nnew;
  Amatrix Acopy;
  double x;

  mnew=*n_size+*m_size;
  nnew=*m_size;
  for (i=0; i< *m_size; i++){
    Acopy[i]=(double *)calloc(*n_size,sizeof x);
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
  for (i=1;i<= *m_size;i++){
    for (j=1;j<= *n_size;j++){
      Acopy[i-1][j-1]=TableauEntry(*m_size, *n_size, A,T,i,j);
    }
  }
  for (i=0;i< *m_size;i++){
    free(A[i]);
  }
  for (i=0; i<mnew; i++){
    A[i]=(double *)calloc(nnew,sizeof x);
  }
  for (i=1;i<= *n_size;i++){
    for (j=1;j<= *m_size;j++){
      A[i-1][j-1]=-Acopy[j-1][i-1];
    }
  }
  for (i=1;i<= *m_size;i++){
    for (j=1;j<= *m_size;j++){
      if (i==j) A[*n_size+i-1][j-1]=1;
      else A[*n_size+i-1][j-1]=0;
    }
  }
  for (i=0; i< *m_size; i++){
    free(Acopy[i]);
  }
  *m_size=mnew;  *n_size=nnew;
  _L99:;
}

void EnlargeAforInteriorFinding(rowrange *m_size, colrange *n_size, Amatrix A)
/* Add an extra column with all minus ones to the matrix AA, 
   add an objective row with (0,...,0,1), and 
   rows & columns, and change m_size and n_size accordingly 
*/
{
  long i,j,mnew=0,nnew=0;
  Amatrix Acopy;
  double x;

  mnew=*m_size+1;
  nnew=*n_size+1;
  for (i=0; i<mnew; i++){
    Acopy[i]=(double *)calloc(nnew,sizeof x);
  }
  for (i=1; i <= *m_size; i++) {
    for (j=1; j <= *n_size; j++) {
      Acopy[i-1][j-1]=A[i-1][j-1];
      if (debug) WriteReal(stdout, A[i-1][j-1]);
    }
    if (debug) fprintf(stdout,"\n");
  }
  for (i=1;i<= *m_size;i++) {
    Acopy[i-1][*n_size]=-1.0;  /* new column with all minus one's */
  }
  for (j=1;j<= *n_size;j++) {
    Acopy[*m_size][j-1]=0.0;  /* new row with (0,...,0,1) */
  }
  Acopy[*m_size][*n_size]=1.0;  /* new row with (0,...,0,1) */
  for (i=0;i<*m_size;i++){
    free(A[i]);
  }
  for (i=0;i< mnew;i++){
    A[i]=Acopy[i];
  }
  *m_size=mnew;  *n_size=nnew;
}

void ComputeRank(rowrange m_size, colrange n_size, 
    Amatrix A1, rowset TargetRows, rowindex ordervec, long *rank)
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
  set_initialize(&NoPivotRow, m_size);
  set_compl(NoPivotRow,TargetRows);
  set_initialize(&ColSelected, n_size);
  InitializeBmatrix(n_size, Btemp);
  SetToIdentity(n_size, Btemp);
  if (debug) WriteBmatrix(stdout, n_size, Btemp);
  do {   /* Find a set of rows for a basis */
      SelectPivot2(m_size, n_size, A1, Btemp, MinIndex, ordervec, 
        m_size, NoPivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(NoPivotRow, r);
        set_addelem(ColSelected, s);
        (*rank)++;
        GausianColumnPivot2(m_size, n_size, A1,Btemp, r, s);
        if (debug) {
          WriteBmatrix(stdout, n_size, Btemp);
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
      } else {
        stop=TRUE;
      }
  } while (!stop);
  set_free(NoPivotRow);
  set_free(ColSelected);
  free_Bmatrix(n_size, Btemp);
}

void ComputeBInverse(rowrange m_size, colrange n_size, Amatrix A1, long lastrow,
       rowindex ordervec, Bmatrix InvA1, long *rank)
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset RowSelected;
  colset ColSelected;

  *rank = 0;
  stop = FALSE;
  SetToIdentity(n_size, InvA1);
  set_initialize(&RowSelected, m_size);
  set_initialize(&ColSelected, n_size);
  do {
    SelectPivot2(m_size, n_size, A1, InvA1, MinIndex, ordervec,
      lastrow, RowSelected, ColSelected, &r, &s, &chosen);
    if (chosen) {
      (*rank)++;
      if (debug)
        printf("%3ldth pivot on%3ld, %3ld\n", *rank, r, s);
      GausianColumnPivot2(m_size, n_size, A1, InvA1, r, s);
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
    } else
      stop = TRUE;
  } while (!stop);
  set_free(RowSelected);
  set_free(ColSelected);
}


void FindBasis(rowrange m_size, colrange n_size, Amatrix A1, 
    HyperplaneOrderType roword, rowindex ordervec,
    rowset RowSelected,colindex ColInd,
    Bmatrix BasisInverse, long *rank)
{
  boolean stop, chosen;
  rowset NopivotRow;
  colset ColSelected;
  rowrange r;
  colrange j,s;

  *rank = 0;
  stop = FALSE;
  for (j=0;j<=n_size;j++) ColInd[j]=0;
  set_emptyset(RowSelected);
  set_initialize(&ColSelected, n_size);
  set_initialize(&NopivotRow, m_size);
  set_copy(NopivotRow,NonequalitySet);
  SetToIdentity(n_size, BasisInverse);
  if (debug) WriteBmatrix(stdout, n_size, BasisInverse);
  if (DynamicWriteOn && !debug){
    printf("*Initial set of rows:");
    fprintf(writing,"*Initial set of rows:");
  }
  do {   /* Find a set of rows for a basis */
      SelectPivot2(m_size, n_size, A1, BasisInverse, roword, ordervec, m_size, 
        NopivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) 
        printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(NopivotRow, r);
        set_addelem(ColSelected, s);
        ColInd[s]=r;    /* ColInd[s] stores the corr. row index */
        (*rank)++;
        GausianColumnPivot2(m_size, n_size, A1,BasisInverse, r, s);
        if (debug) {
          WriteBmatrix(stdout,n_size, BasisInverse);
          WriteTableau(stdout, m_size, n_size, A1,BasisInverse,NonzeroRHS),
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
	    if (DynamicWriteOn && !debug){
	      printf(" %ld",  r);
	      fprintf(writing, " %ld",  r);
	    }
      } else {
        stop=TRUE;
      }
      if (*rank==n_size) stop = TRUE;
  } while (!stop);
  if (DynamicWriteOn && !debug){
    printf("\n");
    fprintf(writing,"\n");
  }
  set_free(ColSelected);
}



void FindInitialRays(rowrange m_size, colrange n_size, Amatrix A, rowindex ordervec,
    rowset InitHyperplanes, Bmatrix InitRays, colindex PivRow, boolean *found)
{
  Bmatrix BInverse;
  rowset CandidateRows;
  long i, j, rank;
  HyperplaneOrderType roworder=LexMin;

  *found = FALSE;
  set_initialize(&CandidateRows, m_size);
  if (InitBasisAtBottom==TRUE) {
    roworder=MaxIndex;
    PreOrderedRun=FALSE;
  }
  else PreOrderedRun=TRUE;
  for (i = 1; i <= m_size; i++)
    if (!set_member(i,NonequalitySet)) set_addelem(CandidateRows, i);
    /*all rows not in NonequalitySet are candidates for initial cone*/
  if (DynamicWriteOn)
    printf("*Computing an initial set of rays\n");
  InitializeBmatrix(n_size, BInverse);
  FindBasis(m_size, n_size, A, roworder, ordervec, InitHyperplanes, 
    PivRow, BInverse, &rank);
  if (debug) {
    printf("FindInitialBasis: InitHyperplanes=");
    set_fwrite(stdout,InitHyperplanes);
    for (j = 1; j <= n_size; j++)
      printf("Pivotrow[%ld] = %ld \n", j, PivRow[j]);
    printf("n_size = %ld, rank = %ld\n",n_size,rank);
  }
  if (rank < n_size-1) {
    if (debug) WriteBmatrix(stdout, n_size, BInverse);
    Error = LowColumnRank;
    return;
  }
  if (!set_subset(EqualitySet,InitHyperplanes)) {
    Error = DependentMarkedSet;
    return;
  }
  *found = TRUE;
  if (debug) {
    WriteBmatrix(stdout, n_size, BInverse);
  }
  /* free_Bmatrix(n_size, InitRays);  */
  CopyBmatrix(n_size, BInverse,InitRays);
  if (debug) WriteBmatrix(stdout, n_size, InitRays);
  set_free(CandidateRows);
  if (HyperplaneOrder==MaxCutoff||HyperplaneOrder==MinCutoff||HyperplaneOrder==MixCutoff){
    PreOrderedRun=FALSE;
  } else PreOrderedRun=TRUE;
}

void CheckEquality(colrange n_size, RayRecord **RP1, RayRecord **RP2, boolean *equal)
{
  long j;

  if (debug)
    printf("Check equality of two rays\n");
  *equal = TRUE;
  j = 1;
  while (j <= n_size && *equal) {
    if (fabs((*RP1)->Ray[j - 1] - (*RP2)->Ray[j - 1]) > 2.0 * zero)
      *equal = FALSE;
    j++;
  }
  if (*equal)
    printf("Equal records found !!!!\n");
}

void CreateNewRay(rowrange m_size, colrange n_size, Amatrix A, 
    RayRecord *Ptr1, RayRecord *Ptr2, rowrange ii, rowindex ordervec)
{
  /*Create a new ray by taking a linear combination of two rays*/
  colrange j;
  Arow NewRay;
  double v1, v2;

  v1 = fabs(AValue(n_size, A, Ptr1->Ray, ii));
  v2 = fabs(AValue(n_size, A, Ptr2->Ray, ii));
  for (j = 0; j < n_size; j++)
    NewRay[j] = Ptr1->Ray[j] * v2 + Ptr2->Ray[j] * v1;
  Normalize(n_size, NewRay);
  AddRay(m_size, n_size, A, NewRay, ordervec);
  if (debug){
    WriteRayRecord(stdout, n_size, Ptr1);
    WriteRayRecord(stdout, n_size, Ptr2);  
    printf("create a new ray by eliminating %ld:\n",ii);
    WriteRayRecord(stdout, n_size, LastRay);
  }
}

void EvaluateARay1(rowrange i, colrange n_size, Amatrix A)
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
    for (j = 0; j < n_size; j++)
      temp += A[i - 1][j] * Ptr->Ray[j];
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

void EvaluateARay2(rowrange i, colrange n_size, Amatrix A)
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

  PosHead=NULL;ZeroHead=NULL;NegHead=NULL;
  PosLast=NULL;ZeroLast=NULL;NegLast=NULL;
  Ptr = FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=NULL;     /* then clear the Next pointer */
    temp = 0.0;
    for (j = 0; j < n_size; j++)
      temp += A[i - 1][j] * Ptr->Ray[j];
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

void DeleteNegativeRays(colrange n_size)
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
  if (PrevPtr->Next != Ptr) 
    printf("Error at DeleteNegativeRays: ArtificialRay does not point the FirstRay.\n");
  completed=FALSE;
  while (Ptr != NULL && !completed){
    if ( (Ptr->ARay) < -zero ){
      Eliminate(n_size, &PrevPtr);
      Ptr=PrevPtr->Next;
    }
    else{
      completed=TRUE;
    }
  }
  
  /* Sort the zero rays */
  Ptr = FirstRay;
  ZeroRayCount=0;
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

void FeasibilityIndices(long *fnum, long *infnum, rowrange i, colrange n_size, Amatrix A)
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
    for (j = 0; j < n_size; j++)
      temp += A[i - 1][j] * Ptr->Ray[j];
    if (temp >= 0)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}

boolean LexSmaller(double *v1, double *v2, long nmax)
{ /* nmax is the size of vectors v1,v2 */
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
  } while (!(determined) && (j <= nmax));
  return smaller;
}


boolean LexLarger(double *v1, double *v2, long nmax)
{
  Arow u1, u2;
  colrange j;

  for (j = 1; j <= nmax; j++) {
    u1[j-1] = -v1[j-1];
    u2[j-1] = -v2[j-1];
  }
  return (LexSmaller(u1, u2, nmax));
}

void CopyArow(double *vcopy, double *v, long nmax)
{
 colrange j;

  for (j = 1; j <= nmax; j++) {
    vcopy[j-1] = v[j-1];
  }
}

void AddNewHyperplane1(rowrange m_size, colrange n_size, 
     Amatrix A, rowrange hnew, rowrange iter, rowindex ordervec)
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

  EvaluateARay1(hnew, n_size, A);  
   /*Check feasibility of rays w.r.t. hnew 
     and put all infeasible ones consecutively */

  RayPtr0 = ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  value1 = FirstRay->ARay;
  if (value1 > -zero) {
    if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
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
    RayCount=0;
    VertexCount=0;
    if (Inequality==ZeroRHS && Conversion==IneToExt){
      CompStatus=AllFound;
    } else {
      CompStatus=RegionEmpty;
    }
     goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;   /*2nd feasible ray to scan and compare with 1st*/
  RayPtr3 = LastRay;    /*Last feasible for scanning*/
  prevprogress=-10.0;
  pos1 = 1;
  completed=FALSE;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    value1 = RayPtr1->ARay;
    value2 = RayPtr2->ARay;
    if (debug) {
      WriteRayRecord2(stdout, n_size, RayPtr1);
      WriteRayRecord2(stdout, n_size, RayPtr2);
      printf("check feasibility%8.3f%8.3f\n", value1, value2);
    }
    CheckEquality(n_size, &RayPtr1, &RayPtr2, &equal);
    if ((value1 >= zero && value2 <= -zero) || (value2 >= zero && value1 <= -zero)) {
      switch (AdjacencyTest) {

      case Algebraic:
		CheckAdjacency1(m_size, n_size, A, ordervec, &RayPtr1, &RayPtr2, &adj);
		break;

      case Combinatorial:
		CheckAdjacency2(m_size, n_size, A, &RayPtr1, &RayPtr2, &adj);
		break;
      }
      if (adj)
		CreateNewRay(m_size, n_size, A, RayPtr1, RayPtr2, hnew, ordervec);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
    if (value1 <= -zero || equal) {
      Eliminate(n_size, &RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed=TRUE;
    }
    pos1++;
    progress = 100.0 * ((double)pos1 / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (progress-prevprogress>=10 && pos1%10==0 && DynamicWriteOn) {
      printf("*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	     iter, m_size, pos1, pos2, progress);
      fprintf(writing,
	  "*Progress of iteration %5ld(/%ld):   %4ld/%4ld => %4.1f%% done\n",
	  iter, m_size, pos1, pos2, progress);
      prevprogress=progress;
      fflush(writing);
      if (writing_icd != NULL) fflush(writing_icd);
    }
  }
  if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
  _L99:;
}

void AddNewHyperplane2(rowrange m_size, colrange n_size, Amatrix A, 
    rowrange hnew, rowrange iter, rowindex ordervec)
/* This procedure must be used under PreOrderedRun mode */
{
  RayRecord *RayPtr0,*RayPtr1, *RayPtr2;
  AdjacencyRecord *EdgePtr, *EdgePtr0;
  long pos1;
  rowrange fii1, fii2;
  boolean localdebug=FALSE;

  EvaluateARay2(hnew, n_size, A);
   /* Check feasibility of rays w.r.t. hnew 
      and sort them. ( -rays, +rays, 0rays)*/

  if (PosHead==NULL && ZeroHead==NULL) {
    FirstRay=NULL;
    ArtificialRay->Next=FirstRay;
    RayCount=0;
    VertexCount=0;
    if (Inequality==ZeroRHS && Conversion==IneToExt){
      CompStatus=AllFound;
    } else {
      CompStatus=RegionEmpty;
    }
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

  EdgePtr=Edges[iter];
  while (EdgePtr!=NULL){
    RayPtr1=EdgePtr->Ray1;
    RayPtr2=EdgePtr->Ray2;
    fii1=RayPtr1->FirstInfeasIndex;   
    CreateNewRay(m_size, n_size, A, RayPtr1, RayPtr2, hnew, ordervec);
    fii2=LastRay->FirstInfeasIndex;
    if (fii1 != fii2) 
      ConditionalAddEdge(m_size, n_size, RayPtr1,LastRay,PosHead,iter,ordervec);
    EdgePtr0=EdgePtr;
    EdgePtr=EdgePtr->Next;
    free(EdgePtr0);
    EdgeCount--;
  }
  Edges[iter]=NULL;
  
  DeleteNegativeRays(n_size);
    
  set_addelem(AddedHyperplanes, hnew);

  if (iter<m_size){
    if (ZeroHead!=NULL && ZeroHead!=LastRay){
      if (ZeroRayCount>200 && DynamicWriteOn) printf("*New edges being scanned...\n");
      UpdateEdges(m_size, n_size, ZeroHead, LastRay, iter, ordervec);
    }
    if (localdebug) printf("*Edges currently stored = %ld\n", EdgeCount);
  }

  if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
_L99:;
}


void SelectNextHyperplane0(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*A natural way to choose the next hyperplane.  Simply the largest index*/
  long i;
  boolean determined;

  i = m_size;
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

void SelectNextHyperplane1(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
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
  } while (!determined && i<=m_size);
  if (determined) 
    *hnext = i;
  else 
    *hnext=0;
}

void SelectNextHyperplane2(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmin, fi=0;   /*feasibility and infeasibility numbers*/

  infmin = RayCount + 1;
  for (i = 1; i <= m_size; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, n_size, A);
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

void SelectNextHyperplane3(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax, fi=0;   /*feasibility and infeasibility numbers*/

  infmax = -1;
  for (i = 1; i <= m_size; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, n_size, A);
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

void SelectNextHyperplane4(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with the most unbalanced cut*/
  long i, fea, inf, max, tmax, fi=0, infi=0;
      /*feasibility and infeasibility numbers*/

  max = -1;
  for (i = 1; i <= m_size; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i, n_size, A);
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

void SelectNextHyperplane5(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-min*/
  long i, minindex;
  double *v1, *v2;

  minindex = 0;
  v1 = NULL;
  for (i = 1; i <= m_size; i++) {
    if (!set_member(i, excluded)) {
	  v2 = A[i - 1];
      if (minindex == 0) {
	    minindex = i;
	    v1=v2;
      } else if (LexSmaller(v2,v1,n_size)) {
        minindex = i;
	    v1=v2;
      }
    }
  }
  *hnext = minindex;
}


void SelectNextHyperplane6(rowrange m_size, colrange n_size, Amatrix A,
    rowset excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  double *v1, *v2;

  maxindex = 0;
  v1 = NULL;
  for (i = 1; i <= m_size; i++) {
    if (!set_member(i, excluded)) {
      v2= A[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1=v2;
      } else if (LexLarger(v2, v1, n_size)) {
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

void LineShellingOrder(rowrange m_size, colrange n_size, Amatrix A, rowindex OV, double *z, double *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  double temp1,temp2,infinity=10.0e+20;
  static double *beta;
  static long mlast=0;
  boolean localdebug=FALSE;
  
  if ( mlast<m_size ){
    if (beta!=NULL) free(beta);
    beta=(double *)calloc(m_size, sizeof *beta);
    /* initialize only for the first time or when last m_size is smaller */
    mlast=m_size;
  }
  for (i=1; i<= m_size; i++) beta[i-1]=A[i-1][0]; /* store the first column in beta */
  for (i=1; i<= m_size; i++){
    temp1 = 0.0;
    temp2 = 0.0;
    for (j = 1; j <= n_size; j++){
      temp1 += A[i - 1][j-1] * z[j-1];
      temp2 += A[i - 1][j-1] * d[j-1];
    }
    if (abs(temp1)>zero) A[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > 0) A[i-1][0]= infinity;
    else A[i-1][0]= -infinity;
     /* use the first column of A tentatively */
  }
  if (localdebug) 
    for (i=1; i<= m_size; i++){
      printf("set A[%ld] = %g\n", i, A[i-1][0]);
    }
  QuickSort(OV, 1, m_size, A, 1);
  for (i=1; i<= m_size; i++) {
    A[i-1][0]=beta[i-1]; 
     /* restore the first column of A */ 
    if (localdebug) printf("restore A[%ld] with %g\n", i, A[i-1][0]);
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
    if (localdebug) printf("u=%g, k=%ld, r=%g, randmax= %g\n",u,k,r,rand_max);
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
    if (localdebug) printf("row %ld is exchanged with %ld\n",j,k); 
  }
}

void ComputeRowOrderVector(rowrange m_size, colrange n_size, Amatrix A,
    rowindex OV, HyperplaneOrderType ho)
{
  long i,itemp,j;
  Arow zvec, dvec;
  
  OV[0]=0;
  switch (ho){
  case MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case MinIndex: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  case LexMin: case MinCutoff: case MixCutoff: case MaxCutoff:
    for(i=1; i<=m_size; i++) OV[i]=i;
    QuickSort(OV, 1, m_size, A, n_size);
    break;

  case LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    QuickSort(OV, 1, m_size, A, n_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    RandomPermutation(OV, m_size, rseed);
    break;

  case LineShelling:
    for(i=1; i<=m_size; i++) OV[i]=i;
    zvec[0]=1;
    dvec[0]=0;
    if (rseed<=0) rseed=1;
    srand(rseed);
    for(j=2; j<=n_size; j++){
      zvec[j-1]=0;
      dvec[j-1]=n_size-j+1;
      /* dvec[j-1]=rand(); */
    }
    LineShellingOrder(m_size, n_size, A, OV, zvec, dvec);
    break;
  }
}

void UpdateRowOrderVector(rowrange m_size, colrange n_size,
    rowset PriorityRows, rowindex ordervec)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  rowrange i,j,k,j1=0,oj=0;
  long rr;
  boolean found, localdebug=FALSE;
  
  if (debug) localdebug=TRUE;
  found=TRUE;
  rr=set_card(PriorityRows);
  if (localdebug) set_write(PriorityRows);
  for (i=1; i<=rr; i++){
    found=FALSE;
    for (j=i; j<=m_size && !found; j++){
      oj=ordervec[j];
      if (set_member(oj, PriorityRows)){
        found=TRUE;
        if (localdebug) printf("%ldth in sorted list (row %ld) is in PriorityRows\n", j, oj);
        j1=j;
      }
    }
    if (found){
      if (j1>i) {
        /* shift everything lower: OV[i]->OV[i+1]..OV[j1-1]->OV[j1] */
        for (k=j1; k>=i; k--) ordervec[k]=ordervec[k-1];
        ordervec[i]=oj;
        if (localdebug){
          printf("OrderVector updated to:\n");
          for (j = 1; j <= m_size; j++) printf(" %2ld", ordervec[j]);
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

void SelectPreorderedNext(rowrange m_size, colrange n_size, 
    rowset excluded, rowindex OV, rowrange *hnext)
{
  rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k, excluded)) *hnext=k ;
  }
}

void SelectNextHyperplane(rowrange m_size, colrange n_size, Amatrix A, HyperplaneOrderType ho, 
         rowset excluded, rowrange *hh, boolean *RefreshOrderVector, rowindex ordervec)
{
  if (PreOrderedRun){
    if (debug) {
      printf("debug SelectNextHyperplane: Use PreorderNext\n");
    }
    SelectPreorderedNext(m_size, n_size, excluded, ordervec, hh);
  }
  else {
    if (debug) {
      printf("debug SelectNextHyperplane: Use DynamicOrderedNext\n");
    }

    switch (ho) {

    case MaxIndex:
      SelectNextHyperplane0(m_size, n_size, A, excluded, hh);
      break;

    case MinIndex:
      SelectNextHyperplane1(m_size, n_size, A, excluded, hh);
      break;

    case MinCutoff:
      SelectNextHyperplane2(m_size, n_size, A, excluded, hh);
      break;

    case MaxCutoff:
      SelectNextHyperplane3(m_size, n_size, A, excluded, hh);
      break;

    case MixCutoff:
      SelectNextHyperplane4(m_size, n_size, A, excluded, hh);
      break;

    default:
      SelectNextHyperplane0(m_size, n_size, A, excluded, hh);
      break;
    }
  }
}

/* end of cddarith.c */


