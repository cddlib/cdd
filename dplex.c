/* dplex.c:  dual simplex method c-code
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.61, December 1, 1997
*/

/* dplex.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  c^T x subject to  x in P, where
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "setoper.h"  /* set operation library header (Ver. March 16,1995 or later) */
#include "dplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define COPYRIGHT   "Copyright (C) 1997, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DPLEXVERSION   "Version 0.61 (December 1, 1997)"

#define dp_FALSE 0
#define dp_TRUE 1

typedef enum {
  dp_MaxIndex, dp_MinIndex, dp_LexMin, dp_LexMax, dp_RandomRow, dp_LineShelling
} dp_HyperplaneOrderType;

void dp_LPInput(FILE **f, dp_FilenameType, rowrange *m, colrange *n, Amatrix A, 
    dp_LPConversionType *Conversion, rowrange *objrow, colrange *rhscol, 
    dp_ErrorType *err);
void dp_InitializeBmatrix(colrange n_size, Bmatrix T);
void CrissCrossMinimize(rowrange, colrange, Amatrix, Bmatrix T,
  rowrange, colrange, int,
  dp_LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *, dp_ErrorType *);
void CrissCrossMaximize(rowrange, colrange, Amatrix,Bmatrix T,
  rowrange, colrange, int,
  dp_LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *, dp_ErrorType *);
void DualSimplexMinimize(rowrange, colrange,  Amatrix,Bmatrix,
  rowrange, colrange, int,
  dp_LPStatusType *, double *, Arow, Arow, colindex,
  rowrange *, colrange *, long *, dp_ErrorType *);
void DualSimplexMaximize(rowrange, colrange, Amatrix,Bmatrix,
  rowrange, colrange, int,
  dp_LPStatusType *, double *, Arow, Arow, colindex,
  rowrange *, colrange *, long *, dp_ErrorType *);
void FindLPBasis(rowrange, colrange, Amatrix, Bmatrix, rowindex, 
    colindex, rowindex, rowrange, colrange, int,
    colrange *, int *, dp_LPStatusType *, long *);
void FindDualFeasibleBasis(rowrange, colrange, Amatrix, Bmatrix, rowindex, 
    colindex, long *, rowrange, colrange,
    colrange *, int *, dp_LPStatusType *, long *);

void dp_WriteBmatrix(FILE *f, colrange n_size, Bmatrix T);
void dp_SetNumberType(char *line, dp_NumberType *number, dp_ErrorType *Error);
void dp_ComputeRowOrderVector(rowrange m_size, colrange n_size, Amatrix A,
    rowindex OV, dp_HyperplaneOrderType ho, unsigned int rseed);
void dp_SelectPreorderedNext(rowrange m_size, colrange n_size, 
    rowset excluded, rowindex OV, rowrange *hnext);
void dp_WriteReal(FILE *f, double x);
void dp_SetInputFile(FILE **f, dp_FilenameType inputfile,  dp_ErrorType *);
void SetSolutions(rowrange, colrange,
   Amatrix, Bmatrix, rowrange, colrange, dp_LPStatusType,
   double *, Arow, Arow, colindex, rowrange, colrange);

void dp_SetInputFile(FILE **f, dp_FilenameType inputfile, dp_ErrorType *Error)
{
  int opened=0,stop, quit=0;
  int i,dotpos=0,trial=0;
  char ch;
  char *tempname;
  
  
  *Error=dp_None;
  while (!opened && !quit) {
    printf("\n>> Input file (*.ine) : ");
    scanf("%s",inputfile);
    ch=getchar();
    stop=dp_FALSE;
    for (i=0; i<dp_filenamelen && !stop; i++){
      ch=inputfile[i];
      switch (ch) {
        case '.': 
          dotpos=i+1;
          break;
        case ';':  case ' ':  case '\0':  case '\n':  case '\t':     
          stop=dp_TRUE;
          tempname=(char*)calloc(dp_filenamelen,sizeof ch);
          strncpy(tempname, inputfile, i);
          strcpy(inputfile,tempname);
          break;
      }
    }
    if ( ( *f = fopen(inputfile, "r") )!= NULL) {
      printf("input file %s is open\n", inputfile);
      opened=1;
      *Error=dp_None;
    }
    else{
      printf("The file %s not found\n",inputfile);
      trial++;
      if (trial>5) {
        *Error=dp_FileNotFound;
        quit=1;
      }
    }
  }
}

void dp_LPInput(FILE **f, dp_FilenameType inputfile,
    rowrange *m_size, colrange *n_size, Amatrix A, 
    dp_LPConversionType *Conversion, rowrange *objrow, colrange *rhscol, 
    dp_ErrorType *err)
{
  long i,j;
  rowrange m_input;
  colrange n_input;
  double value, cost;
  dp_NumberType Number;
  int found=0, localdebug=0;
  char command[dp_wordlenmax], numbtype[dp_wordlenmax], line[dp_linelenmax];

  *err=dp_None;

  dp_SetInputFile(f, inputfile, err);

  if (*err!=dp_None){
    goto _L99;
  }

  while (!found)
  {
    if (fscanf(*f,"%s",command)==EOF) {
      *err=dp_ImproperInputFormat;
      goto _L99;
    }
    else if (strncmp(command, "begin", 5)==0) {
      found=dp_TRUE;
    }
  }
  fscanf(*f, "%ld %ld %s", &m_input, &n_input, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", m_input, n_input, numbtype);
  dp_SetNumberType(numbtype, &Number, err);
  if (Number==dp_Unknown || Number == dp_Rational) {
      goto _L99;
  }
  *n_size=n_input; *m_size=m_input+1;
  if (*n_size > dp_NMAX || *m_size+1 > dp_MMAX) {
    *err = dp_DimensionTooLarge;
    goto _L99;
  }
  
  for (i = 1; i <= m_input; i++) {
    A[i-1]=(double *) calloc(n_input, sizeof value);
    for (j = 1; j <= n_input; j++) {
      fscanf(*f, "%lf", &value);
      A[i-1][j-1] = value;
      if (localdebug) printf("a(%3ld,%5ld) = %10.4f\n",i,j,value);
    }  /*of j*/
    fgets(line,dp_linelenmax,*f);
    if (localdebug) printf("comments to be skipped: %s\n",line);
    if (localdebug) putchar('\n');
  }  /*of i*/
  if (fscanf(*f,"%s",command)==EOF) {
   	 *err=dp_ImproperInputFormat;
  	 goto _L99;
  }
  else if (strncmp(command, "end", 3)!=0) {
     if (localdebug) printf("'end' missing or illegal extra data: %s\n",command);
   	 *err=dp_ImproperInputFormat;
  	 goto _L99;
  }
  
  *objrow=*m_size;
  *rhscol=1L;
  A[*m_size-1]=(double *) calloc(n_input, sizeof value);
  found=0;
  while (!found)
  {
    if (fscanf(*f,"%s",command)==EOF) {
      *err=dp_ImproperInputFormat;
      goto _L99;
    }
    if (strncmp(command, "maximize", 8)==0) {
      *Conversion=dp_LPmax;
      found=dp_TRUE;
    }
    if (strncmp(command, "minimize", 8)==0) {
      *Conversion=dp_LPmin;
      found=dp_TRUE;
    }
  }
  for (j=0;j<n_input;j++) {
    fscanf(*f,"%lf",&cost);
    A[*m_size-1][j]=cost;
    if (localdebug) printf(" cost[%ld] = %.9E\n",j,A[*m_size-1][j]);
  }
_L99: ;
  if (*f!=NULL) fclose(*f);
}

int dp_Nonnegative(double val)
{
  if (val>=-dp_zero) return dp_TRUE;
  else return dp_FALSE;
}

int dp_Nonpositive(double val)
{
  if (val<=dp_zero) return dp_TRUE;
  else return dp_FALSE;
}

int dp_Positive(double val)
{
  return !dp_Nonpositive(val);
}

int dp_Negative(double val)
{
  return !dp_Nonnegative(val);
}

int dp_Zero(double val)
{
  return (dp_Nonnegative(val) && dp_Nonpositive(val));
}

int dp_Nonzero(double val)
{
  return (dp_Positive(val) || dp_Negative(val));
}


void dp_SetNumberType(char *line, dp_NumberType *number, dp_ErrorType *Error)
{
  if (strncmp(line, "integer", 7)==0) {
    *number = dp_Integer;
    return;
  }
  else if (strncmp(line, "rational", 8)==0) {
    *number = dp_Rational;
    *Error=dp_ImproperInputFormat;  /* Rational Input not supported */
    return;
  }
  else if (strncmp(line, "real", 4)==0) {
    *number = dp_Real;
    return;
  }
  else { 
    *number=dp_Unknown;
    *Error=dp_ImproperInputFormat;
  }
}

void dp_WriteErrorMessages(FILE *f, dp_ErrorType Error)
{
  switch (Error) {
    case dp_DimensionTooLarge:
      fprintf(f, "dp_Error: LP size is too large.  Modify dp_NMAX and/or dp_MMAX.\n");
      break;
      
    case dp_LowColumnRank:
      fprintf(f, "dp_Error: The matrix A is not column full rank.\n");
      break;
 
    case dp_ImproperInputFormat:
      fprintf(f, "dp_Error: Input file format is not correct.\n");
      break;

    case dp_FileNotFound:
      fprintf(f, "dp_Error: The input file does not exist.\n");
      break;
    
    case dp_None:
      fprintf(f, "dp_Error: No error found.\n");
      break;

    default:
      fprintf(f, "dp_Error: Unknown error found.\n");
      break;
  }
}


double dp_TableauEntry(rowrange m_size, colrange n_size, Amatrix A, Bmatrix T,
				rowrange r, colrange s)
/* Compute the (r,s) entry of A.T   */
{
  colrange j;
  double temp;
  
  temp=0;
  for (j=0; j< n_size; j++) {
    temp = temp + A[r-1][j] * T[j][s-1];
  }
  return temp;
}

void dp_WriteTableau(FILE *f, rowrange m_size, colrange n_size, Amatrix A, Bmatrix T,
  colindex nbindex, rowindex bflag)
/* Write the tableau  A.T   */
{
  colrange j;
  rowrange i;
  
  fprintf(f, "  %ld   %ld    real\n",m_size, n_size);
  fprintf(f,"          |");
  for (j=1; j<= n_size; j++) {
    fprintf(f," %12ld", nbindex[j]);
  } fprintf(f,"\n");
  for (j=1; j<= n_size+1; j++) {
    fprintf(f," ------------");
  } fprintf(f,"\n");
  for (i=1; i<= m_size; i++) {
    fprintf(f," %3ld(%3ld) |", i, bflag[i]);  
    for (j=1; j<= n_size; j++) {
      fprintf(f," %12.3f",dp_TableauEntry(m_size, n_size, A,T,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
}


void SelectDualSimplexPivot(rowrange m_size, colrange n_size, 
    int Phase1, Amatrix A, Bmatrix T, rowindex OV, 
    colindex nbindex, rowindex bflag,
    rowrange objrow, colrange rhscol,
    rowrange *r, colrange *s, int *selected, dp_LPStatusType *lps)
{ /* selects a dual simplex pivot (*r, *s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=dp_FALSE and *lps=LPSundecided.
     If Phase1=dp_TRUE, the RHS column will be considered as the negative
     of the column of the largest variable (==m_size).  For this case, it is assumed
     that the caller used the auxiliary row (with variable m_size) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for m_size corresponds to the auxiliary variable.
  */
  int colselected=dp_FALSE, rowselected=dp_FALSE, 
    dualfeasible=dp_TRUE,localdebug=dp_FALSE;
  rowrange i;
  colrange j;
  double val=0, minval=0,rat=0, minrat=0;
  static long lastnn=0;
  static Arow rcost;

  lastnn=n_size;
  *r=0; *s=0;
  *selected=dp_FALSE;
  *lps=dp_LPSundecided;
  for (j=1; j<=n_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n", j, nbindex[j]); 
      rcost[j-1]=dp_TableauEntry(m_size, n_size, A,T,objrow,j);
      if (localdebug) printf("reduced cost =  %f\n", rcost[j-1]); 
      if (dp_Positive(rcost[j-1])) { 
        dualfeasible=dp_FALSE;
      }
    }
  }
  if (dualfeasible){
    while ((*lps==dp_LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=m_size; i++) {
        if (localdebug) printf("checking the row var %ld\n",i); 
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            val=-dp_TableauEntry(m_size, n_size, A,T,i,bflag[m_size]); /* for dual Phase I */
          } 
          else {val=dp_TableauEntry(m_size, n_size, A,T,i,rhscol);}
          if (localdebug) printf("RHS val =  %f\n", val); 
          if (val < minval) {
            *r=i;
            minval=val;
            if (localdebug) printf("update minval with = %f  *r = %ld\n",minval, *r);
          }
        }
      }
      if (dp_Nonnegative(minval)) {
        *lps=dp_Optimal;
        if (localdebug) printf("Select DualSimplexPivot: Optimal solution found.\n");
      }
      else {
        rowselected=dp_TRUE;
        for (j=1; j<=n_size; j++){
          val=dp_TableauEntry(m_size, n_size, A,T,*r,j);
          if (j!=rhscol && dp_Positive(val)) {
            rat=-rcost[j-1]/val;
            if (localdebug) printf("new ratio = %f at col %ld\n",rat, j);
            if (*s==0 || rat < minrat){
              minrat=rat;
              *s=j;
              if (localdebug) printf("update minrat = %f  *s = %ld\n",minrat, *s);
            }
          }
        }
        if (localdebug) printf("*s is %ld\n",*s);
        if (*s>0) {colselected=dp_TRUE; *selected=dp_TRUE;}
        else *lps=dp_Inconsistent;
      }
    } /* end of while */
  }
  if (localdebug) {
     if (Phase1) printf("Phase 1 : select %ld, %ld\n", *r, *s);
     else printf("Phase 2 : select %ld, %ld\n", *r, *s);
  }
}

void dp_SelectPivot2(rowrange m_size, colrange n_size, Amatrix A, Bmatrix T,
            dp_HyperplaneOrderType roworder, rowindex ordervec,
            rowrange rowmax, rowset NopivotRow,
            colset NopivotCol, rowrange *r, colrange *s,
            int *selected)
/* Select a position (*r,*s) in the matrix A.T such that (A.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  int stop;
  rowrange rtemp;
  rowset rowexcluded;
  double Xtemp;

  stop = dp_FALSE;
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=m_size;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = dp_FALSE;
  do {
    rtemp=0;
    dp_SelectPreorderedNext(m_size, n_size, rowexcluded, ordervec, &rtemp);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= n_size && !*selected) {
        Xtemp=dp_TableauEntry(m_size, n_size,A,T,*r,*s);
        if (!set_member(*s,NopivotCol) && dp_Nonzero(Xtemp)) {
          *selected = dp_TRUE;
          stop = dp_TRUE;
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
      stop = dp_TRUE;
    }
  } while (!stop);
  set_free(rowexcluded);
}

void dp_GaussianColumnPivot(rowrange m_size, colrange n_size, 
    Amatrix A, Bmatrix T, colindex nbindex, rowindex bflag, rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix A by
   updating the dual basis inverse  T.
 */
{
  int localdebug=dp_FALSE;
  long j, j1, entering;
  Arow Rtemp;
  double Xtemp0, Xtemp;

  if (localdebug) {
	 fprintf(stdout, "Column pivot: (leaving, entering) = (%ld, %ld)\n",  r, entering);
     dp_WriteBmatrix(stdout, n_size, T);
     dp_WriteTableau(stdout, m_size, n_size, A, T, nbindex, bflag);
  }

  for (j=1; j<=n_size; j++) Rtemp[j-1]=dp_TableauEntry(m_size, n_size, A, T, r,j);
  Xtemp0 = Rtemp[s-1];
  for (j = 1; j <= n_size; j++) {
    if (j != s) {
      Xtemp = Rtemp[j-1];
      for (j1 = 1; j1 <= n_size; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= n_size; j++) T[j-1][s - 1] /= Xtemp0;

  entering=nbindex[s];
  bflag[r]=s;     /* the nonbasic variable r corresponds to column s */
  nbindex[s]=r;   /* the nonbasic variable on s column is r */
  if (entering>0) bflag[entering]=-1;
     /* original variables have negative index and should not affect the row index */

}


void dp_InitializeBmatrix(colrange n_size, Bmatrix T)
{
  colrange j;
  double x;

  for (j = 0; j < n_size; j++) {
    T[j]=(double *)calloc(n_size, sizeof x);
  }
}

void dp_free_Bmatrix(colrange n_size, Bmatrix T)
{
  colrange j;

  for (j = 0; j < n_size; j++) {
    free(T[j]);
  }
}

void dp_SetToIdentity(colrange n_size, Bmatrix T)
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

void ResetTableau(rowrange m_size, colrange n_size, Bmatrix T,
    colindex nbindex, rowindex bflag, rowrange objrow, colrange rhscol,
   int UsePrevBasis)
{
  rowrange i;
  colrange j;
  
  if (!UsePrevBasis) {  /* Initialize T and nbindex */
    for (j=1; j<=n_size; j++) nbindex[j]=-j;
    nbindex[rhscol]=0; 
      /* RHS is already in nonbasis and is considered to be associated
         with the zero-th row of input. */
     dp_SetToIdentity(n_size, T);
  }
  
  /* Set the bflag according to nbindex */
  for (i=1; i<=m_size; i++) bflag[i]=-1;  
    /* all basic variables have index -1 */
  bflag[objrow]= 0; 
    /* bflag of the objective variable is 0, 
       different from other basic variables which have -1 */
  for (j=1; j<=n_size; j++) if (nbindex[j]>0) bflag[nbindex[j]]=j;
    /* bflag of a nonbasic variable is its column number */

}

void dp_SelectCrissCrossPivot(rowrange m_size, colrange n_size, Amatrix A, Bmatrix T,
    rowindex bflag, rowrange objrow, colrange rhscol,
    rowrange *r, colrange *s,
    int *selected, dp_LPStatusType *lps)
{
  int colselected=dp_FALSE, rowselected=dp_FALSE;
  rowrange i;
  double val;
  
  *selected=dp_FALSE;
  *lps=dp_LPSundecided;
  while ((*lps==dp_LPSundecided) && (!rowselected) && (!colselected)) {
    for (i=1; i<=m_size; i++) {
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=dp_TableauEntry(m_size, n_size, A,T,i,rhscol);
        if (dp_Negative(val)) {
          rowselected=dp_TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=dp_TableauEntry(m_size, n_size, A,T,objrow,bflag[i]);
        if (dp_Positive(val)) {
          colselected=dp_TRUE;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=dp_Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=m_size; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          val=dp_TableauEntry(m_size, n_size, A,T,*r,bflag[i]);
          if (dp_Positive(val)) {
            colselected=dp_TRUE;
            *s=bflag[i];
            *selected=dp_TRUE;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=m_size; i++) {
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=dp_TableauEntry(m_size, n_size, A,T,i,*s);
          if (dp_Negative(val)) {
            rowselected=dp_TRUE;
            *r=i;
            *selected=dp_TRUE;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=dp_DualInconsistent;
    }
    else if (!colselected) {
      *lps=dp_Inconsistent;
    }
  }
}

void CrissCrossMinimize(rowrange m_size, colrange n_size, 
    Amatrix A,Bmatrix T, 
    rowrange objrow, colrange rhscol, int UsePrevBasis, dp_LPStatusType *LPS,
    double *optvalue, Arow sol, Arow dsol, colindex nbindex,
    rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
{
   colrange j;

   *err=dp_None;
   for (j=1; j<=n_size; j++)
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   CrissCrossMaximize(m_size, n_size, A,T, objrow, rhscol, 
     UsePrevBasis, LPS, optvalue, sol, dsol, nbindex, re,  se, iter, err);
   *optvalue=-*optvalue;
   for (j=1; j<=n_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   }
}

void CrissCrossMaximize(rowrange m_size, colrange n_size,
    Amatrix A,Bmatrix T, 
    rowrange objrow, colrange rhscol, int UsePrevBasis, dp_LPStatusType *LPS,
    double *optvalue, Arow sol, Arow dsol, colindex nbindex,
    rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop, chosen, found;
  long rank, pivots_p0;
  rowrange i,r;
  colrange j,s;
  static rowindex bflag;
  static long mlast=0;
  static rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  int localdebug=dp_FALSE;
  unsigned int rseed=1;

  *err=dp_None;
  if (bflag==NULL || mlast!=m_size){
     if (mlast!=m_size) {
       free(bflag);   /* called previously with different m_size */
       free(OrderVector);
     }
     bflag=(long *) calloc(m_size+1, sizeof *bflag);
     OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector); 
     /* initialize only for the first time or when a larger space is needed */
     mlast=m_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size, n_size, A, OrderVector, dp_MinIndex, rseed);

  *re=0; *se=0; *iter=0;

  ResetTableau(m_size,n_size, T, nbindex, bflag, 
    objrow, rhscol, UsePrevBasis);

  FindLPBasis(m_size, n_size, A, T, 
      OrderVector, nbindex, bflag, 
      objrow, rhscol, UsePrevBasis, &s, &found, LPS, &pivots_p0);
  *iter+=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  stop=dp_FALSE;
  do {   /* Criss-Cross Method */
    dp_SelectCrissCrossPivot(m_size, n_size, A, T, bflag,
       objrow, rhscol, &r, &s, &chosen, LPS);
    if (chosen) {
      dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, r, s);
      (*iter)++;
    } else {
      switch (*LPS){
        case dp_Inconsistent: *re=r;
        case dp_DualInconsistent: *se=s;
        default: break;
      }
      stop=dp_TRUE;
    }
  } while(!stop);
  
_L99:

  SetSolutions(m_size, n_size, A, T, 
   objrow, rhscol, *LPS, optvalue, sol, dsol, nbindex, *re, *se);

}


int dp_LexSmaller(double *v1, double *v2, long nmax)
{ /* nmax is the size of vectors v1,v2 */
  int determined, smaller;
  colrange j;

  smaller = dp_FALSE;
  determined = dp_FALSE;
  j = 1;
  do {
    if (dp_Nonzero(v1[j - 1]-v2[j - 1])) {
      if (v1[j - 1] < v2[j - 1]) {
	    smaller = dp_TRUE;
	  }
      determined = dp_TRUE;
    } else
      j++;
  } while (!(determined) && (j <= nmax));
  return smaller;
}

int dp_LexLarger(double *v1, double *v2, long nmax)
{
  Arow u1, u2;
  colrange j;

  for (j = 1; j <= nmax; j++) {
    u1[j-1] = -v1[j-1];
    u2[j-1] = -v2[j-1];
  }
  return (dp_LexSmaller(u1, u2, nmax));
}

void FindLPBasis(rowrange m_size, colrange n_size,
    Amatrix A, Bmatrix T, rowindex OV, colindex nbindex, 
    rowindex bflag, rowrange objrow, colrange rhscol,int useprevbasis,
    colrange *cs, int *found, dp_LPStatusType *lps, long *pivot_no)
{ /* Find a LP basis using Gaussian pivots.
     If the problem has an LP basis,
     the procedure returns *found=dp_TRUE, *lps=LPSundecided and an LP basis.
     If the constraint matrix A (excluding the rhs and objective) is not
     column indepent, there are two cases.  If the dependency gives a dual
     inconsistency, this returns *found=dp_FALSE, *lps=dp_StrucDualInconsistent and 
     the evidence column *s.  Otherwise, this returns *found=dp_TRUE, 
     *lps=LPSundecided and an LP basis of size less than n_size.  Columns j
     that do not belong to the basis (i.e. cannot be chosen as pivot because
     they are all zero) will be indicated in nbindex vector: nbindex[j] will
     be negative and set to -j.
  */
  int localdebug=dp_FALSE,chosen,stop;
  long pivots_p0=0, rank;
  colset ColSelected;
  rowset RowSelected;
  double val;

  rowrange i,r;
  colrange j,s;
  static Arow rcost;

  *found=dp_FALSE; *cs=0; rank=0;
  *lps=dp_LPSundecided;

  set_initialize(&RowSelected,m_size);
  set_initialize(&ColSelected,n_size);
  set_addelem(RowSelected, objrow);
  set_addelem(ColSelected, rhscol);
 
  stop=dp_FALSE;
  do {   /* Find a LP basis */
    dp_SelectPivot2(m_size, n_size, A, T, dp_MinIndex, OV,
      m_size, RowSelected, ColSelected, &r, &s, &chosen);
    if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
    if (chosen) {
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
      dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, r, s);
      pivots_p0++;
      rank++;
    } else {
      for (j=1;j<=n_size  && *lps==dp_LPSundecided; j++) {
        if (j!=rhscol && nbindex[j]<0){
          if (localdebug) printf("col%ld  %ld\n", j, nbindex[j]);
          val=dp_TableauEntry(m_size, n_size, A,T,objrow,j);
          if (dp_Nonzero(val)){  /* dual inconsistent */
            *lps=dp_StrucDualInconsistent;
            *cs=j;
            if (localdebug)  
              printf("dual inconsistent because the nonzero reduced cost: %lf \n",val);
          }
        }
      }
      if (*lps==dp_LPSundecided) *found=dp_TRUE;  
         /* dependent columns but not dual inconsistent. */
      stop=dp_TRUE;
    }
    if (rank==n_size-1) {
      stop = dp_TRUE;
      *found=dp_TRUE;
    }
  } while (!stop);

/* Check whether the objrow is in the basis in case of UsePrevBasis. */
  if (useprevbasis && (s=bflag[objrow])>0){ /* objrow in the nonbasis. */
    for (j=0;j<=n_size;j++){
      if (j!=s) set_addelem(ColSelected,j);
      if (i=nbindex[j]>0) set_addelem(RowSelected,i);
    }
    if (localdebug) printf("UsePrevBasis but the current basis does not contain objrow.\n");
    dp_SelectPivot2(m_size, n_size, A, T, dp_MinIndex, OV,
      m_size, RowSelected, ColSelected, &r, &s, &chosen);
    if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
    if (chosen) {
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
      dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, r, s);
      pivots_p0++;
    } else {
      /* objrow is nonbasic and the corresponding column is zero. */
      *found=dp_FALSE;
      *lps=dp_StrucDualInconsistent;
      *cs=s;
    }
  }

  *pivot_no=pivots_p0;
  set_free(RowSelected);
  set_free(ColSelected);
}

void FindDualFeasibleBasis(rowrange m_size, colrange n_size,
    Amatrix A, Bmatrix T, rowindex OV, 
    colindex nbindex, rowindex bflag, rowrange objrow, colrange rhscol,
    colrange *s, int *found, dp_LPStatusType *lps, long *pivot_no)
{ /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *found=dp_TRUE, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *found=dp_FALSE, *lps=DualInconsistent and the evidence column *s.
  */
  int phase1, dualfeasible=dp_TRUE,localdebug=dp_FALSE,chosen,stop;
  dp_LPStatusType LPSphase1;
  long pivots_p1=0;
  rowrange i,rtemp;
  colrange j,l,mmsave,ms=0,stemp;
  double val=0,purezero=0,maxcost=-1;
  static long lastnn=0;
  static Arow rcost;

  *found=dp_TRUE; *lps=dp_LPSundecided; *s=0;
  mmsave=m_size;
  m_size=m_size+1;  /* increase m_size by 1 temporally  */
  A[m_size-1]= (double *) calloc(n_size, sizeof val);   /* create an auxiliary row  */
  lastnn=n_size;

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=n_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n", j,nbindex[j]); 
      rcost[j-1]=dp_TableauEntry(m_size, n_size, A,T,objrow,j);
      if (localdebug) printf("reduced cost =  %f\n",rcost[j-1]); 
      if (rcost[j-1] > maxcost) {maxcost=rcost[j-1]; ms = j;}
    }
  }
  if (dp_Positive(maxcost)) dualfeasible=dp_FALSE;

  if (!dualfeasible){
    for (j=1; j<=n_size; j++){
      A[m_size-1][j-1]=purezero;
      for (l=1; l<=n_size; l++){
        if (nbindex[l]>0) {
          A[m_size-1][j-1]-=A[nbindex[l]-1][j-1]; 
          /* To make the auxiliary row (0,-1,-1,...,-1).  */
        }
      }
    }
    if (localdebug){
      printf("Auxiliary row =");
      for (j=1; j<=n_size; j++){
        printf(" ( %ld):%f",j,dp_TableauEntry(m_size, n_size, A,T,m_size,j)); 
      }
      printf("\n");
    }

    if (localdebug){
      printf("FindDualFeasibleBasis: curruent basis is not dual feasible.\n");
      printf("because of the column %ld assoc. with var %ld   dual cost =%f\n",
       ms,nbindex[ms],maxcost);
    }

    /* Pivot on (m_size, ms) so that the dual basic solution becomes feasible */
    dp_GaussianColumnPivot(m_size, n_size, A,T, nbindex, bflag, m_size, ms);
    pivots_p1=pivots_p1+1;

    phase1=dp_TRUE; stop=dp_FALSE;
    do {   /* Dual Simplex Phase I */
      chosen=dp_FALSE; LPSphase1=dp_LPSundecided;
      SelectDualSimplexPivot(m_size, n_size, phase1, A, T, OV, nbindex, bflag,
        objrow, rhscol, &rtemp, &stemp, &chosen, &LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dp_TableauEntry(m_size, n_size, A,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but m_size is still in nonbasis.
           We must pivot in the auxiliary variable m_size. */

        double minval=0;
        rtemp=0;
        for (i=1; i<=m_size; i++){
          if (bflag[i]<0) { 
             /* i is basic and not the objective variable */
            val=dp_TableauEntry(m_size, n_size, A,T,i,ms);  /* auxiliary column*/
            if (val < minval) {
              rtemp=i;
              minval=val;
              if (localdebug) printf("update minval with = %f  rtemp = %ld\n",minval, rtemp);
            }
          }
        }

        dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, rtemp, ms);
        pivots_p1=pivots_p1+1;

        if (dp_Negative(dp_TableauEntry(m_size, n_size, A, T, objrow, ms))){
          if (localdebug){
            printf("Dual infeasible.\n");
            printf("obj-ms: %f  dp_zero = %f\n", 
              dp_TableauEntry(m_size, n_size, A,T,objrow,ms), dp_zero);
          }
          *found=dp_FALSE; *lps=dp_DualInconsistent;  *s=ms;
        }
        stop=dp_TRUE;
      } else {
        dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, rtemp, stemp);
        pivots_p1=pivots_p1+1;
        if (bflag[m_size]<0) {
          stop=dp_TRUE; 
          if (localdebug) printf("Dual Phase I: the auxiliary variable entered the basis, go to phase II\n");
        }
      }
    } while(!stop);
  }
  free(A[m_size-1]);
  m_size=mmsave;
  *pivot_no=pivots_p1;
}

void DualSimplexMinimize(rowrange m_size, colrange n_size,
   Amatrix A,Bmatrix T, 
   rowrange objrow, colrange rhscol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex nbindex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
{
   colrange j;

   *err=dp_None;
   for (j=1; j<=n_size; j++)
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   DualSimplexMaximize(m_size, n_size, A,T, objrow, rhscol, UsePrevBasis, 
     LPS, optvalue, sol, dsol, nbindex, re,  se, iter, err);
   *optvalue=-*optvalue;
   for (j=1; j<=n_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[objrow-1][j-1]=-A[objrow-1][j-1];
   }
}

void DualSimplexMaximize(rowrange m_size, colrange n_size,
   Amatrix A, Bmatrix T, 
   rowrange objrow, colrange rhscol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex nbindex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  double sw;
  int stop, chosen, phase1, found;
  long rank;
  long pivots_ds=0, pivots_p0=0, pivots_p1=0, pivots_pc=0, maxpivots, maxpivfactor=70;
  rowrange i,r;
  colrange j,s;
  static rowindex bflag;
  static long mlast=0,nlast=0;
  int localdebug=dp_FALSE;
  static rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  
  *err=dp_None;
  maxpivots=maxpivfactor*n_size;  /* maximum pivots to be performed before cc pivot is applied. */
  if (mlast!=m_size || nlast!=n_size){
     if (mlast>0) { /* called previously with different m_size */
       if (localdebug) printf("DualSimplex: deleting the old memory space with mlast = %ld\n",mlast);
       free(OrderVector);
       free(bflag);
     }
     if (localdebug) printf("DualSimplex: allocating a new memory space with m_size = %ld\n",m_size);
     OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);
     bflag=(long *) calloc(m_size+2, sizeof *bflag);  /* one more element for an auxiliary variable  */
     mlast=m_size;nlast=n_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size, n_size, A, OrderVector, dp_MinIndex, rseed);

  *re=0; *se=0; *iter=0;
  
  ResetTableau(m_size,n_size, T, nbindex, bflag, objrow, rhscol, UsePrevBasis);
   
  FindLPBasis(m_size, n_size, A, T, OrderVector, nbindex, bflag, 
      objrow, rhscol, UsePrevBasis, &s, &found, LPS, &pivots_p0);
  *iter+=pivots_p0;

  if (!found){
     *se=s;
     goto _L99;
     /* No LP basis is found, and thus Inconsistent.  
     Output the evidence column. */
  }

  FindDualFeasibleBasis(m_size, n_size, A, T, 
      OrderVector, nbindex, bflag, 
      objrow, rhscol, &s, &found, LPS, &pivots_p1);
  *iter+=pivots_p1;

  if (!found){
     *se=s;
     goto _L99;
     /* No dual feasible basis is found, and thus DualInconsistent.  
     Output the evidence column. */
  }
  
  /* Dual Simplex Method */
  stop=dp_FALSE;
  do {
    chosen=dp_FALSE; *LPS=dp_LPSundecided; phase1=dp_FALSE;
    if (pivots_ds<maxpivots) {
      SelectDualSimplexPivot(m_size, n_size, 
        phase1, A, T, OrderVector, nbindex, bflag,
        objrow, rhscol, &r, &s, &chosen, LPS);
    }
    if (chosen) pivots_ds=pivots_ds+1;
    if (!chosen && *LPS==dp_LPSundecided) {  
      /* In principle this should not be executed because we already have dual feasibility
         attained and dual simplex pivot should have been chosen.  This might occur
         under floating point computation, or the case of cycling.
      */
      dp_SelectCrissCrossPivot(m_size, n_size, A, T, bflag,
        objrow, rhscol, &r, &s, &chosen, LPS);
      if (chosen) pivots_pc=pivots_pc+1;
    }
    if (chosen) {
      dp_GaussianColumnPivot(m_size, n_size, A, T, nbindex, bflag, r, s);
      (*iter)++;
    } else {
      switch (*LPS){
        case dp_Inconsistent: *re=r;
        case dp_DualInconsistent: *se=s;
        default: break;
      }
      stop=dp_TRUE;
    }
  } while(!stop);

_L99:

  if (localdebug){
     printf("LP solved with %ld pivots. (ds pivt#= %ld,  p1 piv#= %ld",*iter,pivots_ds, pivots_p1);
     if (pivots_pc > 0) printf(", cc piv#= %ld",pivots_pc);
     printf(")\n");
  }
  
  SetSolutions(m_size, n_size, A, T, 
   objrow, rhscol, *LPS, optvalue, sol, dsol, nbindex, *re, *se);

}

void SetSolutions(rowrange m_size, colrange n_size,
   Amatrix A, Bmatrix T, 
   rowrange objrow, colrange rhscol, dp_LPStatusType LPS,
   double *optvalue, Arow sol, Arow dsol, colindex nbindex,
   rowrange re, colrange se)
/* 
Assign the solution vectors to sol, dsol, *optvalue after solving
the LP.
*/
{
  colrange j;
  double sw;
  int localdebug=dp_FALSE;
  
  switch (LPS){
  case dp_Optimal:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=T[j-1][rhscol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, T,objrow,j);
      *optvalue=dp_TableauEntry(m_size, n_size, A, T,objrow,rhscol);
      if (localdebug) printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    break;
  case dp_Inconsistent:
    if (localdebug) printf("DualSimplexSolve: LP is inconsistent.\n");
    for (j=1;j<=n_size; j++) {
      sol[j-1]=T[j-1][rhscol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, T, re,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    break;
  case dp_DualInconsistent:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=T[j-1][se-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, T,objrow,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  case dp_StrucDualInconsistent:
    if (dp_Positive(dp_TableauEntry(m_size, n_size, A, T,objrow, se))) sw=1;
    else sw=-1;
    for (j=1;j<=n_size; j++) {
      sol[j-1]=sw*T[j-1][se-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, T,objrow,j);
      if (localdebug)  printf("dsol[%ld]= %f\n",nbindex[j],dsol[j-1]);
    }
    if (localdebug) printf( "DualSimplexSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }
}


long dp_Partition(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  double *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  i=p-1;
  j=r+1;
  while (dp_TRUE){
    do{
      j--;
    } while (dp_LexLarger(A[OV[j]-1],x,nmax));
    do{
      i++;
    } while (dp_LexSmaller(A[OV[i]-1],x,nmax));
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

void dp_QuickSort(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  long q;
  
  if (p < r){
    q = dp_Partition(OV, p, r, A, nmax);
    dp_QuickSort(OV, p, q, A, nmax);
    dp_QuickSort(OV, q+1, r, A, nmax);
  }
}

void dp_LineShellingOrder(rowrange m_size, colrange n_size, Amatrix A, rowindex OV, double *z, double *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  double temp1,temp2,infinity=10.0e+20;
  static double *beta;
  static long mlast=0;
  int localdebug=dp_FALSE;
  
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
    if (dp_Nonzero(temp1)) A[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > 0) A[i-1][0]= infinity;
    else A[i-1][0]= -infinity;
     /* use the first column of A tentatively */
  }
  if (localdebug) 
    for (i=1; i<= m_size; i++){
      printf("set A[%ld] = %g\n", i, A[i-1][0]);
    }
  dp_QuickSort(OV, 1, m_size, A, 1);
  for (i=1; i<= m_size; i++) {
    A[i-1][0]=beta[i-1]; 
     /* restore the first column of A */ 
    if (localdebug) printf("restore A[%ld] with %g\n", i, A[i-1][0]);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void dp_RandomPermutation(rowindex OV, long t, unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  int localdebug=dp_FALSE;

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

void dp_ComputeRowOrderVector(rowrange m_size, colrange n_size, Amatrix A,
    rowindex OV, dp_HyperplaneOrderType ho, unsigned int rseed)
{
  long i,itemp,j;
  Arow zvec, dvec;
  
  OV[0]=0;
  switch (ho){
  case dp_MaxIndex:
    for(i=1; i<=m_size; i++) OV[i]=m_size-i+1;
    break;

  case dp_MinIndex: 
    for(i=1; i<=m_size; i++) OV[i]=i;
    break;

  case dp_LexMin:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dp_QuickSort(OV, 1, m_size, A, n_size);
    break;

  case dp_LexMax:
    for(i=1; i<=m_size; i++) OV[i]=i;
    dp_QuickSort(OV, 1, m_size, A, n_size);
    for(i=1; i<=m_size/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[m_size-i+1];
      OV[m_size-i+1]=itemp;
    }
    break;

  case dp_RandomRow:
    for(i=1; i<=m_size; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    dp_RandomPermutation(OV, m_size, rseed);
    break;

  case dp_LineShelling:
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
    dp_LineShellingOrder(m_size, n_size, A, OV, zvec, dvec);
    break;
  }
}


void dp_SelectPreorderedNext(rowrange m_size, colrange n_size, 
    rowset excluded, rowindex OV, rowrange *hnext)
{
  rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=m_size && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k, excluded)) *hnext=k ;
  }
}


void dp_LPSolve(dp_LPConversionType lpconv, dp_LPSolverType solver, 
   rowrange m_size, colrange n_size,
   Amatrix A, Bmatrix T, 
   rowrange objrow, colrange rhscol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex nbindex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  static Amatrix Acopy;
  rowrange i;

  for (i=1; i<=m_size; i++) Acopy[i-1]=A[i-1];

  *err=dp_None;
  switch (lpconv) {
    case dp_LPmax:
      if (solver==dp_CrissCross)
         CrissCrossMaximize(m_size, n_size, Acopy, T, objrow, rhscol, UsePrevBasis,
           LPS, optvalue, sol, dsol,nbindex, re, se, iter, err);
      else
         DualSimplexMaximize(m_size, n_size, Acopy, T, objrow, rhscol, UsePrevBasis,
           LPS, optvalue, sol, dsol,nbindex, re, se, iter, err);
      break;
      
    case dp_LPmin:
      if (solver==dp_CrissCross)
         CrissCrossMinimize(m_size, n_size, Acopy, T, objrow, rhscol, UsePrevBasis,
           LPS, optvalue, sol, dsol,nbindex, re, se, iter, err);
      else
         DualSimplexMinimize(m_size, n_size, Acopy, T, objrow, rhscol, UsePrevBasis,
           LPS, optvalue, sol, dsol,nbindex, re, se, iter, err);
      break;
  }
}

void MakeAforInteriorFinding(rowrange m_size, colrange n_size, Amatrix A, 
  rowrange OBJrow, colrange RHScol,
  rowrange *new_m_size, colrange *new_n_size, Amatrix new_A, rowrange *new_OBJrow)
/* Delete the objective row,
   add an extra column with -1's to the matrix A,
   add an extra row with (bceil, 0,...,0,-1),
   add an objective row with (0,...,0,1), and 
   rows & columns, and change m_size and n_size accordingly, to output new_A.
  This sets up the LP:
  maximize      x_{d+1}
  s.t.    A x + x_{d+1}  <=  b
                x_{d+1}  <=  bm * bmax,
  where bm is set to 2 by default, and bmax=max{1, b[1],...,b[m_size]}.
*/
{
  rowrange i; colrange j;
  double x;
  double bm=2.0, bmax=1, bceil;
  int localdebug=dp_FALSE;

  *new_m_size=m_size+1;
  *new_n_size=n_size+1;
  *new_OBJrow=OBJrow;
  for (i=1; i<=m_size; i++) {
    if (A[i-1][RHScol-1]>bmax) bmax = A[i-1][RHScol-1];
  }
  bceil=bm*bmax;
  if (localdebug) printf("bceil is set to %lg\n",bceil);
  
  for (i=1; i<=*new_m_size; i++){
    new_A[i-1]=(double *)calloc(*new_n_size,sizeof x);
  }
  for (i=1; i <= m_size; i++) {
    for (j=1; j <= n_size; j++) {
      new_A[i-1][j-1]=A[i-1][j-1];
      if (localdebug) dp_WriteReal(stdout, new_A[i-1][j-1]);
    }
    if (localdebug) fprintf(stdout,"\n");
  }
  for (i=1;i<= *new_m_size;i++) {
    if (i!=OBJrow) new_A[i-1][n_size]=-1.0;  /* new column with all minus one's */
  }
  for (j=1;j<= n_size;j++) {
    new_A[*new_m_size-1][j-1]=0.0;  /* new row (bceil, 0,...,0,-1) */
  }
  new_A[*new_m_size-1][RHScol-1]=bceil;    /* new row (bceil, 0,...,0,-1) */
  for (j=1;j<= n_size;j++) {
    new_A[OBJrow-1][j-1]=0.0;  /* new obj row with (0,...,0,1) */
  }
  new_A[OBJrow-1][n_size]=1.0; /* new obj row with (0,...,0,1) */
}

void dp_FindInteriorPoint(dp_LPSolverType solver, 
   rowrange m_size, colrange n_size,
   Amatrix A, rowrange objrow, colrange rhscol,  dp_LPStatusType *LPS,
   double *optvalue, Arow sol, long *iter, dp_ErrorType *err)
/* 
  This solvs the LP:
  maximize      x_{d+1}
  s.t.    A x + x_{d+1}  <=  b
                x_{d+1}  <=  bm * bmax,
  where bm is set to 2 and bmax=max{1, b[0]...,b[m_size-1]}.
  Thus, the optimum value is zero     if the polyhedron has no interior point.
        the optimum value is negative if the polyhedron is empty.
        the optimum value is positive if the polyhedron admits an interior point.
  For the last case, sol will be set to the interior point found.
*/
{
  dp_LPConversionType lpconv;
  rowrange m_new, re;
  colrange n_new, se;
  colindex nbindex;
  Arow dsol;
  Amatrix A_new;
  Bmatrix T;
  int UsePrevBasis=0;
  int localdebug=dp_FALSE;

  lpconv=dp_LPmax;
  MakeAforInteriorFinding(m_size, n_size, A, objrow, rhscol, &m_new, &n_new, A_new, &objrow);
  dp_InitializeBmatrix(n_new, T);

  *err=dp_None;
  dp_LPSolve(lpconv, solver, m_new, n_new, A_new, T, objrow, rhscol, UsePrevBasis, 
     LPS, optvalue, sol, dsol, nbindex, &re, &se, iter, err);

  if (localdebug) dp_WriteLPResult(stdout, lpconv, solver, m_new, n_new, A_new, objrow, rhscol,
    *LPS, *optvalue, sol, dsol, nbindex, re, se, *iter, *err);
}

void dp_WriteLPResult(FILE *f, dp_LPConversionType Conversion, dp_LPSolverType LPSolver,
    rowrange m_size, colrange n_size, 
    Amatrix A, rowrange objrow, colrange rhscol,
    dp_LPStatusType LPS, double optval,
    Arow sol, Arow dsol, colindex nbindex, rowrange re, colrange se,
    long iter, dp_ErrorType err)
{
  long j;

  fprintf(f,"\n*dplex LP result\n");
  
  if (err!=dp_None) {
    dp_WriteErrorMessages(f, err);
    goto _L99;
  }

  fprintf(f,"* #constraints = %ld\n", m_size-1);
  fprintf(f,"* #variables   = %ld\n", n_size-1);

  switch (LPSolver) {
    case dp_DualSimplex:
      fprintf(f,"*Algorithm: dual simplex algorithm\n");break; 
    case dp_CrissCross:
      fprintf(f,"*Algorithm: criss-cross method\n");break;
  }

  switch (Conversion) {
    case dp_LPmax:
      fprintf(f,"*maximization is chosen\n");break; 
    case dp_LPmin:
      fprintf(f,"*minimization is chosen\n");break;
  }
  
  if (Conversion==dp_LPmax||Conversion==dp_LPmin){
    fprintf(f,"*Objective function is\n");  
    for (j=0; j<n_size; j++){
      if (j>0 && A[objrow-1][j]>=0 ) fprintf(f," +");
      if (j>0 && (j % 5) == 0) fprintf(f,"\n");
      dp_WriteReal(f, A[objrow-1][j]);
      if (j>0) fprintf(f," X[%3ld]",j);
    }
    fprintf(f,"\n");
  }

  switch (LPS){
  case dp_Optimal:
    fprintf(f,"*LP status: a dual pair (x, y) of optimal solutions found.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_solution\n");
    for (j=1; j<n_size; j++) {
      fprintf(f,"  %3ld : ",j);
      dp_WriteReal(f,sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"  dual_solution\n");
    for (j=1; j<n_size; j++){
      if (nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",nbindex[j+1]);
        dp_WriteReal(f,dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"  optimal_value : % .9E\n", optval);
    fprintf(f,"end\n");
    break;

  case dp_Inconsistent:
    fprintf(f,"*LP status: LP is inconsistent.\n");
    fprintf(f,"*The positive combination of original inequalities with\n");
    fprintf(f,"*the following coefficients will prove the inconsistency.\n");
    fprintf(f,"begin\n");
    fprintf(f,"  dual_direction\n");
    fprintf(f,"  %3ld : ",re);
    dp_WriteReal(f,1.0);  fprintf(f,"\n");
    for (j=1; j<n_size; j++){
      if (nbindex[j+1]>0) {
        fprintf(f,"  %3ld : ",nbindex[j+1]);
        dp_WriteReal(f,dsol[j]); fprintf(f,"\n");
      }
    }
    fprintf(f,"end\n");
    break;

  case dp_DualInconsistent: case dp_StrucDualInconsistent:
    fprintf(f,"*LP status: LP is dual inconsistent.\n");
    fprintf(f,"*The linear combination of columns with\n");
    fprintf(f,"*the following coefficients will prove the dual inconsistency.\n");
    fprintf(f,"*(It is also an unbounded direction for the primal LP.)\n");
    fprintf(f,"begin\n");
    fprintf(f,"  primal_direction\n");
    for (j=1; j<n_size; j++) {
      fprintf(f,"  %3ld : ",j);
      dp_WriteReal(f,sol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  default:
    break;
  }
  fprintf(f,"*number of pivot operations = %ld\n", iter);
_L99:;
}

void dp_WriteBmatrix(FILE *f, colrange n_size, Bmatrix T)
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

void dp_WriteReal(FILE *f, double x)
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


/* end of dplex.c */


