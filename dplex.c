/* dplex.c:  dual simplex method c-code
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.60, August 21, 1996
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

#define COPYRIGHT   "Copyright (C) 1996, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DPLEXVERSION   "Version 0.60 (August 21, 1996)"

#define dp_FALSE 0
#define dp_TRUE 1

typedef enum {
  dp_MaxIndex, dp_MinIndex, dp_LexMin, dp_LexMax, dp_RandomRow, dp_LineShelling
} dp_HyperplaneOrderType;

void dp_LPInput(FILE **f, dp_FilenameType, rowrange *m, colrange *n, Amatrix A, 
    dp_LPConversionType *Conversion, rowrange *objrow, colrange *rhscol, 
    dp_ErrorType *err);
void dp_InitializeBmatrix(colrange n_size, Bmatrix T);
void CrissCrossMinimize(rowrange, colrange, Amatrix, Bmatrix BasisInverse,
  rowrange, colrange, 
  dp_LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *, dp_ErrorType *);
void CrissCrossMaximize(rowrange, colrange, Amatrix,Bmatrix BasisInverse,
  rowrange, colrange, 
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
void dp_WriteBmatrix(FILE *f, colrange n_size, Bmatrix T);
void dp_SetNumberType(char *line, dp_NumberType *number, dp_ErrorType *Error);
void dp_ComputeRowOrderVector(rowrange m_size, colrange n_size, Amatrix A,
    rowindex OV, dp_HyperplaneOrderType ho, unsigned int rseed);
void dp_SelectPreorderedNext(rowrange m_size, colrange n_size, 
    rowset excluded, rowindex OV, rowrange *hnext);
void dp_WriteReal(FILE *f, double x);
void dp_SetInputFile(FILE **f, dp_FilenameType inputfile,  dp_ErrorType *);

void dp_SetInputFile(FILE **f, dp_FilenameType inputfile, dp_ErrorType *Error)
{
  int opened=0,stop, quit=0;
  int i,dotpos=0, semipos=0,trial=0;
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
  long value1,value2;
  int fileopened=0, found=0, localdebug=0;
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
      A[i-1][j - 1] = value;
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


double dp_TableauEntry(rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
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

void dp_WriteTableau(FILE *f, rowrange m_size, colrange n_size, Amatrix X, Bmatrix T)
/* Write the tableau  X.T   */
{
  colrange j;
  rowrange i;
  
  fprintf(f, "begin\n");
  fprintf(f, "  %ld   %ld    real\n",m_size, n_size);
  for (i=1; i<= m_size; i++) {
    for (j=1; j<= n_size; j++) {
      fprintf(f," %12.6f",dp_TableauEntry(m_size, n_size, X,T,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"end\n");
}


void SelectDualSimplexPivot(rowrange m_size, colrange n_size, 
    int Phase1, Amatrix X, Bmatrix T, rowindex OV, 
    colindex NBIndex, long bflag[], rowrange objrow, colrange rhscol,
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
  int colselected=dp_FALSE, rowselected=dp_FALSE, dualfeasible=dp_TRUE,localdebug=dp_FALSE;
  rowrange i,k;
  colrange j;
  double val=0, minval=0,rat=0, minrat=0;
  static long lastnn=0;
  static double purezero=0;
  static Arow rcost;

  lastnn=n_size;
  *r=0; *s=0;
  *selected=dp_FALSE;
  *lps=dp_LPSundecided;
  for (j=1; j<=n_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n", j, NBIndex[j]); 
      rcost[j-1]=dp_TableauEntry(m_size, n_size, X,T,objrow,j);
      if (localdebug) printf("reduced cost =  %lf\n", rcost[j-1]); 
      if (rcost[j-1] > dp_zero) {
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
            val=-dp_TableauEntry(m_size, n_size, X,T,i,bflag[m_size]); /* for dual Phase I */
          } 
          else {val=dp_TableauEntry(m_size, n_size, X,T,i,rhscol);}
          if (localdebug) printf("RHS val =  %lf\n", val); 
          if (val < minval) {
            *r=i;
            minval=val;
            if (localdebug) printf("update minval with = %lf  *r = \n",minval, *r);
          }
        }
      }
      if (minval>=-dp_zero) {
        *lps=dp_Optimal;
        if (localdebug) printf("Select DualSimplexPivot: Optimal solution found.\n");
      }
      else {
        rowselected=dp_TRUE;
        for (j=1; j<=n_size; j++){
          val=dp_TableauEntry(m_size, n_size, X,T,*r,j);
          if (j!=rhscol && val > dp_zero) {
            rat=-rcost[j-1]/val;
            if (localdebug) printf("new ratio = %lf at col %ld\n",rat, j);
            if (*s==0 || rat < minrat){
              minrat=rat;
              *s=j;
              if (localdebug) printf("update minrat = %lf  *s = \n",minrat, *s);
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

void dp_SelectPivot2(rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
            dp_HyperplaneOrderType roworder, rowindex ordervec,
            rowrange rowmax, rowset NopivotRow,
            colset NopivotCol, rowrange *r, colrange *s,
            int *selected)
/* Select a position (*r,*s) in the matrix X.T such that (X.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  int stop;
  rowrange i,rtemp;
  rowset rowexcluded;
  double Xtemp;
  int localdebug=dp_FALSE;

  stop = dp_FALSE;
  set_initialize(&rowexcluded,m_size);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=m_size;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = dp_FALSE;
  do {
    rtemp=0;
/*    dp_SelectNextHyperplane(m_size, n_size, X,
      roworder, rowexcluded, &rtemp, ordervec);
*/
    dp_SelectPreorderedNext(m_size, n_size, rowexcluded, ordervec, &rtemp);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= n_size && !*selected) {
        Xtemp=dp_TableauEntry(m_size, n_size,X,T,*r,*s);
        if (!set_member(*s,NopivotCol) && fabs(Xtemp) > dp_zero) {
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

void dp_GausianColumnPivot2(rowrange m_size, colrange n_size, 
    Amatrix X, Bmatrix T, rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  long j, j1;
  Arow Rtemp;
  double Xtemp0, Xtemp;

  for (j=1; j<=n_size; j++) Rtemp[j-1]=dp_TableauEntry(m_size, n_size, X, T, r,j);
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


void dp_DualizeA(rowrange *m_size, colrange *n_size, Amatrix A, Bmatrix T, 
    dp_ErrorType dp_Error)
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
  if (mnew > dp_MMAX) {
    printf("MMAX  is too small for ray computation. MMAX must be >= %ld.\n",mnew);
    dp_Error = dp_DimensionTooLarge;
    goto _L99;
  }
  if (nnew > dp_NMAX) {
    printf("NMAX  is too small for ray computation. NMAX must be >= %ld.\n",nnew);
    dp_Error = dp_DimensionTooLarge;
    goto _L99;
  }
  for (i=1;i<= *m_size;i++){
    for (j=1;j<= *n_size;j++){
      Acopy[i-1][j-1]=dp_TableauEntry(*m_size, *n_size, A,T,i,j);
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

void dp_EnlargeAforInteriorFinding(rowrange *m_size, colrange *n_size, Amatrix A)
/* Add an extra column with all minus ones to the matrix AA, 
   add an objective row with (0,...,0,1), and 
   rows & columns, and change m_size and n_size accordingly 
*/
{
  long i,j,mnew=0,nnew=0;
  Amatrix Acopy;
  double x;
  int localdebug=0;

  mnew=*m_size+1;
  nnew=*n_size+1;
  for (i=0; i<mnew; i++){
    Acopy[i]=(double *)calloc(nnew,sizeof x);
  }
  for (i=1; i <= *m_size; i++) {
    for (j=1; j <= *n_size; j++) {
      Acopy[i-1][j-1]=A[i-1][j-1];
      if (localdebug) dp_WriteReal(stdout, A[i-1][j-1]);
    }
    if (localdebug) fprintf(stdout,"\n");
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

void dp_ComputeRank(rowrange m_size, colrange n_size, 
    Amatrix A1, rowset TargetRows, rowindex ordervec, long *rank)
/* Compute the rank of the submatrix of a Amatrix A1 indexed by TargetRows.
   This procedure does not change the matrix A1.
 */
{
  int stop, chosen, localdebug=0;
  rowrange r;
  colrange s;
  rowset NoPivotRow;
  colset ColSelected;
  Bmatrix Btemp;   /* dual basis inverse */
  
  *rank = 0;
  stop = dp_FALSE;
  set_initialize(&NoPivotRow, m_size);
  set_compl(NoPivotRow,TargetRows);
  set_initialize(&ColSelected, n_size);
  dp_InitializeBmatrix(n_size, Btemp);
  dp_SetToIdentity(n_size, Btemp);
  if (localdebug) dp_WriteBmatrix(stdout, n_size, Btemp);
  do {   /* Find a set of rows for a basis */
      dp_SelectPivot2(m_size, n_size, A1, Btemp, dp_MinIndex, ordervec, 
        m_size, NoPivotRow, ColSelected, &r, &s, &chosen);
      if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(NoPivotRow, r);
        set_addelem(ColSelected, s);
        (*rank)++;
        dp_GausianColumnPivot2(m_size, n_size, A1,Btemp, r, s);
        if (localdebug) {
          dp_WriteBmatrix(stdout, n_size, Btemp);
	      printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	    }
      } else {
        stop=dp_TRUE;
      }
  } while (!stop);
  set_free(NoPivotRow);
  set_free(ColSelected);
  dp_free_Bmatrix(n_size, Btemp);
}

void dp_ComputeBInverse(rowrange m_size, colrange n_size, Amatrix A1, long lastrow,
       rowindex ordervec, Bmatrix InvA1, long *rank)
{
  int stop, chosen, localdebug=0;
  rowrange r;
  colrange s;
  rowset RowSelected;
  colset ColSelected;

  *rank = 0;
  stop = dp_FALSE;
  dp_SetToIdentity(n_size, InvA1);
  set_initialize(&RowSelected, m_size);
  set_initialize(&ColSelected, n_size);
  do {
    dp_SelectPivot2(m_size, n_size, A1, InvA1, dp_MinIndex, ordervec,
      lastrow, RowSelected, ColSelected, &r, &s, &chosen);
    if (chosen) {
      (*rank)++;
      if (localdebug)
        printf("%3ldth pivot on%3ld, %3ld\n", *rank, r, s);
      dp_GausianColumnPivot2(m_size, n_size, A1, InvA1, r, s);
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
    } else
      stop = dp_TRUE;
  } while (!stop);
  set_free(RowSelected);
  set_free(ColSelected);
}



void dp_SelectCrissCrossPivot(rowrange m_size, colrange n_size, Amatrix X, Bmatrix T,
    long bflag[], rowrange objrow, colrange rhscol,
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
      if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=dp_TableauEntry(m_size, n_size, X,T,i,rhscol);
        if (val < -dp_zero) {
          rowselected=dp_TRUE;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=dp_TableauEntry(m_size, n_size, X,T,objrow,bflag[i]);
        if (val > dp_zero) {
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
          val=dp_TableauEntry(m_size, n_size, X,T,*r,bflag[i]);
          if (val > dp_zero) {
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
        if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=dp_TableauEntry(m_size, n_size, X,T,i,*s);
          if (val < -dp_zero) {
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
    Amatrix A,Bmatrix BasisInverse, 
    rowrange OBJrow, colrange RHScol, dp_LPStatusType *LPS,
    double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
    rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
{
   colrange j;

   *err=dp_None;
   for (j=1; j<=n_size; j++)
     A[OBJrow-1][j-1]=-A[OBJrow-1][j-1];
   CrissCrossMaximize(m_size, n_size, A,BasisInverse, OBJrow, RHScol, 
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter, err);
   *optvalue=-*optvalue;
   for (j=1; j<=n_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[OBJrow-1][j-1]=-A[OBJrow-1][j-1];
   }
}

void CrissCrossMaximize(rowrange m_size, colrange n_size,
    Amatrix A,Bmatrix BasisInverse, 
    rowrange OBJrow, colrange RHScol, dp_LPStatusType *LPS,
    double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
    rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop, chosen;
  long rank;
  rowrange i,r,entering,leaving;
  colrange j,s;
  colset ColSelected;
  rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0;
  static rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  double adet; /* abs value of the determinant of a basis */
  int localdebug=dp_FALSE,Q_debug=dp_FALSE;
  unsigned int rseed=1;

   *err=dp_None;
  if (BasisFlag==NULL || mlast!=m_size){
     if (mlast!=m_size) {
       free(BasisFlag);   /* called previously with different m_size */
       free(OrderVector);
     }
     BasisFlag=(long *) calloc(m_size+1, sizeof *BasisFlag);
     OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector); 
     /* initialize only for the first time or when a larger space is needed */
     mlast=m_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size, n_size, A, OrderVector, dp_MinIndex, rseed);

  *re=0; *se=0; *iter=0;
  rank = 0;
  stop = dp_FALSE;
  adet=1.0;
  set_initialize(&Cobasis,m_size);
  set_initialize(&Basis,m_size);
  set_initialize(&RowSelected, m_size);
  set_initialize(&ColSelected, n_size);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  for (i=0; i<=m_size; i++) BasisFlag[i]=0;
  for (j=0; j<=dp_NMAX; j++) NBIndex[j]=0;
  for (i=1; i<=m_size; i++) {
    set_addelem(Basis,i);
    BasisFlag[i]=-1;    /* basic variable has index -1 */
  }
  BasisFlag[OBJrow]= 0;
  dp_SetToIdentity(n_size, BasisInverse);
  if (localdebug) dp_WriteBmatrix(stdout, n_size, BasisInverse);
  do {   /* Find a LP basis */
      dp_SelectPivot2(m_size, n_size, A, BasisInverse, dp_MinIndex, OrderVector,
        m_size, RowSelected, ColSelected, &r, &s, &chosen);
      if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        set_addelem(Cobasis, r);
        set_delelem(Basis,r);
        if (localdebug) {
          fprintf(stdout, "CC: find initial basis: nonbasis = ");   
          set_fwrite(stdout,Cobasis); fprintf(stdout,"\n");
          fprintf(stdout, "CC: find initial basis: basis = ");
          set_fwrite(stdout,Basis);fprintf(stdout,"\n");
        }
        BasisFlag[r]=s;   /* the nonbasic variable r corresponds to column s */
        NBIndex[s]=r;     /* the nonbasic variable on s column is r */
        if (localdebug) fprintf(stdout,"nonbasic variable %ld has index %ld\n", r, BasisFlag[r]);
        rank++;
        dp_GausianColumnPivot2(m_size, n_size, A, BasisInverse, r, s);
        if (localdebug) {
          dp_WriteBmatrix(stdout, n_size, BasisInverse);
          dp_WriteTableau(stdout, m_size, n_size, A, BasisInverse),
	      fprintf(stdout, "%3ldth row added to the initial set (%ldth elem)\n",  r, rank);
	    }
      } else {
        stop=dp_TRUE;
      }
      if (rank==n_size-1) stop = dp_TRUE;
  } while (!stop);
  
  stop=dp_FALSE;
  do {   /* Criss-Cross Method */
    dp_SelectCrissCrossPivot(m_size, n_size, A, BasisInverse, BasisFlag,
       OBJrow, RHScol, &r, &s, &chosen, LPS);
    if (localdebug && chosen) fprintf(stdout, "Procedure Criss-Cross: pivot on (r,s) =(%ld, %ld).\n", r, s);
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
      if (localdebug) {
        fprintf(stdout, "nonbasis = "); set_write(Cobasis);
        fprintf(stdout, "basis = "); set_write(Basis);
        fprintf(stdout, "new nonbasic variable %ld has index %ld\n", leaving, BasisFlag[leaving]);
      }
      dp_GausianColumnPivot2(m_size, n_size, A, BasisInverse, r, s);
      (*iter)++;
      if (localdebug) {
        dp_WriteBmatrix(stdout, n_size, BasisInverse);
        dp_WriteTableau(stdout, m_size, n_size, A, BasisInverse);
	    fprintf(stdout, "%3ldth row added to the initial set (%ldth elem)\n",  r, rank);
	  }
    } else {
      switch (*LPS){
        case dp_Inconsistent: *re=r;
        case dp_DualInconsistent: *se=s;
        default: break;
      }
      stop=dp_TRUE;
    }
  } while(!stop);
  switch (*LPS){
  case dp_Optimal:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,j);
      *optvalue=dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,RHScol);
      if (localdebug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    break;
  case dp_Inconsistent:
    if (localdebug) printf("CrissCrossSolve: LP is inconsistent.\n");
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,*re,j);
      if (localdebug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    break;
  case dp_DualInconsistent:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][*se-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,j);
      if (localdebug) printf("dsol %ld  %5.2f \n",NBIndex[j],dsol[j-1]);
    }
    if (localdebug) printf("CrissCrossSolve: LP is dual inconsistent.\n");
    break;

  default:break;
  }

  set_free(ColSelected);
  set_free(RowSelected);
  set_free(Basis);
  set_free(Cobasis);
}


int dp_LexSmaller(double *v1, double *v2, long nmax)
{ /* nmax is the size of vectors v1,v2 */
  int determined, smaller;
  colrange j;

  smaller = dp_FALSE;
  determined = dp_FALSE;
  j = 1;
  do {
    if (fabs(v1[j - 1]-v2[j - 1])>dp_zero) {
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

void FindDualFeasibleBasis(rowrange m_size, colrange n_size,
    Amatrix X, Bmatrix T, rowset basis, rowset cobasis, rowindex OV, 
    colindex NBIndex, long bflag[], rowrange objrow, colrange rhscol,
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
  rowrange i,k,rtemp,entering,leaving,msvar;
  colrange j,l,mmsave,ms=0,stemp;
  double val=0,purezero=0,maxcost=-1;
  static long lastnn=0;
  static Arow rcost;

  *found=dp_TRUE; *lps=dp_LPSundecided; *s=0;
  mmsave=m_size;
  m_size=m_size+1;  /* increase m_size by 1 temporally  */
  if (lastnn != n_size){
    if (lastnn>0){
      free(X[m_size-1]);
    }
    X[m_size-1]= (double *) calloc(n_size, sizeof val);   /* create an auxiliary row  */
    lastnn=n_size;
  }

  ms=0;  /* ms will be the index of column which has the largest reduced cost */
  for (j=1; j<=n_size; j++){
    if (j!=rhscol){
      if (localdebug) printf("checking the column %ld var %ld\n", j,NBIndex[j]); 
      rcost[j-1]=dp_TableauEntry(m_size, n_size, X,T,objrow,j);
      if (localdebug) printf("reduced cost =  %lf\n",rcost[j-1]); 
      if (rcost[j-1] > maxcost) {maxcost=rcost[j-1]; ms = j; msvar=NBIndex[ms];}
    }
  }
  if (maxcost > dp_zero) dualfeasible=dp_FALSE;

  if (!dualfeasible){
    for (j=1; j<=n_size; j++){
      X[m_size-1][j-1]=purezero;
      for (l=1; l<=n_size; l++){
        if (NBIndex[l]>0) {
          X[m_size-1][j-1]-=X[NBIndex[l]-1][j-1]; 
          /* To make the auxiliary row (0,-1,-1,...,-1).  */
        }
      }
    }
    if (localdebug){
      printf("Auxiliary row =");
      for (j=1; j<=n_size; j++){
        printf(" ( %ld):%lf",j,dp_TableauEntry(m_size, n_size, X,T,m_size,j)); 
      }
      printf("\n");
    }

    if (localdebug){
      printf("FindDualFeasibleBasis: curruent basis is not dual feasible.\n");
      printf("because of the column %ld assoc. with var %ld   dual cost =%lf\n",
       ms,NBIndex[ms],maxcost);
    }

    /* Pivot on (m_size, ms) so that the dual basic solution becomes feasible */
    leaving=m_size; entering=msvar;
    set_addelem(cobasis, leaving);
    set_delelem(cobasis, entering);
    set_delelem(basis,leaving);
    set_addelem(basis,entering);
    bflag[leaving]=ms;
    bflag[entering]=-1;
    NBIndex[ms]=leaving;
    if (localdebug) {
      printf("Phase I ini: nonbasis = "); set_write(cobasis);
      printf("\nPhase I ini: basis = "); set_write(basis);
      printf("\nPhase I ini: new nonbasic variable %ld  has index %ld\n",leaving, bflag[leaving]);
    }
    dp_GausianColumnPivot2(m_size, n_size, X,T, m_size, ms);
    pivots_p1=pivots_p1+1;

    phase1=dp_TRUE; stop=dp_FALSE;
    do {   /* Dual Simplex Phase I */
      chosen=dp_FALSE; LPSphase1=dp_LPSundecided;
      SelectDualSimplexPivot(m_size, n_size, phase1, X, T, OV, NBIndex, bflag,
        objrow, rhscol, &rtemp, &stemp, &chosen, &LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           dp_TableauEntry(m_size, n_size, X,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but m_size is still in nonbasis.
           We must pivot in the auxiliary variable m_size. */

        double minval=0;
        rtemp=0;
        for (i=1; i<=m_size; i++){
          if (bflag[i]<0) { 
             /* i is basic and not the objective variable */
            val=dp_TableauEntry(m_size, n_size, X,T,i,ms);  /* auxiliary column*/
            if (val < minval) {
              rtemp=i;
              minval=val;
              if (localdebug) printf("update minval with = %lf  rtemp = %ld\n",minval, rtemp);
            }
          }
        }

        leaving=rtemp; entering=m_size;
        set_addelem(cobasis, leaving);
        set_delelem(cobasis, entering);
        set_delelem(basis,leaving);
        set_addelem(basis,entering);
        bflag[leaving]=ms;
        bflag[entering]=-1;
        NBIndex[ms]=leaving;
        if (localdebug) {
          printf("Phase I fin: nonbasis = "); set_write(cobasis);
          printf("\nPhase I fin: basis = "); set_write(basis);
          printf("\nPhase I ini: new nonbasic variable %ld  has index %ld\n",leaving, bflag[leaving]);
        }
        dp_GausianColumnPivot2(m_size, n_size, X, T, rtemp, ms);
        pivots_p1=pivots_p1+1;

        if (dp_TableauEntry(m_size, n_size, X, T, objrow, ms)<-dp_zero){
          if (localdebug){
            printf("Dual infeasible.\n");
            printf("obj-ms: %lf  dp_zero = %lf\n", 
              dp_TableauEntry(m_size, n_size, X,T,objrow,ms), dp_zero);
          }
          *found=dp_FALSE; *lps=dp_DualInconsistent;  *s=ms;
        }
        stop=dp_TRUE;
      } else {
        entering=NBIndex[stemp];
        leaving=rtemp;
        if (localdebug) {
          printf("Dual Phase I: pivot on (r,s) = (%ld, %ld)\n", 
            rtemp,stemp);
          printf("Dual Phase I: (leaving, entering) = (%ld, %ld)\n", 
            leaving,entering);
        }
        set_addelem(cobasis, leaving);
        set_delelem(cobasis, entering);
        set_delelem(basis,leaving);
        set_addelem(basis,entering);
        bflag[leaving]=stemp;
        bflag[entering]=-1;
        NBIndex[stemp]=leaving;
        if (localdebug) {
          printf("nonbasis = "); set_write(cobasis);
          printf("\nbasis = "); set_write(basis);
          printf("\nnew nonbasic variable %ld  has index %ld\n",leaving, bflag[leaving]);
        }
        dp_GausianColumnPivot2(m_size, n_size, X, T, rtemp, stemp);
        pivots_p1=pivots_p1+1;
        if (entering==m_size) {
          stop=dp_TRUE; 
          if (localdebug) printf("Dual Phase I: the auxiliary variable enter the basis, go to phase II\n");
        }
      }
    } while(!stop);
  }
  m_size=mmsave;
  *pivot_no=pivots_p1;
}

void DualSimplexMinimize(rowrange m_size, colrange n_size,
   Amatrix A,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
{
   colrange j;

   *err=dp_None;
   for (j=1; j<=n_size; j++)
     A[OBJrow-1][j-1]=-A[OBJrow-1][j-1];
   DualSimplexMaximize(m_size, n_size, A,BasisInverse, OBJrow, RHScol, UsePrevBasis, 
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter, err);
   *optvalue=-*optvalue;
   for (j=1; j<=n_size; j++){
     dsol[j-1]=-dsol[j-1];
     A[OBJrow-1][j-1]=-A[OBJrow-1][j-1];
   }
}

void DualSimplexMaximize(rowrange m_size, colrange n_size,
   Amatrix A, Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  int stop, chosen, phase1, found;
  long rank;
  long pivots_ds=0, pivots_p1=0, pivots_pc=0, maxpivots, maxpivfactor=70;
  rowrange i,r,entering,leaving;
  colrange j,s;
  static colset ColSelected;
  static rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0,nlast=0;
  double adet=1; /* abs value of the determinant of a basis */
  int localdebug=dp_FALSE;
  static rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
  unsigned int rseed=1;
  
  *err=dp_None;
  maxpivots=maxpivfactor*n_size;  /* maximum pivots to be performed before cc pivot is applied. */
  if (mlast!=m_size || nlast!=n_size){
     if (mlast>0) { /* called previously with different m_size */
       if (localdebug) printf("DualSimplex: deleting the old memory space with mlast = %ld\n",mlast);
       free(OrderVector);
       free(BasisFlag);
       set_free(ColSelected);
       set_free(RowSelected);
       set_free(Basis);
       set_free(Cobasis);
     }
     if (localdebug) printf("DualSimplex: allocating a new memory space with m_size = %ld\n",m_size);
     OrderVector=(long *)calloc(m_size+1, sizeof *OrderVector);
     BasisFlag=(long *) calloc(m_size+2, sizeof entering);  /* one more element for an auxiliary variable  */
     set_initialize(&Cobasis,m_size+1);
     set_initialize(&Basis,m_size+1);
     set_initialize(&RowSelected, m_size+1);
     set_initialize(&ColSelected, n_size);
     /* initialize only for the first time or when a larger space is needed */
     mlast=m_size;nlast=n_size;
  }
  /* Initializing control variables. */
  dp_ComputeRowOrderVector(m_size, n_size, A, OrderVector, dp_MinIndex, rseed);

  BasisFlag[OBJrow]= 0; /*  BasisFlag of the objective variable is 0, 
    different from other basic variables which have -1 */
  set_emptyset(RowSelected);
  set_emptyset(ColSelected);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  *re=0; *se=0; *iter=0;

  if (UsePrevBasis){
    if (s=BasisFlag[OBJrow]>0){ /* OBJrow is in cobasis, and we must pivot it in. */
      for (j=0;j<=m_size+1;j++){if (j!=s) set_addelem(ColSelected,j);}
      rank=n_size-2;
      if (localdebug) printf("UsePrevBasis but the current basis does not contain OBJrow.\n");
      set_copy(RowSelected,Cobasis);
    } else {
      if (localdebug){
        printf("UsePrevBasis. The current basis contains OBJrow and thus an LP basis.\n");
        printf(" Current nonbasis = ");   
        set_write(Cobasis); printf("\n");
      }
      goto _L88;
    }
  } else {
    if (localdebug) printf("Do not UsePrevBasis.\n");
    rank = 0;
    set_emptyset(Cobasis);
    set_emptyset(Basis);
    for (j=0; j<=n_size; j++) NBIndex[j]=0;
    for (i=1; i<=m_size+1; i++) {
      BasisFlag[i]=0;
      set_addelem(Basis,i);
      BasisFlag[i]=-1;    /* basic variable has index -1 */
    }
    dp_SetToIdentity(n_size, BasisInverse);
  }

  stop = dp_FALSE;
  do {   /* Find a LP basis */
    dp_SelectPivot2(m_size, n_size, A, BasisInverse, dp_MinIndex, OrderVector
      , m_size, RowSelected, ColSelected, &r, &s, &chosen);
    if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
    if (chosen) {
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
      set_addelem(Cobasis, r);
      set_delelem(Basis,r);
      if (localdebug) {
        printf("DSimplex: find initial basis: nonbasis = ");   
        set_write(Cobasis); printf("\n");
      }
      BasisFlag[r]=s;   /* the nonbasic variable r corresponds to column s */
      NBIndex[s]=r;     /* the nonbasic variable on s column is r */
      if (localdebug) printf("nonbasic variable %ld has index %ld\n", r, BasisFlag[r]);
      rank++;
      dp_GausianColumnPivot2(m_size, n_size, A, BasisInverse, r, s);
      if (localdebug) {
        dp_WriteBmatrix(stdout, n_size, BasisInverse);
        dp_WriteTableau(stdout,m_size, n_size, A, BasisInverse);
	    printf("%ld th row added to the initial set (%ld elem)\n", r,rank);
      }
    } else {
      stop=dp_TRUE;
    }
    if (rank==n_size-1) {
       stop = dp_TRUE;
    }
  } while (!stop);

_L88:
  stop=dp_FALSE;

  FindDualFeasibleBasis(m_size, n_size, A, BasisInverse, 
      Basis, Cobasis, OrderVector, NBIndex,BasisFlag, 
      OBJrow, RHScol, &s, &found, LPS, &pivots_p1);
  *iter+=pivots_p1;
  if (!found){
     *se=s;   
     /* No dual feasible basis is found, and thus DualInconsistent.  
     Output the evidence column. */
  }
  else {
    do {   /* Dual Simplex Method */
      chosen=dp_FALSE; *LPS=dp_LPSundecided; phase1=dp_FALSE;
      if (pivots_ds<maxpivots) {
        SelectDualSimplexPivot(m_size, n_size, 
          phase1, A, BasisInverse, OrderVector, NBIndex, BasisFlag,
          OBJrow, RHScol, &r, &s, &chosen, LPS);
      }
      if (chosen) pivots_ds=pivots_ds+1;
      if (!chosen && *LPS==dp_LPSundecided) {  
        /* In principle this should not be executed because we already have dual feasibility
           attained and dual simplex pivot should have been chosen.  This might occur
           under floating point computation, or the case of cycling.
        */
        dp_SelectCrissCrossPivot(m_size, n_size, A, BasisInverse, BasisFlag,
          OBJrow, RHScol, &r, &s, &chosen, LPS);
        if (chosen) pivots_pc=pivots_pc+1;
      }
      if (chosen) {
        entering=NBIndex[s];
        leaving=r;
        if (localdebug) {
          printf("Dual Phase I: pivot on (r,s) = (%ld, %ld)\n", 
            r,s);
          printf("Dual Phase I: (leaving, entering) = (%ld, %ld)\n", 
            leaving,entering);
        }
        set_addelem(Cobasis, leaving);
        set_delelem(Cobasis, entering);
        set_delelem(Basis,leaving);
        set_addelem(Basis,entering);
        BasisFlag[leaving]=s;
        BasisFlag[entering]=-1;
        NBIndex[s]=leaving;
        if (localdebug) {
          printf("nonbasis = "); set_write(Cobasis);
          printf("\nbasis = "); set_write(Basis);
          printf("\nnew nonbasic variable %ld  has index %ld\n",leaving, BasisFlag[leaving]);
        }
        dp_GausianColumnPivot2(m_size, n_size, A, BasisInverse, r, s);
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
  }
  if (localdebug){
     printf("LP solved with %ld pivots. (ds pivt#= %ld,  p1 piv#= %ld",*iter,pivots_ds, pivots_p1);
     if (pivots_pc > 0) printf(", cc piv#= %ld",pivots_pc);
     printf(")\n");
  }
  switch (*LPS){
  case dp_Optimal:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,j);
      *optvalue=dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,RHScol);
      if (localdebug) printf("dsol[%ld]= %lf\n",NBIndex[j],dsol[j-1]);
    }
    break;
  case dp_Inconsistent:
    if (localdebug) printf("DualSimplexSolve: LP is inconsistent.\n");
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,*re,j);
      if (localdebug)  printf("dsol[%ld]= %lf\n",NBIndex[j],dsol[j-1]);
    }
    break;
  case dp_DualInconsistent:
    for (j=1;j<=n_size; j++) {
      sol[j-1]=BasisInverse[j-1][*se-1];
      dsol[j-1]=-dp_TableauEntry(m_size, n_size, A, BasisInverse,OBJrow,j);
      if (localdebug)  printf("dsol[%ld]= %lf\n",NBIndex[j],dsol[j-1]);
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
    if (abs(temp1)>dp_zero) A[i-1][0]=temp2/temp1;  
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
   Amatrix A, Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *err)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  *err=dp_None;
  switch (lpconv) {
    case dp_LPmax:
      if (solver==dp_DualSimplex)
         DualSimplexMaximize(m_size, n_size, A, BasisInverse, OBJrow, RHScol, UsePrevBasis,
           LPS, optvalue, sol, dsol,NBIndex, re, se, iter, err);
      else
         CrissCrossMaximize(m_size, n_size, A, BasisInverse, OBJrow, RHScol, 
           LPS, optvalue, sol, dsol,NBIndex, re, se, iter, err);
      break;
      
    case dp_LPmin:
      if (solver==dp_DualSimplex)
         DualSimplexMinimize(m_size, n_size, A, BasisInverse, OBJrow, RHScol, UsePrevBasis,
           LPS, optvalue, sol, dsol,NBIndex, re, se, iter, err);
      else
         CrissCrossMinimize(m_size, n_size, A, BasisInverse, OBJrow, RHScol, 
           LPS, optvalue, sol, dsol,NBIndex, re, se, iter, err);
      break;
  }
}

void dp_WriteLPResult(FILE *f, dp_LPConversionType Conversion, dp_LPSolverType LPSolver,
    rowrange m_size, colrange n_size, 
    Amatrix A, rowrange objrow, colrange rhscol,
    dp_LPStatusType LPS, double optval,
    Arow sol, Arow dsol, colindex NBIndex, rowrange re, colrange se,
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
      fprintf(f,"  %3ld : ",NBIndex[j+1]);
      dp_WriteReal(f,dsol[j]);
      fprintf(f,"\n");
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
      fprintf(f,"  %3ld : ",NBIndex[j+1]);
      dp_WriteReal(f,dsol[j]);
      fprintf(f,"\n");
    }
    fprintf(f,"end\n");
    break;

  case dp_DualInconsistent:
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


