/* dplex_test.c: Main test program to call the dplex library
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.61, December 1, 1997
   Standard ftp site: ftp.ifor.math.ethz.ch, Directory: pub/fukuda/cdd
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

#include "setoper.h"
#include "dplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

FILE *reading, *writing;


void SetWriteFile(FILE **f)
{
  char *fname;
  fname="dplex_test.out";
  *f = fopen(fname, "w");
  printf("file %s is open\n",fname);
}


void main(int argc, char *argv[])
{
 /* Variables one must declare before calling dp_LPSolve  */
  rowrange mm;   /* row size of LP to be solved by cdd = #constraints + 1 */
  colrange nn;   /* column size of LP to be solved by cdd = #variables + 1 */
  Amatrix AA;
  /* The original LP data  mmxnn matrix 
     = | b   -A  |
       | c0  c^T |,
   
  where the LP to be solved is to
  maximize  c^T x  +   c0
  subj. to
            A   x  <=  b.
  */
        
  colindex NBIndex;    /*  NBIndex[j] is to store the jth nonbasic variable, j=0,1,...,nn-1*/ 
  Arow LPsol, LPdsol;  /*  LP solution x* and the dual solution y* (basic vars only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se,j;  /* evidence col when LP is dual-inconsistent */
  double ov;    /* LP optimum value */
  long LPiter;  /* iteration (=number of pivots) number */
  Bmatrix BasisInverse;    /* dual basis inverse matrix. */ 
  int UsePrevBasis=0;      /* set this variable to be 0  */
  rowrange OBJrow;         /* objective row = the last row of A (= mm-1 st row in C) */
  colrange RHScol;         /* rhs column = the first column of A (= 0th column in C) */
  dp_LPSolverType solver;      /* to be used to specify algorithm */
  dp_LPConversionType lpconv;  /* to be used to specify maximization or minimization */
  dp_ErrorType error=dp_None;  /* In case of errors, this will return the clue */
  dp_LPStatusType LPStatus;    /* to be used to notify the status of current solution */
  dp_FilenameType inputfile;   
    /* This filename (array of char's) is necessary if one uses dp_LPInput to input an LP */
   
  while (error==dp_None) {
    /* Input an LP using the dplex library  */
    dp_LPInput(&reading, inputfile, &mm, &nn, AA, &lpconv, &OBJrow, &RHScol, &error);

    if (error!=dp_None) {
      dp_WriteErrorMessages(stdout, error);
      goto _L99;
    }

    SetWriteFile(&writing);

    dp_InitializeBmatrix(nn, BasisInverse);  
      /* One must initialize BasisInverse before calling dp_LPSolve */

    solver=dp_DualSimplex;   /* either dp_DualSimplex or dp_CrissCross  */
  
    /* Find an interior point with dplex */
    printf("\n--- Running dp_FindInteriorPoint ---\n");
    dp_FindInteriorPoint(solver, mm, nn, AA, OBJrow, RHScol, 
      &LPStatus, &ov, LPsol, &LPiter, &error);
    if (dp_Positive(ov)){
      printf("An interior point found: (");
      for (j=1; j<nn; j++) dp_WriteReal(stdout,LPsol[j]);
      printf(")\n");
    }
    if (dp_Negative(ov)) printf("The feasible region is empty.\n");
    if (dp_Zero(ov)) printf("The feasible region is nonempty but has no interior point.\n");

    /* Solve the LP by dplex LP solver  */
    printf("\n--- Running dp_LPSolve ---\n");
    dp_LPSolve(lpconv, solver, mm, nn, AA, BasisInverse, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol, NBIndex, &re, &se, &LPiter, &error);

    /* Write the LP solutions by dplex LP reporter  */
    dp_WriteLPResult(stdout, lpconv, solver, mm, nn, AA, OBJrow, RHScol,
      LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
    dp_WriteLPResult(writing, lpconv, solver, mm, nn, AA, OBJrow, RHScol,
      LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter, error);
  }
_L99:;
}

/* end of dplex_test.c */
