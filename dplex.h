/* dplex.h: Header file for dplex.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.6, August 21, 1996
*/

/* dplex.c : C-Implementation of the dual simplex method for
   solving an LP: max/min  c^T x subject to  x in P, where
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/


#include <time.h>
#include "dplexdef.h"

typedef long rowrange;
typedef long colrange;
typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;
typedef long *rowindex;   
    /* rowindex should be intialized to be an array of [mm+1] components */
typedef long colindex[dp_NMAX+1];
typedef double *Amatrix[dp_MMAX];
typedef double Arow[dp_NMAX];
typedef double *Bmatrix[dp_NMAX];

typedef char dp_FilenameType[dp_filenamelen];

typedef enum {
  dp_DimensionTooLarge, dp_LowColumnRank, dp_ImproperInputFormat, 
  dp_FileNotFound, dp_None
} dp_ErrorType;

typedef enum {
  dp_Real, dp_Rational, dp_Integer, dp_Unknown
} dp_NumberType;

typedef enum {
  dp_LPmax, dp_LPmin
} dp_LPConversionType;

typedef enum {
  dp_CrissCross, dp_DualSimplex
} dp_LPSolverType;

typedef enum {
  dp_LPSundecided, dp_Optimal, dp_Inconsistent, dp_DualInconsistent, 
  dp_Unbounded, dp_DualUnbounded
} dp_LPStatusType;

void dp_LPInput(FILE **f, dp_FilenameType, rowrange *m, colrange *n, Amatrix A, 
    dp_LPConversionType *Conversion, rowrange *objrow, colrange *rhscol, 
    dp_ErrorType *err);

void dp_LPSolve(dp_LPConversionType, dp_LPSolverType, 
   rowrange, colrange, Amatrix, Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, int UsePrevBasis, dp_LPStatusType *LPS,
   double *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter, dp_ErrorType *);

void dp_WriteLPResult(FILE *, dp_LPConversionType, dp_LPSolverType,
  rowrange m_size, colrange n_size, Amatrix A, rowrange objrow, colrange rhscol,
  dp_LPStatusType, double, Arow, Arow, colindex, rowrange, colrange, long, dp_ErrorType);

void dp_WriteErrorMessages(FILE *, dp_ErrorType);

/* end of dplex.h */
