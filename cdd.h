/* cdd.h: Header file for cdd.c 
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.33, Jan. 16, 1994 
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <time.h>

typedef char boolean;
typedef long rowrange;
typedef long colrange;
typedef long rowset[rowsetblocks];
typedef long colset[colsetblocks];
typedef long rowindex[MMAX+1];
typedef long colindex[NMAX+1];
typedef double Amatrix[MMAX][NMAX];
typedef double Arow[NMAX];
typedef double Bmatrix[NMAX][NMAX];
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
  IneToExt, ExtToIne, Projection, LPmax, LPmin
} ConversionType;
typedef enum {
  IncOff, IncCardinality, IncSet
} IncidenceOutputType;
typedef enum {
  DimensionTooLarge, LowColumnRank, ImproperInputFormat, DependentMarkedSet, None
} ErrorType;
typedef enum {
  InProgress, AllFound, RegionEmpty
} CompStatusType;
typedef enum {
  LPSundecided, Optimal, Inconsistent, DualInconsistent, Unbounded, DualUnbounded
} LPStatusType;

typedef char LineType[linelenmax];
typedef char WordType[wordlenmax];

extern long minput, ninput;   /*size of input data [b -A] */
extern long mm, nn;   /*size of the homogenous system to be solved by dd*/
extern rowrange OBJrow;
extern colrange RHScol;
extern long projdim;  /*dimension of orthogonal preprojection */
extern colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
extern rowset MarkedSet, GroundSet;
extern rowrange Iteration, hh;
extern rowset AddedHyperplanes, InitialHyperplanes;
extern long RayCount, FeasibleRayCount, TotalRayCount, VertexCount;
extern boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
extern Amatrix AA;
extern Bmatrix InitialRays;
extern Arow LPcost;  /* LP cost vector to be maximized  */
extern RayRecord *ArtificialRay, *FirstRay, *LastRay;
extern boolean inputsuccessful;
extern HyperplaneOrderType HyperplaneOrder;
extern AdjacencyTestType AdjacencyTest;
extern NumberType Number;
extern InequalityType Inequality;
extern boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
extern boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
extern boolean PartialEnumeration; /* Partial enumeration Switch (TRUE if it is restricted on the intersection of MarkedSet hyperplanes) */
extern CompStatusType CompStatus;     /* Computation Status */
extern ConversionType Conversion;
extern IncidenceOutputType IncidenceOutput;
extern ErrorType Error;
extern DataFileType inputfile,outputfile,icdfile,logfile;
extern FILE *reading, *writing, *writing_icd, *writing_log;
extern time_t starttime, endtime;

void SetInputFile(FILE **);
void SetWriteFile(FILE **);
void SetIncidenceFile(FILE **);
void SetLogFile(FILE **);
void SetNumberType(char *);
void ProcessCommandLine(char *);
void AmatrixInput(boolean *);

long Cardinality(long *);
void WriteReal(FILE *, double);
void WriteSetElements(FILE *, long *);
void WriteRayRecord(FILE *, RayRecord *);
void WriteRayRecord2(FILE *, RayRecord *);
double AValue(double *, rowrange );
void WriteIncidence(FILE *, RayRecord *);
void StoreRay(double *, RayRecord *, boolean *);
void AddRay(double *);
void AddArtificialRay(void);
void Normalize(double *);
void ZeroIndexSet(double *, long *);
void CopyBmatrix(double (*T)[NMAX], double (*TCOPY)[NMAX]);
void SelectPivot1(double (*X)[NMAX], HyperplaneOrderType,
   rowrange, long *,long *, rowrange *, colrange *,boolean *);
double TableauEntry(double (*X)[NMAX], double (*T)[NMAX], rowrange, colrange);
void WriteTableau(FILE *,double (*X)[NMAX], double (*T)[NMAX],InequalityType);
void SelectPivot2(double (*X)[NMAX], double (*T)[NMAX],
   HyperplaneOrderType,rowrange, long *,long *, rowrange *, colrange *,boolean *);
void GausianColumnPivot1(double (*X)[NMAX], rowrange, colrange,rowrange);
void GausianColumnPivot2(double (*X)[NMAX], double (*T)[NMAX],rowrange, colrange);
void InitializeBmatrix(double (*T)[NMAX]);
void WriteBmatrix(FILE *, double (*T)[NMAX]);
void ReduceAA(long *, long *);
void DualizeAA(double (*T)[NMAX]);
void WriteSubMatrixOfAA(FILE *,long *, long *, InequalityType);
void WriteAmatrix(FILE *, double (*A)[NMAX], long, long, InequalityType);
void ComputeRank(double (*A1)[NMAX], long *, long *);
void ComputeBInverse(double (*A1)[NMAX], long, double (*InvA1)[NMAX], long *);
void FindBasis(double (*A1)[NMAX], HyperplaneOrderType, long *, long *,
   double (*BasisInverse)[NMAX], long *);
void SelectCrissCrossPivot(double (*X)[NMAX], double (*T)[NMAX],
  long bflag[], rowrange,colrange,rowrange *,colrange *,
  boolean *, LPStatusType *);
void CrissCrossSolve(double (*A1)[NMAX],double (*BasisInverse)[NMAX],
   rowrange, colrange, double *);
void CheckAdjacency1(RayRecord **, RayRecord **,boolean *);
void CheckAdjacency2(RayRecord **, RayRecord **,boolean *);
void CheckEquality(RayRecord **, RayRecord **, boolean *);
void Eliminate(RayRecord **);
void CreateNewRay(RayRecord *, RayRecord *, rowrange);
void EvaluateARay(rowrange);
void FeasibilityIndices(long *, long *, rowrange);
boolean LexSmaller(double *, double *);
boolean LexLarger(double *, double *);
void CopyArow(double *, double *);
void SelectNextHyperplane(HyperplaneOrderType,long *, rowrange *);
void AddNewHyperplane(rowrange);
void WriteRunningMode(FILE *);
void WriteCompletionStatus(FILE *);
void WriteTimes(FILE *);

/* End. */
