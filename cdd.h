/* cdd.h: Header file for cdd.c 
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.38, Jan. 31, 1994 
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
typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;
typedef long rowindex[MMAX+1];
typedef long colindex[NMAX+1];
typedef double *Amatrix[MMAX];
typedef double Arow[NMAX];
typedef double *Bmatrix[NMAX];
typedef char DataFileType[filenamelen];
typedef char LineType[linelenmax];
typedef char WordType[wordlenmax];

typedef struct RayRecord {
  double *Ray;
  rowset ZeroSet;
  double ARay;   /* temporary area to store some row of A*Ray */
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
  IneToExt, ExtToIne, Projection, LPmax, LPmin, InteriorFind
} ConversionType;
typedef enum {
  IncOff=0, IncCardinality, IncSet
} IncidenceOutputType;
typedef enum {
  AdjOff=0, OutputAdjacency, InputAdjacency, IOAdjacency
} AdjacencyOutputType;
typedef enum {
  Auto, SemiAuto, Manual
} FileInputModeType;   
   /* Auto if a input filename is specified by command arguments */
typedef enum {
  DimensionTooLarge, LowColumnRank, ImproperInputFormat, DependentMarkedSet, 
  FileNotFound, None
} ErrorType;
typedef enum {
  InProgress, AllFound, RegionEmpty
} CompStatusType;
typedef enum {
  LPSundecided, Optimal, Inconsistent, DualInconsistent, Unbounded, DualUnbounded
} LPStatusType;

extern long minput, ninput;   /*size of input data [b -A] */
extern long mm, nn;   /*size of the homogenous system to be solved by dd*/
extern long projdim;  /*dimension of orthogonal preprojection */
extern colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
extern rowset MarkedSet, GroundSet, Face, Face1;
extern rowrange Iteration, hh;
extern rowset AddedHyperplanes, InitialHyperplanes;
extern long RayCount, FeasibleRayCount, TotalRayCount, VertexCount;
extern boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
extern Amatrix AA;
extern Bmatrix InitialRays;
extern LPStatusType LPStatus;
extern colrange RHScol;   /* LP RHS column */
extern rowrange OBJrow;   /* LP OBJ row */
extern Arow LPcost;
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
extern AdjacencyOutputType AdjacencyOutput;
extern ErrorType Error;
extern FileInputModeType FileInputMode;
extern DataFileType inputfile,ifilehead,ifiletail,
     outputfile,projfile, icdfile,adjfile,logfile;

extern FILE *reading, *writing, *writing_icd,*writing_adj, *writing_log;
extern time_t starttime, endtime;

void SetInputFile(FILE **, boolean *);
void SetWriteFile(FILE **, DataFileType, char, char *);


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
void CopyBmatrix(Bmatrix T, Bmatrix TCOPY);
void SelectPivot1(Amatrix, HyperplaneOrderType,
   rowrange, long *,long *, rowrange *, colrange *,boolean *);
double TableauEntry(Amatrix, Bmatrix T, rowrange, colrange);
void WriteTableau(FILE *,Amatrix, Bmatrix T,InequalityType);
void SelectPivot2(Amatrix, Bmatrix T,
   HyperplaneOrderType,rowrange, long *,long *, rowrange *, colrange *,boolean *);
void GausianColumnPivot1(Amatrix, rowrange, colrange,rowrange);
void GausianColumnPivot2(Amatrix, Bmatrix T,rowrange, colrange);
void InitializeBmatrix(Bmatrix T);
void SetToIdentity(Bmatrix T);
void free_Bmatrix(Bmatrix T);
void WriteBmatrix(FILE *, Bmatrix T);
void ReduceAA(long *, long *);
void DualizeAA(Bmatrix T);
void EnlargeAAforInteriorFinding(void);
void WriteSubMatrixOfAA(FILE *,long *, long *, InequalityType);
void WriteAmatrix(FILE *, Amatrix, long, long, InequalityType);
void ComputeRank(Amatrix, long *, long *);
void ComputeBInverse(Amatrix, long, Bmatrix InvA1, long *);
void FindBasis(Amatrix, HyperplaneOrderType, long *, long *,
   Bmatrix BasisInverse, long *);
void SelectCrissCrossPivot(Amatrix, Bmatrix T,
  long bflag[], rowrange,colrange,rowrange *,colrange *,
  boolean *, LPStatusType *);
void CrissCrossSolve(Amatrix,Bmatrix BasisInverse,
  rowrange, colrange, 
  LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void WriteLPResult(FILE *, LPStatusType, double, 
  Arow, Arow, colindex, rowrange, colrange, long);
void FindInitialRays(rowset InitHyperplanes,
			    Bmatrix InitRays, boolean *found);
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
void WriteAdjacency(FILE *);
void WriteRunningMode(FILE *);
void WriteCompletionStatus(FILE *);
void WriteTimes(FILE *);

/* End. */
