/* cdd.h: Header file for cdd.c 
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.55a, December 18, 1994 
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define COPYRIGHT   "Copyright (C) 1994, Komei Fukuda, fukuda@dma.epfl.ch"
#define DDVERSION   "Version C0.55a (December 18, 1994)"
#include <time.h>

typedef char boolean;
typedef long rowrange;
typedef long colrange;
typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;
typedef long *rowindex;   
    /* rowindex should be intialized to be an array of [mm+1] components */
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
  rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  boolean feasible;  /* flag to store the feasibility */
  double ARay;   /* temporary area to store some row of A*Ray */
  struct RayRecord *Next;
} RayRecord;

typedef struct AdjacencyRecord {
  RayRecord *Ray1, *Ray2;
  struct AdjacencyRecord *Next;
} AdjacencyRecord;

typedef struct node {long key; struct node *next;} node;

typedef enum {
  Combinatorial, Algebraic
} AdjacencyTestType;

typedef enum {
  MaxIndex, MinIndex, MinCutoff, MaxCutoff, MixCutoff,
   LexMin, LexMax, RandomRow, LineShelling
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
extern rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1;
extern rowrange Iteration, hh;
extern rowindex OrderVector;
extern rowindex EqualityIndex;  
extern rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
extern long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
 TotalRayCount, VertexCount,ZeroRayCount;
extern long EdgeCount,TotalEdgeCount;
extern long count_int,count_int_good,count_int_bad;
extern boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
extern Amatrix AA;
extern Bmatrix InitialRays;
extern colindex InitialRayIndex;
extern LPStatusType LPStatus;
extern colrange RHScol;   /* LP RHS column */
extern rowrange OBJrow;   /* LP OBJ row */
extern Arow LPcost;
extern RayRecord *ArtificialRay, *FirstRay, *LastRay;
extern RayRecord *PosHead, *ZeroHead, *NegHead, *PosLast, *ZeroLast, *NegLast;
extern AdjacencyRecord *Edges[MMAX];  /* adjacency relation storage for iteration k */
extern boolean RecomputeRowOrder, inputsuccessful;
extern HyperplaneOrderType HyperplaneOrder;
extern AdjacencyTestType AdjacencyTest;
extern NumberType Number;
extern InequalityType Inequality;
extern boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
extern boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
extern boolean RestrictedEnumeration; /* Restricted enumeration Switch (TRUE if it is restricted on the intersection of EqualitySet hyperplanes) */
extern boolean RelaxedEnumeration; /* Relaxed enumeration Switch (TRUE if NonequalitySet inequalities must be satisfied with strict inequality) */
extern boolean RowDecomposition; /* Row decomposition enumeration Switch */
extern boolean VerifyInput; /* Verification switch for the input data */
extern boolean PreOrderedRun; 
extern boolean QPivotOn; /* QPivot Switch (TRUE if Q-Pivot scheme of Jack Edmonds is chosen) */
extern CompStatusType CompStatus;     /* Computation Status */
extern ConversionType Conversion;
extern IncidenceOutputType IncidenceOutput;
extern AdjacencyOutputType AdjacencyOutput;
extern ErrorType Error;
extern FileInputModeType FileInputMode;
extern DataFileType inputfile,ifilehead,ifiletail,
     outputfile,projfile, icdfile,adjfile,logfile,dexfile,verfile;
extern FILE *reading, *writing, *writing_icd,
     *writing_adj, *writing_log, *writing_dex, *writing_ver, *reading_dex;
extern time_t starttime, endtime;
extern unsigned int rseed;


void SetInputFile(FILE **, boolean *);
void SetWriteFile(FILE **, DataFileType, char, char *);


void SetNumberType(char *);
void ProcessCommandLine(char *);
void AmatrixInput(boolean *);
void SetInequalitySets(rowindex);

void WriteReal(FILE *, double);
void WriteSetElements(FILE *, long *);
void WriteRayRecord(FILE *, RayRecord *);
void WriteRayRecord2(FILE *, RayRecord *);
double AValue(double *, rowrange );
void WriteIncidence(FILE *, RayRecord *);
void StoreRay1(double *, RayRecord *, boolean *);
void StoreRay2(double *, RayRecord *, boolean *, boolean *);
void AddRay(double *);
void AddArtificialRay(void);
void ConditionalAddEdge(RayRecord *Ray1, RayRecord *Ray2, RayRecord *ValidFirstRay);
void CreateInitialEdges(void);
void UpdateEdges(RayRecord *RRbegin, RayRecord *RRend);
void FreeDDMemory(void);
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
void GausianColumnPivot2(Amatrix, Bmatrix,  rowrange, colrange);
void GausianColumnQPivot2(Amatrix, Bmatrix, double *,rowrange, colrange);
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
void CrissCrossMinimize(Amatrix,Bmatrix BasisInverse,
  rowrange, colrange, 
  LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void CrissCrossMaximize(Amatrix,Bmatrix BasisInverse,
  rowrange, colrange, 
  LPStatusType *, double *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void WriteLPResult(FILE *, LPStatusType, double, 
  Arow, Arow, colindex, rowrange, colrange, long);
void FindInitialRays(rowset, Bmatrix, colindex, boolean *);
void CheckAdjacency1(RayRecord **, RayRecord **,boolean *);
void CheckAdjacency2(RayRecord **, RayRecord **,boolean *);
void CheckEquality(RayRecord **, RayRecord **, boolean *);
void Eliminate(RayRecord **);
void CreateNewRay(RayRecord *, RayRecord *, rowrange);
void EvaluateARay1(rowrange);
void EvaluateARay2(rowrange);
void FeasibilityIndices(long *, long *, rowrange);
boolean LexSmaller(double *, double *, long);
boolean LexLarger(double *, double *, long);
void CopyArow(double *, double *, long);
void ComputeRowOrderVector(rowindex OV, HyperplaneOrderType ho);
void UpdateRowOrderVector(rowset PriorityRows);
void SelectNextHyperplane(HyperplaneOrderType,long *, rowrange *, boolean *);
void SelectPreorderedNext(long *excluded, rowindex, rowrange *hnext);
void AddNewHyperplane1(rowrange);
void AddNewHyperplane2(rowrange);
void WriteAdjacency(FILE *);
void WriteRunningMode(FILE *);
void WriteRunningMode2(FILE *);
void WriteCompletionStatus(FILE *);
void WriteTimes(FILE *);
void WriteProgramDescription(FILE *f);
void WriteSolvedProblem(FILE *f);


/* end of cdd.h */
