 /* cdd.h: Header file for cdd.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.61b, November 29, 1997
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define COPYRIGHT   "Copyright (C) 1996, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DDVERSION   "Version 0.61b (November 29, 1997)"
#include <time.h>
#include "dplex.h"
#include "dplexdef.h"

typedef char boolean;

/* 
typedef long rowrange;
typedef long colrange;
typedef set_type rowset;
typedef set_type colset;
typedef long *rowindex;   
typedef long colindex[NMAX+1];
typedef double *Amatrix[MMAX];
typedef double Arow[NMAX];
typedef double *Bmatrix[NMAX];
*/


typedef char DataFileType[filenamelen];

/*
typedef char LineType[linelenmax];
typedef char WordType[wordlenmax];
*/

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
  CrissCross,DualSimplex,CombMaxImprove
} LPsolverType;

typedef enum {
  IncOff=0, IncCardinality, IncSet
} IncidenceOutputType;

typedef enum {
  AdjOff=0, AdjacencyList,  AdjacencyDegree
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

extern boolean OutputReordered;  /* if TRUE, output the reordered problem  */
extern long projdim;  /*dimension of orthogonal preprojection */
extern colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
extern rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1, CheckPoints;
extern rowindex EqualityIndex;  
extern rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
extern long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
 TotalRayCount, VertexCount,ZeroRayCount;
extern long EdgeCount,TotalEdgeCount;
extern long count_int,count_int_good,count_int_bad;
extern boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, debug;
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
extern CompStatusType CompStatus;     /* Computation Status */
extern ConversionType Conversion;
extern LPsolverType LPsolver;
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
void ProcessCommandLine(rowrange, colrange, char *);
void AmatrixInput(rowrange *, colrange *, rowrange *, colrange *, Amatrix, boolean *);
void SetInequalitySets(rowrange, rowindex);

void WriteReal(FILE *, double);
void WriteRayRecord(FILE *, colrange, RayRecord *);
void WriteRayRecord2(FILE *, colrange, RayRecord *);
double AValue(colrange, Amatrix, double *, rowrange );
void WriteIncidence(FILE *, rowrange, colrange, RayRecord *);
void StoreRay1(rowrange, colrange, Amatrix, double *, RayRecord *, rowindex, boolean *);
void StoreRay2(rowrange, colrange, Amatrix, double *, RayRecord *, rowindex, boolean *, boolean *);
void AddRay(rowrange, colrange, Amatrix, double *, rowindex);
void AddArtificialRay(rowrange, colrange, Amatrix, rowindex);
void ConditionalAddEdge(rowrange, colrange, RayRecord *Ray1, RayRecord *Ray2, 
    RayRecord *ValidFirstRay, rowrange, rowindex);
void CreateInitialEdges(rowrange, colrange, rowindex);
void UpdateEdges(rowrange, colrange, RayRecord *RRbegin, RayRecord *RRend, rowrange, rowindex);
void FreeDDMemory(rowindex);
void Normalize(colrange, double *);
void ZeroIndexSet(rowrange, colrange, Amatrix, double *, rowset);
void CopyBmatrix(colrange, Bmatrix T, Bmatrix TCOPY);
void SelectPivot1(rowrange, colrange, Amatrix, HyperplaneOrderType, rowindex,
    rowrange, rowset, colset, rowrange *, colrange *,boolean *);
double TableauEntry(rowrange, colrange, Amatrix, Bmatrix T, rowrange, colrange);
void WriteTableau(FILE *, rowrange, colrange, Amatrix, Bmatrix T,InequalityType);
void SelectPivot2(rowrange, colrange, Amatrix, Bmatrix T, HyperplaneOrderType, rowindex,
    rowrange, rowset, colset, rowrange *, colrange *,boolean *);
void GausianColumnPivot1(rowrange, colrange, Amatrix, rowrange, colrange,rowrange);
void GausianColumnPivot2(rowrange, colrange, Amatrix, Bmatrix,  rowrange, colrange);
void InitializeBmatrix(colrange, Bmatrix T);
void SetToIdentity(colrange, Bmatrix T);
void free_Bmatrix(colrange, Bmatrix T);
void WriteBmatrix(FILE *, colrange, Bmatrix T);
void ReduceA(rowrange *, colrange *, Amatrix, rowset, colset);
void DualizeA(rowrange *, colrange *, Amatrix, Bmatrix);
void EnlargeAforInteriorFinding(rowrange *, colrange *, Amatrix);
void WriteSubMatrixOfA(FILE *, rowrange, colrange, Amatrix, rowset, colset, InequalityType);
void WriteAmatrix(FILE *, Amatrix, long, long, InequalityType);
void WriteAmatrix2(FILE *f, Amatrix A, long rowmax, long colmax,
      InequalityType ineq, rowindex OV, rowrange iteration);
void ComputeRank(rowrange, colrange, Amatrix, rowset, rowindex, long *);
void ComputeBInverse(rowrange, colrange, Amatrix, long, rowindex, Bmatrix InvA1, long *);
void FindBasis(rowrange, colrange, Amatrix, HyperplaneOrderType, rowindex,rowset,colindex,
    Bmatrix, long *rank);
void SelectCrissCrossPivot(rowrange, colrange, Amatrix, Bmatrix T,
  long bflag[], rowrange,colrange,rowrange *,colrange *,
  boolean *, dp_LPStatusType *);
void SelectDualSimplexPivot(rowrange, colrange, boolean, Amatrix, Bmatrix, rowindex,
    colindex, long *, rowrange, colrange,
    rowrange *, colrange *, boolean *, dp_LPStatusType *);
void FindInitialRays(rowrange, colrange, Amatrix, rowindex, rowset, Bmatrix, 
    colindex, boolean *);
void CheckAdjacency1(rowrange, colrange, Amatrix, rowindex, RayRecord **, RayRecord **, boolean *);
void CheckAdjacency2(rowrange, colrange, Amatrix, RayRecord **, RayRecord **, boolean *);
void CheckEquality(colrange, RayRecord **, RayRecord **, boolean *);
void Eliminate(colrange, RayRecord **);
void CreateNewRay(rowrange, colrange, Amatrix, RayRecord *, RayRecord *, rowrange, rowindex);
void EvaluateARay1(rowrange, colrange n_size, Amatrix);
void EvaluateARay2(rowrange, colrange n_size, Amatrix);
void FeasibilityIndices(long *, long *, rowrange, colrange n_size, Amatrix);
boolean LexSmaller(double *, double *, long);
boolean LexLarger(double *, double *, long);
void CopyArow(double *, double *, long);
void ComputeRowOrderVector(rowrange, colrange, Amatrix, rowindex OV, HyperplaneOrderType ho);
void UpdateRowOrderVector(rowrange, colrange, rowset PriorityRows, rowindex);
void SelectNextHyperplane(rowrange, colrange, Amatrix, HyperplaneOrderType,
    rowset, rowrange *, boolean *, rowindex);
void SelectPreorderedNext(rowrange, colrange, rowset, rowindex, rowrange *hnext);
void AddNewHyperplane1(rowrange, colrange, Amatrix, rowrange, rowrange iter, rowindex);
void AddNewHyperplane2(rowrange, colrange, Amatrix, rowrange, rowrange iter, rowindex);
void WriteAdjacency(FILE *, rowrange, colrange, rowrange, colrange, Amatrix);
void WriteAdjacencyDegree(FILE *, rowrange, colrange, rowrange, colrange, Amatrix, rowrange);
void WriteRunningMode0(FILE *);
void WriteRunningMode(FILE *);
void WriteRunningMode2(FILE *, rowrange, colrange);
void WriteCompletionStatus(FILE *,  rowrange, colrange, rowrange);
void WriteTimes(FILE *);
void WriteProgramDescription(FILE *);
void WriteSolvedProblem(FILE *, rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix);
void WriteSolvedSubProblem(FILE *, rowrange m_input, colrange n_input, 
    rowrange m_size, colrange n_size, Amatrix, rowindex OV, 
    rowrange iteration);
void InitialWriting(rowrange, colrange, rowrange, colrange);
void WriteDDResult(rowrange, colrange, rowrange, colrange, Amatrix, rowrange);
void WriteErrorMessages(FILE *);
void WriteDecompResult(rowrange, colrange, rowrange, colrange, rowrange iter);
void WriteProjResult(rowrange, colrange, rowrange, colrange, Amatrix, long *, rowrange iter);
void WriteHeading(void);

/* end of cdd.h */
