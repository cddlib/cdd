#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randlp.h"

void RandomRowVector(Arow row, long nn, long minvalue, long maxvalue)
{
  long j;

  for (j=1; j<=nn ;j++) {
    row[j-1] =  -( rand()%(maxvalue-minvalue+1) + minvalue);
  }
}

void GenerateKuhnQuandt(Amatrix A, long m, long n)
{  /* generate a dual Kuhn-Quandt instance */
  long i, j;
  long maxv=0,minv=-1000;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-10000;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=-10000;
  }
}

void GenerateSignedKuhnQuandt(Amatrix A, long m, long n)
{
  long i, j;
  long maxv=1000,minv=-1000;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-10000;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=10000;
  }
}

void GenerateZeroOne(Amatrix A, long m, long n)
{
  long i, j;
  long maxv=1, minv=-1;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-1;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=1;
  }
}

void GenerateDegenerate(Amatrix A, long m, long n)
{
  long i, j;
  long maxv=1000, minv=-1000;

  for (i=1; i<=m; i++){
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=0;
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=10000;
  }
}


void WriteAmatrix(FILE *f, Amatrix A, long m, long n, unsigned seed, int lptype)
{
  long i, j;

  fprintf(f, "LP type =  %d   Seed = %u\n", lptype, seed);
  fprintf(f, "begin\n  %ld  %ld  integer\n", m, n);
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      fprintf(f, " %ld", A[i-1][j-1]);
    }
    fprintf(f, "\n");
  }
  fprintf(f,"end\nmaximize ");
  for (j=1; j<=n; j++){
    fprintf(f," %ld", A[m][j-1]);
  }
  fprintf(f, "\n");
  fprintf(f, "!criss-cross\n");
}


void main()
{
  long mm,nn,i,j;
  unsigned sd=1;
  int lptype=1;
  Amatrix AA;
  
  FILE *writing;
  char outputfile[20];
  char *s, *numbtype;
  long value;


  printf("\nINPUT m n seed LPtype\n");
  printf("  LPtype 1 (dual Kuhn-Quandt), 2 (Signed KQ), 3 (-1,0,1) or 0 (degenerate).\n");
  scanf("%ld %ld %d %d", &mm, &nn, &sd, &lptype);
  if (nn > NMAX-1 || mm > MMAX-1) {
    printf("Size too large. Recompile with larger MMAX or NMAX\n");
    goto _L99;
  } else {
    srand(sd);
    for (i=0; i<=mm; i++){
      AA[i]=(long *) calloc(nn, sizeof value);
    }
    switch (lptype) {
        case 1:
          GenerateKuhnQuandt(AA, mm, nn);
          break;

        case 2:
          GenerateSignedKuhnQuandt(AA, mm, nn);
          break;

        case 3:
          GenerateZeroOne(AA, mm, nn);
          break;

        case 0:
          GenerateDegenerate(AA, mm, nn);
          break;
    }
    
    printf("Output (*.ine) file for a random lp problem: ");
    scanf("%s",outputfile);

    writing=fopen(outputfile,"w+");
    WriteAmatrix(writing, AA,mm,nn,sd,lptype);
  }

_L99:;
}

