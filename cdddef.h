/* cdddef.h:  Definition file for cdd.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.61b, November 29, 1997
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include "dplexdef.h"

#define MMAX      dp_MMAX  /* USER'S CHOICE: max row size of A plus one */
#define NMAX      dp_NMAX   /* USER'S CHOICE: max column size of A plus one */

#define rowsetsize MMAX   /* The size of the column index set */
#define colsetsize NMAX   /* The size of the row index set */

#define datawidth       10
#define filenamelen     256 
#define wordlenmax      128 
#define linelenmax      256

#define FALSE 0
#define TRUE 1

#define zero            dp_zero   /*real zero*/

/* end of cdddef.h */

