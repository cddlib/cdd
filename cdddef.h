/* cdddef.h:  Definition file for cdd.c 
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.38,  Jan. 31, 1994
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define MMAX      4001  /* USER'S CHOICE: max row size of A plus one */
#define NMAX      101   /* USER'S CHOICE: max column size of A plus one */

#define rowsetsize MMAX   /* The size of the column index set */
#define colsetsize NMAX   /* The size of the row index set */

#define datawidth       10
#define filenamelen     255 
#define wordlenmax      128 
#define linelenmax      255

#define FALSE 0
#define TRUE 1

#define zero            1.0e-5   /*real zero*/

/* End of cdddef.h */
