/* cdddef.h:  Definition file for cdd.c 
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.33,  Jan. 16, 1994 
*/

/* cdd.c : C-Implementation of the double description method for
   computing all vertices and extremal rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define MMAX           101  /* USER'S CHOICE: max row size of A plus one */
#define NMAX            33   /* USER'S CHOICE: max column size of A plus one */

#define SETBITS 32        /* Important Constant: Number of bits in a long integer    */
#define rowsetsize MMAX   /* The size of the column index set */
#define colsetsize NMAX   /* The size of the row index set */
#define rowsetblocks (rowsetsize-1)/SETBITS+2   /* Number of (long) blocks for row set */
#define colsetblocks (colsetsize-1)/SETBITS+2   /* Number of (long) blocks for column set */

#define datawidth       10
#define filenamelen     255 
#define wordlenmax      128 
#define linelenmax      255

#define FALSE 0
#define TRUE 1

#define zero            1.0e-5   /*real zero*/

/* End of cdddef.h */
