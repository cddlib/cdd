/* dplexdef.h:  Definition file for dplex.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.60, August 21, 1996
*/

/* LP to be solved is of form
  maximize  c^T x  +   c0
  subj. to
            A   x  <=  b.
*/

#define dp_MMAX    5002 /* USER'S CHOICE: max row size of A plus two */
#define dp_NMAX    51   /* USER'S CHOICE: max column size of A plus one */
#define dp_zero      1.0e-5   /*real zero*/
#define dp_wordlenmax  32
#define dp_linelenmax  255
#define dp_filenamelen 32

/* end of dplexdef.h */

