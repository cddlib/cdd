/* Header file for setoper.c  */

/* setoper.c: 
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * last modified on March 16, 1995
 */

#define SETBITS 32       /* Important Constant: Number of bits in a long integer */
#include <stdio.h>
#include <stdlib.h>

typedef long *set_type;  /* set type definition */

/* Definitions for optimized set_card function 
   by David Bremner bremner@cs.mcgill.ca  
 */

typedef unsigned char set_card_lut_t;

#define LUTBLOCKS(set) (((set[0]-1)/SETBITS+1)*(sizeof(long)/sizeof(set_card_lut_t)))

long set_blocks(long len);
void set_initialize(set_type *setp,long len);
void set_free(set_type set);
void set_emptyset(set_type set);
void set_copy(set_type setcopy,set_type set);
void set_addelem(set_type set, long elem);
void set_delelem(set_type set, long elem);
void set_int(set_type set,set_type set1,set_type set2);
void set_uni(set_type set,set_type set1,set_type set2);
void set_diff(set_type set,set_type set1,set_type set2);
void set_compl(set_type set,set_type set1);
int set_subset(set_type set1,set_type set2);
int set_member(long elem, set_type set);
long set_card(set_type set);
void set_write(set_type set);
void set_fwrite(FILE *f,set_type set);
void set_binwrite(set_type set);
void set_fbinwrite(FILE *f,set_type set);

/* End of File: setoper.h */
/* End of File: setoper.h */

