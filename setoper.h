/* Header file for setoper.c  */

/* setoper.c: 
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * modified on Dec. 8, 1993
 */

#define SETBITS 32   /* Important Constant: Number of bits in a long integer */
#include <stdio.h>

long set_blocks(long len);
void set_initialize(long set[],long len);
void set_copy(long setcopy[],long set[]);
void set_addelem(long set[], long elem);
void set_delelem(long set[], long elem);
void set_int(long set[],long set1[],long set2[]);
void set_uni(long set[],long set1[],long set2[]);
void set_diff(long set[],long set1[],long set2[]);
void set_compl(long set[],long set1[]);
int set_subset(long set1[],long set2[]);
int set_member(long elem, long set[]);
void set_write(long set[]);
void set_binwrite(long set[]);

/* End of File: setoper.h */