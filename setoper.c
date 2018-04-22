/* setoper.c:
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * modified on December 5, 1994 
   (set_card function replaced with a better code by David Bremner) 
 * last modified on March 16, 1995 
 */

#include "setoper.h"

static unsigned char set_card_lut[]={
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
/* End of Definitions for optimized set_card */

long set_blocks(long len)
{
	long blocks;
	
	blocks=(len-1)/SETBITS+2;
	return blocks;
}

void set_initialize(set_type *setp, long len)
/* Make a set with a given bit lengths  */
{
	long i,forlim1;
	
	forlim1=set_blocks(len);
	*setp=(long *) calloc(forlim1, sizeof i);
	(*setp)[0]=len;  /* size of the ground set */
	for (i=1; i<forlim1; i++)
		(*setp)[i]=0;
}

void set_free(set_type set)
/* Free the space created with the set pointer set*/
{
    free(set);
}

void set_emptyset(set_type set)
/* Set set to be the emptyset  */
{
	long i,forlim;
	
	forlim=set_blocks(set[0])-1;
	for (i=1; i<=forlim; i++)
		set[i]=0;
}


void set_copy(set_type setcopy,set_type set)
/* Copy the set set[] to setcopy[] with setcopy[] length */
{
	long i,forlim;

	forlim=set_blocks(setcopy[0])-1;
	for (i=1; i<=forlim; i++)
		setcopy[i]=set[i];
}

void set_addelem(set_type set, long elem)
/* add elem only if it is within the set[] range */
{
	long i,j;
	long change;
	long one=1;
	
	if (elem<=set[0])    
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change= one << j;  /* put 1 in jth position */
		set[i]=set[i] | change;
	}
}

void set_delelem(set_type set, long elem)
/* delete elem only if it is within the set[] range */
{
	long  i,j;
	long change;
	long one=1;	 
	
	if (elem<=set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change=one << j;  /* put 1 in jth position */
		set[i]=(set[i] | change) ^ change;
	}
}

void set_int(set_type set,set_type set1,set_type set2)
/* Set intersection, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;
	
	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]=(set1[i] & set2[i]);
}

void set_uni(set_type set,set_type set1,set_type set2)
/* Set union,assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] | set2[i];
}

void set_diff(set_type set,set_type set1,set_type set2)
/* Set difference se1/set2, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] & (~set2[i]);
}

void set_compl(set_type set,set_type set1)
/* set[] will be set to the complement of set1[] */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]= ~set1[i];
}

int set_subset(set_type set1,set_type set2)
/* Set containment check, set1 <= set2 */
{
	int  yes=1;
	long i,forlim;
	
	forlim=set_blocks(set2[0])-1;
	for (i=1;i<=forlim && yes;i++)
		if ((set1[i] | set2[i])!=set2[i])
			yes=0;
	return yes;
}

int set_member(long elem, set_type set)
/* Set membership check, elem in set */
{
	int  yes=0;
	long  i,j;
	long testset;
	long one=1;	 
	
	if (elem<=set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		testset=set[i] | (one<<j);   /* add elem to set[i] */
		if (testset==set[i])
			yes=1;
	}
	return yes;
}

long set_card(set_type set)
/* set cardinality, modified by David Bremner bremner@cs.mcgill.ca
   to optimize for speed. */
{
  long block,car=0;
  set_card_lut_t *p;
  
  p=(set_card_lut_t *)&set[1];
  for (block=0; block< LUTBLOCKS(set);block++) {
    car+=set_card_lut[p[block]];
  }
  return car;
}

/* old cardinality code 
long set_card(set_type set)
{
	long elem,car=0;
	
	for (elem=1; elem<=set[0]; elem++) {
		if (set_member(elem,set)) car++;
    }
	return car;
}
*/ 

void set_write(set_type set)
{
	long elem;
	
	for (elem=1;elem<=set[0];elem++)
	{
		if (set_member(elem,set))
			printf("%ld ",elem);
	}
	printf("\n");
}

void set_fwrite(FILE *f,set_type set)
{
	long elem;
	
	for (elem=1;elem<=set[0];elem++)
	{
		if (set_member(elem,set))
			fprintf(f,"%ld ",elem);
	}
	fprintf(f,"\n");
}

void set_binwrite(set_type set)
{
	int i,j;
	long forlim;
	long e1,e2;
	
	printf("max element = %ld,\n",set[0]);
	forlim=set_blocks(set[0])-1;
	for (i=forlim;i>=1;i--)
	{
		e1=e2=set[i];
		for (j=SETBITS-1;j>=0;j--)
		{
			e1=(e1>>j);
			printf("%1ld",e1);
			e1=e2-(e1<<j);
			e2=e1;
		}
		printf(" ");
	}
	printf("\n");
}


void set_fbinwrite(FILE *f,set_type set)
{
	int i,j;
	long forlim;
	long e1,e2;
	
	printf("max element = %ld,\n",set[0]);
	forlim=set_blocks(set[0])-1;
	for (i=forlim;i>=1;i--)
	{
		e1=e2=set[i];
		for (j=SETBITS-1;j>=0;j--)
		{
			e1=(e1>>j);
			fprintf(f,"%1ld",e1);
			e1=e2-(e1<<j);
			e2=e1;
		}
		printf(" ");
	}
	fprintf(f,"\n");
}

/* End of the library:  setoper.c  */
