/* setoper.c:
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * modified on Dec. 20, 1993 
 */

#include "setoper.h"

long set_blocks(long len)
{
	long blocks;
	
	blocks=(len-1)/SETBITS+2;
	return blocks;
}

void set_initialize(long set[],long len)
/* Make a set with a given bit lengths  */
{
	long i,forlim;
	
	forlim=set_blocks(len)-1;
	set[0]=len;  /* size of the ground set */
	for (i=1; i<=forlim; i++)
		set[i]=0;
}

void set_copy(long setcopy[],long set[])
/* Copy the set set[] to setcopy[] with setcopy[] length */
{
	long i,forlim;

	forlim=set_blocks(setcopy[0])-1;
	for (i=1; i<=forlim; i++)
		setcopy[i]=set[i];
}

void set_addelem(long set[], long elem)
/* add elem only if it is within the set[] range */
{
	long i,j;
	long change;
	unsigned long one=1;
	
	if (elem<=set[0])    
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change= one << j;  /* put 1 in jth position */
		set[i]=set[i] | change;
	}
}

void set_delelem(long set[], long elem)
/* delete elem only if it is within the set[] range */
{
	long  i,j;
	long change;
	unsigned long one=1;	 
	
	if (elem<=set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change=one << j;  /* put 1 in jth position */
		set[i]=(set[i] | change) ^ change;
	}
}

void set_int(long set[],long set1[],long set2[])
/* Set intersection, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;
	
	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]=(set1[i] & set2[i]);
}

void set_uni(long set[],long set1[],long set2[])
/* Set union,assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] | set2[i];
}

void set_diff(long set[],long set1[],long set2[])
/* Set difference se1/set2, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] & (~set2[i]);
}

void set_compl(long set[],long set1[])
/* set[] will be set to the complement of set1[] */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]= ~set1[i];
}

int set_subset(long set1[],long set2[])
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

int set_member(long elem, long set[])
/* Set membership check, elem in set */
{
	int  yes=0;
	long  i,j;
	long testset;
	unsigned long one=1;	 
	
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

long set_card(long set[])
/* set cardinality  */
{
	long elem,car=0;
	
	for (elem=1; elem<=set[0]; elem++) {
		if (set_member(elem,set)) car++;
    }
	return car;
}

void set_write(long set[])
{
	long elem;
	
	for (elem=1;elem<=set[0];elem++)
	{
		if (set_member(elem,set))
			printf("%ld ",elem);
	}
	printf("\n");
}

void set_binwrite(long set[])
{
	int i,j;
	long forlim;
	unsigned long e1,e2;
	
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

/* End of the library:  setoper.c  */