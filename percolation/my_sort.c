#include <stdio.h>
#include <math.h>
#include "my_sort.h"
/* *************************************************************** */
/*       Qsort is based on the book written by Kruse et al         */
/* *************************************************************** */

void q_sort_double_with_indx(double *ary,int *indx_lst,int low,int high)
{
	int pivotloc;

	if(low<high){
		pivotloc=q_partition_double_with_indx(ary,indx_lst,low,high);
		q_sort_double_with_indx(ary,indx_lst,low,pivotloc-1);
		q_sort_double_with_indx(ary,indx_lst,pivotloc+1,high);
	}
}

int q_partition_double_with_indx(double *ary,int *indx_lst,int low,int high)
{
	int i,pivotloc;
	double pivotkey;

	swap_double_with_indx(ary,indx_lst,low,(low+high)/2);
	pivotkey=ary[low];
	pivotloc=low;

	for(i=low+1;i<=high;i++){
		if(ary[i]<pivotkey){
			swap_double_with_indx(ary,indx_lst,++pivotloc,i);
		}
	}
	swap_double_with_indx(ary,indx_lst,low,pivotloc);
	return pivotloc;
}

void swap_double_with_indx(double *ary,int *indx_lst,int i,int j)
{
	double tmp_ary;
	int tmp_lst;


	tmp_ary=ary[i]; tmp_lst=indx_lst[i];
	ary[i]=ary[j]; indx_lst[i]=indx_lst[j];
	ary[j]=tmp_ary; indx_lst[j]=tmp_lst;
}
