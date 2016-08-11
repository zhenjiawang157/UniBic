
#include "struct.h" 

/****************************************************/
void GetLCSLen(int *s1,int *s2,Matrix pc,Matrix pb,int nrow,int ncolumn)
{
	int i,j;
	/************initial the edge***************/
	for(i=0; i<nrow; i++)
	{
		pc[i][0] = 0;
		pb[i][0] = 0;
	}
	for(j=0; j<ncolumn; j++)
	{
		pc[0][j] = 0;
		pb[0][j] = 0;
	}
	/************DP*****************************/
	for(i=1; i<nrow; i++)
	{
		for(j=1; j<ncolumn; j++)
		{
			if(s1[i-1] == s2[j-1])
			{
				pc[i][j] = pc[i-1][j-1] + 1;
				pb[i][j] = 1;
			}
			else if(pc[i-1][j] >= pc[i][j-1])
			{
				pc[i][j] = pc[i-1][j];
				pb[i][j] = 2;
			}
			else
			{
				pc[i][j] = pc[i][j-1];
				pb[i][j] = 3;
			}
		}
	}
}
/*create matrix*/
Matrix CreateMatrix(int nrow,int ncolumn)
{
	Matrix p;
	int i;
	p = (Matrix)malloc(nrow*sizeof(short *));
	if(p==NULL)
		return NULL;
	for(i=0; i<nrow; i++)
	{
		p[i] = (short *)malloc(ncolumn*sizeof(short));
	}
	return p;
}
/*free matrix*/
void DeleteMatrix(Matrix p,int nrow,int ncolumn)
{
	int n;
	if(p!=NULL && nrow>0 && ncolumn>0)
	{
		for(n=0; n<nrow; n++)
		{
			free(p[n]);
		}
		free(p);
	}
}
/*track back the matrix*/
void TrackBack(Matrix pc,Matrix pb,int nrow,int ncolumn)
{
	int ntemp;
	if(nrow == 0 || ncolumn == 0)
		return;
	ntemp = pb[nrow-1][ncolumn-1];
	pc[nrow-1][ncolumn-1] = -1;
	switch(ntemp)
	{
	case 1:		
		TrackBack(pc, pb, nrow-1, ncolumn-1);
		break;
	case 2:
		TrackBack(pc, pb, nrow-1, ncolumn);
		break;
	case 3:
		TrackBack(pc, pb, nrow, ncolumn-1);
		break;
	default:
		break;
	}
}



