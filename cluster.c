/************************************************************************/
/* Author:  Zhenjia Wang<zhenjia.sdu@gmail.com> Dec. 22, 2014
 * Biclustering procedure, greedy heuristic by picking an edge with highest
 * score and then dynamically adding vertices into the block and see if 
 * the block score can be improved.
 */

#include "cluster.h"

/************************************************************************/

static void update_colcand(bool *colcand, const char *s1)
{	
	int i=0;

	for (i=0; i< cols; i++)
		if (colcand[i] &&(s1[i]==0))
			colcand[i] = FALSE;
}


/**********************************************************/
static int intersect_row_LCS(const bool *colcand, char *colcand2)
/*caculate the weight of the edge*/
{
	int i;
	int cn = 0;
	for (i=0; i< cols; i++)
	{	if (colcand[i] && colcand2[i]!=0) 
			cn++;
	}
	return cn;
}

/***********************************************************/
static int seed_current(struct dyStack *genes, bool *colcand, int components,int *colsStat)
/* calculate the coverage of any row to the current consensus
 * cnt = # of valid consensus columns
 */
{
	int i,j;
	int cnt =0;
	int threshold = floor(components * 0.7)-1;
	if(threshold <1)
		threshold=1;
	/*get the statistical results of each column produced by seed*/
	char *temptag;
	AllocArray(temptag,cols);
	for(i=0;i<cols;i++)
	{	
		colsStat[i] = 0;
		temptag[i] = 0;
	}
	for(i=1;i<components;i++)
	{
		get_Genes_LCS(arr_c[dsItem(genes,0)], arr_c[dsItem(genes,i)],temptag);
		for(j=0;j<cols;j++)
		{
			if(temptag[j]!=0)
				colsStat[j]++;
			temptag[j]=0;
		}	
	}
/*	printf("\n");*/
	for(i=0;i<cols;i++)
	{	
		if (colsStat[i] >= threshold)
		{
			colcand[i] = TRUE; 
			cnt++;
		}
	}
	free(temptag);
	return cnt;
}

static bool check_seed(Edge *e, Block **bb, const int block_id)
/*check whether current edge can be treat as a seed*/
{
	int profiles[rows];
	int i,b1,b2,b3;
	bool fg = FALSE;
	b1 = b2 = -1;
	for (i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes,e->gene_one) && isInStack(bb[i]->genes, e->gene_two) ) 
			return FALSE; 

	for ( i = 0; i < rows; i++) profiles[i] = 0;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_one) ) 
		{ 
			fg = TRUE;
		       	break; 
		}
	if (fg) 
		b1 = i;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_two) ) 
		{ 
			fg = TRUE; 
			break; 
		}
	if (fg) 
		b2 = i;
	if ( (b1 == -1)||(b2 == -1) ) 
		return TRUE;
	else
	{
		for ( i = 0; i < bb[b1]->block_rows; i++)
			profiles[dsItem(bb[b1]->genes,i)]++;
		for ( i = 0; i < bb[b2]->block_rows; i++)
			profiles[dsItem(bb[b2]->genes,i)]++;
		for ( i = 0; i < rows; i++)
 			if (profiles[i] > 1) 
				return FALSE;
		b3 = MAX(bb[b1]->block_cols, bb[b2]->block_cols);
		if ( e->score <b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/ ) 
			return FALSE;
		else 
			return TRUE;
	}
	err("never see this message\n");
	return FALSE;
}

   
/************************************************/
static void block_init(Edge *e, Block *b, 
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster, long double *pvalues)
{
	int i,score,j;
	int cnt = 0, cnt_all=0, pid=0;
	continuous cnt_ave=0, row_all = rows;
	long double pvalue;
	int max_cnt, max_i;
	int t0,t1;
	int *arr_rows, *arr_rows_b;
	AllocArray(arr_rows, rows);
	AllocArray(arr_rows_b, rows);	
	bool *colcand;
	AllocArray(colcand, cols);
	for (i=0; i< cols; i++) 
		colcand[i] = FALSE;
	discrete *g1, *g2;
	t0=dsItem(genes,0);
	t1=dsItem(genes,1);
	g1 = arr_c[t0];
	g2 = arr_c[t1];
	
	for(i=0;i<rows;i++)	
	{
		lcs_length[i]=0;
		for(j=0;j<cols;j++)
			lcs_tags[i][j]=0;
	}
	
	/*************************calculate the lcs*********************************/
	
	lcs_length[t1]=get_Genes_LCS(g1,g2,lcs_tags[t1]); 
	for(i=0;i<cols;i++)
	{
		if(lcs_tags[t1][i]!=0)
			colcand[i]=TRUE;
	}

	for(j=0;j<rows;j++)
	{	
		if (j==t1 || j==t0) continue;
		lcs_length[j]= get_Genes_LCS_with_lcs(g1,arr_c[j],lcs_tags[j],lcs_tags[t1]); 
		
	}
/*	printf("\n==chose genes from MAXT:");*/
	while (*components < rows)
	{	
		max_cnt = -1;
		max_i = -1;
		(*components)++;
		cnt_all =0;
		cnt_ave = 0;
		/******************************************************/
		/*add a function of controling the bicluster by pvalue*/
		/******************************************************/
		for (i=0; i< rows; i++)
		{
			if (!candidates[i]) continue;
			if (po->IS_list && !sublist[i]) continue;
			cnt = intersect_row_LCS(colcand,lcs_tags[i]);
			cnt_all += cnt;
			if (cnt < cand_threshold) 
				candidates[i] = FALSE;
			if (cnt > max_cnt)
			{
				max_cnt = cnt;
				max_i = i;
			}
		}
	/**	printf("\t%d,%d",lcs_rowNo[max_i],max_cnt);*******************/

		cnt_ave = cnt_all/row_all;
		pvalue = get_pvalue (cnt_ave, max_cnt);
		if (po->IS_cond)
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0|| max_cnt < b->cond_low_bound) break;
		}
		else
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0) break;
		}


		if (po->IS_area)
			score = *components*max_cnt;
		else
			score = MIN(*components, max_cnt);
		if (score > b->score)
			b->score = score;
		if (pvalue < b->pvalue)
			b->pvalue = pvalue;
		dsPush(genes, max_i);
		dsPush(scores,score);
		pvalues[pid++] = pvalue;
		update_colcand(colcand, lcs_tags[max_i]);
		candidates[max_i] = FALSE;
	}
	/*be sure to free a pointer when you finish using it*/
	free(colcand);
	free(arr_rows);
	free(arr_rows_b);
}
/************************************************************************/
/* Core algorithm */
int cluster (FILE *fw, Edge **el, int n)
{
	int block_id = 0;
	Block **bb;
	int cnt = 0;
	int allocated = po->SCH_BLOCK;
	AllocArray(bb, allocated);
	Edge *e;
	Block *b;
	struct dyStack *genes, *scores, *b_genes, *allincluster;
	
	int i, j, k, components;
	
	int *colsStat;
	AllocArray(colsStat,cols);

	genes = dsNew(rows);
	scores = dsNew(rows);
	allincluster = dsNew(rows);

    
	long double *pvalues;
	AllocArray(pvalues, rows);

	bool *candidates;
	AllocArray(candidates, rows);
	e = *el; 
	

	AllocArray(lcs_length,rows);
	AllocArray(lcs_tags,rows);

	for(i=0;i<rows;i++)
	{	
		AllocArray(lcs_tags[i],cols);
	}	

	i = 0;
	while (i++ < n)
	{	
		e = *el++;
		/* check if both genes already enumerated in previous blocks */
		bool flag = TRUE;
		/* speed up the program if the rows bigger than 200 */
	        if (rows > 250)
		{ 
			if ( isInStack(allincluster,e->gene_one) && isInStack(allincluster,e->gene_two) )
				flag = FALSE;
			else if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			else if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag =FALSE;
		}
		else   
		{
			flag = check_seed(e, bb, block_id);
			if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag = FALSE;
		}
		/*if(e->gene_one > rows || e->gene_two > rows)
			flag=FALSE;*/
		if (!flag) continue;

		/*you must allocate a struct if you want to use the pointers related to it*/
		AllocVar(b);
		/*initial the b->score*/
                b->score = MIN(2, e->score);
		/*initial the b->pvalue*/
		b->pvalue = 1;
		/* initialize the stacks genes and scores */		
		int ii;		
		dsClear(genes);
		dsClear(scores);		
		for(ii = 0; ii < rows; ii ++)
		{
			dsPush(genes,-1);
			dsPush(scores,-1);
		}		
		dsClear(genes);
		dsClear(scores);
		
		dsPush(genes, e->gene_one);
		dsPush(genes, e->gene_two);

		dsPush(scores, 1);
		dsPush(scores, b->score);

		/* branch-and-cut condition for seed expansion */
		int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);
                if (cand_threshold < 2) 
			cand_threshold = 2;

		/* maintain a candidate list to avoid looping through all rows */		
		for (j = 0; j < rows; j++) 
			candidates[j] = TRUE;
		candidates[e->gene_one] = candidates[e->gene_two] = FALSE;
		components = 2;
		/* expansion step, generate a bicluster without noise */
		block_init(e, b, genes, scores, candidates, cand_threshold, &components, allincluster, pvalues);
		/* track back to find the genes by which we get the best score*/
		for(k = 0; k < components; k++)
		{
			if (po->IS_pvalue)
				if ((pvalues[k] == b->pvalue) &&(k >= 2) &&(dsItem(scores,k)!=dsItem(scores,k+1))) break;
			if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) break;
		}
		
		components = k + 1;
		int ki;
		for (ki=0; ki < rows; ki++)
		{	
			candidates[ki] = TRUE;
		}

		for (ki=0; ki < components - 1 ; ki++)
		{
			candidates[dsItem(genes,ki)] = FALSE;
		}

		candidates[dsItem(genes,k)] = FALSE;
		genes->top = k;
		
		bool *colcand;
		AllocArray(colcand, cols);
		for(ki = 0; ki < cols; ki++) 
			colcand[ki] = FALSE;             
    
		/* get init block */ 
		cnt = seed_current(genes, colcand, components,colsStat);
		/* add some new possible genes */
		int m_ct=0;
		bool colChose = TRUE;
		for(ki=0;ki < rows;ki++)
		{
			colChose=TRUE;
			if ((po->IS_list && !sublist[ki]) || !candidates[ki])
				continue;
			m_ct = intersect_row_LCS(colcand,lcs_tags[ki]);
			if (candidates[ki]&& (m_ct >= floor(cnt * po->TOLERANCE)-1))
			{
				int temp;
				for(temp=0;temp<cols;temp++)
				{
					if(colcand[temp])
					{	int tmpcount = colsStat[temp];
						if(lcs_tags[ki][temp]!=0)
							tmpcount++;
						if(tmpcount < floor(components * 0.1)-1)
						{	
							colChose = FALSE;
							break;
						}
					}
				}
				if(colChose==TRUE)			
				{	
					dsPush(genes,ki);
					components++;
					candidates[ki] = FALSE;
					for(temp=0;temp<cols;temp++)
					{
						if(lcs_tags[ki][temp]!=0 && colcand[ki])
						{
							colsStat[temp]++;
						}
					}
				}
			}
		}
                b->block_rows_pre = components;
		/* add genes that negative regulated to the consensus */
		char * reve_tag;
		for ( ki = 0; ki < rows; ki++)
		{
			colChose=TRUE;
			if (po->IS_list && !sublist[ki] && !candidates[ki]) continue;
				
			AllocArray(reve_tag,cols);
			int kk=0;
			for(kk=0;kk<cols;kk++)
				reve_tag[kk]=0;
			if(str_intersect_r(arr_c[dsItem(genes,0)],arr_c[ki],'N') < floor(cnt * po->TOLERANCE))
			{
				candidates[ki] = FALSE;
				continue;
			}
			get_Genes_LCS_R_with_lcs(arr_c[dsItem(genes,0)],arr_c[ki],reve_tag,lcs_tags[dsItem(genes,1)]);
			m_ct = intersect_row_LCS(colcand,reve_tag);
			if (candidates[ki] && (m_ct >= floor(cnt * po->TOLERANCE)-1))
			{
				int temp;
				for(temp=0;temp<cols;temp++)
				{
					if(colcand[temp])
					{	int tmpcount = colsStat[temp];
						if(reve_tag[temp]!=0)
							tmpcount++;
						if(tmpcount < floor(components * 0.1)-1)
						{	
							colChose = FALSE;
							break;
						}
					}
				}
				if(colChose == TRUE)			
				{	
					dsPush(genes,ki);
					components++;
					candidates[ki] = FALSE;
					for(temp=0;temp<cols;temp++)
					{
						if(reve_tag[temp]!=0 && colcand[ki])
						{
							colsStat[temp]++;
						}
					}
				}
			}
			free(reve_tag);
			reve_tag=NULL;
		}
		/* save the current cluster*/
		b_genes = dsNew(b->block_rows_pre);
		for (ki = 0; ki < b->block_rows_pre; ki++)
			dsPush(b_genes, dsItem(genes,ki));
		/* store gene arrays inside block */
		b->genes = dsNew(components);
		b->conds = dsNew(cols);
	
		scan_block(colcand, b);
		free(colcand);
		colcand=NULL;
		if (b->block_cols < 4 || components < 5) continue;
		b->block_rows = components;
                if (po->IS_pvalue)
			b->score = -(100*log(b->pvalue));
		else
			b->score = b->block_rows * b->block_cols;		

		dsClear(b->genes);
		for ( ki=0; ki < components; ki++)
			dsPush(b->genes,dsItem(genes,ki));
		for(ki = 0; ki < components; ki++)
			if(!isInStack(allincluster, dsItem(genes,ki))) 
				dsPush(allincluster,dsItem(genes,ki));	
		/*save the current block b to the block list bb so that we can sort the blocks by their score*/
		bb[block_id++] = b;

		/* reaching the results number limit */
		if (block_id == po->SCH_BLOCK) break;
		verboseDot();	
	}
	/* writes character to the current position in the standard output (stdout) and advances the internal file position indicator to the next position.
	 * It is equivalent to putc(character,stdout).*/
	putchar('\n');
	/* free-up the candidate list */
	free(candidates);
	free(allincluster);
	free (pvalues);
	free(colsStat);
	return report_blocks(fw, bb, block_id);
}

/************************************************************************/
static void print_params(FILE *fw)
{
	char filedesc[LABEL_LEN];
	strcpy(filedesc, "continuous");
	if (po->IS_DISCRETE) 
		strcpy(filedesc, "discrete");
	fprintf(fw, "# Unibic version %.1f output\n", VER);
	fprintf(fw, "# Datafile %s: %s type\n", po->FN, filedesc);
	fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
			po->COL_WIDTH, po->FILTER, po->TOLERANCE, po->RPT_BLOCK);
	if (!po->IS_DISCRETE) 
		fprintf(fw, " -q %.2f -r %d", po->QUANTILE, po->DIVIDED);
	fprintf(fw, "\n\n");
}

/************************************************************************/
static int report_blocks(FILE* fw, Block** bb, int num)
{
	print_params(fw);
	sort_block_list(bb, num);
	
	int i, j,k;
	/*MIN MAX et al functions can be accessed in struct.h*/
        int n = MIN(num, po->RPT_BLOCK);
	bool flag;

	Block **output;
	AllocArray(output, n);

	Block **bb_ptr = output;
	Block *b_ptr;
	double cur_rows, cur_cols;
	double inter_rows, inter_cols;
        /*double proportion;*/
	
	/* the major post-processing here, filter overlapping blocks*/
	i = 0; j = 0;
	while (i < num && j < n)
	{
		b_ptr = bb[i];
		cur_rows = b_ptr->block_rows;
		cur_cols = b_ptr->block_cols;

		flag = TRUE;
		k = 0;
		while (k < j)
		{
			inter_rows = dsIntersect(output[k]->genes, b_ptr->genes);
			inter_cols = dsIntersect(output[k]->conds, b_ptr->conds);
			
			if (inter_rows*inter_cols > po->FILTER*cur_rows*cur_cols)
			{
				flag = FALSE; 
				break;
			}
                        k++;
		}
	        i++;
		if (flag)
		{
			print_bc(fw, b_ptr, j++);
			*bb_ptr++ = b_ptr;
		}
	}
	return j;
}
/************************************************************************/

static int block_cmpr(const void *a, const void *b)
/* compare function for qsort, descending by score */
{
	return ((*(Block **)b)->score - (*(Block **)a)->score);
}

static void sort_block_list(Block **el, int n)
{
	qsort(el, n, sizeof *el, block_cmpr);
}
/************************************************************************/

long double get_pvalue (continuous a, int b)
{
	int i =0;
	long double one = 1, pvalue=0;
	long double poisson=one/exp(a);
	for (i=0;i<b+300;i++)
	{
		if (i>(b-1)) 
			pvalue=pvalue+poisson;
		else 
			poisson=poisson*a/(i+1);
	}
	return pvalue;
}

/**************************************************************************/
int get_Genes_LCS(const discrete *s1, const discrete *s2,char *lcs_tg)
{
	struct dyStack *maxRecord;/*record the max value of matrix*/
	int maxvalue,rank,i,j,length1,length2;
	int *temp1,*temp2;
	Matrix C,B;

	maxvalue = 0;
	length1=length2=0;
	rank = po ->DIVIDED;
	maxRecord = dsNew(cols);
	AllocArray(temp1,cols);
	AllocArray(temp2,cols);
	for(i=0;i<cols;i++)
	{
		temp2[i]=0;
		temp1[i]=0;
	}
	/*get the sorted sequence*/
	for(i=1;i<=rank;i++)
		for(j=0;j<cols;j++)
		{
			if(s1[j] == i)
				temp1[length1++]=j+1;
			/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
			if(s2[j] == i)
				temp2[length2++]=j+1;
	
		}
	if(po->dataMode==0 && po->QUANTILE < 0.5)
	{
		for(i=rank*(-1);i<=-1;i++)
			for(j=0;j<cols;j++)
			{
				if(s1[j] == i)
					temp1[length1++]=j+1;
				/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
				if(s2[j] == i)
					temp2[length2++]=j+1;
		
			}
	}
	/*create matrix for lcs*/
	C = CreateMatrix(length1 + 1,length2 + 1);
	B = CreateMatrix(length1 + 1,length2 + 1);
	
	/*start to get the LCS of temp1 and temp2*/
	GetLCSLen(temp1,temp2,C,B,length1 + 1,length2 + 1);

	/*put the max value of Matrix C into stack*/
	maxvalue = C[length1][length2];

	
	for (j=1;j<length2+1;j++)
	{
		if (C[length1][j] == maxvalue)
			dsPush(maxRecord,j);
	}
	
	/*find all the columns of all LCSs*/
	for (i=maxRecord->top;i>=0;i--)
	{
		TrackBack(C,B, length1+1, dsItem(maxRecord,i)+1);break;
	}
	for (i=1;i<length1+1;i++)
	{
		for (j=1;j<length2+1;j++)
		{
			if (C[i][j] == -1 && B[i][j]==1)
			{
				/*printf("lcs_tg:%d,%d\n",temp1[i-1],lcs_tg[temp1[i-1]-1]);*/
				lcs_tg[temp1[i-1]-1] = 1;
			}
		}
	}
	
	dsFree(maxRecord);
	free(temp1);
	free(temp2);
	for(i=0;i<=length1;i++)
	{
		free(C[i]);
		free(B[i]);
	}
	free(C);
	free(B);
	
	return maxvalue;

}
/**************************************************************************/
int get_Genes_LCS_length(const discrete *s1, const discrete *s2)
{
	int maxvalue,rank,i,j,length1,length2;
	int *temp1,*temp2;
	Matrix C,B;

	maxvalue = 0;
	length1=length2=0;
	rank = po ->DIVIDED;
	AllocArray(temp1,cols);
	AllocArray(temp2,cols);
	for(i=0;i<cols;i++)
	{
		temp2[i]=0;
		temp1[i]=0;
	}
	/*get the sorted sequence*/
	for(i=1;i<=rank;i++)
		for(j=0;j<cols;j++)
		{
			if(s1[j] == i)
				temp1[length1++]=j+1;
			/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
			if(s2[j] == i)
				temp2[length2++]=j+1;
	
		}
	if(po->dataMode==0 && po->QUANTILE < 0.5)
	{
		for(i=rank*(-1);i<0;i++)
			for(j=0;j<cols;j++)
			{
				if(s1[j] == i)
					temp1[length1++]=j+1;
				if(s2[j] == i)
					temp2[length2++]=j+1;		
			}
	}
	/*create matrix for lcs*/
	C = CreateMatrix(length1 + 1,length2 + 1);
	B = CreateMatrix(length1 + 1,length2 + 1);
	
	/*start to get the LCS of temp1 and temp2*/
	GetLCSLen(temp1,temp2,C,B,length1 + 1,length2 + 1);

	/*put the max value of Matrix C into stack*/
	maxvalue = C[length1][length2];

	free(temp1);
	free(temp2);
	for(i=0;i<=length1;i++)
	{
		free(C[i]);
		free(B[i]);
	}
	free(C);
	free(B);
	
	return maxvalue;

}

/**************************************************************************/
int get_Genes_LCS_with_lcs(const discrete *s1, const discrete *s2,char *lcs_tg,char *lcs_seed)
{
	struct dyStack *maxRecord;/*record the max value of matrix*/
	int maxvalue,rank,i,j,length1,length2;
	int *temp1,*temp2;
	Matrix C,B;

	maxvalue = 0;
	length1=length2=0;
	rank = po ->DIVIDED;
	maxRecord = dsNew(cols);
	AllocArray(temp1,cols);
	AllocArray(temp2,cols);
	for(i=0;i<cols;i++)
	{
		temp2[i]=0;
		temp1[i]=0;
	}
	/*get the sorted sequence*/
	for(i=1;i<=rank;i++)
		for(j=0;j<cols;j++)
		{	
			if(lcs_seed[j] != 0)
			{
				if(s1[j] == i)
					temp1[length1++]=j+1;
				/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
				if(s2[j] == i)
					temp2[length2++]=j+1;
			}
	
		}
	if(po->dataMode==0 && po->QUANTILE < 0.5)
	{
		for(i=rank*(-1);i<=-1;i++)
			for(j=0;j<cols;j++)
			{	
				if(lcs_seed[j] != 0)
				{
					if(s1[j] == i)
						temp1[length1++]=j+1;
					/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
					if(s2[j] == i)
						temp2[length2++]=j+1;
				}
			}
	}
	/*create matrix for lcs*/
	C = CreateMatrix(length1 + 1,length2 + 1);
	B = CreateMatrix(length1 + 1,length2 + 1);
	
	/*start to get the LCS of temp1 and temp2*/
	GetLCSLen(temp1,temp2,C,B,length1 + 1,length2 + 1);

	/*put the max value of Matrix C into stack*/
	maxvalue = C[length1][length2];

	
	for (j=1;j<length2+1;j++)
	{
		if (C[length1][j] == maxvalue)
			dsPush(maxRecord,j);
	}
	

	/*find all the columns of all LCSs*/
	for (i=maxRecord->top;i>=0;i--)
	{
		TrackBack(C,B, length1+1, dsItem(maxRecord,i)+1);break;
	}
	for (i=1;i<length1+1;i++)
	{
		for (j=1;j<length2+1;j++)
		{
			if (C[i][j] == -1 && B[i][j]==1)
			{
				/*printf("lcs_tg:%d,%d\n",temp1[i-1],lcs_tg[temp1[i-1]-1]);*/
				lcs_tg[temp1[i-1]-1] = 1;
			}
		}
	}
	
	dsFree(maxRecord);
	free(temp1);
	free(temp2);
	for(i=0;i<=length1;i++)
	{
		free(C[i]);
		free(B[i]);
	}
	free(C);
	free(B);
	
	return maxvalue;

}

/**************************************************************************/
int get_Genes_LCS_R_with_lcs(const discrete *s1,const discrete *s2,char *lcs_tagR,char *lcs_seed)
{
	struct dyStack *maxRecord;/*record the max value of matrix*/
	int maxvalue,rank,i,j,length1,length2;
	int *temp1,*temp2;
	Matrix C,B;

	maxvalue = 0;
	length1=length2=0;
	rank = po ->DIVIDED;
	maxRecord = dsNew(cols);
	AllocArray(temp1,cols);
	AllocArray(temp2,cols);

	/*get the sorted sequence*/
	if(po->dataMode==0 && po->QUANTILE < 0.5)
	{
		for(i=1;i<=rank;i++)
			for(j=0;j<cols;j++)
			{	
				if(lcs_seed[j] != 0)
				{
					if(s1[j] == i)
						temp1[length1++]=j+1;
					/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
					if(s2[j] == i*(-1))
						temp2[length2++]=j+1;
				}
			}
		
		for(i=rank*(-1);i<=-1;i++)
			for(j=0;j<cols;j++)
			{	
				if(lcs_seed[j] != 0)
				{
					if(s1[j] == i)
						temp1[length1++]=j+1;
					/*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
					if(s2[j] == i*(-1))
						temp2[length2++]=j+1;
				}
			}
	}
	else
	{
		for(i=1;i<=rank;i++)
			for(j=0;j<cols;j++)
			{
				if(lcs_seed[j]!=0)
				{
					if(s1[j] == i)
						temp1[length1++]=j+1;
					if(s2[j] == rank-i+1)
						temp2[length2++]=j+1;
				}
			}
	}
	/*create matrix for lcs*/
	C = CreateMatrix(length1 + 1,length2 + 1);
	B = CreateMatrix(length1 + 1,length2 + 1);
	/*start to get the LCS of temp1 and temp2*/
	GetLCSLen(temp1,temp2,C,B,length1 + 1,length2 + 1);
	/*put the max value of Matrix C into stack*/
	maxvalue = C[length1][length2];
	
	for (j=1;j<length2+1;j++)
	{
		if (C[length1][j] == maxvalue)
			dsPush(maxRecord,j);
	}
	/*find all the columns of all LCSs*/
	for (i=maxRecord->top;i>=0;i--)
	{
		TrackBack(C,B, length1+1, dsItem(maxRecord,i)+1);break;
	}

	for (i=1;i<length1+1;i++)
	{
		for (j=1;j<length2+1;j++)
		{
			if (C[i][j] == -1 && B[i][j]==1)
			{
				lcs_tagR[temp1[i-1]-1] = 1;
			}
		}
	}
	dsFree(maxRecord);
	free(temp1);
	free(temp2);
	for(i=0;i<=length1;i++)
	{
		free(C[i]);
		free(B[i]);
	}
	free(C);
	free(B);
	return maxvalue;

}
