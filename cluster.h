#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "struct.h"

/* this structure holds the matching score for each row */
struct rowMatch
{
    int row_id;
    int matches;
};

/* struct */
extern void verboseDot();

/* make_graph */
extern void seed_update (const discrete *s);

/* write_block */
extern void scan_block (bool * cand, Block *b_ptr);
extern void print_bc (FILE *fw, Block *b, int num);
extern int str_intersect_r(const discrete *s1, const discrete *s2,char tag);

/* prototypes */
static void update_colcand(bool *colcand, const char *s1);
static int intersect_row_LCS(const bool *colcand, char *colcand2);
static int seed_current(struct dyStack *genes, bool *colcand, int components,int *colsStat);
static bool check_seed(Edge *e, Block **bb, const int block_id);
static void print_params(FILE *fw);
int cluster (FILE *fw, Edge **el, int n);
static int report_blocks(FILE *fw, Block **bb, int num);
static void sort_block_list(Block **el, int n);
static void block_init(Edge *e, Block *b, struct dyStack *genes, struct dyStack *scores,
bool *candidates, const int cand_threshold,int *components, struct dyStack *allincluster, long double *pvalues);
long double get_pvalue (continuous a, int b);
/*get positive lcs of each seed,return the cosine similarity*/
int get_Genes_LCS(const discrete *s1, const discrete *s2,char *lcs_tg);
int get_Genes_LCS_with_lcs(const discrete *s1, const discrete *s2,char *lcs_tg,char *lcs_seed);
/*get the CSvalue of each seed*/
int get_Genes_LCS_R(const discrete *s1,const discrete *s2,char *lcs_ptrR);
int get_Genes_LCS_R_with_lcs(const discrete *s1,const discrete *s2,char *lcs_ptrR,char *lcs_seed);
int get_Genes_LCS_length(const discrete *s1, const discrete *s2);
bits16 **profile;

#endif
