#ifndef _LCS_H_
#define _LCS_H_

#ifndef NULL
#define NULL 0
#endif

typedef short **Matrix;

/*struct dyStack *maxRecord;*/
void GetLCSLen(int *s1,int *s2,Matrix pc,Matrix pb,int nrow,int ncolumn);
Matrix CreateMatrix(int nrow,int ncolumn);
void DeleteMatrix(Matrix p,int nrow,int ncolumn);
void TrackBack(Matrix pc,Matrix pb,int nrow,int ncolumn);

#endif
