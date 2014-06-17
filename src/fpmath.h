#ifndef __FPMATH_H__
#define __FPMATH_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_randist.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_permute.h"
#include "fp.h"

#define NA (54321)
#define invsr2pi 0.3989423
#define TINY 1e-100
//#define isnan(x) ((x) != (x))   //for some compiler use this macro for isnan
using namespace std;

typedef float real;
//globle random number generator; 
extern const gsl_rng_type * gslType;
extern gsl_rng * gsl_r;

class Pair {
	public:	
		real bf; 
		int pos;
		int label; 
};
double sumlog(double *, int); 
int geometric(int len, int meanrank); 
double pgeometric(int pick, int len, int meanrank); 
int imap(int, int, int); 
string getline(streambuf * ); 
void DistinctIntArray(int, int, int, int *); 		// generate arg3 many distinct intergers on [arg1, arg2); 
//void lu_decomp(real ** a, int n, int *indx, int *d);
//void lu_back_sub(real ** a, int n, int *indx, real b[]);
void safe_exit(); 
int compare_pair(const void * a, const void * b); 
int compare(const void * a, const void * b);

void * Allocate1D(size_t, int);
void ** Allocate2D(size_t, int, int);  
void *** Allocate3D(size_t, int, int, int);  
void Free1D(void *);
void Free2D(void **); 
void Free3D(void ***); 

real ** 	Allocate2DMatrix(int, int);
int ** 		Allocate2DIntMatrix(int, int);
short int** Allocate2DShortMatrix(int, int);
char ** 	Allocate2DCharMatrix(int, int);
real *** 	Allocate3DMatrix(int, int, int);
void 	Free2DCharMatrix(char **);
void 	Free2DIntMatrix(int **);
void 	Free2DShortMatrix(short int **);
void 	Free2DMatrix(real **);
void 	Free3DMatrix(real ***);
 
void center_col(double ** xx, int nrow, int ncol); 
void center_col(real ** xx, int nrow, int ncol); 
int nchoosek(int,int); 

#endif
