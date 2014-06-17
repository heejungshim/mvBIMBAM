#include "fpmath.h"
#include "string.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

using namespace std;

double sumlog(double * ta, int n)
{
	double dmax = ta[0]; 
	for (int i = 1; i < n; i++)
	{
        if(dmax < ta[i]) 
			dmax = ta[i]; 
	}

	double dmul = 0; 
	for (int i = 0; i < n; i++)
	{
		double diff = ta[i] - dmax; 
        dmul += exp(diff); 
	}
	return(dmax + log(dmul)); 
}

int geometric(int len, int meanrank)
{
	if(meanrank > len) 
		meanrank = len;
	double tmin = exp(-(double) len / (double) meanrank); 
	double x = gsl_rng_uniform(gsl_r) * (1.0 - tmin) + tmin; 
	return ((int) floor(-log(x) * meanrank));  
}

double pgeometric(int pick, int len, int meanrank)
{
	if(meanrank > len) 
		meanrank = len;   
	double p = 1.0 - exp(-1.0 / (double) meanrank);
	double q = 1.0 - p;

	double temp = len * log(q); 
	if (temp < -100) 
		temp = 0; 
	else 
		temp = exp(temp); 
	double logp = pick * log(q) + log(p) - log(1.0 - temp); 
	return(logp); 
}

string getline(streambuf * pbuf)
{
    char ch;
	string str;
	size_t pos; 
	while((ch = pbuf->sgetc()) != EOF)
	{
		if(ch != '\n' && ch != '\r')
		{
			str.push_back(ch);
			ch = pbuf->snextc();
		}
		else {
			pbuf->sbumpc();  //chomp;
			pos = str.find_first_not_of(";, \t", 0); 
			if(str.empty() || pos  == string::npos || str.at(pos) == '#')
			{
				str.clear();
				continue;
			}
			else
				break;
		}
	}
	return str;
}   //this getline use ;, \t as delimit, and ignore the lines either full of delimit or starting with #. 

int imap(int nK, int x, int y)
{
	int tmax = max(x,y); 
	int tmin = min(x,y); 
	int res = (tmin) * nK - (tmin * tmin - tmin)/2; 
	res += (tmax - tmin);  
	return (res); 
}

void safe_exit() 
{
#if defined (MPI_ENABLED)
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
#endif
	exit(0); 
}

int compare_pair(const void * a, const void * b)
{
	Pair * ta = (Pair *) a;
	Pair * tb = (Pair *) b; 
	if((tb->bf) > (ta->bf)) return 1;
	else if((tb->bf) == (ta->bf)) return 0;
	else return -1;
}

int compare(const void * a, const void * b)
{
	int *ta = (int *) a;
	int *tb = (int *) b;
	if ((*ta) > (*tb)) return 1;
	else if ((*ta) == (*tb)) return 0;
	else return -1;
}

void * Allocate1D(size_t us, int dim)
{
	size_t size = dim * us; 
	char * m = (char *) malloc(size);
	memset(m, 0, size); 
	return (void *) m; 
}

void Free1D(void * m)
{
	if(m == NULL) return; 
	free(m); 
	m = NULL; 
}

void ** Allocate2D(size_t us,  int dim1, int dim2)
{
	char ** m;
	m = (char **) malloc((size_t)(dim1 * us));
	m[0] = (char *) malloc((size_t)(dim1 *dim2 * us));
   	memset(m[0], 0, (size_t) (dim1 * dim2 * us));  
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D real matrix. \n");
		safe_exit();
	}
	for(int i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + (size_t) (dim2 * us);
	}
	return ((void **) (m));
}

void Free2D(void ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

real ** Allocate2DMatrix(int dim1, int dim2)
{
	int i;
	real ** m;
	
	m = (real **) malloc((size_t)((dim1)*sizeof(real*)));
	m[0] = (real *) malloc((size_t)((dim1*dim2)*sizeof(real)));
	memset(m[0], 0, (dim1*dim2)*sizeof(real)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D real matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

//
//void *** Allocate3DMatrix(size_t us, int dim1, int dim2, int dim3)
//{
//	void *** m; 
//	m = (void ***) malloc((size_t)(dim1 * us));
//	m[0] = (void **) malloc((size_t)(dim1 * dim2 * us));
//	m[0][0] = (void *) malloc((size_t)(dim1 * dim2 * dim3 * us));
//	if (!(m && m[0] &&  m[0][0]))
//	{
//		printf("Error: Problem allocating a 3D real matrix. \n");
//		safe_exit();
//	}
//
//	for (int j = 1; j < dim2; j++)
//	{
//		m[0][j] = m[0][j-1] + dim3;
//	}
//	for (int i = 1; i < dim1; i++)
//	{
//		m[i] = m[i-1] + dim2;
//		m[i][0] = m[i-1][dim2 - 1] + dim3;
//		for(int j = 1; j < dim2; j++)
//		{
//			m[i][j] = m[i][j-1] + dim3;
//		}
//	}
//	return (m);
//}
//
//void Free3DMatrix(void *** m)
//{
//	if(m == NULL) return;
//	free(m[0][0]);
//	free(m[0]);
//	free(m);
//	m = NULL; 
//}

real *** Allocate3DMatrix(int dim1, int dim2, int dim3)
{
	int i, j;
	real *** m; 

	m = (real ***) malloc((size_t)((dim1)*sizeof(real**)));
	m[0] = (real **) malloc((size_t)((dim1) * (dim2) * sizeof(real *)));
	m[0][0] = (real *) malloc((size_t)((dim1) * (dim2) * (dim3) * sizeof(real)));
	memset(m[0][0], 0, (dim1) * (dim2) * (dim3) * sizeof(real));  
	if (!(m && m[0] &&  m[0][0]))
	{
		printf("Error: Problem allocating a 3D real matrix. \n");
		safe_exit();
	}

	for (j = 1; j < dim2; j++)
	{
		m[0][j] = m[0][j-1] + dim3;
	}
	for (i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
		m[i][0] = m[i-1][dim2 - 1] + dim3;
		for(j = 1; j < dim2; j++)
		{
			m[i][j] = m[i][j-1] + dim3;
		}
	}
	return (m);
}

void Free2DMatrix(real ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free3DMatrix(real *** m)
{
	if(m == NULL) return;
	free(m[0][0]);
	free(m[0]);
	free(m);
	m = NULL; 
}
 

int ** Allocate2DIntMatrix(int dim1, int dim2)
{
	int i;
	int ** m;
	
	m = (int **) malloc((size_t)((dim1)*sizeof(int*)));
	m[0] = (int *) malloc((size_t)((dim1*dim2)*sizeof(int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(int)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

short int ** Allocate2DShortMatrix(int dim1, int dim2)
{
	int i;
	short int ** m;
	
	m = (short int **) malloc((size_t)((dim1)*sizeof(short int*)));
	m[0] = (short int *) malloc((size_t)((dim1*dim2)*sizeof(short int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(short)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    

void Free2DShortMatrix(short int ** m)
{
	if(m == NULL) return;
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free2DIntMatrix(int ** m)
{                        
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

char  ** Allocate2DCharMatrix(int dim1, int dim2)
{
	int i;
	char ** m;
	
	m = (char **) malloc((size_t)((dim1)*sizeof(char*)));
	m[0] = (char *) malloc((size_t)((dim1*dim2)*sizeof(char)));
	memset(m[0], 0, (dim1*dim2)*sizeof(char)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D char matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    


void Free2DCharMatrix(char ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

void DistinctIntArray(int low, int high, int dim, int * A)
{
    if (dim >= high - low)
		return;  
	
	int i, j, r; 
	int bingle; 
	for (i=0; i<dim; i++)
		A[i] = -1;
	
	int howmany = 0;
	for (i=high - dim; i<high; i++)
	{
		bingle = 0;
		r = low + gsl_rng_uniform_int(gsl_r, i-low);
		for (j = 0; j < howmany; j++)
		{
			if (r == A[j])
			{
				bingle = 1; 
				break;
			}
		}

		if (bingle) A[howmany] = i;
		else A[howmany] = r; 
		howmany++;
	}
}

//void lu_decomp(real ** a, int n, int *indx, int *d)
//{
//	int i, j, k;
//	int imax = -1;
//	real big, dum, sum, temp;
//	real *vv = new real[n];
//	*d = 1;
//
//	for (i = 0; i < n; i++)
//	{
//		big = 0.0; 
//		for (j = 0; j < n; j++)
//			if((temp = fabs(a[i][j])) > big) big = temp;
//		if (big == 0.0)
//		{
//			cout << "singular matrix in routine ludcmp" << endl;
//			exit(0); 
//		}
//		vv[i] = 1.0 / big; 
//	}
//
//	for (j = 0; j < n; j++)
//	{
//		for (i = 0; i < j; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//		}
//
//		big = 0.0;
//		for (i = j; i < n; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//			if ((dum = vv[i] * fabs(sum)) >= big)
//			{
//				big = dum; 
//				imax = i; 
//			}
//		}
//
//		if ( j != imax)
//		{
//			for (k = 0; k < n; k++)
//			{
//				dum = a[imax][k];
//				a[imax][k] = a[j][k];
//				a[j][k] = dum; 
//			}
//			*d = -(*d); 
//			vv[imax] = vv[j];
//		}
//		indx[j] = imax;
//		if(a[j][j] == 0.0) a[j][j] = TINY;
//		if (j != n)
//		{
//			dum = 1.0 / a[j][j];
//			for (i=j+1; i < n; i++) a[i][j] *= dum; 
//		}
//	}
//	delete[] vv; 
//}
//
//void lu_back_sub(real ** a, int n, int *indx, real b[])
//{
//	int i, ii = -1, ip, j; 
//	real sum;
//
//	for (i = 0; i < n; i++)
//	{
//		ip = indx[i];
//		sum = b[ip];
//		b[ip] = b[i];
//		if(ii>=0)
//			for (j=ii;j <i;j++) sum -= a[i][j] * b[j];
//		else if(sum) ii = i; 
//		b[i] = sum; 
//	}
//
//	for (i = n-1; i >= 0; i--) 
//	{
//		sum = b[i];
//		for (j=i+1; j < n; j++) sum -= a[i][j] * b[j];
//		b[i] = sum/a[i][i];
//	}
//} 

void center_col(double ** xx, int nrow, int ncol)
{
    for (int j = 0; j< ncol; j++)
    {
        double sum = 0;
        for (int i = 0; i < nrow; i++)
            sum += xx[i][j];
        sum /= (double) nrow;
        for (int i = 0; i < nrow; i++)
            xx[i][j] -= sum;
    }

}

void center_col(real ** xx, int nrow, int ncol)
{
    for (int j = 0; j< ncol; j++)
    {
        double sum = 0;
        for (int i = 0; i < nrow; i++)
            sum += xx[i][j];
        sum /= (double) nrow;
        for (int i = 0; i < nrow; i++)
            xx[i][j] -= sum;
    }

}

int nchoosek(int n, int k) 
{
  	if (k <= 0 || k >= n)
  		return 1;
  	else
  		return nchoosek(n-1, k-1) + nchoosek(n-1, k);
}

