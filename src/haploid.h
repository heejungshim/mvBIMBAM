#ifndef __HAPLOID_H__
#define __HAPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class HapInd : public Individual
{
private:
	real ** m_phi; 						// forward MxK
	real ** m_beta; 					// backward MxK
	real ** m_top; 
	real ** m_bot; 
	real ** m_prZ; 
	real *  snp_dstr; 	             //nx2; 
	real *  mgt; 
public:
	HapInd();
	~HapInd();
	
	real get_snpmgt(int m) {if(mgt) return mgt[m]; else if(snp_dstr) return snp_dstr[2*m+1]; else return (real) (snpGT[m]-'0'); }; 
	void get_snp_dstr(int m, real * r) { r[0] = snp_dstr[2*m];r[1] = snp_dstr[2*m+1];}
	void allocate_snp_mgt(int n) {mgt = (real *) Allocate1D(sizeof(real),  n); }
	void allocate_snp_dstr(int n) {snp_dstr = (real *) Allocate1D(sizeof(real),  n * 2); }
	void set_snp_dstr(int m, real p0, real p1) {snp_dstr[2*m] = p0; snp_dstr[2*m+1]=p1; }
	void set_snp_mgt(int m, real s) {mgt[m] = s; }
	real * get_snp_dstr(void) {return snp_dstr;}
	real * get_snp_mgt(void) {return mgt;}
	void norm_snp_dstr(int, int);
	void calc_snp_dstr(int nLoci, int nK, real ** theta); 

	virtual inline int** Getzpair(void) {return NULL;}
	virtual real ** Gettopptr(void) {return m_top;}
	virtual real ** Getbotptr(void) {return m_bot;}
//	virtual real *** GetprZmptr(void) {return &m_prZ;}   
	
	virtual void CalcAll(int, int, class ModelParam *, int);
	virtual void FreeMemAfterEM(void); 
	virtual void joint_imputation(class ModelParam *, int, short *, int, int, int, int *);
	
#if defined (IMPUTATION)
	virtual void MaskSNPs(int, int, int *);      
	virtual void ImputeMaskedCounting(int, int*, int*); 
	virtual real * Getmaf(void) {return maf;}
#endif 
};

#endif 
