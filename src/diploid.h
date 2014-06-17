#ifndef __DIPLOID_H__
#define __DIPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class DipInd : public Individual
{
private:
    real ** 	m_phi;  		// forward probability, for each invidual, this is an MxKxK matrix; 
	real ** 	m_beta; 		// backward probability.  
	real ** 	m_top; 			// MxK, the top block of C1 formula in appendix C;  
	real ** 	m_bot; 
	real **     m_prZ; 
	int ** 	    zpair;
	real ** 	phiz;
	real *     snp_dstr;       // snp distribution;   nx2; 
	real *     mgt;            // snp mean genotype;  nx1; 
public:
	DipInd();
	~DipInd();

	real get_snpmgt(int m);  
	void get_snp_dstr(int m, real * r) {r[0] = snp_dstr[2*m]; r[1] = snp_dstr[2*m+1];}
	void allocate_snp_mgt(int n) {mgt = (real *) Allocate1D(sizeof(real),  n); }
	void allocate_snp_dstr(int n) {snp_dstr = (real *) Allocate1D(sizeof(real),  n * 2); }
	void set_snp_dstr(int m, real p0, real p1) {snp_dstr[2*m] = p0; snp_dstr[2*m+1]=p1; }
	void set_snp_mgt(int m, real s) {mgt[m] = s;}
	real * get_snp_dstr(void) {return snp_dstr;}
	real * get_snp_mgt(void) {return mgt;}
	void norm_snp_dstr(int, int);
	void calc_snp_dstr(int nLoci, int nK, real ** theta); 
	
	//for cluster-based association; working with joint_impute to get cluster membership; 
	virtual real ** Gettopptr(void) {return m_top;}
	virtual real ** Getbotptr(void) {return m_bot;}
//	virtual real ** GetprZmptr(void) {return m_prZ;}       //m_phi = m_prZ; 
	virtual int ** Getzpair(void) {return zpair;}

	virtual void joint_imputation(class ModelParam *, int, short *, int, int, int, int *);
	virtual void CalcAll(int, int, class ModelParam *, int);	
	virtual void FreeMemAfterEM(void);
	
#if defined (IMPUTATION)
	virtual void MaskSNPs(int, int, int *);      
	virtual void ImputeMaskedCounting(int, int*, int*); 
	virtual real * Getmaf(void) {return maf;}
#endif 
};

#endif
