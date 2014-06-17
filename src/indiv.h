#ifndef __INDIV_H__
#define __INDIV_H__

#include <iostream>
#include <string>
#include <vector> 
#include "fpmath.h"
using namespace std;

#define 	QQ		(9)

#if defined (IMPUTATION)
class Mask
{
public:
	int pos;                            //masked position;
	char snp;                   	//original genotype;
	Mask() {pos = -1;}; 
//	friend class Individual;
};
#endif

//the parent class of diploid and haploid; 
class Individual
{
protected:
	int   m_ploid; 
	int   m_isPanel; 
	char * snpGT; 					    // 1 dimension array;
	string		huskyID;           		// name or ID of the individual;
    double * 		phiScale;               // the scale to prevent phi from underflow; 
 	double * 		betaScale;              // the scale to prevent beta from underflow; 
 	real ** 	m_expectJmk; 			// MxK, expect jumps at locus m of cluster k; 
	real 		logLikelihood; 			// likelihood of genotype = log(m_phi(nLoci-1)); 
	int 		popLabel; 
	int         nMissingGT; 
#if defined (IMPUTATION)
	class Mask * pMask; 				// 100 * 20%; 
	real * maf;                 		//monir allele freq. 
#endif 
public:
	Individual();					   	// construct function;
	virtual ~Individual();              // destruct function; 
	
    inline void SetID(char * str){huskyID.assign(str);}
	// set husky ID for each individual if available;
	inline void AllocatesnpGT(int len)
		{snpGT = (char*) Allocate1D(sizeof(char), len);}
	// allocate snpGT for pointer; 
    inline void SetsnpGT(int pos, char val) {snpGT[pos] = val;}
   	// for readdipdata;
 	inline string GetID(void){return huskyID;}
	inline short GetsnpGT(int pos) {return short(snpGT[pos]-'0');}
	inline void GetsnpGT(int nLoci, short * snp_gt) 
	{ 
		for (int m = 0; m < nLoci; m++)
			snp_gt[m] = (short) (snpGT[m] - '0'); 
	}
	//		memcpy(snp_gt, snpGT, sizeof(char) * nLoci);}
	inline real GetLikelihood(void){return logLikelihood;}

	inline real ** GetexpectJmkptr(void) {return m_expectJmk;}
	inline void SetpopLabel(int s){popLabel = s;}
	inline int  GetpopLabel(void){return popLabel;}
	inline void SetnMissingGTpp(void) {nMissingGT++;}
	inline int  GetnMissingGT(void) {return nMissingGT;}

	inline void Setploid(int s) {m_ploid = s;}
	inline void SetisPanel(int s) {m_isPanel = s;}
	inline int Getploid(void) {return m_ploid;}
	inline int GetisPanel(void) {return m_isPanel;}
	// interface; 
public: 
	virtual int ** Getzpair(void) = 0; 
	//with joint_impute to get cluster membership; 

	virtual real get_snpmgt(int m) = 0; 
	
	virtual real ** Gettopptr(void) = 0;
	virtual real ** Getbotptr(void) = 0;
//	virtual real *** GetprZmptr(void) = 0;

	virtual void allocate_snp_mgt(int n) = 0; 
	virtual void allocate_snp_dstr(int n) = 0; 
	virtual void set_snp_dstr(int m, real p0, real p1) = 0; 
	virtual void set_snp_mgt(int m, real mg) = 0; 
	virtual void get_snp_dstr(int m, real * r) = 0;
	virtual real * get_snp_dstr(void) = 0;
	virtual real * get_snp_mgt(void) = 0;
	virtual void norm_snp_dstr(int, int) = 0;
	virtual void calc_snp_dstr(int nLoci, int nK, real ** theta) = 0; 
	
	virtual void CalcAll(int, int, class ModelParam *, int) = 0;	
	virtual void FreeMemAfterEM(void) = 0;
	virtual void joint_imputation(class ModelParam *, int, short *, int, int, int, int *) = 0; 
	
#if defined (IMPUTATION)
	virtual void MaskSNPs(int, int, int *) = 0; 
	virtual void ImputeMaskedCounting(int, int*, int*) = 0; 
	virtual real * Getmaf(void) = 0; 
#endif
};
#endif
