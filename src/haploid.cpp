#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "haploid.h"

#include <iostream>
using namespace std;

                 
HapInd::HapInd(void)
{
	m_phi = NULL;
	m_beta = NULL;
	m_top = NULL;
	m_bot = NULL; 
	m_prZ = NULL; 
	snp_dstr = NULL; 
	mgt = NULL; 
	m_ploid = 1; 
}

HapInd::~HapInd(void)
{
	if(m_phi) {Free2DMatrix(m_phi);    m_phi = NULL; }
	if(m_beta) {Free2DMatrix(m_beta);       m_beta = NULL;}
	if(m_top) {Free2DMatrix(m_top); m_top = NULL;}
	if(m_bot) {Free2DMatrix(m_bot); m_bot = NULL;}
}

real prG(real t, short snp)
{
	real res = 0.0; 
	switch(snp)
	{
		case 0:
			res = 1.0 - t;
			break;
		case 1:
			res = t;
			break;
		case QQ: 
			res = 1.0;
			break;
		default:
			cout << "error in prG" << endl;
			break;
	}
	return res; 
}

void HapInd::CalcAll(const int nLoci, const int nK, ModelParam * pMD, const int ball)
{
	real * r = pMD->Getr(popLabel); 
	real ** theta = pMD->Gettheta();
	real ** alpha = pMD->Getalpha(popLabel);
	
	int k, m;

	/***************************************************************************  
	| calculate backward prob. m_beta for each individual. 
	| input: alpha.
	| output: m_beta. 
	***************************************************************************/ 
	if (phiScale == NULL)
		phiScale = new double[nLoci]; //(int *) malloc(sizeof(int) * nLoci); 
	if (betaScale == NULL)
		betaScale = new double[nLoci]; //(int *) malloc(sizeof(int) * nLoci); 

	real * tSum = new real[nLoci]; //(real *) malloc(sizeof(real) * nLoci);  
	
	if(m_beta == NULL)
		m_beta = Allocate2DMatrix(nLoci, nK); 
	
	//at the rightmost locus M; 
	for (k = 0; k < nK; k++)
		(m_beta)[nLoci-1][k] = 1.0;
	betaScale[nLoci-1] = 0;
	
	tSum[nLoci-1] = 0.0;
	for (k = 0; k < nK; k++)
		tSum[nLoci-1] += prG(theta[nLoci-1][k], GetsnpGT(nLoci-1)) * alpha[nLoci-1][k]; //beta[nLoci-1] = 1.0
	
	//do the recursive calc backwards; 
	for (m = nLoci - 1; m > 0; m--)
	{
		for (k = 0; k < nK; k++)
			{
				double temp = (1.0 - r[m]) * prG(theta[m][k], GetsnpGT(m)) * m_beta[m][k]; //1-r[m] = probJ[m][0]
				m_beta[m-1][k] = temp + r[m] * tSum[m];
			}
		
		//margin sum and real sum for locus m-1; 
		tSum[m-1] = 0.0;
		for (k = 0; k < nK; k++)
			tSum[m-1] += prG(theta[m-1][k], GetsnpGT(m-1)) * m_beta[m-1][k] * alpha[m-1][k];
		
		double tscale = 0;
		if (tSum[m-1] <= 1e-100)
		{
			cout << "sum of beta is 0!" << endl;
			exit(0);
		}
		else
			tscale = -log10(tSum[m-1]); 

		if(tscale > 2) 
		{
			betaScale[m-1] = betaScale[m] + tscale; 
			double dummy = pow(10., (double)tscale);
			tSum[m-1] *= dummy;
			for (k = 0; k < nK; k++)
				 m_beta[m-1][k] *= dummy;
			//rescale the partial sums; 
		}
		else 
			betaScale[m-1] = betaScale[m]; 

	}
	/* end calc beta, real checked. */ 
	
	/***************************************************************************  
	| calc forward prob. m_phi 
	| input: alpha, probJ  
	| output: m_phi 
	***************************************************************************/ 
		
	//reuse tSumk1, tSumk2, and tDoubleSum;
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nK); 
	
	//at locus 0 -- leftmost; 
	for (k = 0; k < nK; k++)
		m_phi[0][k] = alpha[0][k] * prG(theta[0][k], GetsnpGT(0));
	
	//cal the marginal sum at locus 0 as in appendix A;
	tSum[0] = 0.0;
	for (k = 0; k < nK; k++)
		tSum[0] += m_phi[0][k];
	
	phiScale[0] = 0;
		
	//do the recursive calculations;
	for (m = 0; m < nLoci-1; m++)
	{
		for (k = 0; k < nK; k++)
		{
			double temp = (1.0 - r[m+1]) * m_phi[m][k] + r[m+1] * alpha[m+1][k] * tSum[m]; 
			m_phi[m+1][k] = temp * prG(theta[m+1][k], GetsnpGT(m+1));
		}   //note in above recursion we use probJ[m+1], instead of m as in appendix. 
		
		//sum for locus m+1. 
		tSum[m+1] = 0.0;
		for (k = 0; k < nK; k++)
			tSum[m+1] += m_phi[m+1][k];
		 
		double tscale = 0; 
		if(tSum[m+1] < 1e-100) 
		{
			cout << "sum of phi is 0!" << endl; 
			exit(0); 
		}
		else 
			tscale = -log10(tSum[m+1]); 
		if(tscale > 2)
		{
			phiScale[m+1] = phiScale[m] + tscale;  
			double dummy = pow(10., tscale);
			tSum[m+1] *= dummy; 
			for (k = 0; k < nK; k++)
				m_phi[m+1][k] *= dummy;
			//rescale; 
		}
		else
			phiScale[m+1] = phiScale[m]+tscale;  
	 }                                                          
	/*end calculate Phi for each individual; */

	logLikelihood = log10(tSum[nLoci-1]) - phiScale[nLoci-1];
//	cout << "logLikelihood " << logLikelihood << endl; 
	//tSum[nLoci-1] is the likilihood of g_i given model parameters. 
	
	/***************************************************************************  
	| calc probz given G and Nu and normalize, see appendix A at the end. 
	| input: m_phi, m_beta. 
	| output: m_prZmk = m_phi. 
	***************************************************************************/ 
		
	
	if(ball == 0) 
	{
		for(m = 0; m < nLoci; m++)
		{
			real tPZtotal = 0.0; //normalization factor. 
			for (k = 0; k < nK; k++)
			{		
				(m_beta)[m][k] *= (m_phi)[m][k]; 
				tPZtotal += (m_beta)[m][k];
			} 
					
			for (k = 0; k < nK; k++)
				(m_beta)[m][k] /= (tPZtotal); 	
		}
		m_prZ = m_beta; 
	}

	else if (ball == 1)
	{
		for(m = 0; m < nLoci; m++)
		{
			real tPZtotal = 0.0; //normalization factor. 
			for (k = 0; k < nK; k++)
			{		
				(m_phi)[m][k] *= (m_beta)[m][k]; 
				tPZtotal += (m_phi)[m][k];
			} 
					
			for (k = 0; k < nK; k++)
				(m_phi)[m][k] /= (tPZtotal); 	
		}
		m_prZ = m_phi; 
	}
	
	/* m_phi is  probZ. */

	if (ball == 0) 
	{
		delete[] tSum;
		return; 
	}
	/***************************************************************************  
	| calc expected J for eack {m,k} for each individual;  as in appendix c.
	| input: m_phi, m_beta, and alpha. 
	| output: m_expectJmk, and m_expectJm. 
	| reuse the caculation did in calc Phi;
	| with prG[m]=p(g_i(<=m)|nu) <=> tDoubleSum. tSumk1 is the marginal phi. 
	***************************************************************************/ 
	if (m_expectJmk == NULL)
		m_expectJmk = Allocate2DMatrix(nLoci, nK);

    //locus 0;
	for (k = 0; k < nK; k++)
		m_expectJmk[0][k] = log(alpha[0][k]);  
	
	for (m = 1; m < nLoci; m++)
	{
		double dummy = log(10.0) * (phiScale[m-1] + betaScale[m] - phiScale[nLoci-1]);
		for (k = 0; k < nK; k++)
		{
			m_expectJmk[m][k] = log(tSum[m-1]) + log(r[m]) + log(prG(theta[m][k], GetsnpGT(m))) + log(m_beta[m][k]);  
			//tsum is for m_phi 
			m_expectJmk[m][k] += log(alpha[m][k]);
			m_expectJmk[m][k] -= tSum[nLoci-1];  
			m_expectJmk[m][k] -= dummy;  
		}
	} 

	for (int m = 0; m < nLoci; m++)
	{
		for (int k = 0; k < nK; k++)
		{
			if(m_expectJmk[m][k] < -100) m_expectJmk[m][k] = 0; 
			else m_expectJmk[m][k] = exp(m_expectJmk[m][k]); 
		}
	}
	delete[] tSum;

    //locus 0;
//	for (k = 0; k < nK; k++)
//		m_expectJmk[0][k] = 1.0 * (alpha)[0][k];  
//    /* for the first locus, the recursion as in the end of appendix no longer work. 
//	   here we use the "fact" or assumption that conditional on jump, where to jump is 
//	   independent of where it comes from. 
//	*/
//	
//	for (m = 1; m < nLoci; m++)
//	{
//		dummy = pow(10.0, (double)(phiScale[m-1] + betaScale[m] - phiScale[nLoci-1]));
//		
//		for (k = 0; k < nK; k++)
//		{
//			m_expectJmk[m][k] = tSum[m-1] * r[m] * prG(theta[m][k], GetsnpGT(m)) * m_beta[m][k];  
//			//tSum is for m_phi 
//			m_expectJmk[m][k] *= (alpha)[m][k];
//			m_expectJmk[m][k] /= tSum[nLoci-1];  
//			m_expectJmk[m][k] /= dummy;
//		}
//	} 
//
//	delete[] tSum;
	/* end calc expectJmk. */

	/***************************************************************************  
	| calc thetablock and its inner product with m_prZ, as (c1) in appendix c. 
	| input: snpGT,  m_prZ[m][k1][k2], theta[m][k]. 
	| output: m_top. 
	***************************************************************************/ 
	if (m_top == NULL)
		m_top = Allocate2DMatrix(nLoci, nK);
	if (m_bot == NULL)
		m_bot = Allocate2DMatrix(nLoci, nK);
			
	for (m = 0; m < nLoci; m++)
	{
		for (k = 0; k < nK; k++)
			m_bot[m][k] = m_prZ[m][k];
		
		switch(GetsnpGT(m))
		{
			case 0: 
				for (k = 0; k < nK; k++)
					m_top[m][k] = 0.0;
				break;
			case 1:
			    for (k = 0; k < nK; k++)
					m_top[m][k] = m_prZ[m][k];      
				break;
		  
			case QQ:
			   	for (k = 0; k < nK; k++)
				{
//					m_top[m][k] = m_prZ[m][k] * theta[m][k];        
					m_top[m][k] = 0.0; 
					m_bot[m][k] = 0.0; 
				}
				break;    
			default:
				fprintf(stderr, "wrong genotype encountered in calc theta block \n");
				exit(0);
				break;
		}
	}	/* end calc top. */
}
//{
//	real * r = pMD->Getr(popLabel); 
//	real ** theta = pMD->Gettheta();
//	real ** alpha = pMD->Getalpha(popLabel);
//	
//	int k, m;
//
//	/***************************************************************************  
//	| calculate backward prob. m_beta for each individual. 
//	| input: alpha.
//	| output: m_beta. 
//	***************************************************************************/ 
//
//	real * tSum = new real[nLoci]; //(real *) malloc(sizeof(real) * nLoci);  
//	
//	if(m_beta == NULL)
//		m_beta = Allocate2DMatrix(nLoci, nK); 
//	
//	//at the rightmost locus M; 
//	for (k = 0; k < nK; k++)
//		(m_beta)[nLoci-1][k] = 0.0; //log(1.0); 
//	
//	double * ta = new double[nK]; 
//	tSum[nLoci-1] = 0.0;
//	for (k = 0; k < nK; k++)
//		ta[k] = log(prG(theta[nLoci-1][k], GetsnpGT(nLoci-1))) + log(alpha[nLoci-1][k]); //beta[nLoci-1] = 1.0
//	tSum[nLoci-1] = sumlog(ta, nK); 
//	
//	//do the recursive calc backwards; 
//	for (m = nLoci - 1; m > 0; m--)
//	{
//		for (k = 0; k < nK; k++)
//		{
//			double tb[2]; 
//			tb[0] = log(1.0 - r[m]) + log(prG(theta[m][k], GetsnpGT(m))) + m_beta[m][k]; //1-r[m] = probJ[m][0]
//			tb[1] = log(r[m]) + tSum[m]; 
//			m_beta[m-1][k] = sumlog(tb, 2);
//		}
//		
//		//margin sum and real sum for locus m-1; 
//		for (k = 0; k < nK; k++)
//			ta[k] = log(prG(theta[m-1][k], GetsnpGT(m-1))) + m_beta[m-1][k] + log(alpha[m-1][k]);
//		tSum[m-1] = sumlog(ta, nK); 
//	}
//	/* end calc beta, real checked. */ 
//	
//	/***************************************************************************  
//	| calc forward prob. m_phi 
//	| input: alpha, probJ  
//	| output: m_phi 
//	***************************************************************************/ 
//		
//	//reuse tSumk1, tSumk2, and tDoubleSum;
//	if(m_phi == NULL)
//		m_phi = Allocate2DMatrix(nLoci, nK); 
//	
//	//at locus 0 -- leftmost; 
//	for (k = 0; k < nK; k++)
//		m_phi[0][k] = log(alpha[0][k]) + log(prG(theta[0][k], GetsnpGT(0)));
//	
//	//cal the marginal sum at locus 0 as in appendix A;
//	for (k = 0; k < nK; k++)
//		ta[k] = m_phi[0][k];
//	tSum[0] = sumlog(ta, nK);
//	
//	//do the recursive calculations;
//	for (m = 0; m < nLoci-1; m++)
//	{
//		for (k = 0; k < nK; k++)
//		{
//			double tt = log(prG(theta[m+1][k], GetsnpGT(m+1)));
//			double tb[2]; 
//			tb[0] = tt + log(1.0 - r[m+1]) + m_phi[m][k]; 
//			tb[1] =  tt + log(r[m+1]) + log(alpha[m+1][k]) + tSum[m]; 
//			m_phi[m+1][k] = sumlog(tb, 2); 
//		}   //note in above recursion we use probJ[m+1], instead of m as in appendix. 
//		
//		//sum for locus m+1. 
//		for (k = 0; k < nK; k++)
//			ta[k] = m_phi[m+1][k];
//		tSum[m+1] = sumlog(ta, nK);  
//	 }                                                          
//	/*end calculate Phi for each individual; */
//
//	logLikelihood = tSum[nLoci-1];
////	cout << "logLikelihood " << logLikelihood << endl; 
//	//tSum[nLoci-1] is the likilihood of g_i given model parameters. 
//	
//	/***************************************************************************  
//	| calc probz given G and Nu and normalize, see appendix A at the end. 
//	| input: m_phi, m_beta. 
//	| output: m_prZmk = m_phi. 
//	***************************************************************************/ 
//		
//	
//	if(ball == 0) 
//	{
//		for(m = 0; m < nLoci; m++)
//		{
//			real tPZtotal = 0.0; //normalization factor. 
//			for (k = 0; k < nK; k++)
//			{		
//				(m_beta)[m][k] += (m_phi)[m][k]; 
//				ta[k] = (m_beta)[m][k];
//			}    
//			tPZtotal = sumlog(ta, nK); 
//					
//			for (k = 0; k < nK; k++)
//				(m_beta)[m][k] -= (tPZtotal); 	
//		}
//		m_prZ = m_beta; 
//	}
//
//	else if (ball == 1)
//	{
//		for(m = 0; m < nLoci; m++)
//		{
//			real tPZtotal = 0.0; //normalization factor. 
//			for (k = 0; k < nK; k++)
//			{		
//				(m_phi)[m][k] += (m_beta)[m][k]; 
//				ta[k] = (m_phi)[m][k];
//			} 
//			tPZtotal = sumlog(ta, nK); 
//					
//			for (k = 0; k < nK; k++)
//				(m_phi)[m][k] -= (tPZtotal); 	
//		}
//		m_prZ = m_phi; 
//	}
//	
//	/* m_phi is  probZ. */
//
//	if (ball == 0) 
//	{
//		delete[] tSum;
//		return; 
//	}
//	/***************************************************************************  
//	| calc expected J for eack {m,k} for each individual;  as in appendix c.
//	| input: m_phi, m_beta, and alpha. 
//	| output: m_expectJmk, and m_expectJm. 
//	| reuse the caculation did in calc Phi;
//	| with prG[m]=p(g_i(<=m)|nu) <=> tDoubleSum. tSumk1 is the marginal phi. 
//	***************************************************************************/ 
//	if (m_expectJmk == NULL)
//		m_expectJmk = Allocate2DMatrix(nLoci, nK);
//
//    //locus 0;
//	for (k = 0; k < nK; k++)
//		m_expectJmk[0][k] = log(alpha[0][k]);  
//	
//	for (m = 1; m < nLoci; m++)
//	{
//		for (k = 0; k < nK; k++)
//		{
//			m_expectJmk[m][k] = tSum[m-1] + log(r[m]) + log(prG(theta[m][k], GetsnpGT(m))) + m_beta[m][k];  
//			//tSum is for m_phi 
//			m_expectJmk[m][k] += log(alpha[m][k]);
//			m_expectJmk[m][k] -= tSum[nLoci-1];  
//		}
//	} 
//
//	for (int m = 0; m < nLoci; m++)
//	{
//		for (int k = 0; k < nK; k++)
//		{
//			if(m_expectJmk[m][k] < -100) m_expectJmk[m][k] = 0; 
//			else m_expectJmk[m][k] = exp(m_expectJmk[m][k]); 
//		}
//	}
//	delete[] tSum;
//	/* end calc expectJmk. */
//
//	/***************************************************************************  
//	| calc thetablock and its inner product with m_prZ, as (c1) in appendix c. 
//	| input: snpGT,  m_prZ[m][k1][k2], theta[m][k]. 
//	| output: m_top. 
//	***************************************************************************/ 
//	if (m_top == NULL)
//		m_top = Allocate2DMatrix(nLoci, nK);
//	if (m_bot == NULL)
//		m_bot = Allocate2DMatrix(nLoci, nK);
//			
//	for (m = 0; m < nLoci; m++)
//	{
//		for (k = 0; k < nK; k++)
//			m_bot[m][k] = exp(m_prZ[m][k]);
//		
//		switch(GetsnpGT(m))
//		{
//			case 0: 
//				for (k = 0; k < nK; k++)
//					m_top[m][k] = 0;
//				break;
//			case 1:
//			    for (k = 0; k < nK; k++)
//					m_top[m][k] = exp(m_prZ[m][k]);      
//				break;
//		  
//			case QQ:
//			   	for (k = 0; k < nK; k++)
//				{
////					m_top[m][k] = m_prZ[m][k] * theta[m][k];        
//					m_top[m][k] = 0; 
//					m_bot[m][k] = 0; 
//				}
//				break;    
//			default:
//				fprintf(stderr, "wrong genotype encountered in calc theta block \n");
//				exit(0);
//				break;
//		}
//	}	/* end calc top. */
//	delete[] ta; 
//}

void HapInd::FreeMemAfterEM(void)
{
	m_prZ = NULL; 
	if(phiScale) {delete[] phiScale; phiScale = NULL;}
	if(betaScale) {delete[] betaScale; betaScale = NULL;}
	Free2DMatrix(m_expectJmk); 	m_expectJmk = NULL;
	Free2DMatrix(m_phi);    	m_phi = NULL; 
	Free2DMatrix(m_beta);       m_beta = NULL;
    Free2DMatrix(m_top); 		m_top = NULL;
    Free2DMatrix(m_bot); 		m_bot = NULL;
}

void HapInd::calc_snp_dstr(int nLoci, int nK, real ** theta)
{                                
	if(snp_dstr == NULL) 
	{
		snp_dstr = new real[nLoci * 2]; 
		for (int m = 0; m < (2 * nLoci); m++)
			snp_dstr[m] = 0; 
	}
	for (int m = 0; m < nLoci; m++)
	{
		switch(GetsnpGT(m))
		{
			case 0: 
 				snp_dstr[2*m] += 1.0;
				break;
			case 1:
  				snp_dstr[2*m] += 0.0;
				break;
			case QQ:
				pair<double, double> pp;
				pp.first = pp.second = 0.0; 
				for (int k = 0; k < nK; k++)
				{
					pp.first += (1.0 - theta[m][k]) * (m_prZ[m][k]);	   	
				}
				snp_dstr[2*m] += pp.first;
				break;
		}
	}
} 

void HapInd::norm_snp_dstr(int nLoci, int nEM)
{
	if(snp_dstr == NULL)
	{
		cout << " illegal operation" << endl; 
		return; 
	}
	
	for (int m = 0; m < nLoci; m++)
	{
		snp_dstr[2*m] /= nEM; 
		snp_dstr[2*m+1] = 1.0 - snp_dstr[2*m]; 
	}
}

void HapInd::joint_imputation(ModelParam * pMP, int runs, short * snpImputed, int nLoci, int nK, int ns, int * index)
{
	snpImputed = NULL;  
}

#if defined (IMPUTATION)
void HapInd::MaskSNPs(int nLoci, int nMasked, int * A)
{
	int i;
    if (pMask == NULL)
   		pMask = new class Mask[nMasked];
	for (i = 0; i < nMasked; i++)
	{
		pMask[i].pos = A[i];
		pMask[i].snp = GetsnpGT(pMask[i].pos);
		SetsnpGT(pMask[i].pos, (char)('0'+QQ));
	}
}

void HapInd::ImputeMaskedCounting(int nMasked, int * err_count, int * na_count)
{
	if(maf == NULL) maf = new real[nMasked];
	for (int i = 0; i < nMasked; i++)
	{
	    int m = pMask[i].pos;
		int original = (pMask[i].snp - '0'); 
		if(original == QQ)
		{
			na_count[i]++; 
			continue;
		}
		double p0 = snp_dstr[2 * m]; 
		double p1 = snp_dstr[2*m+1]; 
		double max = p0;
		int gt = 0; 
		if(p1 > max) gt = 1; 
	   
		if(gt != original) err_count[i]++; 
		maf[i] += p1; 
	}
}
#endif 
