#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "indiv.h"
#include "fpmath.h"
#include "diploid.h"
#include <math.h>

using namespace std;


DipInd::DipInd(void)
{
	m_phi = NULL;
	m_beta = NULL;
	m_top = NULL;
	m_bot = NULL; 
    m_prZ = NULL; 

	zpair = NULL;
	phiz = NULL; 
	snp_dstr = NULL; 
	mgt = NULL; 
	m_ploid = 2; 
}

DipInd::~DipInd(void)
{
 	Free2DMatrix(m_phi);        m_phi = NULL;
	Free2DMatrix(m_beta);      	m_beta = NULL;
 	Free2DMatrix(m_top); 	   	m_top = NULL;
 	Free2DMatrix(m_bot); 	   	m_bot = NULL;

	Free2DIntMatrix(zpair);    	zpair = NULL;
	Free2DMatrix(phiz); 	   	phiz = NULL; 
} 

real DipInd::get_snpmgt(int m) 
{
	real res = 0; 
	if(mgt != NULL) 
		res = mgt[m]; 
	else if(snp_dstr != NULL) 
		res = snp_dstr[2*m+1] + 2 * (1 - snp_dstr[2*m] - snp_dstr[2*m+1]); 
	else if(snpGT[m] == QQ) 
	    res = NA; 
	else 
		res = (real) (snpGT[m]-'0'); 
	return res; 
}

real probJ(int m, int s, real * r)
{
	real res = 0; 
	real temp; 
	switch(s) {
    	case 0: 
			 temp = 1.0 - r[m]; 
			 res = temp * temp; 
			 break; 
		 case 1:
			 res = 2.0 * (1.0 - r[m]) * r[m]; 
			 break;
		 case 2:
			 res = r[m] * r[m];
			 break;
		 default:
			 cout << "error in probJ" << endl; 
			 break;
	}
	return res; 
}

real prG(real t1, real t2, short gt)
{
    real res = 0.0; 
	switch(gt)
	{
		case 0:
			res = (1.0 - t1) * (1.0 - t2);
			break;
		case 1:
			res = t1 + t2 - 2.0 * t1 * t2;
			break;
		case 2:
			res = t1 * t2;
			break;
		case QQ: 
			res = 1.0;
			break;
		default:
			cout << "error in prG " << gt << endl;
			break;
	}
	return res; 
}

//if ball = 0, then run upto probZ; otherwise run all.  
void DipInd::CalcAll(const int nLoci, const int nK,  ModelParam * pMP, const int ball) 
{
    real * r = pMP->Getr(popLabel);
	real ** theta = pMP->Gettheta();
	real ** alpha = pMP->Getalpha(popLabel);
	
	/***************************************************************************  
	| calculate backward prob. m_beta for each individual. 
	***************************************************************************/ 
	if (phiScale == NULL)
		phiScale = new double[nLoci]; 
	if (betaScale == NULL)
		betaScale = new double[nLoci]; // (int *) malloc(sizeof(int) * nLoci); 

	real ** tSumk = Allocate2DMatrix(nLoci, nK);
	real * tDoubleSum = (real*) Allocate1D(sizeof(real), nLoci); //(real *) malloc(sizeof(real) * nLoci);  
	
	if(m_beta == NULL)
		m_beta = Allocate2DMatrix(nLoci, nK*(nK+1)/2); 
	
	//at the rightmost locus M; 
	for (int k1 = 0; k1 < nK; k1++)             
		for (int k2 = k1; k2 < nK; k2++)
			(m_beta)[nLoci-1][imap(nK,k1,k2)] = 1.0;
	betaScale[nLoci-1] = 0;
	
	//cal the marginal sum at locus M as in appendix A;
	tDoubleSum[nLoci-1] = 0.0;
	for (int k1 = 0; k1 < nK; k1++)
	{	
		tSumk[nLoci-1][k1] = 0;
		for (int k = 0; k < nK; k++)
			tSumk[nLoci-1][k1] += (prG(theta[nLoci-1][k1], theta[nLoci-1][k], GetsnpGT(nLoci-1)) * m_beta[nLoci-1][imap(nK,k1,k)] * alpha[nLoci-1][k]);
		tDoubleSum[nLoci-1] += tSumk[nLoci-1][k1] * alpha[nLoci-1][k1];
	}                                                                  

	{
		double tscale = -log10(tDoubleSum[nLoci-1]); 
		betaScale[nLoci-1] = tscale;  
		double dummy = pow(10., (double)tscale);
		tDoubleSum[nLoci-1] *= dummy;
		for (int k1 = 0; k1 < nK; k1++)
		{
			tSumk[nLoci-1][k1] *= dummy;
			for (int k2 = k1; k2 < nK; k2++)
				m_beta[nLoci-1][imap(nK,k1,k2)] *= dummy;
		}   //rescale the partial sums; 
	}

	//do the recursive calc backwards; 
	for (int m = nLoci - 1; m > 0; m--)
	{
		for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = k1; k2 < nK; k2++)
			{
				double temp = 0.5 * probJ(m, 1, r) * (tSumk[m][k1] + tSumk[m][k2]);
				temp += probJ(m, 0, r) * prG(theta[m][k1], theta[m][k2], GetsnpGT(m)) * m_beta[m][imap(nK,k1,k2)];
				m_beta[m-1][imap(nK,k1,k2)] = temp + probJ(m, 2, r) * tDoubleSum[m];
			}   //upper
		
		//margin sum and real sum for locus m-1; 
		tDoubleSum[m-1] = 0.0;
		for (int k1 = 0; k1 < nK; k1++)
		{	
			tSumk[m-1][k1] = 0;
			for (int k = 0; k < nK; k++)
				tSumk[m-1][k1] += prG(theta[m-1][k1], theta[m-1][k], GetsnpGT(m-1)) * m_beta[m-1][imap(nK,k1,k)] * alpha[m-1][k];
			tDoubleSum[m-1] += tSumk[m-1][k1] * alpha[m-1][k1];  
		}
			
								
		if (tDoubleSum[m-1]< 1e-100)
		{
			cout << "sum of beta is 0!" << endl;
			string sfn("output/");
			sfn.append("beta.err");
			fstream outfile;
			outfile.open(sfn.c_str(), ios::out);
			if(!outfile.is_open()) 
			{
				cout << "-bimbam: cannot open file to write:" << sfn << endl;
				return;
			}               

			outfile << "#### "; 
			for (int i = 0; i < nLoci; i++)
				outfile << r[i] << " "; 
			outfile << endl; 


			for (int k = 0; k < nK; k++)
			{
				for (int j = 0; j < nK; j++)
					outfile << m_beta[m][imap(nK, k, j)] * pow(10.0, -betaScale[m]) << " "; 
				outfile << endl; 
			}
			outfile << endl; 

			for (int k = 0; k < nK; k++)
			{
				for (int j = 0; j < nK; j++)
					outfile << m_beta[m-1][imap(nK, k, j)] * pow(10.0, -betaScale[m-1]) << " "; 
				outfile << endl; 
			}
			    
//			outfile << m << " " << m << " " << m << endl; 
//
//			for (int m = 0; m < nLoci; m++)
//			{
//				outfile << tDoubleSum[m] << " "; 
//				outfile << phiScale[m] << " "; 
//				outfile << GetsnpGT(m) << " "; 
//				for (int k = 0; k < nK; k++)
//					outfile << alpha[m][k] << " "; 
//				for (int k = 0; k < nK; k++)
//					outfile << theta[m][k] << " "; 
//				outfile << endl; 
//			}

			outfile.close();
//			return; 
  		exit(0);
		}
		else
		{
			double tscale = -log10(tDoubleSum[m-1]); 
			if (1)
			{
				betaScale[m-1] = betaScale[m] + tscale;  
				double dummy = pow(10., (double)tscale);
				tDoubleSum[m-1] *= dummy;
				for (int k1 = 0; k1 < nK; k1++)
				{
					tSumk[m-1][k1] *= dummy;
					for (int k2 = k1; k2 < nK; k2++)
						m_beta[m-1][imap(nK,k1,k2)] *= dummy;
				}   //rescale the partial sums; 
			}
			else 
				betaScale[m-1] = betaScale[m];  

		}
	}
	/* end calc beta, real checked. */ 
	
	/***************************************************************************  
	| calc forward prob. m_phi 
	| input: alpha,   
	| output: m_phi 
	***************************************************************************/ 
		
	//reuse tSumk and tDoubleSum;
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nK*(nK+1)/2); 
	
	//at locus 0 -- leftmost; 
	for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = k1; k2 < nK; k2++)
			m_phi[0][imap(nK,k1,k2)] = alpha[0][k1] * alpha[0][k2] * prG(theta[0][k1], theta[0][k2], GetsnpGT(0));
	//cal the marginal sum at locus 0 as in appendix A;
	tDoubleSum[0] = 0;
	for (int k1 = 0; k1 < nK; k1++)
	{	
		tSumk[0][k1] = 0;
		for (int k = 0; k < nK; k++)
			tSumk[0][k1] += m_phi[0][imap(nK,k,k1)];
		tDoubleSum[0] += tSumk[0][k1];
	}
	phiScale[0] = 0; 
		
	//do the recursive calculations;
	for (int m = 0; m < nLoci-1; m++)
	{
		for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = k1; k2 < nK; k2++)
			{ 
				double temp = alpha[m+1][k1] * tSumk[m][k2] + alpha[m+1][k2] * tSumk[m][k1];
			    temp *= (.5 * probJ(m+1, 1, r));
				temp += (probJ(m+1, 0, r) * m_phi[m][imap(nK,k1,k2)]); 
				temp += (probJ(m+1, 2, r) * alpha[m+1][k1] * alpha[m+1][k2] * tDoubleSum[m]);
				m_phi[m+1][imap(nK,k1,k2)] = temp * prG(theta[m+1][k1], theta[m+1][k2], GetsnpGT(m+1));
			} //note in above recursion we use probJ[m+1], instead of m as in appendix. 
		//marginal sum for locus m+1. 
		tDoubleSum[m+1] = 0.0;
		for (int k1 = 0; k1 < nK; k1++)
		{	
			tSumk[m+1][k1] = 0;
			for (int k = 0; k < nK; k++)
				tSumk[m+1][k1] += m_phi[m+1][imap(nK,k,k1)];
			tDoubleSum[m+1] += tSumk[m+1][k1];
		}
		
		if(tDoubleSum[m+1] < 1e-100) 
		{
			cout << "sum of beta is 0!" << endl;
			string sfn("output/");
			sfn.append("phi.err");
			fstream outfile;
			outfile.open(sfn.c_str(), ios::out);
			if(!outfile.is_open()) 
			{
				cout << "-bimbam: cannot open file to write:" << sfn << endl;
				return;
			}               

			outfile << m << " " << m << " " << m << endl; 

			for (int m = 0; m < nLoci; m++)
			{
				outfile << tDoubleSum[m] << " "; 
				outfile << phiScale[m] << " "; 
				outfile << GetsnpGT(m) << " "; 
				for (int k = 0; k < nK; k++)
					outfile << alpha[m][k] << " "; 
				for (int k = 0; k < nK; k++)
					outfile << theta[m][k] << " "; 
				outfile << endl; 
			}

			outfile.close();
//			return; 
  		exit(0);
		}
		else
		{
			double tscale = (-log10(tDoubleSum[m+1])); 
			if(1)
			{
				phiScale[m+1] = phiScale[m]+tscale;  
				double dummy = pow(10., (double)tscale);
				tDoubleSum[m+1] *= dummy; 
				for (int k1 = 0; k1 < nK; k1++)
				{
					tSumk[m+1][k1] *= dummy;
					for (int k2 = k1; k2 < nK; k2++)
						m_phi[m+1][imap(nK,k1,k2)] *= dummy;
				}   //rescale the partial sums; 
			}
			else
				phiScale[m+1] = phiScale[m];  
		}                         
	 }                                                          
	/*end calculate Phi for each individual; */

	logLikelihood = log10(tDoubleSum[nLoci-1]) - phiScale[nLoci-1];
	//tDoubleSum[[nLoci-1] is the likilihood of g_i given model parameters. 
	
	/***************************************************************************  
	| calc probz given G and Nu and normalize, see appendix A at the end. 
	| input: m_phi, m_beta. 
	| output: m_prZ,  
	***************************************************************************/ 
		
	if (ball == 0) //joint imputation
	{
		for(int m = 0; m < nLoci; m++)
		{
			real tPZtotal = 0.0;     //for normalize; 
			for (int k1 = 0; k1 < nK; k1++)
				for (int k2 = k1; k2 < nK; k2++)
				{
					(m_beta)[m][imap(nK,k1,k2)] *= (m_phi)[m][imap(nK,k1,k2)]; 
					tPZtotal += (k1 == k2 ? 1.0 : 2.0) * (m_beta)[m][imap(nK,k1,k2)];               //m_beta = m_prZ
				} 
					
			for (int k1 = 0; k1 < nK; k1++)
				for (int k2 = k1; k2 < nK; k2++)
					(m_beta)[m][imap(nK,k1,k2)] /= (tPZtotal); 	          //m_beta = m_prZ
			//normalization for each m.
		}
		m_prZ = m_beta; 
	} 
	
	else if (ball == 1)
	{
		for(int m = 0; m < nLoci; m++)
		{
			real tPZtotal = 0.0;     //for normalize; 
			for (int k1 = 0; k1 < nK; k1++)
				for (int k2 = k1; k2 < nK; k2++)
				{
					(m_phi)[m][imap(nK,k1,k2)] *= (m_beta)[m][imap(nK,k1,k2)]; 
					tPZtotal += (k1 == k2 ? 1.0 : 2.0) * (m_phi)[m][imap(nK,k1,k2)];               //m_phi = m_prZ
				} 
					
			for (int k1 = 0; k1 < nK; k1++)
				for (int k2 = k1; k2 < nK; k2++)
					(m_phi)[m][imap(nK,k1,k2)] /= (tPZtotal); 	          //m_phi = m_prZ
			//normalization for each m.
		}
		m_prZ = m_phi; 
	}

	else cout << "illegal options" << endl; 
	// end calc prZ;	

	if (ball == 0) 
	{
//		if(m_phi) {Free3DMatrix(m_phi); m_phi = NULL;} 
//      m_phi was needed for joint imputation; 		
		Free2DMatrix(tSumk);    tSumk = NULL; 
		free(tDoubleSum);       tDoubleSum = NULL; 
		return; 
	}
	/***************************************************************************  
	| calc expected J for eack {m,k} for each individual;  as in appendix c.
	| input: m_phi, m_beta, and alpha. 
	| output: m_expectJmk 
	| reuse the caculation did in calc Phi;
	| with prG[m]=p(g_i(<=m)|nu) <=> tDoubleSum. tSumk is the marginal phi. 
	***************************************************************************/ 
	if (m_expectJmk == NULL)
		m_expectJmk = Allocate2DMatrix(nLoci, nK);

	// for the first locus, the recursion as in the end of appendix no longer work. 
	// here we use the "fact" or assumption that conditional on jump, where to jump is 
	// independent of where it comes from. 
	for (int k = 0; k < nK; k++)
		m_expectJmk[0][k] = 2.0 * (alpha)[0][k];  
	
	for (int m = 1; m < nLoci; m++)
	{
		double dummy = (phiScale[m-1] + betaScale[m] - phiScale[nLoci-1]);
		for (int k = 0; k < nK; k++)
		{
			m_expectJmk[m][k] = 0.0;
			if(alpha[m][k] > 1e-10) 
			{
				for (int k1 = 0; k1 < nK; k1++)
				{
					double temp = tSumk[m-1][k1] * probJ(m, 1, r);  //tSumk, tDoubleSum is for m_phi 
					temp += 2.0 * probJ(m, 2, r) * tDoubleSum[m-1] * (alpha)[m][k1]; // at scale of phiScale[m-1]; 
					m_expectJmk[m][k] += (temp * prG(theta[m][k], theta[m][k1], GetsnpGT(m)) * m_beta[m][imap(nK,k,k1)]);  
					// salce= phiScale[m-1]+betaScale[m]
				}
				if(m_expectJmk[m][k] > 1e-10) 
				{
					double tt = log(m_expectJmk[m][k]) - log(10.0) * dummy -log(tDoubleSum[nLoci-1]) +log(alpha[m][k]); 
				    m_expectJmk[m][k] = exp(tt); 
				}
			}
		}
	} 

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
			
	for (int m = 0; m < nLoci; m++)
	{
		for (int k = 0; k < nK; k++)
		{
			m_bot[m][k] = 0.0; 
			for (int k1 = 0; k1 < nK; k1++)
				m_bot[m][k] += m_prZ[m][imap(nK,k,k1)]; 
		}
		
		switch(GetsnpGT(m))
		{
			case 0: 
				for (int k = 0; k < nK; k++)
					m_top[m][k] = 0.0;
				break;
			case 1:
			    for (int k = 0; k < nK; k++)
				{
					m_top[m][k] = 0.0;
					for (int k1 = 0; k1 < nK; k1++)
					{
						real t1 = (theta)[m][k] * (1. - (theta)[m][k1]);
						real t2 = t1 + (theta)[m][k1] * (1. - (theta)[m][k]);					
						m_top[m][k] += (m_prZ[m][imap(nK,k,k1)] * t1 / t2);    
					}
				}
				break;
			case 2:
				for (int k = 0; k < nK; k++)
				{
					m_top[m][k] = 0.0;  
					for (int k1 = 0; k1 < nK; k1++)
						m_top[m][k] += (m_prZ)[m][imap(nK,k,k1)];              
				}           
				break;

			case QQ:
			   	for (int k = 0; k < nK; k++)
				{
					m_top[m][k] = 0.0; 
					m_bot[m][k] = 0.0; 
				}
				break;    
				
			default:
				cout << "wrong genotype encountered in calc theta block" << endl;
				safe_exit();
				break;
		}
	}	
	Free2DMatrix(m_beta);       m_beta = NULL;
	Free2DMatrix(tSumk);        tSumk = NULL; 
   	free(tDoubleSum);           tDoubleSum = NULL; 
	/* end calc top. */
}   
//{
//    real * r = pMP->Getr(popLabel);
//	real ** theta = pMP->Gettheta();
//	real ** alpha = pMP->Getalpha(popLabel);
//	
//	/***************************************************************************  
//	| calculate backward prob. m_beta for each individual. 
//	***************************************************************************/ 
//
//	real ** tSumk = Allocate2DMatrix(nLoci, nK);
//	real * tDoubleSum = (real*) Allocate1D(sizeof(real), nLoci); //(real *) malloc(sizeof(real) * nLoci);  
//
//	int nK2 = nK * nK; 
//	double * ta = new double[nK]; 
//	double * ta2 = new double[nK2]; 
//
//	if(m_beta == NULL)
//		m_beta = Allocate2DMatrix(nLoci, nK*(nK+1)/2); 
//	//m_beta becomes log(m_beta); 
//
//	//at the rightmost locus M; 
//	for (int k = 0; k < nK*(nK+1)/2; k++)
//		m_beta[nLoci-1][k] = 0.0; //log(1.0); 
//
//	//cal the marginal sum at locus M as in appendix A;
//	for (int k1 = 0; k1 < nK; k1++)
//	{	
//		for (int k = 0; k < nK; k++)
//			ta[k] = log(prG(theta[nLoci-1][k1], theta[nLoci-1][k], GetsnpGT(nLoci-1))) + m_beta[nLoci-1][imap(nK,k1,k)] + log(alpha[nLoci-1][k]);
// 		tSumk[nLoci-1][k1] = sumlog(ta, nK); 
//	}  
//	for (int k = 0; k < nK; k++)
//		ta[k] = tSumk[nLoci-1][k] + log(alpha[nLoci-1][k]); 
//	tDoubleSum[nLoci-1] = sumlog(ta, nK);
//	
//	//do the recursive calc backwards; 
//	for (int m = nLoci - 1; m > 0; m--)
//	{
//		for (int k1 = 0; k1 < nK; k1++)
//			for (int k2 = k1; k2 < nK; k2++)
//			{
//                double tb[4]; 
//				double tt = log(0.5) + log(probJ(m, 1, r));
//				tb[0] = tt + tSumk[m][k1]; 
//				tb[1] = tt + tSumk[m][k2]; 
//				tb[2] = log(probJ(m, 0, r)) + log(prG(theta[m][k1], theta[m][k2], GetsnpGT(m))) + m_beta[m][imap(nK,k1,k2)];
//				tb[3] = log(probJ(m, 2, r)) + tDoubleSum[m];
//				m_beta[m-1][imap(nK,k1,k2)] = sumlog(tb, 4);  
//			}   //upper
//		
//		//margin sum and real sum for locus m-1; 
//		for (int k1 = 0; k1 < nK; k1++)
//		{	
//			for (int k = 0; k < nK; k++)
//				ta[k] = log(prG(theta[m-1][k1], theta[m-1][k], GetsnpGT(m-1))) + m_beta[m-1][imap(nK,k1,k)] + log(alpha[m-1][k]);
//			tSumk[m-1][k1] = sumlog(ta, nK); 
//		}  
//		for (int k = 0; k < nK; k++)
//			ta[k] = tSumk[m-1][k] + log(alpha[m-1][k]); 
//		tDoubleSum[m-1] = sumlog(ta, nK);
//	}
//
//	/* end calc beta, real checked. */ 
//	
//	/***************************************************************************  
//	| calc forward prob. m_phi 
//	| input: alpha,   
//	| output: m_phi 
//	***************************************************************************/ 
//		
//	//reuse tSumk and tDoubleSum;
//	if(m_phi == NULL)
//		m_phi = Allocate2DMatrix(nLoci, nK*(nK+1)/2); 
//
//	//m_phi to contain log(m_phi); 
//	//at locus 0 -- leftmost; 
//	for (int k1 = 0; k1 < nK; k1++)
//		for (int k2 = k1; k2 < nK; k2++)
//			m_phi[0][imap(nK,k1,k2)] = log(alpha[0][k1]) + log(alpha[0][k2]) + log(prG(theta[0][k1], theta[0][k2], GetsnpGT(0)));
//	//cal the marginal sum at locus 0 as in appendix A;
//	for (int k1 = 0; k1 < nK; k1++)
//	{	
//		for (int k = 0; k < nK; k++)
//			ta[k] = m_phi[0][imap(nK,k,k1)];
//		tSumk[0][k1] = sumlog(ta, nK); 
//	}                                  
//	for (int k = 0; k < nK; k++)
//		ta[k] = tSumk[0][k]; 
//	tDoubleSum[0] = sumlog(ta, nK);
//
//	//do the recursive calculations;
//	for (int m = 0; m < nLoci-1; m++)
//	{
//		for (int k1 = 0; k1 < nK; k1++)
//			for (int k2 = k1; k2 < nK; k2++)
//			{ 
//				double tb[4]; 
//				double la1 = log(alpha[m+1][k1]); 
//				double la2 = log(alpha[m+1][k2]); 
//				double tt = log(0.5) + log(probJ(m+1, 1, r)); 
//				tb[0] = tt + tSumk[m][k1] + la1; 
//				tb[1] = tt + tSumk[m][k2] + la2; 
//				tb[2] = log(probJ(m+1, 0, r)) + m_phi[m][imap(nK,k1,k2)]; 
//				tb[3] = log(probJ(m+1, 2, r)) + la1 + la2 + tDoubleSum[m];
//				m_phi[m+1][imap(nK,k1,k2)] = sumlog(tb, 4) + log(prG(theta[m+1][k1], theta[m+1][k2], GetsnpGT(m+1)));  
//			} //note in above recursion we use probJ[m+1], instead of m as in appendix. 
//		//marginal sum for locus m+1. 
//		for (int k1 = 0; k1 < nK; k1++)
//		{	
//			for (int k = 0; k < nK; k++)
//				ta[k] = m_phi[m+1][imap(nK,k,k1)];
//			tSumk[m+1][k1] = sumlog(ta, nK); 
//		}                                  
//		for (int k = 0; k < nK; k++)
//			ta[k] = tSumk[m+1][k]; 
//		tDoubleSum[m+1] = sumlog(ta, nK);
//		
//	}                                                          
//	/*end calculate Phi for each individual; */
//
//	logLikelihood = tDoubleSum[nLoci-1];
//	//tDoubleSum[[nLoci-1] is the likilihood of g_i given model parameters. 
//	
//	/***************************************************************************  
//	| calc probz given G and Nu and normalize, see appendix A at the end. 
//	| input: m_phi, m_beta. 
//	| output: m_prZ,  
//	***************************************************************************/ 
//		
//	if (ball == 0) //joint imputation
//	{
//		for(int m = 0; m < nLoci; m++)
//		{
//			for (int k1 = 0; k1 < nK; k1++)
//				for (int k2 = 0; k2 < nK; k2++)
//				{
//					(m_beta)[m][imap(nK,k1,k2)] += (m_phi)[m][imap(nK,k1,k2)]; 
//					ta2[nK * k1 + k2] = (m_beta)[m][imap(nK,k1,k2)];  
//				} 
//			double tPZtotal = sumlog(ta2, nK2);               
//					
//			for (int k1 = 0; k1 < nK; k1++)
//				for (int k2 = k1; k2 < nK; k2++)
//					(m_beta)[m][imap(nK,k1,k2)] -= (tPZtotal); 	          //m_beta = m_prZ
//			//normalization for each m.
//		}
//		m_prZ = m_beta; 
//	} 
//	
//	else if (ball == 1)
//	{
//		for(int m = 0; m < nLoci; m++)
//		{
//			for (int k1 = 0; k1 < nK; k1++)
//				for (int k2 = 0; k2 < nK; k2++)
//				{
//					(m_phi)[m][imap(nK,k1,k2)] += (m_beta)[m][imap(nK,k1,k2)]; 
//					ta2[nK * k1 + k2] = (m_phi)[m][imap(nK,k1,k2)]; 
//				} 
//			double tPZtotal = sumlog(ta2, nK2);               
//					
//			for (int k1 = 0; k1 < nK; k1++)
//				for (int k2 = k1; k2 < nK; k2++)
//					(m_phi)[m][imap(nK,k1,k2)] -= (tPZtotal); 	          //m_phi = m_prZ
//			//normalization for each m.
//		}
//		m_prZ = m_phi; 
//	}
//
//	else cout << "illegal options" << endl; 
//	// end calc prZ;	
//
//	if (ball == 0) 
//	{
////		if(m_phi) {Free3DMatrix(m_phi); m_phi = NULL;} 
////      m_phi was needed for joint imputation; 		
//		Free2DMatrix(tSumk);    tSumk = NULL; 
//		free(tDoubleSum);       tDoubleSum = NULL; 
//		return; 
//	}
//	/***************************************************************************  
//	| calc expected J for eack {m,k} for each individual;  as in appendix c.
//	| input: m_phi, m_beta, and alpha. 
//	| output: m_expectJmk 
//	| reuse the caculation did in calc Phi;
//	| with prG[m]=p(g_i(<=m)|nu) <=> tDoubleSum. tSumk is the marginal phi. 
//	***************************************************************************/ 
//	if (m_expectJmk == NULL)
//		m_expectJmk = Allocate2DMatrix(nLoci, nK);
//	//m_expectJmk contains log(m_expectJmk); 
//
//	// for the first locus, the recursion as in the end of appendix no longer work. 
//	// here we use the "fact" or assumption that conditional on jump, where to jump is 
//	// independent of where it comes from. 
//	for (int k = 0; k < nK; k++)
//		m_expectJmk[0][k] = log(2.0) + log(alpha[0][k]);  
//	                                 
//	for (int m = 1; m < nLoci; m++)
//	{
//		for (int k = 0; k < nK; k++)
//		{
//			for (int k1 = 0; k1 < nK; k1++)
//			{
//				double tb[2]; 
//				tb[0] = tSumk[m-1][k1] + log(probJ(m, 1, r));  //tSumk, tDoubleSum is for m_phi 
//				tb[1] = log(2.0) + log(probJ(m, 2, r)) + tDoubleSum[m-1] + log(alpha[m][k1]); // at scale of phiScale[m-1]; 
//				ta[k1] = sumlog(tb, 2) + log(prG(theta[m][k], theta[m][k1], GetsnpGT(m))) + m_beta[m][imap(nK,k,k1)];  
//			}
//			m_expectJmk[m][k] = sumlog(ta, nK); 
//
//			m_expectJmk[m][k] += log(alpha[m][k]);
//			m_expectJmk[m][k] -= tDoubleSum[nLoci-1]; 
////			if(isnan(m_expectJmk[m][k])) 
////			cout << m << "\t" << k  << "\t"<< alpha[m][k] << "\t" << tDoubleSum[nLoci-1] << endl; 
//		}
//
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
//
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
//	for (int m = 0; m < nLoci; m++)
//	{
//		for (int k = 0; k < nK; k++)
//		{
//			for (int k1 = 0; k1 < nK; k1++)
//				ta[k] = m_prZ[m][imap(nK,k,k1)]; 
//			m_bot[m][k] = sumlog(ta, nK); 
//		}        
//		
//		switch(GetsnpGT(m))
//		{
//			case 0: 
//				for (int k = 0; k < nK; k++)
//					m_top[m][k] = -1e100;
//				break;
//			case 1:
//			    for (int k = 0; k < nK; k++)
//				{
//					for (int k1 = 0; k1 < nK; k1++)
//					{
//						real t1 = (theta)[m][k] * (1. - (theta)[m][k1]);
//						real t2 = t1 + (theta)[m][k1] * (1. - (theta)[m][k]);					
//						ta[k1] = m_prZ[m][imap(nK,k,k1)] + log(t1) - log(t2);    
//					}    
//					m_top[m][k] = sumlog(ta, nK); 
//				}
//				break;
//			case 2:
//				for (int k = 0; k < nK; k++)
//				{
//					for (int k1 = 0; k1 < nK; k1++)
//						ta[k1] = (m_prZ)[m][imap(nK,k,k1)];              
// 					m_top[m][k] = sumlog(ta, nK); 
//				}           
//				break;
//
//			case QQ:
//			   	for (int k = 0; k < nK; k++)
//				{
//					m_top[m][k] = -1e100; 
//					m_bot[m][k] = -1e100; 
//				}
//				break;    
//				
//			default:
//				cout << "wrong genotype encountered in calc theta block" << endl;
//				safe_exit();
//				break;
//		}                              
//	}	
//
//	for (int m = 0; m < nLoci; m++)
//	{
//		for (int k = 0; k < nK; k++)
//		{
//			if(m_top[m][k] <= -100) m_top[m][k] = 0; 
//			else m_top[m][k] = exp(m_top[m][k]); 
//			if(m_bot[m][k] <= -100) m_bot[m][k] = 0; 
//			else m_bot[m][k] = exp(m_bot[m][k]); 
//		}
//	}
//
//	delete[] ta; 
//	delete[] ta2; 
//	Free2DMatrix(m_beta);       m_beta = NULL;
//	Free2DMatrix(tSumk);        tSumk = NULL; 
//   	free(tDoubleSum);           tDoubleSum = NULL; 
//	/* end calc top. */
//}   

void DipInd::FreeMemAfterEM(void)
{
	m_prZ = NULL; 
	if(phiScale) {delete[] phiScale; phiScale = NULL;}
	if(betaScale) {delete[] betaScale; betaScale = NULL;}
	Free2DMatrix(m_expectJmk); 	m_expectJmk = NULL;
	if(m_phi) {	Free2DMatrix(m_phi);       	m_phi = NULL;  }
	if(m_beta) { Free2DMatrix(m_beta);      	m_beta = NULL; }
 	Free2DMatrix(m_top); 	   	m_top = NULL;
 	Free2DMatrix(m_bot); 	   	m_bot = NULL;
	Free2DIntMatrix(zpair);     zpair = NULL;
	Free2DMatrix(phiz);         phiz = NULL; 
}

                             
#if defined (IMPUTATION)
void DipInd::MaskSNPs(int nLoci, int nMasked, int * A)
{
    if (pMask == NULL)
   		pMask = new class Mask[nMasked];
	for (int i = 0; i < nMasked; i++)
	{
		pMask[i].pos = A[i];
		pMask[i].snp = (char)('0'+GetsnpGT(pMask[i].pos));
		SetsnpGT(pMask[i].pos, (char)('0'+QQ));
	}
}


void DipInd::ImputeMaskedCounting(int nMasked, int * err_count, int * na_count)
{
	if(maf == NULL) 
	{	
		maf = new real[nMasked];
		for (int i = 0; i < nMasked; i++)
			maf[i] = 0; 
	}
	for (int i = 0; i < nMasked; i++)
	{
		int m = pMask[i].pos;
		double p0 = snp_dstr[2 * m]; 
		double p1 = snp_dstr[2*m+1]; 
		double p2 = 1.0 - p0 - p1; 
		maf[i] += p1 + 2.0 * p2; 
		int original = (pMask[i].snp - '0'); 
		if(original == QQ)
			na_count[i]++; 
		else 
		{
			double max = p0;
			int gt = 0; 
			if(p1 > max) { gt = 1; max = p1; }
			if(p2 > max) { gt = 2;}
			if(gt != original) err_count[i]++; 
		}
	}
}
#endif 

void DipInd::calc_snp_dstr(int nLoci, int nK, real ** theta)
{
	if(snp_dstr == NULL) 
	{
		snp_dstr = new real[nLoci * 2]; 
		for (int m = 0; m < (2 * nLoci); m++)
			snp_dstr[m] = 0.0; 
	}
	
	for (int m = 0; m < nLoci; m++)
	{
		switch(GetsnpGT(m))
		{
			case 0: 
 				snp_dstr[2*m] += 1.0;
    	    	snp_dstr[2*m+1] += 0;
				break;
			case 1:
  				snp_dstr[2*m] += 0.0;
    	    	snp_dstr[2*m+1] += 1.0; 
				break;
			case 2:
			 	snp_dstr[2*m] += 0.0;
        		snp_dstr[2*m+1] += 0.0;
				break;
			case QQ:
				pair<double, double>  pp;
				pp.first = pp.second = 0.0; 
				for (int k1 = 0; k1 < nK; k1++)
					for (int k2 = 0; k2 < nK; k2++)
					{
						real t1 = theta[m][k1];
						real t2 = theta[m][k2];
						pp.first += (1.0 - t1) * (1.0 - t2)  * (m_prZ[m][imap(nK,k1,k2)]);		  	
						pp.second += (t1 + t2 - 2.0 * t1 * t2) * (m_prZ[m][imap(nK,k1,k2)]);       
					}
				snp_dstr[2*m] += pp.first;
				snp_dstr[2*m+1] += pp.second; 
				break;
		}
	}
}   //calc density for this EM run and the sum upto this run. 

void DipInd::norm_snp_dstr(int nLoci, int nEM)
{
	if(snp_dstr == NULL) 
		return;
	for (int m = 0; m < (2*nLoci); m++)
		snp_dstr[m] /= (double) nEM; 
}

void DipInd::joint_imputation(ModelParam * pMP, int runs, short * snpImputed, int nLoci, int nK, int ns, int * index)
{
	if (snpImputed != NULL && index != NULL)
	{
		for (int m = 0; m < ns; m++)
		{
			int pos = index[m];
			snpImputed[m] = GetsnpGT(pos);
		}
		if (runs == -1 || pMP == NULL)           
			return; 
	}
	else if (snpImputed != NULL) {
		for (int m = 0; m < nLoci; m++)
			snpImputed[m] = GetsnpGT(m);
		if (runs == -1 || pMP == NULL)           
			return; 
	}
	
	real * r = pMP[runs].Getr(popLabel); 
	real ** alpha = pMP[runs].Getalpha(popLabel);
	real ** theta = pMP[runs].Gettheta(popLabel);
	if(zpair == NULL)
		zpair = Allocate2DIntMatrix(nLoci, 2); 
	if(phiz == NULL)
		phiz = Allocate2DMatrix(nK, nK);
	real prob = gsl_rng_uniform(gsl_r);
	real cur_prob = 0;
	int bingle = 0; 
	for ( int k1 = 0; k1 < nK; k1++)
	{
		for ( int k2 = 0; k2 < nK; k2++)
		{
			cur_prob += (m_prZ[nLoci-1][imap(nK,k1,k2)]);                 
			if(prob < cur_prob) {
				bingle = 1;
				zpair[nLoci-1][1] = k2;
				break;
			}
		}
		if(bingle)
		{
			zpair[nLoci-1][0] = k1;
			break;
		}
	}
	if(bingle == 0)
		zpair[nLoci-1][0] = zpair[nLoci-1][1] = nK-1;
    //sampling Z_M; 
	
	for (int m = nLoci - 2; m >= 0; m--)
	{
		int kp1 = zpair[m+1][0];
		int kp2 = zpair[m+1][1];             
		double total = 0.0; 
		for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = 0; k2 < nK; k2++)
			{
				double pk1_kp1, pk2_kp2;
				if (k1 == kp1)
					pk1_kp1 = 1 - r[m+1] + r[m+1] * alpha[m+1][kp1];
				else
					pk1_kp1 = r[m+1] * alpha[m+1][kp1];
				
                if (k2 == kp2)
					pk2_kp2 = 1 - r[m+1] + r[m+1] * alpha[m+1][kp2];
				else
					pk2_kp2 = r[m+1] * alpha[m+1][kp2];
                
				total += m_phi[m][imap(nK,k1,k2)] * pk1_kp1 * pk2_kp2;
			}    
		for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = 0; k2 < nK; k2++)
			{               
				if(total < 1e-100) phiz[k1][k2] = 0; 
				else 
					phiz[k1][k2] = phiz[k1][k2] / total; 
			}

        real prob = gsl_rng_uniform(gsl_r);
		real cur_prob = 0;
		int bingle = 0; 
		for ( int k1 = 0; k1 < nK; k1++)
		{
			for ( int k2 = 0; k2 < nK; k2++)
			{
				cur_prob += phiz[k1][k2];
				if(prob < cur_prob) {
					bingle = 1;
					zpair[m][1] = k2;
					break;
				}
			}
			if(bingle)
			{
				zpair[m][0] = k1;
				break;
			}
		}       	
		if (bingle == 0)
			zpair[m][0] = zpair[m][1] = nK - 1 ;
	} // finish sample z vector.

	if(snpImputed != NULL && index != NULL) 
	{	
		for (int i = 0; i < ns; i++)
		{
			if(snpImputed[i] != QQ) continue; 
			int m = index[i];
			int k1 = zpair[m][0];
			int k2 = zpair[m][1];
			real theta1 = theta[m][k1];
			real theta2 = theta[m][k2];
			short int geno = 0; 
			if (gsl_rng_uniform(gsl_r) < theta1)
				geno++; 

			if (gsl_rng_uniform(gsl_r) < theta2)
				geno++;
			
			snpImputed[i] = geno;
		}
	}
	else if(snpImputed != NULL) {
		for (int m = 0; m < nLoci; m++)
		{
			if(GetsnpGT(m) != QQ) continue; 
			int k1 = zpair[m][0];
			int k2 = zpair[m][1];
			real theta1 = theta[m][k1];
			real theta2 = theta[m][k2];
			short int geno = 0; 
			if (gsl_rng_uniform(gsl_r) < theta1)
				geno++; 

			if (gsl_rng_uniform(gsl_r) < theta2)
				geno++;
			
			snpImputed[m] = geno;
		} 
	}
}

