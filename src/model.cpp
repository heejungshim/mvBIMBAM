#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include "model.h"
#include "fpmath.h"
#include "control.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_permute.h"
#include "gsl/gsl_cdf.h"
using namespace std;

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif


void ModelnData::print_progress_bar(int last, char* str, int p, int total)
{
	if(m_silence == 1) return; 
	int progress = (int) (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[100];
	memset(bar, '\0', 100); 
	if(!last) {
		for (int i = 0; i < barsize; i++)
			bar[i] = '>'; 
		printf("%s [%-50s] %d%%\r", str, bar, progress); 
		fflush(stdout); 	
	} else {
		for (int i = 0; i < barsize; i++)
			bar[i] = '>'; 
		printf("%s [%-50s] 100%%\n", str, bar); 
	}
}

CaseCtrlData::CaseCtrlData(int row, int col, real * vph) 
{
	ni = row;        //number of indiv;
	np = col;       // dim of beta; 
	mx=gsl_matrix_alloc(ni, np); 
	vy=gsl_vector_alloc(ni); 
	delta = gsl_vector_alloc(np); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(vy, i, vph[i]); 
}

CaseCtrlData::~CaseCtrlData() 
{
	gsl_matrix_free(mx);
	gsl_vector_free(vy); 
	gsl_vector_free(delta); 
}

PosnBF::PosnBF()
{
	bf = 0.0; 
	var = 0.0; 
}

PosnBF::~PosnBF()
{
	pos.clear(); 
}


ModelnData::ModelnData(void)
{
	snpmgt = NULL; 
	snpIn01 = NULL; 
	m_not_snp = 0; 
	m_silence = 0; 
	m_exclude_maf = -1; 
	m_exclude_miss = 1.00; 
	m_exclude_nopos = 1; 
	m_allele_coding_mode = 0; // 0 minor is reference; 1 major is referece; genotype is count of reference alleles. 
	m_sortQ = 0; 
	m_num = 0; 
	m_df = 2; 
	nK = 10; 
	randSeed = 0;
	nEMRuns = 10;
	nMaxSteps = 1;
	nWarmSteps = 20; 
#if defined (IMPUTATION)
	percentMasked = .10;
	vPARCrs.clear(); 
	vPARCpos.clear(); 
	parcLoci = 0; 
#endif
	
	nLoci = 0; 
	nPH = 1; 
	nLevel = 1; 
	nMultiSnp = 0; 
	nImpute = 0; 
	
	fnOutput.assign("\0");
	fnRsPos.assign("\0");
	hasPopStructure.assign("0"); 
	
	nGeneFlank = 0; 
	nDip = 0;
	nHap = 0;
	nSubPop = 1;
	pMP = NULL;
	vnSubHap.clear(); 
	m_pv.clear(); 
	cc = 0; 
#if defined (MPI_ENABLED)
	pvLoadDip = NULL; 
#endif 
}

void ModelnData::process_prior(void)
{

  if(vsigma_a.size() == 0 || vsigma_d.size() == 0 || vsigma_a.size() != vsigma_d.size())      
    {
      fplog << "## BIMBAM: Use default priors" << endl; 
      vsigma_a.clear(); 
      vsigma_d.clear();
      vsigma_a.push_back(0.05);  vsigma_d.push_back(0.0125); 
      vsigma_a.push_back(0.1);   vsigma_d.push_back(0.025); 
      vsigma_a.push_back(0.2);   vsigma_d.push_back(0.05); 
      vsigma_a.push_back(0.4);   vsigma_d.push_back(0.1); 
      //		vsigma_a.push_back(0.8);   vsigma_d.push_back(0.2); 
    }                             
  else 
    fplog << "## BIMBAM: Use user specified priors" << endl; 
}

#if defined (MPI_ENABLED)
void ModelnData::Loading(void)
{
    pvLoadDip = new vector<int>[nProc];
    int i = 0;
    while (i < nIndiv) {
        int which_proc = i % nProc;
        pvLoadDip[which_proc].push_back(i++);
    }

	fplog << "## indiv load on process 0: "; 
    if(procID == 0)
    {
        for (unsigned j = 0; j < pvLoadDip[0].size(); j++)
            fplog << pvLoadDip[0].at(j) << " ";
        fplog << endl;
    }
	
	nMaxLoad = (int) (1.0 * nIndiv / nProc) + 1; 
	nSlice = (int) (1.0 * nLoci / nProc) + 1; 
	if (nSlice > nLoci) nSlice = nLoci; 
	
    int div = nLoci / nProc;
    int rem = nLoci % nProc;
    loadLoci = new int[nProc+1];
    for (int p = 1; p <= nProc; p++)
        loadLoci[p] = div;
    for (int p = 1; p <= rem; p++)
        loadLoci[p+1]++;
    loadLoci[0] = 0;
    for (int p = 1; p <= nProc; p++)
        loadLoci[p] += loadLoci[p-1];
    //this is to slice the column.  
    if (procID == 0) {
        fplog << "## loadLoci: ";
        for (int i = 0; i <= nProc; i++)
            fplog << loadLoci[i] << " ";
        fplog << endl;
    }
}
#endif

ModelnData::~ModelnData()
{
	if(pMP) {delete[] pMP;  pMP = NULL;} 
	vGin.clear();
	vPin.clear(); 
	vFileIndiv.clear(); 
	vFilePloid.clear(); 
	vv_phval.clear();
	vv_phdex.clear();
	vsRsnum.clear(); 
	vPerm.clear();
	vnSubHap.clear(); 
	nLoci = 0; 
	
	nDip = 0;
	nHap = 0;
	nSubPop = 1;
	m_pv.clear(); 
#if defined (MPI_ENABLED)
   if(pvLoadDip) {delete[] pvLoadDip; pvLoadDip = NULL;}
   vvLoad.clear(); 
#endif 
	//cout << "model destructor being called" << endl; 
}

void ModelnData::clean_geno_pheno(void)
{
	if(pMP) {delete[] pMP;  pMP = NULL;} 
	vFileIndiv.clear(); 
	vFilePloid.clear(); 
	vv_phval.clear();
	vv_phdex.clear();
	vsRsnum.clear(); 
	vPerm.clear();
	vnSubHap.clear(); 
	
	for (int i = 0; i < nIndiv; i++)
		if(pIndiv && pIndiv[i]) delete pIndiv[i];
	if(pIndiv) delete[] pIndiv; 

	nDip = 0;
	nHap = 0;
	nSubPop = 1;
	m_pv.clear(); 
#if defined (MPI_ENABLED)
   if(pvLoadDip) {delete[] pvLoadDip; pvLoadDip = NULL;}
   vvLoad.clear(); 
#endif 
}

void ModelnData::InitModelParam(void)
{
	if (pMP == NULL)
		pMP = new ModelParam[nEMRuns];
	
    for (int i = 0; i < nEMRuns; i++)
		pMP[i].Init(nLoci, nK, nSubPop);
  
}

double ModelnData::calc_bf(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, real ** snpInpr, class PosnBF * pt) 
{  
	m_df = 2; 
	int ns = pt->pos.size();
	int ni = vin->size(); 
	int col =  1 + m_df * ns; 
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, vph->at(vin->at(i))); 
	gsl_matrix * gX = gsl_matrix_alloc(ni, col); 
	for (int ind = 0; ind < ni; ind++)
	{
		int i = vin->at(ind); 
		gsl_matrix_set(gX, ind, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			real p0 = snpInpr[i][2 * s];
			real p1 = snpInpr[i][2*s+1];
			real p2 = 1.0 - p0 - p1; 
			if(m_not_snp == 0 && (p2 < -1e-3)) 
			{
				cout << "-bimbam: illegal probablities in calc_bf " << p0 << "\t" << p1 << "\t" << p2 << endl; 
			}
			gsl_matrix_set(gX, ind, 2*j+1, p1 + 2 * p2);
			if(m_df == 2) 
				gsl_matrix_set(gX, ind, 2*j+2, p1);
		}
	}

	real logbf = bf_core(inv_va, inv_vd, ni, ns, gX, gph);  
	gsl_matrix_free(gX);
	gsl_vector_free(gph); 
	return (logbf);    	
}
                  
double ModelnData::calc_bf(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, short ** genotype, class PosnBF * pt) 
{  
	int ns = pt->pos.size();
	vector<int> ind_has_all;
	for (unsigned ind = 0; ind < vin->size(); ind++)
	{
		int i = vin->at(ind); 
		int bingle = 0; 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			if(fabs(genotype[i][s] - NA) < 1e-6)
			{
				bingle = 1; 
				break;
			}
		}
		if(!bingle)  ind_has_all.push_back(i); 
	}
	if(ind_has_all.size() == 0) return 1.0; 
	
	int ni = ind_has_all.size(); 
	int col = 1 + m_df * ns; 
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, vph->at(ind_has_all.at(i))); 
	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int ind = 0; ind < ni; ind++)
	{
		int i = ind_has_all.at(ind); 
		gsl_matrix_set(gXX, ind, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			gsl_matrix_set(gXX, ind, m_df*j+1, genotype[i][s]);
			if(m_df == 2)
				gsl_matrix_set(gXX, ind, m_df*j+2, (genotype[i][s] == 1 ? 1 : 0));
		}
	}
	
	real logbf = bf_core(inv_va, inv_vd, ni, ns, gXX, gph);  
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	vector<int>().swap(ind_has_all); 
	return (logbf);    	
}

double ModelnData::calc_bf_mgt(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, real ** mgt, class PosnBF * pt) 
{  
	m_df = 1; 
	int ns = pt->pos.size();
	int ni = vin->size(); 
	int col =  1 + m_df * ns; 
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, vph->at(vin->at(i))); 
	gsl_matrix * gX = gsl_matrix_alloc(ni, col); 
	for (int ind = 0; ind < ni; ind++)
	{
		int i = vin->at(ind); 
		gsl_matrix_set(gX, ind, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			gsl_matrix_set(gX, ind, m_df*j+1, mgt[i][s]);
		}
	}

	real logbf = bf_core(inv_va, inv_vd, ni, ns, gX, gph);  
	gsl_matrix_free(gX);
	gsl_vector_free(gph); 
	return (logbf);    	
}

double ModelnData::bf_core(real inv_va, real inv_vd, int ni, int ns, gsl_matrix * gXX, gsl_vector * gph)
{
	int col = 1 + m_df * ns; 

	double mean = 0; 
	double var = 0; 
	for (int i = 0; i < ni; i++)
		mean += gsl_vector_get(gph, i); 
	mean /= ni; 
	for (int i = 0; i < ni; i++)
	{
		double temp = gsl_vector_get(gph, i) - mean; 
		var += temp * temp; 
	}
	if (fabs(var) < 1e-10) 
	{
		m_beta.clear();
		for (int i = 0; i < col; i++)
			m_beta.push_back(0); 
		return 0; 
	}
	double yy = 0; 
	gsl_blas_ddot(gph, gph, &yy); 
	
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_matrix * gXtX = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gXX, gXX, 0.0, gXtX); 
//	cout << " gXX %*% gXX" << endl; 
//	for (int i = 0; i < col; i++)
//	{
//		for (int j = 0; j < col; j++)
//			cout << gsl_matrix_get(ginvOmega, i, j) << "\t"; 
//    	cout << endl; 
//	}

	gsl_vector * gyx = gsl_vector_alloc(col); 
	gsl_vector * res = gsl_vector_alloc(col); 
	m_beta.clear(); 
	for (int i = 0; i < col; i++)
		m_beta.push_back(0); 
	
    vector<double> vsa;   
	vector<double> vsd; 
	
	if(inv_va < 1e-10  && inv_vd < 1e-10)
	{
		for (unsigned r = 0; r < vsigma_a.size(); r++)
		{
			vsa.push_back(1.0 / vsigma_a.at(r) / vsigma_a.at(r)); 
			vsd.push_back(1.0 / vsigma_d.at(r) / vsigma_d.at(r)); 
		}
	}
	else 
	{
		vsa.push_back(inv_va);  // (1.0 / 0.2 /0.2); 
		vsd.push_back(inv_vd); // (1.0 / 0.05 / 0.05);    
	}
	int repeat = (int) vsa.size(); 
	gsl_vector * bf = gsl_vector_alloc(repeat); 
	for (unsigned p = 0; p < vsa.size(); p++)
	{
		gsl_matrix_memcpy(ginvOmega, gXtX); 
		double inv_va = vsa.at(p); 
		double inv_vd = vsd.at(p); 
		for (int i = 0; i < ns; i++)
		{
			real tmp = 0; 
			if(m_df == 1)
			{
				tmp = gsl_matrix_get(ginvOmega, i+1, i+1) + inv_va; 
				gsl_matrix_set(ginvOmega, i+1, i+1, tmp);
			}
			if(m_df == 2) 
			{
				tmp = gsl_matrix_get(ginvOmega, 2*i+1, 2*i+1) + inv_va; 
				gsl_matrix_set(ginvOmega, 2*i+1, 2*i+1, tmp);
				tmp = gsl_matrix_get(ginvOmega, 2*i+2, 2*i+2) + inv_vd; 
				gsl_matrix_set(ginvOmega, 2*i+2, 2*i+2, tmp);
			}
		}   // add prior on diagnal. 
		
		gsl_linalg_cholesky_decomp(ginvOmega); 
		double logdet = 0; 
		for (int i = 0; i < col; i++)
			logdet += log(fabs(gsl_matrix_get(ginvOmega, i, i))); 
		logdet *= 2.0; 
		
		//y^t X;
		gsl_blas_dgemv(CblasTrans, 1.0, gXX, gph, 0.0, gyx); 
		gsl_linalg_cholesky_solve(ginvOmega, gyx, res); //res = omega %*% gyx; 
		
//		cout << "logdet: " << logdet << endl; 
		
//		cout << "gyx : " ; 
//		for (int j = 0; j < col; j++)
//			cout << gsl_vector_get(gyx, j) << " "; 
//		cout << endl; 
														  
		double bob = 0;
		gsl_blas_ddot(res, gyx, &bob); 
		
		double tau = sqrt((yy - bob) / ni); 
		for (int i = 0; i < col; i++)
			m_beta.at(i) += (gsl_vector_get(res, i) / tau); 

		double ph_2 = gsl_vector_get(gyx, 0);
		ph_2 *= ph_2; 
		
		double tlast = (yy - bob) / (yy - ph_2/ni);
//		cout << "tlast: " << tlast << endl; 
											  
		double logvavdn = ns * log((double) inv_va) + log((double)ni);    
		if (m_df == 2) 
			logvavdn += ns * log((double) inv_vd); 
		double logbf = -0.5 * logdet + 0.5 * logvavdn - ni * 0.5 * log((double)tlast); 
		gsl_vector_set(bf, p, logbf); 
//		cout << "lobgf: " << logbf << endl; 
	}

	for (int i = 0; i < col; i++)
		m_beta.at(i) /= repeat; 
	double maxbf = gsl_vector_max(bf); 
	gsl_vector_add_constant(bf, -maxbf); 
	double factor = 0; 
	for (int i = 0; i < repeat; i++)
		factor += exp(gsl_vector_get(bf, i)); 
	factor /= (double) repeat; 
	double logbf = maxbf + log(factor); 
		
	gsl_matrix_free(gXtX);
	gsl_matrix_free(ginvOmega);
	gsl_vector_free(gyx);
	gsl_vector_free(res);
	gsl_vector_free(bf); 
	return (logbf);  
}

double ModelnData::calc_bf(real sigma_a, real sigma_d, vector<real>& phval, vector<int>& curindex, real * genotype)
{  
	int ni = (int) curindex.size(); 
	int ns = 1; 
	int col = m_df + 1; 
	
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
	{
		int pos = curindex.at(i); 
		gsl_vector_set(gph, i, phval.at(pos)); 
	}

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		gsl_matrix_set(gXX, i, 0, 1); 
		gsl_matrix_set(gXX, i, 1, genotype[i]);
		if (m_df == 2)
			gsl_matrix_set(gXX, i, 2, (genotype[i] == 1 ? 1 : 0));
	}
	
	real logbf = bf_core(inv_va, inv_vd, ni, ns, gXX, gph);  
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (logbf);    	
}
	

double ModelnData::calc_bf(real sigma_a, real sigma_d, real * phenoval, real * genotype, int ni)
{  
	if (ni == 0) 
	{
		m_beta.clear(); 
		m_beta.push_back(0); 
		m_beta.push_back(0); 
		m_beta.push_back(0); 
		return 0.0;  
	}
	int ns = 1;
	int col = m_df + 1; 
	
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, phenoval[i]); 

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		gsl_matrix_set(gXX, i, 0, 1); 
		gsl_matrix_set(gXX, i, 1, genotype[i]);
		if(m_df == 2)
			gsl_matrix_set(gXX, i, 2, (genotype[i] == 1 ? 1 : 0));
	}
	
	real logbf = bf_core(inv_va, inv_vd, ni, ns, gXX, gph);  
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (logbf);    	
}

double ModelnData::calc_bf_mgt(real sigma_a, real sigma_d, vector<real>* phval, vector<int>* base, vector<int>* phdex, real* mgt)
{
	m_df = 1; 
	int ni = (int) phdex->size();
	if(ni == 0) return 1.0;  
	int ns = 1; 
	int col = m_df + 1;
	
	gsl_vector * gph = gsl_vector_alloc(ni); 
	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		gsl_vector_set(gph, i, phval->at(phdex->at(i))); 
		int j = base->at(i); 
		gsl_matrix_set(gXX, i, 0, 1); 
		gsl_matrix_set(gXX, i, 1, mgt[j]);
	}           

	///////////////////////////////////////////////////////
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);
	
	real logbf = bf_core(inv_va, inv_vd, ni, ns, gXX, gph);  
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (logbf);    	
}

double ModelnData::calc_bf_mean(real sigma_a, real sigma_d, vector<real>* phval, vector<int>* base, vector<int>* phdex, real** prob_snp)
{
	int ni = (int) phdex->size();
	if(ni == 0) return 1.0;  
	int ns = 1; 
	int col = m_df + 1;
	
	gsl_vector * gph = gsl_vector_alloc(ni); 
	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		gsl_vector_set(gph, i, phval->at(phdex->at(i))); 
		int j = base->at(i); 
		gsl_matrix_set(gXX, i, 0, 1); 
		gsl_matrix_set(gXX, i, 1, -prob_snp[j][1] + 2.0 * (1-prob_snp[j][0]));
        if (m_df == 2) 
			gsl_matrix_set(gXX, i, 2, prob_snp[j][1]);
	}           

	int prior = m_num; //0: default, 1: maf; 
	
	///////////////////////////////////////////////////////
	if(prior == 1)
	{      
		real af = m_current_maf; 
		double var = sqrt(af * (1.0 - af)); 
		sigma_a =  sigma_a / var;
		sigma_d =  sigma_d / var; 
	}   // maf prior; 

	if(prior == 2)
	{
	    real mean = 0;
   	    for (int i = 0; i < ni; i++)
  	        mean += gsl_matrix_get(gXX, i, 1);
		mean /= (real) ni;
		real var = 0;
		for (int i = 0; i < ni; i++)
		{
			double temp = gsl_matrix_get(gXX, i, 1) - mean;
			var += temp * temp;
		}
		var = sqrt(var); 
		sigma_a = sigma_a / var; 
		sigma_d = sigma_d / var; 
	}
						
	///////////////////////////////////////////////////////
	real inv_va = 0.0; 
	real inv_vd = 0.0;
	if(sigma_a > 0) 
		inv_va = 1./(sigma_a * sigma_a); 
	if(sigma_d > 0) 
		inv_vd = 1./(sigma_d * sigma_d);
	
	real logbf = bf_core(inv_va, inv_vd, ni, ns, gXX, gph);  
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (logbf);    	
}

real ModelnData::f_stat(real * phenoval, real * genotype, int ni)
{  
	if (ni == 0) return 0.0; 
	int col = m_df + 1; 

	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, phenoval[i]); 

	double yy = 0; 
	gsl_blas_ddot(gph, gph, &yy); 

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		gsl_matrix_set(gXX, i, 0, 1); 
		if(genotype[i] < 0 || genotype[i] > 2)                         
			cout << "-bimbam: illegal genotype" << endl;
		gsl_matrix_set(gXX, i, 1, genotype[i]);
		if(m_df == 2) 
			gsl_matrix_set(gXX, i, 2, (genotype[i] == 1 ? 1 : 0));
	}
	
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gXX, gXX, 0.0, ginvOmega); 
	

	gsl_linalg_cholesky_decomp(ginvOmega); 
	double logdet = 0; 
	for (int i = 0; i < col; i++)
		logdet += log(fabs(gsl_matrix_get(ginvOmega, i, i))); 
	logdet *= 2.0; 
	
	//y^t X;
	gsl_vector * gyx = gsl_vector_alloc(col); 
    gsl_blas_dgemv(CblasTrans, 1.0, gXX, gph, 0.0, gyx); 
//	cout << "gyx : " ; 
//	for (int j = 0; j < col; j++)
//		cout << gsl_vector_get(gyx, j) << " "; 
//	cout << endl; 

	gsl_vector * res = gsl_vector_alloc(col); 
	gsl_linalg_cholesky_solve(ginvOmega, gyx, res); //res = omega %*% gyx; 
    
	gsl_blas_dgemv(CblasNoTrans, 1.0, gXX, res, -1.0, gph); //gph = x %*% res - gph; 
	double rss = 0; 
	gsl_blas_ddot(gph, gph, &rss); 
	if(isnan(rss)) return (1e-100); 
//	cout << "rss = " << rss << endl; 
			
	double ph_2 = gsl_vector_get(gyx, 0);
	ph_2 *= ph_2;
	double rssnull = yy - ph_2 / ni; 
	
//	cout << endl << rss << "\t" << rssnull << endl; 
	double fstat = (rssnull - rss) * m_df / (rss / (ni - 2));

//	gsl_permutation_free(perm); 
	gsl_matrix_free(ginvOmega);
	gsl_vector_free(gyx);
	gsl_vector_free(res);
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (fstat);    	
}

real ModelnData::f_stat_mean(vector<real>* phval, vector<int>* base, vector<int>* curindex, real ** prob_snp)
{
	int ni = (int) curindex->size();
	if (ni == 0) return 0.0; 
	int col = m_df + 1; 
//	cout << ni << endl; 
	
	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
		gsl_vector_set(gph, i, phval->at(curindex->at(i))); 
	
	double yy = 0; 
	gsl_blas_ddot(gph, gph, &yy); 

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
		int j = base->at(i); 
		gsl_matrix_set(gXX, i, 0, 1); 
		gsl_matrix_set(gXX, i, 1, prob_snp[j][1] + 2.0 * prob_snp[j][2]);
        if(m_df == 2)
			gsl_matrix_set(gXX, i, 2, prob_snp[j][1]);
	}

	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gXX, gXX, 0.0, ginvOmega); 
	
	//y^t X;
	gsl_vector * gyx = gsl_vector_alloc(col); 
    gsl_blas_dgemv(CblasTrans, 1.0, gXX, gph, 0.0, gyx); 
	
	gsl_linalg_cholesky_decomp(ginvOmega); 
	gsl_vector * res = gsl_vector_alloc(col); 
	gsl_linalg_cholesky_solve(ginvOmega, gyx, res); //res = omega %*% gyx; 

	gsl_blas_dgemv(CblasNoTrans, 1.0, gXX, res, -1.0, gph); //gph = x %*% res - gph; 
	double rss = 0; 
	gsl_blas_ddot(gph, gph, &rss); 
	if(isnan(rss)) return (1e-100); 
//	cout << "rss = " << rss << endl; 
			
	double ph_2 = gsl_vector_get(gyx, 0);
	ph_2 *= ph_2;
	double rssnull = yy - ph_2 / ni; 
	
//	cout << endl << rss << "\t" << rssnull << endl; 
	double fstat = (rssnull - rss) * m_df / (rss / (ni - 2));

	gsl_matrix_free(ginvOmega);
	gsl_vector_free(gyx);
	gsl_vector_free(res);
	gsl_matrix_free(gXX);
	gsl_vector_free(gph); 
	return (fstat);    	
}

int ModelnData::compatibleQ(pair<char, char> p1, pair<char, char> p2)
{   
	pair<char, char> empty;
	if(p1 == empty || p2 == empty) return 1; 
	map<char, char> cc; 
	cc['A'] = 'T'; cc['T'] = 'A';
	cc['C'] = 'G'; cc['G'] = 'C';
	cc['N'] = 'N'; cc['?'] = '?'; 
	cc['+'] = '-'; cc['-'] = '+'; 

	int ac1 = (p1.first != 'N' ? 1 : 0) + (p1.second != 'N' ? 1 : 0); 
	int ac2 = (p2.first != 'N' ? 1 : 0) + (p2.second != 'N' ? 1 : 0); 
	
	int action;   //0: do nothing; 1:flip the second; -1:error; 
	if(ac1 == 1 && ac2 == 1) 
	{
		if(p1.first == p2.first) 
		{
			action = 0;          //no flip; 
		}
		else 
		{
			action = 1;         //flip second one; 
		}   // don't care about major minor here, will take care of it at 012 stage. 2-gt; 
	}

	else if (ac1 == 2 && ac2 == 2) 
	{
		if(p1.first == p2.first && p1.second == p2.second)  //match;
		{
	        action = 0; 
		}
		else if (p1.first == p2.second && p1.second == p2.first) //match allele type, but not maf; 
		{
			action = 1; 
		}   //should generate a warning; 
		else if (p1.first == cc[p2.first] && p1.second == cc[p2.second]) //flip allele type, maf match; 
		{
			action = 0; 
		}   //should generate a warning; 
		else if (p1.first == cc[p2.second] && p1.second == cc[p2.first]) //flip allele type, maf not match; 
		{
			action = 1; 
		}   //should generate a warning; 
		else
			action = -1; 
	}
	else if (ac1 == 1 && ac2 == 2) 
	{
		if(p1.first == p2.first) //match 
		{
			action = 0; 
		}
		else if (p1.first == p2.second) 
		{
			action = 1; 
		}
		else if (p1.first == cc[p2.first])  //flip
		{
			action = 0; 
		}
		else if (p1.first == cc[p2.second]) 
		{
			action = 1; 
		}
		else 
			action = -1; 
	}
	else if (ac1 == 2 && ac2 == 1) 
	{
		if (p2.first == p1.first) 
		{
			action = 0; 
		}
		else if (p2.first == p1.second) 
		{
			action = 1; 
		}
		else if (cc[p2.first] == p1.first) 
		{
			action = 0; 
		}
		else if (cc[p2.first] == p1.second) 
		{
			action = 1; 
		}               
		else 
			action = -1;  
	}
	else    //one of them is missing; 
	{
		action = -1; 
	}
	return action; 
}

char ModelnData::merge_mmSNP(class mmSNP & ref, class mmSNP & s)
{   
	//grant:
//	cout << ref.A << ref.countA << " | " << ref.B << ref.countB << endl; 
//	cout << s.A << s.countA << " | " << s.B << s.countB << endl; 
	map<char, char> cc; 
	cc['A'] = 'T'; cc['T'] = 'A';
	cc['C'] = 'G'; cc['G'] = 'C';
	cc['N'] = 'N'; cc['?'] = '?'; 
	cc['+'] = '-'; cc['-'] = '+'; 
    
    char coding = 'N';  
	//both mono;   //need to update ref allele type; 
	//both poly;   //match allele type, update allele count; 
	//ref mono, new poly;  //direct match? Yes (update ref allele) : No (flip and update ref allele); flippable? 
	//ref poly, new mono; 
	int ac1 = (ref.A != 'N' ? 1 : 0) + (ref.B != 'N' ? 1 : 0); 
	int ac2 = (s.A != 'N' ? 1 : 0) + (s.B != 'N' ? 1 : 0); 

	//grant:
//	cout << ac1 << " \t" << ac2 << endl; 
	
	if(ac1 == 1 && ac2 == 1) 
	{
		coding = ref.A; 
		if(ref.A == s.A) 
		{
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else 
		{
			ref.B = s.A; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ; 
		}   // don't care about major minor here, will take care of it at 012 stage. 2-gt; 
	}

	else if (ac1 == 2 && ac2 == 2) 
	{
		if(ref.A == s.A && ref.B == s.B)  //match;
		{
			coding = ref.A; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else if (ref.A == s.B && ref.B == s.A) //match allele type, but not maf; 
		{
			coding = ref.A; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}   //should generate a warning; 
		else if (ref.A == cc[s.A] && ref.B == cc[s.B]) //flip allele type, maf match; 
		{
			coding = s.A; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}   //should generate a warning; 
		else if (ref.A == cc[s.B] && ref.B == cc[s.A]) //flip allele type, maf not match; 
		{
			coding = s.B; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}   //should generate a warning; 
		else
			coding = 'X'; 
	}
	else if (ac1 == 1 && ac2 == 2) 
	{
		if(ref.A == s.A) //match 
		{
			coding = s.A; 
			ref.B = s.B;  //update allele type; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else if (ref.A == s.B) 
		{
			coding = s.B; 
			ref.B = s.A; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}
		else if (ref.A == cc[s.A])  //flip
		{
			coding = cc[s.A]; 
			ref.B = cc[s.B]; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else if (ref.A == cc[s.B]) 
		{
			coding = cc[s.B]; 
			ref.B = cc[s.A]; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}
		else 
			coding = 'X'; 
	}
	else if (ac1 == 2 && ac2 == 1) 
	{
		if (s.A == ref.A)         //reference allele match the first allele; 
		{
			coding = s.A; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else if (s.A == ref.B)    //fist allele match the second reference allele; 
		{
			coding = ref.B; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}
		else if (cc[s.A] == ref.A) 
		{
			coding = s.A; 
			ref.countA += s.countA; 
			ref.countB += s.countB; 
			ref.countQ += s.countQ;   
		}
		else if (cc[s.A] == ref.B) 
		{
			coding = ref.B; 
			ref.countA += s.countB; 
			ref.countB += s.countA; 
			ref.countQ += s.countQ;   
		}               
		else 
			coding = 'X'; 
	}
	else if(ac1 == 0)    
	{
		coding = s.A; 
		ref.countA += s.countA; 
		ref.countB += s.countB; 
		ref.countQ += s.countQ;   
	}
	else if(ac2 == 0) 
	{
		coding = ref.A; 
	}

	//grant: 
//	cout << coding << " : " << ref.A << ref.countA << " \t" << ref.B << ref.countB << endl; 

	return coding; 
}

int ModelnData::read_rs2pos(string fn, long start_pos, long end_pos, long flank) 
{
	ifstream infile; 
	char delimit[] = ";, \t";
	streambuf * pbuf;
	int snpcount = 0; 

    map<string, long> :: iterator liter; 
    map<string, int> :: iterator iter; 
	
	//first, read in the rs pos file; and put into vRsPos;
	infile.open(fn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		return 0; 
	}
	pbuf = infile.rdbuf(); 	
	string line; 
	line.assign(getline(pbuf)); 
	while (!line.empty()) 
	{
		line.append(" 0 "); 
		char * res = strtok((char *)line.c_str(), delimit); 
		string rs(res); 
		res = strtok(NULL, delimit); 
		long pos = atol(res);
		res = strtok(NULL, delimit); 
		int chr = atoi(res);
		
		if (end_pos > 0) 
		{
			if(pos > start_pos - flank && pos < end_pos + flank)
			{
				liter = mapRs2pos.find(rs); 
				if(liter == mapRs2pos.end())
					mapRs2pos.insert(pair<string, long> (rs, (long)pos));
				else 
					liter->second = pos; 
				
				iter = mapRs2chr.find(rs); 
				if(iter == mapRs2chr.end())
					mapRs2chr.insert(pair<string, int> (rs, chr)); 
				else 
					iter->second = chr; 
				snpcount++;
			}
		}
		else 
		{
			mapRs2pos.insert(pair<string, long> (rs, pos));
			mapRs2chr.insert(pair<string, int> (rs, chr)); 
			snpcount++;
		}
		line.assign(getline(pbuf)); 
	}
	infile.close();
	return snpcount; 
}

bool rscomp(pair<string, pair<long, int> > a, pair<string, pair<long, int> > b)
{
	if(a.second.second == b.second.second)
		return (a.second.first < b.second.first); 
	else 
		return(a.second.second < b.second.second); 
}

int ModelnData::read_bimbam_genotype(int mcmc, long start_pos, long end_pos)
{   
	if(vGin.size() == 0 || vGin.size() != vPin.size()) 
	{
		fplog << "## BIMBAM: genotype phenotype files not match" << endl; 
		cout << "-bimbam: genotype phenotype files not match" << endl; 
		return 0; 
	}
	fplog << "## BIMBAM: Number of position files = " << vPos.size() << endl; 
	unsigned readpos = 0; 
	for(unsigned f = 0; f < vPos.size(); f++)
	{
		fnRsPos.assign(vPos.at(f)); 
		readpos += read_rs2pos(fnRsPos, start_pos, end_pos, nGeneFlank); 
	}
	if(readpos == 0 && vGin.size() > 1) 
	{
		fplog << "## Quit: BIMBAM requires SNP position file(s) for multiple genotypes" << endl; 
		cout << "-bimbam: requires position file(s) for multiple genotypes" << endl; 
		safe_exit(); 
	}   //initialize mapRs2pos; 
	if(readpos < vPos.size()) 
	{
		fplog << "## WARNING: One of position files failed" << endl; 
		cout << "-bimbam: one of position files failed" << endl; 
	}
	fplog << "## BIMBAM: Position files contain " << mapRs2pos.size() << " records." << endl; 
	if(vPos.size() > 0)
		cout << "-bimbam: position files contain " << mapRs2pos.size() << " records." << endl; 

	//first pass, go over all the files collecting summary for each SNP. 
	//allele type, allele count, etc. treating phased panel as unphased (but get the correct nHap, nDip);  
	//merge SNPs based on summary. generate a hash keyed by rs whose content has reference allele
	//for different files (so that genotypes are the counts of reference alleles). If the reference
	//allele is '?', the whole SNP for that file were to set 'QQ'; 
	//sort the key according to the physical position, to get vsRsnum; 
	//Allocate memory for each individual, i.e., we have nIndiv x nLoci; 
	 
	nHap = 0; 
	nDip = 0; 
	nLoci = 0; 
	map<string, int> mCohortRsQ; //cohort SNPs; 
	map<string, int> mrs2yes; 
	map<int, map<string, class mmSNP> > mfile2map; 
	vector<int> vindploid; 

	for (unsigned nf = 0; nf < vGin.size(); nf++)
	{   
		int filerows = 0; 
		int filecols = 0; 
		int filecols_diff = 0; 
		int isPanel = 0; 
		if (vPin.at(nf).compare("0") == 0)  
			isPanel = 1; 
		map<string, class mmSNP> rs2ss; //rs to snp summary; 
		
		ifstream infile; 
		streambuf * pbuf;
		char delimit[] = ";, \t";
		infile.open(vGin.at(nf).c_str(), ios::in);
		if(!infile.is_open()) 
		{
			fplog << "Quit: BIMBAM cannot open genotype file: " << vGin.at(nf) << endl; 
			cout << "-bimbam: cannot open genotype file: " << vGin.at(nf) << endl; 
			safe_exit(); 
		} 
		pbuf = infile.rdbuf();
		// open file; 

		string line; 
		int ni = 0; 
		{
			line.assign(getline(pbuf)); 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res != NULL) ni = atoi(res);
			
			string::size_type loc; 
			loc = line.find(".", 0);
			if (loc != string::npos || ni <= 1)
			{
				fplog << "## Quit: Confused -p with -g ? " << endl; 
				cout << "-bimbam: confused -p with -g ? " << endl; 
				safe_exit(); 
			}
			
			loc = line.find("=", 0);
			if( loc != string::npos ) 
				vFilePloid.push_back(1); 
			else 
				vFilePloid.push_back(2); 
		}   //number of individuals;
		
		int ns = 0; 
		{
			line.assign(getline(pbuf)); 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res != NULL) ns = atoi(res);
		}	//number of snps;
				
		string tline(getline(pbuf)); 
		line.assign(tline); 
		char * res = strtok((char*)tline.c_str(), delimit); 
		string head(res);
		if(head.compare("IND") == 0 || head.compare("rsID") == 0)
		{
			line.assign(getline(pbuf)); 
		}   //individual ids
			
		while(!line.empty())
		{   
			//printf("%s\n", line.c_str()); 
			filerows ++; 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res == NULL) break; 
			string rs;
			rs.assign(res);  
			mrs2yes[rs] = 1;
			if(isPanel) mPanelRsQ[rs] = 1; 
			else mCohortRsQ[rs] = 1; 
			
			if(end_pos > 0) 
			{
				map<string, long> :: iterator iter; 
				iter = mapRs2pos.find(rs); 
				if(iter == mapRs2pos.end() || iter->second < 0) 
				{
					line.assign(getline(pbuf));
					continue; 
				}
			}   // rs; 
			
			map<char, int> a2c;
			map<char, int> :: iterator iter; 
			vector<char> snpgt; 
			filecols = 0; 

			int tni = ni; 
			for (int j = 0; j < tni; j++)    // while(1); 
			{          
				string tmp;
				res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				tmp.assign(res);
				for (int i = 0; i < 2; i++)
				{
					if (tmp.at(i) == 'N' || tmp.at(i) == '0' || tmp.at(i) == '?') 
					{
						iter = a2c.find('N'); 
						if(iter == a2c.end()) 
							a2c['N'] = 1; 
						else 
							iter->second ++; 
					}
					else 
					{
						iter = a2c.find(tmp.at(i)); 
						if(iter == a2c.end()) 
							a2c[tmp.at(i)] = 1; 
						else 
							iter->second ++; 
					}
				}                   
				filecols++; 
			}   // snp summary in a2c;    //A C G T ? + -
			if(filecols != ni) filecols_diff++; 
			
			class mmSNP ss; 
			iter = a2c.find('N'); 
			if(iter == a2c.end()) ss.countQ = 0; 
			else ss.countQ = iter->second; 
			
			vector<char> type; 
			for (iter = a2c.begin(); iter != a2c.end(); ++iter)
			{
				if(iter->first != 'N') type.push_back(iter->first); 
			}   //type only contains ACGT+-; 
			
			if (type.size() == 0) 
			{
				; //fplog << "## Warning: In file " << vGin.at(nf) << " snp " << rs << " is null || "; 
			}
			else if (type.size() > 2) 
			{
				fplog << "## Warning: In file " << vGin.at(nf) << " snp " << rs << " has more than 2 alleles" << endl; 
				ss.A = 'N'; 
				ss.B = 'N'; 
			}

			else if(type.size() == 1) 
			{
				ss.A = type.at(0); 
				ss.countA = a2c[ss.A]; 
				ss.B = 'N'; 
				ss.countB = 0; 
			}   //note: this is consistant with two allele case; 
			else 
			{
				double n0 = (double) a2c[type.at(0)]; 
				double n1 = (double) a2c[type.at(1)]; 
			
				if(n0 > n1) 
				{
					ss.A = type.at(0); 
					ss.countA = n0; 
					ss.B = type.at(1); 
					ss.countB = n1; 
				}
				else 
				{
					ss.A = type.at(1); 
					ss.countA = n1; 
					ss.B = type.at(0); 
					ss.countB = n0; 
				}
					   
				//note: A is major allele; 
			}   //translate a2c into mmSNP format; 
			a2c.clear(); 
			map<string, class mmSNP> :: iterator mmiter; 
			mmiter = rs2ss.find(rs); 
			if (mmiter != rs2ss.end())
			{
				fplog << "## WARNING: Two SNPs have the same rsnum " << mmiter->first << endl; 
				cout << "-bimbam: two snp have same rsnum " << mmiter->first << endl; 
			}
			
			rs2ss.insert(pair<string, class mmSNP> (rs, ss)); 
//	cout << rs << " " << rs2ss.size() << " " <<  ss.A << ss.countA << " " << ss.B << ss.countB << " " << ss.countQ << endl; 
			
			if(readpos == 0) 
			{
				mapRs2pos.insert( pair<string, long> (rs, (long) rs2ss.size()));  
				mapRs2chr.insert( pair<string, int> (rs, 0)); 
			}
			//if no rs-pos file and only one genotype file, use snp order as pos.
			line.assign(getline(pbuf));
		}

		infile.close();
		mfile2map[nf] = rs2ss; 
		
		if(filerows != ns) 
		{
			fplog << "## WARNING: Number of SNPs in file is different from the second line" << endl; 
			cout << "-bimbam: actual snp number differ from the second line" << endl; 
		}
		if(filecols != ni) 
		{
			fplog << "## WARNING: Number of ind's is different from the first line " << filecols << " " << ni << endl; 
			cout << "-bimbam: actual indiv number differ from the first line " << filecols << " " << ni << endl; 
		}
		if(filecols_diff > 0)
		{
			fplog << "## WARNING: " << filecols_diff << " SNPs has different # of ind's from the first line" << endl; 
			cout << "-bimbam: in " << filecols_diff << " snps, individual differ from the first line" << endl; 
		}
		if (ni > filecols) ni = filecols;
		
		if(vFilePloid.at(nf) == 1)
		{
			vFileIndiv.push_back(2 * ni);
			nHap += 2 * ni; 
			for(int i = 0; i < 2 * ni; i++)
				vindploid.push_back(1); 
		}
		else 
		{
		
			vFileIndiv.push_back(ni);
			nDip += ni; 
			for(int i = 0; i < ni; i++)
				vindploid.push_back(2); 
		}
		fplog << "## BIMBAM: File " << nf << " has " << vFileIndiv.at(nf) << " ind's and " \
			<< rs2ss.size() << " SNPs" << endl; 
		cout << "-bimbam: file " << nf << " has " << vFileIndiv.at(nf) << " individual and " \
			<< rs2ss.size() << " snps" << endl; 
		rs2ss.clear(); 
	}   // SNP summary data;
	
//		{
//			map<string, class mmSNP> :: iterator iter; 
//			for (int nf = 0; nf < mfile2map.size(); nf++)
//			{
//				for (iter = mfile2map[nf].begin(); iter != mfile2map[nf].end(); ++iter)
//				{
//					cout << nf << "\t" <<  iter->first << "\t" << iter->second.countA << "\t" << iter->second.countB;   
//					cout << "\t" << iter->second.countQ << endl; 
//				}
//			}
//		}   //test if the mmSNP build ok; 
	////////////////////////////////////////////////////////////////////////////
	// now perform merger;   //snp by snp; 

	map<string, int> ::iterator m_iter; 
	map<string, vector<char> > mrs2coding;   //rs to reference allele coding for each file;  
	for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); ++m_iter)
	{
		string rs; 
		rs.assign(m_iter->first); 
		vector<char> coding; 
		class mmSNP ref; 
		int bingle = 0; 
		for (unsigned nf = 0; nf < vGin.size(); nf++)
		{
			map<string, class mmSNP> :: iterator iter; 
			iter = mfile2map[nf].find(rs); 
			if(iter == mfile2map[nf].end())
			{
				ref.countQ += vFileIndiv.at(nf) * vFilePloid.at(nf); 
				coding.push_back('X');
				continue; 
			}
			if(ref.A == 'N' && ref.B == 'N')  // first one; 
			{
				ref.assign(iter->second); 
				ref.A = (iter->second).A; 
				ref.B = (iter->second).B; 
				ref.countA = (iter->second).countA; 
				ref.countB = (iter->second).countB; 
				ref.countQ += (iter->second).countQ; 
//					for (int i = 0; i < nf; i++)
//						ref.countQ += vFileIndiv.at(i); 
				// add the missing count in previous files;
				coding.push_back(ref.A);
				continue; 
			}
			
			if(iter->second.A == 'X') 
			{
				bingle = 1; 
				break; 
			}
			char cr = merge_mmSNP(ref, iter->second);
			if(cr == 'X') 
			{
				bingle = 1; 
				break; 
			}
			else 
				coding.push_back(cr); 
		}
		
//   			cout << rs << "\t" << ref.A << ref.B << "\t" << ref.countA << "\t" << ref.countB;  
//   			cout << "\t" << ref.countQ << endl; 
		
		if (bingle == 1) 
			m_iter->second = 0; //remove this SNP; 
		
		mrs2coding[rs] = coding;                  
		mapRs2mm[rs] = pair<char, char> (ref.A, ref.B); 
		double countAB = ref.countA + ref.countB; 
		double missrate = (double ) ref.countQ / (countAB + ref.countQ); 
		//cout << missrate << "\t" << ref.countA << "\t" << ref.countB << "\t" << ref.countQ << endl; 
		if( missrate > m_exclude_miss)        //missing test; 
			mrs2yes[rs] = -1; 
		else if(countAB < 1e-6) 
		{
			mapRs2maf[rs] = 0.0; 
			mapRs2var[rs] = 1e-10;
		}
		else 
		{
			double af = (double) (ref.countB) / countAB;  
			mapRs2maf[rs] = af; 
			mapRs2var[rs] = 2.0 * af * (1.0 - af); 
		}
			
		vector<char> ().swap(coding); 
	}
	
	int count_pos = 0;  //failed to have a position; 
	int count_am = 0;   //failed to match allele type; 
	int count_maf = 0;  //failed due to maf too extreme; 
	int count_miss = 0; 
	vector<pair<string, pair<long, int> > > vp;
	map<string, long> :: iterator pos_iter; 
	map<string, int> :: iterator chr_iter; 
	for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); m_iter++)
	{
		string rs(m_iter->first);
		if(m_iter->second == -1) 
		{
			count_miss ++; 
			continue; 
		}
		if(m_iter->second == 0) 
		{
			count_am ++; 
			continue; 
		}
		if(mapRs2maf[rs] <= m_exclude_maf) 
		{
			count_maf ++;    
//			fplog << rs << "\t" << mapRs2maf[rs] << endl; 
			m_iter->second = 0;  //mark this SNP to be excluded. 
			continue; 
		}
		if(m_exclude_nopos && mapRs2pos[rs] == 0)   //if snp no position; 
		{
			count_pos ++; 
			m_iter ->second = 0; //mark this SNP to be excluded. 
			continue; 
		}
		pair<string, pair<long, int> > tmp;
		tmp.first.assign(m_iter->first);
		pos_iter = mapRs2pos.find(m_iter->first);
		if(pos_iter != mapRs2pos.end())
			tmp.second.first = pos_iter->second;
		else 
			tmp.second.first = 0; 

		chr_iter = mapRs2chr.find(tmp.first); 
		if(chr_iter != mapRs2chr.end())
			tmp.second.second = chr_iter->second; 
		else
			tmp.second.second = 0; 
								   
		vp.push_back(tmp); 
	}
	if(count_miss > 0) 
	{
		fplog << "## BIMBAM: Exclude " << count_miss << " SNPs due to miss proportion > " << m_exclude_miss << endl;
		cout << "-bimbam: exclude " << count_miss << " snps due to miss proportion > " << m_exclude_miss << endl;
	}
	if(count_am > 0) 
	{
		fplog << "## BIMBAM: Exclude " << count_am << " SNPs due to failure to match betweeen files." << endl;
		cout << "-bimbam: exclude " << count_am << " snps due to failure to match betweeen files." << endl;
	}
	if(count_maf > 0) 
	{
		fplog << "## BIMBAM: Exclude " << count_maf << " SNPs due to maf <= " << m_exclude_maf << endl;
		cout << "-bimbam: exclude " << count_maf << " snps due to maf <= " << m_exclude_maf << endl;
	}
	if(count_pos > 0) 
	{
		fplog << "## BIMBAM: Exclude " << count_pos << " SNPs due to no position information" << endl;
		cout << "-bimbam: exclude " << count_pos << " snps due to no position information" << endl;
	}
	stable_sort(vp.begin(), vp.end(), rscomp); 
	//sort the rs pos vector in order of chr && pos;
	map<string, int> mrs2index; //the column of each snp; 
	vsRsnum.clear(); 
	for(int i = 0; i < (int) vp.size(); i++)
	{
		vsRsnum.push_back(vp.at(i).first);
		mrs2index[vp.at(i).first] = i; 
	}
	vector<pair<string, pair<long, int> > >().swap(vp); 
	mrs2yes.clear(); 
	nLoci = (int) vsRsnum.size(); 
	nIndiv = nDip + nHap; 
	nCohort = nDip; 
	if(nLoci == 0) 
	{
		cout << "-bimbam: no valid snp" <<  endl; 
		safe_exit(); 
	}
	//now vsRsnum hold ID of all valid SNPs in order; 
	
	vPARCrs.clear(); 
	vPARCpos.clear(); 
	for (unsigned i = 0; i < vsRsnum.size(); i++)
	{
		m_iter = mCohortRsQ.find(vsRsnum.at(i));
		if(m_iter != mCohortRsQ.end())  
		{
			vPARCrs.push_back(m_iter->first); 
			vPARCpos.push_back(i); 
		}
	}

	mCohortRsQ.clear(); 
	////////////////////////////////////////////////////////////////////////////////////
	
	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".snpdata.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(outfile.is_open()) 
	{
		outfile << "## af is the allele freq for A" << endl;  
		outfile << "rs\t A\t B\t af \t chr \t pos" << endl;  
		for (int m = 0; m < nLoci; m++)
		{
			string rs = vsRsnum.at(m); 

//			cSNP tsnp(rs, mapRs2chr[rs], mapRs2pos[rs], m, mapRs2mm[rs], maf); 
//			vsnp.push_back(tsnp); 
			char buf[100]; 
			sprintf(buf, "%-s\t", rs.c_str()); 
			outfile << buf; 
			int tac = m_allele_coding_mode; 
//				if(tac > 1) 
//				{
//					if(mapRs2mm[rs].first < mapRs2mm[rs].second) tac = 1; 
//					else 
//						tac = 0; 
//				}   //only for gain; ignore for now; 
			if (tac == 0) 
			{
				outfile << mapRs2mm[rs].second << "\t" << mapRs2mm[rs].first << "\t";
				sprintf(buf, "%.3f\t", mapRs2maf[rs]); 
			}
			else //tac == 1;  
			{
				outfile << mapRs2mm[rs].first << "\t" << mapRs2mm[rs].second << "\t";
				sprintf(buf, "%.3f\t", 1.0 - mapRs2maf[rs]); 
			}

			
			outfile << buf;
			if(mapRs2chr[rs] < 0)
				outfile << " NA \t"; 
			else
				outfile << mapRs2chr[rs] << "\t"; 
			outfile << mapRs2pos[rs] << endl; 
		}
		outfile.close(); 
	}
	else 
		cout << "-bimbam: skip writing snpdata" << endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	pIndiv = new Individual * [nIndiv]; 

	for (int i = 0; i < nIndiv; i++)
	{
		if (vindploid.at(i) == 1) {
			pIndiv[i] = new HapInd; 
			pIndiv[i]->AllocatesnpGT(nLoci); 
			for (int m = 0; m < nLoci; m++)
				pIndiv[i]->SetsnpGT(m, char('0' + QQ)); 
		} else {
			pIndiv[i] = new DipInd;
			pIndiv[i]->AllocatesnpGT(nLoci); 
			for (int m = 0; m < nLoci; m++)
				pIndiv[i]->SetsnpGT(m, char('0' + QQ)); 
		}
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////
	// go over genotype files again; 
	// for each diploid file, scan each SNP, fill in the genotype table. 
	int start_ni = 0; 
	for (int nf = 0; nf < (int)vGin.size(); nf++)
	{       
		streambuf * pbuf;
		ifstream infile; 
		char delimit[] = ";, \t";
		infile.open(vGin.at(nf).c_str(), ios::in);
		if(!infile.is_open()) 
		{
			cout << "-bimbam: cannot open genotype file: " << vGin.at(nf) << endl; 
			safe_exit(); 
		} 
		pbuf = infile.rdbuf();
		 
		{
			if (vPin.at(nf).compare("0") == 0)  
			{
				for (int i = 0; i < vFileIndiv.at(nf); i++)
					pIndiv[start_ni + i]->SetisPanel(1); 
			} 
			else 
			{
				for (int i = 0; i < vFileIndiv.at(nf); i++)
					pIndiv[start_ni + i]->SetisPanel(0); 
			}
		}
		string line; 
		line.assign(getline(pbuf));   //first line; 
		line.assign(getline(pbuf));   //second line; 
		
		string tline(getline(pbuf));   //third line; 
		line.assign(tline); 
		char * res = strtok((char*)tline.c_str(), delimit); 
		string head(res);
		if(head.compare("IND") == 0 || head.compare("rsID") == 0)
			line.assign(getline(pbuf)); 
		//possible individual ids
			
		while(line.size() > 0)
		{   
			char * res = strtok((char*)line.c_str(), delimit); 
			if(res == NULL) break; 
			string rs; 
			rs.assign(res);  
			if(end_pos > 0) 
			{
				map<string, long> :: iterator iter; 
				iter = mapRs2pos.find(rs); 
				if(iter == mapRs2pos.end() || iter->second < 0) 
				{
					line.assign(getline(pbuf));
					continue; 
				}
			}    
			int m = -1; 
			map<string, int> :: iterator s2i; 
			s2i = mrs2index.find(rs); 
			if(s2i == mrs2index.end())
			{
				line.assign(getline(pbuf));
				continue; 
			}
			else 
				m = s2i->second; ; 
			// rs number; 
			char ref = mrs2coding[rs].at(nf); //reference allele; 
//				cout << nf << " " << ref << endl; 
//				if(ref == 'X') continue; 
			int ni = start_ni; 
			for (int i = 0; i < vFileIndiv.at(nf); i++)
			{                              
				char * res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				string tmp;
				tmp.assign(res);
				if(vFilePloid.at(nf)  == 2) 
				{
					if (tmp.at(0) == 'N' || tmp.at(0) == '?' || tmp.at(1) == 'N' || tmp.at(1) == '?')
					{
						pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(0) == ref) gt++;
						if (tmp.at(1) == ref) gt++;
						if(m_allele_coding_mode == 0) gt = 2 - gt; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}           
					ni++; 
				}
				else if (vFilePloid.at(nf) == 1) 
				{
					if (tmp.at(0) == 'N' || tmp.at(0) == '?')
					{
						 pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(0) == ref) gt++;
						if(m_allele_coding_mode == 0) gt = 1 - gt; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}                                        
					ni++; 
					if (tmp.at(1) == 'N' || tmp.at(1) == '?')
					{
						pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(1) == ref) gt++;
						if(m_allele_coding_mode == 0) gt = 1 - gt; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}     
					ni++;
				}
			}
			line.assign(getline(pbuf));
		}
		infile.close();
		cout << "-bimbam: read file " << nf << " again " << endl; 
		start_ni += vFileIndiv.at(nf); 
	}

	nPanel = 0; 
	for (int i = 0; i < nIndiv; i++)
		nPanel += pIndiv[i]->GetisPanel(); 
	nCohort = nIndiv - nPanel; 

	//////////////////////////////////////////////////////////////////////////////////
	
	fplog << "## number genotype files = " << vGin.size() << endl; 
	fplog << "## number phenotype files = " << vPin.size() << endl; 
	fplog << "## number of diploid = " << nDip << endl; 
	fplog << "## number of haploid = " << nHap << endl; 
	fplog << "## number of panel individuals = " << nPanel << endl; 
	fplog << "## number of cohort = " << nCohort << endl; 
	fplog << "## number of snp = " << nLoci  << endl;
	fplog << "## number cohort only snp = " << vPARCrs.size() << endl; 

//////////////////////////////////////////////////////////////////////////////////////////

	return 1; 
}  

int ModelnData::read_bimbam_genotype_distribution(int mode, int mcmc, long start_pos, long end_pos)   
	//mode=1: mean geno; mode=2:posterior distribution; 
{   
	if(vGin.size() == 0 || vGin.size() != vPin.size()) 
	{
		fplog << "## BIMBAM: genotype phenotype files not match" << endl; 
		cout << "-bimbam: genotype phenotype files not match" << endl; 
		return 0; 
	}   //check geno-pheno files match; 
	
	fplog << "## BIMBAM: Number of position files = " << vPos.size() << endl; 
	unsigned readpos = 0; 
	for(unsigned f = 0; f < vPos.size(); f++)
	{
		fnRsPos.assign(vPos.at(f)); 
		readpos += read_rs2pos(fnRsPos, start_pos, end_pos, nGeneFlank); 
	}    //initialize mapRs2pos; 
	if(readpos == 0 && vGin.size() > 1) 
	{
		fplog << "## Quit: BIMBAM requires SNP position file(s) for multiple genotypes" << endl; 
		cout << "-bimbam: requires position file(s) for multiple genotypes" << endl; 
		safe_exit(); 
	}   
	if(readpos < vPos.size()) 
	{
		fplog << "## WARNING: One of position files failed" << endl; 
		cout << "-bimbam: one of position files failed" << endl; 
	}   //check position files; 

	fplog << "## BIMBAM: Position files contain " << mapRs2pos.size() << " records." << endl; 
	if(vPos.size() > 0)
		cout << "-bimbam: position files contain " << mapRs2pos.size() << " records." << endl; 

	{

		//first pass, go over all the files collecting summary for each SNP. 
		//allele type, allele count, etc. treating phased panel as unphased (but get the correct nHap, nDip);  
		//merge SNPs based on summary. generate a hash keyed by rs whose content has reference allele
		//for different files (so that genotypes are the counts of reference alleles). If the reference
		//allele is '?', the whole SNP for that file were to set 'QQ'; 
		//sort the key according to the physical position, to get vsRsnum; 
		//Allocate memory for each individual, i.e., we have nIndiv x nLoci; 
		 
			
		nHap = 0; 
		nDip = 0; 
		nLoci = 0; 
		map<string, int> mCohortRsQ; //cohort SNPs; 
		map<string, int> mrs2yes; 
		map<int, map<string, class mmSNP> > mfile2map; 
		vector<int> vindploid; 

		for (unsigned nf = 0; nf < vGin.size(); nf++)
		{   
			int filerows = 0; 
			int filecols = 0; 
			int filecols_diff = 0; 
			int isPanel = 0; 
			if (vPin.at(nf).compare("0") == 0)  
				isPanel = 1; 
			map<string, class mmSNP> rs2ss; //rs to snp summary; 
			
			ifstream infile; 
			streambuf * pbuf;
			char delimit[] = ";, \t";
			infile.open(vGin.at(nf).c_str(), ios::in);
			if(!infile.is_open()) 
			{
				fplog << "Quit: BIMBAM cannot open genotype file: " << vGin.at(nf) << endl; 
				cout << "-bimbam: cannot open genotype file: " << vGin.at(nf) << endl; 
				safe_exit(); 
			} 
			pbuf = infile.rdbuf();
			// open file; 

			int ni = 0; 
			int ns = 0; 
			vFilePloid.push_back(2); 
			        
			string line; 
			line.assign(getline(pbuf)); 
			int len = min(100, (int)line.size()); 
			for (int i = 0; i < len; i++)
			{
				if((int)line.at(i) == 0) 
					line.at(i) = 'N'; 
			}   //this is because in a previous version the null allele was coded "\0"; 
			    //however, this won't hurt for the current mgt files. 

			while(!line.empty())
			{   
				//printf("%s\n", line.c_str()); 
				filerows ++; 
				char *res = strtok((char*)line.c_str(), delimit); 
				if(res == NULL) break; 
				string rs(res);  
				if(end_pos > 0) 
				{
					map<string, long> :: iterator iter; 
					iter = mapRs2pos.find(rs); 
					if(iter == mapRs2pos.end() || iter->second < 0) 
					{
						line.assign(getline(pbuf));
						continue; 
					}
				}   // rs; 

				mrs2yes[rs] = 1;
				if(isPanel) mPanelRsQ[rs] = 1; 
				else mCohortRsQ[rs] = 1; 
				
				pair<char, char> mm; 
				res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				mm.first = res[0]; 
				res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				mm.second = res[0];
					
				mapRs2mm[rs] = mm; 
				
				map<char, double> a2c;
				map<char, double> :: iterator iter; 
				vector<char> snpgt; 
				
				filecols = 0; 
				double ex1 = 0; 
				double ex2 = 0; 
				
				//p0 = pr(major, major); 
				//p1 = pr(minor, major); 
				//p2 = pr(minor, minor); 

				//mm.first == minor; 
				//mm.second == major; 
				//note: this is the same as in the output of mean geotype, 
				//where the first allele is the minor allele of the panel; 
				while(1) 
				{    
					double g0 = 0.0; 
					double g1 = 0.0; 
					res = strtok(NULL, delimit); 
					if(res == NULL) break; 
					g0 = atof(res); 
					if(mode == 2) 
					{
						res = strtok(NULL, delimit); 
						if(res == NULL) 
						{
							cout << "-bimbam: odd number of numerical entries" << endl;
							safe_exit(); 
						}
						g1 = atof(res); 
					}
					
					if(mode == 1) 
					{
						double mg = 0; 
					
						iter = a2c.find(mm.first); 
						if(iter == a2c.end())
							a2c[mm.first] = g0;
						else 
							iter->second += g0; 
						
						iter = a2c.find(mm.second); 
						if(iter == a2c.end())
							a2c[mm.second] = 2.0 - g0;
						else 
							iter->second += 2.0 - g0; 
						
						mg = g0; 
						ex1 += mg; 
						ex2 += mg * mg; 
					}    
					else if (mode == 2) 
					{                                                                  
						double mg = 0; 
						double g2 = 1.0 - g0 - g1; 

						iter = a2c.find(mm.first); 
						if(iter == a2c.end())
							a2c[mm.first] = (2 * g2 + g1);
						else 
							iter->second += (2 * g2 + g1); 
						
						iter = a2c.find(mm.second); 
						if(iter == a2c.end())
							a2c[mm.second] = (2 * g0 + g1);
						else 
							iter->second += (2 * g0 + g1); 
						
						mg = g1 + 2.0 * g2; 
						ex1 += mg; 
						ex2 += mg * mg; 
					}
					else 
					{
						cout << "-bibmam: wrong mode in read genotype distribution" << endl; 
						safe_exit(); 
					}
					
					filecols++; 
				}   // snp summary in a2c;    //A C G T ? + -

				double tmu = ex1 / (double) filecols; 
				double gvar = (ex2 / filecols) - tmu * tmu; 
				mapRs2var[rs] = gvar; 
				
                if(ni == 0) ni = filecols;  
				if(filecols != ni) { 
					filecols_diff++; 
				    cout << rs << "\t" << filecols << "\t" << ni  << endl; 
					for (int i = 1; i < 20; i++)
						fplog << line.at(i) << " " << (int)line.at(i) << "  "; 
					fplog << endl; 
				}
				class mmSNP ss; 

				vector<char> type; 
				for (iter = a2c.begin(); iter != a2c.end(); ++iter)
				{
					if(iter->first != 'N') type.push_back(iter->first); 
				}

				if(type.size() <= 1) 
				{
					type.push_back('N'); 
				    a2c['N'] = 0; 
				}

				//grant: 
//				for (iter = a2c.begin(); iter != a2c.end(); ++iter)
//					cout << iter->first << " " << iter->second << " ||| "; 
//				cout << endl; 
				
				//make the ss.A as major and ss.B as minor; for monomorphic SNP. ss.B = 'N'. 
				//note: there is a complication here that our major/minor in mgt 
				// is referring to the panel alleles. 
				else
				{
					double n0 = a2c[type.at(0)]; 
					double n1 = a2c[type.at(1)]; 
						 
					if(n0 > n1) 
					{
						ss.A = type.at(0); 
						ss.countA = n0; 
						ss.B = type.at(1); 
						ss.countB = n1; 
					}
					else 
					{
						ss.A = type.at(1); 
						ss.countA = n1; 
						ss.B = type.at(0); 
						ss.countB = n0; 
					}

//					if(ss.A == 'N') 
//					{
//						ss.A = ss.B; 
//						ss.B = 'N';
//						double tt = ss.countA; 
//						ss.countA = ss.countB; 
//						ss.countB = tt; 
//					}   //specially handling for monomorphic; 
				}   //translate a2c into mmSNP format; 
				a2c.clear(); 
				
				//grant: 
//				cout << ss.A << " " << ss.countA << " ||| "; 
//				cout << ss.B << " " << ss.countB<< " ||| "; 
//				cout << endl; 

				map<string, class mmSNP> :: iterator mmiter; 
				mmiter = rs2ss.find(rs); 
				if (mmiter != rs2ss.end())
				{
					fplog << "## WARNING: Two SNPs have the same rsnum " << mmiter->first << endl; 
					cout << "-bimbam: two snp have same rsnum " << mmiter->first << endl; 
				}
				
				rs2ss.insert(pair<string, class mmSNP> (rs, ss)); 
//				cout << rs << " " << rs2ss.size() << " "; 
//				cout << ss.A <<  ss.countA << " " << ss.B << ss.countB << " " << ss.countQ << endl; 
				
				if(readpos == 0) 
				{
					mapRs2pos.insert( pair<string, long> (rs, (long) rs2ss.size()));  
					mapRs2chr.insert( pair<string, int> (rs, 0)); 
				}
				//if no rs-pos file and only one genotype file, use snp order as pos.
				line.assign(getline(pbuf));
				int len = min(100, (int)line.size()); 
				for (int i = 0; i < len; i++)
				{
					if((int)line.at(i) == 0) 
						line.at(i) = 'N'; 
				}                       
			}

			infile.close();
			mfile2map[nf] = rs2ss; 
			
			ns = filerows; 
            ni = filecols;  
			
			if(filecols_diff > 0)
			{
				fplog << "## WARNING: " << filecols_diff << " SNPs has different # of ind's from the first line" << endl; 
				cout << "-bimbam: in " << filecols_diff << " snps, individual differ from the first line" << endl; 
//				safe_exit(); 
			}
			
			if(vFilePloid.at(nf) == 1)
			{
				vFileIndiv.push_back(2 * ni);
				nHap += 2 * ni; 
				for(int i = 0; i < 2 * ni; i++)
					vindploid.push_back(1); 
			}
			else 
			{
				vFileIndiv.push_back(ni);
				nDip += ni; 
				for(int i = 0; i < ni; i++)
					vindploid.push_back(2); 
			}
			fplog << "## BIMBAM: File " << nf << " has " << vFileIndiv.at(nf) << " ind's and " \
				<< rs2ss.size() << " SNPs" << endl; 
			cout << "-bimbam: file " << nf << " has " << vFileIndiv.at(nf) << " individual and " \
				<< rs2ss.size() << " snps" << endl; 
			rs2ss.clear(); 
		}   // SNP summary data;
		
//		{
//			map<string, class mmSNP> :: iterator iter; 
//			for (int nf = 0; nf < mfile2map.size(); nf++)
//			{
//				for (iter = mfile2map[nf].begin(); iter != mfile2map[nf].end(); ++iter)
//				{
//					cout << nf << "\t" <<  iter->first << "\t" << iter->second.A << iter->second.B;  
//					cout << "\t" << iter->second.countA << "\t" << iter->second.countB;   
//					cout << "\t" << iter->second.countQ << endl; 
//				}
//			}
//		}   //test if the mmSNP build ok; 

		
		////////////////////////////////////////////////////////////////////////////
		// now perform merger;   //snp by snp; 

		map<string, int> ::iterator m_iter; 
		map<string, vector<char> > mrs2coding;   //rs to reference allele coding for each file;  
		for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); ++m_iter)
		{
			string rs; 
			rs.assign(m_iter->first); 
			vector<char> coding; 
			class mmSNP ref; 
			int bingle = 0; 
			for (unsigned nf = 0; nf < vGin.size(); nf++)
			{
				map<string, class mmSNP> :: iterator iter; 
				iter = mfile2map[nf].find(rs); 
				if(iter == mfile2map[nf].end())
				{
					bingle = 1;  //if missing in one file, remove the snp; 
					break; 
				}
				if(ref.A == 'N' && ref.B == 'N')  // first one; 
				{
					ref.assign(iter->second); 
					ref.A = (iter->second).A; 
					ref.B = (iter->second).B; 
					ref.countA = (iter->second).countA; 
					ref.countB = (iter->second).countB; 
					ref.countQ += (iter->second).countQ; 
					coding.push_back(ref.A);
					continue; 
				}

				char cr = merge_mmSNP(ref, iter->second);
				if(cr == 'X') 
				{
					bingle = 1; 
					break; 
				}
				else 
					coding.push_back(cr); 
			}
			
//   		cout << bingle << " " << ref.A << ref.B << "\t" << ref.countA << "\t" << ref.countB;  
//			cout << "\t" << ref.countQ << endl; 
			
			if (bingle == 1) 
				m_iter->second = 0; //remove this SNP; 
			
			mrs2coding[rs] = coding;                  
			mapRs2mm[rs] = pair<char, char> (ref.A, ref.B); 
			double countAB = ref.countA + ref.countB; 
			double missrate = (double ) ref.countQ / (countAB + ref.countQ); 
			//cout << missrate << "\t" << ref.countA << "\t" << ref.countB << "\t" << ref.countQ << endl; 
			if((m_not_snp == 0) &&  missrate > m_exclude_miss)        //missing test; 
				mrs2yes[rs] = -1; 
			else if(countAB < 1e-6) 
				mapRs2maf[rs] = 0.0; 
			else 
				mapRs2maf[rs] = (double) (ref.countB) / countAB;  
			//reference allele is A, i.e. genotypes  = 2 - counts of reference alleles. 
			vector<char> ().swap(coding); 
		}
		
		int count_pos = 0;  //failed to have a position; 
		int count_am = 0;   //failed to match allele type; 
		int count_maf = 0;  //failed due to maf too extreme; 
		int count_miss = 0; 
		vector<pair<string, pair<long, int> > > vp;
		map<string, long> :: iterator pos_iter; 
		map<string, int> :: iterator chr_iter; 
		for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); m_iter++)
		{
			string rs(m_iter->first);
			if(m_iter->second == -1) 
			{
				count_miss ++; 
				continue; 
			}
			if(m_iter->second == 0) 
			{
				count_am ++; 
				continue; 
			}
			//to accommodate different allele coding; 
			if((m_not_snp == 0) && mapRs2maf[rs] <= m_exclude_maf)
			{
//				cout << mapRs2maf[rs] << endl; 
				count_maf ++;    
				m_iter->second = 0;  //mark this SNP to be excluded. 
				continue; 
			}
			if(m_exclude_nopos && mapRs2pos[rs] == 0)   //if snp no position; 
			{
				count_pos ++; 
				m_iter ->second = 0; //mark this SNP to be excluded. 
				continue; 
			}
			pair<string, pair<long, int> > tmp;
			tmp.first.assign(m_iter->first);
			pos_iter = mapRs2pos.find(m_iter->first);
			if(pos_iter != mapRs2pos.end())
				tmp.second.first = pos_iter->second;
			else 
				tmp.second.first = 0; 

			chr_iter = mapRs2chr.find(tmp.first); 
			if(chr_iter != mapRs2chr.end())
				tmp.second.second = chr_iter->second; 
			else
				tmp.second.second = 0; 
									   
			vp.push_back(tmp); 
		}
		if(count_miss > 0) 
		{
			fplog << "## BIMBAM: Exclude " << count_miss << " SNPs due to miss proportion > " << m_exclude_miss << endl;
			cout << "-bimbam: exclude " << count_miss << " snps due to miss proportion > " << m_exclude_miss << endl;
		}
		if(count_am > 0) 
		{
			fplog << "## BIMBAM: Exclude " << count_am << " SNPs due to failure to match betweeen files." << endl;
			cout << "-bimbam: exclude " << count_am << " snps due to failure to match betweeen files." << endl;
		}
		if(count_maf > 0) 
		{
			fplog << "## BIMBAM: Exclude " << count_maf << " SNPs due to maf < " << m_exclude_maf << endl;
			cout << "-bimbam: exclude " << count_maf << " snps due to maf < " << m_exclude_maf << endl;
		}
		if(count_pos > 0) 
		{
			fplog << "## BIMBAM: Exclude " << count_pos << " SNPs due to no position information" << endl;
			cout << "-bimbam: exclude " << count_pos << " snps due to no position information" << endl;
		}
		stable_sort(vp.begin(), vp.end(), rscomp); 
		//sort the rs pos vector in order of chr && pos;
		map<string, int> mrs2index; //the column of each snp; 
		vsRsnum.clear(); 
		for(int i = 0; i < (int) vp.size(); i++)
		{
			vsRsnum.push_back(vp.at(i).first);
			mrs2index[vp.at(i).first] = i; 
		}
		vector<pair<string, pair<long, int> > >().swap(vp); 
		mrs2yes.clear(); 
		nLoci = (int) vsRsnum.size(); 
		nIndiv = nDip + nHap; 
		nCohort = nDip; 
		if(nLoci == 0) 
		{
			cout << "-bimbam: no valid snp" <<  endl; 
			safe_exit(); 
		}
        //now vsRsnum hold ID of all valid SNPs in order; 
		
		vPARCrs.clear(); 
		vPARCpos.clear(); 
		for (unsigned i = 0; i < vsRsnum.size(); i++)
		{
			m_iter = mCohortRsQ.find(vsRsnum.at(i));
			if(m_iter != mCohortRsQ.end())  
			{
				vPARCrs.push_back(m_iter->first); 
				vPARCpos.push_back(i); 
			}
		}

		mCohortRsQ.clear(); 
		////////////////////////////////////////////////////////////////////////////////////
		
		fstream outfile; 
		string sfn("output/");
		sfn.append(fnOutput);
		sfn.append(".snpdata.txt");
		outfile.open(sfn.c_str(), ios::out);
		if(outfile.is_open()) 
		{
			outfile << "## af is the allele freq for A, the minor allele" << endl;  
			outfile << "rs\t A\t B\t af \t chr \t pos" << endl;  
			for (int m = 0; m < nLoci; m++)
			{
				string rs = vsRsnum.at(m); 

				char buf[100]; 
				sprintf(buf, "%-s\t", rs.c_str()); 
				outfile << buf; 
				
				outfile << " " << mapRs2mm[rs].second << " " << mapRs2mm[rs].first << "\t";
				sprintf(buf, "%.3f\t", mapRs2maf[rs]); 
				
				outfile << buf;
				if(mapRs2chr[rs] < 0)
					outfile << " NA \t"; 
				else
					outfile << mapRs2chr[rs] << "\t"; 
				outfile << mapRs2pos[rs] << endl; 
			}
			outfile.close(); 
		}
		else 
			cout << "-bimbam: skip writing snpdata" << endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////
        
//		cout << " before allocate pIndiv " << endl; 
		pIndiv = new Individual * [nIndiv]; 
	
		for (int i = 0; i < nIndiv; i++)
		{
			if (vindploid.at(i) == 1) {
				pIndiv[i] = new HapInd; 
				if(mode == 1) 
					pIndiv[i]->allocate_snp_mgt(nLoci); 
				else //if (mode == 2) 
					pIndiv[i]->allocate_snp_dstr(nLoci); 
				pIndiv[i]->SetisPanel(0);
			} else {
				pIndiv[i] = new DipInd;
				if(mode == 1) 
					pIndiv[i]->allocate_snp_mgt(nLoci); 
				else //if (mode == 2) 
					pIndiv[i]->allocate_snp_dstr(nLoci); 
				pIndiv[i]->SetisPanel(0);
			}
		}
		
////////////////////////////////////////////////////////////////////////////////////////////////////
		// go over genotype files again; 
		// for each diploid file, scan each SNP, fill in the genotype table. 
		int start_ni = 0; 
		for (int nf = 0; nf < (int)vGin.size(); nf++)
		{       
			streambuf * pbuf;
			ifstream infile; 
			char delimit[] = ";, \t";
			infile.open(vGin.at(nf).c_str(), ios::in);
			if(!infile.is_open()) 
			{
				cout << "-bimbam: cannot open genotype file: " << vGin.at(nf) << endl; 
				safe_exit(); 
			} 
			pbuf = infile.rdbuf();
			 
			string line; 
			line.assign(getline(pbuf)); 
			//possible individual ids
				
			while(line.size() > 0)
			{   
				int len = min(100, (int)line.size()); 
				for (int i = 0; i < len; i++)
				{
					if((int)line.at(i) == 0) 
						line.at(i) = 'N'; 
				}                       
				char * res = strtok((char*)line.c_str(), delimit); 
				if(res == NULL) break; 
				string rs; 
				rs.assign(res);  
				if(end_pos > 0) 
				{
					map<string, long> :: iterator iter; 
					iter = mapRs2pos.find(rs); 
					if(iter == mapRs2pos.end() || iter->second < 0) 
					{
						line.assign(getline(pbuf));
						continue; 
					}
				}    
				int m = -1; 
				map<string, int> :: iterator s2i; 
				s2i = mrs2index.find(rs); 
				if(s2i == mrs2index.end())
				{
					line.assign(getline(pbuf));
					continue; 
				}
				else 
				 	m = s2i->second;  
				// rs number; 
				char aA, aB; 
				res = strtok(NULL, delimit); 
				if(res != NULL) aA = res[0]; 
				else aA = 'N'; 
				res = strtok(NULL, delimit); 
				if(res != NULL) aB = res[0]; 
				else aB = 'N'; 
				
				
				char ref = mrs2coding[rs].at(nf); //reference allele; 
				//grant:
//				cout << ref << endl; 
				int ni = start_ni; 
				for (int i = 0; i < vFileIndiv.at(nf); i++)
				{    
					double g0 = -1.0; 
					double g1 = -1.0; 
					char * res = strtok(NULL, delimit); 
					if(res == NULL) break; 
					g0 = atof(res); 
					if(mode == 2) 
					{
						res = strtok(NULL, delimit); 
						if(res == NULL) 
						{
							cout << "-bimbam: odd number of numerical entries" << endl;
							safe_exit(); 
						}
						g1 = atof(res); 
					}
				
					if(mode == 1) 
					{
						if(ref == aA) 
							pIndiv[ni]->set_snp_mgt(m, g0); 
						else 
							pIndiv[ni]->set_snp_mgt(m, 2.0-g0); 
					}
					else if (mode == 2)
					{
						if(ref == aA)                      //if reference allele match flip; else no flip;
							pIndiv[ni]->set_snp_dstr(m, g0, g1); 
						else 
							pIndiv[ni]->set_snp_dstr(m, 1.0 - g0 - g1, g1); 
					}
					ni++; 
				}
				line.assign(getline(pbuf));
			}
			infile.close();
			cout << "-bimbam: read file " << nf << " again " << endl; 
			start_ni += vFileIndiv.at(nf); 
		}

		nPanel = 0; 
		for (int i = 0; i < nIndiv; i++)
			nPanel += pIndiv[i]->GetisPanel(); 
		nCohort = nIndiv - nPanel; 

		//////////////////////////////////////////////////////////////////////////////////
		
		fplog << "## number genotype files = " << vGin.size() << endl; 
		fplog << "## number phenotype files = " << vPin.size() << endl; 
		fplog << "## number of diploid = " << nDip << endl; 
		fplog << "## number of haploid = " << nHap << endl; 
	   	fplog << "## number of panel individuals = " << nPanel << endl; 
		fplog << "## number of cohort = " << nCohort << endl; 
		fplog << "## number of snp = " << nLoci  << endl;
		fplog << "## number cohort only snp = " << vPARCrs.size() << endl; 
	}

//////////////////////////////////////////////////////////////////////////////////////////
	
	return 1; 
}  

int ModelnData::read_bimbam_phenotype(int mcmc)
{
	if(vPin.size() == 0 || vPin.size() != vGin.size()) 
		return 0; 
	vector<vector<real> > trans_phval; 
#if defined (MPI_ENABLED)
	real * tmp_phval; 
	if(procID == 0)
#endif
	{
		fstream infile; 
		streambuf * pbuf; 
		string sfn; 
		char delimit[] = ",; :\t";
			   
		int cur_ni = 0; 
		for (int nf = 0; nf < (int) vPin.size(); nf++)
		{   
			int ni = vFileIndiv.at(nf);
			sfn.assign(vPin.at(nf)); 
			
			if(sfn.compare("0") == 0)    //panel, put some dummy values NA; 
			{
				for (int i = 0; i < ni; i++)
				{
					vector<real> ind_phval; 
					for (int np = 0; np < nPH; np++)
						ind_phval.push_back(NA); 
					trans_phval.push_back(ind_phval);
				}
				cur_ni += ni; 
				continue; 
			}
			
			else if(sfn.compare("1") == 0)   //case; 
			{
				for (int i = 0; i < ni; i++)
				{
					vector<real> ind_phval; 
					for (int np = 0; np < nPH; np++)
						ind_phval.push_back(1); 
					trans_phval.push_back(ind_phval);
				}
				cur_ni += ni; 
				continue; 
			}
			
			else if(sfn.compare("z") == 0)  //control; 
			{
				for (int i = 0; i < ni; i++)
				{
					vector<real> ind_phval; 
					for (int np = 0; np < nPH; np++)
						ind_phval.push_back(0); 
					trans_phval.push_back(ind_phval);
				}
				cur_ni += ni; 
				continue; 
			}
			
			infile.open(sfn.c_str(), ios::in);
			if(!infile.is_open())
			{
				fplog << "ERROR: BIMBAM cannot open phenotype file " << nf << endl; 
				cout << "-bimbam: cannot open phenotype file " << nf << endl; 
				safe_exit(); 
			}
			pbuf = infile.rdbuf();
			string line;
			
//			cout << ni << endl; 
//			for (int i = 0; i < vFileIndiv.size(); i++)
//				cout << vFileIndiv.at(i) << endl; 
//			cout << endl; 
			for (int i = 0; i < ni; i++)
			{
				vector<real> ind_phval; 
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int np = 0; np < nPH; np++)
				{
					if(res == NULL) break; 
					string sv(res);
					if(sv.compare("NA") == 0) ind_phval.push_back(NA);
					else ind_phval.push_back(atof(sv.data()));
					res = strtok(NULL, delimit); 
				}

				if ((int) ind_phval.size() != nPH)
				{
					nPH = ind_phval.size(); 
					fplog << "WARNING: File " << nf << " phenotypes columns smaller than specified." << endl;
					cout << "-bimbam: file " << nf << " phenotypes columns smaller than specified." << endl;
				}
					
				trans_phval.push_back(ind_phval);
				ind_phval.clear();
			}
			line.assign(getline(pbuf));
			if(v_reduced_num_indiv.size() == 0 && line.size() > 0) 
			{
				fplog << "WARNING: Phenotypes in file " << nf << " larger than number of indiv." << endl;
				cout << "-bimbam: phenotypes in file " << nf << " larger than number of individuals" << endl;
			}             
			if ((int)trans_phval.size() - cur_ni < ni) 
			{
				fplog << "QUIT: Phenotypes in file " << nf << " smaller than number of indiv." << endl;
				cout << "-bimbam: phenotypes in file " << nf << " smaller than number of individuals" << endl;
                safe_exit(); 
			}
			infile.close(); 
			cur_ni += ni; 
		}
	}
	
#if defined (MPI_ENABLED)
   	tmp_phval = new double[nPH * nIndiv];
	if (procID == 0) 
	{
		for (int np = 0; np < nPH; np++)
			for (int i = 0; i < nIndiv; i++)
				tmp_phval[i * nPH + np] = trans_phval.at(i).at(np); 
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Bcast(tmp_phval, nPH * nIndiv, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

	for(int np = 0; np < nPH; np++)
	{
		vector<int> vindex; 
		vector<real> vt;
		int ni = 0; 
		for(int i = 0; i < nIndiv; i++)
		{
			if (fabs(tmp_phval[i*nPH+np] - NA) < 0.1)
			{
				ni++; 
			} 
			else 
			{
				vindex.push_back(ni++); 
				vt.push_back(tmp_phval[i * nPH + np]);
			}
		}
		vv_phval.push_back(vt);
		vv_phdex.push_back(vindex);
		vindex.clear(); 
		vt.clear();
	}
	
#else 
	
	for(int np = 0; np < nPH; np++)
	{
		vector<int> vindex; 
		vector<real> vt; 
		int ni = 0; 
		for(int i = 0; i < nIndiv; i++)
		{
			if (pIndiv[i]->GetisPanel()) continue;  
			if (fabs(trans_phval.at(i).at(np) - NA) < 0.1)
			{
				ni++; 
				vt.push_back(trans_phval.at(i).at(np));
			} 
			else 
			{
				vindex.push_back(ni++); 
				vt.push_back(trans_phval.at(i).at(np));
			}
		}
		vv_phval.push_back(vt);
		vv_phdex.push_back(vindex);
		vindex.clear(); 
		vt.clear();
	}   //vv_phdex contains "index" of indiv that has valid phval. 
#endif	

//	for (int i = 0; i < vv_phval.at(0).size(); i ++)
//		cout << vv_phval.at(0).at(i) << " "; 
//	cout << endl; 
	
#if defined (MPI_ENABLED)
	delete[] tmp_phval; 
	if (procID == 0)
#endif  
	{
		vector<vector<real> >().swap(trans_phval); 
		cout << "-bimbam: number of phenotypes = " << nPH << endl; 
		fplog << "###: number of phenotypes = " << nPH << endl; 
	}
	return 1; 
}

void ModelnData::assign_population_label(void)
{
    if (hasPopStructure.compare("0") == 0)  //no pop structure <==> nSubPop = 1
	{   
		nSubPop = 1;
		vnSubHap.push_back(0);
		for (int i = 0; i < nIndiv; i++)
		{
			pIndiv[i]->SetpopLabel(0);
			vnSubHap.at(0) += pIndiv[i]->Getploid();   //vnSubHap = how many haploids in a subpop. 
		}
	} 
	else 
	{
		vector<char> vIndex;
		fstream infile(hasPopStructure.c_str(), ios::in); 
		if(!infile.good()) 
		{
			fplog << "## ERROR: BIMBAM can't open population label file." << endl;
			cout << "-bimbam: cannot open population label file" << endl;
			safe_exit();  
		} 
		else 
		{
			for (int i = 0; i < nIndiv; i++)
			{
				char tmp; 
				infile >> tmp;
				vIndex.push_back(tmp); 
			}
		}
		
		vector<char> vtype; 
		for (int i = 0; i < nIndiv; i++)
		{
			bool bingle = 0;
			for (unsigned int j = 0; j < vtype.size(); j++)
			{
				if (vtype.at(j) == vIndex.at(i))
				{
					bingle = 1; 
					break;
				}
			}
			if (!bingle)
				vtype.push_back(vIndex.at(i));
		}   //figure how many types in the poplabel. 
		
		nSubPop = vtype.size();
	 	vnSubHap.insert(vnSubHap.end(), nSubPop, 0); 
	    for (int i = 0; i < nIndiv; i++)
		{
			for (int j = 0; j < nSubPop; j++)
			{
				if (vtype.at(j) == vIndex.at(i))
				{
					pIndiv[i]->SetpopLabel(j);
					vnSubHap.at(j) += pIndiv[i]->Getploid();
					break;
				}
			}
		}
		vector<char>().swap(vIndex); 
		//vIndex.clear();
	}
#if defined (MPI_ENABLED)
	if(procID == 0) 
#endif 
		fplog << "## nSubPop = " << nSubPop << endl; 
}

#if defined (IMPUTATION)
void ModelnData::MaskRand(void)
{
	nMasked = (int) (parcLoci * percentMasked); 
	int * pMaskedSNP = new int[nMasked];
	cout << "-bimbam: mask genotype at random each individual masked " << nMasked << endl; 
	
	for (int i = 0; i < nIndiv; i++) 
	{
		if(pIndiv[i]->GetisPanel()) continue; 
		DistinctIntArray(1, parcLoci-1, nMasked, pMaskedSNP); 
		for (int m = 0; m < nMasked; m++)
			pMaskedSNP[m] = vPARCpos.at(pMaskedSNP[m]); 
		pIndiv[i]->MaskSNPs(nLoci, nMasked, pMaskedSNP);
	}
	delete[] pMaskedSNP; 
}

void ModelnData::MaskSNP(void)
{
	nMasked = 0; 
	int * pMaskedSNP = NULL;
	{
		int step = 0; 
		int start = m_num; 
		if(percentMasked > 1) 
			step = (int) percentMasked;
		else 
			step = (int) (1.0 / percentMasked);
		
		fplog << "## step:start:parcLoci " << step << " : " << start << " : " << vPARCpos.size() << endl;

		int temp_count = 0; 
		map<string, int> ::iterator iter; 
		vector<int> vPos;
		for (unsigned i = start; i < vPARCpos.size(); i += step)
		{
			vMaskedSNPrs.push_back(vPARCrs.at(i));
			iter = mPanelRsQ.find(vPARCrs.at(i)); 
			if(iter == mPanelRsQ.end())
				temp_count++; 
			vPos.push_back(vPARCpos.at(i)); 
		}
		cout << "-bimbam: number of masked snp that is not in panel = " << temp_count << endl; 
		fplog << "## number of masked SNP that is not in panel = " << temp_count << endl; 
		
		nMasked = (int) vPos.size(); 
		pMaskedSNP = new int[nMasked];
		fplog << "## masked SNP : " << vPos.size() << " : " ; 
		for (int i = 0; i < (int) vPos.size(); i++)
		{
			pMaskedSNP[i] = vPos.at(i); 
			fplog << pMaskedSNP[i] << " "; 
		}
		fplog << endl; 
		vector<int>().swap(vPos); 
		//vPos.clear(); 
   	}

	for (int i = 0; i < nIndiv; i++) 
	{
		if(pIndiv[i]->GetisPanel()) continue; 
		pIndiv[i]->MaskSNPs(nLoci, nMasked, pMaskedSNP);
	}
}

void ModelnData::ImputeMasked(void)
{
	double * pmaf = new double[nMasked];
	double * sum_maf = new double[nMasked]; 
	for (int i = 0; i < nMasked; i++)
		pmaf[i] = sum_maf[i] = 0.0; 
	int * err_count_all = new int[nMasked];
	int * na_count_all = new int[nMasked];
	int * err_count = new int[nMasked]; 
	int * na_count = new int[nMasked]; 
	for (int m = 0; m < nMasked; m++)
	{
		err_count_all[m] = na_count_all[m] = 0; 
		err_count[m] = na_count[m] = 0; 
	}
	
	for (int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->GetisPanel()) continue; 
		pIndiv[i]->ImputeMaskedCounting(nMasked, err_count_all, na_count_all);
		real * maf = pIndiv[i]->Getmaf();  
		if(maf != NULL)
		{
			for (int m = 0; m < nMasked; m++)
				sum_maf[m] += maf[m]; 
		}
	}	                       
	{
		real err_total = 0; 
		real na_total = 0; 
		int temp_count = 0; 
		map<string, int> ::iterator iter; 
		for (int m = 0; m < nMasked; m++)
		{
			iter = mPanelRsQ.find(vMaskedSNPrs.at(m)); 
			if(iter == mPanelRsQ.end())
			{
				temp_count++; 
				continue; 
			}
			err_total += err_count_all[m];
			na_total += na_count_all[m]; 
		}
		real  total =  nCohort * (nMasked - temp_count);
		fplog << "## masked per indiv. = " << nMasked - temp_count << endl; 
		fplog << "## total masked = " << total << endl;
		fplog << "## countQ = " <<  na_total << endl; 
		fplog << "## total missed = " << err_total << endl;
		fplog << "## the pecentage of miss = " << err_total /(total - na_total) << endl;
		fplog << "## rs \t hapmap \t err \t qq \t maf" << endl; 
		
		int step = 0; 
		int start = 0; 
		if(percentMasked > 1) {
			step = (int) percentMasked;
			start = m_num;
		} else {
			step = (int) (1.0 / percentMasked);
			start = gsl_rng_uniform_int(gsl_r, step); 
		}

		for (int i = 0; i < nMasked; i++) 
		{
		   	int pos = start + i * step; 
			fplog << vMaskedSNPrs.at(i) << "\t";  
			iter = mPanelRsQ.find(vPARCrs.at(pos)); 
			fplog << (iter != mPanelRsQ.end()) << "\t"; 
			fplog << err_count_all[i] << "\t" << na_count_all[i] << "\t";
			double tmp = sum_maf[i] / (2.0 * (nCohort - na_count_all[i]));  
			fplog << tmp << "  " << endl; 
		}
		fplog << endl;                                                
	}
	delete[] sum_maf; 
	delete[] pmaf; 
}
#endif

void ModelnData::EM(int warm)
{
	stringstream lss; 
	
	sumiTop = Allocate2DMatrix(nLoci, nK);
	sumiBot = Allocate2DMatrix(nLoci, nK);
	sumiJm  = Allocate2DMatrix(nSubPop, nLoci);
	sumiJmk = Allocate3DMatrix(nSubPop, nLoci, nK);

	if(pMP == NULL) InitModelParam(); 
#if defined (MPI_ENABLED)
    int size = (nLoci * nK * (2 + nSubPop) + nLoci * nSubPop);
//	int size = nSubPop * nLoci * nK + nSubPop * nLoci  + 2 * nLoci * nK;  
    pParamCollect = new real[size];
    real * pBuffer = new real[size];
	memset(pParamCollect, 0, sizeof(real) * size); 
	memset(pBuffer, 0, sizeof(real) * size); 
#endif
	int total_steps; 
	if(warm) total_steps = nWarmSteps;
	else total_steps = nMaxSteps; 
	int hap_in_panel = 0;  
	for (int i = 0; i < nIndiv; i++)
		hap_in_panel += pIndiv[i]->GetisPanel() * pIndiv[i]->Getploid(); 
//	cout << procID << " has hap_in panel = " << hap_in_panel << endl; 

	int prg_cur = 0; 
	int prg_all = nEMRuns * total_steps;
	if(warm) prg_all *= nPanel; 
	else prg_all *= (nPanel + nCohort); 
	char str[100]; 
	if(warm) sprintf(str, "-bimbam: em on panel: "); 
	else sprintf(str, "-bimbam: em on cohort: "); 
	int prg_mod = (int) (prg_all / 100.0); 
	if(prg_mod == 0) prg_mod = 1; 
	//need more for mpi; 

	for(int runs = 0; runs < nEMRuns; runs++)
	{ 
		for(int steps = 0; steps < total_steps; steps++)
		{
			//initialize sumi's before each steps. 
			total_likelihood = 0.0;
			logLikelihood = 0.0;
			for (int m = 0; m < nLoci; m++) 
				for (int k = 0; k < nK; k++)
				{
					sumiTop[m][k] = 0.0;
					sumiBot[m][k] = 0.0;
    			}    	                                     
			for (int s = 0; s < nSubPop; s++)
				for (int m = 0; m < nLoci; m++)
				{
					sumiJm[s][m] = 0.0;
					for (int k = 0; k < nK; k++)
						sumiJmk[s][m][k] = 0.0;
				}
		
			int count = 0; 
#if defined (MPI_ENABLED)
			for (unsigned int wp = 0; wp < pvLoadDip[procID].size(); wp++)
#else
			for (int i = 0; i < nIndiv; i++)
#endif         
			{
#if defined (MPI_ENABLED)
				int i = pvLoadDip[procID].at(wp);
#endif			
				if (warm == 1 && pIndiv[i]->GetisPanel() == 0) continue;
//				if (warm == 0 && pIndiv[i]->GetisPanel() == 1) continue; 
				count ++; 
				pIndiv[i]->CalcAll(nLoci, nK, pMP+runs, 1);
				logLikelihood += pIndiv[i]->GetLikelihood();
				
				real ** top = pIndiv[i]->Gettopptr();
				for (int m = 0; m < nLoci; m++)
					for (int k = 0; k < nK; k++)
						sumiTop[m][k] += top[m][k];
				
				real ** bot = pIndiv[i]->Getbotptr();
				for (int m = 0; m < nLoci; m++)
					for (int k = 0; k < nK; k++)
						sumiBot[m][k] += bot[m][k];
				
				real ** jmk = pIndiv[i]->GetexpectJmkptr();
				int s = pIndiv[i]->GetpopLabel(); 
				for (int m = 0; m < nLoci; m++)
					for (int k = 0; k < nK; k++)
					{                                                        
						sumiJmk[s][m][k] += jmk[m][k]; 
						sumiJm[s][m] += jmk[m][k]; 
					}
				
				pIndiv[i]->FreeMemAfterEM();
#if defined (MPI_ENABLED)
				if(procID == 0)
#endif 
				{
				    prg_cur++;             
					if(prg_cur % prg_mod == 0) 
						print_progress_bar(0, str, prg_cur, prg_all);
				}
			}

#if defined (MPI_ENABLED)           
            MPI_Barrier(MPI_COMM_WORLD); //wait for all processors finish the job;  
            if(sizeof(real) == sizeof(double))
                MPI_Reduce(&logLikelihood, &total_likelihood, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else
                MPI_Reduce(&logLikelihood, &total_likelihood, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); //wait for all processors finish the job;  
            // do the sumi on each processor; 

            //tar everything into an arrayfor communication; 
            int pos = 0;
            for (int m = 0; m < nLoci; m++)
                for (int k = 0; k < nK; k++)
				{
					pBuffer[pos] = sumiTop[m][k];
					pos++;
				}
            for (int m = 0; m < nLoci; m++)
                for (int k = 0; k < nK; k++)
				{
					pBuffer[pos] = sumiBot[m][k];
					pos++;
				}

            for (int s = 0; s < nSubPop; s++)
                for (int m = 0; m < nLoci; m++)
                    for (int k = 0; k < nK; k++)
					{
						pBuffer[pos] = sumiJmk[s][m][k];
                        pos++;
					}
            for (int s = 0; s < nSubPop; s++)
                for (int m = 0; m < nLoci; m++)
				{
                    pBuffer[pos] = sumiJm[s][m];
                 	pos++; 
				}
            for (int i = 0; i < size; i++)
                pParamCollect[i] = 0.0;

            MPI_Barrier(MPI_COMM_WORLD);
            if(sizeof(real) == sizeof(double))
                MPI_Allreduce(pBuffer, pParamCollect, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            else
                MPI_Allreduce(pBuffer, pParamCollect, size, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            pos = 0;
            for (int m = 0; m < nLoci; m++)
                for (int k = 0; k < nK; k++)
				{
                    sumiTop[m][k] = pParamCollect[pos];
					pos++;
				}

            for (int m = 0; m < nLoci; m++)
                for (int k = 0; k < nK; k++)
				{
                    sumiBot[m][k] = pParamCollect[pos];
                 	pos++;
				}
            for (int s = 0; s < nSubPop; s++)
                for (int m = 0; m < nLoci; m++)
                    for (int k = 0; k < nK; k++)
					{
                        sumiJmk[s][m][k] = pParamCollect[pos];
                        pos++;
					}
            for (int s = 0; s < nSubPop; s++)
                for (int m = 0; m < nLoci; m++)
				{ 
					sumiJm[s][m] = pParamCollect[pos];
					pos++;
				}
#endif
	        // Update Param;
		    for(int m = 0; m < nLoci; m++)
		    {    
				for (int s = 0; s < nSubPop; s++)
				{
					real tr;
					if(warm) tr = sumiJm[s][m] / hap_in_panel; 
					else tr = sumiJm[s][m] / vnSubHap.at(s); 
					if(isnan(tr)) 
					{
//						fplog << "### tr is nan " << sumiJm[s][m] << "\t" << m<< endl; 
                        continue; 
					}
					if (tr > 0.99999) {
//						fplog << "tr > 1.0 " << m << " " << tr << endl; 
						tr = 0.99999;                   
					}
					else if (tr < 1e-5) {
//						fplog << "tr < 0.0 " << m << " " << tr << endl; 
						tr = 0.00001;
					}
					
					pMP[runs].Setr(s, m, tr);
				}   // r 
			   
				for (int s = 0; s < nSubPop; s++)
				{
					if(!(isnan(sumiJm[s][m]) || sumiJm[s][m] < 1e-100 )) 
					{
						double maxj = -1; 
						for (int k = 0; k < nK; k++)
						{
                            if(sumiJmk[s][m][k] > maxj) 
								maxj = sumiJmk[s][m][k];  
						}

						double tj = 0; 
						for (int k = 0; k < nK; k++)
						{
							double ratio = sumiJmk[s][m][k] / maxj; 
                            if(ratio < 1e-6) 
								sumiJmk[s][m][k] = 1e-6; 
							else 
								sumiJmk[s][m][k] = ratio; 
							tj += sumiJmk[s][m][k]; 
						}

						for(int k = 0; k < nK; k++)
						{
							double ta = sumiJmk[s][m][k] / tj;
							pMP[runs].Setalpha(s, m, k, ta);
						}	
					}
				}   //alpha

				for (int k = 0; k < nK; k++)
				{
					real tt = sumiTop[m][k] / sumiBot[m][k];
					if(isnan(tt)) 
					{
//				    	fplog  << "theta is nan " << sumiTop[m][k] << " " << sumiBot[m][k] << endl; 
						continue;
					}
					if (tt > 0.999) tt = 0.999;
					else if (tt < 0.001) tt = 0.001;
					
					pMP[runs].Settheta(m, k, tt);  
				}   //theta
			}
			
			//in a way, this looks like maximization step. 
			//we do so seperately for each individual, then pool them togehter later;
			//in a way, this part looks like expectation. 
			
#if !defined (MPI_ENABLED)
			total_likelihood = logLikelihood; 
#else
			if(procID == 0)
#endif      
			{
				if(warm) fplog << "## warm up em runs = ";
				else fplog << "## em runs = "; 
				fplog << runs << " of " << nEMRuns << "\t steps = " << steps << " of " << total_steps;
				if(warm) fplog << "\t log-likelihood = " << total_likelihood << endl;
				else fplog << "\t log-likelihood = " << total_likelihood << endl;
			}
		} // end total_steps; 
	} //end em loop 
  
#if defined (MPI_ENABLED)
	delete[] pBuffer; 
	delete[] pParamCollect; 
	if(procID== 0) 
#endif 
	    print_progress_bar(1, str, 1, 1);
	Free2DMatrix(sumiTop);
	Free2DMatrix(sumiBot);
	Free2DMatrix(sumiJm);
	Free3DMatrix(sumiJmk);   
	fplog << "##EM -c " << nK << endl; 
}

void ModelnData::SNP_Density()
{
	int prg_cur = 0; 
	int prg_all = nEMRuns * nCohort;
	char str[100] = "-bimbam: calculate snp posterior: "; 
	int prg_mod = (int) (prg_all / 100.0); 
	if(prg_mod == 0) prg_mod = 1; 
// calculate and normalize genotype distributions; 
#if !defined (MPI_ENABLED)		
	for (int runs = 0; runs < nEMRuns; runs++)
	{
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
		   	pIndiv[i]->CalcAll(nLoci, nK, pMP+runs, 0);
			pIndiv[i]->calc_snp_dstr(nLoci, nK, pMP[runs].Gettheta());
			pIndiv[i]->FreeMemAfterEM();   
			prg_cur ++;
			if(prg_cur % prg_mod == 0) 
				print_progress_bar(0, str, prg_cur, prg_all);
		}
	}

	for (int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->GetisPanel()) continue; 
		pIndiv[i]->norm_snp_dstr(nLoci, nEMRuns); 
	}
	print_progress_bar(1, str, 1, 1);
#else
	for (int runs = 0; runs < nEMRuns; runs++)
	{ 
 		for (unsigned int wp = 0; wp < pvLoadDip[procID].size(); wp++)
        {
            int i = pvLoadDip[procID].at(wp);
			if(pIndiv[i]->GetisPanel()) continue; 
			pIndiv[i]->CalcAll(nLoci, nK, pMP+runs, 0);
			pIndiv[i]->calc_snp_dstr(nLoci, nK, pMP[runs].Gettheta());
			pIndiv[i]->FreeMemAfterEM();   
        }   
	}
	for (unsigned int wp = 0; wp < pvLoadDip[procID].size(); wp++)
	{
		int i = pvLoadDip[procID].at(wp);
		if(pIndiv[i]->GetisPanel()) continue; 
		pIndiv[i]->norm_snp_dstr(nLoci, nEMRuns); 
	}
#endif
}  

bool bfcomp(pair<real, int> a, pair<real, int> b)
{
	return (a.first > b.first); 
}

void ModelnData::single_snp_analysis(int nShuffle)
{
	fplog << "### nPH = " << nPH << endl; 
	real * bf1 = (real *) Allocate1D(sizeof(real), nPH * nLoci); 
	if (bf1 == NULL) 
	{
		cout << "-bimbam: the allocation for " << sizeof(real) * nPH * nLoci << " bytes failed."; 
		fplog << "-bimbam: the allocation for " << sizeof(real) * nPH * nLoci << " bytes failed."; 
		safe_exit(); 
	}
	real * bf2 = NULL; 
	
	time_t sec_beg, sec_end; 
	sec_beg = time (NULL);
	if(nImpute == 0)
		single_snp_exact(1, bf1);   // bf exact; 
	else if(nImpute == 21)
		single_snp_exact(2, bf1);   // f-stat exact; 
	else if (nImpute == 22)
		single_snp_exact(3, bf1);   // lrt exact; 
	
	else if (nImpute < 100)
		single_snp_mean_genotype(nImpute, bf1);     // 1, 2, 3, for mean gentype; 11, 12, 13 for best guess genotype; 
	else
	{ 
		bf2 = (real *) Allocate1D(sizeof(real), nPH * nLoci); 
		if (bf2 == NULL) 
		{
			cout << "-bimbam: the allocation for " << sizeof(real) * nPH * nLoci << " bytes failed."; 
			fplog << "-bimbam: the allocation for " << sizeof(real) * nPH * nLoci << " bytes failed."; 
			safe_exit(); 
		}
		single_snp_importance_sampling(bf1, bf2);
	}
	sec_end = time(NULL); 
	fplog << "## calc stat seconds used =" <<  sec_end - sec_beg << endl; 
	
	m_sumbf.clear();               
	m_sumbf.push_back(0); 
	for (int m = 0; m < nLoci; m++)
		m_sumbf.at(0) += bf1[m]; 
	m_sumbf.at(0) /= nLoci;
	if(nImpute > 0) 
		m_sumbf.at(0) /= nImpute; 
	
	if(m_sortQ == 1) 
	{
		vector< pair<real, int> > bfsort;
		for (int m = 0; m < nLoci; m++)
		{
			pair<real, int> tp; 
			tp.first = bf1[m];
			tp.second = m;
			bfsort.push_back(tp); 
		}
		stable_sort(bfsort.begin(), bfsort.end(), bfcomp); 
		
		for(int m = 0; m < nLoci; m++)
			mrank2loc[m] = bfsort.at(m).second; 
		//  m is the rank of bfsort.at(m)
		vector< pair<real, int> >().swap(bfsort);
	}
	else  
	{
		for(int m = 0; m < nLoci; m++)
			mrank2loc[m] = m; 
	}    
		
	if (nPH == 1 && nShuffle >= 100) single_snp_pvalue(nShuffle); 
	//calculate p-value only doable for single phenotype;  
	
	{
		string sfn("output/");
		sfn.append(fnOutput);
		sfn.append(".single.txt");
		fstream outfile;
		outfile.open(sfn.c_str(), ios::out);
		if(!outfile.is_open()) 
		{
			cout << "-bimbam: cannot open file to write:" << sfn << endl;
			return;
		}
		real total = (real) nImpute; 
		if(nImpute < 100) total = 1; 
		
		if(nPH == 1)
		{
			outfile << "## bf=log10_BF\tse=log10_std_err\tpv=p-value\t" << endl; 
			outfile << "rsnum\t\tpos\t chr \tbf\tse\trank\tpv\t mu\t a "; 
			if(m_df == 2) outfile << " \t d" << endl; 
			else outfile << endl; 

			map<string, vector<real> > :: iterator miter; 
			map<string, int> :: iterator chr_iter; 
			for (int r = 0; r < nLoci; r++) 
			{
				char buf[1000];
				int loc = mrank2loc[r]; 
				double t1, t2; 
				t1 = bf1[loc];
				if(bf2 == NULL) t2 = 0.0; 
				else t2 = bf2[loc];

				if(total > 1) 
				{
					double tmax = max(t1 * 2.0, t2); 
					double diff1 = t1 * 2.0 - tmax; 
					double diff2 = t2 - tmax; 
				 //   if(diff2 < diff1) cout << diff2 << "\t " << diff1 << endl; 
					t2 = 0.5 * (tmax + log(exp(diff2) - exp(diff1)) - log(total)); 
				}
				t1 /= log(10.0); 
				t2 /= log(10.0); 
				
				string rs;
				rs.assign(vsRsnum.at(loc)); 
				chr_iter = mapRs2chr.find(rs); 
				int chr = 0; 
				if(chr_iter != mapRs2chr.end())
					chr = chr_iter->second; 
				
				sprintf(buf, "%-15s\t%ld\t", rs.data(), mapRs2pos[rs]); 
				outfile << buf; 
				if(chr < 0) outfile << " NA \t"; 
				else outfile << chr << "\t"; 
				if (total >= 2) 
					sprintf(buf, "%+5.3f\t%+5.3f\t %d", t1, t2, r+1);
				else 
					sprintf(buf, "%+5.3f\t NA \t %d",  t1, r+1);
				outfile << buf; 
				if(nShuffle >= 100) 
					sprintf(buf, "\t\t %.2e", m_pv.at(loc)); 
				else 
					sprintf(buf, "\t\t NA"); 
				outfile << buf; 
				
				miter = mapRs2beta.find(rs); 
				if (miter == mapRs2beta.end())
					outfile << "\t\t NA \t NA \t NA";
				else
				{
					outfile << "\t"; 
					for (unsigned i = 0; i < miter->second.size(); i++)
					{
						sprintf(buf, "\t%5.3f", miter->second.at(i)); 
						outfile << buf; 
					}
				}
				outfile << endl; 
					
			}
		}
		else
		{
			outfile << "## bf = log10_BF\t se = log10_std_err\t" << endl; 
			outfile << "rsnum\tpos"; 
			for (int p = 0; p < nPH; p++)
			{
				outfile << "\tbf" << p; 
				if(total > 1) 
					outfile << "\tse" << p; 
			}
			outfile << endl; 
			for (int m = 0; m < nLoci; m++)
			{
				char buf[1000];
				sprintf(buf, "%-15s\t%ld\t", vsRsnum.at(m).data(), mapRs2pos[vsRsnum.at(m)]); 
				outfile << buf; 
				for (int p = 0; p < nPH; p++)
				{  
					double t1, t2; 
					t1 = bf1[p * nLoci + m];
					if(bf2 == NULL)
						t2 = 0.0; 
					else
						t2 = bf2[p * nLoci + m];
					
					if(total > 1) 
					{
						double tmax = max(t1 * 2.0, t2); 
						double diff1 = t1 * 2.0 - tmax; 
						double diff2 = t2 - tmax; 
						t2 = 0.5 * (tmax + log(exp(diff2) - exp(diff1)) - log(total)); 
					}
					t1 /= log(10.0); 
					t2 /= log(10.0); 
					
					if(total > 100) 
						sprintf(buf, "%+5.3f %+5.3f\t", t1, t2);
					else
						sprintf(buf, "%+5.3f  NA   \t", t1);
					outfile << buf;
				}
				outfile << endl; 
			}
		}
		outfile.close(); 

		if(nLevel < 2) 
		{
			if(nPH == 1) 
			{
				double logbfmax = 0; 
				for (int m = 0; m < nLoci; m ++)
				{
					if(logbfmax < bf1[m]) 
						logbfmax = bf1[m];
				} 

				double factor = 0; 
				for (int m = 0; m < nLoci; m++)
				{
					factor += exp(bf1[m] - logbfmax); 
				}
				double logbfmean = logbfmax + log(factor / nLoci); 
				write_summary(logbfmean, NULL);
			}
			else 
				write_summary(-1, NULL); 
		}
	}
	
	Free1D((void *) bf1); 
	Free1D((void *) bf2); 
}

void ModelnData::single_snp_pvalue(int nShuffle)  //p-value for single phenotype; 
{
	if(nLoci == 0) return; 

	m_pv.clear(); 
	m_pv.insert(m_pv.end(), nLoci+1, 0); 
	vector<int> perm; 
	perm = vv_phdex.at(0); 
	
	real * bfv = new real[nLoci]; 
 	real * bf0 = new real[nLoci];    
	char str[100];
	sprintf(str, "-bimbam: calculate p value of %d permutations: ", nShuffle); 
	
	int prg_cur = 0;
	int prg_all = nLoci * (nShuffle + 1); 
	int prg_mod = (int) (prg_all / 100.0); 
	if(prg_mod == 0) prg_mod = 1; 

	if(nImpute == 0)   //p-val based on exact genotype; 
	{
		real * snpgt = new real[nCohort];   //single snp genotype; 
		real * snpph = new real[nCohort];   //single snp phenotype; 
		for (int p = 0; p < nShuffle+1; p++)
		{
			vector<int>* base = &vv_phdex.at(0);   
			vector<real>* phval = &vv_phval.at(0);
			
			if(p > 0) 
			{
				random_shuffle(perm.begin(), perm.end());
			}
			for (int m = 0; m < nLoci; m++)
			{   
				int ni = 0; 
				for (unsigned i = 0; i < base->size(); i++)
				{
					int g = base->at(i); 
					int t = perm.at(i); 
					short sgt = pIndiv[g]->GetsnpGT(m); 
					if(sgt == QQ) continue;
					else 
					{
						snpgt[ni] = (real) sgt; 
						snpph[ni] = phval->at(t); 
					}
					ni++;  
				}

				real tbf = 0; 
				if(cc) tbf = cc_bf(0.2, 0.05, snpph, snpgt, ni); 
				else tbf = calc_bf(0.2, 0.05, snpph, snpgt, ni); 
				bfv[m] = tbf; 
				
				prg_cur++; 
				if(prg_cur % prg_mod == 0) 
					print_progress_bar(0, str, prg_cur, prg_all); 
			}   // finish bf for one permute;   
		
			real sumbf0 = 0; 
			real sumbf = 0; 
			if(p == 0) 
			{
				for (int m = 0; m < nLoci; m++)
				{
					bf0[m] = bfv[m]; 
					sumbf0 += bf0[m];
				}
			}
			else 
			{
				for(int m = 0; m < nLoci; m++)
				{
					sumbf += bfv[m]; 
					if(bfv[m] > bf0[m]) m_pv.at(m) += 1.0; 
				}
				if(sumbf > sumbf0) m_pv.at(nLoci) += 1.0; 
			}   //count; 
			
		} // permute loop;   
		delete[] snpgt; 
		delete[] snpph; 
	}
	else   //p-val based on mean genotype; 
	{
		real ** snppr = Allocate2DMatrix(nCohort, 2); 	
		for (int p = 0; p < nShuffle+1; p++)
		{
			vector<int>* base = &vv_phdex.at(0);   
			vector<real>* phval = &vv_phval.at(0);
			
			if(p > 0) random_shuffle(perm.begin(), perm.end());
			
			for (int m = 0; m < nLoci; m++)
			{   
				int ni = 0; 
				for (int i = 0; i < nIndiv; i++)
				{
					if(pIndiv[i]->GetisPanel()) continue;
					real res[2];
					pIndiv[i]->get_snp_dstr(m, res); 
					snppr[ni][0] = res[0];
					snppr[ni][1] = res[1];
					real p2 = 1 - res[0] - res[1]; 
					if(m_not_snp == 0 && (res[0] < -1e-3 || res[1] < -1e-3 || p2 < -1e-3)) 
					{
						cout << "-bimbam: illegal prob. in single_snp_pvalue" << res[0] << "\t" << res[1] << endl; 
						safe_exit(); 
					}
					ni++;
				}   // prepare genotype distribution; 

				real tbf = 0; 
				if(cc) tbf = cc_bf_mean(0.2, 0.05, phval, base, &perm, snppr); 
				else tbf = calc_bf_mean(0.2, 0.05, phval, base, &perm, snppr); 
				bfv[m] = tbf; 
				prg_cur++; 
				if(prg_cur % prg_mod == 0) 
					print_progress_bar(0, str, prg_cur, prg_all); 
			}   // snp loop; 
		
			
			real sumbf0 = 0; 
			real sumbf = 0; 
			if(p == 0) 
			{
				for (int m = 0; m < nLoci; m++)
				{
					bf0[m] = bfv[m]; 
					sumbf0 += bf0[m];
				}
			}
			else 
			{
				for(int m = 0; m < nLoci; m++)
				{
					sumbf += bfv[m]; 
					if(bfv[m] > bf0[m]) m_pv.at(m) += 1.0; 
				}
				if(sumbf > sumbf0) m_pv.at(nLoci) += 1.0; 
				
			}
		} // permute loop;   
		Free2DMatrix(snppr);
	}
	print_progress_bar(1, str, 1, 1);
			
	for (int m = 0; m < nLoci+1; m++)
   	{
		m_pv.at(m) = (m_pv.at(m) + 1.) / (nShuffle+1); 
	}
    delete[] bf0; 
   	delete[] bfv; 
}

void ModelnData::write_summary(double total_weighted_bf, real ** marginal_bf2)
{
	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".summary.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write:" << sfn << endl;
		return;
	}

	if (marginal_bf2 == NULL)
	{
		outfile << "## Summary of overall evidence for association with region:" << endl; 
		if(total_weighted_bf < 0) 
			outfile << "## log10(BF for region) = NA" << endl;  
		else
			outfile << "## log10(BF for region) = " << total_weighted_bf / log(10.0) << endl;  
		if (m_pv.size() == 0) 
			outfile << "## p-value for region = NA" << endl;
		else 
			outfile << "## p-value for region = " << m_pv.at(nLoci) << endl;
	}
	else {
	
		if(nMultiSnp < nLoci) 
			outfile << "## Warning: summary statistics for Multi-SNP model are no longer valid "
				<< "because nMultiSnp != nLoci" << endl;
		fplog << "## Write posterior probability matrix to a file." << endl;
		outfile << "## Summary of overall evidence for association with region:" << endl; 
		outfile << "## log10(BF for region) = " << total_weighted_bf / log(10.0) << endl;  
		if (m_pv.size() == 0) 
			outfile << "## p-value for region = NA" << endl;
		else 
			outfile << "## p-value for region = " << m_pv.at(nLoci) << endl;
		outfile << endl; 
		outfile << "## BFs for choosing among multi-SNP models:" << endl;  

		double maxsumbf = 0; 
		for (int i = 0; i < nLevel; i++)
		{
			double temp =  m_sumbf.at(i) + double(i-1) * log(0.5);
			if (maxsumbf < temp) maxsumbf = temp; 
		}
		
		double grandsum = 0;
		for (int i = 0; i < nLevel; i++)
		{
			grandsum += exp(m_sumbf.at(i) + double(i-1) * log(0.5) - maxsumbf);
			if(nMultiSnp == nLoci) 
				outfile << "## log10(average " << i+1 << "-SNP BF) = " << m_sumbf.at(i) / log(10.0) << endl;  
			else
				outfile << "## log10(average " << i+1 << "-SNP BF) = NA" << endl;   
		}
		outfile << endl;
		
		outfile << "## Example posterior probabilities for multi-SNP models, "
			<< "based on default prior: [p(l) propto 0.5^l for l=1,2,...,L]" << endl;
		for (int i = 0; i < nLevel; i++)
			outfile << "## posterior probability of " << i+1 << "-SNP = " 
				<< exp(m_sumbf.at(i) + double(i-1) * log(0.5) - maxsumbf) / grandsum << endl; 
		
		real mchoose2 = 0.5 * nMultiSnp * (nMultiSnp - 1); 
		outfile <<  endl << "## log10(BF) for single SNP and pair of SNPs:" << endl; 
		double maxbf = 0; 
		for (int m2 = 0; m2 < nMultiSnp; m2++) 
		{
			char buf[100];
			int pos = mrank2loc[m2]; 
			sprintf(buf, "%-15s", vsRsnum.at(pos).c_str()); 
			outfile << buf; 
			for (int m1 = 0; m1 < nMultiSnp; m1++) 
			{
				if(maxbf < marginal_bf2[m2][m1]) maxbf = marginal_bf2[m2][m1]; 
				if(m1 >= m2) sprintf(buf, "%+4.3f\t", marginal_bf2[m2][m1] / log(10.0));
				else sprintf(buf, "%+4.3f\t", marginal_bf2[m1][m2] / log(10.0));
				outfile << buf; 
			}
			outfile << endl;;
		}            
		outfile << endl << "## posterior probabilities for single SNP and pair of SNPs:" << endl; 
		
		real dsum = 0.0; 
		for (int m2 = 0; m2 < nMultiSnp; m2++) 
			for (int m1 = 0; m1 < nMultiSnp; m1++) 
			{
				if(m1 == m2) 
				{
					marginal_bf2[m1][m2] -= log((double)nMultiSnp);
					dsum += exp(marginal_bf2[m1][m2] - maxbf);
				}
				else if(m1 > m2) 
				{
					marginal_bf2[m2][m1] -= log((double)mchoose2); 
					dsum += exp(marginal_bf2[m2][m1] - maxbf);
				}
			}
		for (int m1 = 0; m1 < nMultiSnp; m1++)
			for (int m2 = m1; m2 < nMultiSnp; m2++)
				marginal_bf2[m1][m2] = exp(marginal_bf2[m1][m2] - maxbf) / dsum;
		
		for (int m2 = 0; m2 < nMultiSnp; m2++) 
		{
			char buf[100];
			int pos = mrank2loc[m2]; 
			sprintf(buf, "%-15s", vsRsnum.at(pos).c_str()); 
			outfile << buf; 
			for (int m1 = 0; m1 < nMultiSnp; m1++) {
				if (m1 >= m2) 
					sprintf(buf, "%+4.3f\t", marginal_bf2[m2][m1]);
				else	
					sprintf(buf, "%+4.3f\t", marginal_bf2[m1][m2]);
				outfile << buf; 
			}
			outfile << endl;
		}
	}
	outfile.close(); 
}

void ModelnData::write_genotype(int type, int all)   
//type = 0 exact; type = 1 mean; type=2 best guess; type=3 distribution;
//all=0 cohort-only; all = 1 all-SNPs; 
{
	fstream outfile; 
	string fn("output/");
	if (type == 0)
	{
		fn.append(fnOutput);
		fn.append(".exact.genotype.txt");
		outfile.open(fn.c_str(), ios::out);
		if(!outfile.good()) 
		{
			fplog << "ERROR: BIMBAM failed to open file to write." << endl;
			cout << "-bimbam: bimbam failed to open file to write" << endl;
			return; 
		}
		
		if (all == 0)   //write cohort in numeric.   
		{
    		outfile << "### 2 minor allele homozygous, NA missing genotype " << endl; 
			for (unsigned n = 0; n < vPARCpos.size(); n++)
			{   
				int m = vPARCpos.at(n);
				string rs(vsRsnum.at(m)); 
				outfile << rs << " "; 
				if(m_allele_coding_mode == 0)
					outfile << mapRs2mm[rs].second << " " << mapRs2mm[rs].first; 
				else //(m_allele_coding_mode == 1) 
					outfile << mapRs2mm[rs].first << " " << mapRs2mm[rs].second; 
				for(int i = 0; i < nIndiv; i++)
				{
					short gt = 0; 
					if(pIndiv[i]->GetisPanel()) continue; 
					gt = pIndiv[i]->GetsnpGT(m);
					
					if(gt == QQ)
						outfile << " NA ";
					else 
						outfile << " " << gt; 
				}
				outfile << endl; 
			}					
		}
		else           //write cohort in bimbam format. 
		{   
			outfile << nDip << endl; 
			outfile << vPARCpos.size() << endl; 
			for (unsigned n = 0; n < vPARCpos.size(); n++)
			{
				int m = vPARCpos.at(n);
				string rs(vsRsnum.at(m)); 
				outfile << rs << " "; 
				for (int i = 0; i < nIndiv; i++)
				{
					short gt = 0; 
					if(pIndiv[i]->GetisPanel()) continue; 
					gt = pIndiv[i]->GetsnpGT(m);
					if(gt == QQ) 
						outfile << " NN"; 
					else 
					{
						if(m_allele_coding_mode == 0) 
						{
							if(gt == 0)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
							else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
						else
						{
							if(gt == 2)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
							else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
						
					}
				}
				outfile << endl; 
			}
		}
		outfile.close(); 
		return ; 
	}

#if defined (MPI_ENABLED)
	if (procID == 0)
#endif
	{
		if(type == 1)
		{
			fn.append(fnOutput);
			fn.append(".mean.genotype.txt");
			cout <<  "-bimbam: write mean genotype" << endl;
		}
		else if(type == 2) 
		{
			fn.append(fnOutput);
			fn.append(".best.guess.genotype.txt");
			cout <<  "-bimbam: write best guess genotype" << endl;
		}
		else if(type == 3)
		{
			fn.append(fnOutput);
			fn.append(".genotype.distribution.txt");
			cout << "-bimbam: write genotype distribution" << endl;
		}
		
		outfile.open(fn.c_str(), ios::out);
		if(!outfile.good()) 
		{
			fplog << "ERROR: BIMBAM failed to open file to write." << fn << endl;
			cout << "-bimbam: failed to open file to write" << fn << endl;
			return; 
		}
			
	}
	
#if defined (MPI_ENABLED)
	real ** snpInpr = NULL;
    long bufsize = nMaxLoad * nLoci * 2;
	real ** sendbuf = Allocate2DMatrix(nMaxLoad, nLoci * 2); 
    MPI_Status *stat = new MPI_Status; 
	if (procID == 0)
		snpInpr = Allocate2DMatrix(nCohort, nLoci * 2);
#endif 
	
#if defined (MPI_ENABLED)		
	int np = 0; 
	for (unsigned int wp = 0; wp < pvLoadDip[procID].size(); wp++)
	{
		int i = pvLoadDip[procID].at(wp);
		if(pIndiv[i]->GetisPanel()) continue; 
		real * snp_dens = pIndiv[i]->get_snp_dstr();
		for (int m = 0; m < nLoci * 2; m++)
			sendbuf[np][m] = snp_dens[m];
		np++;                
	}   //parepare the send buffer. 
	
    if (procID > 0)   //should look into how to use Isend Irecv to speed this up; 
    {
        if (sizeof(real) == sizeof(double))
            MPI_Send(sendbuf[0], bufsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        else
            MPI_Send(sendbuf[0], bufsize, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
    }
    else if (procID == 0)
    {
		for (unsigned j = 0; j < vvLoad.at(0).size(); j++)
		{
            int i = vvLoad.at(0).at(j);
			for (int m = 0; m < nLoci * 2; m++)
				snpInpr[i][m] = sendbuf[j][m]; 
		}
        for (int p = 1; p < nProc; p++)
        {
            if (sizeof(real) == sizeof(double))
                MPI_Recv(sendbuf[0], bufsize, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, stat);
            else
                MPI_Recv(sendbuf[0], bufsize, MPI_FLOAT, p, 1, MPI_COMM_WORLD, stat);

			for (unsigned j = 0; j < vvLoad.at(p).size(); j++)
			{
				int i = vvLoad.at(p).at(j);
				for (int m = 0; m < nLoci * 2; m++)
					snpInpr[i][m] = sendbuf[j][m]; 
			}
        }
    }

	Free2DMatrix(sendbuf); 

    if(procID == 0)
#endif
	{
		if(all == 1)        //write all SNPs (panel and cohort) not all individual!
		{
			for (int m = 0; m < nLoci; m++) 
			{
				string rs(vsRsnum.at(m)); 
				outfile << rs;             
				if(type != 2)              
				{
					if(m_allele_coding_mode == 0) 
						outfile << " " << mapRs2mm[rs].second << " " << mapRs2mm[rs].first;
					else
						outfile << " " << mapRs2mm[rs].first << " " << mapRs2mm[rs].second;
				}
				//always output rs, if output numerical genotype, output allele coding as well. 
				char buf[100];
				for(int i = 0; i < nIndiv; i++)
				{
					real pr2[2]; 
					if(pIndiv[i]->GetisPanel()) continue; 
					pIndiv[i]->get_snp_dstr(m, pr2); 
					real p0 = pr2[0];
					real p1 = pr2[1];
					real p2 = 1.0 - p0 - p1;
					if(type == 1)
					{
						double mean = 0; 
						mean = p1 + 2.0 * p2;   
						if(mean >= 0 && mean <= 2.0)
							sprintf(buf, " %.3f", mean); 
						else
							outfile << " NA "; 
						outfile << buf; 
					}
					else if (type == 3) 
					{
						sprintf(buf, " %.3f %.3f ", p0, p1); 
						outfile << buf; 
					}
					else // (type == 2)
					{
						int gt = 0; 
						real max = p0; 
						if(p1 > max) {gt = 1; max = p1;}
						if(p2 > max) gt = 2; 
						
						if(m_allele_coding_mode == 0) 
						{
							if(gt == 0) outfile << " " << mapRs2mm[rs].first <<  mapRs2mm[rs].first; 
							else if (gt == 1) outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
						else 
						{
							if(gt == 2) outfile << " " << mapRs2mm[rs].first <<  mapRs2mm[rs].first; 
							else if (gt == 1) outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
					}
				}
				outfile << endl; 
			}
		}
		else  //all == 0; write cohort snps only; 
		{
        	for (unsigned n = 0; n < vPARCpos.size(); n++)
			{
	        	int m = vPARCpos.at(n);
				string rs(vsRsnum.at(m)); 
				outfile << rs;             
				if(type != 2) 
				{
					if(m_allele_coding_mode == 0) 
						outfile << " " << mapRs2mm[rs].second << " " << mapRs2mm[rs].first;
					else
						outfile << " " << mapRs2mm[rs].first << " " << mapRs2mm[rs].second;
				}
				//always output rs, if output numerical genotype, output allele coding as well. 
				char buf[100];
				for(int i = 0; i < nIndiv; i++)
				{
					real pr2[2]; 
					if(pIndiv[i]->GetisPanel()) continue; 
					pIndiv[i]->get_snp_dstr(m, pr2); 
					real p0 = pr2[0];
					real p1 = pr2[1];
					real p2 = 1.0 - p0 - p1;
					if(type == 1)
					{
						double mean = 0; 
						mean = p1 + 2.0 * p2;   
						if(mean >= 0 && mean <= 2.0)
							sprintf(buf, " %.3f", mean); 
						else
							outfile << " NA "; 
						outfile << buf; 
					}
					else if (type == 3) 
					{
						sprintf(buf, " %.3f %.3f ", p0, p1); 
						outfile << buf; 
					}
					else // (type == 2)
					{
						int gt = 0; 
						real max = p0; 
						if(p1 > max) {gt = 1; max = p1;}
						if(p2 > max) gt = 2; 
						if(m_allele_coding_mode == 0) 
						{
							if(gt == 0) outfile << " " << mapRs2mm[rs].first <<  mapRs2mm[rs].first; 
							else if (gt == 1) outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
						else 
						{
							if(gt == 2) outfile << " " << mapRs2mm[rs].first <<  mapRs2mm[rs].first; 
							else if (gt == 1) outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
							else outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
						}
					}
				}
				outfile << endl; 
			}
		}
		outfile.close(); 
	}
#if defined (MPI_ENABLED)
	Free2DMatrix(snpInpr);
#endif
}


void ModelnData::multi_snp_analysis(int mode)
{
	if(nMultiSnp < 2 || nMultiSnp > nLoci) 
	{
		nMultiSnp = nLoci; 
	}
	if(nLevel < 2 || nLevel > nMultiSnp) 
	{
		return; 
		fplog << "## BIMBAM: Parameter in multi_snp_analysis out of range." << endl; 
		cout << "-bimbam: parameter in multi-snp analysis out of range." << endl; 
	}
	
 	vector<long> nChoosek; 
	nChoosek.push_back(1);
	for (int m = 1; m < nMultiSnp+1; m++) 
	{
		long temp = nChoosek.at(m-1) * (nMultiSnp-m+1) / m;
		nChoosek.push_back(temp); 
	}   //nMultiSnp choose k; 
	
	long snp_combo  = 0; 
	for (int i = 1; i <= nLevel; i++)
		snp_combo += nChoosek.at(i);
	if (snp_combo > (long) vPerm.max_size()) 
	{
	   	cout << "-bimbam: number of combintations over limit: " << snp_combo << " vs. " << vPerm.max_size() << endl;
		safe_exit(); 
	}
	
	vPerm.clear(); 
	vPerm.reserve(snp_combo);
	for (int i = 0; i < nMultiSnp; i++)
	{
		PosnBF tmp;
		tmp.pos.push_back(i);
		tmp.bf = 0;
		vPerm.push_back(tmp);
	} //single snp;
	
	for (int i = 0; i < (int)vPerm.size(); i++)
	{
		if ((long)vPerm.size() >= snp_combo)
			break; 
		PosnBF tmp = vPerm.at(i);
		if((int)tmp.pos.size() > nLevel)
			continue;
		for (int j = *(tmp.pos.end()-1) + 1; j < nMultiSnp; j++)
		{
			tmp.pos.push_back(j);
			tmp.bf = 0;
			vPerm.push_back(tmp);
			if((int)vPerm.size() > snp_combo) break;
			tmp.pos.pop_back(); //since we changed tmp in last step, we recover it;
		}
	}  // finish preparing vPerm;
	if(vPerm.size() != (unsigned) snp_combo)
	{
		cout << "vPerm dimension wrong! " << vPerm.size() << "\t" << snp_combo <<  endl; 
		exit(0); 
	}
	
	vector<real>* vph = &vv_phval.at(0); 
	vector<int>* phindex = &vv_phdex.at(0); 
	
	
	int * index = new int[nMultiSnp];
	//which snp in multi_snp study is missing.  
	{
		for (int r = 0; r < nMultiSnp; r++)
			index[r] = mrank2loc[r]; 
	}
	
	if (mode == 1)  //mean genotypes; 
	{   
		real ** snpmgt = Allocate2DMatrix(nCohort, nMultiSnp);
		char str[100] = "-bimbam: multi-snp analysis of exact genotypes: "; 
		int ni = 0; 
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
			for (int m = 0; m < nMultiSnp; m++)
			{
				real tt = pIndiv[i]->get_snpmgt(index[m]); 
				snpmgt[ni][m] = tt; 
			}
			ni++; 
		}
		
		int prg_all = (int) vPerm.size(); 
		int prg_mod = (int) (prg_all / 100.0); 
		if(prg_mod == 0) prg_mod = 1; 
		for (unsigned i = 0; i < vPerm.size(); i++)
		{
			PosnBF * pt = &vPerm.at(i);
			double bf_val = 0.0;
			if(cc) bf_val = cc_bf_mgt(0, 0, vph, phindex, snpmgt, pt);
			else bf_val = calc_bf_mgt(0, 0, vph, phindex, snpmgt, pt);
			pt->bf = bf_val; 
			if(i % prg_mod == 0) print_progress_bar(0, str, i, prg_all); 
		}   
		print_progress_bar(1, str, 1, 1); 
		Free2DMatrix(snpmgt); 
	} 
	
	else if (mode == 2) //genotype distribution;  
	{
		real ** snpInpr = Allocate2DMatrix(nCohort, nMultiSnp * 2);
		char str[100] = "-bimbam: multi-snp analysis of mean genotypes: ";  
		{
			int ni = 0; 
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->GetisPanel()) continue; 
	    		real * snp_dens = pIndiv[i]->get_snp_dstr();
				for (int m = 0; m < nMultiSnp; m++)
				{  
					snpInpr[ni][m * 2] = snp_dens[index[m] * 2];
					snpInpr[ni][m*2+1] = snp_dens[index[m]*2+1];
				}
				ni++; 
			}
		}
		
		int prg_all = (int) vPerm.size(); 
		int prg_mod = (int) (prg_all / 100.0); 
		if(prg_mod == 0) prg_mod = 1; 
		for (unsigned i = 0; i < vPerm.size(); i++)
		{
			PosnBF * pt = &vPerm.at(i);
			double bf_val = 0.0;
			if(cc) bf_val = cc_bf(0, 0, vph, phindex, snpInpr, pt);
			else bf_val = calc_bf(0, 0, vph, phindex, snpInpr, pt);
			pt->bf = bf_val;
			if(i % prg_mod == 0) print_progress_bar(0, str, i, prg_all); 
		}   
		print_progress_bar(1, str, 1, 1); 
		Free2DMatrix(snpInpr); 
	}
	
	else if (nImpute == 0) 
	{
        short ** snpIn01 = Allocate2DShortMatrix(nCohort, nMultiSnp);
        short * snp_imp = new short[nMultiSnp];
        char str[100] = "-bimbam: multi-snp analysis of exact genotypes: ";
        int ni = 0;
        for (int i = 0; i < nIndiv; i++)
        {
            if(pIndiv[i]->GetisPanel()) continue;
            pIndiv[i]->joint_imputation(pMP, -1, snp_imp, nLoci, nK, nMultiSnp, index);
            for (int m = 0; m < nMultiSnp; m++)
                snpIn01[ni][m] = snp_imp[m];
            ni++;
        }
        delete[] snp_imp;

        int prg_all = (int) vPerm.size();
        int prg_mod = (int) (prg_all / 100.0);
        if(prg_mod == 0) prg_mod = 1;
        for (unsigned i = 0; i < vPerm.size(); i++)
        {
            PosnBF * pt = &vPerm.at(i);
            double bf_val = 0.0;
            if(cc) bf_val = cc_bf(0, 0, vph, phindex, snpIn01, pt);
            else bf_val = calc_bf(0, 0, vph, phindex, snpIn01, pt);
            pt->bf = bf_val;
            if(i % prg_mod == 0) print_progress_bar(0, str, i, prg_all);
        }
        print_progress_bar(1, str, 1, 1);
        Free2DShortMatrix(snpIn01); 

	}
	
	else if (nImpute > 1)     //jont imputation; 
	{
		short ** snpIn01 = Allocate2DShortMatrix(nCohort, nMultiSnp);
		int unitrepeat = (int)((real) nImpute / nEMRuns);
		if(unitrepeat < 1) unitrepeat = 1; 
		nImpute = unitrepeat * nEMRuns; 
		char str[100] = "-bimbam: multi-snp analysis with joint imputation: ";  
		int prg_cur = 0; 
		int prg_all = nImpute; 
		int prg_mod = (int) (nImpute / 100.0); 
		if(prg_mod == 0) prg_mod = 1; 
		double count = 0.0; 
		for (int runs = 0; runs < nEMRuns; runs++)
		{    
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->GetisPanel()) continue; 
				pIndiv[i]->CalcAll(nLoci, nK, pMP+runs, 0);
			}
			for (int repeat = 0; repeat < unitrepeat; repeat++)
			{
				int ni = 0; 
				for (int i = 0; i < nIndiv; i++)
				{
					if(pIndiv[i]->GetisPanel()) continue;
					 pIndiv[i]->joint_imputation(pMP, runs, snpIn01[ni++], nLoci, nK, nMultiSnp, index);
				}

				for (int i = 0; i < (int)vPerm.size(); i++)
				{
					PosnBF * pt = &vPerm.at(i);
					double bf_val = 0.0;
					if(cc) bf_val = cc_bf(0, 0, vph, phindex, snpIn01, pt);
					else bf_val = calc_bf(0, 0, vph, phindex, snpIn01, pt);
					double maxlog = max((double)pt->bf, bf_val); 
					double diff0 = pt->bf - maxlog; 
					double diff1 = bf_val -maxlog; 
					pt->bf = maxlog + log(count / (count + 1.0) * exp(diff0) + 1.0 / (count + 1.0) * exp(diff1)); 
					maxlog = max((double)pt->var, bf_val*2.0); 
					diff0 = pt->var - maxlog; 
					diff1 = bf_val*2.0 -maxlog; 
					pt->var = maxlog + log(count / (count + 1.0) * exp(diff0) + 1.0 / (count + 1.0) * exp(diff1)); 
				}  
				count += 1.0; 
				prg_cur++;
				if(prg_cur % prg_mod == 0) print_progress_bar(0, str, prg_cur, prg_all); 
			}
		}
		print_progress_bar(1, str, 1, 1); 
		Free2DShortMatrix(snpIn01); 
	}
	
	{
		double total = (double) nImpute; 
		if (total <= 1) total = 1; 
		
	   	//wirte multi-snp bfs. 
		fstream outfile; 
		string sfn("output/");
		sfn.append(fnOutput);
		sfn.append(".multi.txt");
		outfile.open(sfn.c_str(), ios::out);
		if(!outfile.is_open()) 
		{
			cout << "-bimbam: cannot open file to write:" << sfn << endl;
			return;
		}
		outfile << "## note:bf=log10(Bayes Factor), mapping between SNP ID and its label(rank) are in .single file." << endl;
		outfile << " bf \t se"; 
		for (int m = 0; m < nLevel; m++) 
			outfile << "\t\tsnp" << m+1; 
		outfile << endl; 
		for (unsigned i = 0; i < vPerm.size(); i++)
		{
			PosnBF * pt = &vPerm.at(i);
			char str[100]; 
			double t1 = pt->bf; 
			double t2 = pt->var; 
			
			if(total <= 1)
				sprintf(str, "%+5.3f  NA   ", t1 / log(10.0)); 
			else 
			{
				double tmax = max(t1*2.0, t2); 
				double diff1 = t1*2.0 - tmax; 
				double diff2 = t2 - tmax; 
				t2 = tmax + log(exp(diff2) - exp(diff1)) - 0.5 * log(total); 
				sprintf(str, "%+5.3f  %+5.3f", t1 / log(10.0), t2 / log(10.0)); 
			}
			outfile << str; 
			for (unsigned m = 0; m < pt->pos.size(); m++)
				outfile << "\t\t" << pt->pos.at(m)+1; 
			for (int m = (int)pt->pos.size(); m < nLevel; m++)
				outfile << "\t\tNA"; 
			outfile << endl; 
		}
		outfile.close(); 
		
		double total_weighted_bf = 0;
		m_sumbf.clear(); 
		m_sumbf.insert(m_sumbf.begin(), nLevel, 0); 
		double geometric = 0.0; 
		for (int i = 1; i <= nLevel; i++)
			geometric += pow(0.5, (double) i);
		
		double maxbf = 0.0; 
		for (int i = 0; i < (int) vPerm.size(); i++)
		{
			if(maxbf < vPerm.at(i).bf) maxbf = vPerm.at(i).bf; 
		}
		
		
		for (int i = 0; i < (int) vPerm.size(); i++)
		{
			PosnBF * pt = &vPerm.at(i);
			int size = pt->pos.size();
			m_sumbf.at(size-1) += exp(pt->bf - maxbf); 
		}   // nChoosek.at(size)
		
		for (int i = 0; i < nLevel; i++)
		{
			m_sumbf.at(i) /= nChoosek.at(i+1); 
			total_weighted_bf += m_sumbf.at(i) * pow(0.5, double(i+1)) / geometric;
			m_sumbf.at(i) = maxbf + log(m_sumbf.at(i)); 
		}
		total_weighted_bf = maxbf + log(total_weighted_bf); 
		
		real ** marginal_bf2 = Allocate2DMatrix(nMultiSnp, nMultiSnp);
		for (int i = 0; i < (int) vPerm.size(); i++)
		{
			PosnBF * pt = &vPerm.at(i);
			int size = (int) pt->pos.size(); 
			if(size > 2) continue; 
			if(size == 1) marginal_bf2[pt->pos.at(0)][pt->pos.at(0)] = pt->bf;
			else marginal_bf2[pt->pos.at(0)][pt->pos.at(1)] = pt->bf; 
		} 
	   
		write_summary(total_weighted_bf, marginal_bf2); 
		Free2DMatrix(marginal_bf2);   
	}
	return; 
}    

void ModelnData::single_snp_exact(int mode, real * bf1)
{
	char str[100]; 
	if(mode == 1) 
		sprintf(str, "-bimbam: single snp analysis of typed snps: "); 
	else if (mode == 2) 
		sprintf(str, "-bimbam: single snp  f-stat of typed snps: "); 
	else if (mode == 3) 
		sprintf(str, "-bimbam: single snp lrt statistics of typed snps: "); 
		
	short ** snpIn01 = NULL;
	real ** bf = NULL; 
	int numPH = nPH; 
	snpIn01 = Allocate2DShortMatrix(nCohort, nLoci);
	bf = Allocate2DMatrix(numPH, nLoci); 

	int ni = 0;
	for (int i = 0; i < nIndiv; i++)
	{
		if (pIndiv[i]->GetisPanel()) continue; 
		pIndiv[i]->GetsnpGT(nLoci, snpIn01[ni++]); 
	}
  	int len = nLoci; 
	
	real * snpgt = new real[nCohort];   //single snp genotype; 
	real * snpph = new real[nCohort];   //single snp phenotype; 
	
	int prg_cur = 0; 
	int prg_all = numPH * len; 
	int prg_mod = (int) (numPH * len / 100.0); 
	if (prg_mod == 0) prg_mod = 1;
	for (int p = 0; p < numPH; p++)
	{   
		vector<int> * phdex = &vv_phdex.at(p);
		vector<real> * phval = &vv_phval.at(p); 
		
		for (int m = 0; m < len; m++)
		{
			int ni = 0; 
			for (unsigned i = 0; i < phdex->size(); i++)
			{
				int g = phdex->at(i); 
				if (snpIn01[g][m] >= 0 && snpIn01[g][m] <=2)
				{
					snpgt[ni] = (real) snpIn01[g][m];
					snpph[ni] = phval->at(g);
					ni++; 
				}
			}
			real tbf = 0; 
			if (mode == 1)
			{
				if(cc) tbf = cc_bf(0, 0, snpph, snpgt, ni);
				else tbf = calc_bf(0, 0, snpph, snpgt, ni);
				if(numPH == 1) 
					mapRs2beta[vsRsnum.at(m)] = m_beta; 
			}
			else if(mode == 2) tbf = f_stat(snpph, snpgt, ni); 
			else if(mode == 3) tbf = snp_lrt(snpph, snpgt, ni); 
			bf[p][m] = tbf;
			prg_cur++; 
			if(prg_cur % prg_mod == 0) print_progress_bar(0, str, prg_cur, prg_all); 
		} //loci loop
	}  //pheno loop
	delete[] snpph;
	delete[] snpgt; 
	Free2DShortMatrix(snpIn01); 

	for (int r = 0; r < numPH; r++)
		for (int m = 0; m < nLoci; m++)
			bf1[r * nLoci + m] = bf[r][m]; 
   	print_progress_bar(1, str, 1, 1); 
	
    Free2DMatrix(bf); 
}

void ModelnData::single_snp_importance_sampling(real * bf1, real * bf2)
{    
	real ** bf = NULL;   //mean and variance. 
	real ** snpInpr = NULL; 
	real * phenoVal = new real[nCohort]; 
	real * one_snp = new real[nCohort];
	real ** snppr = Allocate2DMatrix(nCohort, 3);  
	gsl_matrix * XtX = gsl_matrix_alloc(3, 3); 
	
	real ** q = Allocate2DMatrix(nCohort, 3);
	snpInpr = Allocate2DMatrix(nCohort, nLoci * 2);
	bf = Allocate2DMatrix(nCohort, nLoci * 2); 
	gsl_vector * xy = gsl_vector_alloc(3); 
    gsl_vector * beta = gsl_vector_alloc(3); 
	
	char str[100] = "-bimbam: single snp analysis with importance sampling: ";

	int prg_cur = 0; 
	int prg_all = nImpute * nLoci;
	int prg_mod = (int)(prg_all / 100.0);
	if(prg_mod == 0) prg_mod = 1 ;
	
	int ni = 0; 
	for (int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->GetisPanel()) continue; 
		real * snp_dens = pIndiv[i]->get_snp_dstr();
		for (int m = 0; m < nLoci * 2; m++)
			snpInpr[ni][m] = snp_dens[m]; 
		ni++; 
	}
	int loadlen = nLoci; 

	gsl_vector * logbf = gsl_vector_alloc(nImpute); 
	for (int m = 0; m < loadlen; m++)
	{
		for (int i = 0; i < nCohort; i++)
		{
			snppr[i][0] = snpInpr[i][2 * m];
			snppr[i][1] = snpInpr[i][2*m+1];
			snppr[i][2] = 1.0 - snppr[i][0] - snppr[i][1];
			if(m_not_snp == 0 && (snppr[i][0] < -1e-3 || snppr[i][1] < -1e-3 || snppr[i][2] < -1e-3))
			{
				safe_exit(); 
				cout << "-bimbam: illegal probabilities in importance sampling. " << endl; 
			}
		}

		for (int np = 0; np < nPH; np++)
		{
			vector<real>* phval = &vv_phval.at(np); 
			vector<int>* phdex = &vv_phdex.at(np);   
			int nCohort = (int) phdex->size();

			gsl_matrix * X = gsl_matrix_alloc(nCohort, 3); 
			gsl_vector * gph = gsl_vector_alloc(nCohort); 
			
			for (int i = 0; i < nCohort; i++)
			{
				int j = phdex->at(i);
				phenoVal[i] = phval->at(j);
				gsl_vector_set(gph, i, phval->at(j)); 
				gsl_matrix_set(X, i, 0, 1.0); 
				gsl_matrix_set(X, i, 1, snppr[j][1] + 2.0 * snppr[j][2]);
				gsl_matrix_set(X, i, 2, snppr[j][1]); 
			}
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XtX); 
			gsl_blas_dgemv(CblasTrans, 1.0, X, gph, 0.0, xy);  
			
			real var = 0.2 * 0.2; 
			gsl_matrix_set(XtX, 1, 1, gsl_matrix_get(XtX, 1, 1) + 1.0 / var);
			gsl_matrix_set(XtX, 2, 2, gsl_matrix_get(XtX, 2, 2) + 16.0 / var);

			gsl_linalg_cholesky_decomp(XtX); 
			gsl_linalg_cholesky_solve(XtX, xy, beta); //beta = (XtX)^-1 %*% xy; 
			
			for (int i = 0; i < nCohort; i++)
			{        
				double b0 = gsl_vector_get(beta, 0); 
				double b1 = gsl_vector_get(beta, 1); 
				double b2 = gsl_vector_get(beta, 2); 
				double ph = gsl_vector_get(gph, i); 
				
				int j = phdex->at(i); 
				q[i][0] = snppr[j][0] * gsl_ran_gaussian_pdf(ph - b0, 1.0);
				q[i][1] = snppr[j][1] * gsl_ran_gaussian_pdf(ph - (b0 + b1 + b2), 1.0);
				q[i][2] = snppr[j][2] * gsl_ran_gaussian_pdf(ph - (b0 + b1 + b1), 1.0);   	
				//generate q;
				real sum = q[i][0] + q[i][1] + q[i][2];
				q[i][0] /= sum; q[i][1] /= sum; q[i][2] /= sum; 
			}   //normalization;       		   
			
			//now generate samples based on q.
			for (int repeat = 0; repeat < nImpute; repeat++)
			{    
				double logpwght = 0;
				double logqwght = 0;                       
				double logratio = 0;    
				if (gsl_rng_uniform(gsl_r) < 0.9)
				{ 
					for (int i = 0; i < nCohort; i++)
					{
						int j = phdex->at(i); 
						real unif = gsl_rng_uniform(gsl_r); 
						if(unif < q[i][0])
						{
							one_snp[i] = 0; 
							logpwght += log(snppr[j][0]);
							logqwght += log(q[i][0]);
						}
						else if(unif < (q[i][0] + q[i][1]))
						{
							one_snp[i] = 1; 
							logpwght += log(snppr[j][1]);
							logqwght += log(q[i][1]);    						
						}
						else 
						{
							one_snp[i] = 2;
							logpwght += log(snppr[j][2]);
							logqwght += log(q[i][2]);
						}
					}
					logratio = (logpwght - logqwght);    
				}
				else
				{   
					for (int i = 0; i < nCohort; i++)
					{
						int j = phdex->at(i); 
						real unif = gsl_rng_uniform(gsl_r); 
						if(unif < snppr[j][0])
							one_snp[i] = 0; 
						else if(unif < (snppr[j][0] + snppr[j][1]))
							one_snp[i] = 1; 
						else 
							one_snp[i] = 2;
					}
				}
				
				real tbf = 0; 
				if(cc) tbf = cc_bf(0, 0, phenoVal, one_snp, nCohort);
				else tbf = calc_bf(0, 0, phenoVal, one_snp, nCohort);
				tbf += logratio;        
				gsl_vector_set(logbf, repeat, tbf); 
				prg_cur++; 
			   	if(prg_cur % prg_mod == 0) print_progress_bar(0, str, prg_cur, prg_all);
			}                                
			double maxbf = gsl_vector_max(logbf); 
			gsl_vector_add_constant(logbf, -maxbf); 
			double f1 = 0; 
			double f2 = 0; 
			for (int r = 0; r < nImpute; r++)
			{
				double temp = exp(gsl_vector_get(logbf, r)); 
				f1 += temp; 
				f2 += temp * temp; 
			}
			f1 /= (double) nImpute; 
			f2 /= (double) nImpute; 
			bf[np][2 * m] = maxbf + log(f1);        //log first-moment; 
			bf[np][2*m+1] = 2.0 * maxbf + log(f2);      //log second-moment; 
			gsl_vector_free(gph); 
			gsl_matrix_free(X); 
		} //end pheno. 
	} // end for(m=0;m<nLoci;m++)
  
	print_progress_bar(1, str, 1, 1);
    Free2DMatrix(snpInpr); 
	delete[] phenoVal;
	delete[] one_snp;
	Free2DMatrix(snppr);
	Free2DMatrix(q);
	gsl_vector_free(xy); 
	gsl_vector_free(beta); 
	gsl_matrix_free(XtX); 
	
	
	for (int r = 0; r < nPH; r++)
		for (int m = 0; m < nLoci; m++)
		{              
			bf1[r * nLoci + m] = bf[r][2 * m]; 
			bf2[r * nLoci + m] = bf[r][2*m+1]; 
		}

	Free2DMatrix(bf);  
}  

void ModelnData::single_snp_mean_genotype(int mode, real * bf1)
{    
	int bgQ = 0; 
	mode = mode % 10;          //nImpute = 1, 2, 3, 4, 5;  based on mean genotype; 
	char str[100]; 
	if(bgQ == 0) 
	{
		if(mode == 2) sprintf(str, "-bimbam: single snp f-stat of mean genotypes: "); 
		else if(mode == 3) sprintf(str, "-bimbam: single snp lrt of genotype distributions: "); 
		else  sprintf(str, "-bimbam: single snp bf of mean genotypes: "); 
	}
	else 
	{
		if(mode == 2) sprintf(str, "-bimbam: single snp f-stat of best guess genotypes: "); 
		else if(mode == 3) sprintf(str, "-bimbam: single snp lrt of best guess genotypes: "); 
		else  sprintf(str, "-bimbam: single snp bf of best guess genotypes: "); 
	}
	
	real ** bf = NULL;   //mean bf. 
	real ** snpInpr = NULL; 
	int numPH = nPH; 
	
	bf = Allocate2DMatrix(numPH, nLoci); 

	int len = nLoci; 
	real * snpmgt = new real[nCohort]; 	
	int prg_cur = 0;
	int prg_all = len * numPH; 
	int prg_mod = (int) (len * numPH / 100.0); 
	if(prg_mod == 0) prg_mod = 1; 
	for (int m = 0; m < len; m++)
	{   
		int ni = 0; 
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
	    	snpmgt[ni] = pIndiv[i]->get_snpmgt(m);
			ni++; 
		}
		
		for (int p = 0; p < numPH; p++)
		{
			vector<int> * phdex = &vv_phdex.at(p);   
			vector<real> * phval = &vv_phval.at(p); 
			
			real tbf = 0; 
			real mu, add, dom; 
			mu = add = dom = 0; 
			if (mode == 1) 
			{
				if(cc) tbf = cc_bf_mgt(0, 0, phval, phdex, phdex, snpmgt); 
				else tbf = calc_bf_mgt(0, 0, phval, phdex, phdex, snpmgt); 
				if(numPH == 1) 
					mapRs2beta[vsRsnum.at(m)] = m_beta; 
			}
			else cout << "-bimbam: wrong mode" << endl; 
			bf[p][m] = tbf; 
			prg_cur++; 
			if(prg_cur % prg_mod == 0) 
				print_progress_bar(0, str, prg_cur, prg_all); 
		} //end pheno. 
	} // end for(m=0;m<nLoci;m++)
	delete[] snpmgt; 

	if(snpInpr) Free2DMatrix(snpInpr); 
	for (int r = 0; r < numPH; r++)
		for (int m = 0; m < nLoci; m++)
			bf1[r * nLoci + m] = bf[r][m]; 
	print_progress_bar(1, str, 1, 1); 
	
	Free2DMatrix(bf); 
}  

void ModelnData::open_log(void)
{
#if defined (MPI_ENABLED)
	if(procID == 0)
#endif
	{
		string sfn("output/");
		sfn.append(fnOutput);
		sfn.append(".log.txt");
		fplog.open(sfn.c_str(), ios::out);
		if(!fplog.is_open()) 
		{
			cout << "-bimbam: cannot open log file" << endl;
			safe_exit(); 
		}
	}
}

void ModelnData::close_log(void)
{
#if defined (MPI_ENABLED)
	if(procID == 0)
#endif
		fplog.close();
}

double func_f(const gsl_vector * beta, void * params)
{
    CaseCtrlData * pD = (CaseCtrlData *) params;
    int ni = pD->ni;
    int np = pD->np;

    gsl_vector * xb = gsl_vector_alloc(ni);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pD->mx, beta, 0.0, xb);

    double xyb = 0.0;
    gsl_blas_ddot(pD->vy, xb, &xyb);

    double sum = 0;
    for (int i = 0; i < ni; i++)
        sum += logl(1.0 + expl(gsl_vector_get(xb, i)));

    double prior = 0; //-np / 2.0 * log(6.2832);   //const move to prDM1; 
    for (int k = 0; k < np; k++)
    {
        double tmp = gsl_vector_get(beta, k);
        double dlt = gsl_vector_get(pD->delta, k);
//		prior -= log(dlt);  //this steps been took out to the prDM1; 
        prior -= tmp * tmp / (2.0 * dlt * dlt);
    }
    gsl_vector_free(xb);
	
    return -(xyb - sum + prior);
}

void func_df(const gsl_vector * beta, void * params, gsl_vector * df)
{
    CaseCtrlData * pD = (CaseCtrlData *) params;
    int ni = pD->ni;
    int np = pD->np;

    gsl_vector * xb = gsl_vector_alloc(ni);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pD->mx, beta, 0.0, xb);

    gsl_vector * exb = gsl_vector_alloc(ni);
    for (int i = 0; i < ni; i++)
    {
        double tt = expl(gsl_vector_get(xb, i));
        tt /= (1.0 + tt);
        gsl_vector_set(exb, i, tt);
    }   // exb is inverse logit;

    gsl_vector_memcpy(df, beta);
    gsl_vector_scale(df, -1.0);
    gsl_vector_div(df, pD->delta);
    gsl_vector_div(df, pD->delta);
    // -beta / delta^2; 
   	gsl_vector * xty = gsl_vector_alloc(np);
    gsl_blas_dgemv(CblasTrans, 1.0, pD->mx, pD->vy, 0.0, xty);
    gsl_vector_add(df, xty);
    gsl_blas_dgemv(CblasTrans, 1.0, pD->mx, exb, 0.0, xty);
    gsl_vector_sub(df, xty);
    gsl_vector_scale(df, -1.0);
    //this scale comes from we put a minus sign before likelihood. 
    gsl_vector_free(xb);
    gsl_vector_free(exb);
    gsl_vector_free(xty);
}

void func_fdf(const gsl_vector * beta, void * params, double * f, gsl_vector * df)
{
    *f = func_f(beta,  params);
    func_df(beta, params, df);
}

double Hessian(const gsl_vector * beta, void * params)
{
    CaseCtrlData * pD = (CaseCtrlData *) params;
    int ni = pD->ni;
    int np = pD->np;

    gsl_vector * xb = gsl_vector_alloc(ni);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pD->mx, beta, 0.0, xb);
    for (int i = 0; i < ni; i++)
    {
        double tmp = expl(gsl_vector_get(xb, i));
        gsl_vector_set(xb, i, tmp/ (1.0+tmp) / (1.0 + tmp));
    }

    gsl_matrix * mx = gsl_matrix_alloc(ni, np);
    gsl_matrix * wmx = gsl_matrix_alloc(ni, np);

    gsl_matrix_memcpy(mx, pD->mx);
    gsl_matrix_memcpy(wmx, pD->mx);

    for (int row = 0; row < ni; row++)
    {
        double w = gsl_vector_get(xb, row);
        for (int col = 0; col < np; col++)
        {
            double * pt = gsl_matrix_ptr(wmx, row, col);
            (*pt) *= w;
        }
    }
	// this is for fomular (5) on the wall. 
	// \sum_{i}{x_{ij}x_{ik}w_i}

    gsl_matrix * h2 = gsl_matrix_alloc(np, np);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1.0, mx, wmx, 0.0, h2);

    for (int i = 0; i < np; i++)
    {
        double * pt = gsl_matrix_ptr(h2, i, i);
        double delta = gsl_vector_get(pD->delta, i);
        (*pt) -= (1.0 / delta / delta);
    }

    gsl_permutation * perm = gsl_permutation_alloc(np);
    int signum;
    gsl_linalg_LU_decomp(h2, perm, &signum);
	double logdet = 0; 
	for (int i = 0; i < np; i++)
		logdet += log(1e-10+fabs(gsl_matrix_get(h2, i, i))); 
//    double logdet = log(fabs(gsl_linalg_LU_det(h2, signum)));
    gsl_vector_free(xb);
    gsl_matrix_free(mx);
    gsl_matrix_free(wmx);
    gsl_permutation_free(perm);
    gsl_matrix_free(h2);
    return (logdet);
}

double prDM1(CaseCtrlData * pd, gsl_vector * x)
{
    int np = pd->np;

    const gsl_multimin_fdfminimizer_type * t;
    gsl_multimin_fdfminimizer * s;
    gsl_multimin_function_fdf logit;

    logit.f = &func_f;
    logit.df = &func_df;
    logit.fdf = &func_fdf;
    logit.n = np;
    logit.params = pd;

    t = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (t, np);
    gsl_multimin_fdfminimizer_set(s, &logit, x, 0.01, 1e-2);

    int iter = 0;
	int max_iter = 100 * np; 
    int status;
    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);
        if (status)
	        break;
        status = gsl_multimin_test_gradient(s->gradient, 1e-2);
    } while (status == GSL_CONTINUE && iter < max_iter);

    double logdet = -0.5 *  Hessian(s->x, pd);

    double logmdel = 0.0;
    for (int i = 0; i < np; i++)
        logmdel += log(gsl_vector_get(pd->delta, i));
    double fval = s->f;
	for (int i = 0; i < np; i++)
		gsl_vector_set(x, i, gsl_vector_get(s->x, i)); 
//	cout << iter << "\t" << s->f << "\t" << logdet << endl; 
	//grant; 
    gsl_multimin_fdfminimizer_free(s);
    return (-fval + logdet - logmdel);
}                                       

double ModelnData::cc_bf(real sigma_a,  real sigma_d,  real * phenotype, real * genotype, int ni)
{
	if(ni == 0) return 0; 
	int ns = 1; 
	int np = m_df * ns + 1; 
    CaseCtrlData * palt = new CaseCtrlData(ni, np, phenotype);
    CaseCtrlData * pnul = new CaseCtrlData(ni, 1, phenotype);

    for (int i = 0; i < ni; i++)
    {
        gsl_matrix_set(palt->mx, i, 0, 1.0);
        gsl_matrix_set(palt->mx, i, 1, genotype[i]);
		if(m_df == 2) 
        	gsl_matrix_set(palt->mx, i, 2, ((int) genotype[i] == 1 ? 1.0 : 0.0));
    }
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pnul, palt);  
	delete palt; 
	delete pnul;
    return (logbf);
}

double ModelnData::cc_bf_mgt(real sigma_a, real sigma_d, vector<real>* phval, vector<int>* base, vector<int>* phdex, real * mgt)
{
	m_df = 1; 
	int ni = (int) phdex->size();
	int ns = 1; 
	int np = m_df * ns + 1; 
	real * ph = new real[ni]; 
	for (int i = 0; i < ni; i++)
		ph[i] = phval->at(phdex->at(i)); 
    CaseCtrlData * palt = new CaseCtrlData(ni, np, ph);
    CaseCtrlData * pnul = new CaseCtrlData(ni, 1, ph);
	
	for (int i = 0; i < ni; i++)
	{
		int j = base->at(i); 
		gsl_matrix_set(palt->mx, i, 0, 1); 
		gsl_matrix_set(palt->mx, i, 1, mgt[j]);
	}
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pnul, palt);  
	delete palt; 
	delete pnul;
	delete[] ph; 
    return (logbf);
}

double ModelnData::cc_bf_mean(real sigma_a, real sigma_d, vector<real>* phval, vector<int>* base, vector<int>* phdex, real ** prob_snp)
{
	int ni = (int) phdex->size();
	int ns = 1; 
	int np = m_df * ns + 1; 
	real * ph = new real[ni]; 
	for (int i = 0; i < ni; i++)
		ph[i] = phval->at(phdex->at(i)); 
    CaseCtrlData * palt = new CaseCtrlData(ni, np, ph);
    CaseCtrlData * pnul = new CaseCtrlData(ni, 1, ph);
	
	for (int i = 0; i < ni; i++)
	{
		int j = base->at(i); 
		gsl_matrix_set(palt->mx, i, 0, 1); 
		gsl_matrix_set(palt->mx, i, 1, prob_snp[j][1] + 2.0 * prob_snp[j][2]);
		if(m_df == 2) 
			gsl_matrix_set(palt->mx, i, 2, prob_snp[j][1]);
	}
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pnul, palt);  
	delete palt; 
	delete pnul;
	delete[] ph; 
    return (logbf);
}

double ModelnData::cc_bf_core(double sigma_a, double sigma_d, int ni, int ns, class CaseCtrlData* pNul, class CaseCtrlData* pAlt)
{
    //////////////////////////////////////////////////////////////////////////////////
	double mu = 0; 
	for (int i = 0; i < ni; i++)
		mu += gsl_vector_get(pAlt->vy, i); 
	mu /= (double) ni; 
	double scale = mu * (1.0 - mu); 
    gsl_vector * nul_start = gsl_vector_alloc(1);
    gsl_vector_set(nul_start, 0, log(mu) - log(1.0 - mu));

	for (int i = 0; i < ni; i++)
		gsl_matrix_set(pNul->mx, i, 0, 1.0);
	gsl_vector_set(pNul->delta, 0, 1000);
	double pr0 = prDM1(pNul, nul_start);
	//for null mode, calc pr(D|M0); 
	
	int col = m_df * ns + 1; 
	////////////////////////////////////////////////////////////////////////////////
	gsl_matrix * gXtX = gsl_matrix_alloc(col, col); 
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pAlt->mx, pAlt->mx, 0.0, gXtX); 
	gsl_vector * gyx = gsl_vector_alloc(col); 
	gsl_vector * alt_start = gsl_vector_alloc(col); 
	
    vector<real> vsa;   
	vector<real> vsd; 
	
	if(sigma_a == 0 && sigma_d == 0)
	{
		for (unsigned r = 0; r < vsigma_a.size(); r++)
		{
			vsa = vsigma_a; 
			vsd = vsigma_d; 
		}
	}
	else 
	{
		vsa.push_back(sigma_a);  
		vsd.push_back(sigma_d);    
	}
	int repeat = (int) vsa.size(); 
	gsl_vector * logbf = gsl_vector_alloc(repeat); 
	for (int r = 0; r < repeat; r++)
	{
		gsl_matrix_memcpy(ginvOmega, gXtX); 
		gsl_vector_set(pAlt->delta, 0, 1000);
		for(int j = 0; j < ns; j++)
		{
			gsl_vector_set(pAlt->delta, m_df*j+1, vsa.at(r));
			if(m_df == 2) 
				gsl_vector_set(pAlt->delta, m_df*j+2, vsd.at(r));
		}
		
		double inv_va = 1.0 / vsa.at(r) / vsa.at(r); 
		double inv_vd = 1.0 / vsd.at(r) / vsd.at(r); 
		
		double tmp = 0; 
		for (int s = 0; s < ns; s++)
		{
			tmp = gsl_matrix_get(ginvOmega, m_df*s+1, m_df*s+1) + inv_va; 
			gsl_matrix_set(ginvOmega, m_df*s+1, m_df*s+1, tmp);
			if(m_df == 2) 
			{
				tmp = gsl_matrix_get(ginvOmega, 2*s+2, 2*s+2) + inv_vd; 
				gsl_matrix_set(ginvOmega, 2*s+2, 2*s+2, tmp);
			}
		}   // add prior on diagnal. 

		//y^t X;
		gsl_blas_dgemv(CblasTrans, 1.0, pAlt->mx, pAlt->vy, 0.0, gyx); 
		
		gsl_linalg_cholesky_decomp(ginvOmega); 
		gsl_linalg_cholesky_solve(ginvOmega, gyx, alt_start); //res = omega %*% gyx; 
		real beta0 = gsl_vector_get(alt_start, 0); 
		gsl_vector_set(alt_start, 0, log(beta0) - log(1. - beta0)); 
		
		for (int s = 0; s < ns; s++)
		{
			gsl_vector_set(alt_start, m_df*s+1, gsl_vector_get(alt_start, m_df*s+1) / scale); 
			if(m_df == 2) 
				gsl_vector_set(alt_start, m_df*s+2, gsl_vector_get(alt_start, m_df*s+2) / scale); 
		}

		double pr1 = prDM1(pAlt, alt_start);
		// for alternative model, setup design matrix, prior, starting point, calc pr(D|M1); 
		// use mle to approximate the starting the point; 
		gsl_vector_set(logbf, r, pr1 - pr0); 
	}
	
	double maxbf = gsl_vector_max(logbf); 
	gsl_vector_add_constant(logbf, -maxbf); 
	double factor = 0; 
	for (int r = 0; r < repeat; r++)
		factor += exp(gsl_vector_get(logbf, r)); 
	double meanlogbf = maxbf + log(factor / (double) repeat); 
				
	gsl_matrix_free(gXtX);
	gsl_matrix_free(ginvOmega);
	gsl_vector_free(gyx);
	gsl_vector_free(alt_start); 
	gsl_vector_free(nul_start); 
	gsl_vector_free(logbf); 
	
	return (meanlogbf); 
}

double ModelnData::cc_bf(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, short ** genotype, class PosnBF * pt) 
{  
	int ns = pt->pos.size();
	int col = 1 + m_df * ns; 
	vector<int> ind_has_all;
	for (unsigned ind = 0; ind < vin->size(); ind++)
	{
		int i = vin->at(ind); 
		int bingle = 0; 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			if(genotype[i][s] < 0 || genotype[i][s] > 2) {
				bingle = 1; 
				break;
			}
		}
		if(!bingle)  ind_has_all.push_back(i); 
	}
	if(ind_has_all.size() == 0) return 0.0; 
	
	int ni = ind_has_all.size(); 

	real * ph = new real[ni];
	for (int i = 0; i < ni; i++)
		ph[i] = vph->at(ind_has_all.at(i)); 
	
    CaseCtrlData * pAlt = new CaseCtrlData(ni, col, ph);
    CaseCtrlData * pNul = new CaseCtrlData(ni, 1, ph);

	for (int ind = 0; ind < ni; ind++)
	{
		int i = ind_has_all.at(ind); 
		gsl_matrix_set(pAlt->mx, ind, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			gsl_matrix_set(pAlt->mx, ind, m_df*j+1, genotype[i][s]);
			if(m_df == 2) 
				gsl_matrix_set(pAlt->mx, ind, 2*j+2, (genotype[i][s] == 1 ? 1 : 0));
		}
	}

	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pNul, pAlt);  
	
	delete pAlt; 
	delete pNul;
	delete[] ph; 
    return (logbf);
}
	
double ModelnData::cc_bf(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, real ** snpInpr, class PosnBF * pt) 
{  
	int ni = vin->size(); 
	int ns = pt->pos.size();
	int np = 1 + m_df * ns; 

	real * ph = new real[ni];
	for (int i = 0; i < ni; i++)
		ph[i] = vph->at(vin->at(i)); 
	
    CaseCtrlData * pAlt = new CaseCtrlData(ni, np, ph);
    CaseCtrlData * pNul = new CaseCtrlData(ni, 1, ph);

	for (int i = 0; i < ni; i++)
	{
		gsl_matrix_set(pAlt->mx, i, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			real p0 = snpInpr[i][2 * s];
			real p1 = snpInpr[i][2*s+1];
			real p2 = 1.0 - p0 - p1; 
			if (p0 < 0 || p1 < 0 || p2 < 0) 
				cout << "-bimbam: illegal probabilities in cc_bf" << endl; 
			gsl_matrix_set(pAlt->mx, i, m_df*j+1, p1 + 2 * p2);
			if(m_df == 2) 
				gsl_matrix_set(pAlt->mx, i, 2*j+2, p1);
		}
	}
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pNul, pAlt);  

	delete pAlt; 
	delete pNul;
	delete[] ph; 
    return (logbf);
}

double ModelnData::cc_bf_mgt(real sigma_a, real sigma_d, vector<real>* vph, vector<int>* vin, real ** mgt, class PosnBF * pt) 
{  
	m_df = 1; 
	int ni = vin->size(); 
	int ns = pt->pos.size();
	int np = 1 + m_df * ns; 

	real * ph = new real[ni];
	for (int i = 0; i < ni; i++)
		ph[i] = vph->at(vin->at(i)); 
	
    CaseCtrlData * pAlt = new CaseCtrlData(ni, np, ph);
    CaseCtrlData * pNul = new CaseCtrlData(ni, 1, ph);

	for (int i = 0; i < ni; i++)
	{
		gsl_matrix_set(pAlt->mx, i, 0, 1); 
		for (int j = 0; j < ns; j++)
		{
			int s = pt->pos.at(j);
			gsl_matrix_set(pAlt->mx, i, m_df*j+1, mgt[i][s]);
		}
	}
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pNul, pAlt);  

	delete pAlt; 
	delete pNul;
	delete[] ph; 
    return (logbf);
}

double ModelnData::cc_bf(real sigma_a, real sigma_d, vector<real>& phval, vector<int>& curindex, real * genotype)
{  
	int ni = (int) curindex.size();
	int ns = 1; 
	int col = ns * m_df + 1; 
	real * ph = new real[ni]; 
	for (int i = 0; i < ni; i++)
		ph[i] = phval.at(curindex.at(i)); 
    CaseCtrlData * pAlt = new CaseCtrlData(ni, col, ph);
    CaseCtrlData * pNul = new CaseCtrlData(ni, 1, ph);
	
	for (int i = 0; i < ni; i++)                                                 
	{
		gsl_matrix_set(pAlt->mx, i, 0, 1); 
		gsl_matrix_set(pAlt->mx, i, 1, genotype[i]);
		if(m_df == 2) 
			gsl_matrix_set(pAlt->mx, i, 2, (genotype[i] == 1 ? 1 : 0));
	}
	
	double logbf = cc_bf_core(sigma_a, sigma_d, ni, ns, pNul, pAlt);  
	delete pAlt; 
	delete pNul;
	delete[] ph;                                          
    return (logbf);
}

double mylike_f(const gsl_vector * beta, void * param)
{                       
	double mu, add, sigma; 
	if(beta->size == 3)
	{
		mu = gsl_vector_get(beta, 0); 
		add = gsl_vector_get(beta, 1); 
		sigma = gsl_vector_get(beta, 2); 
	}
	else 
	{
		mu = gsl_vector_get(beta, 0); 
		add = 0; 
		sigma = gsl_vector_get(beta, 1); 
	}

	double * vpg = (double *) param; 
	int ni = (int) vpg[0]; 

	double loglike = 0; 

//	for (int i = 0; i < ni; i++)
//	{
//		double sum = 0; 
//		for (int j = 0; j < 3; j++)
//			sum += vpg[ni+1+3*i+j] * gsl_ran_gaussian_pdf(vpg[i+1] - (mu + add * j), sigma); 
//		loglike += log(sum); 
//	}   //full likelihood; 
	
	for (int i = 0; i < ni; i++)
	{
		double mg = vpg[ni+1+3*i+1] + 2.0 * vpg[ni+1+3*i+2]; 
		double	sum = gsl_ran_gaussian_pdf(vpg[i+1] - (mu + add * mg), sigma); 
		loglike += log(sum); 
	}   // likelihood based on mean genotypes; 

	return (-loglike); 
}

real ModelnData::snp_lrt(real * phval, real * genotype, int ni)
{   
	if(ni == 0) return 1; 
	double ssn = 0; 
	double mean = 0; 
	double * vpg = new double[ni * 4 + 1];
	vpg[0] = ni;  
	for (int i = 0; i < ni; i++)
	{
		vpg[i+1]  = phval[i]; 
		mean += vpg[i+1]; 
	}
	mean /= ni; 
	for (int i = 0; i < ni; i++)
	{
		double temp = (vpg[i+1] - mean); 
		ssn += temp * temp; 
	}
	ssn /= ni; 
	ssn = sqrt(ssn); 
	
	for (int i = 0; i < ni; i++)
	{
		int g=(int)genotype[i]; 
		vpg[ni+1+3*i] = 0;  
		vpg[ni+1+3*i+1] = 0;  
		vpg[ni+1+3*i+2] = 0; 
		vpg[ni+1+3*i+g] = 1; 
		
	}

///////////////////////////////////////////////////////////////////////////
//    const gsl_multimin_fdfminimizer_type * t;
//    gsl_multimin_fdfminimizer * s;
//    gsl_multimin_function_fdf minlike;
//
//    minlike.f = &mylike_f;
//    minlike.df = &mylike_df;
//    minlike.fdf = &mylike_fdf;
//    minlike.n = 3;
//    minlike.params = (void *)vpg;
//	
//	/* Starting point */
//	gsl_vector * x = gsl_vector_alloc (3); //mu, beta, sigma,
//	gsl_vector_set (x, 0, mean); //gsl_vector_get(res, 0));
//	gsl_vector_set (x, 1, 0); //gsl_vector_get(res, 1));
//	gsl_vector_set (x, 2, ssn); 
//
//    t = gsl_multimin_fdfminimizer_conjugate_fr;
//    s = gsl_multimin_fdfminimizer_alloc (t, 3);
//    gsl_multimin_fdfminimizer_set(s, &minlike, x, 0.01, 1e-4);
//
//    int iter = 0;
//	int max_iter = 100; 
//    int status;
//    do {
//        iter++;
//        status = gsl_multimin_fdfminimizer_iterate(s);
//        if (status)
//	        break;
//        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
//    } while (status == GSL_CONTINUE && iter < max_iter);
//
//	double loglike1 = s->f; 
//    double mu1 = gsl_vector_get(s->x, 0); 
//    double add1 = gsl_vector_get(s->x, 1); 
//    double sigma1 = gsl_vector_get(s->x,2); 
//    cout << endl << loglike1 << " " << mu1 << " " << add1 << " " << sigma1 << " " << iter << endl;
//    gsl_multimin_fdfminimizer_free(s);

/////////////////////////////////////////////////////////////////////////////////////////////////////

	int col = 2; 
	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
	  gsl_vector_set(gph, i, phval[i]); 

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
	  gsl_matrix_set(gXX, i, 0, 1); 
	  gsl_matrix_set(gXX, i, 1, genotype[i]);
	}
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gXX, gXX, 0.0, ginvOmega); 

	//y^t X;
	gsl_vector * gyx = gsl_vector_alloc(col); 
	gsl_blas_dgemv(CblasTrans, 1.0, gXX, gph, 0.0, gyx); 

	int sign; 
	gsl_permutation * perm = gsl_permutation_alloc(col); 
	gsl_linalg_LU_decomp(ginvOmega, perm, &sign); 
	//after LU_decomp, ginvOmega is LU; 
	real det = gsl_linalg_LU_det(ginvOmega, sign); 
	gsl_vector * res = gsl_vector_alloc(col); 
	if(det != 0) {
	  gsl_linalg_LU_solve(ginvOmega, perm, gyx, res); //res = omega %*% gyx; 
	}
	else {
	  gsl_vector_set(res, 0, mean); 
	  gsl_vector_set(res, 1, 0); 
	}

	gsl_blas_dgemv(CblasNoTrans, 1.0, gXX, res, -1.0, gph); //now gph is the residue; 
	double gsse;
	gsl_blas_ddot(gph, gph, &gsse);
    double gse = sqrt(gsse / ni); 
	
///////////////////////////////////////////////////////////////////////////////////////////
	
	double loglike1 = 0; 
	{
		gsl_vector *x; 
		x = gsl_vector_alloc (3); //mu, beta, sigma,
		gsl_vector_set (x, 0, gsl_vector_get(res, 0));
		gsl_vector_set (x, 1, gsl_vector_get(res, 1));
		gsl_vector_set (x, 2, gse); 
	   
		loglike1 = mylike_f(x,(void*)vpg); 
		gsl_vector_free(x); 
	}
		
/////////////////////////////////////////////////////////////////////////////////////////
	   
       /* Starting point */
	   gsl_vector *x; 
       x = gsl_vector_alloc (2); //mu, beta, sigma,
       gsl_vector_set (x, 0, mean);
	   gsl_vector_set (x, 1, ssn); 
	   
	   double loglike0 = mylike_f(x,(void*)vpg); 
       gsl_vector_free(x);
	   
/////////////////////////////////////////////////////////////////////////////////////////
    gsl_vector_free(gph);
	gsl_matrix_free(gXX); 
    gsl_matrix_free(ginvOmega); 
    gsl_vector_free(gyx); 
    gsl_permutation_free(perm); 
    gsl_vector_free(res); 


	delete[] vpg; 
    return (expl(-loglike1 + loglike0));
}

real ModelnData::snp_lrt(vector<real>* phval, vector<int>* base, vector<int>* curindex, real ** prob_snp)
{   
	int ni = (int) curindex->size();
	double ssn = 0; 
	double mean = 0; 
	double * vpg = new double[ni * 4 + 1];
	vpg[0] = ni;  
	for (int i = 0; i < ni; i++)
	{
		vpg[i+1]  = phval->at(curindex->at(i)); 
		mean += vpg[i+1]; 
	}
	mean /= ni; 
	for (int i = 0; i < ni; i++)
	{
		double temp = (vpg[i+1] - mean); 
		ssn += temp * temp; 
	}
	ssn /= ni; 
	ssn = sqrt(ssn); 
	
	for (int i = 0; i < ni; i++)
	{
		int j = base->at(i); 
		vpg[ni+1+3*i] =  prob_snp[j][0];  
		vpg[ni+1+3*i+1] = prob_snp[j][1];  
		vpg[ni+1+3*i+2] =  prob_snp[j][2];  
	}

///////////////////////////////////////////////////////////////////////////
//    const gsl_multimin_fdfminimizer_type * t;
//    gsl_multimin_fdfminimizer * s;
//    gsl_multimin_function_fdf minlike;
//
//    minlike.f = &mylike_f;
//    minlike.df = &mylike_df;
//    minlike.fdf = &mylike_fdf;
//    minlike.n = 3;
//    minlike.params = (void *)vpg;
//	
//	/* Starting point */
//	gsl_vector * x = gsl_vector_alloc (3); //mu, beta, sigma,
//	gsl_vector_set (x, 0, mean); //gsl_vector_get(res, 0));
//	gsl_vector_set (x, 1, 0); //gsl_vector_get(res, 1));
//	gsl_vector_set (x, 2, ssn); 
//
//    t = gsl_multimin_fdfminimizer_conjugate_fr;
//    s = gsl_multimin_fdfminimizer_alloc (t, 3);
//    gsl_multimin_fdfminimizer_set(s, &minlike, x, 0.01, 1e-4);
//
//    int iter = 0;
//	int max_iter = 100; 
//    int status;
//    do {
//        iter++;
//        status = gsl_multimin_fdfminimizer_iterate(s);
//        if (status)
//	        break;
//        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
//    } while (status == GSL_CONTINUE && iter < max_iter);
//
//	double loglike1 = s->f; 
//    double mu1 = gsl_vector_get(s->x, 0); 
//    double add1 = gsl_vector_get(s->x, 1); 
//    double sigma1 = gsl_vector_get(s->x,2); 
//    cout << endl << loglike1 << " " << mu1 << " " << add1 << " " << sigma1 << " " << iter << endl;
//    gsl_multimin_fdfminimizer_free(s);

/////////////////////////////////////////////////////////////////////////////////////////////////////

	int col = 2; 
	gsl_vector * gph = gsl_vector_alloc(ni); 
	for (int i = 0; i < ni; i++)
	  gsl_vector_set(gph, i, phval->at(curindex->at(i))); 

	gsl_matrix * gXX = gsl_matrix_alloc(ni, col); 
	for (int i = 0; i < ni; i++)
	{
	  int j = base->at(i); 
	  gsl_matrix_set(gXX, i, 0, 1); 
	  gsl_matrix_set(gXX, i, 1, prob_snp[j][1] + 2.0 * prob_snp[j][2]);
	}
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gXX, gXX, 0.0, ginvOmega); 
//	gsl_matrix_set(ginvOmega, 1, 1, gsl_matrix_get(ginvOmega, 1, 1) + 25); //1/0.2^2, to stablize it;

	//y^t X;
	gsl_vector * gyx = gsl_vector_alloc(col); 
	gsl_blas_dgemv(CblasTrans, 1.0, gXX, gph, 0.0, gyx); 

	int sign; 
	gsl_permutation * perm = gsl_permutation_alloc(col); 
	gsl_linalg_LU_decomp(ginvOmega, perm, &sign); 
	//after LU_decomp, ginvOmega is LU; 
	real det = gsl_linalg_LU_det(ginvOmega, sign); 
	gsl_vector * res = gsl_vector_alloc(col); 
	if(det != 0) {
	  gsl_linalg_LU_solve(ginvOmega, perm, gyx, res); //res = omega %*% gyx; 
	}
	else {
	  gsl_vector_set(res, 0, mean); 
	  gsl_vector_set(res, 1, 0); 
	}

	gsl_blas_dgemv(CblasNoTrans, 1.0, gXX, res, -1.0, gph); //now gph is the residue; 
	double gsse;
	gsl_blas_ddot(gph, gph, &gsse);
    double gse = sqrt(gsse / ni); 
	
///////////////////////////////////////////////////////////////////////////////////////////
	
	double loglike1 = 0; 
	unsigned max_iter = 1;
	if(max_iter == 1)
	{
		gsl_vector *x; 
		x = gsl_vector_alloc (3); //mu, beta, sigma,
		gsl_vector_set (x, 0, gsl_vector_get(res, 0));
		gsl_vector_set (x, 1, gsl_vector_get(res, 1));
		gsl_vector_set (x, 2, gse); 
	   
	    loglike1 = mylike_f(x,(void*)vpg); 
	    gsl_vector_free(x); 
	}
    else
	{
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		gsl_multimin_function minex_func;

		size_t iter = 0;
		int status;
		double size;

		/* Starting point */
		x = gsl_vector_alloc (3); //mu, beta, sigma,
		gsl_vector_set (x, 0, gsl_vector_get(res, 0));
		gsl_vector_set (x, 1, gsl_vector_get(res, 1));
		gsl_vector_set (x, 2, gse); 

		// cout << endl << "starting pts: " << gsl_vector_get(x,0) << "\t" << gsl_vector_get(x,1) << "\t";
		//cout << gsl_vector_get(x,2);
		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc (3);
		gsl_vector_set_all (ss, .2);

		/* Initialize method and iterate */
		minex_func.n = 3;
		minex_func.f = &mylike_f;
		minex_func.params = (void *)vpg;

		s = gsl_multimin_fminimizer_alloc (T, 3);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

		do {
		   iter++;
		   status = gsl_multimin_fminimizer_iterate(s);
		   
		   if (status) 
			 break;

		   size = gsl_multimin_fminimizer_size (s);
		   status = gsl_multimin_test_size (size, 1e-2);
//		double mu1 = gsl_vector_get(s->x, 0); 
//		double add1 = gsl_vector_get(s->x, 1); 
//		double sigma1 = gsl_vector_get(s->x,2); 
//		cout << endl << s->fval << " " << mu1 << " " << add1 << " " << sigma1 << " " << iter ;

		} while (status == GSL_CONTINUE && iter < max_iter);

		loglike1 = s->fval; 

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	}

		
/////////////////////////////////////////////////////////////////////////////////////////
	   
       /* Starting point */
	   gsl_vector *x; 
       x = gsl_vector_alloc (2); //mu, beta, sigma,
       gsl_vector_set (x, 0, mean);
	   gsl_vector_set (x, 1, ssn); 
	   
	   double loglike0 = mylike_f(x,(void*)vpg); 
       gsl_vector_free(x);
//	   cout << endl << loglike0 << endl; 
	   
/////////////////////////////////////////////////////////////////////////////////////////
    gsl_vector_free(gph);
	gsl_matrix_free(gXX); 
    gsl_matrix_free(ginvOmega); 
    gsl_vector_free(gyx); 
    gsl_permutation_free(perm); 
    gsl_vector_free(res); 


	delete[] vpg; 
    return (expl(-loglike1 + loglike0));
}

//use -sem str prefix.str.em
//    -rem prefix.str.em 
void ModelnData::write_em_param(void)  
{
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".em"); 
	fstream outfile;
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write:" << sfn << endl;
		return;
	}

    outfile << nEMRuns <<  " #number of EM runs" << endl;
    outfile << nLoci <<  " #number of loci" << endl;
    outfile << nK <<  " #number of clusters" << endl;
    outfile << nSubPop << " #number of subpopulations" << endl;

	for (int runs = 0; runs < nEMRuns; runs++)
	{
		outfile << "###theta" << endl; 
		real ** theta = pMP[runs].Gettheta();
		for (int k = 0; k < nK; k++)
		{
			char str[100]; 
			for (int m = 0; m < nLoci; m++)
			{
				sprintf(str, " %.3f",theta[m][k]);
				outfile << str; 
			}
			outfile << endl; 
		}
		outfile << "###alpha" << endl; 
		for (int p = 0; p < nSubPop; p++)
		{
			real ** alpha = pMP[runs].Getalpha(p); 
			for (int k = 0; k < nK; k++)
			{
				char str[100]; 
				for (int m = 0; m < nLoci; m++)
				{
					sprintf(str, " %.3f",alpha[m][k]);
					outfile << str; 
				}
				outfile << endl; 
			}
			outfile << endl; 
		}
		outfile << "###r" << endl; 
		for (int p = 0; p < nSubPop; p++)
		{
			char str[100]; 
			real * r = pMP[runs].Getr(p); 
			for (int m = 0; m < nLoci; m++)
			{
				sprintf(str, " %.6f", r[m]);
				outfile << str; 
			}
			outfile << endl; 
		}
	}
	
    outfile.close();
}

int ModelnData::read_em_param(string sfn)   
{
	real * par = NULL; 
#if defined (MPI_ENABLED)
	int * pconst = NULL; 
	if(procID == 0)
#endif
	{
		fstream infile;
		streambuf * pbuf;
		infile.open(sfn.c_str(), ios::in);
		if(!infile.is_open()) 
		{
			cout << "-bimbam: cannot open file to read:" << sfn << endl;
			return 0;
		}
		pbuf = infile.rdbuf(); 	
		char delimit[] = ",; \t:"; 
		string line; 
		line.assign(getline(pbuf));
		nEMRuns = atoi(line.data()); 
		line.assign(getline(pbuf));
		int loci = atoi(line.data()); 
		if(loci != nLoci) 
		{
			cout << "-bimbam: em parameters not consistent" << endl; 
			return 0;
		}
		line.assign(getline(pbuf));
		nK = atoi(line.data()); 
		line.assign(getline(pbuf));
		nSubPop = atoi(line.data()); 

		if (pMP) delete[] pMP; 
		pMP = new ModelParam[nEMRuns]; 
		InitModelParam(); 

		int bufsize = nEMRuns * (nLoci * nK + nLoci * nK * nSubPop + nLoci * nSubPop);
		par = (real*) Allocate1D(sizeof(real), bufsize);
		int p = 0; 
		for (int runs = 0; runs < nEMRuns; runs++)
		{
			for (int k = 0; k < nK; k++)
			{
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int m = 0; m < nLoci; m++)
				{
					if(res == NULL) 
					{
						cout << "-bimbam: bad file in read em" << endl; 
						safe_exit(); 
					}
					par[p] = atof(res);
					p++; 
					res = strtok(NULL, delimit); 
				}
			}
			
			for (int s = 0; s < nSubPop; s++)
				for (int k = 0; k < nK; k++)
				{
					line.assign(getline(pbuf)); 
					char * res = strtok((char*)line.c_str(), delimit); 
					for (int m = 0; m < nLoci; m++)
					{
						if(res == NULL) 
						{
							cout << "-bimbam: bad file in read em" << endl; 
							safe_exit(); 
						}
						par[p] = atof(res);
						p++; 
						res = strtok(NULL, delimit); 
					}
				}
			
			for (int s = 0; s < nSubPop; s++)
			{
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int m = 0; m < nLoci; m++)
				{
					if(res == NULL) 
					{
						cout << "-bimbam: bad file in read em" << endl; 
						safe_exit(); 
					}
					par[p] = atof(res);
					p++; 
					res = strtok(NULL, delimit); 
				}
			}
		}
	}
	
#if defined (MPI_ENABLED)
	pconst = new int[4]; 
	if (procID == 0) 
	{
		pconst[0] = nEMRuns;
	    pconst[1] = nLoci; 
		pconst[2] = nK; 
		pconst[3] = nSubPop; 
	} 
	MPI_Bcast(pconst, 4, MPI_INT, 0, MPI_COMM_WORLD); 
	if (procID > 0) 
	{
		nEMRuns = pconst[0];
		nLoci = pconst[1];
		nK = pconst[2];
		nSubPop = pconst[3]; 
		if(pMP) delete[] pMP; 
		pMP = new ModelParam[nEMRuns]; 
		InitModelParam(); 
	}
	
//	for (int p = 0; p < nProc; p++)
//	{
//		if(procID == p)
//			cout << nEMRuns << "\t" << nLoci << "\t" << nK << "\t" << nSubPop << endl; 
//		MPI_Barrier(MPI_COMM_WORLD); 
//	}

	delete[] pconst; 
	int bufsize = nEMRuns * (nLoci * nK + nLoci * nK * nSubPop + nLoci * nSubPop);
	if (procID > 0) 
		par = (real *) Allocate1D(sizeof(real), bufsize); 
	if(sizeof(real) == sizeof(double))
		MPI_Bcast(par, bufsize, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	else
		MPI_Bcast(par, bufsize, MPI_FLOAT, 0, MPI_COMM_WORLD); 
#endif
	
	int p = 0; 
	for (int runs = 0; runs < nEMRuns; runs++)
	{
		for (int k = 0; k < nK; k++)
			for (int m = 0; m < nLoci; m++)
				pMP[runs].Settheta(m, k, par[p++]); 
		for (int s = 0; s < nSubPop; s++)
			for (int k = 0; k < nK; k++)
				for (int m = 0; m < nLoci; m++)
					pMP[runs].Setalpha(s, m, k, par[p++]); 
		for (int s = 0; s < nSubPop; s++)
			for (int m = 0; m < nLoci; m++)
				pMP[runs].Setr(s, m, par[p++]); 
	}
	free(par); 
	return 1; 
}

real max(real a, real b)
{
	return (a > b ? a:b); 
}

real ModelnData::calc_bf_of_snp_summary(real sigma_a, real sigma_d, class SNPSummaryData * ssd)
{
	int ni = ssd->ni; 
	if (ni == 0) return 1.0;  
	int ns = 1;
	int col = m_df + 1; 
	gsl_matrix * gXtX = gsl_matrix_alloc(col, col); 
	gsl_matrix * ginvOmega = gsl_matrix_alloc(col, col); 
	gsl_matrix_set(gXtX, 0, 0, ssd->ni); 
	gsl_matrix_set(gXtX, 0, 1, ssd->sg);
	gsl_matrix_set(gXtX, 1, 0, ssd->sg);
	gsl_matrix_set(gXtX, 1, 1, ssd->sg2);
	if(m_df == 2) 
	{
		gsl_matrix_set(gXtX, 0, 2, ssd->sd);
		gsl_matrix_set(gXtX, 2, 0, ssd->sd);
		gsl_matrix_set(gXtX, 2, 2, ssd->sd2);
		gsl_matrix_set(gXtX, 1, 2, ssd->sgd);
		gsl_matrix_set(gXtX, 2, 1, ssd->sgd);
	}
	
	//y^t X;
	gsl_vector * gyx = gsl_vector_alloc(col); 
	gsl_vector_set(gyx, 0, ssd->sy);
	gsl_vector_set(gyx, 1, ssd->syg);
	if(m_df == 2) 
		gsl_vector_set(gyx, 2, ssd->syd);  
	double yy = ssd->sy2; 

	gsl_vector * res = gsl_vector_alloc(col); 
	
    vector<double> vsa;   
	vector<double> vsd; 
	
	if(sigma_a < 1e-6 && sigma_d < 1e-6)
	{
		for (unsigned r = 0; r < vsigma_a.size(); r++)
		{
			vsa.push_back(1.0 / vsigma_a.at(r) / vsigma_a.at(r)); 
			vsd.push_back(1.0 / vsigma_d.at(r) / vsigma_d.at(r)); 
		}
	}
	else 
	{
		vsa.push_back(1.0 / sigma_a / sigma_a);  // (1.0 / 0.2 /0.2); 
		vsd.push_back(1.0 / sigma_d / sigma_d); // (1.0 / 0.05 / 0.05);    
	}
	
	int repeat = (int) vsa.size(); 
	gsl_vector * bf = gsl_vector_alloc(repeat); 
	m_beta.clear(); 
	for (int i = 0; i < col; i++)
		m_beta.push_back(0); 
	for (unsigned p = 0; p < vsa.size(); p++)
	{
		gsl_matrix_memcpy(ginvOmega, gXtX); 
		double inv_va = vsa.at(p); 
		double inv_vd = vsd.at(p); 
		for (int i = 0; i < ns; i++)
		{
			real tmp = 0; 
			if(m_df == 1)
			{
				tmp = gsl_matrix_get(ginvOmega, i+1, i+1) + inv_va; 
				gsl_matrix_set(ginvOmega, i+1, i+1, tmp);
			}
			if(m_df == 2) 
			{
				tmp = gsl_matrix_get(ginvOmega, 2*i+1, 2*i+1) + inv_va; 
				gsl_matrix_set(ginvOmega, 2*i+1, 2*i+1, tmp);
				tmp = gsl_matrix_get(ginvOmega, 2*i+2, 2*i+2) + inv_vd; 
				gsl_matrix_set(ginvOmega, 2*i+2, 2*i+2, tmp);
			}
		}   // add prior on diagnal. 

//		cout << inv_va << "\t" << inv_vd << endl; 
//		for (int i = 0; i < col; i++)
//		{
//			for (int j = 0; j < col; j++)
//			{
//				cout << gsl_matrix_get(ginvOmega, i, j) << " "; 
//			}
//			cout << endl; 
//		}
//		cout << endl; 
//		//grant; 
		
		gsl_linalg_cholesky_decomp(ginvOmega); 
		double logdet = 0; 
		for (int i = 0; i < col; i++)
			logdet += log(fabs(gsl_matrix_get(ginvOmega, i, i))); 
		logdet *= 2.0; 
		
//		cout << "logdet: " << logdet << endl; 
		

		gsl_linalg_cholesky_solve(ginvOmega, gyx, res); //res = omega %*% gyx; 
														  
		double bob = 0;
		gsl_blas_ddot(res, gyx, &bob); 
		
		double tau = sqrt((yy - bob) / ni); 
		for (int i = 0; i < col; i++)
			m_beta.at(i) += (gsl_vector_get(res, i) / tau); 

		double ph_2 = gsl_vector_get(gyx, 0);
		ph_2 *= ph_2; 
		
		if(yy - ph_2/ni < 1e-100) 
			gsl_vector_set(bf, p, 0); 
		else 
		{
			double tlast = (yy - bob) / (yy - ph_2/ni);
	//		cout << "tlast: " << tlast << endl; 
												  
			double logvavdn = ns * log((double) inv_va) + log((double)ni);    
			if (m_df == 2) 
				logvavdn += ns * log((double) inv_vd); 
			double logbf = -0.5 * logdet + 0.5 * logvavdn - ni * 0.5 * log((double)tlast); 
			gsl_vector_set(bf, p, logbf); 
//			cout << "lobgf: " << logbf << endl; 
		}
	}

	for (int i = 0; i < col; i++)
		m_beta.at(i) /= repeat; 
	double maxbf = gsl_vector_max(bf); 
	gsl_vector_add_constant(bf, -maxbf); 
	double factor = 0; 
	for (int i = 0; i < repeat; i++)
		factor += exp(gsl_vector_get(bf, i)); 
	factor /= (double) repeat; 
	double logbf = maxbf + log(factor); 
	
	gsl_matrix_free(gXtX);
	gsl_matrix_free(ginvOmega);
	gsl_vector_free(gyx);
	gsl_vector_free(res);
	gsl_vector_free(bf);
	return (logbf);  
}

void ModelnData::combine_snp_summary_data(class SNPSummaryData * ssda, class SNPSummaryData * ssdb)
{
    pair<char, char> p1, p2; 
	p1.first = ssda->major; 
	p1.second = ssda->minor; 

	p2.first = ssdb->major; 
	p2.second = ssdb->minor; 
	
	int action = compatibleQ(p1, p2); 
	if(action == -1) 
	{
		cout << "-bimbam: alleles are not compatiable between two studies: " << ssda->rs << endl; 
		fplog << "-bimbam: alleles are not compatiable between two studies: " << ssda->rs << endl; 
		return; 
	} 

	else if(action == 0 ) 
	{
		ssdb->ni += ssda->ni; 
		ssdb->sg += ssda->sg; 
		ssdb->sg2 += ssda->sg2; 
		ssdb->sgd += ssda->sgd; 
		ssdb->sd += ssda->sd; 
		ssdb->sd2 += ssda->sd2; 
		ssdb->sy += ssda->sy; 
		ssdb->sy2 += ssda->sy2; 
		ssdb->syg += ssda->syg; 
		ssdb->syd += ssda->syd; 
	}
	else if(action == 1) 
	{   //flip the ssda allele coding; 
		cout << "-bimbam: flip allele coding in summary data: " << ssda->rs << endl; 
		fplog << "-bimbam: flip allele coding in summary data: " << ssda->rs << endl; 
		double tni = ssda->ni; 
		double tsg = 2 * tni - ssda->sg; 
		double tsg2 = ssda->sg2 + 4.0 * tni - 4.0 * ssda->sg; 
		double tsgd = 2.0 * ssda->sd - ssda->sgd; 
		double tsyg = 2.0 * ssda->sy - ssda->syg;
		
		ssdb->ni += ssda->ni; 
		ssdb->sg += tsg; 
		ssdb->sg2 += tsg2; 
		ssdb->sgd += tsgd; 
		ssdb->sd += ssda->sd; 
		ssdb->sd2 += ssda->sd2; 
		ssdb->sy += ssda->sy; 
		ssdb->sy2 += ssda->sy2; 
		ssdb->syg += tsyg; 
		ssdb->syd += ssda->syd; 
	}
}

void ModelnData::ProduceSummaryData(int gmode, string fn)
{
	fstream outfile; 
	string sfn("output/");
	if(fn.size() == 0)
	{
		sfn.append(fnOutput);
		sfn.append(".snp.summary.data.ssd");
	}
	else
		sfn.append(fn); 
	
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write: " << sfn << endl;
		return;
	}
	
	vector<real> *phval = &vv_phval.at(0); 
	vector<int> *phindex = &vv_phdex.at(0); 

	outfile << "SNP A B STRAND ni sg sg2 sgd sd sd2 sy sy2 syg syd" << endl; 
	for (int m = 0; m < nLoci; m++)
	{
		SNPSummaryData st; 
		string rs(vsRsnum.at(m)); 
//		get_snp_summary_data(m, &st); 
		st.ni = 0; 
		st.sg = 0; 
		st.sg2 = 0;
		st.sgd = 0; 
		st.sd = 0; 
		st.sd2 = 0; 
		st.sy = 0;
		st.sy2 = 0;
		st.syg = 0;
		st.syd = 0; 
		
		pair<char, char> tp = mapRs2mm[rs]; 
		st.major = tp.first; 
		st.minor = tp.second; 
		
		for (unsigned i = 0; i < phindex->size(); i++)
		{
			int ind = phindex->at(i); 
			double gt = -1; 
			real ph = phval->at(ind); 
			real dt = 0;   
			if(gmode == 0) 
			{
				short temp =  pIndiv[ind]->GetsnpGT(m);
				if(temp == QQ) continue; 
				dt = (temp == 1 ? 1 : 0);
				gt = (double) temp;
			}
			else if(gmode == 1) 
				gt = (double) pIndiv[ind]->get_snpmgt(m);
			else 
			{
				real res[2]; 
				pIndiv[ind]->get_snp_dstr(m, res); 
				dt = res[1]; 
				gt = res[1] + (1-res[0]-res[1]) * 2.0; 
			}
			
			st.ni ++; 
			st.sg += gt; 
			st.sg2 += (gt * gt); 
			st.sgd += (gt * dt); 
			st.sd += dt; 
			st.sd2 += (dt * dt); 

			st.sy += ph; 
			st.sy2 += (ph * ph); 
			st.syg += (ph * gt); 
			st.syd += (ph * dt);  
		}
		char buf[500];   
		if(m_allele_coding_mode == 0) 
			sprintf(buf, "%.15s %c %c NA ", rs.c_str(), st.minor, st.major);  
		else 
			sprintf(buf, "%.15s %c %c NA ", rs.c_str(), st.major, st.minor);  
		outfile << buf; 

		sprintf(buf, "%d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n", 
				st.ni, st.sg, st.sg2, st.sgd, st.sd, st.sd2, st.sy, st.sy2, st.syg, st.syd); 
		outfile << buf; 
	} 
	outfile.close(); 
}

void ModelnData::CombineStudy()
{
	map<string, class SNPSummaryData > :: iterator siter; 
	vector<string> vSNPrs; 
	for (unsigned nf = 0; nf < vSSD.size(); nf++)
	{
		ifstream infile; 
		char delimit[] = ";, \t";
		streambuf * pbuf;
		//first, read in the rs pos file; and put into vRsPos;
		infile.open(vSSD.at(nf).c_str(), ios::in); 
		if(!infile.is_open()) 
		{
			cout << "-bimbam: cannot open snp summary data to read" << endl;
			safe_exit(); 
		}
		pbuf = infile.rdbuf(); 	
		string line(getline(pbuf)); //the first line; 
			
		while(pbuf->sgetc() != EOF) 
		{
			class SNPSummaryData objssd; 
			string line(getline(pbuf));
			char * res = strtok((char*)line.c_str(), delimit); 
			objssd.rs.assign(res); 
			
			res = strtok(NULL, delimit); 
			objssd.major = res[0]; 
			
			res = strtok(NULL, delimit); 
			objssd.minor = res[0]; 

			res = strtok(NULL, delimit); 
			//strand NOT in use; 
			
			res = strtok(NULL, delimit); 
			objssd.ni = atoi(res);
			res = strtok(NULL, delimit); 
			objssd.sg = atof(res);
			res = strtok(NULL, delimit); 
			objssd.sg2 = atof(res);
			res = strtok(NULL, delimit); 
			objssd.sgd = atof(res);
			res = strtok(NULL, delimit); 
			objssd.sd = atof(res);
			res = strtok(NULL, delimit); 
			objssd.sd2 = atof(res);
			
			res = strtok(NULL, delimit); 
			objssd.sy = atof(res);
			res = strtok(NULL, delimit); 
			objssd.sy2 = atof(res);
			res = strtok(NULL, delimit); 
			objssd.syg = atof(res);
			res = strtok(NULL, delimit); 
			objssd.syd = atof(res);

			siter = mapRs2SNPSummary.find(objssd.rs); 
			if(siter == mapRs2SNPSummary.end())
			{
				mapRs2SNPSummary[objssd.rs] = objssd; 
				vSNPrs.push_back(objssd.rs); 
			}
			else 
				combine_snp_summary_data(&objssd, &mapRs2SNPSummary[objssd.rs]); 
		}
		infile.close();
	}

	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".combined.study.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		fplog << "ERROR: cannot open file to write: " << sfn << endl;
		cout << "-bimbam: cannot open file to write: " << sfn << endl;
		safe_exit(); 
	}
	
	outfile << " rs  log10bf" << endl; 
	for (unsigned m = 0; m < vSNPrs.size(); m++)
	{
		string rs(vSNPrs.at(m)); 
		real bfs = 0; 
		bfs = calc_bf_of_snp_summary(0, 0, &mapRs2SNPSummary[rs]); 
		char buf[100]; 
		sprintf(buf, "%.15s \t\t %+5.3f\n", rs.c_str(), bfs/log(10.0)); 
		outfile << buf; 
	} 
	outfile.close(); 
}

#if defined (POWER)
void ModelnData::read_is_write()
{
	stringstream lss; 
	char delimit[] = ";, \t";
	streambuf * pbuf;
	ifstream infile; 
	//first, read in the rs pos file; and put into vRsPos;
	infile.open(fnFILE.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		lss << ":< Can't open genotype distribution file, proceed with care." << endl; 
		fplog << lss.str(); 
		cout << lss.str(); 
		lss.str("\0"); 
	}
	
	pbuf = infile.rdbuf(); 	

	string line; 
	int ni = 0; 
	{
		line.assign(getline(pbuf)); 
		char * res = strtok(line.c_str(), delimit); 
		ni = atoi(res);
	}
	int ns = 0; 
	{
		line.assign(getline(pbuf)); 
		char * res = strtok(line.c_str(), delimit); 
		ns = atoi(res);
	}	//number of snps;
	
	vector<string> rs; 
	real ** gd = Allocate2DMatrix(ns, ni * 3); 
	for (int m = 0; m < ns; m++) 
	{
		string line(getline(pbuf));
		string temp; 
		char * res = strtok(line.c_str(), delimit); 
		temp.assign(res); 
		rs.push_back(temp); 

		for (int i = 0; i < ni; i++)
		{
			res = strtok(NULL, delimit); 
			gd[m][i*3+0] = atof(res);
			res = strtok(NULL, delimit); 
			gd[m][i*3+1] = atof(res);
			res = strtok(NULL, delimit); 
			gd[m][i*3+2] = atof(res);
			
		}
	}
	infile.close(); 
	//read in genotype distribution 
	
	cout << "ni = " << ni << endl;
	cout << "ns = " << ns << endl; 
	cout << "nPH = " << nPH << endl; 
	
	vector<vector<real> > trans_phval; 
	string sfn; 
		   
	sfn.assign(vPin.at(0)); 
	infile.open(sfn.c_str(), ios::in);
	if(!infile.is_open())
	{
		cout << "Error: can't open file in read phenotype data!" << endl; 
		vGin.clear();
		vPin.clear(); 
		return; 
	}
	pbuf = infile.rdbuf();
	
	for (int i = 0; i < ni; i++)
	{
		vector<real> ind_phval; 
		line.assign(getline(pbuf)); 
		char * res = strtok(line.c_str(), delimit); 
		for (int s = 0; s < nWarmSteps; s++)     //nWarmSteps as shift;
		{
			res = strtok(NULL, delimit); 
		}
			
		for (int np = 0; np < nPH; np++)
		{
			if(res == NULL) break; 
			string sv(res);
			if(sv.compare("NA") == 0)
				ind_phval.push_back(NA);
			else                       	
				ind_phval.push_back(atof(sv.data()));
			res = strtok(NULL, delimit); 
		}

		if ((int) ind_phval.size() == 0) {
			cout << "Error: number of phenotypes in " << sfn << " is smaller than number of individuals..." << endl;
			vGin.clear(); 
			vPin.clear(); 
			return; 
		}
			
		if ((int) ind_phval.size() < nPH) {
			cout << "Warning: use a smaller number of phenotypes" << endl; 
			nPH = (int)ind_phval.size(); 
		}
		trans_phval.push_back(ind_phval);
		ind_phval.clear();
	}
	line.assign(getline(pbuf));
	if(line.size() > 0) {
			cout << "Error: number of phenotypes in " << sfn << " is larger than number of individuals..." << endl;
			vGin.clear(); 
			vPin.clear(); 
			return; 
	}
	infile.close(); 
	
	
	cout << trans_phval.size() << endl; 
	cout << trans_phval.at(0).size() << endl; 
	
	for(int np = 0; np < nPH; np++)
	{
		vector<int> vindex; 
		vector<real> vt; 
		for(int i = 0; i < ni; i++)
		{
			if (fabs(trans_phval.at(i).at(np) - NA) < 0.1)
			{  
				vt.push_back(trans_phval.at(i).at(np));
			}
			else 
			{
				vindex.push_back(i); 
				vt.push_back(trans_phval.at(i).at(np));
			}
		}
		vv_phval.push_back(vt);
		vv_phdex.push_back(vindex);
		vindex.clear(); 
		vt.clear();
	}
	vector<vector<real> >().swap(trans_phval); 
	//trans_phval.clear();	
	for (int i = 0; i < ni; i++)
		cout << vv_phval.at(0).at(i) << " "; 
	cout << endl; 
	for (int i = 0; i < vv_phdex.at(0).size(); i++)
		cout << vv_phdex.at(0).at(i) << " "; 
	cout << endl; 
	//read in phenotype; 
 	
	int nCohort = ni; 
	int nLoci = ns; 
	{
		cout << "num Phenotype = " << nPH << endl; 
		cout << "num Cohort = " << nCohort << endl; 
	}

	real ** bf = Allocate2DMatrix(nPH, nLoci*2);   //mean and variance. 
	real * phenoVal = new real[nCohort]; 
	real * one_snp = new real[nCohort];
	real ** snppr = Allocate2DMatrix(nCohort, 3);  
	gsl_matrix * XtX = gsl_matrix_alloc(3, 3); 
	
	real ** q = Allocate2DMatrix(nCohort, 3);
	gsl_vector * xy = gsl_vector_alloc(3); 
    gsl_vector * beta = gsl_vector_alloc(3); 
	
	char str[100] = "Importance Sampling and Single-SNP BFs: ";

	int prg_cur = 0; 
	int prg_all = nLoci;
	int prg_mod = (int)(prg_all / 100.0);
	if(prg_mod == 0) prg_mod = 1 ;
	
	gsl_vector * logbf = gsl_vector_alloc(nImpute); 
	for (int m = 0; m < nLoci; m++)
	{
		for (int i = 0; i < nCohort; i++)
		{
			snppr[i][0] = gd[m][3*i+0];
			snppr[i][1] = gd[m][3*i+1];
			snppr[i][2] = gd[m][3*i+2];
			if(snppr[i][0] < 0 || snppr[i][0] < 0 || snppr[i][0] < 0)
				cout << "illegal probabilities in importance sampling. " << endl; 
		}

		for (int np = 0; np < nPH; np++)
		{
			vector<real>* phval = &vv_phval.at(np); 
			vector<int>* phdex = &vv_phdex.at(np);   
			int nCohort = (int) phdex->size();

			gsl_matrix * X = gsl_matrix_alloc(nCohort, 3); 
			gsl_vector * gph = gsl_vector_alloc(nCohort); 
			
			for (int i = 0; i < nCohort; i++)
			{
				int j = phdex->at(i);
				phenoVal[i] = phval->at(j);
				gsl_vector_set(gph, i, phval->at(j)); 
				gsl_matrix_set(X, i, 0, 1.0); 
				gsl_matrix_set(X, i, 1, snppr[j][1] + 2.0 * snppr[j][2]);
				gsl_matrix_set(X, i, 2, snppr[j][1]); 
			}
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XtX); 
			gsl_blas_dgemv(CblasTrans, 1.0, X, gph, 0.0, xy);  
			
			real var = 0.2 * 0.2; 
			gsl_matrix_set(XtX, 1, 1, gsl_matrix_get(XtX, 1, 1) + 1.0 / var);
			gsl_matrix_set(XtX, 2, 2, gsl_matrix_get(XtX, 2, 2) + 16.0 / var);

			gsl_linalg_cholesky_decomp(XtX); 
			gsl_linalg_cholesky_solve(XtX, xy, beta); 
			
			for (int i = 0; i < nCohort; i++)
			{        
				double b0 = gsl_vector_get(beta, 0); 
				double b1 = gsl_vector_get(beta, 1); 
				double b2 = gsl_vector_get(beta, 2); 
				double ph = gsl_vector_get(gph, i); 
				
				int j = phdex->at(i); 
				q[i][0] = snppr[j][0] * gsl_ran_gaussian_pdf(ph - b0, 1.0);
				q[i][1] = snppr[j][1] * gsl_ran_gaussian_pdf(ph - (b0 + b1 + b2), 1.0);
				q[i][2] = snppr[j][2] * gsl_ran_gaussian_pdf(ph - (b0 + b1 + b1), 1.0);   	
				//generate q;
				real sum = q[i][0] + q[i][1] + q[i][2];
				q[i][0] /= sum; q[i][1] /= sum; q[i][2] /= sum; 
			}   //normalization;       		   
			
//			real tbf = 0.0; 
//			if(cc) tbf = cc_bf_mean(0, 0, phval, phdex, phdex, snppr); 
//			else tbf = calc_bf_mean(0, 0, phval, phdex, phdex, snppr); 
//			bf[np][2 * m] = tbf; 
//			bf[np][2*m+1] = 0;
//			//mean genotype bf. 
			
			//now generate samples based on q.
			for (int repeat = 0; repeat < nImpute; repeat++)
			{    
				double logpwght = 0;
				double logqwght = 0;                       
				double ratio = 1;    
				if (gsl_rng_uniform(gsl_r) < 0.9)
				{ 
					for (int i = 0; i < nCohort; i++)
					{
						int j = phdex->at(i); 
						real unif = gsl_rng_uniform(gsl_r); 
						if(unif < q[i][0])
						{
							one_snp[i] = 0; 
							logpwght += log(snppr[j][0]);
							logqwght += log(q[i][0]);
						}
						else if(unif < (q[i][0] + q[i][1]))
						{
							one_snp[i] = 1; 
							logpwght += log(snppr[j][1]);
							logqwght += log(q[i][1]);    						
						}
						else 
						{
							one_snp[i] = 2;
							logpwght += log(snppr[j][2]);
							logqwght += log(q[i][2]);
						}
					}
					ratio = exp(logpwght - logqwght);    
				}
				else
				{   
					for (int i = 0; i < nCohort; i++)
					{
						int j = phdex->at(i); 
						real unif = gsl_rng_uniform(gsl_r); 
						if(unif < snppr[j][0])
							one_snp[i] = 0; 
						else if(unif < (snppr[j][0] + snppr[j][1]))
							one_snp[i] = 1; 
						else 
							one_snp[i] = 2;
					}
				}
				
				real tbf = 0; 
				if(cc) tbf = cc_bf(0, 0, phenoVal, one_snp, nCohort);
				else tbf = calc_bf(0, 0, phenoVal, one_snp, nCohort);
				tbf += log(ratio);        
				gsl_vector_set(logbf, repeat, tbf); 
			}                                
			double maxbf = gsl_vector_max(logbf); 
			gsl_vector_add_constant(logbf, -maxbf); 
			double f1 = 0; 
			double f2 = 0; 
			for (int r = 0; r < nImpute; r++)
			{
				double temp = exp(gsl_vector_get(logbf, r)); 
				f1 += temp; 
				f2 += temp * temp; 
			}
			f1 /= (double) nImpute; 
			f2 /= (double) nImpute; 
			bf[np][2 * m] = maxbf + log(f1);        //log first-moment; 
			bf[np][2*m+1] = 2.0 * maxbf + log(f2);      //log second-moment; 
			gsl_vector_free(gph); 
			gsl_matrix_free(x); 
		} //end pheno. 
		prg_cur++; 
		if(prg_cur % prg_mod == 0) print_progress_bar(0, str, prg_cur, prg_all);
	} // end for(m=0;m<nLoci;m++)
	print_progress_bar(1, str, 1, 1);
		
	delete[] phenoVal;
	delete[] one_snp;
	Free2DMatrix(snppr);
	Free2DMatrix(q);
	gsl_vector_free(xy); 
	gsl_vector_free(beta); 
	gsl_matrix_free(XtX); 
	//importance sampling; 
	
	sfn.assign("output/");
	sfn.append(fnOutput);
	char buf[20];
	sprintf(buf, ".shift%d", nWarmSteps); 
	sfn.append(buf); 
	sprintf(buf, ".nph%d", nPH); 
	sfn.append(buf); 

	sfn.append(".gd.mean.01.txt");
	fstream outfile;
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) {
		cout << "Error: cannot open file to write:" << sfn << endl;
		return;
	}
	
	real total = (real) nImpute; 
	if(total <= 0 ) total = 1; 	
	cout << "nImpute = " << total << endl; 
	outfile << "## bf = log10_BF\t se = log10_std_err\t" << endl; 
	for (int m = 0; m < nLoci; m++)
	{
		char buf[1000];
		sprintf(buf, "%-15s\t", rs.at(m).data()); 
		outfile << buf; 
		for (int p = 0; p < nPH; p++)
		{  
			real t1, t2; 
			t1 = bf[p][2 * m];
			t2 = bf[p][2*m+1];
			if(total >= 2) 
			{
				double tmax = max(t1 * 2.0, t2); 
				double diff1 = t1 * 2.0 - tmax; 
				double diff2 = t2 - tmax; 
				t2 = 0.5 * (tmax + log(exp(diff2) - exp(diff1)) - log(total)); 
			}
			t1 /= log(10.0); 
			t2 /= log(10.0); 
			if(total >= 2) 
				sprintf(buf, "%+5.3f %+5.3f\t", t1, t2);
			else
				sprintf(buf, "%+5.3f  NA   \t", t1);
			outfile << buf;
		}
		outfile << endl; 
	}
	outfile.close(); 
	Free2DMatrix(bf);  
	//output
}  

void ModelnData::read_gd_stat_write()
{
	stringstream lss; 
	char delimit[] = ";, \t";
	streambuf * pbuf;
	ifstream infile; 
	
	infile.open(fnFILE.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		lss << ":< Can't open genotype distribution file, proceed with care." << endl; 
		fplog << lss.str(); 
		cout << lss.str(); 
		lss.str("\0"); 
	}
	
	pbuf = infile.rdbuf(); 	

	string line; 
	int ni = 0; 
	ni = nWarmSteps; 
	int ns = 0; 
	ns = nMaxSteps; 
	
	vector<string> rs; 
	real ** gd = Allocate2DMatrix(ns, ni * 3); 
	for (int m = 0; m < ns; m++) 
	{
		string line(getline(pbuf));
		char * res = strtok(line.c_str(), delimit); 
		rs.push_back(string temp(res)); 

		for (int i = 0; i < ni; i++)
		{
			res = strtok(NULL, delimit); 
			gd[m][i*3+0] = atof(res);
			res = strtok(NULL, delimit); 
			gd[m][i*3+1] = atof(res);
			gd[m][i*3+2] = 1.0 - gd[m][i*3] - gd[m][i*3+1];
		}
	}
	infile.close(); 
	//read in genotype distribution 
	
	cout << "ni = " << ni << endl;
	cout << "ns = " << ns << endl; 
	cout << "nPH = " << nPH << endl; 
	for(int m = 0; m < 1; m++)
	{
		for (int i = 0; i < 10; i++)
		{
			cout << gd[m][3*i] << " "<< gd[m][3*i+1] << " " <<  gd[m][3*i+2]  << " \t"; 
		}
	    cout << endl; 
	}

	string sfn; 
	
	vector<vector<real> > trans_phval; 
		   
	sfn.assign(fnDOC); 
	infile.open(sfn.c_str(), ios::in);
	if(!infile.is_open())
	{
		cout << "Error: can't open file in read phenotype data!" << endl; 
		vGin.clear();
		vPin.clear(); 
		return; 
	}
	pbuf = infile.rdbuf();
	
	for (int i = 0; i < ni; i++)
	{
		vector<real> ind_phval; 
		line.assign(getline(pbuf)); 
		char * res = strtok(line.c_str(), delimit); 	
		for (int np = 0; np < nPH; np++)
		{
			if(res == NULL) break;
			string sv(res);
			if(sv.compare("NA") == 0)
				ind_phval.push_back(NA);
			else                       	
				ind_phval.push_back(atof(sv.data()));
			res = strtok(NULL, delimit); 
		}

		if ((int) ind_phval.size() == 0) {
			cout << "Error: number of phenotypes in " << sfn << " is smaller than number of individuals..." << endl;
			vGin.clear(); 
			vPin.clear(); 
			return; 
		}
			
		if ((int) ind_phval.size() < nPH) {
			cout << "Warning: use a smaller number of phenotypes" << endl; 
			nPH = (int)ind_phval.size(); 
		}
		trans_phval.push_back(ind_phval);
		ind_phval.clear();
	}
	line.assign(getline(pbuf));
	if(line.size() > 0) {
			cout << "Error: number of phenotypes in " << sfn << " is larger than number of individuals..." << endl;
			vGin.clear(); 
			vPin.clear(); 
			return; 
	}
	infile.close(); 
	
	
	cout << trans_phval.size() << endl; 
	cout << trans_phval.at(0).size() << endl; 
	
	for(int np = 0; np < nPH; np++)
	{
		vector<int> vindex; 
		vector<real> vt; 
		for(int i = 0; i < ni; i++)
		{
			if (fabs(trans_phval.at(i).at(np) - NA) < 0.1)
			{  
				vt.push_back(trans_phval.at(i).at(np));
			}
			else 
			{
				vindex.push_back(i); 
				vt.push_back(trans_phval.at(i).at(np));
			}
		}
		vv_phval.push_back(vt);
		vv_phdex.push_back(vindex);
		vindex.clear(); 
		vt.clear();
	}
	vector<vector<real> >().swap(trans_phval); 
	//trans_phval.clear();	
//	for (int i = 0; i < ni; i++)
//		cout << vv_phval.at(0).at(i) << " "; 
//	cout << endl; 
//	for (unsigned i = 0; i < vv_phdex.at(0).size(); i++)
//		cout << vv_phdex.at(0).at(i) << " "; 
//	cout << endl; 
	//note: vv_phval contains all phenotype. 
	//    : vv_phdex contains "index" of legal phenotype that ignore NA when counting. 
	//read in phenotype; 
 	
	int nCohort = ni; 
	int nLoci = ns; 
	{
		cout << "num Phenotype = " << nPH << endl; 
		cout << "num Cohort = " << nCohort << endl; 
	}

	real ** bf = Allocate2DMatrix(nPH, nLoci*2);   //mean and variance. 
	real * phenoVal = new real[nCohort]; 
	real * one_snp = new real[nCohort];
	real ** snppr = Allocate2DMatrix(nCohort, 3);  
	
	char str[100] = "Calc Statistics Using Genotype Distributions : ";

	int prg_cur = 0; 
	int prg_all = nLoci;
	int prg_mod = (int)(prg_all / 100.0);
	if(prg_mod == 0) prg_mod = 1 ;
	
	time_t sec_beg, sec_end; 
	sec_beg = time (NULL);
	for (int m = 0; m < nLoci; m++)
	{
		if(m < (int)vsRsnum.size()) {
			map<string, real> :: iterator iter; 
			iter = mapRs2maf.find(vsRsnum.at(m)); 
			if(iter != mapRs2maf.end())
				m_current_maf = iter->second; 
			else
				m_current_maf = 0.5; 
		} 
//		cout << m_current_maf << endl; 
		for (int i = 0; i < nCohort; i++)
		{
			snppr[i][0] = gd[m][3*i+0];
			snppr[i][1] = gd[m][3*i+1];
			snppr[i][2] = gd[m][3*i+2];
			if(snppr[i][0] < 0 || snppr[i][0] < 0 || snppr[i][0] < 0)
				cout << "illegal probabilities in importance sampling. " << endl; 
		}

		for (int np = 0; np < nPH; np++)
		{
			vector<real>* phval = &vv_phval.at(np); 
			vector<int>* phdex = &vv_phdex.at(np);   
			int nCohort = (int) phdex->size();

			for (int i = 0; i < nCohort; i++)
			{
				int j = phdex->at(i);
				phenoVal[i] = phval->at(j);
			}
			
			real tbf = 0.0; 
			if(nImpute == 1)
			{
				for (unsigned t = 0; t < vsigma_a.size(); t++)
			   		tbf += calc_bf_mean(vsigma_a.at(t), vsigma_d.at(t), phval, phdex, snppr); 
				tbf /= vsigma_a.size(); 
			}   //default prior to calc bfs. use -A -D whne using other priors. 
			else if(nImpute == 33)
				tbf = snp_lrt(phval, phdex, phdex, snppr); 
			
			bf[np][2 * m] = tbf; 
			bf[np][2*m+1] = 0;
			//mean genotype bf. 
		} //end pheno. 
		prg_cur++; 
		if(prg_cur % prg_mod == 0) print_progress_bar(0, str, prg_cur, prg_all);
	} // end for(m=0;m<nLoci;m++)
	print_progress_bar(1, str, 1, 1);
		
	 sec_end = time(NULL); 
	 cout << "## seconds used =" <<  sec_end - sec_beg << endl; 
	 fplog << "## seconds used =" <<  sec_end - sec_beg << endl; 
	delete[] phenoVal;
	delete[] one_snp;
	Free2DMatrix(snppr);
	
	sfn.assign("output/");
	sfn.append(fnOutput);
	sfn.append(".stat.txt");
	
	fstream outfile;
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) {
		cout << "Error: cannot open file to write:" << sfn << endl;
		return;
	}
	
	for (int m = 0; m < ns; m++)
	{
		char buf[1000];
		for (int p = 0; p < nPH; p++)
		{  
			real t1 = log10(bf[p][2 * m]);
			sprintf(buf, "%+5.3f ", t1);
			outfile << buf;
		}
		outfile << endl; 
	}
	outfile.close(); 
	Free2DMatrix(bf);  
	//output
}  
#endif


cGENE::cGENE()
{
	;
}

cGENE::~cGENE()
{
	;
}

cGENE::cGENE(string sgn, int schr, long sbeg, long send) 
{
	 gn.assign(sgn); 
	 chr = schr; 
	 beg = sbeg; 
	 end = send; 
}	

void cGENE::assign(class cGENE ct)
{
	gn.assign(ct.get_name()); 
	chr = ct.get_chr(); 
	beg = ct.get_beg(); 
	end = ct.get_end(); 
}

cSNP::cSNP(string srs, int schr, long spos, long sindex, pair<char, char> smm, real smaf)
{
	rs.assign(srs); 
	chr = schr; 
	pos = spos; 
	index = sindex; 
	major = smm.first; 
	minor = smm.second; 
	maf = smaf; 
	flag = 0; 
}                                                                         

#if defined  (SPECIAL)
void ModelnData::fix_mean_genotype()
{
	fstream outfile; 
	string sfn(vPin.at(0));
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write:" << sfn << endl;
		return;
	}
	
	ifstream infile; 
	string tfn(vGin.at(0)); 
	infile.open(tfn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		cout << " can't open file to read ..." << tfn << endl; 
		safe_exit(); 
	}

	streambuf * pbuf;
	pbuf = infile.rdbuf(); 	
	string line; 
	line.assign(getline(pbuf));
	while (!line.empty()) 
	{
		for (int i = 0; i < 100; i++)
		{
			if((int) line.at(i) == 0)
				line.at(i) = 'N'; 
		}
		outfile << line << endl; 
		line.assign(getline(pbuf)); 
	}

	outfile.close(); 
	infile.close(); 
	
}

void ModelnData::write_for_gain()
{
	read_bimbam_genotype_distribution(2, 0, 0, -1); 
	
	int wchr = m_num; 
	
	ifstream infile; 
	char delimit[] = ";, \t";
	streambuf * pbuf;

	map<string, int> rs2i; 
    map<string, int> :: iterator iter; 
	
	//first, read in the rs pos file; and put into vRsPos;
	string tfn("/Users/yguan/GAIN/snp2impute119.txt"); 
	infile.open(tfn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		cout << " can't open file to read ..." << tfn << endl; 
		safe_exit(); 
	}

	pbuf = infile.rdbuf(); 	
	string line; 
	line.assign(getline(pbuf));
	cout << "read in " << tfn << endl; 
	while (!line.empty()) 
	{
		line.append(" 0 "); 
		char * res = strtok((char *)line.c_str(), delimit); 
		int chr = atoi(res); 
		if(chr != wchr) 
		{
			line.assign(getline(pbuf)); 
			continue; 
		}
		if(chr > wchr) break; 
		res = strtok(NULL, delimit); 
		string rs(res); 
		rs2i[rs] = 1; 
		line.assign(getline(pbuf)); 
	}
	cout << rs2i.size() << endl; 
	//get the vec and harsh for snp2impute; 	
	infile.close(); 	
	tfn.assign("/Users/yguan/GAIN/fam.ind.id"); 
	infile.open(tfn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		cout << " can't open file to read ..." << tfn << endl; 
		safe_exit(); 
	}

	vector<string> vid; 
	pbuf = infile.rdbuf(); 	
	line.assign(getline(pbuf)); 
	while (!line.empty()) 
	{
		line.append(" 0 "); 
		char * res = strtok((char *)line.c_str(), delimit); 
		string id(res); 
	    vid.push_back(id); 
		line.assign(getline(pbuf)); 
	}
    //get the fam and indiv id, they are the same in this case.
	infile.close(); 
	cout << vid.size() << endl; 
	fstream outfile; 
	string sfn("/Users/yguan/GAIN/ggain.cimputed.chr");
	char buf[100]; 
	sprintf(buf, "%d", wchr); 
	sfn.append(buf); 
	sfn.append(".txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write:" << sfn << endl;
		return;
	}
	
	outfile << "FAMID INDID "; 
	for (int m = 0; m < nLoci; m++)
	{
		string rs(vsRsnum.at(m)); 
		iter = rs2i.find(rs); 
		if(iter == rs2i.end()) continue; 
		map<string, pair<char, char> > :: iterator miter; 
		miter = mapRs2mm.find(rs); 
		if(miter == mapRs2mm.end()) 
		{
			cout << " missing snp " << rs << endl; 
		}
		char a = miter->second.first; 
		char b = miter->second.second; 
		if(a == 'N') a = 'N'; 
		if(b == 'N') b = 'N'; 
		if(a < b) 
			outfile << rs << "." << a << a << " " << rs << "." << a << b << " "; 
		else 
			outfile << rs << "." << b << b << " " << rs << "." << b << a << " "; 
	}
	outfile << endl; 
	//build the header line; 

	for (int i = 0; i < nCohort; i++)
	{
		real * gd = pIndiv[i]->get_snp_dstr(); 
		outfile << vid.at(i) << " " << vid.at(i) << " "; 
		for (int m = 0; m < nLoci; m++)
		{
			string rs(vsRsnum.at(m)); 
			iter = rs2i.find(rs); 
			if(iter == rs2i.end()) continue; 
			map<string, pair<char, char> > :: iterator miter; 
			miter = mapRs2mm.find(rs); 
			
			char a = miter->second.first; 
			char b = miter->second.second; 
			char buf[100]; 
			
			if(a < b) 
				sprintf(buf, " %.3f %.3f ", gd[2*m], gd[2*m+1]); 
			else 
			{
				double temp = 1.0 - gd[2*m] - gd[2*m+1];
				if(temp < 0) temp = 0.0; 
				sprintf(buf, " %.3f %.3f ", temp, gd[2*m+1]); 
			}
			outfile << buf; 
		}
		outfile << endl; 
	}
	outfile.close(); 
}

void ModelnData::for_barbara()
{
	char delimit[] = ";, \t";
	streambuf * pbuf;
	ifstream infile;    

	int ni = 1000; 
	int ns = 3533; 
	int nf = 12; 
	string line; 
	
	string fgeno("chr22.geno"); 
	string ffeature("chr22.feature"); 
	string fph("threepheno.phe"); 
	
	infile.open(fgeno.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
	    cout << ":< Can't open file..." << endl; 
		exit(0); 
	}
	
	pbuf = infile.rdbuf(); 	
    real ** gt = Allocate2DMatrix(ns, ni); 
	
	for (int m = 0; m < ns; m++) 
	{
		string line(getline(pbuf));
		char * res = strtok((char*)line.c_str(), delimit); 

		for (int i = 0; i < ni; i++)
		{
			if(res != NULL); 
			gt[m][i] = atof(res);
			res = strtok(NULL, delimit); 
		}
	}
	infile.close(); 
	//read in genotype;  
	
	infile.open(ffeature.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
	    cout << ":< Can't open file..." << endl; 
		exit(0); 
	}
	
	pbuf = infile.rdbuf(); 	
    real ** ft = Allocate2DMatrix(ns, nf); 
	
	for (int m = 0; m < ns; m++) 
	{
		string line(getline(pbuf));
		char * res = strtok((char*)line.c_str(), delimit); 

		for (int i = 0; i < nf; i++)
		{
			if(res != NULL); 
			ft[m][i] = atof(res);
			res = strtok(NULL, delimit); 
		}
	}
	infile.close(); 
	//read in feature;  


	infile.open(fph.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
	    cout << ":< Can't open file..." << endl; 
		exit(0); 
	}
	
	pbuf = infile.rdbuf(); 	
    real ** ph = Allocate2DMatrix(ni, 3); 
	
	for (int i = 0; i < ni; i++) 
	{
		string line(getline(pbuf));
		char * res = strtok((char*)line.c_str(), delimit); 

		for (int j = 0; j < 3; j++)
		{
			if(res != NULL); 
			ph[i][j] = atof(res);
			res = strtok(NULL, delimit); 
		}     
	}
	infile.close(); 
	//read in pheno;  
	

	double curlike = 0; 
	double nextlike = 0; 
	int which = m_num; //which phenotype;
	real * gval = new real[ni]; 
	real * phval = new real[ni]; 
    double *wgts = new double[nf+1]; 
	double *next_wgts = new double[nf+1]; 
	double *sig = new double[ns]; 
    real * bf = new real[ns]; 
	for (int i = 0; i < nf+1; i++)
		wgts[i] = gsl_rng_uniform(gsl_r); 
	for (int m = 0; m < ns; m++)
	{
		sig[m] = wgts[0]; 
		for (int i = 0; i < nf; i++)
			sig[m] += wgts[i+1] * ft[m][i]; 
		if(sig[m] < 0) sig[m] = -sig[m]; 
	}
	
	m_df = 1; 
    double maxbf = -1e100; 
	for (int m = 0; m < ns; m++)
	{
		for (int i = 0; i < ni; i++)
		{
			gval[i] = gt[m][i]; 
			phval[i] = ph[i][which]; 
		}
		
		bf[m] = calc_bf(sig[m], 0.00001, gval, phval, ni);
		if (bf[m] > maxbf) maxbf=bf[m]; 
	}
	for (int m = 0; m < ns; m++)
		bf[m] -= maxbf; 
	double sum = 0; 
	for (int m = 0; m < ns; m++)
		sum += exp(bf[m]); 
	curlike = maxbf + log(sum); 
	
	string sfn; 
	sfn.append(fnOutput);
	sfn.append(".txt");
	
	fstream outfile;
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "Error: cannot open file to write:" << sfn << endl;
		return;
	}
	
	for (int s = 0; s < nMaxSteps; s++)
	{
		for (int i = 0; i < nf+1; i++)
		{
			next_wgts[i] = wgts[i] + (gsl_rng_uniform(gsl_r) - 0.5)*0.2;                                      
			if(next_wgts[i] < 0) next_wgts[i] = -next_wgts[i]; 
			else if(next_wgts[i] > 1) next_wgts[i] = 2-next_wgts[i]; 
		}
		for (int m = 0; m < ns; m++)
		{
			sig[m] = next_wgts[0]; 
			for (int i = 0; i < nf; i++)
				sig[m] += next_wgts[i+1] * ft[m][i]; 
//			if(sig[m] < 0) sig[m] = -sig[m]; 
		}
		//propose new wgts; 

		{
			double maxbf = -1e100; 
			for (int m = 0; m < ns; m++)
			{
				for (int i = 0; i < ni; i++)
				{
					gval[i] = gt[m][i]; 
					phval[i] = ph[i][which]; 
				}
				
				bf[m] = calc_bf(sig[m], 0.00001, gval, phval, ni);
				if (bf[m] > maxbf) maxbf=bf[m]; 
			}
			for (int m = 0; m < ns; m++)
				bf[m] -= maxbf; 
			double sum = 0; 
			for (int m = 0; m < ns; m++)
				sum += exp(bf[m]); 
			nextlike = maxbf + log(sum); 
		}

		double alpha = nextlike - curlike; 
//		alpha *= log(1.0 + s); 
		if(alpha > 100) alpha = 1;
		else if(alpha < -100) alpha = 0; 
		else alpha = exp(alpha); 

		if(gsl_rng_uniform(gsl_r) < alpha) 
		{
			curlike = nextlike; 
			for (int i = 0; i < nf+1; i++)
				wgts[i] = next_wgts[i]; 
		}

		if(s % 10 == 0) 
		{
		    char temp[100]; 
			sprintf(temp, "%3.2f \t", curlike); 
			outfile << temp; 
			
			for (int i = 0; i < nf+1; i++)
			{
				sprintf(temp, "%3.2f ", wgts[i]); 
				outfile << temp; 
			}
			outfile << endl; 
		}
	}
	
	outfile.close(); 

	delete[] gval; 
	delete[] phval; 
	delete[] wgts; 
	delete[] next_wgts; 
	delete[] sig; 
	delete[] bf; 

	Free2DMatrix(gt); 
	Free2DMatrix(ft); 
	Free2DMatrix(ph); 
	
}  

//only for scan; 
void ModelnData::single_snp_cluster(void)
{
	vector<real> phval = vv_phval.at(0); 
	vector<int> phindex = vv_phdex.at(0); 
	int nCohort = (int) phindex.size(); 
	
	real * phenoVal = new real[nCohort];
	for (int i = 0; i < nCohort; i++)
		phenoVal[i] = vv_phval.at(0).at(phindex.at(i));
	
	short ** snpIn01 = Allocate2DShortMatrix(nLoci, nCohort * 2);
	
	real * bf1 = new real[nLoci];
	real * bf2 = new real[nLoci];
	for (int m = 0; m < nLoci; m++)
		bf1[m] = bf2[m] = 0.0; 
	real * onesnp = new real[nCohort];
	int unitrepeat = (int)((real) nImpute / nEMRuns);
	nImpute = unitrepeat * nEMRuns; 
	for (int runs = 0; runs < nEMRuns; runs++)
	{
		int ** zpair = NULL; 
		for (int repeat = 0; repeat < unitrepeat; repeat++)
		{
			for (int i = 0; i < nCohort; i++)
			{
				int ind = phindex.at(i); 
				pIndiv[ind]->CalcAll(nLoci, nK, pMP+runs, 0);
				pIndiv[ind]->joint_imputation(pMP, runs, NULL, nLoci, nK, nLoci, NULL);
				zpair = pIndiv[ind]->Getzpair(); 

				for (int m = 0; m < nLoci; m++)
				{
					snpIn01[m][i * 2] = zpair[m][0];
					snpIn01[m][i*2+1] = zpair[m][1];
				}
				pIndiv[ind]->FreeMemAfterEM();
			}
		
			for (int m = 0; m < nLoci; m++)
			{
				double bf_val = 0.0;
				for (int k = 0; k < nK; k++)
				{
					for (int i = 0; i < nCohort; i++)
					{
						int tmp = 0; 
						if(snpIn01[m][i * 2] == k) tmp++; 
						if(snpIn01[m][i*2+1] == k) tmp++; 
						onesnp[i] = tmp; 
					}
					
					if(cc) bf_val += cc_bf(0, 0, phenoVal, onesnp, nCohort);
					else bf_val += calc_bf(0, 0, phenoVal, onesnp, nCohort);
				}
				bf_val /= (double) nK; 
				bf1[m] += bf_val;
				bf2[m] += bf_val * bf_val; 
			}
		}
	} //end runs

	for (int m = 0; m < nLoci; m++)
	{
		bf1[m] /= (real) nImpute; 
	    bf2[m] /= (real) nImpute;
		bf2[m] -= bf1[m] * bf1[m];
		bf2[m] /= (real) nImpute; 
		bf1[m] = log10(bf1[m]);
		bf2[m] = 0.5 * log10(bf2[m]); 
	}                           
	
	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".cluster.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "-bimbam: cannot open file to write:" << sfn << endl;
		return;
	}
	for (int m = 0; m < nLoci; m++)
		outfile <<  bf1[m] << "\t" << bf2[m] << endl; 

	outfile.close(); 
	
   	Free2DShortMatrix(snpIn01);  
	delete[] onesnp; 
	delete[] phenoVal;  
	delete[] bf1;  
	delete[] bf2;  
}    
#endif 

//private:
//    int nn;   //n; 
//    int np;   //p; 
//    int nd;   //d; 
//    //same notation as in paper. 
//
//    double ** gt;    //nxp; 
//    double ** ph;    //nxd; 
//    //data; 
//
//    int pm;
//    double ** psi;   //dxd; 
//    double sa;  //sigma_a;   K=diag(0, sa, sa, ...); 
//    //prior; 
//
//public:
//
//   logbf_summaries(gsl_matrix * m1phi, gsl_matrix * m1syy, gsl_vector * v1s11, int * zz, int pm, double sa, double * lbf)
//   logbf_rank1(double sa, double pi0, int pm, int * config, double * lbf); 
//    single_snp_multi_phenotype(int mode); 

void ModelnData::single_snp_multi_phenotype(int mode)
{


  //--Heejung start --//
  if(mode > 0){

    if(vsigma_a.size() == 0){
      fplog << "## BIMBAM: Use default priors on additive effects" << endl; 
      vsigma_a.clear(); 
      vsigma_a.push_back(0.05); 
      vsigma_a.push_back(0.1); 
      vsigma_a.push_back(0.2);   
      vsigma_a.push_back(0.4);
    }else{
      fplog << "## BIMBAM: Use user specified priors on additive effects" << endl; 
    }
    
    if(vsigma_d.size() > 0)
      vsigma_d.clear(); 

    //if(vsigma_d.size() == 0){
    //  fplog << "## BIMBAM: Use default value for a prior probability of the null hypothesis." << endl; 
    //  vsigma_d.clear(); 
    //  vsigma_d.push_back(0.5); 
    //}else if(vsigma_d.size() > 1){
    //  fplog << "## BIMBAM: Use user specified value for a prior probability of the null hypothesis." << endl; 
    //}   
   
  }

  //-- Heejung start --//
  //cout << "Heejung : " << vsigma_a.size() << endl;
  //for (unsigned t = 0; t < vsigma_a.size(); t++)
  //  cout << "t " << t << " and  " << vsigma_a.at(t) << endl;
    //tbf += calc_bf_mean(vsigma_a.at(t), vsigma_d.at(t), phval, phdex, snppr); 
  //-- Heejung end --//

  //	double sa = vsigma_a.at(0);  
 //double pi0 = vsigma_d.at(0);
  double sa;
  double pi0 = 0;
  


    //cout << "pi0 : " << pi0 << endl;

    nn = nCohort; 
	np = nLoci; 
	nd = nPH; 
//	cout << " nn, np, nd " << nn << " " << np << " " << nd << endl; 
	
	gt = Allocate2DMatrix(nCohort, nLoci); 
	int ni = 0; 
	for (int i = 0; i < nIndiv; i++)
	  {
		if(pIndiv[i]->GetisPanel()) continue; 
		for (int m = 0; m < nLoci; m++)
		{
			real tt = pIndiv[i]->get_snpmgt(m); 
			gt[ni][m] = tt; 
		}
		ni++; 
	}

	ph = Allocate2DMatrix(nCohort, nPH); 
	ni = 0; 
	for (int i = 0; i < nIndiv; i++)
	{
//		vv_phdex; 
//		vv_phval; 
		if(pIndiv[i]->GetisPanel()) continue; 
		for (int p = 0; p < nPH; p++)
		{
			ph[ni][p] = vv_phval.at(p).at(i); 
		}
		ni++;
	}   //note: here we assume no missing value in both gt and ph; 
	
    //center columns of gt, ph;  only need to be done once in after read in data; 
	center_col(gt, nn, np);
	center_col(ph, nn, nd);


	//-- Heejung start : erase all parts --//

	/*	

	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".mph.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}
	
	if(mode == 1)    //Bfall; 
	{


	  double * lbf = new double[np]; 

	  int * config = new int[nd]; 
	  for (int d = 0; d < nd; d++)
	    {
	      config[d] = 1; 
	    }
	  logbf_rank1(sa, pi0, 0, config, lbf); 
	  outfile << "rs " ; 
	  for (int d = 0; d < nd; d++)
	    outfile << config[d]; 
	  outfile << endl; 
	  for (int m = 0; m < np; m++)
	    {
	      outfile << vsRsnum.at(m) << " "; 
	      char buf[100]; 
	      if(lbf[m] < 1e-3) 
		sprintf(buf, "%.3f ", lbf[m]); 
	      else 
		sprintf(buf, "%+.3f ", lbf[m]); 
	      outfile << buf << endl; 
	    }
	  delete[] lbf; 
	}
	
	else if(mode == 2)  //enumerate all combination of phenotypes; 
	  {

	    int total = (int) pow(3, nd);
        int* zz = new int[nd];
	double * lbf = new double[np]; 
	double * lpr = new double[total]; 
	real ** mbf = Allocate2DMatrix(total, np); 
	
        for (int i = 0; i < total; i++)
	  {
            for (int j = 0; j < nd; j++)
	      {
                int aa = i % (int)pow(3,(j+1));
                zz[j] = aa / (int)pow(3,j);
	      }
	    
            double prior = compute_prior(zz, pi0);
	    lpr[i] = prior; 
            if(prior > 1e-6)
	      logbf_rank1(sa, pi0, 0, zz, lbf); 
            else
	      {
                for (int m = 0; m < np; m++)
		  lbf[m] = 0.0; 
	      }
	      
	    for (int m = 0; m < np; m++)
	      mbf[i][m] = lbf[m]; 
	  }

	*/
	//-- Heejung end --//


	//-- Heejung start : to average over all different sigma_a
	if(mode == 1){    //Bfall; 
       

	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".mph.BFs.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}
	



	  int lenSa = vsigma_a.size();

	  vector<vector<double> > BFs(0);
	  double maxlogBF;
	  BFs.resize(nd + 1);

	  double * lbf = new double[np]; 
	  int * config = new int[nd]; 
	 
	  for(int j = 0; j <= nd; j++){
	  
	    BFs[j].resize(np);

	    if(j <nd){
	      for (int d = 0; d < nd; d++)
		config[d] = 2; 
	      config[j] = 1;
	    }else{
	      for (int d = 0; d < nd; d++)
		config[d] = 1; 
	    }

	    for(int i =0; i < lenSa; i++){

	      sa = vsigma_a.at(i);
	      logbf_rank1(sa, pi0, 0, config, lbf);

	      if(i == 0){
		for(int d = 0; d < np; d++)
		  BFs[j][d] = lbf[d];
	      }else{
		for(int d = 0; d < np; d++){
		  maxlogBF = max(BFs[j][d], lbf[d]);
		  BFs[j][d] = log(exp(BFs[j][d] - maxlogBF) + exp(lbf[d] - maxlogBF)) + maxlogBF;
		}
	      }
	    }

	    for(int d = 0; d < np; d++)
	      BFs[j][d] = BFs[j][d] - log((double)lenSa);

	  }

	  delete[] lbf; 
	  //outfile << "rs " ; 
	  //for (int d = 0; d < nd; d++)
	  //  outfile << config[d]; 
	  //outfile << endl; 
	  for (int m = 0; m < np; m++){
	    outfile << vsRsnum.at(m) << " "; 
	    char buf[100]; 

	    outfile << "NA ";
	    
	    if((BFs[nd][m]/log(10)) < 1e-5) 
	      sprintf(buf, "%.5f ", (BFs[nd][m]/log(10))); 
	    else 
	      sprintf(buf, "%+.5f ", (BFs[nd][m]/log(10))); 
	    outfile << buf; 	    

	    for(int j =0; j < nd; j++){
	      if((BFs[j][m]/log(10)) < 1e-5) 
		sprintf(buf, "%.5f ", (BFs[j][m]/log(10))); 
	      else 
		sprintf(buf, "%+.5f ", (BFs[j][m]/log(10))); 
	      outfile << buf;
	    }
	    outfile << endl;
	  }

	  outfile.close(); 
	  cout << sfn << " has been created." << endl;
	 

	  for(int j =0; j <= nd; j++)
	    BFs[j].resize(0);
	  BFs.resize(0);
	  /*

	  int lenSa = vsigma_a.size();

	  vector<double> BFs(0);
	  double maxlogBF;
	  BFs.resize(np);

	  double * lbf = new double[np]; 
	  int * config = new int[nd]; 

	  for (int d = 0; d < nd; d++){
	    config[d] = 1; 
	  }

	  for(int i =0; i < lenSa; i++){

	    sa = vsigma_a.at(i);
	    logbf_rank1(sa, pi0, 0, config, lbf);

	    if(i == 0){
	      for(int d = 0; d < np; d++)
		BFs[d] = lbf[d];
	    }else{
	      for(int d = 0; d < np; d++){
		maxlogBF = max(BFs[d], lbf[d]);
		BFs[d] = log(exp(BFs[d] - maxlogBF) + exp(lbf[d] - maxlogBF)) + maxlogBF;
	      }
	    }
	  }

	  for(int d = 0; d < np; d++)
	    lbf[d] = BFs[d] - log((double)lenSa);

	  BFs.resize(0);

	  outfile << "rs " ; 
	  for (int d = 0; d < nd; d++)
	    outfile << config[d]; 
	  outfile << endl; 
	  for (int m = 0; m < np; m++){
	    outfile << vsRsnum.at(m) << " "; 
	    char buf[100]; 
	    if(lbf[m] < 1e-3) 
	      sprintf(buf, "%.3f ", lbf[m]); 
	    else 
	      sprintf(buf, "%+.3f ", lbf[m]); 
	    outfile << buf << endl; 
	  }
	  delete[] lbf; 

	  */


	}else if(mode == 2){  //enumerate all combination of phenotypes; 

	  int total = (int) pow(3, nd);
	  int* zz = new int[nd];
	  double * lbf = new double[np]; 
	  double * lpr = new double[total]; 
	  real ** mbf = Allocate2DMatrix(total, np); 
	  
	  int lenSa = vsigma_a.size();

	  //-- numD start --//
	  vector<int> numD(0);
	  vector<int> numA(0);
	  numD.resize(total);
	  numA.resize(total);
	  int Dcount = 0;
	  int Icount = 0;
	  vector<vector<double> > distD(0);
	  vector<vector<double> > distA(0);
	  distD.resize(np);
	  distA.resize(np);
	  //-- numD end --//

	  vector<double> BFs(0);
	  double maxlogBF;
	  BFs.resize(np);


	  for (int i = 0; i < total; i++){

	    //-- numD start --//
	    Dcount = 0;
	    Icount = 0;
	    //-- numD end --//
            for (int j = 0; j < nd; j++){
	      int aa = i % (int)pow(3,(j+1));
	      zz[j] = aa / (int)pow(3,j);
	      //-- numD start --//
	      if(zz[j] == 1){
		Dcount++;
	      }else if(zz[j] == 2){
		Icount++;
	      }
	      //-- numD end --//
	    }
	    
	    //-- numD start --//
	    numD[i] = Dcount;
	    numA[i] = Dcount + Icount;
	    //-- numD end --//

            double prior = compute_prior(zz, pi0);
	    lpr[i] = prior; 
            if(prior > 1e-6){

	      for(int k =0; k < lenSa; k++){

		sa = vsigma_a.at(k);
		logbf_rank1(sa, pi0, 0, zz, lbf);

		if(k == 0){
		  for(int d = 0; d < np; d++)
		    BFs[d] = lbf[d];
		}else{
		  for(int d = 0; d < np; d++){
		    maxlogBF = max(BFs[d], lbf[d]);
		    BFs[d] = log(exp(BFs[d] - maxlogBF) + exp(lbf[d] - maxlogBF)) + maxlogBF;
		  }
		}
	      }

	      for(int d = 0; d < np; d++)
	        BFs[d] -= log((double)lenSa);

            }else{
	      for (int d = 0; d < np; d++)
		BFs[d] = 0.0; 
	    }
	    
	    for (int d = 0; d < np; d++)
	      mbf[i][d] = BFs[d]; 

	  
	  }

	  BFs.resize(0);
	 
       


	  //-- Heejung End --//
	

	//--- Heejung Start --//
	// To calculate average BF and posterior prob
	
	vector<double> avBFs(np);
	vector<vector<double> > lPP(np);
       
	for(int i = 0; i < np; i++){
	  lPP[i].resize(total);
	  avBFs[i] = 0;
        }




	// calculating log of averaged BFs (avBFs) and not scaled log of posterior prob (lPP)
	for (int i = 1; i < total; i++){

	  if(lpr[i] > 1e-6){

	    for(int j = 0; j < np; j++){
	      //lPP[j][i] = mbf[i][j] + log(lpr[i]/pi0);
	      lPP[j][i] = mbf[i][j] + log(lpr[i]);
	      maxlogBF = max(avBFs[j], lPP[j][i]);
	      if(i > 1)
		avBFs[j] = log(exp(avBFs[j] - maxlogBF) + exp(lPP[j][i] - maxlogBF)) + maxlogBF;
	      else
		avBFs[j] = lPP[j][i];
	    }

	  }else{

	    for(int j = 0; j < np; j++)
	      lPP[j][i] = 0.0;

          }
	}


	// should memory of mbf fee here!!
	// scale log of posterior prob and calculate PP for U, I, D
	vector<vector<vector<double> > > ppUID(np);
	for(int i = 0; i < np; i++){
	  ppUID[i].resize(nd);
	  //-- numD start --//
	  distD[i].resize(nd);
	  distA[i].resize(nd);
	  //-- numD end --//
 	  for(int j = 0; j < nd; j++){
	    ppUID[i][j].resize(3);
	    //-- numD start --//
	    distD[i][j] = 0;
	    distA[i][j] = 0;
	    //-- numD end --//
	    for(int k = 0; k < 3; k++)
	      ppUID[i][j][k] = 0;
	  }
	}



	for(int i = 1; i < total; i++){
	  
	  if(lpr[i] > 1e-6){

	    //-- numD start--//
	    for(int k = 0; k < np; k++){
	      maxlogBF = exp(lPP[k][i] - avBFs[k]);
	      distD[k][numD[i]-1] += maxlogBF;
	      distA[k][numA[i]-1] += maxlogBF;
	    }  
	    //-- numD end --//
	    for (int j = 0; j < nd; j++){
	      int aa = i % (int)pow(3,(j+1));
	      zz[j] = aa / (int)pow(3,j);

	      for(int k = 0; k < np; k++)
		ppUID[k][j][zz[j]] += exp(lPP[k][i] - avBFs[k]);

	    }	   

	  }

	}


	//-- numD start --//
	numD.resize(0);
	numA.resize(0);
	//-- numD end --//

	//--- Heejung End --//



	
	//		for (int i = 0; i < total; i++)
	//		{
//            for (int j = 0; j < nd; j++)
//            {
//                int aa = i % (int)pow(3,(j+1));
//                zz[j] = aa / (int)pow(3,j);
//            }
//			for (int j = 0; j < nd; j++)
//				outfile << zz[j]; 
//			outfile << "\t"; 
//			char buf[100]; 
//			sprintf(buf, "%.5f ", lpr[i]); 
//			outfile << buf; 
//			for (int m = 0; m < np; m++)
//			{
//				sprintf(buf, "%.5f ", mbf[i][m]); 
//				outfile << buf; 
//			}
//			outfile << endl; 
//		}


	//-- Heejung start --//

	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".mph.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(!outfile.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}
	


	vector<int> uniposi(0);
	uniposi.resize(nd);
	int count1, count2, posi1, posiAll;
	posi1 = 0;
	posiAll = 0;
	//-- Heejung end --//


	cout << "total = " << total << endl; 
	outfile << "rs "; 
        for (int i = 1; i < total; i++)
        {
	
	  //-- Heejung start --//
	  count1 = 0;
	  count2 = 0;
	  //-- Heejung end --//
	  for (int j = 0; j < nd; j++){
	    int aa = i % (int)pow(3,(j+1));
	    zz[j] = aa / (int)pow(3,j);
	    //-- Heejung start --//
	    if(zz[j] == 1){
	      count1++;
	      posi1 = j;
	    }else if(zz[j] == 2)
	      count2++;
	    //-- Heejung end --//
	    
	  }
	  //-- Heejung start --//
	  if((count1==1) & (count2 == (nd -1)))
	    uniposi[posi1] = i;

	  if(count1 == nd)
	    posiAll = i;
	  //-- Heejung end --//
	  for (int j = 0; j < nd; j++)
	    outfile << zz[j]; 
	  outfile << " ";

	}
	

	

	outfile << endl; 
		
	//for(unsigned i = 0; i < uniposi.size(); i++)
	//  cout << "uniposi : " << i << " " << uniposi[i] << endl;


	outfile << "prior "; 
	for (int i = 1; i < total; i++)
	  {
	    char buf[100]; 
	    if(lpr[i] < 1e-5) 
	      sprintf(buf, "%.5f ", lpr[i]); 
	    else 
	      sprintf(buf, "%+.5f ", lpr[i]); 
	    outfile << buf; 
	  }
	outfile << endl; 
	
	
	for (int m = 0; m < np; m++)
	  {
	    outfile << vsRsnum.at(m) << " "; 
	    for (int i = 1; i < total; i++)
	      {
		char buf[100]; 
		if((mbf[i][m]/log(10)) < 1e-5) 
		  sprintf(buf, "%.5f ", (mbf[i][m]/log(10))); 
		else 
		  sprintf(buf, "%+.5f ", (mbf[i][m]/log(10))); 
		outfile << buf; 
	      }
	    outfile << endl; 
	  }
	


	//--- Heejung start --//

	outfile.close(); 
	cout << sfn << " has been created." << endl;


	fstream outfile_prob; 
	string sfn_prob("output/");
	sfn_prob.append(fnOutput);
	sfn_prob.append(".mph.prob.txt");
	outfile_prob.open(sfn_prob.c_str(), ios::out);
	if(!outfile_prob.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}

	for (int m = 0; m < np; m++){

	  outfile_prob << vsRsnum.at(m) << " "; 
	  for (int i = 0; i < nd; i++){
	    char buf[100]; 
	    if(ppUID[m][i][0] < 1e-5) 
	      sprintf(buf, "%.5f ", ppUID[m][i][0]); 
	    else 
	      sprintf(buf, "%+.5f ", ppUID[m][i][0]); 
	    outfile_prob << buf; 

	    if(ppUID[m][i][1] < 1e-5) 
	      sprintf(buf, "%.5f ", ppUID[m][i][1]); 
	    else 
	      sprintf(buf, "%+.5f ", ppUID[m][i][1]); 
	    outfile_prob << buf; 

	  }
	  outfile_prob << endl; 
	}
	

	outfile_prob.close(); 
	cout << sfn_prob << " has been created." << endl; 
	

	fstream outfile_BFs; 
	string sfn_BFs("output/");
	sfn_BFs.append(fnOutput);
	sfn_BFs.append(".mph.BFs.txt");
	outfile_BFs.open(sfn_BFs.c_str(), ios::out);
	if(!outfile_BFs.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}
	

	 
	for (int m = 0; m < np; m++){

	  outfile_BFs << vsRsnum.at(m) << " "; 
	  char buf[100]; 
	  if((avBFs[m]/log(10)) < 1e-5) 
	    sprintf(buf, "%.5f ", (avBFs[m]/log(10))); 
	  else 
	    sprintf(buf, "%+.5f ", (avBFs[m]/log(10))); 
	  outfile_BFs << buf; 

	  if((mbf[posiAll][m]/log(10)) < 1e-5) 
	    sprintf(buf, "%.5f ", (mbf[posiAll][m]/log(10))); 
	  else 
	    sprintf(buf, "%+.5f ", (mbf[posiAll][m]/log(10))); 
	  outfile_BFs << buf; 
	  
	  for(int i = 0; i < nd; i++){
	    if((mbf[uniposi[i]][m]/log(10)) < 1e-5) 
	      sprintf(buf, "%.5f ", (mbf[uniposi[i]][m]/log(10))); 
	    else 
	      sprintf(buf, "%+.5f ", (mbf[uniposi[i]][m]/log(10))); 
	    outfile_BFs << buf; 
	  }
	  outfile_BFs << endl; 
	}
	

	outfile_BFs.close(); 
	cout << sfn_BFs << " has been created." << endl; 


	uniposi.resize(0);

	avBFs.resize(0);
	for(int i = 0; i < np; i++)
	  lPP[i].resize(0);
	lPP.resize(0);


	//-- numD start --//
	/*
	fstream outfile_distDA; 
	string sfn_distDA("output/");
	sfn_distDA.append(fnOutput);
	sfn_distDA.append(".mph.distNumDA.txt");
	outfile_distDA.open(sfn_distDA.c_str(), ios::out);
	if(!outfile_distDA.is_open()) 
	{
		cout << "can't open file ... " << endl;  
		exit(0); 
	}
	

	 
	for (int m = 0; m < np; m++){

	  outfile_distDA << vsRsnum.at(m) << " "; 
	  char buf[100];
	  for(int j = 0; j < nd; j++){
	    if(distD[m][j] < 1e-5) 
	      sprintf(buf, "%.5f ", distD[m][j]); 
	    else 
	      sprintf(buf, "%+.5f ", distD[m][j]); 
	    outfile_distDA << buf; 
	  }

	  for(int j = 0; j < nd; j++){
	    if(distA[m][j] < 1e-5) 
	      sprintf(buf, "%.5f ", distA[m][j]); 
	    else 
	      sprintf(buf, "%+.5f ", distA[m][j]); 
	    outfile_distDA << buf; 
	  }	  
	  outfile_distDA << endl; 
	}
	
	outfile_distDA.close(); 
	cout << sfn_distDA << " has been created." << endl; 
	//-- numD end --//

	*/


        for(int i = 0; i < np; i++){
	  for(int j = 0; j < nd; j++)
	    ppUID[i][j].resize(0);
	  ppUID[i].resize(0);
	  //-- numD start --//
	  distD[i].resize(0);
	  distA[i].resize(0);
	  //-- numD end --//
	}
	ppUID.resize(0);
	//-- numD start --//
	distD.resize(0);
	distA.resize(0);
	//-- numD end --//

	//--- Heejung end --//
	

	delete[] zz; 
	delete[] lbf; 
	delete[] lpr; 
	Free2DMatrix(mbf); 
    }
	
	//-- Heejung start --//
	//outfile.close(); 
	//cout << sfn << " has been created." << endl;
	//-- Heejung end --//
	Free2DMatrix(gt); 
	Free2DMatrix(ph); 
	
}

void ModelnData::logbf_rank1(double sa, double pi0, int pm, int * config, double * lbf)
{
    if(pm == 0) pm = nd -1;
                                 
//	cout << " nn, np, nd " << nn << " " << np << " " << nd << endl; 

    gsl_matrix * m1g = gsl_matrix_alloc(nn, np);
    for (int i = 0; i < nn; i++)
	{
		for (int j = 0; j < np; j++)
			gsl_matrix_set(m1g, i, j, gt[i][j]); 
	}
	
    gsl_matrix * m1y = gsl_matrix_alloc(nn, nd);
    for (int i = 0; i < nn; i++)
	{
		for (int j = 0; j < nd; j++)
			gsl_matrix_set(m1y, i, j, ph[i][j]); 
	}

    gsl_matrix * m1phi = gsl_matrix_alloc(nd, np);
    gsl_matrix * m1syy = gsl_matrix_alloc(nd, nd);
    gsl_vector * v1s11 = gsl_vector_alloc(np);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, m1y, m1g, 0, m1phi);  //phi=t(y) %*% g; 
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, m1y, m1y, 0, m1syy);  //syy=t(y) %*% y;

    for (int j = 0; j < np; j++)
    {
        double sum = 0;
        for (int i = 0; i < nn; i++)
        {
            double tt = gsl_matrix_get(m1g, i, j);
            sum += tt * tt;
        }
        gsl_vector_set(v1s11, j, sum);
    }
    // s11 = colsum(g*g); 
	
//	for (int i = 0; i < np; i++)
//		cout << gsl_vector_get(v1s11, i) << " "; 
//	cout << endl; 

	double prior = compute_prior(config, pi0);
	if(prior > 1e-6)
		logbf_summaries(m1phi, m1syy, v1s11, config, pm, sa, lbf);
	else
	{
		for (int m = 0; m < np; m++)
			lbf[m] = 0.0; 
	}
	gsl_matrix_free(m1g); 
	gsl_matrix_free(m1y); 
	gsl_matrix_free(m1phi); 
	gsl_matrix_free(m1syy); 
	gsl_vector_free(v1s11); 

}

void  ModelnData::logbf_summaries(gsl_matrix * m1phi, gsl_matrix * m1syy, gsl_vector * v1s11, int * zz, int pm, double sa, double * lbf)
{
	//zz is config; 
    int dd = 0;
    int du = 0;
    vector<int> vd;  //index of direct phenotype; 
    vector<int> vu;
    for (int i = 0; i < nd; i++)
    {
        if(zz[i] == 0) vu.push_back(i);
        if(zz[i] == 1) vd.push_back(i);
    }

    dd = (int) vd.size();
    du = (int) vu.size();

	if(dd == 0) 
	{
        for (int m = 0; m < np; m++)
			lbf[m] = 0; 
		return; 
	}

    if(du > 0)
    {
        gsl_matrix * syy = gsl_matrix_alloc(du, du);
        gsl_matrix * invsyy = gsl_matrix_alloc(du, du);
        for (unsigned i = 0; i < vu.size(); i++)
        {
            for (unsigned j = 0; j < vu.size(); j++)
                gsl_matrix_set(syy, i, j, gsl_matrix_get(m1syy, vu.at(i), vu.at(j)));
        }   //syy is a submatrix of syy=Y^t%*%Y; 

        gsl_matrix * phi0 =gsl_matrix_alloc(du, dd);
        for (unsigned i = 0; i < vu.size(); i++)
        {
            for (unsigned j = 0; j < vd.size(); j++)
                gsl_matrix_set(phi0, i, j, gsl_matrix_get(m1syy, vu.at(i), vd.at(j)));
			//this is confusing notation for phi0 is a submatrix of syy. 
        }

        gsl_matrix * phiu =gsl_matrix_alloc(du, np);
        for (unsigned i = 0; i < vu.size(); i++)
        {
            for (int j = 0; j < np; j++)
                gsl_matrix_set(phiu, i, j, gsl_matrix_get(m1phi, vu.at(i), j));
        }   //phiu is a sub-matrix of phi=Y^t%*%G;

		
		gsl_linalg_cholesky_decomp(syy); 
//		gsl_vector * tb = gsl_vector_alloc(du); 
//		gsl_vector * tx = gsl_vector_alloc(du); 
//
//        gsl_matrix * mc = gsl_matrix_alloc(du, np);
//		for (int col = 0; col < np; col++)
//		{
//			gsl_matrix_get_col(tb, phiu, col); 
//			gsl_linalg_cholesky_solve(syy, tb, tx); 
//            gsl_matrix_set_col(mc, col, tx); 
//		}
//
//        gsl_matrix * mb = gsl_matrix_alloc(du, dd);
//		for (int col = 0; col < dd; col++)
//		{
//			gsl_matrix_get_col(tb, phi0, col); 
//			gsl_linalg_cholesky_solve(syy, tb, tx); 
//			gsl_matrix_set_col(mb, col, tx); 
//		}
//		gsl_vector_free(tx); 
//		gsl_vector_free(tb); 

		for (int i = 0; i < du; i++)
			for (int j = i+1; j < du; j++)
				gsl_matrix_set(syy, i, j, 0); 
		gsl_permutation * perm;
		perm = gsl_permutation_alloc(du);
		int sig; 
		gsl_linalg_LU_decomp(syy, perm, &sig);  
		gsl_linalg_LU_invert(syy, perm, invsyy); 
		gsl_permutation_free(perm); 

        gsl_matrix * mc = gsl_matrix_alloc(du, np);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, invsyy, phiu, 0, mc);  //mc=(luu)^-1 %*% phiu;

        gsl_matrix * mb = gsl_matrix_alloc(du, dd);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, invsyy, phi0, 0, mb);  //mb=(luu)^-1 %*% phi0;

        gsl_matrix_free(syy);
        gsl_matrix_free(invsyy);
        gsl_matrix_free(phi0);
        gsl_matrix_free(phiu);


        gsl_vector * vc = gsl_vector_alloc(np);
        for (int j = 0; j < np; j++)
        {
            double sum = 0;
            for (int i = 0; i < du; i++)
            {
                double tt = gsl_matrix_get(mc, i, j);
                sum += tt * tt;
            }
            gsl_vector_set(vc, j, sum);
        }
        // vc = colsum(mc*mc); 

        gsl_vector_scale(vc, -1);
        gsl_vector_add(vc, v1s11);

        gsl_matrix * mu = gsl_matrix_alloc(dd, np);
        gsl_matrix * phid = gsl_matrix_alloc(dd, np);

        for (unsigned i = 0; i < vd.size(); i++)
        {
            for (int j = 0; j < np; j++)
                gsl_matrix_set(phid, i, j, gsl_matrix_get(m1phi, vd.at(i), j));
        }   //phid is a submatrix of m1phi=Y^t%*%G; 

        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, mb, mc, 0, mu);
        gsl_matrix_scale(mu, -1);
        gsl_matrix_add(mu, phid);

        gsl_matrix * rss0 = gsl_matrix_alloc(dd, dd);
        gsl_matrix * syyd = gsl_matrix_alloc(dd, dd);
        for (unsigned i = 0; i < vd.size(); i++)
        {
            for (unsigned j = 0; j < vd.size(); j++)
                gsl_matrix_set(syyd, i, j, gsl_matrix_get(m1syy, vd.at(i), vd.at(j)));
        }   //syyd is a submatrix of m1syy=Y^t%*%Y; 

        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, mb, mb, 0, rss0);
        gsl_matrix_scale(rss0, -1);
        gsl_matrix_add(rss0, syyd);

		gsl_linalg_cholesky_decomp(rss0); 
//        gsl_matrix * ma = gsl_matrix_alloc(dd, np);
//		tb = gsl_vector_alloc(dd); 
//		tx = gsl_vector_alloc(dd); 
//
//		for (int col = 0; col < np; col++)
//		{
//			gsl_matrix_get_col(tb, mu, col); 
//			gsl_linalg_cholesky_solve(rss0, tb, tx); 
//            gsl_matrix_set_col(ma, col, tx); 
//		}
//		gsl_vector_free(tx); 
//		gsl_vector_free(tb); 

		for (int i = 0; i < dd; i++)
			for (int j = i+1; j < dd; j++)
				gsl_matrix_set(rss0, i, j, 0); 

		perm = gsl_permutation_alloc(dd); 
		sig = 0; 
		gsl_linalg_LU_decomp(rss0, perm, &sig);  
		gsl_linalg_LU_invert(rss0, perm, syyd); 
		gsl_permutation_free(perm); 

        gsl_matrix * ma = gsl_matrix_alloc(dd, np);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, syyd, mu, 0, ma);

//		for (int i = 0; i < dd; i++)
//		{
//			for (int j = 0; j < dd; j++)
//				cout << gsl_matrix_get(syyd, i, j) << " "; 
//			cout << endl; 
//		}

        gsl_vector * va = gsl_vector_alloc(np);
        for (int j = 0; j < np; j++)
        {
            double sum = 0;
            for (int i = 0; i < dd; i++)
            {
                double tt = gsl_matrix_get(ma, i, j);
                sum += tt * tt;
            }
            gsl_vector_set(va, j, sum);
        }
        // va = colsum(ma*ma); 

        //lbf[1:np]; 
        double  df =  (nn + pm - (nd - dd - du));
        for (int m = 0; m < np; m++)
        {
            double tc = gsl_vector_get(vc, m);
            double lam = 1.0 /sa / sa / tc;
            double k = 1.0 / (1.0 + lam);

            lbf[m] = log(1.0-k) * dd / 2.0 - 0.5 * df  * log(1.0 - k / tc * gsl_vector_get(va, m));
        }

        gsl_vector_free(va);
        gsl_vector_free(vc);
        gsl_matrix_free(ma);
        gsl_matrix_free(syyd);
        gsl_matrix_free(rss0);
        gsl_matrix_free(phid);
        gsl_matrix_free(mu);

	    gsl_matrix_free(mb);
   		gsl_matrix_free(mc);
		
    }
    else
    {
        gsl_matrix * phid = gsl_matrix_alloc(dd, np);
        for (unsigned i = 0; i < vd.size(); i++)
        {
            for (int j = 0; j < np; j++)
                gsl_matrix_set(phid, i, j, gsl_matrix_get(m1phi, vd.at(i), j));
        }   //phid is a submatrix of phi=Y^t%*%G

        gsl_matrix * rss0 = gsl_matrix_alloc(dd, dd);
        gsl_matrix * syyd = gsl_matrix_alloc(dd, dd);
        for (unsigned i = 0; i < vd.size(); i++)
        {
            for (unsigned j = 0; j < vd.size(); j++)
                gsl_matrix_set(rss0, i, j, gsl_matrix_get(m1syy, vd.at(i), vd.at(j)));
        }    //syyd is a submatrix of syy=Y^t%*%Y; 
        

		gsl_linalg_cholesky_decomp(rss0); 
//        gsl_matrix * ma = gsl_matrix_alloc(dd, np);
//		gsl_vector * tb = gsl_vector_alloc(dd); 
//		gsl_vector * tx = gsl_vector_alloc(dd); 
//
//		for (int col = 0; col < np; col++)
//		{
//			gsl_matrix_get_col(tb, phid, col); 
//			gsl_linalg_cholesky_solve(rss0, tb, tx); 
//            gsl_matrix_set_col(ma, col, tx); 
//		}
//		gsl_vector_free(tx); 
//		gsl_vector_free(tb); 

		for (int i = 0; i < dd; i++)
			for (int j = i+1; j < dd; j++)
				gsl_matrix_set(rss0, i, j, 0); 


		gsl_permutation * perm = gsl_permutation_alloc(dd); 
		int sig; 
		gsl_linalg_LU_decomp(rss0, perm, &sig);  
		gsl_linalg_LU_invert(rss0, perm, syyd); 

        gsl_matrix * ma = gsl_matrix_alloc(dd, np);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, syyd, phid, 0, ma);
        gsl_permutation_free(perm); 

//		for (int i = 0; i < dd; i++)
//		{
//			for (int j = 0; j < np; j++)
//				cout << gsl_matrix_get(ma, i, j) << " "; 
//			cout << endl; 
//		}

        gsl_vector * va = gsl_vector_alloc(np);
        for (int j = 0; j < np; j++)
        {
            double sum = 0;
            for (int i = 0; i < dd; i++)
            {
                double tt = gsl_matrix_get(ma, i, j);
                sum += tt * tt;
            }
            gsl_vector_set(va, j, sum);
        }
        // va = colsum(ma*ma); 

        //lbf[1:np]; 
        double  df =  (nn + pm - (nd - dd - du));
        for (int m = 0; m < np; m++)
        {
            double tc = gsl_vector_get(v1s11, m);
            double lam = 1.0 /sa / sa / tc;
            double k = 1.0 / (1.0 + lam);

            lbf[m] = log(1.0-k) * dd / 2.0 - 0.5 * df  * log(1.0 - k / tc * gsl_vector_get(va, m));
        }

        gsl_matrix_free(phid); 
        gsl_matrix_free(rss0); 
        gsl_matrix_free(syyd); 
        gsl_matrix_free(ma); 
        gsl_vector_free(va); 
		
    }

}
                                         
double ModelnData::compute_prior(int * zz, double pi0)
{
	int n0, n1, n2; 
	n0 = n1 = n2 = 0; 
	for (int d = 0; d < nd; d++)
	{
		if (zz[d] == 0) n0++; 
		else if (zz[d] == 1) n1++; 
		else n2++;  
	}

	if(n1 == 0 && n2 == 0) return pi0; 
	else if(n1 == 0) return 0.0; 
	else  
		return (1.0 - pi0) / nd / (n1+n2) / nchoosek(nd, n0) / nchoosek(nd-n0, n1);  
	
}

