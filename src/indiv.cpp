#include <stdio.h>
#include <stdlib.h>

#include "indiv.h"
#include <iostream>
using namespace std;

Individual::Individual()
{
	m_expectJmk = NULL;
	snpGT = NULL; 
	phiScale = NULL; 
	betaScale = NULL; 
	nMissingGT = 0; 
#if defined (IMPUTATION)
	pMask = NULL;
	maf = NULL; 
#endif
}



Individual::~Individual()
{  
	huskyID.assign("\0"); 
	if(snpGT) {Free1D((void*) snpGT); snpGT = NULL;}
	if(phiScale) {delete[] phiScale;       phiScale = NULL;}
	if(betaScale) {delete[] betaScale;     betaScale = NULL;} 
	if(m_expectJmk) {Free2DMatrix(m_expectJmk);   m_expectJmk = NULL;} 
#if defined (IMPUTATION)
	if(pMask) {delete[] pMask; pMask = NULL;}
	delete[] maf;  
#endif
//	cout << "indiv destructor being called : " << snpGT << endl;    
}
