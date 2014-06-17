#ifndef __PARAM_H__
#define __PARAM_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include <iostream>
using namespace std;
		 
//differenet populations share theta while having different alpha and r;
class ModelParam
{                                    
private:
	real ** theta;  			// MxK, emission para;
	real *** alpha; 
	real ** r; 
public:
	ModelParam();
	~ModelParam(); 
	void Init(int, int, int); 
	inline void Settheta(int m, int k, real s) {theta[m][k] = s;}
	inline void Setalpha(int p, int m, int k, real s) {alpha[p][m][k] = s;}
	inline void Setr(int p, int m, real s) {r[p][m] = s;}
	//set is used to update parameters, one by one; 
	//while get functions return pointers;
	inline real ** Gettheta(void) {return theta;}            
	inline real ** Gettheta(int p) {return theta;}            
	
	inline real ** Getalpha(void) {return alpha[0];} 
	inline real * Getr(void) {return r[0];}
	//if don't specify which subpoplation, return the first one; 
	inline real ** Getalpha(int p) {return alpha[p];} 
	inline real * Getr(int p) {return r[p];}
};

#endif
