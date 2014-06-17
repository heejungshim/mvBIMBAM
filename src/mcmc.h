#if !define __MCMC__
#define __MCMC__

#include <string>
#include <map>
#include <vector>

using namespace std; 

class singleSNP
{
private:
	string rs; 
	real af; 
public:
	inline real Getaf(void) {return af;}
	inline string Getrs(void) {return rs;}
}

class multiSNP
{
private:
	map<string, int> map_rs01;  //hash if a snp is in or out; 
	vector<int> gamma;  		//position of 1's in gamma vector; 
	real h; 					//heritability;                                                         
	//note: the state space of our MCMC is the product space of vector gamma and scalar h; 
	real sigma_a; 
	vector<class singleSNP> vsnp; 
public:
	void updateh();
	real calc_prior(h); 
	void updategamma(); 
}

#endif


void multiSNP::updateh()
{
	real prop =( gsl_rand(gsl) -0.5) * 0.1;
	real newh = h + prop;
	if (newh > 0.5) newh = 1.0 - newh; 
	if (newh < 0.0) newh = -newh;
	
	real prior_new = calc_prior(newh);
	real prior_old = calc_prior(h); 

	real bf_new = calc_bf(prior_new, gamma); 
	real bf_old = calc_bf(prior_old, gamma); 

	if(gsl_rand(gsl) < exp(bf_new - bf_old))
	{
		h = newh; 
	}
	
}


real multiSNP::cacl_prior(real hh)
{
	real temp = 0.0; 
	for (unsigned i = 0; i < gamma.size(); i++)
	{
		class singleSNP * ps = vsnp.at(gamma.at(i)); 
		real af = ps->Getaf(); 
		temp += af * (1.0 - af); 
	}
	temp *= 2.0; 
	temp = 1.0 / temp * hh / (1.0 - h);
	return (sqrt(temp));
	
}

void updategamma()
{
	//should take into account of r^2; 
	//should take into account of number of SNPs in the study; 
	int max = 20; 
	int min = 2; 
	if(gsl_rand(gsl) < 0.5 && gamma.size() < max) 
	{
		//pick a snp and try to put into gamma list; 
		
	}
	else //decrease; 
	{
		//pick a snp and try to through out of gamma list; 
	}
	
}
