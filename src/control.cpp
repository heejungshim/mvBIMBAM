#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdio.h>    
#include <stdlib.h>
#include "control.h"
#include "model.h"
#include "indiv.h"
#include "diploid.h"
#include "haploid.h"
#if defined (READLINE)
#include "readline/readline.h"
#include "readline/history.h"
#include "ncurses.h"
#endif
#if defined (MPI_ENABLED)
#include "mpi.h"
#endif 
#include "fpmath.h"
#include "fp.h" 

using namespace std;
#define com_empty 	0
#define com_geno	1
#define com_pheno	2
#define com_out		3
#define com_pos		4
#define com_struct  5
#define com_nPH		6
#define com_clean   7
#define com_rand 	8
#define com_em 		11
#define com_step 	12
#define com_warm 	13
#define com_cluster 14
#define com_impute  20
#define com_multi_snp	21
#define com_level 	23
#define com_pval 	24
#define com_msnp 	25
#define com_rem     30
#define com_sem     31
#define com_sigma_a 41
#define com_sigma_d 42
#define com_sort    43
#define com_gene    2000
#define com_flank   2001

#define com_version 99
#define com_help    100
#define com_wmg		101
#define com_wbg     102   //write best guess; 
#define com_wgd     103   //write genotype distribution; 
#define com_weg     104

#define com_casectrl       1001
#define com_cluster_method 1002
#define com_mask_impute    1003
#define com_ssd 	3000
#define com_psd     3001
#define com_file    3002
#define com_doc     3005
#define com_num     3003
#define com_df      3004
#define com_mcmc    3210
#define com_note    3300
#define com_test    3400
#define com_allele_coding 4000
#define com_exclude_maf 4001
#define com_exclude_miss 4005
#define com_exclude_nopos 4006

#define com_nobf 4002
#define com_reduce_indiv 4003
#define com_silence 4004
#define com_read_mode 5000
#define com_not_snp 5001

#define com_multi_ph 6001

CtrlParam::CtrlParam(void)
{
	pMD = NULL;
#if defined (READLINE)
	real rand = gsl_rng_uniform(gsl_r); 
	if(rand < 0.5) bimbam.assign("BimBam $ ");
	else bimbam.assign("BimBam ¥ ");
#endif
	hcom[""] = com_empty;
	hcom["-g"] = com_geno; 
	hcom["-gen"] = com_geno;
	hcom["-p"] = com_pheno; 
	hcom["-phe"] = com_pheno; 
	hcom["-o"] = com_out;
	hcom["-out"] = com_out;
	hcom["-pos"] = com_pos; // rspos file name;
	hcom["-P"] = com_pos; // rspos file name;
	hcom["-S"] = com_struct; 
	hcom["-f"] = com_nPH; 
	hcom["-F"] = com_nPH;
	hcom["-r"] = com_rand;
	hcom["-R"] = com_rand;
	hcom["-rand"] = com_rand;
	hcom["-clean"] = com_clean;
//1* for em related. 	
	hcom["-e"] = com_em;
	hcom["-em"] = com_em;
	hcom["-s"] = com_step;
    hcom["-w"] = com_warm;
	hcom["-c"] = com_cluster;
	hcom["-k"] = com_cluster; 
	hcom["-rem"] = com_rem; 
	hcom["--rem"] = com_rem; 
	hcom["-sem"] = com_sem; 
	hcom["--sem"] = com_sem; 
//2* for snp calc. 
	hcom["-i"] = com_impute;
	hcom["-m"] = com_multi_snp;
	hcom["-mul"] = com_multi_snp; 
	hcom["-l"] = com_level; 
	hcom["-pval"] = com_pval;
	hcom["-msnp"] = com_msnp;
	hcom["-cc"] = com_casectrl; 
	hcom["-cam"] = com_cluster_method;
	hcom["-mask"] = com_mask_impute; 
	hcom["-gene"] = com_gene; 
	hcom["-GF"] = com_flank; 
	hcom["-gf"] = com_flank; 
	hcom["--ssd"] = com_ssd;   //snp_summary_data;
	hcom["-ssd"] = com_ssd;   //snp_summary_data;
	hcom["--psd"] = com_psd; 
	hcom["-psd"] = com_psd; 
//help;     
	hcom["-v"] = com_version; 
	hcom["-ver"] = com_version; 
	hcom["-h"] = com_help;
	hcom["-help"] = com_help; 
	hcom["--help"] = com_help; 
	hcom["-wmg"] = com_wmg; 	//write mean genotype
	hcom["--wmg"] = com_wmg; 	//wriet mean genotype
	hcom["-wbg"] = com_wbg;
	hcom["--wbg"] = com_wbg; 
	hcom["-wgd"] = com_wgd; 
	hcom["--wgd"] = com_wgd; 
	hcom["--weg"] = com_weg; 
	hcom["-weg"] = com_weg; 
	hcom["-file"] = com_file; 
	hcom["-doc"] = com_doc; 
	
	hcom["-num"] = com_num;
	hcom["-df"] = com_df;
	hcom["-mcmc"] = com_mcmc;
	hcom["--mcmc"] = com_mcmc;
	hcom["-note"] = com_note;
	hcom["-test"] = com_test;
	hcom["-ac"] = com_allele_coding; 
	hcom["-allele-coding"] = com_allele_coding; 
	hcom["-exclude-maf"] = com_exclude_maf; 
	hcom["--exclude-maf"] = com_exclude_maf; 
	hcom["-exclude-miss"] = com_exclude_miss; 
	hcom["--exclude-miss"] = com_exclude_miss; 
	hcom["-exclude-nopos"] = com_exclude_nopos; 
	hcom["--exclude-nopos"] = com_exclude_nopos; 
	hcom["--nobf"] = com_nobf;
	hcom["-nobf"] = com_nobf;
	hcom["--reduce-indiv"] = com_reduce_indiv; 
	hcom["--silence"] = com_silence; 
	hcom["-silence"] = com_silence; 
	hcom["--silent"] = com_silence; 
	hcom["-silent"] = com_silence; 
		 
	hcom["-A"] = com_sigma_a; 
	hcom["-a"] = com_sigma_a; 
	hcom["-D"] = com_sigma_d; 
	hcom["-d"] = com_sigma_d; 
	hcom["--sort"] = com_sort;
	hcom["-sort"] = com_sort;
	hcom["-gmode"] = com_read_mode; 
	hcom["--notsnp"] = com_not_snp; 
	hcom["-notsnp"] = com_not_snp; 
	hcom["-mp"] = com_multi_ph; 
	hcom["-mph"] = com_multi_ph; 
} 

CtrlParam::~CtrlParam(void)
{
	;
} 

void CtrlParam::OnHelp(void)
{
	cout << endl; 
	cout << " BIMBAM Version 0.99 released in Sept. 2008 by Y. Guan and M. Stephens" << endl; 
	cout << " Usage: ./bimbam -g panel.txt -p 0 -g genotype.txt -p phenotype.txt -i 1 -o prefix" << endl;
	cout << " " << endl;
	cout << " FILE I/O RELATED OPTIONS." << endl;
	cout << " -g <filename>    " << " repeatable for multiple input files. must pair with -p" << endl;
	cout << " -p <file/0/1/z>  " << " repeatable for multiple input files. must pair with -g" << endl; 
	cout << "    <0> pairing genotypes are panel, <z/1> pairing genotype individuals have phenotype 0/1" << endl;
	cout << " -pos <file>      " << " repeatable for multiple input files" << endl; 
	cout << " -f <num>         " << " specify number of phenotypes (columns in the phenotype files)" << endl;  
    cout << " -o <prefix>      " << " prefix of all output files, random seeds will be used by default" << endl;  
	cout << " -weg <0/1>       " << " write exact genotype, missing denote by NA. default 0" << endl; 
	cout << "                  " << " <0> cohort in numerical format, <1> cohort in bimbam format" << endl; 
	cout << " -wmg <0/1>       " << " write mean genotype. defualt 1" << endl;
	cout << " -wbg <0/1>       " << " write best guess genotype, in ACGT+- form. default 0" << endl;
	cout << " -wgd <0/1>       " << " write genotype distribution, pr(0), pr(1) for each genotype. default 1" << endl;
	cout << "                  " << " <0> write cohort only, <1> write both panel and cohort" << endl; 
	cout << " " << endl; 
	cout << " EM RELATED OPTIONS." << endl;
	cout << " -e(em) <num>     " << " specify number of EM runs, default 10" << endl; 
	cout << " -w(warm) <num>   " << " specify steps of warm up EM run, default 10" << endl; 
	cout << " -s(step) <num>   " << " specify steps of each EM run, default 0" << endl; 
	cout << " -c <num>         " << " specify number of clusters in EM algorithm, default 20" << endl; 
	cout << " -r <num>         " << " specify random seed, system time by default" << endl; 
	cout << " -sem <file>      " << " to save EM results, if <file> missing prefix.em will be used" << endl; 
	cout << " -rem <file>      " << " to read EM results" << endl; 
	cout << " " << endl; 
	cout << " MULTI-SNP CACULATION RELATED OPTIONS."  << endl; 
	cout << " -i <num>         " << " specify number of imputations to calc BF, default <0>" << endl;
	cout << "                  " << " <0> no imputation <1> use mean genotype, >100 importance sampling" << endl;
	cout << " -m <num>         " << " specify number of snps for multiple SNP study" << endl;
	cout << " -l <num>         " << " specify maximum number of SNPs in all combinations" << endl;
	cout << " -gene <file>     " << " to read gene file that specify regions of interests" << endl; 
	cout << " -gf(GF) <num>    " << " pair with -gene to specify gene flanking region in kb" << endl; 
	cout << " " << endl;
	cout << " BAYES FACTOR RELATED OPTIONS." << endl;
	cout << " -a(A) <num>         " << " repeatable, specify priors for additive effect. must pair with -D" << endl; 
	cout << " -d(D) <num>         " << " repeatable, specify priors for dominant effect. must pair with -A" << endl;
	cout << " -df <num>        " << " <1> additive effect, <2> (default) additive and dominance effect" << endl;
	cout << " -pval <num>      " << " calculate p-values via permuations." << endl; 
	cout << " -cc              " << " calc bf of logit regression on binary phenotype." << endl; 
	cout << " " << endl; 
	cout << " COMBINE STUDIES." << endl; 
	cout << " -psd <file_name> " << " to convert genotype and pheotype to summary statitics and save to a file" << endl; 
	cout << " -ssd <file_name> " << " take multiple summary data and calculate BF after combining them" << endl; 
	cout << " " << endl; 
	cout << " OTHER OPTIONS." << endl;
	cout << " -v(ver)          " << "print version and citation" << endl;
	cout << " -h(help)         " << "print this help" << endl;
	cout << " -exclude-maf <num> " << endl; 
	cout << "          exclude SNPs whose maf < num, default 0.01" << endl; 
	cout << " -exclude-miss <num> " << endl; 
	cout << "          exclude SNPs whose missing rate > num, default 1" << endl; 
	cout << " -exclude-nopos <num> " << endl; 
	cout << "          exclude SNPs that has no position information, 1 = yes (default), 0 = no" << endl; 
	cout << " --silence        " << "no terminal output" <<  endl; 
	cout << " " << endl; 
	cout << " Screen running away? try: ./bimbam -h | less" << endl; 
	cout << " " << endl; 
}

void CtrlParam::PrintHeader(void)
{
	cout << endl; 
	cout << " BIMBAM version 0.99a, visit http://stephenslab.uchicago.edu for possible update." << endl;
	cout << " Developed by Yongtao Guan ytguan.at.gmail.com, all rights to be lefted to GNU GPL." << endl; 
	cout << " References: Guan and Stephens (2008), Servin and Stephens (2007), Scheet and Stephens (2006)." << endl;
	cout << " Updated March 14 2009." << endl; 
	cout << endl; 
}

void CtrlParam::BatchRun(int argc, char ** argv)
{
	pMD = new class ModelnData;
                             
	int silence = 0; 
	int needem = 0; 
	int permute = 0;    //how many permutation for pval option; 
	int msnp = 0; 
	int yes_weg = 0; 
	int yes_wmg = 0; 
	int yes_wbg = 0;   
	int yes_wgd = 0; 
	int all_weg = 0; 
	int all_wmg = 1; 
	int all_wbg = 0;   //best guess default only write cohort snps; 
	int all_wgd = 1; 
	int gene = 0; 
	int rem = 0; //read em param; 
	int sem = 0; //save em param; 
	int ssd = 0; //snp summary data; 
	int psd = 0; //produce summary data; 
	int mcmc = 0; 
	int file = 0; 
	int test = 0; 
	int nobf = 0; 
	int mph = 0; 
	string fnEM; 
	string fnGene;  
	string fnssd; 
	string fnpsd;        
	string note; 
	int cluster_method = 0;

	int gmode = 0; 
	
#if defined (IMPUTATION)
	int mask = 0; 
#endif
	for(int i = 1; i < argc; i++) 
	{   
		string str;  
		string opt; 
		if(argv[i][0] != '-') 
			continue;
		
		str.assign(argv[i]);
		opt.assign(str, 0, str.length());

		map<string, int> :: iterator iter; 
		iter = hcom.find(opt); 
		if(iter == hcom.end())
		{
			cout << "-bimbam: unknown option: " << opt << endl; 
			safe_exit(); 
		}
		
		switch (hcom[opt]) 
		{			
			case com_not_snp:
				pMD->Setmnotsnp(1); 
				break;
			case com_read_mode:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				gmode = atoi(argv[i+1]); 
				break;
			case com_silence:
				silence = 1; 
				pMD->Setsilence(1);
				break; 
			case com_reduce_indiv:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->Set_reduced_num_indiv(atoi(argv[i+1]));
				break; 
			case com_nobf:
				nobf = 1; 
				break; 
			case com_exclude_maf:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') 
					continue;
				{
					double exmaf = atof(argv[i+1]); 
					if(exmaf == 0)
						pMD->Setexcludemaf(-100);
					else 
						pMD->Setexcludemaf(exmaf);
				}
				break; 
			case com_exclude_miss:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->Setexcludemiss(atof(argv[i+1]));
				break; 
			case com_exclude_nopos:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') 
					continue;
				
				{
				 	int temp = atoi(argv[i+1]); 
					if(temp != 0 && temp != 1) 
					{
						cout << "-bimbam: bad argument follow option: " << opt << endl; 
						safe_exit(); 
					}
					pMD->Setexcludenopos(temp);
				}
				break; 
			case com_allele_coding:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->Setallelecoding(atoi(argv[i+1]));
				break; 
			case com_test:
				test = 1; 
				break;
			case com_note:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				note.assign(argv[i+1]);
				break;
			case com_file:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
#ifndef WINDOWS
				str.clear();
#else
				str.resize(0); 
#endif
				str.assign(argv[i+1]);
				pMD->SetfnFILE(str);
				file = 1; 
				break;
			case com_doc:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.clear();
				str.assign(argv[i+1]);
				pMD->SetfnDOC(str);
				break;
			case com_geno:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
#ifndef WINDOWS
				str.clear();
#else
				str.resize(0); 
#endif
				str.clear();
				str.assign(argv[i+1]);
				pMD->SetvGin(str);
				break;
			case com_pheno:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.assign(argv[i+1]);
				pMD->SetvPin(str);
				break;
			case com_out:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.assign(argv[i+1]);
				pMD->SetfnOutput(str);
				break;
			case com_pos:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.clear(); 
				str.assign(argv[i+1]);
				pMD->SetvPos(str); 
				break;
			case com_struct:  
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->SethasPopStructure(argv[i+1]);
				break;
			case com_rand: // seed 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetrandSeed(atoi(argv[i+1]));
				break;
			case com_em:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnEMRuns(atoi(argv[i+1]));
				break;
            case com_num:
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    cout << "wrong augument after option." << endl;
                    exit(0);
                }
                pMD->Setnum(atoi(argv[i+1]));
                break;
				
            case com_df:
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    cout << "wrong augument after option." << endl;
                    exit(0);
                }
				if(atoi(argv[i+1]) != 1 && atoi(argv[i+1]) != 2)
				{
					cout << "degree of freedom must be either 1 or 2 !" << endl; 
					exit(0); 
				}
                pMD->Setdf(atoi(argv[i+1]));
                break;

			case com_step:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnMaxSteps(atoi(argv[i+1]));
				break;			
			case com_warm:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnWarmSteps(atoi(argv[i+1]));
				break;
			case com_cluster:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnK(atoi(argv[i+1]));
				break;
			case com_impute:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				{
					int temp = atoi(argv[i+1]); 
				pMD->SetnImpute(temp);
				if(temp > 100 || temp == 1) needem = 1;
				}
				break;				
			case com_multi_snp:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnMultiSnp(atoi(argv[i+1]));
				break;
			case com_nPH:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnPH(atoi(argv[i+1]));
				break;

			case com_help:
				PrintHeader(); 
				OnHelp();
				safe_exit();
				break;

			case com_version:
				PrintHeader(); 
				safe_exit(); 
				break; 
			case com_sigma_a:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}        
				{
					double ta = atof(argv[i+1]); 
					if(fabs(ta) < 1e-10) ta = 1e-10; 
					pMD->Setsigmaa(ta);
				}
				break;
				
			case com_sigma_d:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}        
				{
					double td = atof(argv[i+1]); 
					if(fabs(td) < 1e-10) { pMD->Setdf(1); td = 1e-10;}
					pMD->Setsigmad(td);
				}
				break;

			case com_sort:
				pMD->SetsortQ(1); 
				break; 
				
			case com_level:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnLevel(atoi(argv[i+1]));
				if(atoi(argv[i+1])) msnp = 1; 
				break;
			case com_msnp:
                msnp = 1; 
				break; 
			case com_mcmc:
                mcmc += 1; 
				break; 

			case com_pval:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				permute = atoi(argv[i+1]);   
				break; 
				
			case com_weg:
				needem = 0; 
				yes_weg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_weg = atoi(argv[i+1]); 
				break;
			case com_wmg:
				needem = 1; 
				yes_wmg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wmg = atoi(argv[i+1]); 
				break;
				
			case com_wbg:
				needem = 1; 
				yes_wbg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wbg = atoi(argv[i+1]); 
				break;
				
			case com_wgd:
				needem = 1; 
				yes_wgd = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wgd = atoi(argv[i+1]); 
				break;
   			
			case com_casectrl:
				pMD->SetOnCaseCtrl(); 
				break; 

			case com_cluster_method:
				cluster_method = 1; 
				needem = 1; 
				break; 

#if defined (IMPUTATION)				
			case com_mask_impute:
				needem = 1; 
				mask = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->Setmask(atof(argv[i+1])); 
				break; 
#endif
			case com_gene:
				gene = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fnGene.assign(argv[i+1]);
				break;
				
			case com_flank:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->Setflank(atoi(argv[i+1])); 
				break; 
			case com_sem:
				needem = 1;
				sem = 1; 
				break;
				
			case com_rem:
				rem = 1;
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fnEM.assign(argv[i+1]);
				break;
			case com_ssd:
				ssd = 1;
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->SetvSSD(argv[i+1]);
				break;
			case com_psd:
				psd = 1;
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fnpsd.assign(argv[i+1]);
				break;

			case com_multi_ph:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				mph = atoi(argv[i+1]);
				break;
				
			default:
				fprintf(stderr,"Bad option %s\n", argv[i]);
				OnHelp();
				exit(0);
				break; 
		}              
	}
    
	if(gmode > 0) 
	{
		needem = 0; 
		if(pMD->GetnImpute() == 0)
			pMD->SetnImpute(1); 
	}
	
	streambuf * cout_strbuf(cout.rdbuf()); 
	//save current cout into a streambuf; 
	ostringstream output; 
	if(silence == 1) 
		cout.rdbuf(output.rdbuf()); 
	//redirect cout to a streambuf output; 

	int randseed = 0;
#if defined (MPI_ENABLED)    
    int nProc;
    int procID;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    pMD->SetnProc(nProc);
    pMD->SetprocID(procID);
    if (procID == 0)
    {
        cout << "MPI is up and runing..." << endl;
        cout << "total processors in use = " << pMD->GetnProc() << endl;
	}
    if (pMD->GetprocID() == 0)
#endif
    {
		randseed = pMD->GetrandSeed();  // default returns -1;
		if (randseed <= 0)
			randseed = (unsigned) time(NULL);
		gsl_rng_set(gsl_r, randseed);
    }

	string prefix(pMD->GetfnOutput()); 
	if(prefix.compare("\0") == 0) {
		char tmp[100];
		memset(tmp, '\0', 100);
		sprintf(tmp, "%d", randseed); 
		prefix.assign(tmp); 
		pMD->SetfnOutput(prefix); 
	}
	//if no random seed provided, use system clock;    
	//if no prefix provided, use randseed; 
	
#if defined (MPI_ENABLED)
  	MPI_Bcast(&rs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pMD->SetrandSeed(randseed);
    gsl_rng_set(gsl_r, randseed + procID);
	if(nProc > 1 && procID == 0 )
        cout << "Broadcasted random seed..." << randseed;
	MPI_Barrier(MPI_COMM_WORLD); 
	if(nProc > 1 && procID == 0)
		cout << " ... done." << endl; 
	if (procID == 0)
#endif 
	{
		pMD->open_log(); 
		(pMD->fplog) << "## " << note << endl;
		(pMD->fplog) << "## COMMAND: "; 
		for (int i = 0; i < argc; i++)
			(pMD->fplog) << argv[i] << " "; 
		(pMD->fplog) << endl; 
		(pMD->fplog) << "## randseed = " << randseed << endl; 
	}   // write command line and notes into log; 
	
	//-- Heejung start --//
	if(mph == 0) 
	  pMD->process_prior(); 
	//pMD->process_prior();
	//-- Heejung end --//
	
//	pMD->for_barbara(); 
//   return; 
	
//	pMD->write_for_gain(); 
//	return; 
//  pMD->fix_mean_genotype(); 
//	return; 

	
	vector<class cGENE> vGene;
    if(gene == 0 || read_gene_map(fnGene, vGene) == 0) 
	{
		class cGENE tg("", 0, 0, -1); 
		vGene.push_back(tg); 
	}
	
	for (unsigned g = 0; g < vGene.size(); g++)
	{
	
#if defined (MPI_ENABLED)
		if (procID == 0 ) 
#endif      
		{
			if (vGene.size() > 1)
			{
				pMD->fplog << vGene.at(g).get_name() << "\t" << vGene.at(g).get_chr() << "\t"\
					 << vGene.at(g).get_beg() << "\t" << vGene.at(g).get_end() << endl; 
			}
		}
		string newpref(prefix); 
		if(vGene.size() > 0 && vGene.at(0).get_name().size() > 0)
			newpref.append(".");  
		newpref.append(vGene.at(g).get_name()); 
		pMD->SetfnOutput(newpref); 
		// get the new prefix for each gene; 

#if defined (POWER)
//		pMD->read_is_write(); 
//		return; 
//for IS of causal SNP only; 

// 		pMD->read_gd_stat_write(); 
// 		return; 
//for read in genotype density and calc different statistics only; 
#endif

		int read_dip, read_phe; 
		if(ssd == 1) 
		{
			pMD->CombineStudy();
			return;
		}
		if(gmode > 0) 
		{
		 	read_dip = pMD->read_bimbam_genotype_distribution(gmode, mcmc, vGene.at(g).get_beg(), vGene.at(g).get_end()); 
			read_phe = pMD->read_bimbam_phenotype(mcmc); 
		}
		else 
		{
		    read_dip = pMD->read_bimbam_genotype(mcmc, vGene.at(g).get_beg(), vGene.at(g).get_end()); 
		    read_phe = pMD->read_bimbam_phenotype(mcmc); 
		}
		if(psd == 1) 
		{
			if(read_dip == 1 && read_phe == 1) 
				pMD->ProduceSummaryData(gmode, fnpsd); 
			else 
				continue; 
		}

		if(read_dip == 0) continue;
		else if(needem == 1) pMD->assign_population_label();
		//divide into different subpopulations. 

//////////////////////////////////////
/////////	simulate phenotypes; 
//		pMD->simulate_phenotype(40, 0.02); 
//		return; 
		if(mph > 0) 
		{
			pMD->single_snp_multi_phenotype(mph); 
			return; 
		}
/////////////////////////////////////
		
#if defined (IMPUTATION)	
		if(mask == 1)                         
		{
			if(pMD->Getmask() < 1)
				pMD->MaskRand();
			else 
				pMD->MaskSNP(); 
		}   //for mask and imputation; 
#endif
		 
        time_t sec_beg, sec_end; 
		sec_beg = time(NULL); 
		if(needem == 1) 
		{
			pMD->InitModelParam(); 
			int needwarm = 0; 
			if(rem == 0 || (rem == 1 && pMD->read_em_param(fnEM) == 0)) //no em data or read em data failed; 
			{
				needwarm = 1; 
				if(pMD->GetnPanel() == 0 || pMD->GetnWarmSteps() == 0) 
					needwarm = 0;
			}
			
			if(needwarm) pMD->EM(1);
			if(pMD->GetnCohort() > 0 && pMD->GetnMaxSteps() > 0) pMD->EM(0); 
			if(sem == 1) pMD->write_em_param(); 
			if(pMD->GetnEMRuns() > 0) pMD->SNP_Density(); 
		}
		sec_end = time(NULL); 
		(pMD->fplog) << "## EM seconds used = " << sec_end - sec_beg << endl; 
		
#if defined (IMPUTATION)
		if(mask == 1) 
			pMD->ImputeMasked();
#endif 
		
//#if !defined (MPI_ENABLED)
//		if(cluster_method == 1) 
//		{ 
//			pMD->single_snp_cluster();
//			pMD->close_log(); 
//			cout << "-bimbam: write a file in output subdir: " << newpref << ".cluster.txt" << endl; 
//			return; 
//		}
//#endif
		
		if(yes_weg == 1) pMD->write_genotype(0, all_weg); 
		if(yes_wmg == 1) pMD->write_genotype(1, all_wmg); 
		if(yes_wbg == 1) pMD->write_genotype(2, all_wbg); 
		if(yes_wgd == 1) pMD->write_genotype(3, all_wgd); 
		
		if(nobf == 0) 
		{
			pMD->single_snp_analysis(permute); 
			if(msnp == 1) 
			{
				int action = gmode; 
				pMD->multi_snp_analysis(action);
				//if gmode == 0; then either using exact gt or do imputation; 
			}
		}

#if defined (MPI_ENABLED)
		if (procID == 0)
#endif
		{
			pMD->fplog << "## Bimbam generate following files in the output directory." << endl; 
			pMD->fplog << "## " << newpref << ".snpdata.txt" << endl; 
			if(!nobf) pMD->fplog << "## "<< newpref << ".single.txt" << endl;
			if(!nobf) pMD->fplog << "## "<< newpref << ".summary.txt" << endl; 
			if(msnp == 1) pMD->fplog << "## "<< newpref << ".multi.txt" << endl; 
			if(yes_weg == 1) pMD->fplog << "## "<< newpref << ".exact.genotype.txt" << endl;
			if(yes_wmg == 1) pMD->fplog << "## "<< newpref << ".mean.genotype.txt" << endl;
			if(yes_wbg == 1) pMD->fplog << "## "<< newpref << ".best.guess.genotype.txt" << endl; 
			if(yes_wgd == 1) pMD->fplog << "## "<< newpref << ".genotype.distribution.txt" << endl; 
			if(sem == 1) pMD->fplog << "## "<< newpref << ".em" << endl; 
			if(ssd == 1) pMD->fplog << "## "<< newpref << ".combined.study.txt" << endl;
			if(psd == 1 && fnpsd.size() > 0) pMD->fplog << "## "<< fnpsd << endl; 
			if(psd == 1 && fnpsd.size() == 0) pMD->fplog << "## "<< newpref << ".snp.summary.data.ssd" << endl; 
			
	    	(pMD->fplog) << "## random seed = " << randseed << endl;
			pMD->close_log();
			cout << "-bimbam: finished, for details see log: "<<  newpref << ".log.txt" << endl; 
		} 
		pMD->clean_geno_pheno(); 
	}

	if(silence == 1)
		cout.rdbuf(cout_strbuf); 
	//restore the original cout; 
	delete pMD; 
}

int CtrlParam::read_gene_map(string fnGene, vector<class cGENE>& vGene) 
{                                            
	fstream infile;
	string delimit(";, \t");
	streambuf * pbuf;
	infile.open(fnGene.c_str(), ios::in);
	if(!infile.is_open()) 
	{
		cout << "-bimbam: cannot open gene file " << fnGene << endl; 
		safe_exit();  
	}   // file;
	
	pbuf = infile.rdbuf(); 	
	while(pbuf->sgetc() != EOF) 
	{
		string line(getline(pbuf));
		int beg = line.find_first_not_of(delimit,0);
		int end = line.find_first_of(delimit); 
		string name(line.substr(beg, end-beg)); 
		beg = line.find_first_not_of(delimit, end);
		end = line.find_first_of(delimit, beg); 
		int chr = atoi(line.substr(beg, end-beg).c_str());
		beg = line.find_first_not_of(delimit, end);
		end = line.find_first_of(delimit, beg); 
		long start = (long) atoi(line.substr(beg, end-beg).c_str());
		beg = line.find_first_not_of(delimit, end);
		end = line.find_first_of(delimit, beg); 
		long send = (long) atoi(line.substr(beg, end-beg).c_str());
		class cGENE cg(name, chr, start, send);
		vGene.push_back(cg); 
	}
	infile.close(); 
	return (vGene.size() > 0);  
}
