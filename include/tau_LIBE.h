#ifndef __tau_LIBE__
#define __tau_LIBE__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "T_min.h"
#include "input.h"
#include "g_minus_2_utilities.h"
#include "HVP.h"
#include "highPrec.h"
#include "Spectral.h"
#include "tau_decay_LIBE_ISO.h"



using namespace std;

class FLAV {
  
public:
  FLAV() {}
  FLAV(double c,const distr_t& m,const distr_t& mcr) : q(c)  { dm=m; dmcr=mcr;}
  double q;
  distr_t dm;
  distr_t dmcr;
};



void Compute_tau_LIBE();
void Read_file_QED(string path, VVVfloat &QED, int reim);
void Read_file_SIB(string path, VVVfloat &SIB, int r1, int r2, int reim);
distr_t_list Get_LIBE_correlator(const vector<distr_t_list>& C_QED, const vector<distr_t_list>& C_SIB, const FLAV& U, const FLAV& D, string path);
distr_t_list Get_LIBE_correlator(const vector<distr_t_list>& C_QED, const vector<distr_t_list>& C_SIB, const FLAV& U, const FLAV& D, const distr_t& F,  string path);


#endif
