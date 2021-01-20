#ifndef __corranalysis__
#define __corranalysis__

#include "numerics.h"
#include "stat.h"


using namespace std;

class CorrAnalysis {

 public:
  CorrAnalysis() {Perform_Nt_t_average=true; seed = 4324; UseJack=1; Tmin=0; Tmax=0; Nt=0; Reflection_sign=1;}
 CorrAnalysis(bool UseJack, int Njacks, int nboots) : UseJack(UseJack), Njacks(Njacks), Nboots(nboots) { Perform_Nt_t_average=true; seed=4324; Tmin=0; Tmax=0; Nt=0; Reflection_sign=1;} 
  Pfloat Fit_(const distr_t_list& data_t_distr);
  distr_t Fit_distr(const distr_t_list& data_t_distr);
  distr_t_list effective_mass_t(const VVfloat &corr_A, string Obs); //returns jackknife or bootstrap distribution of the effective mass m(t) "acts on data"
   distr_t_list effective_mass_t(const distr_t_list& corr_A_distr, string Obs); //returns jackknife or bootstrap distribution of the effective mass m(t) "acts on distribution"
  distr_t_list effective_slope_t(const VVfloat &corr_A,const VVfloat &corr_B, string Obs); //returns jackknife or bootstrap distribution of the effective slope dm(t) "acts on  data"
  distr_t_list effective_slope_t(const distr_t_list& corr_A_distr,const distr_t_list& corr_B_distr, string Obs);//returns jackknife or bootstrap distribution of the effective slope dm(t) "acts on distributions"
  distr_t_list effective_slope2_t(const VVfloat &corr_A,const VVfloat &corr_B, string Obs); //returns jackknife or bootstrap distribution of the effective slope squared dm2(t) "acts on  data"
  distr_t_list corr_t(const VVfloat &corr, string Obs); //returns jackknife or bootstrap distribution of a correlator
  distr_t_list residue_t(const VVfloat &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole "acts on data"
  distr_t_list residue_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of the residue of the pole "acts on distributions"
  distr_t_list decay_constant_t(const VVfloat &corr, string Obs); //returns jackknife or boostrap distribution of a decay constant "acts on data"
  distr_t_list decay_constant_t(const distr_t_list &corr_A, string Obs); //returns jackknife or boostrap distribution of a decay constant "acts on distribution"
  Vfloat ASymm(const VVfloat& data, int t );


  bool UseJack;
  int Njacks;
  int Nboots;
  int Tmin;
  int Tmax;
  int Nt;
  bool Perform_Nt_t_average;
  int Reflection_sign;
  int seed;
;
  


};
  

























#endif
