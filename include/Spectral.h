#ifndef __Spectral__
#define __Spectral__

#include "numerics.h"
#include "highPrec.h"
#include "random.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "Corr_analysis.h"
#include "LatInfo.h"

using namespace std;

double Get_exact_gauss(const double &E,const double &m,const double &s,const double &E0);
PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_gaussE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_lego(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_legoE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_cauchy(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_func(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f);
PrecFloat Get_gaussE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0);
PrecFloat Get_legoE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0);
PrecFloat BaseFunc(const PrecFloat& E, int t, int T);
PrecFloat aE0(const PrecFloat &E0, int t, int n);
PrecFloat F_gauss(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_gaussE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_lego(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_legoE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
//PrecFloat F_cauchy(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax) ;
void Get_Rt(PrecVect& Rt, const PrecFloat &E0,   int T, int tmin, int tmax);
void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int  T, int tmin, int tmax, string SMEARING_FUNC,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f) ;
void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax);
PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f);
PrecFloat Get_M2(PrecFloat &m, PrecFloat &s, PrecFloat &E0, const function<PrecFloat(const PrecFloat&,const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f);
void Get_Rt_up_to_N(PrecFloat &E0, int T, int tmin, int tmax, vector<PrecVect> &Rt_n);
void Get_G_matrix(PrecMatr &G,const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n);
void Get_M_N(PrecFloat &m, PrecFloat &s, PrecFloat &E0,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f, PrecVect &M_n);
void Get_M_tilde_N( const PrecVect &ft, const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n,  PrecVect &M_tilde_n);
void Compute_covariance_matrix(PrecMatr &B,const PrecMatr &Atr, const distr_t_list &corr, int tmin, int tmax, string MODE);
void Get_optimal_lambda(const PrecMatr &Atr,const PrecMatr &B,const PrecVect &ft,const PrecVect &Rt,const PrecFloat & M2,const double &mean, const double &sigma, const double &Estart,  double& lambda_opt, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f, vector<PrecVect> Rt_n, const PrecVect &M_n ,const distr_t_list & corr,int T, int tmin, int tmax, const double mult,   string MODE, string curr_type, string SMEARING_FUNC, string CORR_NAME, string FLAV);
distr_t Get_Laplace_transfo(double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f, const distr_t_list& corr, double &syst,const double mult, double& lambda_ret, string MODE, string cur_type, string CORR_NAME, string FLAV);


#endif
