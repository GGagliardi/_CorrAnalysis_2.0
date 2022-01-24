#ifndef __Spectral__
#define __Spectral__

#include "numerics.h"
#include "highPrec.h"
#include "random.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "Corr_analysis.h"

using namespace std;

PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_gaussE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_lego(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_legoE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_cauchy(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0);
PrecFloat Get_exact_func(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0, string SMEARING_FUNC);
PrecFloat Get_gaussE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0);
PrecFloat Get_legoE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0);
PrecFloat BaseFunc(const PrecFloat& E, int t, int T);
PrecFloat aE0(const PrecFloat &E0, int t);
PrecFloat F_gauss(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_gaussE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_lego(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
PrecFloat F_legoE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
//PrecFloat F_cauchy(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t);
void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax) ;
void Get_Rt(PrecVect& Rt, const PrecFloat &E0,   int T, int tmin, int tmax);
void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int  T, int tmin, int tmax, string SMEARING_FUNC) ;
void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax);
PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, string SMEARING_FUNC);
void Get_Laplace_transfo(double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC);


#endif
