#ifndef __Num_integrate_l4_decay_rate__
#define __Num_integrate_l4_decay_rate__



#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "T_min.h"
#include "input.h"
#include "kernel.h"
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_monte.h>
//#include <gsl/gsl_monte_plain.h>
//#include <gsl/gsl_monte_miser.h>
//#include <gsl/gsl_monte_vegas.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/naive_monte_carlo.hpp>

using namespace std;


class Decay_Rate_Integration_Result {

 public:
  Decay_Rate_Integration_Result() {}
  Decay_Rate_Integration_Result(int UseJack) : Jack_Distr_Int_Quad_mumu(UseJack), Jack_Distr_Int_Quad_ee(UseJack), Jack_Distr_Int_MonteCarlo_mumu(UseJack), Jack_Distr_Int_MonteCarlo_ee(UseJack) {}
  

  double Int_Quad_val_mumu, Int_Quad_val_ee, Int_MonteCarlo_val_mumu, Int_MonteCarlo_val_ee;

  double Int_Quad_err_mumu, Int_Quad_err_ee, Int_MonteCarlo_err_mumu, Int_MonteCarlo_err_ee;

  //int Nsteps_MonteCarlo_mumu, Nsteps_MonteCarlo_ee;

  //int Nfunc_eval_Quad_mumu, Nfunc_eval_Quad_elel;

  double eps_rel_mumu, eps_rel_ee;

  Vfloat Stat_err_MonteCarlo_mumu, Stat_err_MonteCarlo_ee;

  bool Exit_status;

  distr_t Jack_Distr_Int_Quad_mumu, Jack_Distr_Int_Quad_ee, Jack_Distr_Int_MonteCarlo_mumu, Jack_Distr_Int_MonteCarlo_ee;

};


Decay_Rate_Integration_Result Num_Integrate_Decay_Rate(const vector<function<double(double, double)>> &H1, const vector<function<double(double, double)>> &H2  , const vector<function<double(double, double)>> &FA, const vector<function<double(double, double)>> &FV, distr_t &m_distr, distr_t & fp_distr, bool UseJack);






#endif
