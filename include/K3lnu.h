#ifndef __K3lnu__
#define __K3lnu__


#include "g_minus_2_utilities.h"
#include "gm2_fits.h"
#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "gm2.h"
#include "stat.h"
#include "binary_io.h"
#include "RC_WI_analysis.h"

using namespace std;

class An_cont_ret {
public:
  An_cont_ret() {}

  distr_t unsub_FF;
  distr_t meff_FF;
  distr_t mGS_FF;
  distr_t mRES_FF;

  distr_t_list unsub_FF_list;
  distr_t_list meff_FF_list;
  distr_t_list mGS_FF_list;
  distr_t_list mRES_FF_list;

  void Print(string out);
  
};

class exc_state_info {
public:
  exc_state_info() {}

  distr_t mGS;
  distr_t mRES;
};


An_cont_ret operator+(const An_cont_ret &A, const An_cont_ret &B);
An_cont_ret operator-(const An_cont_ret &A, const An_cont_ret &B);
An_cont_ret operator+(const An_cont_ret &A, const distr_t &B);
An_cont_ret operator-(const An_cont_ret &A, const distr_t &B);
An_cont_ret operator*(const double &a, const An_cont_ret &A);
An_cont_ret operator*(const An_cont_ret &A, const double &a);
An_cont_ret operator+(const distr_t &A, const An_cont_ret &B);
An_cont_ret operator-(const distr_t &A, const An_cont_ret &B);
An_cont_ret Analytic_continuation(const distr_t_list& C_in, const distr_t& E, const exc_state_info& ex_info, int tw, int TO,  string out) ;



void K3lnu();
void Get_electrounquenching();
void Analyze_Aprime();



#endif
