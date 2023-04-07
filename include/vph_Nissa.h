#ifndef __vph_Nissa__
#define __vph_Nissa__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "pt3_momenta.h"
#include "T_min.h"
#include "input.h"
#include "Num_integrate_l4_decay_rate.h"
#include "virtual_ff_fit.h"
#include "virtual_FF_t_interval_list.h"
#include "ChPT_form_factors.h"

using namespace std;


class rt_FF {

public:
  rt_FF() : UseJack(1), FA(1), FV(1), FA_u(1), FV_u(1), FA_d(1), FV_d(1), Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) { };
  rt_FF(bool x) : UseJack(x), FA(x), FV(x), FA_u(x), FV_u(x), FA_d(x), FV_d(x), Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) {    };
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({FA,FV,FA_u,FV_u,FA_d,FV_d}); return A[i];}
  Vfloat   Get_ch2(int i) { VVfloat A({Ch2_FA, Ch2_FV, Ch2_FA_u, Ch2_FV_u, Ch2_FA_d, Ch2_FV_d}); return A[i];}

  bool UseJack;
  distr_t_list FA, FV, FA_u, FV_u, FA_d, FV_d;
  Vfloat Ch2_FA, Ch2_FV, Ch2_FA_u, Ch2_FV_u, Ch2_FA_d, Ch2_FV_d;
  int Nmeas;
  int Npars;
  int Ndof;
  bool Use_three_finest;
  bool Include_a4;
  int num_xg;
  
  
};

void Get_xg_t_list(int num_xg);
void Get_xg_to_spline();
void Get_xg_to_spline_VMD();
void Get_lattice_spacings_to_print();
void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens);
void Compute_form_factors_Nissa();
rt_FF Get_form_factors_Nissa(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list );


#endif
