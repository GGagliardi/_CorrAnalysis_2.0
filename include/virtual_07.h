#ifndef __virtual_07__
#define __virtual_07__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "LatInfo.h"
#include "pt3_momenta.h"
#include "T_min.h"
#include "input.h"
#include "ChPT_form_factors.h"
#include "highPrec.h"
#include "Spectral.h"
#include "kernel.h"



class rt_07_Bs {

public:
  rt_07_Bs() : UseJack(1),  FA_T_u(1), FA_T_d(1), FV_T_u(1), FV_T_d(1),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) { };
  rt_07_Bs(bool x) : UseJack(x), FA_T_u(x), FA_T_d(x), FV_T_u(x), FV_T_d(x),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) {    };
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({ FA_T_u, FA_T_d, FV_T_u, FV_T_d}); return A[i];}
  Vfloat   Get_ch2(int i) { VVfloat A({ Ch2_FA_T_u, Ch2_FA_T_d, Ch2_FV_T_u, Ch2_FV_T_d}); return A[i];}
  

  bool UseJack;
  distr_t mass;
  distr_t fp;
  distr_t_list FA_T_u, FA_T_d, FV_T_u, FV_T_d;
  Vfloat  Ch2_FA_T_u, Ch2_FA_T_d,  Ch2_FV_T_u, Ch2_FV_T_d;
  int Nmeas;
  int Npars;
  int Ndof;
  bool Use_three_finest;
  bool Include_a4;
  int num_xg;
};


rt_07_Bs Get_virtual_tensor_FF(int num_xg, bool UseJack, int Njacks, string MESON,  string Corr_path, string path_out);


#endif
