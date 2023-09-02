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
  rt_07_Bs() : UseJack(1),  F_T_u(1), F_T_d(1),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) { };
  rt_07_Bs(bool x) : UseJack(x),  F_T_u(x), F_T_d(x),  Nmeas(0), Npars(0), Ndof(0), Use_three_finest(0), Include_a4(0), num_xg(0) {    };
  distr_t_list Get_FF(int i) { vector<distr_t_list> A({ F_T_u, F_T_d}); return A[i];}
  Vfloat   Get_ch2(int i) { VVfloat A({  Ch2_F_T_u, Ch2_F_T_d}); return A[i];}
  

  bool UseJack;
  distr_t mass;
  distr_t fp;
  distr_t_list F_T_u, F_T_d;
  Vfloat  Ch2_F_T_u, Ch2_F_T_d;
  int Nmeas;
  int Npars;
  int Ndof;
  bool Use_three_finest;
  bool Include_a4;
  int num_xg;
};


rt_07_Bs Get_virtual_tensor_FF(int num_xg, bool UseJack, int Njacks, string MESON,  string Corr_path, string path_out);


#endif
