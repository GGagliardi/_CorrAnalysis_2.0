#ifndef __RC_WI_analysis__
#define __RC_WI_analysis__


#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "gm2.h"
#include "stat.h"
#include "binary_io.h"

class RCs_info {

public:
  RCs_info()  {}
  
  distr_t_list Za;
  distr_t_list Zv;
  distr_t_list Zp_ov_Zs;
  vector<string> Ens;

};


RCs_info Get_RCs(string AN_TYPE);
void RC_WI_analysis();


#endif
