#ifndef __scale_setting_main__
#define __scale_setting_main__


#include "g_minus_2_utilities.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "binary_io.h"
#include "stat.h"
#include "numerics.h"




using namespace std;

class scale_setting_info {

  

public:
  scale_setting_info() {}
  distr_t_list ms;
  distr_t_list ms_OS;
  vector<string> Ens;
  vector<string> Ens_l;
  vector<string> Ens_c;
  distr_t a_A;
  distr_t a_B;
  distr_t a_C;
  distr_t a_D;
  distr_t a_E;
  distr_t_list MK1;
  distr_t_list MK2;
  distr_t_list MDs1;
  distr_t_list Mpi;
  distr_t_list fpi;
  
  
  
};

scale_setting_info Get_scale_setting_info();


void Get_scale_setting() ;





#endif
