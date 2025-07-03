#ifndef __tau_decay_LIBE_ISO__
#define __tau_decay_LIBE_ISO__

#include "numerics.h"
#include "highPrec.h"
#include "Spectral.h"
#include "random.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "Corr_analysis.h"
#include "LatInfo.h"
#include "input.h"
#include "Meson_mass_extrapolation.h"
#include "g_minus_2_utilities.h"
#include "binary_io.h"



using namespace std;


class LIBE_tau_ret {

  

public:
  LIBE_tau_ret() {}


  distr_t A0_TM;
  distr_t V0_TM;
  distr_t VK_TM;
  distr_t AK_TM;

  distr_t A0_OS;
  distr_t V0_OS;
  distr_t VK_OS;
  distr_t AK_OS;
    
};


double Customized_plateaux_tau_spectre_LIBE_ISO( double alpha, double Emax, string channel, string reg, double s, string Ens );
LIBE_tau_ret Compute_tau_decay_width_LIBE_ISO(bool Is_Emax_Finite, double Emax, double alpha, double sigma, int Njacks, string Ens);
void Generate_data();

#endif
