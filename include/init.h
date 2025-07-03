#ifndef __init__
#define __init__


#include "3pts_mes_gamma_W.h"
#include "PionMassAnalysis.h"
#include "Axion_l7.h"
#include "PionMassAnalysis_twisted.h"
#include "PionMassAnalysis_twisted_adim.h"
#include "PionMassAnalysis_twisted_ov_X.h"
#include "PionMassAnalysis_ov_X.h"
#include "gm2.h"
#include "R_ratio.h"
#include "vph_Nissa.h"
#include "vph_Nissa_3d.h"
#include "tau_decay.h"
#include "tau_decay_strange.h"
#include "semileptonic.h"
#include "Bs_mumu_gamma.h"
#include "HVP.h"
#include "HVP_blinded_analysis.h"
#include "HVP_strange.h"
#include "RC_analysis.h"
#include "scale_setting_main.h"
#include "l7_Weinberg.h"
#include "tests.h"
#include "weak_annihilation_inclusive.h"
#include "RC_WI_analysis.h"
#include "Bs_phi_gamma.h"
#include "Kl4_HLT.h"
#include "Ds_phi_lnu.h"
#include "LIBE.h"
#include "RC_fits.h"
#include "sea_quark_effects.h"
#include "heavy_radiative.h"
#include "multi_shift_HLT.h"
#include "K3lnu.h"
#include "PINGU.h"
#include "axial_WI_disco.h"
#include "I0_gm2.h"
#include "sphaleron.h"
#include "tau_LIBE.h"
#include "inclusive_Ds_CKM_analysis.h"
#include "hc_ee.h"

using namespace std;

class MasterClass_analysis {

 public:
  MasterClass_analysis(string Input);


 private:
  void Analysis_manager();
  string Analysis_Mode;
  string Meson_to_analyze;
  bool IncludeDisconnected;
  string CURRENT_TYPE;


} ;




#endif
