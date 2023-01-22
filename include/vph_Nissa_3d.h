#ifndef __vph_Nissa_3d__
#define __vph_Nissa_3d__


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
#include "highPrec.h"
#include "Spectral.h"
#include "kernel.h"

using namespace std;

void Get_virt_list();
void Integrate_over_photon_insertion(const distr_t_list &W, vector<distr_t_list> &H, double Eg, int t_weak, double MP, int TO);
void Integrate_over_photon_insertion_w_subtraction(const distr_t_list &W, vector<distr_t_list> &H, double Eg, int t_weak, double MP, int ixg, int Tmin_mass, int Tmax_mass,  string out_path, string obs);
void GET_AXIAL_FORM_FACTORS_FROM_HADRONIC_TENSOR( distr_t_list &FA_distr_list, distr_t_list &H1_distr_list, distr_t_list &H2_distr_list, const distr_t_list &HA_11_distr_list, const distr_t_list &HA_33_distr_list, const distr_t_list &HA_03_distr_list, const distr_t_list &HA_30_distr_list, double kz, double Eg,const distr_t &MP_distr,const distr_t &FP_distr);
void Get_radiative_form_factors_3d();
void Compute_form_factors_Nissa_3d(double alpha, bool Integrate_Up_To_Infinite, double Emax, string SM_TYPE, bool Perform_theta_average, double E0_fact);

#endif
