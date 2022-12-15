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

void Get_xg_t_list();
void Get_xg_to_spline();
void Get_xg_to_spline_VMD();
void Get_lattice_spacings_to_print();
void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens);
void Compute_form_factors_Nissa();


#endif
