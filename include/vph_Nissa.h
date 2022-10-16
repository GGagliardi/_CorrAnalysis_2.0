#ifndef __vph_Nissa__
#define __vph_NissA__


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


void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg);
void Compute_form_factors_Nissa();


#endif
