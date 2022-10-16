#ifndef __R_ratio__
#define __R_ratio__

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


using namespace std;

void Get_Ergs_list();
void Compute_R_ratio();
void Get_exp_smeared_R_ratio(const Vfloat &Ergs_GeV_list_exp, double sigma);


#endif
