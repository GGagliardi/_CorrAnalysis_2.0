#ifndef __I0_gm2__
#define __I0_gm2__


#include "g_minus_2_utilities.h"
#include "gm2_fits.h"
#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "gm2.h"
#include "stat.h"
#include "binary_io.h"
#include "HVP.h"


void I0_gm2();
void Bounding_HVP_isoscalar(distr_t &HVP, int &Tcut_opt, const distr_t_list &V, const distr_t &a, string path,  distr_t lowest_mass);


#endif
