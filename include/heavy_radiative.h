#ifndef __heavy_radiative__
#define __heavy_radiative__

#include "numerics.h"
#include "Corr_analysis.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "stat.h"
#include "binary_io.h"
#include "input.h"

using namespace std;


void Get_plateaux_int_2pt_H_loc(string Ens, CorrAnalysis &Corr);
void Get_plateaux_int_2pt_H_sm(string Ens, CorrAnalysis &Corr);
void Get_plateaux_int_2pt_dH_loc(string Ens, string M, CorrAnalysis &Corr);
void Get_plateaux_int_2pt_dH_sm(string Ens, string M, CorrAnalysis &Corr);
void Get_plateaux_VEV_ratio(string Ens, string M, CorrAnalysis &Corr);

void heavy_radiative();
void Get_chi_c1_decay();


#endif
