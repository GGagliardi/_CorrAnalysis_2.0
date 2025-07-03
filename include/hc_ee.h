#ifndef __hc_ee__
#define __hc_ee__

#include "numerics.h"
#include "Corr_analysis.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "stat.h"
#include "binary_io.h"
#include "input.h"

using namespace std;


void Get_ee_plateaux_int_2pt_H_loc(string Ens, CorrAnalysis &Corr);
void Get_ee_plateaux_int_2pt_H_sm(string Ens, CorrAnalysis &Corr);


void hc_ee();


#endif
