#ifndef __virtual_ff_fit__
#define __virtual_ff_fit__



#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "Bootstrap_fit.h"
#include "header_file_virph.h"
#include "LatInfo.h"
#include "T_min.h"
#include "pt3_momenta.h"


void Fit_virtual_FF_VMD(vector<function<double(double, double)>> &fit_func, const vector<distr_t> &FF, distr_t &f_p, distr_t &m_p, distr_t &Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag,string Meson, bool UseJack, bool ConstFit, int t_fit);
void Fit_virtual_FF_ChPT(vector<function<double(double, double)>> &fit_func,const vector<distr_t> &FF, distr_t &f_p, distr_t &m_p, distr_t &Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag,string Meson, bool UseJack, bool ConstFit, int t_fit);




#endif 
