#ifndef __init__
#define __init__


#include "3pts_mes_gamma_W.h"
#include "PionMassAnalysis.h"
#include "Axion_l7.h"
#include "PionMassAnalysis_twisted.h"
#include "gm2.h"

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
