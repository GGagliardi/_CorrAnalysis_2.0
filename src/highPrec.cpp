#include "../include/highPrec.h"


using namespace std;

PrecFloat ExpEiComplexSum(PrecFloat MOD, PrecFloat PH, PrecFloat s,  bool MODE) {

  PrecFloat gamma= precEuler();
  PH= atan(PH);
  PrecFloat RE_Log= log(MOD);
  PrecFloat IM_Log= PH;

  //evaluate sum_k=1,infty  z^k /( k * k! )

  PrecFloat START_SUM, SUM_OLD;
  if(MODE==0) START_SUM= (RE_Log+gamma)*sin(s) -IM_Log*cos(s);
  else START_SUM= (RE_Log+gamma)*cos(s) +IM_Log*sin(s);

  PrecFloat SUM=START_SUM;
  SUM_OLD=0.0;

  
  
  bool Arg_zero = (MOD==0);
  
  if(!Arg_zero) {
  
    bool converged=false;

    unsigned long int k=1;

    PrecFloat MOD_TERM_k= MOD;
    SUM+=MOD*((MODE==0)?sin(s- PH):cos(s - PH));
    k++;
   
    while( ! converged) {

      if(k%30==0) SUM_OLD= SUM;

      MOD_TERM_k *=  MOD*(k-1)/(PrecFloat(k)*PrecFloat(k));
      SUM += MOD_TERM_k*((MODE==0)?sin(s- k*PH):cos(s - k*PH));
    
  
      //cout<<"SUM: k="<<k<<" : "<< SUM<<endl;

      k++;
      if(SUM_OLD == SUM) converged=true;
      
    }
  }


  return SUM;
  
}
