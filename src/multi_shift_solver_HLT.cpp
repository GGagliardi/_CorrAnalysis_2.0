#include "../include/multi_shift_HLT.h"
#include "highPrec.h"

using namespace std;
const int N = 200;

void conj_grad_solver_multi_shift(const PrecMatr A, const PrecVect y, const Vfloat& m_n, PrecVect y0 ) {

  //solves using conjugate gradient the system: [ A + m_n ] x = y, using y0 as guess

  
  


}


void multi_shift_solver_HLT() {

  //set precision
  PrecFloat::setDefaultPrecision(54*8);
  
  //generate A matrix
  PrecMatr Atr;
  PrecFloat E0=PrecFloat(0.01);
  Get_Atr_std(Atr, E0, 2*N, 1, N);

  //generate fake covariance e^-(x1-x2)*5
  PrecMatr B;
  B.resize(N,N);
  for(int n=0;n<N;n++) {
    PrecFloat C= exp(-n*PrecFloat(0.02));
    for(int m=n;m<N;m++) { B(n,m) = C*exp(-(m-n)*PrecFloat(2.0)) ; B(m,n) = B(n,m) ; }
  }
    
    


   



  return;

}
