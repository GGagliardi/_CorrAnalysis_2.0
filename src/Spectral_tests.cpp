#include "../include/Spectral_tests.h"

using namespace std;


void Spectral_tests() {

  double mean = 1;
  double Estart = 0.1;
  double sigma = 0.03;
  int tmax= 256;
  int T = 512;
  int prec = 4096*2;

  Get_Laplace_transfo(  mean,  sigma, Estart,  T, tmax, prec, "GAUSSIAN_E2") ;





  return;
}
