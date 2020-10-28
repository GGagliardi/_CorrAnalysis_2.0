#ifndef __LatInfo__
#define __LatInfo__

#include "numerics.h"

using namespace std;


class LatticeInfo {

 public:
  LatticeInfo() {};
  LatticeInfo(string CURRENT_TYPE) { this->CURRENT_TYPE = CURRENT_TYPE;}
  void LatInfo(string S);
  double ZP;
  double Zv;
  double Zv_err;
  double Zv_M2;
  double Zv_M2_err;
  double Za;
  double Za_err;
  double Za_M2;
  double Za_M2_err;
  double a;
  double a_err;
  double Beta;
  string CURRENT_TYPE;
  Pfloat Retrieve_Zv(string S, int ibranch) {
    LatInfo(S);
    if(ibranch<4) {return make_pair(Zv,Zv_err);}
    else return make_pair(Zv_M2, Zv_M2_err);
  }
  Pfloat Retrieve_Za(string S, int ibranch) {
    LatInfo(S);
    if(ibranch<4) {return make_pair(Za, Za_err);}
    else return make_pair(Za_M2, Za_M2_err);
  }


};


void Read_pars_from_ensemble_tag(vector<string>&, Vfloat& m_lat, Vfloat &L, Vfloat &Nt);
void Read_pars_from_ensemble_tag(string, double&, double&, double&);
void ReadBranch(int k, Eigen::MatrixXd& CovMatrixInput, Eigen::VectorXd& Ave_input_parameters);

  


















#endif
