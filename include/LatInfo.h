#ifndef __LatInfo__
#define __LatInfo__

#include "numerics.h"

using namespace std;


class LatticeInfo {

 public:
  LatticeInfo() {};
  LatticeInfo(string CURRENT_TYPE) { this->CURRENT_TYPE = CURRENT_TYPE;}
  void LatInfo(string S);
  void LatInfo_new_ens(string Tag);
  double ZP;
  double ZP_err;
  double ZS;
  double ZS_err;
  double Zv;
  double Zv_err;
  double Zv_M2;
  double Zv_M2_err;
  double Za;
  double Za_err;
  double Za_M2;
  double Za_M2_err;
  double a; //fm
  double a_err;  //fm
  double Beta;
  int L,T;
  double a_nucleon;
  double a_nucleon_err;
  double ml;
  string CURRENT_TYPE;
  double ainv; //GeV
  double ainv_err; //GeV
  double Za_fact;
  double Zm_fact;
  double Za_fact_err;
  double Zm_fact_err;
  double Za_fact_M2;
  double Zm_fact_M2;
  double Za_fact_M2_err;
  double Zm_fact_M2_err;
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
