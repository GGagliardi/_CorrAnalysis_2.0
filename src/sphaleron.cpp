#include "../include/sphaleron.h"
#include "Corr_analysis.h"

using namespace std;

const int prec = 128;

void get_axion_therm_rate() {

  int Nboots=1000;
  int T=12;
  int Thalf= T/2  + 1;
  
  VVfloat C_raw(T);
  for(auto & c: C_raw) c.resize(Nboots);

  ifstream read("../sphaleron/correlator.dat");

  int Nrows= Nboots*Thalf;

  for(int irow=0;irow<Nrows;irow++) {

    int dummy;
    double dummy2,dummy3;
    double corr;

    read>>dummy>>dummy2>>corr>>dummy3;
    
    int t= irow%Thalf;
    int boot= irow/Thalf;

    C_raw[t][boot] = corr;
    C_raw[(T-t)%T][boot] = corr;


  }

  read.close();

  CorrAnalysis Corr(0, 100,Nboots);
  Corr.Nt= T;

  boost::filesystem::create_directory("../data/sphaleron");

  distr_t_list C= Corr.corr_t( C_raw, "../data/sphaleron/C");

  for(int t=0;t < Thalf;t++) {
    cout<<"C["<<t<<"] : "<<C.ave(t)<<" +- "<<C.err(t)<<endl;
  }

  Vfloat Cov, Corr_matrix;   

  Vfloat TT, RR; 

  for(int tt=0; tt <T; tt++)
    for(int rr=0;rr <T;rr++) {
      TT.push_back(tt);
      RR.push_back(rr);
      Cov.push_back(  C.distr_list[tt]%C.distr_list[rr] );
      Corr_matrix.push_back(  C.distr_list[tt]%C.distr_list[rr]/(C.err(tt)*C.err(rr)));
    }
  
  
  Print_To_File({}, {TT,RR,Corr_matrix, Cov}, "../data/sphaleron/Cov", "", "");

 

  


  double Ag_target= 1e-4;

  double syst, l_A;


  double k= 0.1308996938995747;
   
  
  double mult= 1e7;
  //important input par
  double aE0 = 0.2*k  ;
  double mean=  k;
  double sigma = 0.145833  ;


  const auto c_norm = [&mean, &sigma](double omega) -> double {

   
    return (omega-mean)/sinh( (omega-mean)/sigma); 
  };

  gsl_function_pp<decltype(c_norm)> Fp(c_norm);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function *G = static_cast<gsl_function*>(&Fp);
  double res_GSL, err_GSL;
  double pr = 1e-7;

  
  gsl_integration_qagiu(G, 0.0, 0.0, pr, 10000, w, &res_GSL, &err_GSL);
  gsl_integration_workspace_free (w);

  const auto func = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int Nj) -> PrecFloat {

    PrecFloat omega=E;
    return (omega-m)/sinh( (omega-m)/s); 
  };

  cout<<"norm of reco func: "<<res_GSL<<endl;

  distr_t ax_therm_rate = Get_Laplace_transfo(  mean,  sigma, aE0,  T, T/2 , prec, "k_1",func, C, syst, mult, l_A, "TANT", "STAGG", "ax_therm_rate", Ag_target,0, (1.0/res_GSL)*Get_id_distr(Nboots,0) , 0.0 , "ax_therm_rate", Cov, fake_func,0, fake_func_d ,  0, 2.0, 0.0); 


  



  return;
}
