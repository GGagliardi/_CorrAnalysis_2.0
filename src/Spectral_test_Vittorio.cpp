#include "../include/Spectral_test_Vittorio.h"




const bool UseJack=1;
const int Njacks=25;
const int Nboots=200;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 256;
const double sm_size= 0.05; // [ GeV ]
const double max_energy = 0.6; // [ GeV ]
const double Emin = 0.0001;
const double step_size= sm_size/3.0; //step size in energy in units of sm_size
using namespace std;




void Spectral_test_Vittorio() {


  boost::filesystem::create_directory("../data/toy_spectral");
  

  const auto f = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) -> PrecFloat {
		   
		 
		   
		   return Get_exact_gauss(E, m, s, E0);
		   
		   
		 };


  GaussianMersenne GM(5431);

  //test
  Vfloat corr_ave, corr_err;
  corr_ave = Read_From_File("../corr_Vittorio/corr.data", 1, 3);
  corr_err = Read_From_File("../corr_Vittorio/corr.data", 2, 3);

  //printV(corr_ave, "corr", 1);
  //printV(corr_err, "corr_err", 1);


  distr_t_list corr(UseJack);

  for(int t=0; t < (signed)corr_ave.size(); t++) {

    distr_t corr_t;

    for(int ijack=0; ijack < Njacks;ijack++) {

      corr_t.distr.push_back( corr_ave[t] + GM()*corr_err[t]/sqrt(Njacks -1.0)); 
    }

    corr.distr_list.push_back( corr_t);

  }

  int Nergs= (max_energy- Emin)/step_size;

  distr_t_list T_list(UseJack);
  distr_t_list S_list(UseJack);
  Vfloat Erg_list;

  
  for(int ie=0; ie<Nergs;ie++) {

    double Erg= Emin + (ie+1)*step_size;

    cout<<"Computing erg: "<<Erg<<endl;
    Erg_list.push_back(Erg);
    double s = sm_size;

    double syst_T, syst_S;
    double lambda_T, lambda_S;
    #pragma omp parallel for
    for(int i=0;i<2;i++) {
      if(i==0) {
	distr_t T;
	if(Erg < 0.20) {
	 T = Get_Laplace_transfo( Erg, s, Emin,  96, 46, prec, "G",f, corr, syst_T, 1 , lambda_T, "TANT", "fake", "toy", "fake" );
	}
	else {
	  T = Get_Laplace_transfo( Erg, s, Emin,  96, 46, prec, "G",f, corr, syst_T, 0.1 , lambda_T, "TANT", "fake", "toy", "fake" );
	}
	cout<<"T("<<Erg<<") : "<<T.ave()<<" +- "<<T.err()<<endl;
	T_list.distr_list.push_back(T);
      }
      else {
	distr_t S;
	if(Erg < 0.20) S = Get_Laplace_transfo( Erg, s, Emin,  96, 46, prec, "G",f, corr, syst_S, 0.15, lambda_S,  "SANF", "fake", "toy", "fake" );
	else S = Get_Laplace_transfo( Erg, s, Emin,  96, 46, prec, "G",f, corr, syst_S, 0.05, lambda_S,  "SANF", "fake", "toy", "fake" );
	cout<<"S("<<Erg<<") : "<<S.ave()<<" +- "<<S.err()<<endl;
	S_list.distr_list.push_back(S);
      }
    }

  }

  Print_To_File({}, {Erg_list, T_list.ave(), T_list.err(), S_list.ave(), S_list.err()}, "../data/toy_spectral/dens.dat","", "#E  T   S");

  

  
  

  return;
};
