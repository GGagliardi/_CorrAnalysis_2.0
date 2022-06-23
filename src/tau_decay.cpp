#include "../include/tau_decay.h"


const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const bool UseJack=1;
const int Njacks=50;
const int Nboots=200;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 256;
const double Nc=3;
const double m_mu= 0.10565837; // [ GeV ]
const string SM_TYPE_0= "TAU_KERNEL_0";
const string SM_TYPE_1= "TAU_KERNEL_1";
const double GF= 1.1663787*1e-5; //[GeV^-2]
//CKM matrix elements
const double Vud = 0.97373;
const double m_tau= 1.77686;
const double m_Jpsi= 3.0969;
const double m_Jpsi_err= 0.001;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_etac= 2.9839;
const double m_etac_err = 0.004;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const double m_pi_plus = 0.13957039;
const string Extrapolation_strange_mode="etas";
const string Extrapolation_charm_mode  = "etac";
const bool Include_light_disconnected=true;
const bool Include_strange_disconnected=true;
const bool Include_charm_disconnected=true;
const bool Include_off_diagonal_disconnected=true;
const int Simps_ord=4;
const double E0_l = m_pi_plus/m_tau;
const Vfloat sigma_list({0.05, 0.1, 0.2, 0.25, 0.3});
const double C_V = 2*M_PI;
const double GAMMA_FACT= pow(Vud*GF,2)*pow(m_tau,3)/(16*pow(M_PI,2))*1e12;
const string MODE="SANF";
using namespace std;


void Compute_tau_decay_width() {


  //create directories

  boost::filesystem::create_directory("../data/tau_decay");
  boost::filesystem::create_directory("../data/tau_decay/light");
  boost::filesystem::create_directory("../data/tau_decay/light/Br");
  boost::filesystem::create_directory("../data/tau_decay/light/corr");




 




  //Read data


  //Custom sorting for V_light to account for the two replica r0 and r1
  auto Sort_light_confs = [](string A, string B) {
		      
			    string rA = A.substr(A.length()-2);
			    string rB = B.substr(B.length()-2);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(A.substr(A.length()-1));
			      int n2 = stoi(B.substr(B.length() -1));
			      if(rA==rB) return A<B;
			      else return n1<n2;
			    }
		      
			    return A<B;
			  };



  data_t Vk_data_tm, V0_data_tm, Ak_data_tm, A0_data_tm;
  data_t Vk_data_OS, V0_data_OS, Ak_data_OS, A0_data_OS;


  //light 
  //tm
  Vk_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  V0_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "V0V0", Sort_light_confs);
  Ak_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "AKAK", Sort_light_confs);
  A0_data_tm.Read("../tau_decay_data/light", "mes_contr_2pts_ll_1", "A0A0", Sort_light_confs);
  //OS
  Vk_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);
  V0_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "V0V0", Sort_light_confs);
  Ak_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "AKAK", Sort_light_confs);
  A0_data_OS.Read("../tau_decay_data/light", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs);









  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
      
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      
    }
  }


  //############################################################################################




  int Nens= Vk_data_tm.size;


  vector<distr_t_list> Br_tau_tm, Br_tau_OS;

  for(int iens=0; iens<Nens; iens++) { Br_tau_tm.emplace_back( UseJack, sigma_list.size()); Br_tau_OS.emplace_back( UseJack, sigma_list.size());}



  //loop over the ensembles

  for(int iens=0; iens<Nens;iens++) {


 
     
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = Vk_data_tm.nrows[iens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<Vk_data_tm.Tag[iens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(Vk_data_tm.Tag[iens].substr(1,1)=="A") {a_distr=a_A; Zv = ZV_A; Za = ZA_A;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+Vk_data_tm.Tag[iens]+" not found");



    // lambda function to be used as a smearing func.
  
    const auto K0 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) -> PrecFloat {

		      PrecFloat X= E/(m_tau*a_distr.ave());
		      PrecFloat sm_theta = 1/(1+ exp(-2*(1-X)/s));
		      if( X < E0) return 0.0;
		      else return (1/X)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    const auto K1 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) -> PrecFloat {

		      PrecFloat X = E/(m_tau*a_distr.ave());
		      PrecFloat sm_theta = 1/(1+ exp(-2*(1-X)/s));
		      if( X < E0) return 0.0;
		      else return (1 + 2*pow(X,2))*(1/(X))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);

    //tm
    distr_t_list Vk_tm_distr, V0_tm_distr, Ak_tm_distr, A0_tm_distr;
    //OS
    distr_t_list Vk_OS_distr, V0_OS_distr, Ak_OS_distr, A0_OS_distr;
   
    //light-tm sector
    Vk_tm_distr = Corr.corr_t(Vk_data_tm.col(0)[iens], "../data/tau_decay/light/corr/Vk_tm");
    V0_tm_distr = Corr.corr_t(V0_data_tm.col(0)[iens], "../data/tau_decay/light/corr/V0_tm");
    Ak_tm_distr = Corr.corr_t(Ak_data_tm.col(0)[iens], "../data/tau_decay/light/corr/Ak_tm");
    A0_tm_distr = Corr.corr_t(A0_data_tm.col(0)[iens], "../data/tau_decay/light/corr/A0_tm");

    //light-OS sector
    Vk_OS_distr = Corr.corr_t(Vk_data_OS.col(0)[iens], "../data/tau_decay/light/corr/Vk_OS");
    V0_OS_distr = Corr.corr_t(V0_data_OS.col(0)[iens], "../data/tau_decay/light/corr/V0_OS");
    Ak_OS_distr = Corr.corr_t(Ak_data_OS.col(0)[iens], "../data/tau_decay/light/corr/Ak_OS");
    A0_OS_distr = Corr.corr_t(A0_data_OS.col(0)[iens], "../data/tau_decay/light/corr/A0_OS");


    distr_t_list C0_tm, Cii_tm, C0_OS, Cii_OS;

    
    //######### DEFINE 0th and ii component of C^munu ###########
    //tm
    C0_tm = C_V*(0.0*Za*Za*V0_tm_distr + Zv*Zv*A0_tm_distr);
    Cii_tm = C_V*(Za*Za*Vk_tm_distr + Zv*Zv*Ak_tm_distr);
    //OS
    C0_OS = C_V*(0.0*Zv*Zv*V0_OS_distr + Za*Za*A0_OS_distr);
    Cii_OS = C_V*(Zv*Zv*Vk_OS_distr + Za*Za*Ak_OS_distr);
    //###########################################################

    bool Found_error_less_x_percent=false;
    double x=5;
    //tm
    int tmax_tm_0=1;
    while(!Found_error_less_x_percent && tmax_tm_0 < Corr.Nt/2 -1 ) {
   
      if( (C0_tm.distr_list[tmax_tm_0]).err()/fabs( (C0_tm.distr_list[tmax_tm_0]).ave()) <  0.01*x) tmax_tm_0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_tm_1=1;
    while(!Found_error_less_x_percent && tmax_tm_1 < Corr.Nt/2 -1 ) {
   
      if( (Cii_tm.distr_list[tmax_tm_1]).err()/fabs( (Cii_tm.distr_list[tmax_tm_1]).ave()) <  0.01*x) tmax_tm_1++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    //OS

    int tmax_OS_0=1;
    while(!Found_error_less_x_percent && tmax_OS_0 < Corr.Nt/2 -1 ) {
   
      if( (C0_OS.distr_list[tmax_OS_0]).err()/fabs( (C0_OS.distr_list[tmax_OS_0]).ave()) <  0.01*x) tmax_OS_0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_OS_1=1;
    while(!Found_error_less_x_percent && tmax_OS_1 < Corr.Nt/2 -1) {
   
      if( (Cii_OS.distr_list[tmax_OS_1]).err()/fabs( (Cii_OS.distr_list[tmax_OS_1]).ave()) <  0.01*x) tmax_OS_1++;
      else Found_error_less_x_percent=true;
    }

   

    //loop over sigma
    #pragma omp parallel for
    for(int is=0; is < (signed)sigma_list.size(); is++) {

      distr_t Br_sigma_0_tm;
      distr_t Br_sigma_1_tm;
      distr_t Br_sigma_0_OS;
      distr_t Br_sigma_1_OS;
      double s= sigma_list[is];
      //int tmax= T/2 -4;
      double l0_tm, l1_tm, l0_OS, l1_OS;
      double syst_0_tm, syst_1_tm, syst_0_OS, syst_1_OS;
      double mult=0.05;
      if(MODE=="SANF") mult= 0.00005;
      distr_t resc_GeV = 1.0/(a_distr*a_distr*a_distr*m_tau);
      cout<<"Print times"<<endl;
      cout<<"tmax_tm_0: "<<tmax_tm_0<<" tmax_tm_1: "<<tmax_tm_1<<" tmax_OS_0: "<<tmax_OS_0<<" tmax_OS_1: "<<tmax_OS_1<<endl;
      Br_sigma_0_tm = GAMMA_FACT*resc_GeV*Get_Laplace_transfo(  0.0,  s, E0_l,  T, tmax_tm_0, prec, SM_TYPE_0,K0, -1*C0_tm, syst_0_tm, mult, l0_tm, MODE, "tm", "TAU_K0_light_tm_"+Vk_data_tm.Tag[iens], "ud" );
      cout<<"Br0_tm["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<l0_tm<<" : "<<Br_sigma_0_tm.ave()<<" +- "<<Br_sigma_0_tm.err()<<endl; 
      Br_sigma_1_tm = GAMMA_FACT*resc_GeV*Get_Laplace_transfo(  0.0,  s, E0_l,  T, tmax_tm_1, prec, SM_TYPE_1,K1, Cii_tm, syst_1_tm, mult, l1_tm, MODE, "tm", "TAU_K1_light_tm_"+Vk_data_tm.Tag[iens], "ud" );
      cout<<"Br1_tm["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<l1_tm<<" : "<<Br_sigma_1_tm.ave()<<" +- "<<Br_sigma_1_tm.err()<<endl; 
      Br_sigma_0_OS = GAMMA_FACT*resc_GeV*Get_Laplace_transfo(  0.0,  s, E0_l,  T, tmax_OS_0, prec, SM_TYPE_0,K0, -1*C0_OS, syst_0_OS, mult, l0_OS, MODE, "OS", "TAU_K0_light_OS_"+Vk_data_tm.Tag[iens], "ud" );
      cout<<"Br0_OS["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<l0_OS<<" : "<<Br_sigma_0_OS.ave()<<" +- "<<Br_sigma_0_OS.err()<<endl; 
      Br_sigma_1_OS = GAMMA_FACT*resc_GeV*Get_Laplace_transfo(  0.0,  s, E0_l,  T, tmax_OS_1, prec, SM_TYPE_1,K1, Cii_OS, syst_1_OS, mult, l1_tm, MODE, "OS", "TAU_K1_light_OS_"+Vk_data_tm.Tag[iens], "ud" );
      cout<<"Br1_OS["<<Vk_data_tm.Tag[iens]<<"] , s= "<<s<<", l= "<<l1_OS<<" : "<<Br_sigma_1_OS.ave()<<" +- "<<Br_sigma_1_OS.err()<<endl;

      
      distr_t Br_sigma_tm = Br_sigma_1_tm + Br_sigma_0_tm;
      distr_t Br_sigma_OS = Br_sigma_1_OS + Br_sigma_0_OS;
      Br_tau_tm[iens].distr_list[is] = Br_sigma_tm;
      Br_tau_OS[iens].distr_list[is] = Br_sigma_OS;
    }
    #pragma omp barrier
  }

    
   





  //Print to File
  for(int iens=0; iens<Nens;iens++) Print_To_File({}, {sigma_list, Br_tau_tm[iens].ave(), Br_tau_tm[iens].err(), Br_tau_OS[iens].ave(), Br_tau_OS[iens].err()}, "../data/tau_decay/light/Br/br_"+MODE+"_"+Vk_data_tm.Tag[iens]+".dat", "", "#sigma Br[tm] Br[OS]  [GeV]");


    
 
  


  return;

}
