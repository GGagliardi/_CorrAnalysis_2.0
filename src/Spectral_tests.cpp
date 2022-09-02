#include "../include/Spectral_tests.h"



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
const int prec = 512;
const int prec_charm=512;
const double rho_R= 12*M_PI*M_PI;
const double Nc=3;
const double Rpert= Nc*( qu*qu + qd*qd);
const double Rpert_strange = Nc*(qs*qs);
const double Rpert_charm = Nc*(qc*qc);
const double m_mu= 0.10565837; // [ GeV ]
const string SM_TYPE= "GAUSSIAN";
const bool add_pert_corr_strange_up_to = 1.0;
const bool add_pert_corr_charm_up_to = 1.0;
const double sm_size= 1.0; // [ GeV ]
const double max_energy = 3.2; // [ GeV ]
const double Emin = 2.0*MPiPhys;
const double step_size= 0.2; //step size in energy in units of sm_size
const double m_Jpsi= 3.0969;
const double m_Jpsi_err= 0.001;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_etac= 2.9839;
const double m_etac_err = 0.004;
const double m_etas = 0.68989;
const double m_etas_err= 0.00050;
const string Extrapolation_strange_mode="etas";
const string Extrapolation_charm_mode  = "etac";
const int pert_corr_strange_on_off = 0;
const int pert_corr_charm_on_off = 0;
const int pert_corr_light_on_off = 0;
const bool Include_light_disconnected=false;
const bool Include_strange_disconnected=false;
const bool Include_charm_disconnected=false;
const bool Include_off_diagonal_disconnected=false;
const bool Compute_experimental_smeared_R_ratio=false;
const int Sim_ord=4;
using namespace std;





void Get_exp_smeared_R_ratio(const Vfloat& Ergs_GeV_list, double sigma) {


  const auto f_R = [](const double &E, const double &m, const double &s, const double &E0) -> double {
		     
		     if(SM_TYPE == "WIN")  return (1.0/(1 + exp(-2*(E-(m-s))/s)) - 1.0/(1 + exp(-2*(E-(m+s))/s)))/(2*s);
		     
		     else if (SM_TYPE == "GAUSSIAN") return Get_exact_gauss(E, m, s, E0);
		   
		     else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
		   
		     return 0.0;
		   
		   };  

  //load File
  Vfloat E = Read_From_File("../dispersive_data/ratio.dati", 0,5);
  Vfloat RE_exp = Read_From_File("../dispersive_data/ratio.dati", 1,5);
  Vfloat RE_exp_stat = Read_From_File("../dispersive_data/ratio.dati", 2,5);
  Vfloat RE_exp_syst = Read_From_File("../dispersive_data/ratio.dati", 3,5);

  //define intervals around resonances where to sum absolute systematic error
  Vfloat Res_MIN({ 0.40, 3.09, 3.74, 9.45, 10.021, 10.352, 10.54}); 
  Vfloat Res_MAX({ 1.06, 3.10, 3.81, 9.47, 10.026, 10.358, 10.62});

  //assign E values to interval  [ Res_MIN[i]: Res_MAX[i] ]
  int Nsteps= E.size();
  vector<int> E_intervals(Nsteps, -1);
  for(int istep=0;istep<Nsteps;istep++) {

    for(int j=0;j<(signed)Res_MIN.size();j++) {
      if( E[istep] >= Res_MIN[j] && E[istep] <= Res_MAX[j]) {
	E_intervals[istep] = j;
	break;
      }
      if( E[istep] < Res_MIN[j]) break;
      if( E[istep] > Res_MAX[Res_MAX.size() -1]) break;
    }

  }


  //start integrating using Simpson Rule 8/3
  
  double dE = 0.00025; //GeV
  Vfloat RE_sm_exp_ave(Ergs_GeV_list.size(), 0.0);
  Vfloat RE_sm_exp_stat(Ergs_GeV_list.size(),0.0);
  Vfloat RE_sm_exp_syst(Ergs_GeV_list.size(),0.0);
  Vfloat RE_sm_exp_err(Ergs_GeV_list.size(), 0.0);

  for(int ierg=0;ierg<(signed)Ergs_GeV_list.size();ierg++) {

    Vfloat DE_corr(Res_MIN.size(),0.0);

    //Simpson integration
    for(int istep=0;istep<Nsteps;istep++) {

      double sm_coeff = f_R(E[istep], Ergs_GeV_list[ierg],sigma, Emin);
      RE_sm_exp_ave[ierg] += w(istep,Sim_ord)*dE*RE_exp[istep]*sm_coeff;
      RE_sm_exp_stat[ierg] += pow(w(istep,Sim_ord)*dE*RE_exp_stat[istep]*sm_coeff,2);
      if( E_intervals[istep] != -1) DE_corr[ E_intervals[istep]] += fabs(w(istep,Sim_ord)*dE*RE_exp_syst[istep]*sm_coeff);
      else RE_sm_exp_syst[ierg] += pow( w(istep,Sim_ord)*dE*RE_exp_syst[istep]*sm_coeff,2);
    }

    //add sum of absolute systematic error around resonances
    for(int ires=0;ires< (signed)Res_MIN.size();ires++) RE_sm_exp_syst[ierg] += pow( DE_corr[ires],2);

       
    RE_sm_exp_err[ierg] = RE_sm_exp_stat[ierg] + RE_sm_exp_syst[ierg];

    //take sqrt of stat^2 syst^2 and tot^2
    RE_sm_exp_err[ierg] = sqrt( RE_sm_exp_err[ierg]);
    RE_sm_exp_stat[ierg] = sqrt( RE_sm_exp_stat[ierg]);
    RE_sm_exp_syst[ierg] = sqrt( RE_sm_exp_syst[ierg]);

  }

     
     
  cout<<"Smeared experimental R(e+e- -> hadrons) with sigma: "<<sigma<<" computed!"<<endl;
  //Print To File
  Print_To_File({}, {Ergs_GeV_list ,RE_sm_exp_ave, RE_sm_exp_stat, RE_sm_exp_syst, RE_sm_exp_err}, "../dispersive_data/smeared_data/smeared_R_sigma_"+to_string_with_precision(sigma,3)+".dat", "", "#E[GeV]  R[E]    stat. syst.  total ");
}
   






void Spectral_tests() {


  int Np = (int)( max_energy - Emin)/(step_size*sm_size);
 
  //create output directories
  boost::filesystem::create_directory("../dispersive_data");
  boost::filesystem::create_directory("../dispersive_data/smeared_data");
  boost::filesystem::create_directory("../data/R_ratio");
  boost::filesystem::create_directory("../data/R_ratio/light");
  boost::filesystem::create_directory("../data/R_ratio/light_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/strange");
  boost::filesystem::create_directory("../data/R_ratio/strange_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/charm");
  boost::filesystem::create_directory("../data/R_ratio/charm_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/corr");
  boost::filesystem::create_directory("../data/R_ratio/corr/light");
  boost::filesystem::create_directory("../data/R_ratio/corr/light_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/corr/strange");
  boost::filesystem::create_directory("../data/R_ratio/corr/strange_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/corr/charm");
  boost::filesystem::create_directory("../data/R_ratio/corr/charm_w_disco");
  boost::filesystem::create_directory("../data/R_ratio/total");

  //LOAD GM2 DATA
  GaussianMersenne GM(981832);

  //light
  data_t  V_light_1, V_light_2, V_light_3;
  data_t  V_light_OS_1, V_light_OS_2, V_light_OS_3;
  data_t disco_light, disco_strange, disco_charm, disco_light_strange, disco_light_charm, disco_strange_charm;

  //strange and charm
  
  //L
  data_t  V_strange_1_L, V_strange_2_L, V_strange_3_L, V_strange_OS_1_L, V_strange_OS_2_L, V_strange_OS_3_L;
  data_t  V_charm_1_L, V_charm_2_L, V_charm_3_L, V_charm_OS_1_L, V_charm_OS_2_L, V_charm_OS_3_L;
  data_t  pt2_etaC_L, pt2_etaC_OS_L, pt2_etaS_L, pt2_etaS_OS_L;
  //M
  data_t  V_strange_1_M, V_strange_2_M, V_strange_3_M, V_strange_OS_1_M, V_strange_OS_2_M, V_strange_OS_3_M;
  data_t  V_charm_1_M, V_charm_2_M, V_charm_3_M, V_charm_OS_1_M, V_charm_OS_2_M, V_charm_OS_3_M;
  data_t  pt2_etaC_M, pt2_etaC_OS_M, pt2_etaS_M, pt2_etaS_OS_M;
  //H
  data_t  V_charm_1_H, V_charm_2_H, V_charm_3_H, V_charm_OS_1_H, V_charm_OS_2_H, V_charm_OS_3_H;
  data_t  pt2_etaC_H, pt2_etaC_OS_H;



  

 

  //Ergs_list
  Vfloat Ergs_GeV_list;
  int Nps = (int)((max_energy - Emin)/(sm_size*step_size));
  for(int ip=0;ip<Nps;ip++) Ergs_GeV_list.push_back(  Emin + ip*((max_energy- Emin)/((double)Nps) ));

  Vfloat Ergs_GeV_list_exp;
  double Dx_e= 0.01; //10 MeV
  int Nps_exp= (int)((5.0 - Emin)/Dx_e); //up to 5 GeV
  for(int ip=0;ip<Nps_exp;ip++) Ergs_GeV_list_exp.push_back(  Emin + ip*((5.0- Emin)/((double)Nps_exp)));
  for( auto &sigma: {0.03, 0.05, 0.20, 0.50, 1.0})
    if(Compute_experimental_smeared_R_ratio) Get_exp_smeared_R_ratio(Ergs_GeV_list_exp, sigma);
  


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


  // lambda function to be used as a smearing func.

  const auto f = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) -> PrecFloat {
		   
		   if(SM_TYPE == "WIN")  return (1.0/(1 + exp(-2*(E-(m-s))/s)) - 1.0/(1 + exp(-2*(E-(m+s))/s)))/(2*s*sqr(E));
		   
		   else if (SM_TYPE == "GAUSSIAN") return Get_exact_gauss(E, m, s, E0)/sqr(E);
		   
		   else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
		   
		   return PrecFloat(0);
		   
		 };

  





  //light 
  //#################################END CUSTOM SORTING#################
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);

  V_light_OS_1.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs);
  V_light_OS_2.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_light_OS_3.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);

  //disco_light
  disco_light.Read("../gm2_data_disc_Simone/loops/data/light_light_D", "disco", "", Sort_light_confs);
  //disco_strange
  disco_strange.Read("../gm2_data_disc_Simone/loops/data/strange_strange", "disco", "", Sort_light_confs);
  //disco_charm
  disco_charm.Read("../gm2_data_disc_Simone/loops/data/charm_charm", "disco", "", Sort_light_confs);
  //disco light-strange
  disco_light_strange.Read("../gm2_data_disc_Simone/loops/data/light_D_strange", "disco", "", Sort_light_confs);
  //disco light-charm
  disco_light_charm.Read("../gm2_data_disc_Simone/loops/data/light_D_charm", "disco", "", Sort_light_confs);
  //disco strange-charm
  disco_strange_charm.Read("../gm2_data_disc_Simone/loops/data/strange_charm", "disco", "", Sort_light_confs);
  //#################################END READING LIGHT########################

  //strange
  //L
  V_strange_1_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_strange_2_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  pt2_etaS_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  V_strange_OS_1_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  pt2_etaS_OS_L.Read("../gm2_data/strange_Nhits64_spectral/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  //M
  V_strange_1_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); 
  V_strange_2_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  pt2_etaS_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  V_strange_OS_1_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  pt2_etaS_OS_M.Read("../gm2_data/strange_Nhits64_spectral/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  //#################################END READING STRANGE########################



  //charm
  //L
  V_charm_1_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_L.Read( "../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_L.Read("../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_L.Read( "../gm2_data/charm_Nhits20_spectral/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //M
  V_charm_1_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs); 
  V_charm_2_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_M.Read( "../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_M.Read("../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_M.Read( "../gm2_data/charm_Nhits20_spectral/medium", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //H
  V_charm_1_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_1", "V1V1",  Sort_light_confs);
  V_charm_2_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_1", "V2V2",  Sort_light_confs);
  V_charm_3_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_1", "V3V3",  Sort_light_confs);
  pt2_etaC_H.Read( "../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs); 
  V_charm_OS_1_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_2", "V1V1",  Sort_light_confs); 
  V_charm_OS_2_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_2", "V2V2",  Sort_light_confs);
  V_charm_OS_3_H.Read("../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_2", "V3V3",  Sort_light_confs);
  pt2_etaC_OS_H.Read( "../gm2_data/charm_Nhits20_spectral/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 



  //###########################################################################################


  //Ensemble Tags

  //light Tags
  vector<string> Light_Ens_Tags;
  vector<string> Strange_Ens_Tags;
  vector<string> Charm_Ens_Tags;
  vector<string> Disco_Ens_Tags;


  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
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


  //############################################################################################
  //generate fake jack_distr for M_etas and M_phi phys
  distr_t Mphi_phys_distr(UseJack), Metas_phys_distr(UseJack), MJpsi_phys_distr(UseJack), Metac_phys_distr(UseJack);
  for(int ijack=0; ijack<Njacks;ijack++) {
    Mphi_phys_distr.distr.push_back( m_phi + GM()*m_phi_err*(1.0/sqrt(Njacks -1.0)));
    Metas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err*(1.0/sqrt(Njacks -1.0)));
    MJpsi_phys_distr.distr.push_back( m_Jpsi + GM()*m_Jpsi_err*(1.0/sqrt(Njacks -1.0)));
    Metac_phys_distr.distr.push_back( m_etac + GM()*m_etac_err*(1.0/sqrt(Njacks -1.0)));
  }



 

  int Nens_light =  V_light_1.size;
  int Nens_disco_light = disco_light.size;
  int Nens_strange = V_strange_1_L.size;
  int Nens_charm = V_charm_1_L.size;


  //vector of distr_t_list to store the R(E) results for all ensemble (different vectors for different contributions)

  //light
  vector<distr_t_list> RE_light_TANT_OS(Nens_light);
  vector<distr_t_list> RE_light_TANT_tm(Nens_light);
  vector<distr_t_list> RE_light_SANF_OS(Nens_light);
  vector<distr_t_list> RE_light_SANF_tm(Nens_light);

  //strange
  vector<distr_t_list> RE_strange_TANT_OS(Nens_strange);
  vector<distr_t_list> RE_strange_TANT_tm(Nens_strange);
  vector<distr_t_list> RE_strange_SANF_OS(Nens_strange);
  vector<distr_t_list> RE_strange_SANF_tm(Nens_strange);

  //charm
  vector<distr_t_list> RE_charm_TANT_OS(Nens_charm);
  vector<distr_t_list> RE_charm_TANT_tm(Nens_charm);
  vector<distr_t_list> RE_charm_SANF_OS(Nens_charm);
  vector<distr_t_list> RE_charm_SANF_tm(Nens_charm);

  //disconnected
  vector<distr_t_list> RE_disco_TANT(Nens_charm);
  vector<distr_t_list> RE_disco_SANF(Nens_charm);

  //total
  vector<distr_t_list> RE_tot_TANT_OS(Nens_light);
  vector<distr_t_list> RE_tot_TANT_tm(Nens_light);
  vector<distr_t_list> RE_tot_SANF_OS(Nens_light);
  vector<distr_t_list> RE_tot_SANF_tm(Nens_light);

  //total optimized
  vector<distr_t_list> RE_tot_opt_OS(Nens_light);
  vector<distr_t_list> RE_tot_opt_tm(Nens_light);

  Light_Ens_Tags.resize(Nens_light);
  Strange_Ens_Tags.resize(Nens_strange);
  Charm_Ens_Tags.resize(Nens_charm);










  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE CHARM CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

  //get charm correlators only on B,C,D physical point ensemble

    
  #pragma omp parallel for
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) {

      
    Charm_Ens_Tags[i_ens] = V_charm_1_L.Tag[i_ens];


    string FLAV= "charm";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_charm_1_L.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_charm_1_L.Tag[i_ens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(V_charm_1_L.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv = ZV_A; Za= ZA_A;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za= ZA_B;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za= ZA_C;}
    else if(V_charm_1_L.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za= ZA_D;}
    else crash("lattice spacing distribution for Ens: "+V_charm_1_L.Tag[i_ens]+" not found");



    //set Tmin and Tmax for the eta_C
    int Tmin_etaC, Tmax_etaC;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_etaC=40; Tmax_etaC=65;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_etaC=35; Tmax_etaC=50;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_etaC=32;Tmax_etaC=45;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_etaC=41; Tmax_etaC= 57;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_etaC=28; Tmax_etaC=45;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_etaC=20; Tmax_etaC=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_etaC=20; Tmax_etaC=23;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_etaC=25;Tmax_etaC=31;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens]== "cD211a.054.96") {Tmin_etaC=50; Tmax_etaC=71;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    
    else crash("Ensemble tag not valid");

    Corr.Tmin = Tmin_etaC; Corr.Tmax= Tmax_etaC;

    int Tmin_VV, Tmax_VV;
 
    //set time intervals for vector obs
    if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=43; Tmax_VV=58;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=31; Tmax_VV=43;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=32;Tmax_VV=42;}
      else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=28; Tmax_VV= 40;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=27; Tmax_VV=42;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=19; Tmax_VV=24;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=20; Tmax_VV=24;}
      else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=24;Tmax_VV=31;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_charm_1_L.Tag[i_ens]== "cD211a.054.96") {Tmin_VV=49; Tmax_VV=65;}
      else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
    }
    
    else crash("Ensemble tag not valid");

    

      
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_charm_1_L.Tag[i_ens]);
 
    //L
    distr_t_list  V_charm_1_L_distr, V_charm_2_L_distr, V_charm_3_L_distr, V_charm_L_distr;
    distr_t_list  V_charm_OS_1_L_distr, V_charm_OS_2_L_distr, V_charm_OS_3_L_distr, V_charm_OS_L_distr;
    distr_t  Metac_L_distr, Metac_OS_L_distr, MV_L_distr, MV_OS_L_distr;
    //M
    distr_t_list  V_charm_1_M_distr, V_charm_2_M_distr, V_charm_3_M_distr, V_charm_M_distr;
    distr_t_list  V_charm_OS_1_M_distr, V_charm_OS_2_M_distr, V_charm_OS_3_M_distr, V_charm_OS_M_distr;
    distr_t Metac_M_distr, Metac_OS_M_distr, MV_M_distr, MV_OS_M_distr;
    //H
    distr_t_list  V_charm_1_H_distr, V_charm_2_H_distr, V_charm_3_H_distr, V_charm_H_distr;
    distr_t_list  V_charm_OS_1_H_distr, V_charm_OS_2_H_distr, V_charm_OS_3_H_distr, V_charm_OS_H_distr;
    distr_t Metac_H_distr, Metac_OS_H_distr, MV_H_distr, MV_OS_H_distr;


    //disco
    distr_t_list disco_distr, disco_light_charm_distr;
  
  

    //vector light sector
    //L
    V_charm_1_L_distr = Corr.corr_t(V_charm_1_L.col(0)[i_ens], "");
    V_charm_2_L_distr = Corr.corr_t(V_charm_2_L.col(0)[i_ens], "");
    V_charm_3_L_distr = Corr.corr_t(V_charm_3_L.col(0)[i_ens], "");
    V_charm_OS_1_L_distr = Corr.corr_t(V_charm_OS_1_L.col(0)[i_ens], "");
    V_charm_OS_2_L_distr = Corr.corr_t(V_charm_OS_2_L.col(0)[i_ens], "");
    V_charm_OS_3_L_distr = Corr.corr_t(V_charm_OS_3_L.col(0)[i_ens], "");
    //M
    V_charm_1_M_distr = Corr.corr_t(V_charm_1_M.col(0)[i_ens], "");
    V_charm_2_M_distr = Corr.corr_t(V_charm_2_M.col(0)[i_ens], "");
    V_charm_3_M_distr = Corr.corr_t(V_charm_3_M.col(0)[i_ens], "");
    V_charm_OS_1_M_distr = Corr.corr_t(V_charm_OS_1_M.col(0)[i_ens], "");
    V_charm_OS_2_M_distr = Corr.corr_t(V_charm_OS_2_M.col(0)[i_ens], "");
    V_charm_OS_3_M_distr = Corr.corr_t(V_charm_OS_3_M.col(0)[i_ens], "");

    //L
    V_charm_1_H_distr = Corr.corr_t(V_charm_1_H.col(0)[i_ens], "");
    V_charm_2_H_distr = Corr.corr_t(V_charm_2_H.col(0)[i_ens], "");
    V_charm_3_H_distr = Corr.corr_t(V_charm_3_H.col(0)[i_ens], "");
    V_charm_OS_1_H_distr = Corr.corr_t(V_charm_OS_1_H.col(0)[i_ens], "");
    V_charm_OS_2_H_distr = Corr.corr_t(V_charm_OS_2_H.col(0)[i_ens], "");
    V_charm_OS_3_H_distr = Corr.corr_t(V_charm_OS_3_H.col(0)[i_ens], "");


    Corr.Tmin = Tmin_etaC; Corr.Tmax = Tmax_etaC;
    //P5P5
    //L
    Metac_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_L.col(0)[i_ens], ""));
    Metac_OS_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_L.col(0)[i_ens], ""));
    //M
    Metac_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_M.col(0)[i_ens], ""));
    Metac_OS_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_M.col(0)[i_ens], ""));
    //H
    Metac_H_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_H.col(0)[i_ens], ""));
    Metac_OS_H_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaC_OS_H.col(0)[i_ens], ""));



     
      
    
      
    //sum over the Lorentz indices of the e.m. current
    //L
    V_charm_L_distr= (pow(qc,2)/3.0)*(V_charm_1_L_distr+ V_charm_2_L_distr + V_charm_3_L_distr);
    V_charm_OS_L_distr =  (pow(qc,2)/3.0)*(V_charm_OS_1_L_distr+ V_charm_OS_2_L_distr + V_charm_OS_3_L_distr);
    //M
    V_charm_M_distr= (pow(qc,2)/3.0)*(V_charm_1_M_distr+ V_charm_2_M_distr + V_charm_3_M_distr);
    V_charm_OS_M_distr =  (pow(qc,2)/3.0)*(V_charm_OS_1_M_distr+ V_charm_OS_2_M_distr + V_charm_OS_3_M_distr);
    //M
    V_charm_H_distr= (pow(qc,2)/3.0)*(V_charm_1_H_distr+ V_charm_2_H_distr + V_charm_3_H_distr);
    V_charm_OS_H_distr =  (pow(qc,2)/3.0)*(V_charm_OS_1_H_distr+ V_charm_OS_2_H_distr + V_charm_OS_3_H_distr);


    


    string charm_tag = "charm";

    


    //get MV from vector correlator
    Corr.Tmin= Tmin_VV; Corr.Tmax = Tmax_VV;
    //L
    MV_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_L_distr, ""));
    MV_OS_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_L_distr, ""));
    //M
    MV_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_M_distr, ""));
    MV_OS_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_M_distr, ""));
    //H
    MV_H_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_H_distr, ""));
    MV_OS_H_distr = Corr.Fit_distr( Corr.effective_mass_t(V_charm_OS_H_distr, ""));


    vector<distr_t> Metac_vec, MJpsi_vec;
    if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
      Metac_vec = { Metac_L_distr/a_distr, Metac_M_distr/a_distr, Metac_H_distr/a_distr};
      MJpsi_vec = { MV_L_distr/a_distr, MV_M_distr/a_distr, MV_H_distr/a_distr};
    }
    else {
      Metac_vec = { Metac_L_distr/a_distr, Metac_M_distr/a_distr};
      MJpsi_vec = { MV_L_distr/a_distr, MV_M_distr/a_distr};
    }
    vector<distr_t> X_fit;
    if(Extrapolation_charm_mode == "Jpsi") X_fit = MJpsi_vec;
    else if(Extrapolation_charm_mode == "etac") X_fit = Metac_vec;
    else crash("Extrapolation charm mode: "+Extrapolation_charm_mode+" is invalid");
    distr_t X_phys = (Extrapolation_charm_mode=="Jpsi")?MJpsi_phys_distr:Metac_phys_distr;


   


    //############################################################################################


    //perturbative substraction of lattice artifacts

    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_L,3)+"/OPPOR";
    string Pt_free_oppor_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_M,3)+"/OPPOR";
    string Pt_free_oppor_H= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_H,3)+"/OPPOR";
    Vfloat VV_free_oppor_L= Read_From_File(Pt_free_oppor_L, 1, 4);
    Vfloat VV_free_oppor_M= Read_From_File(Pt_free_oppor_M, 1, 4);
    Vfloat VV_free_oppor_H= Read_From_File(Pt_free_oppor_H, 1, 4);
    if(VV_free_oppor_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L w opposite r");
    if(VV_free_oppor_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M w opposite r");
    if(VV_free_oppor_H.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_H w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_L,3)+"/SAMER";
    string Pt_free_samer_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_M,3)+"/SAMER";
    string Pt_free_samer_H= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.mc_H,3)+"/SAMER";
    Vfloat VV_free_samer_L= Read_From_File(Pt_free_samer_L, 1, 4);
    Vfloat VV_free_samer_M= Read_From_File(Pt_free_samer_M, 1, 4);
    Vfloat VV_free_samer_H= Read_From_File(Pt_free_samer_H, 1, 4);
    if(VV_free_samer_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L  w same r");
    if(VV_free_samer_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M  w same r");
    if(VV_free_samer_H.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_H  w same r");
    //Insert electric charges
    for( auto & OP:VV_free_oppor_L) OP *= pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_L) SA *= pert_corr_charm_on_off*qc*qc;
    for( auto & OP:VV_free_oppor_M) OP *= pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_M) SA *= pert_corr_charm_on_off*qc*qc;
    for( auto & OP:VV_free_oppor_H) OP *= pert_corr_charm_on_off*qc*qc;
    for( auto & SA:VV_free_samer_H) SA *= pert_corr_charm_on_off*qc*qc;
  
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) {  if( t*a_distr.ave() < 1.0*fm_to_inv_Gev && t != 0) {
	free_corr_log_art[t] = -1.0*pert_corr_charm_on_off*(qc*qc)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));
      }
      if(t==0 || t*a_distr.ave() > add_pert_corr_charm_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0; VV_free_samer_H[t] =0; VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0; VV_free_oppor_H[t]=0;}
    }


    //multiply corr using Zv and Za
    //L
    V_charm_L_distr = (V_charm_L_distr)*Za*Za + VV_free_oppor_L ;
    V_charm_OS_L_distr = (V_charm_OS_L_distr)*Zv*Zv + VV_free_samer_L;
    //M
    V_charm_M_distr = (V_charm_M_distr)*Za*Za + VV_free_oppor_M ;
    V_charm_OS_M_distr = (V_charm_OS_M_distr)*Zv*Zv + VV_free_samer_M;
    //H
    V_charm_H_distr = (V_charm_H_distr)*Za*Za + VV_free_oppor_H ;
    V_charm_OS_H_distr = (V_charm_OS_H_distr)*Zv*Zv + VV_free_samer_H;


      


    //print correlator

    //L
    Print_To_File({}, {V_charm_L_distr.ave(), V_charm_L_distr.err(), V_charm_OS_L_distr.ave(), V_charm_OS_L_distr.err()}, "../data/R_ratio/corr/"+charm_tag+"/corr_L_"+V_charm_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //M
    Print_To_File({}, {V_charm_M_distr.ave(), V_charm_M_distr.err(), V_charm_OS_M_distr.ave(), V_charm_OS_M_distr.err()}, "../data/R_ratio/corr/"+charm_tag+"/corr_M_"+V_charm_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //H
    Print_To_File({}, {V_charm_H_distr.ave(), V_charm_H_distr.err(), V_charm_OS_H_distr.ave(), V_charm_OS_H_distr.err()}, "../data/R_ratio/corr/"+charm_tag+"/corr_H_"+V_charm_1_H.Tag[i_ens]+".dat", "", "# t  tm  OS");


    cout<<"Perturbative corrections added"<<endl;





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    //L
    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax_L=1;
    while(!Found_error_less_x_percent && tmax_L < Corr.Nt/2) {
   
      if( (V_charm_L_distr.distr_list[tmax_L]).err()/fabs( (V_charm_L_distr.distr_list[tmax_L]).ave()) <  0.01*x) tmax_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_L=1;
    while(!Found_error_less_x_percent && tmax_OS_L < Corr.Nt/2) {
   
      if( (V_charm_OS_L_distr.distr_list[tmax_OS_L]).err()/fabs( (V_charm_OS_L_distr.distr_list[tmax_OS_L]).ave()) <  0.01*x) tmax_OS_L++;
      else Found_error_less_x_percent=true;
    }

    //M
    Found_error_less_x_percent=false;
    double tmax_M=1;
    while(!Found_error_less_x_percent && tmax_M < Corr.Nt/2) {
   
      if( (V_charm_M_distr.distr_list[tmax_M]).err()/fabs( (V_charm_M_distr.distr_list[tmax_M]).ave()) <  0.01*x) tmax_M++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_M=1;
    while(!Found_error_less_x_percent && tmax_OS_M < Corr.Nt/2) {
   
      if( (V_charm_OS_M_distr.distr_list[tmax_OS_M]).err()/fabs( (V_charm_OS_M_distr.distr_list[tmax_OS_M]).ave()) <  0.01*x) tmax_OS_M++;
      else Found_error_less_x_percent=true;
    }

    //H
    Found_error_less_x_percent=false;
    double tmax_H=1;
    while(!Found_error_less_x_percent && tmax_H < Corr.Nt/2) {
   
      if( (V_charm_H_distr.distr_list[tmax_H]).err()/fabs( (V_charm_H_distr.distr_list[tmax_H]).ave()) <  0.01*x) tmax_H++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_H=1;
    while(!Found_error_less_x_percent && tmax_OS_H < Corr.Nt/2) {
   
      if( (V_charm_OS_H_distr.distr_list[tmax_OS_H]).err()/fabs( (V_charm_OS_H_distr.distr_list[tmax_OS_H]).ave()) <  0.01*x) tmax_OS_H++;
      else Found_error_less_x_percent=true;
    }
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double Emax= max_energy*a_distr.ave();  //[x GeV]*a
    double E0 =  Emin*a_distr.ave(); // two pion treshold ~ 280 MeV
    double Estart= Emin*a_distr.ave(); // start from 300 MeV
    double sigma = sm_size*a_distr.ave(); // smearing
    int Npoints = (int)((Emax - Estart)/(sigma*step_size));
    tmax_L= (T/2 -4);
    tmax_OS_L = T/2 -4;
    tmax_M= (T/2 -4);
    tmax_OS_M = T/2 -4;
    tmax_H = (T/2 -4);
    tmax_OS_H = (T/2 - 4);
    double mult_SANF = 0.00000002;
    double mult_TANT=0.9;
    //#########################################################################################


    cout<<"Reconstructing R(E)....."<<endl;
    cout<<"aEmax: "<<Emax<<" -> Emax: "<<Emax/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"asigma: "<<sigma<<" -> sigma: "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"Npoints: "<<Npoints<<endl;
    cout<<"Zv: "<<L_info.Zv_WI<< " +- "<<L_info.Zv_WI_err<<endl;
    cout<<"Za: "<<L_info.Za_WI<< " +- "<<L_info.Za_WI_err<<endl;
    cout<<"SM_TYPE: "<<SM_TYPE<<endl;
    

    distr_t_list Spectral_dens_L;
    distr_t_list Spectral_dens_OS_L;
    distr_t_list Spectral_dens_SANF_L;
    distr_t_list Spectral_dens_OS_SANF_L;
    distr_t_list Spectral_dens_M;
    distr_t_list Spectral_dens_OS_M;
    distr_t_list Spectral_dens_SANF_M;
    distr_t_list Spectral_dens_OS_SANF_M;
    distr_t_list Spectral_dens_H;
    distr_t_list Spectral_dens_OS_H;
    distr_t_list Spectral_dens_SANF_H;
    distr_t_list Spectral_dens_OS_SANF_H;
    distr_t_list Spectral_dens_Extr;
    distr_t_list Spectral_dens_OS_Extr;
    distr_t_list Spectral_dens_SANF_Extr;
    distr_t_list Spectral_dens_OS_SANF_Extr;
    Vfloat Ergs;
    Vfloat Ergs_GeV;
    Vfloat syst_tm_L(Npoints);
    Vfloat syst_OS_L(Npoints);
    Vfloat syst_tm_SANF_L(Npoints);
    Vfloat syst_OS_SANF_L(Npoints);
    Vfloat syst_tm_M(Npoints);
    Vfloat syst_OS_M(Npoints);
    Vfloat syst_tm_SANF_M(Npoints);
    Vfloat syst_OS_SANF_M(Npoints);
    Vfloat syst_tm_H(Npoints);
    Vfloat syst_OS_H(Npoints);
    Vfloat syst_tm_SANF_H(Npoints);
    Vfloat syst_OS_SANF_H(Npoints);


    //COMPUTE THE SMEARED R-RATIO

    for(int ip=0; ip<Npoints;ip++) {

      double mean = Estart + ip*((Emax-Estart)/((double)Npoints));
      Ergs.push_back(mean);
      double mean_GeV = mean/a_distr.ave();
      Ergs_GeV.push_back(mean_GeV);
      double sigma_resc=sigma;
      double lambda_Estar;
      double lambda_Estar_SANF;
 
      cout<<"Evaluating R(E*) at aE*: "<<mean<<" -> "<<"E*: "<<mean/a_distr.ave()<<" [GeV] ..."<<endl;
      cout<<"Rescaled asigma(E*): "<<sigma_resc<<" -> sigma(E*): "<<sigma_resc/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"[tmin, tmax] (tm L): ["<<1<<", "<<tmax_L<<"]"<<endl;
      cout<<"[tmin, tmax] (OS L): ["<<1<<", "<<tmax_OS_L<<"]"<<endl;
      cout<<"[tmin, tmax] (tm M): ["<<1<<", "<<tmax_M<<"]"<<endl;
      cout<<"[tmin, tmax] (OS M): ["<<1<<", "<<tmax_OS_M<<"]"<<endl;
      cout<<"[tmin, tmax] (tm H): ["<<1<<", "<<tmax_H<<"]"<<endl;
      cout<<"[tmin, tmax] (OS H): ["<<1<<", "<<tmax_OS_H<<"]"<<endl;

      //L
      Spectral_dens_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_L, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_L_distr, syst_tm_L[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+charm_tag+"_L_"+V_charm_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm L (TANTALO): "<<Spectral_dens_L.ave(ip)<<" +- ("<<Spectral_dens_L.err(ip)<<")_stat ("<<syst_tm_L[ip]<<")_syst "<<endl;
      Spectral_dens_SANF_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_L, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_L_distr, syst_tm_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+charm_tag+"_L_"+V_charm_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm L (SANFILIPPO): "<<Spectral_dens_SANF_L.ave(ip)<<" +- ("<<Spectral_dens_SANF_L.err(ip)<<")_stat ("<<syst_tm_SANF_L[ip]<<")_syst "<<endl;
      Spectral_dens_OS_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_L, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_L_distr, syst_OS_L[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+charm_tag+"_L_"+V_charm_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS L (TANTALO): "<<Spectral_dens_OS_L.ave(ip)<<" +- ("<<Spectral_dens_OS_L.err(ip)<<")_stat ("<<syst_OS_L[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_L, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_L_distr, syst_OS_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+charm_tag+"_L_"+V_charm_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS L (SANFILIPPO): "<<Spectral_dens_OS_SANF_L.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF_L.err(ip)<<")_stat ("<<syst_OS_SANF_L[ip]<<")_syst "<<endl;
      //M
      Spectral_dens_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_M, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_M_distr, syst_tm_M[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+charm_tag+"_M_"+V_charm_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm M (TANTALO): "<<Spectral_dens_M.ave(ip)<<" +- ("<<Spectral_dens_M.err(ip)<<")_stat ("<<syst_tm_M[ip]<<")_syst "<<endl;
      Spectral_dens_SANF_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_M, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_M_distr, syst_tm_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+charm_tag+"_M_"+V_charm_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm M (SANFILIPPO): "<<Spectral_dens_SANF_M.ave(ip)<<" +- ("<<Spectral_dens_SANF_M.err(ip)<<")_stat ("<<syst_tm_SANF_M[ip]<<")_syst "<<endl;
      Spectral_dens_OS_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_M, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_M_distr, syst_OS_M[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+charm_tag+"_M_"+V_charm_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS M (TANTALO): "<<Spectral_dens_OS_M.ave(ip)<<" +- ("<<Spectral_dens_OS_M.err(ip)<<")_stat ("<<syst_OS_M[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_M, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_M_distr, syst_OS_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+charm_tag+"_M_"+V_charm_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS M (SANFILIPPO): "<<Spectral_dens_OS_SANF_M.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF_M.err(ip)<<")_stat ("<<syst_OS_SANF_M[ip]<<")_syst "<<endl;
      //H
      Spectral_dens_H.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_H, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_H_distr, syst_tm_H[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+charm_tag+"_H_"+V_charm_1_H.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm H (TANTALO): "<<Spectral_dens_H.ave(ip)<<" +- ("<<Spectral_dens_H.err(ip)<<")_stat ("<<syst_tm_H[ip]<<")_syst "<<endl;
      Spectral_dens_SANF_H.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_H, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_H_distr, syst_tm_SANF_H[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+charm_tag+"_H_"+V_charm_1_H.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm H (SANFILIPPO): "<<Spectral_dens_SANF_H.ave(ip)<<" +- ("<<Spectral_dens_SANF_H.err(ip)<<")_stat ("<<syst_tm_SANF_H[ip]<<")_syst "<<endl;
      Spectral_dens_OS_H.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_H, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_H_distr, syst_OS_H[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+charm_tag+"_H_"+V_charm_1_H.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS H (TANTALO): "<<Spectral_dens_OS_H.ave(ip)<<" +- ("<<Spectral_dens_OS_H.err(ip)<<")_stat ("<<syst_OS_H[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF_H.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_H, prec_charm, SM_TYPE+"_ov_E2",f, rho_R*V_charm_OS_H_distr, syst_OS_SANF_H[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+charm_tag+"_H_"+V_charm_1_H.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS H (SANFILIPPO): "<<Spectral_dens_OS_SANF_H.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF_H.err(ip)<<")_stat ("<<syst_OS_SANF_H[ip]<<")_syst "<<endl;
      


      //Extrapolate to the physical charm point
      vector<distr_t> Y_fit_tm, Y_fit_OS, Y_fit_SANF_tm, Y_fit_SANF_OS;
      if(V_charm_1_L.Tag[i_ens].substr(1,1) != "D") {
	Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip], Spectral_dens_H.distr_list[ip]};
	Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip], Spectral_dens_OS_H.distr_list[ip]};
	Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip], Spectral_dens_SANF_H.distr_list[ip]};
	Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip], Spectral_dens_OS_SANF_H.distr_list[ip]};
      }
      else {
	Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip]};
	Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip]};
	Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip]};
	Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip]};
      }

      Spectral_dens_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_tm, X_fit, X_phys ,  "../data/R_ratio/"+charm_tag  , "R_extr_tm_TANT_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_OS_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_OS, X_fit, X_phys ,  "../data/R_ratio/"+charm_tag  , "R_extr_OS_TANT_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_SANF_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_SANF_tm, X_fit, X_phys ,  "../data/R_ratio/"+charm_tag  , "R_extr_tm_SANF_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_OS_SANF_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_SANF_OS, X_fit, X_phys ,  "../data/R_ratio/"+charm_tag  , "R_extr_OS_SANF_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
	
	

     
      cout<<"...done"<<endl;
      //############################################################################################



    }

    RE_charm_TANT_tm[i_ens] = Spectral_dens_Extr;
    RE_charm_TANT_OS[i_ens] = Spectral_dens_OS_Extr;
    RE_charm_SANF_tm[i_ens] = Spectral_dens_SANF_Extr;
    RE_charm_SANF_OS[i_ens] = Spectral_dens_OS_SANF_Extr;
   

      
    cout<<"printing output charm..."<<endl;

    //print to file
    //L
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_L.ave(), Spectral_dens_L.err(),  Spectral_dens_SANF_L.ave(), Spectral_dens_SANF_L.err(), Spectral_dens_OS_L.ave(), Spectral_dens_OS_L.err(),  Spectral_dens_OS_SANF_L.ave(), Spectral_dens_OS_SANF_L.err()}, "../data/R_ratio/"+charm_tag+"/L_"+SM_TYPE+"_"+V_charm_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat. ) [SANFILIPPO]   R(E)_OS (<> stat. ) [TANTALO]  R(E)_OS (<> stat. ) [SANFILIPPO]");
    //M
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_M.ave(), Spectral_dens_M.err(), Spectral_dens_SANF_M.ave(), Spectral_dens_SANF_M.err(), Spectral_dens_OS_M.ave(), Spectral_dens_OS_M.err(),  Spectral_dens_OS_SANF_M.ave(), Spectral_dens_OS_SANF_M.err()}, "../data/R_ratio/"+charm_tag+"/M_"+SM_TYPE+"_"+V_charm_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat ) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");
    //H
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_H.ave(), Spectral_dens_H.err(), Spectral_dens_SANF_H.ave(), Spectral_dens_SANF_H.err(), Spectral_dens_OS_H.ave(), Spectral_dens_OS_H.err(), Spectral_dens_OS_SANF_H.ave(), Spectral_dens_OS_SANF_H.err()}, "../data/R_ratio/"+charm_tag+"/H_"+SM_TYPE+"_"+V_charm_1_H.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");
    //Extr
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_Extr.ave(), Spectral_dens_Extr.err(), Spectral_dens_SANF_Extr.ave(), Spectral_dens_SANF_Extr.err(), Spectral_dens_OS_Extr.ave(), Spectral_dens_OS_Extr.err(), Spectral_dens_OS_SANF_Extr.ave(), Spectral_dens_OS_SANF_Extr.err() }, "../data/R_ratio/"+charm_tag+"/Extr_"+SM_TYPE+"_"+V_charm_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");

    cout<<"done!"<<endl;


  }
  
  #pragma omp barrier

  cout<<"Calculation of charm contribution COMPLETED"<<endl;


  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE LIGHT CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

  
  #pragma omp parallel for
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {


    
    Light_Ens_Tags[i_ens] = V_light_1.Tag[i_ens];

    string FLAV= "light";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_light_1.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_light_1.Tag[i_ens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(V_light_1.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv = ZV_A; Za = ZA_A;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+V_light_1.Tag[i_ens]+" not found");
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
 
   
    distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr, disco_distr;
    distr_t_list  V_light_OS_1_distr, V_light_OS_2_distr, V_light_OS_3_distr, V_light_OS_distr;  
    distr_t_list  V_light_distr;
    
  

    //vector light sector
    V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "");
    V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "");
    V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "");
    V_light_OS_1_distr = Corr.corr_t(V_light_OS_1.col(0)[i_ens], "");
    V_light_OS_2_distr = Corr.corr_t(V_light_OS_2.col(0)[i_ens], "");
    V_light_OS_3_distr = Corr.corr_t(V_light_OS_3.col(0)[i_ens], "");



    //sum over the Lorentz indices of the e.m. current
    V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
    V_light_OS_distr =  ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_OS_1_distr+ V_light_OS_2_distr + V_light_OS_3_distr);
    

      
    string light_tag = "light";


    //############################################################################################


    //leading order perturbative substraction of a^2/t^2 lattice artifacts
  
    //free corr LO artifacts
    Vfloat free_corr_log_art(Corr.Nt);
    for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave()/fm_to_inv_Gev < 1.0 && t != 0) {   free_corr_log_art[t] = -1.0*(qu*qu + qd*qd)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;}

    
    //multiply corr using Zv and Za
    V_light_distr = (V_light_distr)*Za*Za;
    V_light_OS_distr = (V_light_OS_distr)*Zv*Zv;

    if(pert_corr_light_on_off) {

      V_light_distr = V_light_distr + free_corr_log_art;
      V_light_OS_distr = V_light_OS_distr + free_corr_log_art;

    }

 


    //print correlator

    Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), V_light_OS_distr.ave(), V_light_OS_distr.err()}, "../data/R_ratio/corr/"+light_tag+"/corr_"+V_light_1.Tag[i_ens]+".dat", "", "# t  tm  OS");


    cout<<"Perturbative corrections added"<<endl;





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (V_light_distr.distr_list[tmax]).err()/fabs( (V_light_distr.distr_list[tmax]).ave()) <  0.01*x) tmax++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS=1;
    while(!Found_error_less_x_percent && tmax_OS < Corr.Nt/2) {
   
      if( (V_light_OS_distr.distr_list[tmax_OS]).err()/fabs( (V_light_OS_distr.distr_list[tmax_OS]).ave()) <  0.01*x) tmax_OS++;
      else Found_error_less_x_percent=true;
    }
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double Emax= max_energy*a_distr.ave();  //[x GeV]*a
    double E0 = Emin*a_distr.ave(); // two pion treshold ~ 280 MeV
    double Estart= Emin*a_distr.ave(); // start from 300 MeV
    double sigma = sm_size*a_distr.ave(); // smearing
    
    int Npoints = (int)((Emax - Estart)/(sigma*step_size));
    tmax= (T/2 -4);
    tmax_OS = T/2 -4;
    double mult_TANT= 1.0;
    double mult_SANF = 0.00001;
    //#########################################################################################


    cout<<"Reconstructing R(E)....."<<endl;
    cout<<"aEmax: "<<Emax<<" -> Emax: "<<Emax/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"asigma: "<<sigma<<" -> sigma: "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"Npoints: "<<Npoints<<endl;
    cout<<"Corr(0) tm: "<<(rho_R*V_light_distr).ave(0)<<endl;
    cout<<"Corr(0) OS: "<<(rho_R*V_light_OS_distr).ave(0)<<endl;
    cout<<"Zv: "<<L_info.Zv_WI<< " +- "<<L_info.Zv_WI_err<<endl;
    cout<<"Za: "<<L_info.Za_WI<< " +- "<<L_info.Za_WI_err<<endl;
    cout<<"SM_TYPE: "<<SM_TYPE<<endl;
    

  
 

    distr_t_list Spectral_dens;
    distr_t_list Spectral_dens_OS;
    distr_t_list Spectral_dens_SANF;
    distr_t_list Spectral_dens_OS_SANF;
    Vfloat Ergs;
    Vfloat Ergs_GeV;
    Vfloat lambdas_tm;
    Vfloat lambdas_OS;
    Vfloat syst_tm(Npoints);
    Vfloat syst_OS(Npoints);
    Vfloat syst_tm_SANF(Npoints);
    Vfloat syst_OS_SANF(Npoints);


    //COMPUTE THE SMEARED R-RATIO

    for(int ip=0; ip<Npoints;ip++) {

      double mean = Estart + ip*((Emax-Estart)/((double)Npoints));
      Ergs.push_back(mean);
      Ergs_GeV.push_back(mean/a_distr.ave());
      double sigma_resc=sigma;
      double lambda_Estar;
      double lambda_Estar_SANF;
   
      cout<<"Evaluating R(E*) on Ens "<<V_light_1.Tag[i_ens]<<" at aE*: "<<mean<<" -> "<<"E*: "<<mean/a_distr.ave()<<" [GeV] ..."<<endl;
      cout<<"Rescaled asigma(E*): "<<sigma_resc<<" -> sigma(E*): "<<sigma_resc/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"[tmin, tmax] (tm): ["<<1<<", "<<tmax<<"]"<<endl;
      cout<<"[tmin, tmax] (OS): ["<<1<<", "<<tmax_OS<<"]"<<endl;
      Spectral_dens.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_distr, syst_tm[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+light_tag+"_"+V_light_1.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm (TANTALO): "<<Spectral_dens.ave(ip)<<" +- ("<<Spectral_dens.err(ip)<<")_stat ("<<syst_tm[ip]<<")_syst "<<endl;
      Spectral_dens_SANF.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_distr, syst_tm_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+light_tag+"_"+V_light_1.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm (SANFILIPPO): "<<Spectral_dens_SANF.ave(ip)<<" +- ("<<Spectral_dens_SANF.err(ip)<<")_stat ("<<syst_tm_SANF[ip]<<")_syst "<<endl;
      lambdas_tm.push_back(lambda_Estar);
    
      Spectral_dens_OS.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_OS_distr, syst_OS[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+light_tag+"_"+V_light_1.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS (TANTALO): "<<Spectral_dens_OS.ave(ip)<<" +- ("<<Spectral_dens_OS.err(ip)<<")_stat ("<<syst_OS[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_OS_distr, syst_OS_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+light_tag+"_"+V_light_1.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS (SANFILIPPO): "<<Spectral_dens_OS_SANF.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF.err(ip)<<")_stat ("<<syst_OS_SANF[ip]<<")_syst "<<endl;
      lambdas_OS.push_back(lambda_Estar);
     
      cout<<"...done"<<endl;
      //############################################################################################



    }

    


    RE_light_TANT_tm[i_ens] = Spectral_dens;
    RE_light_TANT_OS[i_ens] = Spectral_dens_OS;
    RE_light_SANF_tm[i_ens] = Spectral_dens_SANF;
    RE_light_SANF_OS[i_ens] = Spectral_dens_OS_SANF;

    
    cout<<"printing output light..."<<endl;

    //light
    //print to file
    Print_To_File({}, {Ergs, Ergs_GeV, Spectral_dens.ave(), Spectral_dens.err(), Spectral_dens_SANF.ave(), Spectral_dens_SANF.err(), Spectral_dens_OS.ave(), Spectral_dens_OS.err(), Spectral_dens_OS_SANF.ave(), Spectral_dens_OS_SANF.err()}, "../data/R_ratio/"+light_tag+"/"+SM_TYPE+"_"+V_light_1.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");

    cout<<"done!"<<endl;


  }

  #pragma omp barrier

  cout<<"Calculation of light contribution COMPLETED"<<endl;

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE STRANGE CONTRIBUTION ##################################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

 
  #pragma omp parallel for
  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {


    
    
    Strange_Ens_Tags[i_ens] = V_strange_1_L.Tag[i_ens];

    string FLAV= "strange";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_strange_1_L.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_strange_1_L.Tag[i_ens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(V_strange_1_L.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv= ZV_A; Za = ZA_A;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv= ZV_B; Za = ZA_B;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv= ZV_C; Za = ZA_C;}
    else if(V_strange_1_L.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv= ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+V_strange_1_L.Tag[i_ens]+" not found");



    //set time intervals for pseudoscalar obs
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=40; Corr.Tmax=70;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=31; Corr.Tmax=58;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=23;Corr.Tmax=44;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=36; Corr.Tmax= 57;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Corr.Tmin=40; Corr.Tmax= 80;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=19; Corr.Tmax=33;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=18; Corr.Tmax=23;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=16; Corr.Tmax=22;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Corr.Tmin=55; Corr.Tmax=88;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    int Tmin_P5P5 = Corr.Tmin;
    int Tmax_P5P5 = Corr.Tmax;

    int Tmin_VV, Tmax_VV;
    int Tmin_VV_OS, Tmax_VV_OS;

    //set time intervals for vector obs
    if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
      if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=26; Tmax_VV=35; Tmin_VV_OS=24; Tmax_VV_OS= 33;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
      if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=18; Tmax_VV=28; Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=19;Tmax_VV=25;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=25; Tmax_VV= 34; Tmin_VV_OS= 26; Tmax_VV_OS=34;} //25:32, 25:32
      else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_VV=23; Tmax_VV= 31; Tmin_VV_OS=25; Tmax_VV_OS=34;} //25:32 , 25:32
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
      if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=18; Tmax_VV=28;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=13; Tmax_VV=20;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=15; Tmax_VV=20;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=17;Tmax_VV=23;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
      if(V_strange_1_L.Tag[i_ens]=="cD211a.054.96") { Tmin_VV=32; Tmax_VV=40;  Tmin_VV_OS=36; Tmax_VV_OS=45;}
      else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
    }
    else crash("Ensemble tag not valid");

    

      
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_strange_1_L.Tag[i_ens]);
 
    //L
    distr_t_list  V_strange_1_L_distr, V_strange_2_L_distr, V_strange_3_L_distr, V_strange_L_distr;
    distr_t_list  V_strange_OS_1_L_distr, V_strange_OS_2_L_distr, V_strange_OS_3_L_distr, V_strange_OS_L_distr;
    distr_t  Metas_L_distr, Metas_OS_L_distr, MV_L_distr, MV_OS_L_distr;
    //M
    distr_t_list  V_strange_1_M_distr, V_strange_2_M_distr, V_strange_3_M_distr, V_strange_M_distr;
    distr_t_list  V_strange_OS_1_M_distr, V_strange_OS_2_M_distr, V_strange_OS_3_M_distr, V_strange_OS_M_distr;
    distr_t Metas_M_distr, Metas_OS_M_distr, MV_M_distr, MV_OS_M_distr;

     
    

    //vector light sector
    //L
    V_strange_1_L_distr = Corr.corr_t(V_strange_1_L.col(0)[i_ens], "");
    V_strange_2_L_distr = Corr.corr_t(V_strange_2_L.col(0)[i_ens], "");
    V_strange_3_L_distr = Corr.corr_t(V_strange_3_L.col(0)[i_ens], "");
    V_strange_OS_1_L_distr = Corr.corr_t(V_strange_OS_1_L.col(0)[i_ens], "");
    V_strange_OS_2_L_distr = Corr.corr_t(V_strange_OS_2_L.col(0)[i_ens], "");
    V_strange_OS_3_L_distr = Corr.corr_t(V_strange_OS_3_L.col(0)[i_ens], "");
    //M
    V_strange_1_M_distr = Corr.corr_t(V_strange_1_M.col(0)[i_ens], "");
    V_strange_2_M_distr = Corr.corr_t(V_strange_2_M.col(0)[i_ens], "");
    V_strange_3_M_distr = Corr.corr_t(V_strange_3_M.col(0)[i_ens], "");
    V_strange_OS_1_M_distr = Corr.corr_t(V_strange_OS_1_M.col(0)[i_ens], "");
    V_strange_OS_2_M_distr = Corr.corr_t(V_strange_OS_2_M.col(0)[i_ens], "");
    V_strange_OS_3_M_distr = Corr.corr_t(V_strange_OS_3_M.col(0)[i_ens], "");

   



    string strange_tag = "strange";

    

    //P5P5
    //L
    Metas_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_L.col(0)[i_ens], ""));
    Metas_OS_L_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_OS_L.col(0)[i_ens], ""));
    //M
    Metas_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_M.col(0)[i_ens], ""));
    Metas_OS_M_distr = Corr.Fit_distr(Corr.effective_mass_t(pt2_etaS_OS_M.col(0)[i_ens], ""));


     
      
      


      
    //sum over the Lorentz indices of the e.m. current
    //L
    V_strange_L_distr= (pow(qs,2)/3.0)*(V_strange_1_L_distr+ V_strange_2_L_distr + V_strange_3_L_distr);
    V_strange_OS_L_distr =  (pow(qs,2)/3.0)*(V_strange_OS_1_L_distr+ V_strange_OS_2_L_distr + V_strange_OS_3_L_distr);
    //M
    V_strange_M_distr= (pow(qs,2)/3.0)*(V_strange_1_M_distr+ V_strange_2_M_distr + V_strange_3_M_distr);
    V_strange_OS_M_distr =  (pow(qs,2)/3.0)*(V_strange_OS_1_M_distr+ V_strange_OS_2_M_distr + V_strange_OS_3_M_distr);



    //get MV from vector correlator
    //L
    Corr.Tmin = Tmin_VV; Corr.Tmax= Tmax_VV;
    MV_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_L_distr, ""));
    Corr.Tmin = Tmin_VV_OS; Corr.Tmax = Tmax_VV_OS;
    MV_OS_L_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_OS_L_distr, ""));
    //M
    Corr.Tmin = Tmin_VV; Corr.Tmax = Tmax_VV;
    MV_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_M_distr, ""));
    Corr.Tmin = Tmin_VV_OS; Corr.Tmax = Tmax_VV_OS;
    MV_OS_M_distr = Corr.Fit_distr( Corr.effective_mass_t(V_strange_OS_M_distr, ""));


    vector<distr_t> Metas_vec({ Metas_L_distr/a_distr, Metas_M_distr/a_distr});
    vector<distr_t> Mphi_vec({ MV_L_distr/a_distr, MV_M_distr/a_distr});
    vector<distr_t> X_fit;
    if(Extrapolation_strange_mode == "phi") X_fit = Mphi_vec;
    else if(Extrapolation_strange_mode == "etas") X_fit = Metas_vec;
    else crash("Extrapolation strange mode: "+Extrapolation_strange_mode+" is invalid");
    distr_t X_phys = (Extrapolation_strange_mode=="phi")?Mphi_phys_distr:Metas_phys_distr;



    //############################################################################################


    //perturbative substraction of lattice artifacts

    //free corr LO artifacts
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
    string Pt_free_oppor_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_L,3)+"/OPPOR";
    string Pt_free_oppor_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_M,3)+"/OPPOR";
    Vfloat VV_free_oppor_L= Read_From_File(Pt_free_oppor_L, 1, 4);
    Vfloat VV_free_oppor_M= Read_From_File(Pt_free_oppor_M, 1, 4);
    if(VV_free_oppor_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L w opposite r");
    if(VV_free_oppor_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M w opposite r");
    //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
    string Pt_free_samer_L= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_L,3)+"/SAMER";
    string Pt_free_samer_M= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(2*L_info.ms_M,3)+"/SAMER";
    Vfloat VV_free_samer_L= Read_From_File(Pt_free_samer_L, 1, 4);
    Vfloat VV_free_samer_M= Read_From_File(Pt_free_samer_M, 1, 4);
    if(VV_free_samer_L.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_L  w same r");
    if(VV_free_samer_M.size() != Corr.Nt) crash("Failed to read properly free VV correlator mc_M  w same r");
    //Insert electric charges
    for( auto & OP:VV_free_oppor_L) OP *= pert_corr_strange_on_off*qs*qs;
    for( auto & SA:VV_free_samer_L) SA *= pert_corr_strange_on_off*qs*qs;
    for( auto & OP:VV_free_oppor_M) OP *= pert_corr_strange_on_off*qs*qs;
    for( auto & SA:VV_free_samer_M) SA *= pert_corr_strange_on_off*qs*qs;
  
    //free corr LO artifacts
    Vfloat free_corr_log_art(Corr.Nt, 0.0);
    for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave() < add_pert_corr_strange_up_to*fm_to_inv_Gev && t != 0) {   free_corr_log_art[t] = -1.0*pert_corr_strange_on_off*(qs*qs)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;

      if(t==0 || t*a_distr.ave() > add_pert_corr_strange_up_to*fm_to_inv_Gev) { VV_free_samer_L[t] =0; VV_free_samer_M[t] = 0;  VV_free_oppor_L[t] = 0; VV_free_oppor_M[t]=0;}


    }


    //multiply corr using Zv and Za
    //L
    V_strange_L_distr = (V_strange_L_distr)*Za*Za + VV_free_oppor_L ;
    V_strange_OS_L_distr = (V_strange_OS_L_distr)*Zv*Zv + VV_free_samer_L;
    //M
    V_strange_M_distr = (V_strange_M_distr)*Za*Za + VV_free_oppor_M ;
    V_strange_OS_M_distr = (V_strange_OS_M_distr)*Zv*Zv + VV_free_samer_M;



     
      

    
      
      


    //print correlator

    //L
    Print_To_File({}, {V_strange_L_distr.ave(), V_strange_L_distr.err(), V_strange_OS_L_distr.ave(), V_strange_OS_L_distr.err()}, "../data/R_ratio/corr/"+strange_tag+"/corr_L_"+V_strange_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");
    //M
    Print_To_File({}, {V_strange_M_distr.ave(), V_strange_M_distr.err(), V_strange_OS_M_distr.ave(), V_strange_OS_M_distr.err()}, "../data/R_ratio/corr/"+strange_tag+"/corr_M_"+V_strange_1_L.Tag[i_ens]+".dat", "", "# t  tm  OS");


    cout<<"Perturbative corrections added"<<endl;





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    //L
    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax_L=1;
    while(!Found_error_less_x_percent && tmax_L < Corr.Nt/2) {
   
      if( (V_strange_L_distr.distr_list[tmax_L]).err()/fabs( (V_strange_L_distr.distr_list[tmax_L]).ave()) <  0.01*x) tmax_L++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_L=1;
    while(!Found_error_less_x_percent && tmax_OS_L < Corr.Nt/2) {
   
      if( (V_strange_OS_L_distr.distr_list[tmax_OS_L]).err()/fabs( (V_strange_OS_L_distr.distr_list[tmax_OS_L]).ave()) <  0.01*x) tmax_OS_L++;
      else Found_error_less_x_percent=true;
    }

    //M
    Found_error_less_x_percent=false;
    double tmax_M=1;
    while(!Found_error_less_x_percent && tmax_M < Corr.Nt/2) {
   
      if( (V_strange_M_distr.distr_list[tmax_M]).err()/fabs( (V_strange_M_distr.distr_list[tmax_M]).ave()) <  0.01*x) tmax_M++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS_M=1;
    while(!Found_error_less_x_percent && tmax_OS_M < Corr.Nt/2) {
   
      if( (V_strange_OS_M_distr.distr_list[tmax_OS_M]).err()/fabs( (V_strange_OS_M_distr.distr_list[tmax_OS_M]).ave()) <  0.01*x) tmax_OS_M++;
      else Found_error_less_x_percent=true;
    }
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double Emax= max_energy*a_distr.ave();  //[x GeV]*a
    double E0 =  Emin*a_distr.ave(); // two pion treshold ~ 280 MeV
    double Estart= Emin*a_distr.ave(); // start from 300 MeV
    double sigma = sm_size*a_distr.ave(); // smearing
    int Npoints = (int)((Emax - Estart)/(sigma*step_size));
    tmax_L= (T/2 -4);
    tmax_OS_L = T/2 -4;
    tmax_M= (T/2 -4);
    tmax_OS_M = T/2 -4;
    double mult_SANF=0.000001;
    double mult_TANT=0.6;
    //#########################################################################################


    cout<<"Reconstructing R(E)....."<<endl;
    cout<<"aEmax: "<<Emax<<" -> Emax: "<<Emax/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"asigma: "<<sigma<<" -> sigma: "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"Npoints: "<<Npoints<<endl;
    cout<<"Zv: "<<L_info.Zv_WI<< " +- "<<L_info.Zv_WI_err<<endl;
    cout<<"Za: "<<L_info.Za_WI<< " +- "<<L_info.Za_WI_err<<endl;
    cout<<"SM_TYPE: "<<SM_TYPE<<endl;
    

  
    distr_t_list Spectral_dens_L;
    distr_t_list Spectral_dens_OS_L;
    distr_t_list Spectral_dens_SANF_L;
    distr_t_list Spectral_dens_OS_SANF_L;
    distr_t_list Spectral_dens_M;
    distr_t_list Spectral_dens_OS_M;
    distr_t_list Spectral_dens_SANF_M;
    distr_t_list Spectral_dens_OS_SANF_M;
    distr_t_list Spectral_dens_Extr;
    distr_t_list Spectral_dens_OS_Extr;
    distr_t_list Spectral_dens_SANF_Extr;
    distr_t_list Spectral_dens_OS_SANF_Extr;
    Vfloat Ergs;
    Vfloat Ergs_GeV;
    Vfloat syst_tm_L(Npoints);
    Vfloat syst_OS_L(Npoints);
    Vfloat syst_tm_SANF_L(Npoints);
    Vfloat syst_OS_SANF_L(Npoints);
    Vfloat syst_tm_M(Npoints);
    Vfloat syst_OS_M(Npoints);
    Vfloat syst_tm_SANF_M(Npoints);
    Vfloat syst_OS_SANF_M(Npoints);


    //COMPUTE THE SMEARED R-RATIO

    for(int ip=0; ip<Npoints;ip++) {

      double mean = Estart + ip*((Emax-Estart)/((double)Npoints));
      Ergs.push_back(mean);
      double mean_GeV= mean/a_distr.ave();
      Ergs_GeV.push_back(mean_GeV);
      double sigma_resc=sigma;
      double lambda_Estar;
      double lambda_Estar_SANF;
 
      cout<<"Evaluating R(E*) at aE*: "<<mean<<" -> "<<"E*: "<<mean/a_distr.ave()<<" [GeV] ..."<<endl;
      cout<<"Rescaled asigma(E*): "<<sigma_resc<<" -> sigma(E*): "<<sigma_resc/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"[tmin, tmax] (tm L): ["<<1<<", "<<tmax_L<<"]"<<endl;
      cout<<"[tmin, tmax] (OS L): ["<<1<<", "<<tmax_OS_L<<"]"<<endl;
      cout<<"[tmin, tmax] (tm M): ["<<1<<", "<<tmax_M<<"]"<<endl;
      cout<<"[tmin, tmax] (OS M): ["<<1<<", "<<tmax_OS_M<<"]"<<endl;

      //L
      Spectral_dens_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_L, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_L_distr, syst_tm_L[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+strange_tag+"_L_"+V_strange_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm L (TANTALO): "<<Spectral_dens_L.ave(ip)<<" +- ("<<Spectral_dens_L.err(ip)<<")_stat ("<<syst_tm_L[ip]<<")_syst "<<endl;
      Spectral_dens_SANF_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_L, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_L_distr, syst_tm_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+strange_tag+"_L_"+V_strange_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm L (SANFILIPPO): "<<Spectral_dens_SANF_L.ave(ip)<<" +- ("<<Spectral_dens_SANF_L.err(ip)<<")_stat ("<<syst_tm_SANF_L[ip]<<")_syst "<<endl;
      Spectral_dens_OS_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_L, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_OS_L_distr, syst_OS_L[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+strange_tag+"_L_"+V_strange_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS L (TANTALO): "<<Spectral_dens_OS_L.ave(ip)<<" +- ("<<Spectral_dens_OS_L.err(ip)<<")_stat ("<<syst_OS_L[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF_L.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_L, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_OS_L_distr, syst_OS_SANF_L[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+strange_tag+"_L_"+V_strange_1_L.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS L (SANFILIPPO): "<<Spectral_dens_OS_SANF_L.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF_L.err(ip)<<")_stat ("<<syst_OS_SANF_L[ip]<<")_syst "<<endl;
      //M
      Spectral_dens_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_M, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_M_distr, syst_tm_M[ip], mult_TANT, lambda_Estar, "TANT", "tm", "VV_"+strange_tag+"_M_"+V_strange_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm M (TANTALO): "<<Spectral_dens_M.ave(ip)<<" +- ("<<Spectral_dens_M.err(ip)<<")_stat ("<<syst_tm_M[ip]<<")_syst "<<endl;
      Spectral_dens_SANF_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_M, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_M_distr, syst_tm_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "tm", "VV_"+strange_tag+"_M_"+V_strange_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) tm M (SANFILIPPO): "<<Spectral_dens_SANF_M.ave(ip)<<" +- ("<<Spectral_dens_SANF_M.err(ip)<<")_stat ("<<syst_tm_SANF_M[ip]<<")_syst "<<endl;
      Spectral_dens_OS_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_M, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_OS_M_distr, syst_OS_M[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+strange_tag+"_M_"+V_strange_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS M (TANTALO): "<<Spectral_dens_OS_M.ave(ip)<<" +- ("<<Spectral_dens_OS_M.err(ip)<<")_stat ("<<syst_OS_M[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF_M.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS_M, prec, SM_TYPE+"_ov_E2",f, rho_R*V_strange_OS_M_distr, syst_OS_SANF_M[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+strange_tag+"_M_"+V_strange_1_M.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) OS M (SANFILIPPO): "<<Spectral_dens_OS_SANF_M.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF_M.err(ip)<<")_stat ("<<syst_OS_SANF_M[ip]<<")_syst "<<endl;


      //Extrapolate to the physical strange mass
      vector<distr_t> Y_fit_tm = {Spectral_dens_L.distr_list[ip], Spectral_dens_M.distr_list[ip]};
      vector<distr_t> Y_fit_OS = {Spectral_dens_OS_L.distr_list[ip], Spectral_dens_OS_M.distr_list[ip]};
      vector<distr_t> Y_fit_SANF_tm = {Spectral_dens_SANF_L.distr_list[ip], Spectral_dens_SANF_M.distr_list[ip]};
      vector<distr_t> Y_fit_SANF_OS = {Spectral_dens_OS_SANF_L.distr_list[ip], Spectral_dens_OS_SANF_M.distr_list[ip]};

      Spectral_dens_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_tm, X_fit, X_phys ,  "../data/R_ratio/"+strange_tag  , "R_extr_tm_TANT_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_OS_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_OS, X_fit, X_phys ,  "../data/R_ratio/"+strange_tag  , "R_extr_OS_TANT_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_SANF_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_SANF_tm, X_fit, X_phys ,  "../data/R_ratio/"+strange_tag  , "R_extr_tm_SANF_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
      Spectral_dens_OS_SANF_Extr.distr_list.push_back(  Obs_extrapolation_meson_mass( Y_fit_SANF_OS, X_fit, X_phys ,  "../data/R_ratio/"+strange_tag  , "R_extr_OS_SANF_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+"_Estar_"+to_string_with_precision(mean_GeV,3)+".dat",  UseJack, "SPLINE" ));
	
	

     
      cout<<"...done"<<endl;
      //############################################################################################



    }


    RE_strange_TANT_tm[i_ens] = Spectral_dens_Extr;
    RE_strange_TANT_OS[i_ens] = Spectral_dens_OS_Extr;
    RE_strange_SANF_tm[i_ens] = Spectral_dens_SANF_Extr;
    RE_strange_SANF_OS[i_ens] = Spectral_dens_OS_SANF_Extr;
   

      
    cout<<"printing output strange..."<<endl;

    //print to file
    //L
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_L.ave(), Spectral_dens_L.err(),  Spectral_dens_SANF_L.ave(), Spectral_dens_SANF_L.err(), Spectral_dens_OS_L.ave(), Spectral_dens_OS_L.err(), Spectral_dens_OS_SANF_L.ave(), Spectral_dens_OS_SANF_L.err()}, "../data/R_ratio/"+strange_tag+"/L_"+SM_TYPE+"_"+V_strange_1_L.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");
    //M
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_M.ave(), Spectral_dens_M.err(), Spectral_dens_SANF_M.ave(), Spectral_dens_SANF_M.err(), Spectral_dens_OS_M.ave(), Spectral_dens_OS_M.err(), Spectral_dens_OS_SANF_M.ave(), Spectral_dens_OS_SANF_M.err()}, "../data/R_ratio/"+strange_tag+"/M_"+SM_TYPE+"_"+V_strange_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");
    //Extr
    Print_To_File({}, {Ergs, Ergs_GeV,  Spectral_dens_Extr.ave(), Spectral_dens_Extr.err(), Spectral_dens_SANF_Extr.ave(), Spectral_dens_SANF_Extr.err(), Spectral_dens_OS_Extr.ave(), Spectral_dens_OS_Extr.err(), Spectral_dens_OS_SANF_Extr.ave(), Spectral_dens_OS_SANF_Extr.err() }, "../data/R_ratio/"+strange_tag+"/Extr_"+SM_TYPE+"_"+V_strange_1_M.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");

    cout<<"done!"<<endl;


  }



 


  cout<<"Calculation of strange contribution COMPLETED"<<endl;


  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  
  //################## CALCULATION OF THE DISCONNECTED CONTRIBUTIONS ############################

  
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################
  //#############################################################################################

 
  #pragma omp parallel for
  for(int i_ens=0;i_ens<Nens_disco_light;i_ens++) {

    Disco_Ens_Tags.push_back(disco_light.Tag[i_ens]);

    string FLAV= "disco";
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = disco_light.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<disco_light.Tag[i_ens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    if(disco_light.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A; Zv= ZV_A; Za = ZA_A;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B; Zv= ZV_B; Za = ZA_B;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C; Zv= ZV_C; Za = ZA_C;}
    else if(disco_light.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D; Zv= ZV_D; Za = ZA_D;}
    else crash("lattice spacing distribution for Ens: "+disco_light.Tag[i_ens]+" not found");


    LatticeInfo L_info;
    L_info.LatInfo_new_ens(disco_light.Tag[i_ens]);
 
    distr_t_list  light_light_distr, strange_strange_distr, charm_charm_distr;
    distr_t_list  light_strange_distr, light_charm_distr, strange_charm_distr;

    

    if( (disco_light.Tag[i_ens] != disco_strange.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_charm.Tag[i_ens]) ||  (disco_light.Tag[i_ens] != disco_light_strange.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_light_charm.Tag[i_ens]) || (disco_light.Tag[i_ens] != disco_strange_charm.Tag[i_ens])) crash("disconnected ensembles do not match"); 
    

    //diagonal disco
    light_light_distr = Corr.corr_t(disco_light.col(0)[i_ens], "");
    strange_strange_distr = Corr.corr_t(disco_strange.col(0)[i_ens], "");
    charm_charm_distr = Corr.corr_t(disco_charm.col(0)[i_ens], "");
    //off-diagonal disco
    light_strange_distr = Corr.corr_t(disco_light_strange.col(0)[i_ens], "");
    light_charm_distr = Corr.corr_t(disco_light_charm.col(0)[i_ens], "");
    strange_charm_distr = Corr.corr_t(disco_strange_charm.col(0)[i_ens], "");


    distr_t_list disco_distr = Zv*Zv*( pow(qu+qd,2)*light_light_distr + pow(qs,2)*strange_strange_distr + pow(qc,2)*charm_charm_distr + 2.0*(qu+qd)*qs*light_strange_distr + 2.0*(qu+qd)*qc*light_charm_distr + 2.0*qs*qd*strange_charm_distr);

    string disco_tag= "disco";



     //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (disco_distr.distr_list[tmax]).err()/fabs( (disco_distr.distr_list[tmax]).ave()) <  0.01*x) tmax++;
      else Found_error_less_x_percent=true;
    }

  
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double Emax= max_energy*a_distr.ave();  //[x GeV]*a
    double E0 = Emin*a_distr.ave(); // two pion treshold ~ 280 MeV
    double Estart= Emin*a_distr.ave(); // start from 300 MeV
    double sigma = sm_size*a_distr.ave(); // smearing
    
    int Npoints = (int)((Emax - Estart)/(sigma*step_size));
    tmax= (T/2 -4);
    double mult_TANT= 1.0;
    double mult_SANF = 0.00001;
    //#########################################################################################
    
    
    cout<<"Reconstructing R(E)....."<<endl;
    cout<<"aEmax: "<<Emax<<" -> Emax: "<<Emax/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"asigma: "<<sigma<<" -> sigma: "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"Npoints: "<<Npoints<<endl;
    cout<<"Corr(0): "<<(rho_R*disco_distr).ave(0)<<endl;
    cout<<"Zv: "<<L_info.Zv_WI<< " +- "<<L_info.Zv_WI_err<<endl;
    cout<<"SM_TYPE: "<<SM_TYPE<<endl;
    
    
  
    

    distr_t_list Spectral_dens;
    distr_t_list Spectral_dens_SANF;
    Vfloat Ergs;
    Vfloat Ergs_GeV;
    Vfloat lambdas;
    Vfloat syst(Npoints);
    Vfloat syst_SANF(Npoints);
    
    
    //COMPUTE THE SMEARED R-RATIO

    for(int ip=0; ip<Npoints;ip++) {
      
      double mean = Estart + ip*((Emax-Estart)/((double)Npoints));
      Ergs.push_back(mean);
      Ergs_GeV.push_back(mean/a_distr.ave());
      double sigma_resc=sigma;
      double lambda_Estar;
      double lambda_Estar_SANF;
      
      cout<<"Evaluating R(E*) on Ens "<<disco_light.Tag[i_ens]<<" at aE*: "<<mean<<" -> "<<"E*: "<<mean/a_distr.ave()<<" [GeV] ..."<<endl;
      cout<<"Rescaled asigma(E*): "<<sigma_resc<<" -> sigma(E*): "<<sigma_resc/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"[tmin, tmax] (tm): ["<<1<<", "<<tmax<<"]"<<endl;
      Spectral_dens.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*disco_distr, syst[ip], mult_TANT, lambda_Estar, "TANT", "OS", "VV_"+disco_tag+"_"+disco_light.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) (TANTALO): "<<Spectral_dens.ave(ip)<<" +- ("<<Spectral_dens.err(ip)<<")_stat ("<<syst[ip]<<")_syst "<<endl;
      Spectral_dens_SANF.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*disco_distr, syst_SANF[ip], mult_SANF, lambda_Estar_SANF, "SANF", "OS", "VV_"+disco_tag+"_"+disco_light.Tag[i_ens], FLAV )) ;
      cout<<"R(E*) (SANFILIPPO): "<<Spectral_dens_SANF.ave(ip)<<" +- ("<<Spectral_dens_SANF.err(ip)<<")_stat ("<<syst_SANF[ip]<<")_syst "<<endl;
      lambdas.push_back(lambda_Estar);
         
      cout<<"...done"<<endl;
      //############################################################################################



    }


    RE_disco_TANT[i_ens] = Spectral_dens;
    RE_disco_SANF[i_ens] = Spectral_dens_SANF;
    

    
    cout<<"printing output light..."<<endl;

    //light
    //print to file
    Print_To_File({}, {Ergs, Ergs_GeV, Spectral_dens.ave(), Spectral_dens.err(), Spectral_dens_SANF.ave(), Spectral_dens_SANF.err()}, "../data/R_ratio/"+disco_tag+"/"+SM_TYPE+"_"+disco_light.Tag[i_ens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# aE* E*(GeV)   R(E) (stat.) [TANTALO]  R(E) (<> stat.) [SANFILIPPO]");

    cout<<"done!"<<endl;
        





  }
  #pragma omp barrier

  cout<<"Putting all contributions together..."<<endl;

   //summ the light, strange, charm and disconnected contributions
   //optimal result is SANF on light TANT on strange charm and discos


   for(int iens=0; iens<Nens_light;iens++) {

     RE_tot_TANT_tm[iens] = RE_light_TANT_tm[iens];
     RE_tot_TANT_OS[iens] = RE_light_TANT_OS[iens];
     RE_tot_SANF_tm[iens] = RE_light_SANF_tm[iens];
     RE_tot_SANF_OS[iens] = RE_light_SANF_OS[iens];
     RE_tot_opt_tm[iens] =  RE_light_SANF_tm[iens];
     RE_tot_opt_OS[iens] =  RE_light_SANF_OS[iens];
     //find strange
     bool Strange_found=false;
     for(int js_ens=0; js_ens<Nens_strange;js_ens++) {
       if( Light_Ens_Tags[iens] == Strange_Ens_Tags[js_ens]) {
	   RE_tot_TANT_tm[iens] = RE_tot_TANT_tm[iens] + RE_strange_TANT_tm[js_ens];
	   RE_tot_TANT_OS[iens] = RE_tot_TANT_OS[iens] + RE_strange_TANT_OS[js_ens];
	   RE_tot_SANF_tm[iens] = RE_tot_SANF_tm[iens] + RE_strange_SANF_tm[js_ens];
	   RE_tot_SANF_OS[iens] = RE_tot_SANF_OS[iens] + RE_strange_SANF_OS[js_ens];
	   RE_tot_opt_tm[iens] = RE_tot_opt_tm[iens] + RE_strange_SANF_tm[js_ens];
	   RE_tot_opt_OS[iens] = RE_tot_opt_OS[iens] + RE_strange_SANF_OS[js_ens];
	   Strange_found=true;
       }
     }

     //find charm
     bool Charm_found=false;
     int B64_pos = -1;
     for(int jc_ens=0; jc_ens<Nens_charm;jc_ens++) {
       if(Charm_Ens_Tags[jc_ens] == "cB211b.072.64") B64_pos = jc_ens;
       
       if( Light_Ens_Tags[iens] == Charm_Ens_Tags[jc_ens]) {
	   RE_tot_TANT_tm[iens] = RE_tot_TANT_tm[iens] + RE_charm_TANT_tm[jc_ens];
	   RE_tot_TANT_OS[iens] = RE_tot_TANT_OS[iens] + RE_charm_TANT_OS[jc_ens];
	   RE_tot_SANF_tm[iens] = RE_tot_SANF_tm[iens] + RE_charm_SANF_tm[jc_ens];
	   RE_tot_SANF_OS[iens] = RE_tot_SANF_OS[iens] + RE_charm_SANF_OS[jc_ens];
	   RE_tot_opt_tm[iens] = RE_tot_opt_tm[iens] + RE_charm_TANT_tm[jc_ens];
	   RE_tot_opt_OS[iens] = RE_tot_opt_OS[iens] + RE_charm_TANT_OS[jc_ens];
	   Charm_found=true;
       }
     }

     if(Charm_found==false && Light_Ens_Tags[iens] == "cB211b.072.96" && B64_pos != -1 ) {
       
       RE_tot_TANT_tm[iens] = RE_tot_TANT_tm[iens] + RE_charm_TANT_tm[B64_pos];
       RE_tot_TANT_OS[iens] = RE_tot_TANT_OS[iens] + RE_charm_TANT_OS[B64_pos];
       RE_tot_SANF_tm[iens] = RE_tot_SANF_tm[iens] + RE_charm_SANF_tm[B64_pos];
       RE_tot_SANF_OS[iens] = RE_tot_SANF_OS[iens] + RE_charm_SANF_OS[B64_pos];
       RE_tot_opt_tm[iens] = RE_tot_opt_tm[iens] + RE_charm_TANT_tm[B64_pos];
       RE_tot_opt_OS[iens] = RE_tot_opt_OS[iens] + RE_charm_TANT_OS[B64_pos];
       Charm_found=true;
     }

     //find disco
     bool Disco_found=false;
     int B64_pos_disco = -1;
     for(int jd_ens=0; jd_ens<Nens_disco_light;jd_ens++) {
       if(Disco_Ens_Tags[jd_ens] == "cB211b.072.64") B64_pos_disco = jd_ens;
       
       if( Light_Ens_Tags[iens] == Disco_Ens_Tags[jd_ens]) {
	   RE_tot_TANT_tm[iens] = RE_tot_TANT_tm[iens] + RE_disco_TANT[jd_ens];
	   RE_tot_TANT_OS[iens] = RE_tot_TANT_OS[iens] + RE_disco_TANT[jd_ens];
	   RE_tot_SANF_tm[iens] = RE_tot_SANF_tm[iens] + RE_disco_SANF[jd_ens];
	   RE_tot_SANF_OS[iens] = RE_tot_SANF_OS[iens] + RE_disco_SANF[jd_ens];
	   RE_tot_opt_tm[iens] = RE_tot_opt_tm[iens] + RE_disco_SANF[jd_ens];
	   RE_tot_opt_OS[iens] = RE_tot_opt_OS[iens] + RE_disco_SANF[jd_ens];
	   Disco_found=true;
       }
     }

     if(Disco_found==false && Light_Ens_Tags[iens] == "cB211b.072.96" && B64_pos_disco != -1 ) {
       
       RE_tot_TANT_tm[iens] = RE_tot_TANT_tm[iens] + RE_disco_TANT[B64_pos_disco];
       RE_tot_TANT_OS[iens] = RE_tot_TANT_OS[iens] + RE_disco_TANT[B64_pos_disco];
       RE_tot_SANF_tm[iens] = RE_tot_SANF_tm[iens] + RE_disco_SANF[B64_pos_disco];
       RE_tot_SANF_OS[iens] = RE_tot_SANF_OS[iens] + RE_disco_SANF[B64_pos_disco];
       RE_tot_opt_tm[iens] = RE_tot_opt_tm[iens] + RE_disco_SANF[B64_pos_disco];
       RE_tot_opt_OS[iens] = RE_tot_opt_OS[iens] + RE_disco_SANF[B64_pos_disco];
       Charm_found=true;
     }

     if(!Strange_found || !Charm_found  || !Disco_found) crash("Cannot sum light-strange-charm-disconnected contributions. Ensembles do not match");
   }


   cout<<"done!...Printing total..."<<endl;

   //Print the total
   for(int iens=0; iens<(signed)RE_tot_TANT_tm.size();iens++) {

      Print_To_File({}, {Ergs_GeV_list,  RE_tot_TANT_tm[iens].ave(), RE_tot_TANT_tm[iens].err(), RE_tot_SANF_tm[iens].ave(), RE_tot_SANF_tm[iens].err(), RE_tot_TANT_OS[iens].ave(), RE_tot_TANT_OS[iens].err(), RE_tot_SANF_OS[iens].ave(), RE_tot_SANF_OS[iens].err() }, "../data/R_ratio/total/"+SM_TYPE+"_"+Light_Ens_Tags[iens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# E*(GeV)   R(E)_tm (<> stat.) [TANTALO]  R(E)_tm (<> stat.) [SANFILIPPO]   R(E)_OS (<> stat.) [TANTALO]  R(E)_OS (<> stat.) [SANFILIPPO]");

       Print_To_File({}, {Ergs_GeV_list,  RE_tot_opt_tm[iens].ave(), RE_tot_opt_tm[iens].err(),  RE_tot_opt_OS[iens].ave(), RE_tot_opt_OS[iens].err() }, "../data/R_ratio/total/optimal_"+SM_TYPE+"_"+Light_Ens_Tags[iens]+"_sigma_"+to_string_with_precision(sm_size,3)+".dat", "", "# E*(GeV)   R(E)_tm (<> stat.)  R(E)_OS (<> stat.)");

   }



 

   
   cout<<"Bye!"<<endl;

    
  
  return;
}
