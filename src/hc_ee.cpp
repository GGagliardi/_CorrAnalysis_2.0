#include "../include/hc_ee.h"
#include "numerics.h"
#include "stat.h"
#include <boost/optional/detail/optional_reference_spec.hpp>

const double alpha = 1.0 / 137.035999;
const bool UseJack = true;
const int Njacks =  100; //1000;// 93;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const double lambda = 1.24283;
const double Qb = -1.0 / 3.0;
const double Qc = 2.0 / 3.0;
const double m_etac = 2.9841;
const double m_hc=3.52537;
const double m_etab = 9.3987;
const double m_hb = 9.8993;
const string TYPE = "WEAK"; // OR TYPE=="STRONG"
const double mc_3GeV_MS = 0.990;
const double syst_mc_A=0.0022;
const double syst_mc_B=0.0019;
const double syst_mc_C= 0.0022;
const double syst_mc_D = 0.0025;
const double syst_mc_E = 0.002;
const double k2= (m_hc*m_hc - m_etac*m_etac)/(m_hc);

using namespace std;


void Get_ee_plateaux_int_2pt_H_loc(string Ens, CorrAnalysis &Corr) {

  if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
  else if(Ens=="cD211a.054.96") {Corr.Tmin = 22; Corr.Tmax = 26;}
  else if(Ens=="cE211a.044.112") {Corr.Tmin = 25; Corr.Tmax = 29;}
  else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
  else if(Ens=="cA211a.12.48") {Corr.Tmin = 16; Corr.Tmax = 19;}
  else crash("Ens: "+Ens+" not found");
  
  return;
}

void Get_ee_plateaux_int_2pt_H_sm(string Ens, CorrAnalysis &Corr) {

    if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
    else if(Ens=="cD211a.054.96") {Corr.Tmin = 21; Corr.Tmax = 26;}
    else if(Ens=="cE211a.044.112") { Corr.Tmin=25; Corr.Tmax=30;}  //24-30
    else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
    else if(Ens=="cA211a.12.48") {Corr.Tmin = 16; Corr.Tmax = 19;}
    else crash("Ens: "+Ens+" not found");

    return;
}








void hc_ee() {



  


  if(TYPE != "WEAK" && TYPE != "STRONG") crash("TYPE must be WEAK or STRONG, chosen: "+TYPE);



  //#############################################################################################################

  auto Sort_light_confs = [](string A, string B) {

			   

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			     
			    string conf_num_A = A_bis.substr(0,4);
			    string conf_num_B = B_bis.substr(0,4);
							       
		      
			    string rA = A_bis.substr(A_bis.length()-2);
			    string rB = B_bis.substr(B_bis.length()-2);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(A_bis.substr(A_bis.length()-1));
			      int n2 = stoi(B_bis.substr(B_bis.length()-1));
			      if(rA == rB) {
			      if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
			      else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
			      else crash("stream not recognized");
			      }
			      else return n1<n2;
			    }
			    return A_bis<B_bis;
  };




    //generate lattice spacings and RCs
   //############################################################################################
   //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
   GaussianMersenne GM(36551294);
   LatticeInfo a_info;
   distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
   distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack), ZV_E(UseJack);
   distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack), ZA_E(UseJack);
   double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
   double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err, ZV_E_ave, ZV_E_err;
   double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err, ZA_E_ave, ZA_E_err;
   a_info.LatInfo_new_ens("cA211a.12.48");
   a_A_ave= a_info.a_from_afp_FLAG;
   a_A_err= a_info.a_from_afp_FLAG_err;
   ZA_A_ave = a_info.Za_WI_strange;
   ZA_A_err = a_info.Za_WI_strange_err;
   ZV_A_ave = a_info.Zv_WI_strange;
   ZV_A_err = a_info.Zv_WI_strange_err;
   a_info.LatInfo_new_ens("cB211b.072.64");
   a_B_ave= a_info.a_from_afp_FLAG;
   a_B_err= a_info.a_from_afp_FLAG_err;
   ZA_B_ave = a_info.Za_WI_FLAG;
   ZA_B_err = a_info.Za_WI_FLAG_err;
   ZV_B_ave = a_info.Zv_WI_FLAG;
   ZV_B_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cC211a.06.80");
   a_C_ave= a_info.a_from_afp_FLAG;
   a_C_err= a_info.a_from_afp_FLAG_err;
   ZA_C_ave = a_info.Za_WI_FLAG;
   ZA_C_err = a_info.Za_WI_FLAG_err;
   ZV_C_ave = a_info.Zv_WI_FLAG;
   ZV_C_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cD211a.054.96");
   a_D_ave= a_info.a_from_afp_FLAG;
   a_D_err= a_info.a_from_afp_FLAG_err;
   ZA_D_ave = a_info.Za_WI_FLAG;
   ZA_D_err = a_info.Za_WI_FLAG_err;
   ZV_D_ave = a_info.Zv_WI_FLAG;
   ZV_D_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cZ211a.077.64");
   a_Z_ave= a_info.a_from_afp_FLAG;
   a_Z_err= a_info.a_from_afp_FLAG_err;
   a_info.LatInfo_new_ens("cE211a.044.112");
   a_E_ave= a_info.a_from_afp_FLAG;
   a_E_err= a_info.a_from_afp_FLAG_err;
   ZA_E_ave = a_info.Za_WI_FLAG;
   ZA_E_err = a_info.Za_WI_FLAG_err;
   ZV_E_ave = a_info.Zv_WI_FLAG;
   ZV_E_err = a_info.Zv_WI_FLAG_err;


   if(UseJack) {
   
     for(int ijack=0;ijack<Njacks;ijack++) {
     
       a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
       a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
       a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
       a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
       a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
       a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
       ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
       ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
       ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
       ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
       ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
       ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
       ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
       ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
       ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err*(1.0/sqrt(Njacks -1.0)));
       ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err*(1.0/sqrt(Njacks -1.0)));
       
     }
   }

   else {

      for(int ijack=0;ijack<Njacks;ijack++) {
     
	a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
	a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
	a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
	a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
	a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err));
	a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err));
	ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
	ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
	ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
	ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
	ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
	ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
	ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
	ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
	ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err);
	ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err);
	
      }
   }

    




  //###########################################################################################################à


  boost::filesystem::create_directory("../data/hc_ee");

  distr_t_list ZV_M1(UseJack);



  vector<distr_t_list> F1;
  vector<distr_t_list> F2;
  distr_t_list metac(UseJack);
  distr_t_list mhc(UseJack);
  distr_t_list a_distr_list(UseJack);

   
  vector<vector<string>> Ens_list(6);
  vector<string> Ens_list_Hc;

  string m="Hc";
  string mm="M1";

  data_t P5P5_sm_loc, A0P5_sm_loc;
  data_t P5P5_sm, P5P5_sm_REST;
  data_t B1B1_sm, B3B3_sm, B2B2_sm, B2B2_MOT_sm, V1V1_sm;
  data_t PT3_B1P5_tw1, PT3_B1P5_tw2;
  data_t B1B1_sm_loc, B2B2_sm_loc, B3B3_sm_loc;
  data_t V1V1_sm_REST;
  data_t V2V2_sm_REST;
  data_t V3V3_sm_REST;

  P5P5_sm_loc.Read("../heavy_radiative_Hc", "mes_contr_2PT_LOC_"+mm, "P5P5", Sort_light_confs);
  A0P5_sm_loc.Read("../heavy_radiative_Hc", "mes_contr_2PT_LOC_"+mm, "A0P5", Sort_light_confs);
  P5P5_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_MOT_"+m, "P5P5", Sort_light_confs);
  P5P5_sm_REST.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "P5P5", Sort_light_confs);
  V1V1_sm_REST.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
  V2V2_sm_REST.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "V2V2", Sort_light_confs);
  V3V3_sm_REST.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "V3V3", Sort_light_confs);
  B2B2_MOT_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_MOT_"+m, "B1B1", Sort_light_confs);
  B2B2_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "B1B1", Sort_light_confs);
  B1B1_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "B2B2", Sort_light_confs);
  B3B3_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "B3B3", Sort_light_confs);
  V1V1_sm.Read("../heavy_radiative_Hc", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
  PT3_B1P5_tw1.Read("../heavy_radiative_Hc", "mes_contr_3PT_"+m+"_tw1", "B1P5", Sort_light_confs);
  PT3_B1P5_tw2.Read("../heavy_radiative_Hc_tw2", "mes_contr_3PT_"+m+"_tw2", "B1P5", Sort_light_confs);
  
  B2B2_sm_loc.Read("../heavy_radiative_Hc", "mes_contr_2PT_LOC_"+mm, "B1B1", Sort_light_confs);
  B1B1_sm_loc.Read("../heavy_radiative_Hc", "mes_contr_2PT_LOC_"+mm, "B2B2", Sort_light_confs);
  B3B3_sm_loc.Read("../heavy_radiative_Hc", "mes_contr_2PT_LOC_"+mm, "B3B3", Sort_light_confs);

  //READ DIFFERENT THETAS

  data_t P5P5_sm_th1, P5P5_sm_th2, P5P5_sm_th3, P5P5_sm_th4, P5P5_sm_th5;
  data_t PT3_B1P5_G1_th1_tw1, PT3_B1P5_G1_th2_tw1, PT3_B1P5_G1_th3_tw1, PT3_B1P5_G1_th4_tw1, PT3_B1P5_G1_th5_tw1;
  data_t PT3_B3P5_G3_th1_tw1, PT3_B3P5_G3_th2_tw1, PT3_B3P5_G3_th3_tw1, PT3_B3P5_G3_th4_tw1, PT3_B3P5_G3_th5_tw1;

  //tw2
  data_t PT3_B1P5_G1_th1_tw2, PT3_B1P5_G1_th2_tw2, PT3_B1P5_G1_th3_tw2, PT3_B1P5_G1_th4_tw2, PT3_B1P5_G1_th5_tw2;
  data_t PT3_B3P5_G3_th1_tw2, PT3_B3P5_G3_th2_tw2, PT3_B3P5_G3_th3_tw2, PT3_B3P5_G3_th4_tw2, PT3_B3P5_G3_th5_tw2;

  P5P5_sm_th1.Read("../hc_ee", "mes_contr_2PT_MOT_Hc_TH1", "P5P5", Sort_light_confs);
  P5P5_sm_th2.Read("../hc_ee", "mes_contr_2PT_MOT_Hc_TH2", "P5P5", Sort_light_confs);
  P5P5_sm_th3.Read("../hc_ee", "mes_contr_2PT_MOT_Hc_TH3", "P5P5", Sort_light_confs);
  P5P5_sm_th4.Read("../hc_ee", "mes_contr_2PT_MOT_Hc_TH4", "P5P5", Sort_light_confs);
  P5P5_sm_th5.Read("../hc_ee", "mes_contr_2PT_MOT_Hc_TH5", "P5P5", Sort_light_confs);

  
  PT3_B1P5_G1_th1_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G1_TH1_tw1", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th2_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G1_TH2_tw1", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th3_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G1_TH3_tw1", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th4_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G1_TH4_tw1", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th5_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G1_TH5_tw1", "B1P5", Sort_light_confs);

  PT3_B3P5_G3_th1_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G3_TH1_tw1", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th2_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G3_TH2_tw1", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th3_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G3_TH3_tw1", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th4_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G3_TH4_tw1", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th5_tw1.Read("../hc_ee", "mes_contr_3PT_Hc_G3_TH5_tw1", "B3P5", Sort_light_confs);

  
  PT3_B1P5_G1_th1_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G1_TH1_tw2", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th2_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G1_TH2_tw2", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th3_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G1_TH3_tw2", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th4_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G1_TH4_tw2", "B1P5", Sort_light_confs);
  PT3_B1P5_G1_th5_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G1_TH5_tw2", "B1P5", Sort_light_confs);

  PT3_B3P5_G3_th1_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G3_TH1_tw2", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th2_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G3_TH2_tw2", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th3_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G3_TH3_tw2", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th4_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G3_TH4_tw2", "B3P5", Sort_light_confs);
  PT3_B3P5_G3_th5_tw2.Read("../hc_ee_tw2", "mes_contr_3PT_Hc_G3_TH5_tw2", "B3P5", Sort_light_confs);
  
  

  //#################################################################################################à

  
  
  
    
    int Nens=P5P5_sm_loc.size;

    for(int iens=0;iens<Nens;iens++) {

      string Ens =P5P5_sm_loc.Tag[iens];

      bool compute_2nd_tw= (Ens=="cB211b.072.64");

      boost::filesystem::create_directory("../data/hc_ee/"+Ens);

      cout<<"########### ANALYZING ENSEMBLE "<<Ens<<" #############"<<endl;


      //########################################## ENSEMBLE INFO ############################################
      int tw1, tw2;
      int L;
      double amc, amh;
      distr_t a_distr(UseJack);
      distr_t ZV_had(UseJack);
      double syst_mc;
      Vfloat moms;
      if(Ens=="cB211b.072.64") {
	tw1=20; tw2=30;
	a_distr=a_B;
	amc= 0.231567;
	ZV_had= ZV_B;
	L=64;
	syst_mc= syst_mc_B;
	moms = {3.076584*M_PI/L, 2.051056*M_PI/L, 1.025528*M_PI/L, (1e-5)*M_PI/L};
      }
      else if(Ens=="cA211a.12.48") {
	tw1=19;
	a_distr=a_A;
	amc= 0.262;
	ZV_had= ZV_A;
	L=48;
	syst_mc= syst_mc_A;
      }
      else if(Ens=="cD211a.054.96") {
	tw1=30; //tw2=30;
	a_distr=a_D;
	amc= 0.164898;
	ZV_had= ZV_D;
	L=96;
	syst_mc= syst_mc_D;
      }
      else if(Ens=="cC211a.06.80") {
	tw1=25; //tw2=30;
	a_distr=a_C;
	amc= 0.1984;
	ZV_had=ZV_C;
	L=80;
	syst_mc= syst_mc_C;
	moms = {3.2999*M_PI/L, 2.199959*M_PI/L, 1.0999797*M_PI/L, (1e-5)*M_PI/L};
      }
      else if(Ens=="cE211a.044.112") {
	tw1=35; //tw2=30;
	a_distr=a_E;
	amc= 0.141255;
	ZV_had=ZV_E;
	L=112;
	syst_mc= syst_mc_E;
      }
      else crash("Error");
      amh= amc*pow(lambda, stoi(mm.substr(1,1))-1);
      //#################################################################################################

      
      //set up statistical analysis
      CorrAnalysis Corr(UseJack,Njacks,Njacks);
      Corr.Nt= P5P5_sm_loc.nrows[iens];
      Corr.Reflection_sign=1;
      Corr.Perform_Nt_t_average=1;

      //print moms

      cout<<"momenta considered except for q^2=0: "<<endl;
      for(int i=0;i<(signed)moms.size();i++) cout<<"n="<<i+1<<": "<<1e3*moms[i]/a_distr.ave()<<" [MeV]"<<endl;



      
      //########################## LOAD 2PT and 3PT FUNCS ################################################
      Corr.Reflection_sign=-1;
      distr_t_list A0P5_sm_loc_distr= Corr.corr_t(A0P5_sm_loc.col(0)[iens], "../data/hc_ee/"+Ens+"/A0P5_sm_loc");
      Corr.Reflection_sign=1;
      distr_t_list P5P5_sm_th0_distr= Corr.corr_t(P5P5_sm.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm_th0");
      distr_t_list P5P5_sm_th1_distr= Corr.corr_t(P5P5_sm_th1.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm_th1");
      distr_t_list P5P5_sm_th2_distr= Corr.corr_t(P5P5_sm_th2.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm_th2");
      distr_t_list P5P5_sm_th3_distr= Corr.corr_t(P5P5_sm_th3.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm_th2");
      distr_t_list P5P5_sm_th4_distr= Corr.corr_t(P5P5_sm_th4.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm_th3");
      distr_t_list P5P5_sm_distr= Corr.corr_t(P5P5_sm_th5.col(0)[iens], "../data/hc_ee/"+Ens+"/P5P5_sm");
      distr_t_list VKVK_sm_REST_distr= Corr.corr_t(summ_master(V1V1_sm_REST.col(0)[iens], V2V2_sm_REST.col(0)[iens], V3V3_sm_REST.col(0)[iens]), "../data/hc_ee/"+Ens+"/VKVK_sm_REST");
      distr_t_list B2B2_sm_distr= Corr.corr_t(B2B2_sm.col(0)[iens], "../data/hc_ee/"+Ens+"/B2B2_sm");
      distr_t_list BKBK_sm_distr= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm.col(0)[iens], B1B1_sm.col(0)[iens], B3B3_sm.col(0)[iens]), "../data/hc_ee/"+Ens+"/BKBK_sm");
      distr_t_list BKBK_sm_loc_distr= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc.col(0)[iens], B1B1_sm_loc.col(0)[iens], B3B3_sm_loc.col(0)[iens]), "../data/hc_ee/"+Ens+"/BKBK_sm_loc");
      distr_t_list V1V1_sm_distr= Corr.corr_t(V1V1_sm.col(0)[iens], "../data/hc_ee/"+Ens+"/V1V1_sm");
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_G1_th0_tw1=Corr.corr_t(PT3_B1P5_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_th0_tw1");
      distr_t_list PT3_G1_th1_tw1=Corr.corr_t(PT3_B1P5_G1_th1_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_th1_tw1");
      distr_t_list PT3_G1_th2_tw1=Corr.corr_t(PT3_B1P5_G1_th2_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_th2_tw1");
      distr_t_list PT3_G1_th3_tw1=Corr.corr_t(PT3_B1P5_G1_th3_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_th3_tw1");
      distr_t_list PT3_G1_th4_tw1=Corr.corr_t(PT3_B1P5_G1_th4_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_th4_tw1");
      distr_t_list PT3_G1_tw1=Corr.corr_t(PT3_B1P5_G1_th5_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G1_tw1");

      distr_t_list PT3_G3_th1_tw1=Corr.corr_t(PT3_B3P5_G3_th1_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G3_th1_tw1");
      distr_t_list PT3_G3_th2_tw1=Corr.corr_t(PT3_B3P5_G3_th2_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G3_th2_tw1");
      distr_t_list PT3_G3_th3_tw1=Corr.corr_t(PT3_B3P5_G3_th3_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G3_th3_tw1");
      distr_t_list PT3_G3_th4_tw1=Corr.corr_t(PT3_B3P5_G3_th4_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G3_th4_tw1");
      distr_t_list PT3_G3_tw1=Corr.corr_t(PT3_B3P5_G3_th5_tw1.col(0)[iens], "../data/hc_ee/"+Ens+"/PT3_G3_tw1");

      distr_t_list PT3_G3_th1_sub_tw1=Corr.corr_t(summ_master(PT3_B3P5_G3_th1_tw1.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw1.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw1.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th1_sub_tw1");
      distr_t_list PT3_G3_th2_sub_tw1=Corr.corr_t(summ_master(PT3_B3P5_G3_th2_tw1.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw1.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw1.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th2_sub_tw1");
      distr_t_list PT3_G3_th3_sub_tw1=Corr.corr_t(summ_master(PT3_B3P5_G3_th3_tw1.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw1.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw1.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th3_sub_tw1");
      distr_t_list PT3_G3_th4_sub_tw1=Corr.corr_t(summ_master(PT3_B3P5_G3_th4_tw1.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw1.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw1.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th4_sub_tw1");
            
      
      distr_t_list PT3_G1_th0_tw2(UseJack), PT3_G1_th1_tw2(UseJack), PT3_G1_th2_tw2(UseJack), PT3_G1_th3_tw2(UseJack), PT3_G1_th4_tw2(UseJack), PT3_G1_tw2(UseJack);
      distr_t_list PT3_G3_th0_tw2(UseJack), PT3_G3_th1_tw2(UseJack), PT3_G3_th2_tw2(UseJack), PT3_G3_th3_tw2(UseJack), PT3_G3_th4_tw2(UseJack), PT3_G3_tw2(UseJack);
      distr_t_list PT3_G3_th1_sub_tw2(UseJack), PT3_G3_th2_sub_tw2(UseJack), PT3_G3_th3_sub_tw2(UseJack), PT3_G3_th4_sub_tw2(UseJack);
      
      
      if(compute_2nd_tw) {
	
	PT3_G1_th0_tw2=Corr.corr_t(PT3_B1P5_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_th0_tw2");
	PT3_G1_th1_tw2=Corr.corr_t(PT3_B1P5_G1_th1_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_th1_tw2");
	PT3_G1_th2_tw2=Corr.corr_t(PT3_B1P5_G1_th2_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_th2_tw2");
	PT3_G1_th3_tw2=Corr.corr_t(PT3_B1P5_G1_th3_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_th3_tw2");
	PT3_G1_th4_tw2=Corr.corr_t(PT3_B1P5_G1_th4_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_th4_tw2");
	PT3_G1_tw2=Corr.corr_t(PT3_B1P5_G1_th5_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G1_tw2");
	
	PT3_G3_th1_tw2=Corr.corr_t(PT3_B3P5_G3_th1_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G3_th1_tw2");
	PT3_G3_th2_tw2=Corr.corr_t(PT3_B3P5_G3_th2_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G3_th2_tw2");
	PT3_G3_th3_tw2=Corr.corr_t(PT3_B3P5_G3_th3_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G3_th3_tw2");
	PT3_G3_th4_tw2=Corr.corr_t(PT3_B3P5_G3_th4_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G3_th4_tw2");
	PT3_G3_tw2=Corr.corr_t(PT3_B3P5_G3_th5_tw2.col(0)[0], "../data/hc_ee/"+Ens+"/PT3_G3_tw2");
	
	PT3_G3_th1_sub_tw2=Corr.corr_t(summ_master(PT3_B3P5_G3_th1_tw2.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw2.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw2.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th1_sub_tw2");
	PT3_G3_th2_sub_tw2=Corr.corr_t(summ_master(PT3_B3P5_G3_th2_tw2.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw2.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw2.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th2_sub_tw2");
	PT3_G3_th3_sub_tw2=Corr.corr_t(summ_master(PT3_B3P5_G3_th3_tw2.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw2.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw2.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th3_sub_tw2");
	PT3_G3_th4_sub_tw2=Corr.corr_t(summ_master(PT3_B3P5_G3_th4_tw2.col(0)[iens], Multiply_Vvector_by_scalar(PT3_B3P5_G3_th5_tw2.col(0)[iens],-1.0), PT3_B1P5_G1_th5_tw2.col(0)[iens]), "../data/hc_ee/"+Ens+"/PT3_G3_th4_sub_tw2");
	
      }
      Corr.Perform_Nt_t_average=1;

     


      //####################    DETERMINE   EFFECTIVE MASSES ############################################
      distr_t_list H_mass_distr= Corr.effective_mass_t(B2B2_sm_distr, "../data/hc_ee/"+Ens+"/H_mass");
      distr_t_list H_mass_averaged_distr= Corr.effective_mass_t(BKBK_sm_distr, "../data/hc_ee/"+Ens+"/H_mass_averaged");
      distr_t_list H_mass_sm_loc_averaged_distr= Corr.effective_mass_t(BKBK_sm_loc_distr, "../data/hc_ee/"+Ens+"/H_mass_sm_loc_averaged");
      distr_t_list eta_mass_th0_distr= Corr.effective_mass_t(P5P5_sm_th0_distr, "../data/hc_ee/"+Ens+"/eta_en_th0");
      distr_t_list eta_mass_th1_distr= Corr.effective_mass_t(P5P5_sm_th1_distr, "../data/hc_ee/"+Ens+"/eta_en_th1");
      distr_t_list eta_mass_th2_distr= Corr.effective_mass_t(P5P5_sm_th2_distr, "../data/hc_ee/"+Ens+"/eta_en_th2");
      distr_t_list eta_mass_th3_distr= Corr.effective_mass_t(P5P5_sm_th3_distr, "../data/hc_ee/"+Ens+"/eta_en_th3");
      distr_t_list eta_mass_th4_distr= Corr.effective_mass_t(P5P5_sm_th4_distr, "../data/hc_ee/"+Ens+"/eta_en_th4");
      distr_t_list eta_mass_rest_distr= Corr.effective_mass_t(P5P5_sm_distr, "../data/hc_ee/"+Ens+"/eta_mass_REST");
      distr_t_list Jpsi_mass_rest_distr= Corr.effective_mass_t(VKVK_sm_REST_distr, "../data/hc_ee/"+Ens+"/Jpsi_mass_REST");
      Print_To_File({} , { (eta_mass_rest_distr/a_distr).ave(), (eta_mass_rest_distr/a_distr).err()},  "../data/hc_ee/"+Ens+"/eta_mass_REST_pu", "", "");
      Print_To_File({} , { (H_mass_averaged_distr/a_distr).ave(), (H_mass_averaged_distr/a_distr).err()},  "../data/hc_ee/"+Ens+"/H_mass_averaged_pu", "", "");
      
           
      //#################################################################################################

      
 
      //determine Hb mass
      Get_ee_plateaux_int_2pt_H_loc(Ens,Corr);
      distr_t H_mass_sm_loc = Corr.Fit_distr(H_mass_sm_loc_averaged_distr);
      Get_ee_plateaux_int_2pt_H_sm(Ens,Corr);
      distr_t H_mass =  Corr.Fit_distr(H_mass_averaged_distr);
      distr_t H_mass_ave= 0.5*(H_mass+H_mass_sm_loc);
      H_mass= H_mass_ave.ave() + (H_mass_ave - H_mass_ave.ave())*sqrt( 1 + pow( 0.5*(H_mass.ave() - H_mass_sm_loc.ave())/H_mass_ave.err(),2));

      
    
      //Print H_mass_fit
      distr_t_list H_mass_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks, UseJack) ;
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	H_mass_fit.distr_list[t] = H_mass;
      }

      Print_To_File({}, { (H_mass_fit/a_distr).ave(), (H_mass_fit/a_distr).err()}, "../data/hc_ee/"+Ens+"/H_mass_fit_pu", "", "");
      Print_To_File({}, { (H_mass_fit).ave(), (H_mass_fit).err()}, "../data/hc_ee/"+Ens+"/H_mass_fit", "", "");

      distr_t_list VV_H = Corr.matrix_element_t(BKBK_sm_distr, "");



          
      if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=22;}
      else if(Ens=="cA211a.12.48") {Corr.Tmin=14; Corr.Tmax=20;}
      else if(Ens=="cD211a.054.96") {Corr.Tmin = 22; Corr.Tmax = 30;}
      else if(Ens=="cC211a.06.80") {Corr.Tmin = 21; Corr.Tmax = 26;}
      else if(Ens=="cE211a.044.112") {Corr.Tmin = 24; Corr.Tmax = 31;}
      else crash("Ens: "+Ens+" not found"); 
   
      distr_t VV_H_fit= Corr.Fit_distr(VV_H);
      distr_t_list H_MEM= EXPT_D(-1.0*H_mass, Corr.Nt)*Corr.Fit_distr(VV_H)/(2*H_mass);
      distr_t_list H_MEM_tw1 = EXP_D(H_mass*tw1)*H_MEM;
      distr_t_list H_MEM_tw2(UseJack);
      if(compute_2nd_tw) {
	H_MEM_tw2 = EXP_D(H_mass*tw2)*H_MEM;
      }


      
      //set time interval and determine eta masses
      if(Ens=="cB211b.072.64") {Corr.Tmin=35; Corr.Tmax=55; }
      else if(Ens=="cA211a.12.48") {Corr.Tmin=26; Corr.Tmax=35; }
      else if(Ens=="cD211a.054.96") {Corr.Tmin = 50; Corr.Tmax = 58;}
      else if(Ens=="cC211a.06.80") {Corr.Tmin = 35; Corr.Tmax = 56;}
      else if(Ens=="cE211a.044.112") {Corr.Tmin = 51; Corr.Tmax = 63;}
      else crash("Ens "+Ens+" not found");
      distr_t eta_mass_th0 = Corr.Fit_distr(eta_mass_th0_distr);
      distr_t eta_mass_th1 = Corr.Fit_distr(eta_mass_th1_distr);
      distr_t eta_mass_th2 = Corr.Fit_distr(eta_mass_th2_distr);
      distr_t eta_mass_th3 = Corr.Fit_distr(eta_mass_th3_distr);
      distr_t eta_mass_th4 = Corr.Fit_distr(eta_mass_th4_distr);
      distr_t eta_mass_rest= Corr.Fit_distr(eta_mass_rest_distr);
      Print_To_File({}, { ((Jpsi_mass_rest_distr- eta_mass_rest)/(a_distr)).ave(), ((Jpsi_mass_rest_distr - eta_mass_rest)/(a_distr)).err() },  "../data/hc_ee/"+Ens+"/hyperfine_splitting_pu", "", "");
            
      //Print eta_mass_rest_fit
      distr_t_list eta_mass_rest_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks,UseJack) ;
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	eta_mass_rest_fit.distr_list[t] = eta_mass_rest;
      }

      Print_To_File({}, { (eta_mass_rest_fit/a_distr).ave(), (eta_mass_rest_fit/a_distr).err()}, "../data/hc_ee/"+Ens+"/eta_mass_rest_fit", "", "");
     
      
      distr_t_list VV_ETA_th0= Corr.matrix_element_t(P5P5_sm_th0_distr, "");
      distr_t ETA_MEM_th0_tw1 = EXP_D(-1.0*eta_mass_th0*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th0_distr, ""))/(2.0*eta_mass_th0);

      distr_t_list VV_ETA_th1= Corr.matrix_element_t(P5P5_sm_th1_distr, "");
      distr_t ETA_MEM_th1_tw1 = EXP_D(-1.0*eta_mass_th1*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th1_distr, ""))/(2.0*eta_mass_th1);
      
      distr_t_list VV_ETA_th2= Corr.matrix_element_t(P5P5_sm_th2_distr, "");
      distr_t ETA_MEM_th2_tw1 = EXP_D(-1.0*eta_mass_th2*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th2_distr, ""))/(2.0*eta_mass_th2);

      distr_t_list VV_ETA_th3= Corr.matrix_element_t(P5P5_sm_th3_distr, "");
      distr_t ETA_MEM_th3_tw1 = EXP_D(-1.0*eta_mass_th3*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th3_distr, ""))/(2.0*eta_mass_th3);

      distr_t_list VV_ETA_th4= Corr.matrix_element_t(P5P5_sm_th4_distr, "");
      distr_t ETA_MEM_th4_tw1 = EXP_D(-1.0*eta_mass_th4*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th4_distr, ""))/(2.0*eta_mass_th4);
      
                  
      
      distr_t ETA_MEM_th0_tw2(UseJack), ETA_MEM_th1_tw2(UseJack),ETA_MEM_th2_tw2(UseJack),ETA_MEM_th3_tw2(UseJack),ETA_MEM_th4_tw2(UseJack);
      if(compute_2nd_tw) {
	ETA_MEM_th0_tw2 = EXP_D(-1.0*eta_mass_th0*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th0_distr, ""))/(2.0*eta_mass_th0);
	ETA_MEM_th1_tw2 = EXP_D(-1.0*eta_mass_th1*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th1_distr, ""))/(2.0*eta_mass_th1);
	ETA_MEM_th2_tw2 = EXP_D(-1.0*eta_mass_th2*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th2_distr, ""))/(2.0*eta_mass_th2);
	ETA_MEM_th3_tw2 = EXP_D(-1.0*eta_mass_th3*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th3_distr, ""))/(2.0*eta_mass_th3);
	ETA_MEM_th4_tw2 = EXP_D(-1.0*eta_mass_th4*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_th4_distr, ""))/(2.0*eta_mass_th4);
      }
      Print_To_File({}, { VV_H.ave(), VV_H.err(), VV_ETA_th0.ave(), VV_ETA_th0.err(),  VV_ETA_th1.ave(), VV_ETA_th1.err(),  VV_ETA_th2.ave(), VV_ETA_th2.err(), VV_ETA_th3.ave(), VV_ETA_th3.err(), VV_ETA_th4.ave(), VV_ETA_th4.err()} ,  "../data/hc_ee/"+Ens+"/VV", "", "");


      
      //#####################       EVALUATE FORM FACTORS    ###################################
      distr_t_list F1_th0_tw1 = ZV_had*PT3_G1_th0_tw1/(ETA_MEM_th0_tw1*H_MEM_tw1*H_mass);
      distr_t_list F1_th1_tw1 = ZV_had*PT3_G1_th1_tw1/(ETA_MEM_th1_tw1*H_MEM_tw1*H_mass);
      distr_t_list F1_th2_tw1 = ZV_had*PT3_G1_th2_tw1/(ETA_MEM_th2_tw1*H_MEM_tw1*H_mass);
      distr_t_list F1_th3_tw1 = ZV_had*PT3_G1_th3_tw1/(ETA_MEM_th3_tw1*H_MEM_tw1*H_mass);
      distr_t_list F1_th4_tw1 = ZV_had*PT3_G1_th4_tw1/(ETA_MEM_th4_tw1*H_MEM_tw1*H_mass);

      distr_t_list G3_MEM_th1_tw1= ZV_had*PT3_G3_th1_sub_tw1/(ETA_MEM_th1_tw1*H_MEM_tw1);
      distr_t_list G3_MEM_th2_tw1= ZV_had*PT3_G3_th2_sub_tw1/(ETA_MEM_th2_tw1*H_MEM_tw1);
      distr_t_list G3_MEM_th3_tw1= ZV_had*PT3_G3_th3_sub_tw1/(ETA_MEM_th3_tw1*H_MEM_tw1);
      distr_t_list G3_MEM_th4_tw1= ZV_had*PT3_G3_th4_sub_tw1/(ETA_MEM_th4_tw1*H_MEM_tw1);
     
           
             
      distr_t_list F1_th0_tw2(UseJack);
      distr_t_list F1_th1_tw2(UseJack);
      distr_t_list F1_th2_tw2(UseJack);
      distr_t_list F1_th3_tw2(UseJack);
      distr_t_list F1_th4_tw2(UseJack);

      distr_t_list G3_MEM_th1_tw2(UseJack), G3_MEM_th2_tw2(UseJack), G3_MEM_th3_tw2(UseJack), G3_MEM_th4_tw2(UseJack);
               
      if(compute_2nd_tw) {
	F1_th0_tw2 = ZV_had*PT3_G1_th0_tw2/(ETA_MEM_th0_tw2*H_MEM_tw2*H_mass);
	F1_th1_tw2 = ZV_had*PT3_G1_th1_tw2/(ETA_MEM_th1_tw2*H_MEM_tw2*H_mass);
	F1_th2_tw2 = ZV_had*PT3_G1_th2_tw2/(ETA_MEM_th2_tw2*H_MEM_tw2*H_mass);
	F1_th3_tw2 = ZV_had*PT3_G1_th3_tw2/(ETA_MEM_th3_tw2*H_MEM_tw2*H_mass);
	F1_th4_tw2 = ZV_had*PT3_G1_th4_tw2/(ETA_MEM_th4_tw2*H_MEM_tw2*H_mass);
	
	G3_MEM_th1_tw2= ZV_had*PT3_G3_th1_sub_tw2/(ETA_MEM_th1_tw2*H_MEM_tw2);
	G3_MEM_th2_tw2= ZV_had*PT3_G3_th2_sub_tw2/(ETA_MEM_th2_tw2*H_MEM_tw2);
	G3_MEM_th3_tw2= ZV_had*PT3_G3_th3_sub_tw2/(ETA_MEM_th3_tw2*H_MEM_tw2);
	G3_MEM_th4_tw2= ZV_had*PT3_G3_th4_sub_tw2/(ETA_MEM_th4_tw2*H_MEM_tw2);
	
      }
      //################################################################################################


      

      //################################# PRINT FORM FACTORS ###########################################
      //reduce
      distr_t_list F1_th0_tw1_RED(UseJack), F1_th1_tw1_RED(UseJack), F1_th2_tw1_RED(UseJack), F1_th3_tw1_RED(UseJack), F1_th4_tw1_RED(UseJack);
      distr_t_list F1_th0_tw2_RED(UseJack), F1_th1_tw2_RED(UseJack), F1_th2_tw2_RED(UseJack), F1_th3_tw2_RED(UseJack), F1_th4_tw2_RED(UseJack);

          
      distr_t_list G3_MEM_th1_tw1_RED(UseJack), G3_MEM_th2_tw1_RED(UseJack), G3_MEM_th3_tw1_RED(UseJack), G3_MEM_th4_tw1_RED(UseJack);
      distr_t_list G3_MEM_th1_tw2_RED(UseJack), G3_MEM_th2_tw2_RED(UseJack), G3_MEM_th3_tw2_RED(UseJack), G3_MEM_th4_tw2_RED(UseJack);

      // < eta | G0 | hc > =  Af*F1 + Bf*F2, we determine now Af and Bf
      // < eta | G3 | hc > =  Cf*F1 + Df*F2, we determine now Cf and Df

      distr_t_list Af(UseJack), Bf(UseJack), Cf(UseJack), Df(UseJack);

      for(int i=0;i<moms.size();i++) {

	assert(i < 4);
	distr_t eta_en(UseJack);
	if(i==0) eta_en = eta_mass_th1;
	else if(i==1) eta_en = eta_mass_th2;
	else if(i==2) eta_en = eta_mass_th3;
	else if(i==3) eta_en = eta_mass_th4;
	
	
	double k = moms[i];
	distr_t m_eta= eta_mass_rest;
	distr_t m_h = H_mass;
	distr_t q0= H_mass - eta_en;
	distr_t tq0= H_mass + eta_en;
	distr_t q2= POW_D(H_mass,2) + POW_D(m_eta,2) - 2*H_mass*eta_en;

	distr_t af= -1.0*m_h*q0/q2 ;
	distr_t bf= (m_h*m_h - m_eta*m_eta)*q0/q2 -tq0;
	distr_t cf = m_h + m_h*k*k/q2;
	distr_t df= (m_h*m_h - m_eta*m_eta)*(-k*k/q2) + k*k;

	Af.distr_list.push_back(1.0*af);
	Bf.distr_list.push_back(1.0*bf);
	Cf.distr_list.push_back(1.0*cf);
	Df.distr_list.push_back(1.0*df);

	cout<<"cf: "<<cf.ave()<<" df: "<<df.ave()<<endl;
	
	
      }

          
      distr_t_list F23_th0_tw1_RED(UseJack), F23_th1_tw1_RED(UseJack), F23_th2_tw1_RED(UseJack), F23_th3_tw1_RED(UseJack), F23_th4_tw1_RED(UseJack);
      distr_t_list F23_th0_tw2_RED(UseJack), F23_th1_tw2_RED(UseJack), F23_th2_tw2_RED(UseJack), F23_th3_tw2_RED(UseJack), F23_th4_tw2_RED(UseJack);
     
      distr_t fact= a_distr*k2/H_mass;
      
      for(int t=tw1;t<Corr.Nt;t++)  {
	F1_th0_tw1_RED.distr_list.push_back( F1_th0_tw1[t]);
	F1_th1_tw1_RED.distr_list.push_back( F1_th1_tw1[t]);
	F1_th2_tw1_RED.distr_list.push_back( F1_th2_tw1[t]);
	F1_th3_tw1_RED.distr_list.push_back( F1_th3_tw1[t]);
	F1_th4_tw1_RED.distr_list.push_back( F1_th4_tw1[t]);




	F23_th0_tw1_RED.distr_list.push_back( F1_th0_tw1[t] );
	F23_th1_tw1_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th1_tw1[t] - Cf[0]*F1_th1_tw1[t])/Df[0]);
	F23_th2_tw1_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th2_tw1[t] - Cf[1]*F1_th2_tw1[t])/Df[1]);
	F23_th3_tw1_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th3_tw1[t] - Cf[2]*F1_th3_tw1[t])/Df[2]);
	F23_th4_tw1_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th4_tw1[t] - Cf[3]*F1_th4_tw1[t])/Df[3]);
	

      }
      if(compute_2nd_tw) {
	for(int t=tw2;t<Corr.Nt;t++) {
	  F1_th0_tw2_RED.distr_list.push_back( F1_th0_tw2[t]);
	  F1_th1_tw2_RED.distr_list.push_back( F1_th1_tw2[t]);
	  F1_th2_tw2_RED.distr_list.push_back( F1_th2_tw2[t]);
	  F1_th3_tw2_RED.distr_list.push_back( F1_th3_tw2[t]);
	  F1_th4_tw2_RED.distr_list.push_back( F1_th4_tw2[t]);


	  F23_th0_tw2_RED.distr_list.push_back( F1_th0_tw2[t] );
	  F23_th1_tw2_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th1_tw2[t] - Cf[0]*F1_th1_tw2[t])/Df[0]);
	  F23_th2_tw2_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th2_tw2[t] - Cf[1]*F1_th2_tw2[t])/Df[1]);
	  F23_th3_tw2_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th3_tw2[t] - Cf[2]*F1_th3_tw2[t])/Df[2]);
	  F23_th4_tw2_RED.distr_list.push_back( fact*H_mass*(G3_MEM_th4_tw2[t] - Cf[3]*F1_th4_tw2[t])/Df[3]);
	  
	}
      }
      
      Print_To_File({}, {F1_th0_tw1_RED.ave(), F1_th0_tw1_RED.err(),F1_th1_tw1_RED.ave(), F1_th1_tw1_RED.err(),F1_th2_tw1_RED.ave(), F1_th2_tw1_RED.err(),F1_th3_tw1_RED.ave(), F1_th3_tw1_RED.err(), F1_th4_tw1_RED.ave(), F1_th4_tw1_RED.err()  }, "../data/hc_ee/"+Ens+"/F1_tw1", "", "");


      Print_To_File({}, {F23_th0_tw1_RED.ave(), F23_th0_tw1_RED.err(), F23_th1_tw1_RED.ave(), F23_th1_tw1_RED.err(),F23_th2_tw1_RED.ave(), F23_th2_tw1_RED.err(),F23_th3_tw1_RED.ave(), F23_th3_tw1_RED.err(), F23_th4_tw1_RED.ave(), F23_th4_tw1_RED.err()  }, "../data/hc_ee/"+Ens+"/F23_tw1", "", "");

      
      if(compute_2nd_tw) {

        Print_To_File({}, { F1_th0_tw2_RED.ave(), F1_th0_tw2_RED.err(),F1_th1_tw2_RED.ave(), F1_th1_tw2_RED.err(),F1_th2_tw2_RED.ave(), F1_th2_tw2_RED.err(),F1_th3_tw2_RED.ave(), F1_th3_tw2_RED.err(), F1_th4_tw2_RED.ave(), F1_th4_tw2_RED.err()}, "../data/hc_ee/"+Ens+"/F1_tw2", "", "");


	Print_To_File({}, {F23_th0_tw2_RED.ave(), F23_th0_tw2_RED.err(), F23_th1_tw2_RED.ave(), F23_th1_tw2_RED.err(),F23_th2_tw2_RED.ave(), F23_th2_tw2_RED.err(),F23_th3_tw2_RED.ave(), F23_th3_tw2_RED.err(), F23_th4_tw2_RED.ave(), F23_th4_tw2_RED.err()  }, "../data/hc_ee/"+Ens+"/F23_tw2", "", "");

      }
      
      //################################################################################################


      Corr.Tmin= (int)(1.1/(a_distr.ave()/fm_to_inv_Gev));
      Corr.Tmax= (int)(1.4/(a_distr.ave()/fm_to_inv_Gev));

      distr_t F1_th0=Corr.Fit_distr(F1_th0_tw1_RED);
      distr_t F1_th1=Corr.Fit_distr(F1_th1_tw1_RED);
      distr_t F1_th2=Corr.Fit_distr(F1_th2_tw1_RED);
      distr_t F1_th3=Corr.Fit_distr(F1_th3_tw1_RED);
      distr_t F1_th4=Corr.Fit_distr(F1_th4_tw1_RED);


      distr_t F2_th0=Corr.Fit_distr(F23_th0_tw1_RED);
      distr_t F2_th1=Corr.Fit_distr(F23_th1_tw1_RED);
      distr_t F2_th2=Corr.Fit_distr(F23_th2_tw1_RED);
      distr_t F2_th3=Corr.Fit_distr(F23_th3_tw1_RED);
      distr_t F2_th4=Corr.Fit_distr(F23_th4_tw1_RED);

      //store in a vector and print

      distr_t_list F1_list(UseJack), F2_list(UseJack);
      F1_list.distr_list.push_back( F1_th0);
      F1_list.distr_list.push_back( F1_th1);
      F1_list.distr_list.push_back( F1_th2);
      F1_list.distr_list.push_back( F1_th3);
      F1_list.distr_list.push_back( F1_th4);
      F2_list.distr_list.push_back( F2_th0);
      F2_list.distr_list.push_back( F2_th1);
      F2_list.distr_list.push_back( F2_th2);
      F2_list.distr_list.push_back( F2_th3);
      F2_list.distr_list.push_back( F2_th4);

     

      cout<<"Ens: "<<Ens<<endl;
      cout<<"F1(th0): "<<F1_th0.ave()<<" +- "<<F1_th0.err()<<endl;
      cout<<"F1(th1): "<<F1_th1.ave()<<" +- "<<F1_th1.err()<<endl;
      cout<<"F1(th2): "<<F1_th2.ave()<<" +- "<<F1_th2.err()<<endl;
      cout<<"F1(th3): "<<F1_th3.ave()<<" +- "<<F1_th3.err()<<endl;
      cout<<"F1(th4): "<<F1_th4.ave()<<" +- "<<F1_th4.err()<<endl;


      Vfloat q2_sqrt_list;
      q2_sqrt_list.push_back( 0.0);

      

      for(int it=0;it<(signed)moms.size();it++) {
	double mom= moms[it]/a_distr.ave();
	double mom_half= mom*1.25;
	double E_etac= sqrt( m_etac*m_etac + mom*mom);
	double E_etac_half= sqrt( m_etac*m_etac + mom_half*mom_half);
	q2_sqrt_list.push_back(sqrt( m_hc*m_hc + m_etac*m_etac - 2*m_hc*E_etac));
	cout<<"it: "<<it<<" sqrt(q2):"<<q2_sqrt_list[it+1]<<" half mom: "<< sqrt(m_hc*m_hc + m_etac*m_etac - 2*m_hc*E_etac_half)<<endl;
      }

      Print_To_File({}, {q2_sqrt_list, F1_list.ave(), F1_list.err(), F2_list.ave(), F2_list.err() },  "../data/hc_ee/"+Ens+"/F_list_q2", "", "");
	       

  
    }

    //print info on FFs 
}
  
