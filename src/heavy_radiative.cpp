#include "../include/heavy_radiative.h"
#include "numerics.h"
#include "stat.h"
#include <boost/optional/detail/optional_reference_spec.hpp>

const double alpha = 1.0 / 137.035999;
const bool UseJack = true;
const int Njacks = 50;
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
const double l= 1.165;
//const Vfloat Mh_masses({m_etac, m_etac, m_etac *l, m_etac*l *l, m_etac *l *l *l, m_etac *l *l *l, m_etac *l *l *l *l});
const Vfloat Mh_masses({m_etac, m_etac, 3.45, 3.99, 4.64, 5.41, 6.32});


using namespace std;


void Get_plateaux_int_2pt_H_loc(string Ens, CorrAnalysis &Corr) {

  if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
  else if(Ens=="cD211a.054.96") {Corr.Tmin = 22; Corr.Tmax = 26;}
  else if(Ens=="cE211a.044.112") {Corr.Tmin = 25; Corr.Tmax = 29;}
  else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
  else if(Ens=="cA211a.12.48") {Corr.Tmin = 15; Corr.Tmax = 21;}
  else crash("Ens: "+Ens+" not found");
  
  return;
}

void Get_plateaux_int_2pt_H_sm(string Ens, CorrAnalysis &Corr) {

    if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
    else if(Ens=="cD211a.054.96") {Corr.Tmin = 21; Corr.Tmax = 26;}
    else if(Ens=="cE211a.044.112") { Corr.Tmin=25; Corr.Tmax=30;}  //24-30
    else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
    else if(Ens=="cA211a.12.48") {Corr.Tmin = 14; Corr.Tmax = 19;}
    else crash("Ens: "+Ens+" not found");

    return;
}

void Get_plateaux_dmc(string Ens, string M,  CorrAnalysis &Corr) {

  if(Ens=="cA211a.12.48") {
    Corr.Tmin=10; Corr.Tmax=14;
  }
  
  else if(Ens=="cB211b.072.64") {
    Corr.Tmin=11; Corr.Tmax=14;
  }
  else if(Ens=="cC211a.06.80") {
    if(M=="M1" || M=="Hc") { Corr.Tmin=9; Corr.Tmax=13 ;}
    if(M=="M2") { Corr.Tmin=9; Corr.Tmax=13;  }
    if(M=="M3") { Corr.Tmin=9; Corr.Tmax=12;  }
    if(M=="M4") { Corr.Tmin=12; Corr.Tmax=14;   }
    if(M=="M5" || M=="M6") { Corr.Tmin=12; Corr.Tmax=14;  }
    
  }
  else if(Ens=="cD211a.054.96") {

    if(M=="M1" || M=="Hc") { Corr.Tmin=13; Corr.Tmax=15 ;}
    if(M=="M2") {  Corr.Tmin=12; Corr.Tmax=15;  }
    if(M=="M3") {  Corr.Tmin=12; Corr.Tmax=15; }
    if(M=="M4") {  Corr.Tmin=12; Corr.Tmax=15; }
    if(M=="M5" || M=="M6") { Corr.Tmin=9; Corr.Tmax=14;   }

  }
  else if(Ens=="cE211a.044.112") {

    if(M=="M1" || M=="Hc") { Corr.Tmin=14; Corr.Tmax=18 ;}
    if(M=="M2") { Corr.Tmin=14; Corr.Tmax=18;  }
    if(M=="M3") { Corr.Tmin=14; Corr.Tmax=18;  }
    if(M=="M4") { Corr.Tmin=14; Corr.Tmax=18;  }
    if(M=="M5" || M=="M6") { Corr.Tmin=14; Corr.Tmax=18;  }

  }
  else crash("CRASH Get_plateaux_dmc");
  
  

    return;
}

void Get_plateaux_int_2pt_dH_loc(string Ens, string M, CorrAnalysis &Corr) {

  if(M=="Hc" || M == "M1") return;

  if(Ens=="cB211b.072.64") {
    if(M=="M2") {Corr.Tmin=19; Corr.Tmax=23;}
    else if(M=="M3") {Corr.Tmin=12; Corr.Tmax=18;}
    else if(M=="M4") {Corr.Tmin=17; Corr.Tmax=19;}
    else if(M=="M5")  {Corr.Tmin=18; Corr.Tmax=22;}
    else if(M=="M6")  {Corr.Tmin=13; Corr.Tmax=20;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cA211a.12.48") {
    if(M=="M2") {Corr.Tmin=12; Corr.Tmax=14;}
    else if(M=="M3") {Corr.Tmin=11; Corr.Tmax=14;}
    else if(M=="M4") {Corr.Tmin=13; Corr.Tmax=18;}
    else if(M=="M5")  {Corr.Tmin=16; Corr.Tmax=20;}
    else if(M=="M6")  {Corr.Tmin=12; Corr.Tmax=18;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cC211a.06.80") {
    if(M=="M2") {Corr.Tmin=21; Corr.Tmax=27;}
    else if(M=="M3") {Corr.Tmin=19; Corr.Tmax=25;}
    else if(M=="M4") {Corr.Tmin=21; Corr.Tmax=25;}
    else if(M=="M5")  {Corr.Tmin=19; Corr.Tmax=22;}
    else if(M=="M6")  {Corr.Tmin=19; Corr.Tmax=26;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cD211a.054.96") {
    if(M=="M2") {Corr.Tmin=21; Corr.Tmax=25; }
    else if(M=="M3") {Corr.Tmin=24; Corr.Tmax=30;}

    else if(M=="M4") {Corr.Tmin=21; Corr.Tmax=23;}
    else if(M=="M5")  {Corr.Tmin=21; Corr.Tmax=26;}
    else if(M=="M6")  {Corr.Tmin=24; Corr.Tmax=30;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cE211a.044.112") {
    if(M=="M2") {Corr.Tmin=24; Corr.Tmax=31;}
    else if(M=="M3") {Corr.Tmin=25; Corr.Tmax=32;}
    else if(M=="M4") {Corr.Tmin=26; Corr.Tmax=29;}
    else if(M=="M5")  {Corr.Tmin=22; Corr.Tmax=23;}
    else if(M=="M6")  {Corr.Tmin=27; Corr.Tmax=29;}
    else crash("ERROR IN PLATEAUX");
  }
  else crash("Ens not found");


  return;
}

void Get_plateaux_int_2pt_dH_sm(string Ens, string M, CorrAnalysis &Corr) {

  if(M=="Hc" || M == "M1") return;

  if(Ens=="cB211b.072.64") {
    if(M=="M2") {Corr.Tmin=15; Corr.Tmax=21;}
    else if(M=="M3") {Corr.Tmin=14; Corr.Tmax=18;}
    else if(M=="M4") {Corr.Tmin=14; Corr.Tmax=18;}
    else if(M=="M5")  {Corr.Tmin=14; Corr.Tmax=22;}
    else if(M=="M6")  {Corr.Tmin=18; Corr.Tmax=22;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cA211a.12.48") {
    if(M=="M2") {Corr.Tmin=10; Corr.Tmax=13;}
    else if(M=="M3") {Corr.Tmin=9; Corr.Tmax=12;}
    else if(M=="M4") {Corr.Tmin=13; Corr.Tmax=16;}
    else if(M=="M5")  {Corr.Tmin=12; Corr.Tmax=20;}
    else if(M=="M6")  {Corr.Tmin=15; Corr.Tmax=20;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cC211a.06.80") {
    if(M=="M2") {Corr.Tmin=17; Corr.Tmax=26;}
    else if(M=="M3") {Corr.Tmin=15; Corr.Tmax=20;}
    else if(M=="M4") {Corr.Tmin=18; Corr.Tmax=25;}
    else if(M=="M5")  {Corr.Tmin=21; Corr.Tmax=23;}
    else if(M=="M6")  {Corr.Tmin=19; Corr.Tmax=26;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cD211a.054.96") {
    if(M=="M2") {Corr.Tmin=24; Corr.Tmax=28; }
    else if(M=="M3") {Corr.Tmin=20; Corr.Tmax=25;}
    else if(M=="M4") {Corr.Tmin=26; Corr.Tmax=30;}
    else if(M=="M5")  {Corr.Tmin=20; Corr.Tmax=28;}
    else if(M=="M6")  {Corr.Tmin=22; Corr.Tmax=30;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cE211a.044.112") {
    if(M=="M2") {Corr.Tmin=21; Corr.Tmax=26;}
    else if(M=="M3") {Corr.Tmin=22; Corr.Tmax=27;}
    else if(M=="M4") {Corr.Tmin=22; Corr.Tmax=27;}
    else if(M=="M5")  {Corr.Tmin=24; Corr.Tmax=28;}
    else if(M=="M6")  {Corr.Tmin=28; Corr.Tmax=32;}
    else crash("ERROR IN PLATEAUX");
  }
  else crash("Ens not found");


  return;
  
}

void Get_plateaux_VEV_ratio(string Ens, string M, CorrAnalysis &Corr) {
  
  if(Ens=="cB211b.072.64") {
    Corr.Tmin=15; Corr.Tmax=25;

    if(M=="M4") {Corr.Tmin=13; Corr.Tmax=19;}
  }
  else if (Ens=="cA211a.12.48") {
    Corr.Tmin=13; Corr.Tmax=22;
  }
  else if(Ens=="cC211a.06.80") {
    Corr.Tmin=17; Corr.Tmax=22;
  }
  else if(Ens=="cD211a.054.96") {
    if(M=="M2") {Corr.Tmin=24;Corr.Tmax=34;}
    else if(M=="M3")  { Corr.Tmin=19; Corr.Tmax=23; }
    else if(M=="M5")  { Corr.Tmin=23; Corr.Tmax=28;}
    else if(M=="M6")  { Corr.Tmin=19; Corr.Tmax=23;}
    else {Corr.Tmin=24;Corr.Tmax=30;}
  }
  else if(Ens=="cE211a.044.112") {
    if(M=="M2") {Corr.Tmin=21; Corr.Tmax=26;}
    else if(M=="M6") { Corr.Tmin=27; Corr.Tmax=32;}
    else if(M=="M5") {Corr.Tmin=24; Corr.Tmax=28;}
    else {Corr.Tmin=23; Corr.Tmax=27;}
  }
  else crash("Ens not found");
  
  
  return;
}





void heavy_radiative() {


  //heavy_radiative_test();


  if(TYPE != "WEAK" && TYPE != "STRONG") crash("TYPE must be WEAK or STRONG, chosen: "+TYPE);

  auto Get_F_R = [](int ir) {


     double L_QCD=Get_Lambda_MS_bar(4);
     double M= pow(lambda,ir)*mc_3GeV_MS; //mc(3GeV) ~ 0.99 GeV FLAG'24
     //solve v= alphas(mv)
     double MRS= MS_bar_to_MRS_mass(3.0, 4, L_QCD, M,  mc_3GeV_MS);
     cout<<"MRS("<<M<<"): "<<MRS<<endl;   
     auto lambda_func= [&MRS,&L_QCD](double x) {

       return Get_4l_alpha_s(MRS*x, 4 , L_QCD) - x;

     };
     double Z=R_brent( lambda_func, 3.0*Get_4l_alpha_s(MRS,4,L_QCD), 0.7*Get_4l_alpha_s(MRS,4,L_QCD));
     cout<<"Z: "<<Z<<endl;
     //double Mm=  m_MS_bar_m(3,4, L_QCD, M);
     double Mm=M;
     cout<<"LQCD: "<<L_QCD<<endl;
     cout<<"m(m,"<<ir<<"): "<<M<<endl;
     cout<<"alpha(mc): "<<Get_4l_alpha_s(1.0, 4 , L_QCD)<<endl;
     cout<<"alpha(mb): "<<Get_4l_alpha_s(4.8, 4 , L_QCD)<<endl;
     double F_R = MRS/M;
     F_R = MRS/M;
     if(TYPE=="WEAK") return F_R*Z; //F_R = Get_4l_alpha_s( Mm, 4,  L_QCD);

     return F_R;
    
  };

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















  //###########################################################################################################à

  vector<string> M({"Hc","M1", "M2", "M3", "M4", "M5", "M6"});

  boost::filesystem::create_directory("../data/heavy_radiative");

  distr_t_list ZV_M1(UseJack);



  distr_t_list Fc(UseJack);
  distr_t_list Fm1(UseJack);
  distr_t_list metac(UseJack);
  distr_t_list mhc(UseJack);
  distr_t_list a_distr_list(UseJack);

  vector<distr_t_list> mom_k_list;
  for(int i=0;i<6;i++) mom_k_list.emplace_back(UseJack);
  vector<distr_t_list> m_etah_list;
  for(int i=0;i<6;i++) m_etah_list.emplace_back(UseJack);
  vector<distr_t_list> m_etah_list_glb;
  for(int i=0;i<6;i++) m_etah_list_glb.emplace_back(UseJack);
  vector<distr_t_list> aH_glb_list;
  for(int i=0;i<6;i++) aH_glb_list.emplace_back(UseJack);
  vector<distr_t_list> aH_list;
  for(int i=0;i<6;i++) aH_list.emplace_back(UseJack);


  //b-quark mass extrapolation

  distr_t_list R1(UseJack);
  distr_t_list R2(UseJack);
  distr_t_list R3(UseJack);
  distr_t_list R4(UseJack);
  distr_t_list R5(UseJack);

  vector<distr_t_list> R;
  for(int i=0;i<5;i++) R.emplace_back(UseJack);
  vector<distr_t_list> aR;
  for(int i=0;i<5;i++) aR.emplace_back(UseJack);


  vector<distr_t_list> C_beta;
  for(int i=0;i<6;i++) C_beta.emplace_back(UseJack);


  vector<distr_t_list> DH_mass_list;
  for(int i=0;i<6;i++) DH_mass_list.emplace_back(UseJack);

  vector<distr_t_list> etah_mass_list;
  for(int i=0;i<6;i++) etah_mass_list.emplace_back(UseJack);

  vector<vector<string>> Ens_list(6);

  vector<distr_t_list> AP_splitting;
  for(int i=0;i<6;i++) AP_splitting.emplace_back(UseJack);

  
  distr_t Metac_phys(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++) Metac_phys.distr.push_back( 2.984 + 0.005*GM()/sqrt(Njacks-1.0)         );

  distr_t_list amc_etac_list(UseJack);
    
  

  for(int im=0;im<(signed)M.size();im++) {

    

    string m=M[im];
    string mm=M[im];
    if (m=="Hc") mm="M1";
    string m_minus;
    if(m=="Hc" || m=="M1") m_minus="M1";
    else m_minus=M[im-1];
    
    boost::filesystem::create_directory("../data/heavy_radiative/"+m);
    
    
    data_t P5P5_sm_loc, A0P5_sm_loc;
    data_t P5P5_sm, P5P5_sm_REST;
    data_t B1B1_sm, B3B3_sm, B2B2_sm, B2B2_MOT_sm, V1V1_sm;
    data_t PT3_B2P5_tw1, PT3_B2P5_tw2;
    data_t B1B1_sm_loc, B2B2_sm_loc, B3B3_sm_loc;
    data_t V1V1_sm_REST;
    data_t V2V2_sm_REST;
    data_t V3V3_sm_REST;

    P5P5_sm_loc.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+mm, "P5P5", Sort_light_confs);
    A0P5_sm_loc.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+mm, "A0P5", Sort_light_confs);
    P5P5_sm.Read("../heavy_radiative", "mes_contr_2PT_MOT_"+m, "P5P5", Sort_light_confs);
    P5P5_sm_REST.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "P5P5", Sort_light_confs);
    V1V1_sm_REST.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    V2V2_sm_REST.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "V2V2", Sort_light_confs);
    V3V3_sm_REST.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "V3V3", Sort_light_confs);
    B2B2_MOT_sm.Read("../heavy_radiative", "mes_contr_2PT_MOT_"+m, "B1B1", Sort_light_confs);
    B2B2_sm.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "B3B3", Sort_light_confs);
    V1V1_sm.Read("../heavy_radiative", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1.Read("../heavy_radiative", "mes_contr_3PT_"+m+"_tw1", "B1P5", Sort_light_confs);
    PT3_B2P5_tw2.Read("../heavy_radiative_tw2", "mes_contr_3PT_"+m+"_tw2", "B1P5", Sort_light_confs);

    B2B2_sm_loc.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm_loc.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm_loc.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+mm, "B3B3", Sort_light_confs);


    //#################################################################################################à
    //############################ mc1 and mc2 ##########################################################


    
    data_t P5P5_sm_loc_mc1, A0P5_sm_loc_mc1;
    data_t P5P5_sm_mc1, P5P5_sm_REST_mc1;
    data_t B1B1_sm_mc1, B3B3_sm_mc1, B2B2_sm_mc1, B2B2_MOT_sm_mc1, V1V1_sm_mc1;
    data_t PT3_B2P5_tw1_mc1;
    data_t B1B1_sm_loc_mc1, B2B2_sm_loc_mc1, B3B3_sm_loc_mc1;
    data_t V1V1_sm_REST_mc1;
    data_t V2V2_sm_REST_mc1;
    data_t V3V3_sm_REST_mc1;

    P5P5_sm_loc_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+mm, "P5P5", Sort_light_confs);
    A0P5_sm_loc_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+mm, "A0P5", Sort_light_confs);
    P5P5_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_MOT_"+m, "P5P5", Sort_light_confs);
    P5P5_sm_REST_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "P5P5", Sort_light_confs);
    V1V1_sm_REST_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    V2V2_sm_REST_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "V2V2", Sort_light_confs);
    V3V3_sm_REST_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "V3V3", Sort_light_confs);
    B2B2_MOT_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_MOT_"+m, "B1B1", Sort_light_confs);
    B2B2_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "B3B3", Sort_light_confs);
    V1V1_sm_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_3PT_"+m+"_tw1", "B1P5", Sort_light_confs);
   
    B2B2_sm_loc_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm_loc_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm_loc_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+mm, "B3B3", Sort_light_confs);



       
    data_t P5P5_sm_loc_mc2, A0P5_sm_loc_mc2;
    data_t P5P5_sm_mc2, P5P5_sm_REST_mc2;
    data_t B1B1_sm_mc2, B3B3_sm_mc2, B2B2_sm_mc2, B2B2_MOT_sm_mc2, V1V1_sm_mc2;
    data_t PT3_B2P5_tw1_mc2;
    data_t B1B1_sm_loc_mc2, B2B2_sm_loc_mc2, B3B3_sm_loc_mc2;
    data_t V1V1_sm_REST_mc2;
    data_t V2V2_sm_REST_mc2;
    data_t V3V3_sm_REST_mc2;

    P5P5_sm_loc_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+mm, "P5P5", Sort_light_confs);
    A0P5_sm_loc_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+mm, "A0P5", Sort_light_confs);
    P5P5_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_MOT_"+m, "P5P5", Sort_light_confs);
    P5P5_sm_REST_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "P5P5", Sort_light_confs);
    V1V1_sm_REST_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    V2V2_sm_REST_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "V2V2", Sort_light_confs);
    V3V3_sm_REST_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "V3V3", Sort_light_confs);
    B2B2_MOT_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_MOT_"+m, "B1B1", Sort_light_confs);
    B2B2_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "B3B3", Sort_light_confs);
    V1V1_sm_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+mm, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_3PT_"+m+"_tw1", "B1P5", Sort_light_confs);
  

    B2B2_sm_loc_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+mm, "B1B1", Sort_light_confs);
    B1B1_sm_loc_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+mm, "B2B2", Sort_light_confs);
    B3B3_sm_loc_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+mm, "B3B3", Sort_light_confs);




    //##################################################################################################
    

    data_t P5P5_sm_loc_M_minus, A0P5_sm_loc_M_minus;
    data_t P5P5_sm_M_minus, P5P5_sm_REST_M_minus;
    
    data_t B1B1_sm_M_minus, B3B3_sm_M_minus, B2B2_sm_M_minus, B2B2_MOT_sm_M_minus, V1V1_sm_M_minus;
    data_t B1B1_sm_loc_M_minus, B2B2_sm_loc_M_minus, B3B3_sm_loc_M_minus;
    data_t PT3_B2P5_tw1_M_minus;
    
    P5P5_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "P5P5", Sort_light_confs);
    A0P5_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "A0P5", Sort_light_confs);
    P5P5_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_MOT_"+m_minus, "P5P5", Sort_light_confs);
    P5P5_sm_REST_M_minus.Read("../heavy_radiative", "mes_contr_2PT_"+m_minus, "P5P5", Sort_light_confs);
    B2B2_MOT_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_MOT_"+m_minus, "B1B1", Sort_light_confs);
    B2B2_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_"+m_minus, "B3B3", Sort_light_confs);
    V1V1_sm_M_minus.Read("../heavy_radiative", "mes_contr_2PT_"+m_minus, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1_M_minus.Read("../heavy_radiative", "mes_contr_3PT_"+m_minus+"_tw1", "B1P5", Sort_light_confs);

    B2B2_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B3B3", Sort_light_confs);


    //#####################################################################################################
    //############################# mc1 and mc2 ###########################################################

    data_t P5P5_sm_loc_M_minus_mc1, A0P5_sm_loc_M_minus_mc1;
    data_t P5P5_sm_M_minus_mc1, P5P5_sm_REST_M_minus_mc1;
    
    data_t B1B1_sm_M_minus_mc1, B3B3_sm_M_minus_mc1, B2B2_sm_M_minus_mc1, B2B2_MOT_sm_M_minus_mc1, V1V1_sm_M_minus_mc1;
    data_t B1B1_sm_loc_M_minus_mc1, B2B2_sm_loc_M_minus_mc1, B3B3_sm_loc_M_minus_mc1;
    data_t PT3_B2P5_tw1_M_minus_mc1;
    
    
    P5P5_sm_loc_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+m_minus, "P5P5", Sort_light_confs);
    A0P5_sm_loc_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+m_minus, "A0P5", Sort_light_confs);
    P5P5_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_MOT_"+m_minus, "P5P5", Sort_light_confs);
    P5P5_sm_REST_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+m_minus, "P5P5", Sort_light_confs);
    B2B2_MOT_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_MOT_"+m_minus, "B1B1", Sort_light_confs);
    B2B2_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+m_minus, "B3B3", Sort_light_confs);
    V1V1_sm_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_"+m_minus, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_3PT_"+m_minus+"_tw1", "B1P5", Sort_light_confs);
   

    B2B2_sm_loc_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_loc_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_loc_M_minus_mc1.Read("../heavy_radiative_etac/mc1", "mes_contr_2PT_LOC_"+m_minus, "B3B3", Sort_light_confs);



    
    data_t P5P5_sm_loc_M_minus_mc2, A0P5_sm_loc_M_minus_mc2;
    data_t P5P5_sm_M_minus_mc2, P5P5_sm_REST_M_minus_mc2;
    
    data_t B1B1_sm_M_minus_mc2, B3B3_sm_M_minus_mc2, B2B2_sm_M_minus_mc2, B2B2_MOT_sm_M_minus_mc2, V1V1_sm_M_minus_mc2;
    data_t B1B1_sm_loc_M_minus_mc2, B2B2_sm_loc_M_minus_mc2, B3B3_sm_loc_M_minus_mc2;
    data_t PT3_B2P5_tw1_M_minus_mc2;
    
    P5P5_sm_loc_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+m_minus, "P5P5", Sort_light_confs);
    A0P5_sm_loc_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+m_minus, "A0P5", Sort_light_confs);
    P5P5_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_MOT_"+m_minus, "P5P5", Sort_light_confs);
    P5P5_sm_REST_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+m_minus, "P5P5", Sort_light_confs);
    B2B2_MOT_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_MOT_"+m_minus, "B1B1", Sort_light_confs);
    B2B2_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+m_minus, "B3B3", Sort_light_confs);
    V1V1_sm_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_"+m_minus, "V1V1", Sort_light_confs);
    PT3_B2P5_tw1_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_3PT_"+m_minus+"_tw1", "B1P5", Sort_light_confs);

    B2B2_sm_loc_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_loc_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_loc_M_minus_mc2.Read("../heavy_radiative_etac/mc2", "mes_contr_2PT_LOC_"+m_minus, "B3B3", Sort_light_confs);








    //#####################################################################################################

    

    
    int Nens=P5P5_sm_loc.size;

    for(int iens=0;iens<Nens;iens++) {

      string Ens =P5P5_sm_loc.Tag[iens];

      bool compute_2nd_tw= (Ens=="cB211b.072.64");

      boost::filesystem::create_directory("../data/heavy_radiative/"+m+"/"+Ens);

      cout<<"########### ANALYZING ENSEMBLE "<<Ens<<" #############"<<endl;


      //########################################## ENSEMBLE INFO ############################################
      int tw1, tw2;
      double amc, amc2, amh;
      distr_t a_distr(UseJack);
      distr_t ZV_had(UseJack);
      if(Ens=="cB211b.072.64") {
	tw1=20; tw2=30;
	a_distr=a_B;
	amc= 0.231567;
	amc2=0.239;
	ZV_had= ZV_B;
      }
      else if(Ens=="cA211a.12.48") {
	tw1=19;
	a_distr=a_A;
	amc= 0.264;
	amc2=0.2727;
	ZV_had= ZV_A;
      }
      else if(Ens=="cD211a.054.96") {
	tw1=30; //tw2=30;
	a_distr=a_D;
	amc= 0.164898;
	amc2= 0.167;
	ZV_had= ZV_D;
      }
      else if(Ens=="cC211a.06.80") {
	tw1=25; //tw2=30;
	a_distr=a_C;
	amc= 0.1984;
	amc2=0.203;
	ZV_had=ZV_C;
      }
      else if(Ens=="cE211a.044.112") {
	tw1=35; //tw2=30;
	a_distr=a_E;
	amc= 0.141255;
	amc2= 0.1423;
	ZV_had=ZV_E;
      }
      else crash("Error");
      amh= amc*pow(lambda, stoi(mm.substr(1,1))-1);
      double amh2 = amc2*pow(lambda, stoi(mm.substr(1,1))-1);
      //#################################################################################################

      
      //set up statistical analysis
      CorrAnalysis Corr(UseJack,Njacks,100);
      Corr.Nt= P5P5_sm_loc.nrows[iens];
      Corr.Reflection_sign=1;
      Corr.Perform_Nt_t_average=1;



      
      //########################## LOAD 2PT and 3PT FUNCS ################################################
      distr_t_list P5P5_sm_loc_distr= Corr.corr_t(P5P5_sm_loc.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/P5P5_sm_loc");
      distr_t_list P5P5_sm_loc_distr_M_minus= Corr.corr_t(P5P5_sm_loc_M_minus.col(0)[iens], "");
      Corr.Reflection_sign=-1;
      distr_t_list A0P5_sm_loc_distr= Corr.corr_t(A0P5_sm_loc.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/A0P5_sm_loc");
      distr_t_list A0P5_sm_loc_distr_M_minus= Corr.corr_t(A0P5_sm_loc_M_minus.col(0)[iens], "");
      Corr.Reflection_sign=1;
      distr_t_list P5P5_sm_distr= Corr.corr_t(P5P5_sm.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/P5P5_sm");
      distr_t_list P5P5_sm_distr_M_minus= Corr.corr_t(P5P5_sm_M_minus.col(0)[iens], "");
      distr_t_list P5P5_sm_REST_distr= Corr.corr_t(P5P5_sm_REST.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/P5P5_sm_REST");
      distr_t_list VKVK_sm_REST_distr= Corr.corr_t(summ_master(V1V1_sm_REST.col(0)[iens], V2V2_sm_REST.col(0)[iens], V3V3_sm_REST.col(0)[iens]), "../data/heavy_radiative/"+m+"/"+Ens+"/VKVK_sm_REST");
      distr_t_list P5P5_sm_REST_M_minus_distr= Corr.corr_t(P5P5_sm_REST_M_minus.col(0)[iens], "");
      distr_t_list B2B2_sm_distr= Corr.corr_t(B2B2_sm.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/B2B2_sm");
      distr_t_list BKBK_sm_distr= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm.col(0)[iens], B1B1_sm.col(0)[iens], B3B3_sm.col(0)[iens]), "../data/heavy_radiative/"+m+"/"+Ens+"/BKBK_sm");
      distr_t_list BKBK_sm_distr_M_minus= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_M_minus.col(0)[iens], B1B1_sm_M_minus.col(0)[iens], B3B3_sm_M_minus.col(0)[iens]), "");
      distr_t_list BKBK_sm_loc_distr= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc.col(0)[iens], B1B1_sm_loc.col(0)[iens], B3B3_sm_loc.col(0)[iens]), "../data/heavy_radiative/"+m+"/"+Ens+"/BKBK_sm_loc");
      distr_t_list BKBK_sm_loc_distr_M_minus= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc_M_minus.col(0)[iens], B1B1_sm_loc_M_minus.col(0)[iens], B3B3_sm_loc_M_minus.col(0)[iens]), "");
      distr_t_list V1V1_sm_distr= Corr.corr_t(V1V1_sm.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/V1V1_sm");
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_tw1=Corr.corr_t(PT3_B2P5_tw1.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_tw1");
      distr_t_list PT3_tw1_IM= Corr.corr_t(PT3_B2P5_tw1.col(1)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_tw1_IM");
      distr_t_list PT3_tw1_M_minus=Corr.corr_t(PT3_B2P5_tw1_M_minus.col(0)[iens], "");
      distr_t_list PT3_tw2(UseJack);
      if(compute_2nd_tw) {
	PT3_tw2=Corr.corr_t(PT3_B2P5_tw2.col(0)[0], "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_tw2");
      }
      Corr.Perform_Nt_t_average=1;



      //##################### mc1 and mc2 ##################################################################
      
      distr_t_list P5P5_sm_loc_distr_mc1= Corr.corr_t(P5P5_sm_loc_mc1.col(0)[iens], "");
      distr_t_list P5P5_sm_loc_distr_M_minus_mc1= Corr.corr_t(P5P5_sm_loc_M_minus_mc1.col(0)[iens], "");
      Corr.Reflection_sign=-1;
      distr_t_list A0P5_sm_loc_distr_mc1= Corr.corr_t(A0P5_sm_loc_mc1.col(0)[iens], "");
      distr_t_list A0P5_sm_loc_distr_M_minus_mc1= Corr.corr_t(A0P5_sm_loc_M_minus_mc1.col(0)[iens], "");
      Corr.Reflection_sign=1;
      distr_t_list P5P5_sm_distr_mc1= Corr.corr_t(P5P5_sm_mc1.col(0)[iens], "");
      distr_t_list P5P5_sm_distr_M_minus_mc1= Corr.corr_t(P5P5_sm_M_minus_mc1.col(0)[iens], "");
      distr_t_list P5P5_sm_REST_distr_mc1= Corr.corr_t(P5P5_sm_REST_mc1.col(0)[iens],  "");
      distr_t_list VKVK_sm_REST_distr_mc1= Corr.corr_t(summ_master(V1V1_sm_REST_mc1.col(0)[iens], V2V2_sm_REST_mc1.col(0)[iens], V3V3_sm_REST_mc1.col(0)[iens]), "");
      distr_t_list P5P5_sm_REST_M_minus_distr_mc1= Corr.corr_t(P5P5_sm_REST_M_minus_mc1.col(0)[iens], "");
      distr_t_list B2B2_sm_distr_mc1= Corr.corr_t(B2B2_sm_mc1.col(0)[iens],  "");
      distr_t_list BKBK_sm_distr_mc1= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_mc1.col(0)[iens], B1B1_sm_mc1.col(0)[iens], B3B3_sm_mc1.col(0)[iens]),  "");
      distr_t_list BKBK_sm_distr_M_minus_mc1= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_M_minus_mc1.col(0)[iens], B1B1_sm_M_minus_mc1.col(0)[iens], B3B3_sm_M_minus_mc1.col(0)[iens]), "");
      distr_t_list BKBK_sm_loc_distr_mc1= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc_mc1.col(0)[iens], B1B1_sm_loc_mc1.col(0)[iens], B3B3_sm_loc_mc1.col(0)[iens]),  "");
      distr_t_list BKBK_sm_loc_distr_M_minus_mc1= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc_M_minus_mc1.col(0)[iens], B1B1_sm_loc_M_minus_mc1.col(0)[iens], B3B3_sm_loc_M_minus_mc1.col(0)[iens]), "");
      distr_t_list V1V1_sm_distr_mc1= Corr.corr_t(V1V1_sm_mc1.col(0)[iens], "");
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_tw1_mc1=Corr.corr_t(PT3_B2P5_tw1_mc1.col(0)[iens],  "");
      distr_t_list PT3_tw1_M_minus_mc1= Corr.corr_t(PT3_B2P5_tw1_M_minus_mc1.col(0)[iens], "");


      distr_t_list P5P5_sm_loc_distr_mc2= Corr.corr_t(P5P5_sm_loc_mc2.col(0)[iens], "");
      distr_t_list P5P5_sm_loc_distr_M_minus_mc2= Corr.corr_t(P5P5_sm_loc_M_minus_mc2.col(0)[iens], "");
      Corr.Reflection_sign=-1;
      distr_t_list A0P5_sm_loc_distr_mc2= Corr.corr_t(A0P5_sm_loc_mc2.col(0)[iens], "");
      distr_t_list A0P5_sm_loc_distr_M_minus_mc2= Corr.corr_t(A0P5_sm_loc_M_minus_mc2.col(0)[iens], "");
      Corr.Reflection_sign=1;
      distr_t_list P5P5_sm_distr_mc2= Corr.corr_t(P5P5_sm_mc2.col(0)[iens], "");
      distr_t_list P5P5_sm_distr_M_minus_mc2= Corr.corr_t(P5P5_sm_M_minus_mc2.col(0)[iens], "");
      distr_t_list P5P5_sm_REST_distr_mc2= Corr.corr_t(P5P5_sm_REST_mc2.col(0)[iens],  "");
      distr_t_list VKVK_sm_REST_distr_mc2= Corr.corr_t(summ_master(V1V1_sm_REST_mc2.col(0)[iens], V2V2_sm_REST_mc2.col(0)[iens], V3V3_sm_REST_mc2.col(0)[iens]), "");
      distr_t_list P5P5_sm_REST_M_minus_distr_mc2= Corr.corr_t(P5P5_sm_REST_M_minus_mc2.col(0)[iens], "");
      distr_t_list B2B2_sm_distr_mc2= Corr.corr_t(B2B2_sm_mc2.col(0)[iens],  "");
      distr_t_list BKBK_sm_distr_mc2= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_mc2.col(0)[iens], B1B1_sm_mc2.col(0)[iens], B3B3_sm_mc2.col(0)[iens]),  "");
      distr_t_list BKBK_sm_distr_M_minus_mc2= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_M_minus_mc2.col(0)[iens], B1B1_sm_M_minus_mc2.col(0)[iens], B3B3_sm_M_minus_mc2.col(0)[iens]), "");
      distr_t_list BKBK_sm_loc_distr_mc2= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc_mc2.col(0)[iens], B1B1_sm_loc_mc2.col(0)[iens], B3B3_sm_loc_mc2.col(0)[iens]),  "");
      distr_t_list BKBK_sm_loc_distr_M_minus_mc2= (1.0/3.0)*Corr.corr_t(summ_master(B2B2_sm_loc_M_minus_mc2.col(0)[iens], B1B1_sm_loc_M_minus_mc2.col(0)[iens], B3B3_sm_loc_M_minus_mc2.col(0)[iens]), "");
      distr_t_list V1V1_sm_distr_mc2= Corr.corr_t(V1V1_sm_mc2.col(0)[iens], "");
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_tw1_mc2=Corr.corr_t(PT3_B2P5_tw1_mc2.col(0)[iens],  "");
      distr_t_list PT3_tw1_M_minus_mc2= Corr.corr_t(PT3_B2P5_tw1_M_minus_mc2.col(0)[iens], "");




      
      //######################################################################################################



      


      
      //#################################################################################################
     


      //####################    DETERMINE   EFFECTIVE MASSES ############################################
      distr_t_list H_mass_distr= Corr.effective_mass_t(B2B2_sm_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass");
      distr_t_list H_mass_averaged_distr= Corr.effective_mass_t(BKBK_sm_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_averaged");
      distr_t_list H_mass_averaged_M_minus_distr= Corr.effective_mass_t(BKBK_sm_distr_M_minus, "");
      distr_t_list H_mass_sm_loc_averaged_distr= Corr.effective_mass_t(BKBK_sm_loc_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_sm_loc_averaged");
      distr_t_list H_mass_sm_loc_averaged_M_minus_distr= Corr.effective_mass_t(BKBK_sm_loc_distr_M_minus, "");
      distr_t_list H_mass_MOT_distr = Corr.effective_mass_t( B2B2_MOT_sm.col(0)[iens], "../data/heavy_radiative/"+m+"/"+Ens+"/H_MOT_mass");
      //distr_t_list JPsi_mass_distr= Corr.effective_mass_t(VKVK_sm_REST_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/Jpsi_mass");
      distr_t_list eta_mass_distr= Corr.effective_mass_t(P5P5_sm_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass");
      distr_t_list eta_mass_distr_M_minus= Corr.effective_mass_t(P5P5_sm_distr_M_minus, "");
      distr_t_list eta_mass_rest_distr_M_minus = Corr.effective_mass_t(P5P5_sm_REST_M_minus_distr, "");
      distr_t_list eta_mass_rest_distr= Corr.effective_mass_t(P5P5_sm_REST_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST");
      distr_t_list Jpsi_mass_rest_distr= Corr.effective_mass_t(VKVK_sm_REST_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/Jpsi_mass_REST");
      distr_t_list eta_mass_rest_loc_distr = Corr.effective_mass_t(P5P5_sm_loc_distr, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST_loc");
      Print_To_File({} , { (eta_mass_rest_distr/a_distr).ave(), (eta_mass_rest_distr/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST_pu", "", "");
      Print_To_File({} , { (H_mass_averaged_distr/a_distr).ave(), (H_mass_averaged_distr/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_averaged_pu", "", "");
      distr_t_list DH_distr= H_mass_averaged_distr - H_mass_averaged_M_minus_distr;
      distr_t_list DH_sm_loc_distr = H_mass_sm_loc_averaged_distr - H_mass_sm_loc_averaged_M_minus_distr;
      Print_To_File({}, { DH_distr.ave(), DH_distr.err(), DH_sm_loc_distr.ave(), DH_sm_loc_distr.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/DH_mass", "", "");
      Print_To_File({}, { (DH_distr/a_distr).ave(), (DH_distr/a_distr).err(), (DH_sm_loc_distr/a_distr).ave(), (DH_sm_loc_distr/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/DH_mass_pu", "", "");


      //############################## mc1 and mc2 ######################################################


      distr_t_list H_mass_averaged_distr_mc1= Corr.effective_mass_t(BKBK_sm_distr_mc1, "");
      distr_t_list H_mass_averaged_M_minus_distr_mc1= Corr.effective_mass_t(BKBK_sm_distr_M_minus_mc1, "");
      distr_t_list H_mass_sm_loc_averaged_distr_mc1= Corr.effective_mass_t(BKBK_sm_loc_distr_mc1, "");
      distr_t_list H_mass_sm_loc_averaged_M_minus_distr_mc1= Corr.effective_mass_t(BKBK_sm_loc_distr_M_minus_mc1, "");
      
      distr_t_list eta_mass_distr_mc1= Corr.effective_mass_t(P5P5_sm_distr_mc1, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass");
      distr_t_list eta_mass_distr_M_minus_mc1= Corr.effective_mass_t(P5P5_sm_distr_M_minus_mc1, "");
      distr_t_list eta_mass_rest_distr_M_minus_mc1 = Corr.effective_mass_t(P5P5_sm_REST_M_minus_distr_mc1, "");
      distr_t_list eta_mass_rest_distr_mc1= Corr.effective_mass_t(P5P5_sm_REST_distr_mc1, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST");
      Print_To_File({} , { (eta_mass_rest_distr_mc1/a_distr).ave(), (eta_mass_rest_distr_mc1/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST_pu_mc1", "", "");
      Print_To_File({} , { (H_mass_averaged_distr_mc1/a_distr).ave(), (H_mass_averaged_distr_mc1/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_averaged_pu_mc1", "", "");

      distr_t_list H_mass_averaged_distr_mc2= Corr.effective_mass_t(BKBK_sm_distr_mc2, "");
      distr_t_list H_mass_averaged_M_minus_distr_mc2= Corr.effective_mass_t(BKBK_sm_distr_M_minus_mc2, "");
      distr_t_list H_mass_sm_loc_averaged_distr_mc2= Corr.effective_mass_t(BKBK_sm_loc_distr_mc2, "");
      distr_t_list H_mass_sm_loc_averaged_M_minus_distr_mc2= Corr.effective_mass_t(BKBK_sm_loc_distr_M_minus_mc2, "");
      
      distr_t_list eta_mass_distr_mc2= Corr.effective_mass_t(P5P5_sm_distr_mc2, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass");
      distr_t_list eta_mass_distr_M_minus_mc2= Corr.effective_mass_t(P5P5_sm_distr_M_minus_mc2, "");
      distr_t_list eta_mass_rest_distr_M_minus_mc2 = Corr.effective_mass_t(P5P5_sm_REST_M_minus_distr_mc2, "");
      distr_t_list eta_mass_rest_distr_mc2= Corr.effective_mass_t(P5P5_sm_REST_distr_mc2, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST");
      Print_To_File({} , { (eta_mass_rest_distr_mc2/a_distr).ave(), (eta_mass_rest_distr_mc2/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_REST_pu_mc2", "", "");
      Print_To_File({} , { (H_mass_averaged_distr_mc2/a_distr).ave(), (H_mass_averaged_distr_mc2/a_distr).err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_averaged_pu_mc2", "", "");
      

      Print_To_File({}, { (eta_mass_distr_mc2 -eta_mass_distr_mc1).ave() , (eta_mass_distr_mc2-eta_mass_distr_mc1).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/deta_mc12", "", "");
   
      

      distr_t_list deta_mass_distr = eta_mass_distr_mc2 -eta_mass_distr_mc1;
      distr_t_list dh_mass_distr = H_mass_averaged_distr_mc2 -H_mass_averaged_distr_mc1;
      distr_t_list deta_mass_rest_distr = eta_mass_rest_distr_mc2 -eta_mass_rest_distr_mc1;

      distr_t_list deta_mass_minus_distr =  eta_mass_distr_M_minus_mc2 -eta_mass_distr_M_minus_mc1;
      distr_t_list dh_mass_minus_distr = H_mass_averaged_M_minus_distr_mc2 -H_mass_averaged_M_minus_distr_mc1;
      distr_t_list deta_mass_rest_minus_distr= eta_mass_rest_distr_M_minus_mc2 - eta_mass_rest_distr_M_minus_mc1;
      

      distr_t_list dh_mass_rel_distr = H_mass_averaged_distr_mc1/H_mass_averaged_distr_mc2;

     
      distr_t_list dh_mass_minus_rel_distr = H_mass_averaged_M_minus_distr_mc1/H_mass_averaged_M_minus_distr_mc2;

      distr_t_list dh_mm_mass_rel_distr = dh_mass_rel_distr/dh_mass_minus_rel_distr;

      Print_To_File({}, { (dh_mm_mass_rel_distr).ave() , (dh_mm_mass_rel_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/dh_mm_rel_mc12", "", "");

      Print_To_File({}, { dh_mass_rel_distr.ave(), dh_mass_rel_distr.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/dh_rel_mc12", "","");

      Print_To_File({}, { dh_mass_distr.ave(), dh_mass_distr.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/dh_mc12", "","");
          
      //#################################################################################################

      
      //determine DH mass
      Get_plateaux_int_2pt_dH_loc(Ens,m,Corr);
      distr_t DH_sm_loc = Corr.Fit_distr(DH_sm_loc_distr);
      Get_plateaux_int_2pt_dH_sm(Ens,m,Corr);
      distr_t DH_sm = Corr.Fit_distr(DH_distr);
      //Corr.Tmin=12; Corr.Tmax=15;
      Get_plateaux_dmc(Ens, m, Corr);
      distr_t DH_dmc= Corr.Fit_distr(dh_mass_distr);
      //DH_dmc= DH_dmc.ave() + (DH_dmc - DH_dmc.ave())*2;
      Get_plateaux_dmc(Ens, (m=="M1" || m=="Hc")?m:M[im-1], Corr);
      distr_t DH_minus_dmc = Corr.Fit_distr(dh_mass_minus_distr);
      //DH_minus_dmc = DH_minus_dmc.ave() + (DH_minus_dmc - DH_minus_dmc.ave())*2.0;
      Get_plateaux_dmc(Ens,m, Corr);
      distr_t DH_mm_mass_rel = Corr.Fit_distr(dh_mm_mass_rel_distr);
      distr_t DH_mass_rel = Corr.Fit_distr(dh_mass_rel_distr);
      distr_t DH_mass_minus_rel = Corr.Fit_distr(dh_mass_minus_rel_distr);
     
      double w1= 0.5; // 1/pow(DH_sm.err(),2);
      double w2= 0.5; // 1/pow(DH_sm_loc.err(),2);
      w1= w1/(w1+w2); w2 = 1- w1;
      distr_t DH_mass = DH_sm;
      double syst = fabs(DH_sm.ave() - DH_sm_loc.ave());  //w1*sqrt(pow(DH_sm.ave() - DH_mass.ave(),2)) + w2*sqrt(pow(DH_sm_loc.ave() - DH_mass.ave(),2));
      //syst = sqrt(syst);
      DH_mass= DH_mass.ave() + (DH_mass-DH_mass.ave())*sqrt( 1 + pow( syst/DH_mass.err(),2));
      if(Ens=="cA211a.14.48") DH_mass = DH_mass.ave() - DH_mass.err() +  (DH_mass-DH_mass.ave())*sqrt( 1 + pow( 2.0*syst/DH_mass.err(),2));
      distr_t DH_mass_mc2= DH_mass + DH_dmc - DH_minus_dmc;

      
      Get_plateaux_int_2pt_dH_sm(Ens,m,Corr);
     

      //Print DH_mass
      distr_t_list DH_mass_fit = 0.0*Get_id_jack_distr_list(Corr.Nt, Njacks) ;
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	DH_mass_fit.distr_list[t] = DH_mass;
      }

      Print_To_File({}, {DH_mass_fit.ave(), DH_mass_fit.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/DH_mass_fit", "", "");
      Print_To_File({}, {(DH_mass_fit/a_distr).ave(), (DH_mass_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/DH_mass_fit_pu", "", "");
      
      //determine amputating factor for ratio
      distr_t_list FACT_ratio = EXPT_D(-1.0*DH_mass, Corr.Nt);
      distr_t_list FACT_ratio_II = EXPT_DL( -1.0*DH_distr);
      distr_t_list PT2_ratio= BKBK_sm_distr/(BKBK_sm_distr_M_minus);
      distr_t_list shift_PT2_ratio=PT2_ratio;
      for(int t=0;t<Corr.Nt;t++) shift_PT2_ratio.distr_list[(t+tw1)%Corr.Nt] = PT2_ratio[t];
      distr_t_list VEV_ratio_distr= SQRT_DL( PT2_ratio/FACT_ratio_II);
      Print_To_File({}, {VEV_ratio_distr.ave(), VEV_ratio_distr.err()} , "../data/heavy_radiative/"+m+"/"+Ens+"/VEV_ratio", "", "");
      Get_plateaux_VEV_ratio(Ens, m, Corr);
      distr_t VEV_ratio= Corr.Fit_distr(VEV_ratio_distr);

      //determine amputating factor for mc1 and mc2
      distr_t_list FACT_ratio_II_dmc = EXPT_D( -1.0*DH_dmc, Corr.Nt);
      distr_t_list PT2_ratio_dmc= BKBK_sm_distr_mc2/(BKBK_sm_distr_mc1);
      distr_t_list VEV_ratio_distr_dmc= SQRT_DL( PT2_ratio_dmc/FACT_ratio_II_dmc);
      Print_To_File({}, {VEV_ratio_distr_dmc.ave(), VEV_ratio_distr_dmc.err()} , "../data/heavy_radiative/"+m+"/"+Ens+"/VEV_ratio_dmc", "", "");
      Get_plateaux_VEV_ratio(Ens, m, Corr);
      distr_t VEV_ratio_dmc= Corr.Fit_distr(VEV_ratio_distr_dmc);

      distr_t_list FACT_ratio_II_minus_dmc = EXPT_D( -1.0*DH_minus_dmc, Corr.Nt);
      distr_t_list PT2_ratio_minus_dmc= BKBK_sm_distr_M_minus_mc2/(BKBK_sm_distr_M_minus_mc1);
      distr_t_list VEV_ratio_distr_minus_dmc= SQRT_DL( PT2_ratio_minus_dmc/FACT_ratio_II_minus_dmc);
      Print_To_File({}, {VEV_ratio_distr_minus_dmc.ave(), VEV_ratio_distr_minus_dmc.err()} , "../data/heavy_radiative/"+m+"/"+Ens+"/VEV_ratio_minus_dmc", "", "");
      Get_plateaux_VEV_ratio(Ens, m, Corr);
      distr_t VEV_ratio_minus_dmc= Corr.Fit_distr(VEV_ratio_distr_minus_dmc);




      


      //determine Hb mass
      Get_plateaux_int_2pt_H_loc(Ens,Corr);
      distr_t H_mass_sm_loc = Corr.Fit_distr(H_mass_sm_loc_averaged_distr);
      Get_plateaux_int_2pt_H_sm(Ens,Corr);
      distr_t H_mass =  Corr.Fit_distr(H_mass_averaged_distr);
      distr_t H_mass_ave= 0.5*(H_mass+H_mass_sm_loc);
      H_mass= H_mass_ave.ave() + (H_mass_ave - H_mass_ave.ave())*sqrt( 1 + pow( 0.5*(H_mass.ave() - H_mass_sm_loc.ave())/H_mass_ave.err(),2));

      distr_t H_mass_mc2= H_mass + DH_dmc;
      
      

    
      //Print H_mass_fit
      distr_t_list H_mass_fit = 0.0*Get_id_jack_distr_list(Corr.Nt, Njacks) ;
      distr_t_list H_mass_fit_mc2= 0.0*Get_id_jack_distr_list(Corr.Nt, Njacks);
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	H_mass_fit.distr_list[t] = H_mass;
	H_mass_fit_mc2.distr_list[t] = H_mass_mc2;
      }

      Print_To_File({}, { (H_mass_fit/a_distr).ave(), (H_mass_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_fit", "", "");
      Print_To_File({}, { (H_mass_fit_mc2/a_distr).ave(), (H_mass_fit_mc2/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_fit_mc2", "", "");

      distr_t H_mass_minus= H_mass - DH_mass;
      
      distr_t_list MEM_ratio= VEV_ratio*SQRT_D( H_mass_minus/H_mass)*EXPT_D(-1.0*DH_mass, Corr.Nt)*EXP_D(DH_mass*tw1);
      distr_t_list VV_H = Corr.matrix_element_t(BKBK_sm_distr, "");



      distr_t_list MEM_ratio_dmc= VEV_ratio_dmc*SQRT_D( DH_mass_rel)*EXPT_D(-1.0*DH_dmc, Corr.Nt)*EXP_D(DH_dmc*tw1);
      distr_t_list MEM_ratio_minus_dmc= VEV_ratio_minus_dmc*SQRT_D( DH_mass_minus_rel)*EXPT_D(-1.0*DH_minus_dmc, Corr.Nt)*EXP_D(DH_minus_dmc*tw1);

      

      //distr_t_list MEM_ratio_bis=shift_PT2_ratio/(VEV_ratio*SQRT_D(H_mass/H_mass_minus));
      
           

      
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
      distr_t eta_mass = Corr.Fit_distr(eta_mass_distr);
      distr_t eta_mass_M_minus = Corr.Fit_distr(eta_mass_distr_M_minus);
      distr_t eta_mass_rest= Corr.Fit_distr(eta_mass_rest_distr);
      distr_t eta_mass_M_minus_rest=Corr.Fit_distr(eta_mass_rest_distr_M_minus);
      Print_To_File({}, { ((eta_mass_rest_distr-eta_mass_rest_distr_M_minus)/eta_mass_rest_distr_M_minus).ave(), ((eta_mass_rest_distr-eta_mass_rest_distr_M_minus)/eta_mass_rest_distr_M_minus).err() },  "../data/heavy_radiative/"+m+"/"+Ens+"/Deta_eff_mass_rest_pu", "", "");
      Print_To_File({}, { ((Jpsi_mass_rest_distr- eta_mass_rest)/(a_distr)).ave(), ((Jpsi_mass_rest_distr - eta_mass_rest)/(a_distr)).err() },  "../data/heavy_radiative/"+m+"/"+Ens+"/hyperfine_splitting_pu", "", "");
      //distr_t mom_k = SQRT_D( (eta_mass*eta_mass - eta_mass_rest*eta_mass_rest)/(a_distr*a_distr));

      distr_t eta_mass_dmc = Corr.Fit_distr(deta_mass_distr);
      distr_t eta_mass_M_minus_dmc = Corr.Fit_distr(deta_mass_minus_distr);
      distr_t eta_mass_rest_dmc= Corr.Fit_distr(deta_mass_rest_distr);
      distr_t eta_mass_rest_M_minus_dmc= Corr.Fit_distr(deta_mass_rest_minus_distr);

    

      vector<distr_t> MM({eta_mass_rest/a_distr, (eta_mass_rest + eta_mass_rest_dmc)/a_distr});
      vector<distr_t> MC({ amh*Get_id_jack_distr(Njacks), amh2*Get_id_jack_distr(Njacks)});

      distr_t amh_etah=  Obs_extrapolation_meson_mass(MC, MM, Mh_masses[im] ,  "../data/heavy_radiative"  , "amh_"+m+"_extrapolation_etah_"+Ens+".dat",  UseJack, "SPLINE" );

      
      //extrapolate DH
    

      distr_t rn1= eta_mass_rest/eta_mass_M_minus_rest;
      distr_t rn2 = (eta_mass_rest+eta_mass_rest_dmc)/(eta_mass_M_minus_rest);

      //the ratios rn1 and rn2 must be interpolated to "l"

      cout<<"rn1: "<<rn1.ave()<<" rn2: "<<rn2.ave()<<" l:"<<l<<endl;
    

      cout<<"amh(etah): "<<amh_etah.ave()<<" "<<amh_etah.err()<<endl;

      //extrapolate eta_mass to the physical point

      vector<distr_t> METAC({ eta_mass_rest, (eta_mass_rest + eta_mass_rest_dmc)});
      vector<distr_t> METAC_minus({ eta_mass_M_minus_rest, (eta_mass_M_minus_rest + eta_mass_rest_M_minus_dmc)});
    

      distr_t eta_mass_mc2= eta_mass+ eta_mass_dmc;
      distr_t eta_mass_M_minus_mc2= eta_mass_M_minus + eta_mass_M_minus_dmc;

      
      
      
      

      cout<<"Ens: "<<Ens<<" a: "<<(a_distr/fm_to_inv_Gev).ave()<<" fm"<<endl;
      distr_t mom_k = (2.0/a_distr)*ASIN_D(SQRT_D(SINH_D(eta_mass/2.0)*SINH_D(eta_mass/2.0) -SINH_D(eta_mass_rest/2.0)*SINH_D(eta_mass_rest/2.0) ));
      cout<<"K: "<<mom_k.ave()<<" "<<mom_k.err()<<endl;
      mom_k = eta_mass/a_distr;
      distr_t k_cont= 0.488*a_distr;
      //distr_t C = (POW_D(SINH_D(eta_mass/2.0),2) - POW_D( SIN_D( k_cont/2),2))/POW_D( SINH_D(eta_mass_rest/2.0),2);
      distr_t C = (POW_D(eta_mass,2) - POW_D(k_cont,2))/POW_D(eta_mass_rest,2);
      cout<<"C: "<<C.ave()<<" "<<C.err()<<endl;


      //Print eta_mass_rest_fit
      distr_t_list eta_mass_rest_fit = 0.0*Get_id_jack_distr_list(Corr.Nt, Njacks) ;
      distr_t_list Deta_mass_rest_fit = 0.0*Get_id_jack_distr_list(Corr.Nt,Njacks);
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	eta_mass_rest_fit.distr_list[t] = eta_mass_rest;
	Deta_mass_rest_fit.distr_list[t] = eta_mass_rest - eta_mass_M_minus_rest;
      }

      Print_To_File({}, { (eta_mass_rest_fit/a_distr).ave(), (eta_mass_rest_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_rest_fit", "", "");
      Print_To_File({}, { (Deta_mass_rest_fit/a_distr).ave(), (Deta_mass_rest_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/Deta_mass_rest_fit", "", "");

      
      distr_t_list VV_ETA= Corr.matrix_element_t(P5P5_sm_distr, "");
      distr_t ETA_MEM_tw1 = EXP_D(-1.0*eta_mass*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr, ""))/(2.0*eta_mass);
      distr_t ETA_MEM_tw1_M_minus = EXP_D(-1.0*eta_mass_M_minus*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr_M_minus, ""))/(2.0*eta_mass_M_minus);


      distr_t ETA_MEM_tw1_mc2_ov_mc1= EXP_D(-1.0*(eta_mass_dmc)*tw1)*Corr.Fit_distr( Corr.matrix_element_t(P5P5_sm_distr_mc2,"")/Corr.matrix_element_t(P5P5_sm_distr_mc1,""))/((1.0+eta_mass_dmc/eta_mass));
      distr_t ETA_MEM_M_minus_tw1_mc2_ov_mc1= EXP_D(-1.0*(eta_mass_M_minus_dmc)*tw1)*Corr.Fit_distr( Corr.matrix_element_t(P5P5_sm_distr_M_minus_mc2,"")/Corr.matrix_element_t(P5P5_sm_distr_M_minus_mc1,""));

      

      
      
      distr_t ETA_MEM_tw2(UseJack);
      if(compute_2nd_tw) {
	ETA_MEM_tw2 = EXP_D(-1.0*eta_mass*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr, ""))/(2.0*eta_mass);
      }
      Print_To_File({}, { VV_H.ave(), VV_H.err(), VV_ETA.ave(), VV_ETA.err()} ,  "../data/heavy_radiative/"+m+"/"+Ens+"/VV", "", "");



      //#################### J/Psi - eta_c SPLITTING ##########################

      distr_t_list hyperfine_distr = Jpsi_mass_rest_distr -eta_mass_rest;

      if(Ens=="cB211b.072.64") {Corr.Tmin=26; Corr.Tmax=35;}
      else if(Ens=="cA211a.12.48") {Corr.Tmin=23; Corr.Tmax=30;}
      else if(Ens=="cC211a.06.80") { Corr.Tmin = 37; Corr.Tmax= 43; }
      else if(Ens=="cD211a.054.96") { Corr.Tmin = 50  ; Corr.Tmax= 70  ; }
      else if(Ens=="cE211a.044.112") { Corr.Tmin = 60  ; Corr.Tmax= 80 ; }

      distr_t hyperfine_split = Corr.Fit_distr(hyperfine_distr/a_distr);


      
      
   

      
      
      
      
      
      
      



      //#######################################################################
      

      //###################################     DETERMINE Z FACTORS  ###################################
      distr_t_list ZV_distr =  2*amh*P5P5_sm_loc_distr/(distr_t_list::derivative(A0P5_sm_loc_distr, 0));

      distr_t_list ZV_distr_mc1 =  2*amh*P5P5_sm_loc_distr_mc1/(distr_t_list::derivative(A0P5_sm_loc_distr_mc1, 0));
      distr_t_list ZV_distr_mc2 =  2*amh2*P5P5_sm_loc_distr_mc2/(distr_t_list::derivative(A0P5_sm_loc_distr_mc2, 0));

      distr_t_list dZ_distr_dmc= ZV_distr_mc2/ZV_distr_mc1;

      
      
      distr_t_list ZV_distr_M_minus(UseJack);
      distr_t_list dZ_distr_M_minus_dmc(UseJack);
      
      if(mm != "M1") {
	ZV_distr_M_minus= 2*(amh/lambda)*P5P5_sm_loc_distr_M_minus/(distr_t_list::derivative(A0P5_sm_loc_distr_M_minus, 0));
	distr_t_list ZV_distr_M_minus_mc1= 2*(amh/lambda)*P5P5_sm_loc_distr_M_minus_mc1/(distr_t_list::derivative(A0P5_sm_loc_distr_M_minus_mc1,0));
	distr_t_list ZV_distr_M_minus_mc2= 2*(amh2/lambda)*P5P5_sm_loc_distr_M_minus_mc2/(distr_t_list::derivative(A0P5_sm_loc_distr_M_minus_mc2,0));
	dZ_distr_M_minus_dmc = ZV_distr_M_minus_mc2/ZV_distr_M_minus_mc1;
      }
      else {  ZV_distr_M_minus= ZV_distr; dZ_distr_M_minus_dmc = ZV_distr/ZV_distr;}
      Print_To_File({}, {ZV_distr.ave(), ZV_distr.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/ZV", "", "");
      Print_To_File({}, {ZV_distr_mc1.ave(), ZV_distr_mc1.err(), ZV_distr_mc2.ave(), ZV_distr_mc2.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/ZV_dmc", "", "");
      if(Ens=="cB211b.072.64") {Corr.Tmin=30; Corr.Tmax=50; }
      else if(Ens=="cA211a.12.48") {Corr.Tmin=26; Corr.Tmax=38; }
      else if(Ens=="cD211a.054.96") {Corr.Tmin = 35; Corr.Tmax = 60;}
      else if(Ens=="cC211a.06.80") {Corr.Tmin = 40; Corr.Tmax = 55;}
      else if(Ens=="cE211a.044.112") {Corr.Tmin = 45; Corr.Tmax = 66;}
      else crash("Ens "+Ens+" not found");
      distr_t ZV= Corr.Fit_distr(ZV_distr);
      distr_t ZV_M_minus= Corr.Fit_distr(ZV_distr_M_minus);
      distr_t dZV_dmc= Corr.Fit_distr(dZ_distr_dmc);
      distr_t dZV_M_minus_dmc= Corr.Fit_distr(dZ_distr_M_minus_dmc);
      if(m=="Hc") ZV_M1.distr_list.push_back(ZV);
      cout<<"eta_fact_ratio: "<<(ETA_MEM_tw1/ETA_MEM_tw1_M_minus).ave()<<" +- "<<(ETA_MEM_tw1/ETA_MEM_tw1_M_minus).err()<<endl;
      cout<<"Z[m]/Z[mh-1]: "<<(ZV/ZV_M_minus).ave()<<" +- "<<(ZV/ZV_M_minus).err()<<endl;
      //################################################################################################

      
      //#####################       EVALUATE FORM FACTOR AND RATIOS  ###################################
      distr_t_list PT3_tw1_ratio = PT3_tw1/PT3_tw1_M_minus;
      Print_To_File({}, { PT3_tw1_ratio.ave(), PT3_tw1_ratio.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_ratio_tw1", "", "");
      distr_t_list F0_tw1 = 2.0*ZV*PT3_tw1/(ETA_MEM_tw1*H_MEM_tw1*H_mass);
      distr_t_list F0_tw1_II = 2.0*ZV_had*(ZV/ZV_M1[iens])*PT3_tw1/(ETA_MEM_tw1*H_MEM_tw1*H_mass);
      //distr_t SCALE_FACT = SQRT_D( eta_mass*H_mass/(eta_mass_M_minus*H_mass_minus));
      distr_t SCALE_FACT= SQRT_D(H_mass/H_mass_minus);
      
     
      if(TYPE=="WEAK") SCALE_FACT= 1.0*Get_id_jack_distr(Njacks);  // SCALE_FACT= 1.0/(eta_mass_rest/eta_mass_M_minus_rest);  // SCALE_FACT=1.0*Get_id_jack_distr(Njacks);
      distr_t_list ratio_corr= ((ZV/ZV_M_minus)*PT3_tw1_ratio/( (ETA_MEM_tw1/ETA_MEM_tw1_M_minus)*(H_mass/H_mass_minus)*MEM_ratio))/SCALE_FACT;
      //distr_t_list ratio_corr= (ZV/ZV_M_minus)*PT3_tw1_ratio/( (ETA_MEM_tw1/ETA_MEM_tw1_M_minus)*(H_mass_ave/H_mass_minus)*MEM_ratio_bis);
      distr_t_list ratio_corr_red(UseJack);
      distr_t_list F0_tw2(UseJack);
      distr_t_list F0_tw2_II(UseJack);
      if(compute_2nd_tw) {
	F0_tw2 = 2.0*ZV*PT3_tw2/(ETA_MEM_tw2*H_MEM_tw2*H_mass);
	F0_tw2_II = 2.0*ZV_had*(ZV/ZV_M1[iens])*PT3_tw2/(ETA_MEM_tw2*H_MEM_tw2*H_mass);
      }
      //################################################################################################



      //####################    EVALUATE FORM FACTOR AND RATIO BETWEEN mh^n and mh2^n  and mh^n-1 and mh2^n-1 ################################


      distr_t_list PT3_N_tw1_ratio = PT3_tw1_mc2/PT3_tw1_mc1;

      distr_t_list PT3_Nm1_tw1_ratio = PT3_tw1_M_minus_mc2/PT3_tw1_M_minus_mc1;

      //evaluate change in ZV

   
      distr_t_list F0_tw1_II_CORR = PT3_N_tw1_ratio/(ETA_MEM_tw1_mc2_ov_mc1*(1.0/DH_mass_rel)*MEM_ratio_dmc);

      distr_t ETA_FACT= ETA_MEM_M_minus_tw1_mc2_ov_mc1;

      distr_t_list MEM_ratio_CORR= MEM_ratio_dmc;

      
      cout<<"ETA FACT: "<<ETA_FACT.ave()<<" +- "<<ETA_FACT.err()<<endl;

      cout<<"MEM_ratio_CORR: "<<MEM_ratio_CORR.ave(20)<<" "<<MEM_ratio_CORR.err(20)<<endl;

      cout<<"PT3 N_ratio: "<<(PT3_N_tw1_ratio/PT3_Nm1_tw1_ratio).ave(20)<<" "<<(PT3_N_tw1_ratio/PT3_Nm1_tw1_ratio).err(20)<<endl;

      distr_t_list F0_tw1_M_minus_II_CORR= PT3_Nm1_tw1_ratio/(ETA_MEM_M_minus_tw1_mc2_ov_mc1*(1.0/DH_mass_minus_rel)*MEM_ratio_minus_dmc);

      distr_t_list ratio_corr_CORR= (dZV_dmc)*F0_tw1_II_CORR;

      distr_t_list ratio_corr_1_2= 1.0/ ((dZV_M_minus_dmc)*F0_tw1_M_minus_II_CORR);

      distr_t_list ratio_corr_2_2= ratio_corr_CORR*ratio_corr_1_2;

      //interpolate F0_tw1_II_CORR and ratio_corr_CORR to the reference point

      cout<<"Z_ratio: "<<(dZV_dmc).ave()<<" +- "<<(dZV_dmc).err()<<endl;
      cout<<"Z_ratio(M-minus): "<<(dZV_M_minus_dmc).ave()<<" +- "<<(dZV_M_minus_dmc).err()<<endl;

      cout<<"RATIO_CORR[t=10a]: "<<ratio_corr_CORR.ave(20)<<" +- "<<ratio_corr_CORR.err(20)<<endl;


      distr_t_list CORR_F0_interpol(UseJack), CORR_ratio_1_interpol(UseJack), CORR_ratio_2_interpol(UseJack), CORR_ratio_interpol(UseJack);

      vector<distr_t> RS({ eta_mass_rest/a_distr , (eta_mass_rest + eta_mass_rest_dmc)/a_distr});
      vector<distr_t> RS_MINUS({ eta_mass_M_minus_rest/a_distr , (eta_mass_M_minus_rest + eta_mass_rest_M_minus_dmc)/a_distr});

      vector<distr_t> RATIO_ABS_MASS({ H_mass/a_distr, (H_mass + DH_dmc)/a_distr});
      
      vector<distr_t> RATIOS_MASS({ DH_mass/a_distr, (DH_mass+ DH_dmc)/a_distr});
      
      vector<distr_t> RATIOS_2_MASS({ (DH_mass -DH_minus_dmc)/a_distr, (DH_mass + DH_dmc -DH_minus_dmc)/a_distr});

      distr_t DH1_interpol= Obs_extrapolation_meson_mass( RATIOS_MASS, RS, Mh_masses[im], "../data/heavy_radiative", "amh_"+m+"_extrapolation_DH1_"+Ens+".dat", UseJack, "SPLINE");

      distr_t DH2_interpol= Obs_extrapolation_meson_mass( RATIOS_2_MASS, RS, Mh_masses[im], "../data/heavy_radiative", "amh_"+m+"_extrapolation_DH2_"+Ens+".dat", UseJack, "SPLINE");

      vector<distr_t> RATIOS_3_MASS({ DH1_interpol, DH2_interpol});

      distr_t DH_mass_interpol = Obs_extrapolation_meson_mass( RATIOS_3_MASS, RS_MINUS,  Mh_masses[(im==0)?0:(im-1)],  "../data/heavy_radiative", "amh_"+m+"_extrapolation_DH_"+Ens+".dat", UseJack, "SPLINE");

      distr_t H_mass_interpol =  Obs_extrapolation_meson_mass( RATIO_ABS_MASS, RS, Mh_masses[im], "../data/heavy_radiative", "amh_"+m+"_extrapolation_H_"+Ens+".dat", UseJack, "SPLINE");

      for(int t=0;t<Corr.Nt;t++) {

	vector<distr_t> F0S({ 1.0*Get_id_jack_distr(Njacks), 1.0/F0_tw1_II_CORR.distr_list[t] });
	vector<distr_t> RATIOS({ 1.0*Get_id_jack_distr(Njacks), 1.0/ratio_corr_CORR.distr_list[t]});
	vector<distr_t> RATIOS_2({ 1.0/ratio_corr_1_2.distr_list[t], 1.0/ratio_corr_2_2.distr_list[t]});




      
	CORR_F0_interpol.distr_list.push_back( 1.0/Obs_extrapolation_meson_mass( F0S, RS, Metac_phys*Mh_masses[im]/Mh_masses[0], "../data/heavy_radiative", "amh_"+m+"_extrapolation_F0_t_"+to_string(t)+"_"+Ens+".dat", UseJack, "SPLINE"));

	CORR_ratio_1_interpol.distr_list.push_back( 1.0/Obs_extrapolation_meson_mass( RATIOS, RS, Mh_masses[im], "../data/heavy_radiative", "amh_"+m+"_extrapolation_ratio_t_"+to_string(t)+"_"+Ens+".dat", UseJack, "SPLINE"));

	

	CORR_ratio_2_interpol.distr_list.push_back( 1.0/Obs_extrapolation_meson_mass( RATIOS_2, RS, Mh_masses[im], "../data/heavy_radiative", "amh_"+m+"_extrapolation_ratio_2_t_"+to_string(t)+"_"+Ens+".dat", UseJack, "SPLINE"));

	vector<distr_t> RATIOS_3({ CORR_ratio_1_interpol.distr_list[t], CORR_ratio_2_interpol.distr_list[t]});

	CORR_ratio_interpol.distr_list.push_back( Obs_extrapolation_meson_mass( RATIOS_3, RS_MINUS, Mh_masses[(im==0)?0:(im-1)],  "../data/heavy_radiative", "amh_"+m+"_extrapolation_ratio_3_t_"+to_string(t)+"_"+Ens+".dat", UseJack, "SPLINE"));
	

      }
      
      



      //######################################################################################################################################



      //################################# PRINT FORM FACTORS ###########################################
      //reduce
      distr_t_list F0_tw1_RED(UseJack);
      distr_t_list F0_tw2_RED(UseJack);
      distr_t_list F0_tw1_II_RED(UseJack);
      distr_t_list F0_tw2_II_RED(UseJack);

      distr_t_list F0_tw1_RED_interpol(UseJack);
      distr_t_list ratio_corr_RED_interpol(UseJack);
      
      for(int t=tw1;t<Corr.Nt;t++)  {	F0_tw1_RED.distr_list.push_back( F0_tw1[t]);  F0_tw1_II_RED.distr_list.push_back( F0_tw1_II[t]); F0_tw1_RED_interpol.distr_list.push_back( F0_tw1_II.distr_list[t]*CORR_F0_interpol.distr_list[t]);  ratio_corr_red.distr_list.push_back( ratio_corr[t]); ratio_corr_RED_interpol.distr_list.push_back( ratio_corr.distr_list[t]*CORR_ratio_interpol.distr_list[t]);  }
      if(compute_2nd_tw) {
	for(int t=tw2;t<Corr.Nt;t++) { F0_tw2_RED.distr_list.push_back( F0_tw2[t]);  F0_tw2_II_RED.distr_list.push_back( F0_tw2_II[t]); }
      }
      Print_To_File({}, {ratio_corr_red.ave(), ratio_corr_red.err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/F_ratio", "", "");
      Print_To_File({}, {ratio_corr_RED_interpol.ave(), ratio_corr_RED_interpol.err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/F_ratio_interpol", "", "");
      
      Print_To_File({}, {F0_tw1_RED.ave(), F0_tw1_RED.err(), F0_tw1_II_RED.ave(), F0_tw1_II_RED.err() }, "../data/heavy_radiative/"+m+"/"+Ens+"/F0_tw1", "", "");
      Print_To_File({}, {F0_tw1_RED_interpol.ave(), F0_tw1_RED_interpol.err() }, "../data/heavy_radiative/"+m+"/"+Ens+"/F0_tw1_interpol", "", "");
      if(compute_2nd_tw) Print_To_File({}, {F0_tw2_RED.ave(), F0_tw2_RED.err(), F0_tw2_II_RED.ave(), F0_tw2_II_RED.err() }, "../data/heavy_radiative/"+m+"/"+Ens+"/F0_tw2", "", "");
      //################################################################################################


      //fit Fc

      Corr.Tmin= (int)(0.9/(a_distr.ave()/fm_to_inv_Gev));
      Corr.Tmax= (int)(1.3/(a_distr.ave()/fm_to_inv_Gev));

      distr_t F1=Corr.Fit_distr(F0_tw1_RED_interpol);

      
      //push back results for Fc
      if(m=="Hc") {
	Fc.distr_list.push_back(F1);
	metac.distr_list.push_back( eta_mass_rest/a_distr);
	mhc.distr_list.push_back( H_mass/a_distr);
	a_distr_list.distr_list.push_back(a_distr);
      }

      //push back result for Fm1
      if(m=="M1") {
	Fm1.distr_list.push_back(F1);
      }


      //fit ratios

      //set time intervals

      //test
      
      //H_mass = H_mass.ave() + (H_mass - H_mass.ave())*2;
      //DH_mass= DH_mass.ave() + (DH_mass - DH_mass.ave())*2;


      distr_t a_distr_new= a_distr.ave() + (a_distr - a_distr.ave())*2;
      
      distr_t to_m_etah= eta_mass_rest/a_distr ;

      if(m != "Hc") {

	m_etah_list_glb[im-1].distr_list.push_back(to_m_etah) ;
	aH_glb_list[im-1].distr_list.push_back( a_distr);
      }

      if(m=="M1") {
	mom_k_list[0].distr_list.push_back(mom_k);
	m_etah_list[0].distr_list.push_back(to_m_etah);
	aH_list[0].distr_list.push_back( a_distr);
	C_beta[0].distr_list.push_back( C);
	DH_mass_list[0].distr_list.push_back( H_mass_interpol);
	etah_mass_list[0].distr_list.push_back( eta_mass_rest/a_distr);

	Ens_list[0].push_back(Ens);
	
      }

      if(m=="M2") {
	C_beta[1].distr_list.push_back( C);
	Corr.Tmin = (int)(0.6/(a_distr.ave()/fm_to_inv_Gev));
	Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	distr_t R1= Corr.Fit_distr(ratio_corr_RED_interpol);
	R[0].distr_list.push_back(R1);
	aR[0].distr_list.push_back(a_distr);

	mom_k_list[1].distr_list.push_back(mom_k);
	m_etah_list[1].distr_list.push_back(to_m_etah);
	aH_list[1].distr_list.push_back( a_distr);

	DH_mass_list[1].distr_list.push_back( DH_mass_interpol);
	etah_mass_list[1].distr_list.push_back( eta_mass_rest/a_distr);

	Ens_list[1].push_back(Ens);
	
      }
      else if(m=="M3") {
	C_beta[2].distr_list.push_back( C);
	if(Ens != "cE211a.044.114") {
	  Corr.Tmin = (int)(0.6/(a_distr.ave()/fm_to_inv_Gev));
	  Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	}
	else {
	  Corr.Tmin = (int)(0.45/(a_distr.ave()/fm_to_inv_Gev));
	  Corr.Tmax = (int)(0.8/(a_distr.ave()/fm_to_inv_Gev));
	}
	distr_t R2 = Corr.Fit_distr(ratio_corr_RED_interpol);
	R[1].distr_list.push_back(R2);
	aR[1].distr_list.push_back(a_distr);

	mom_k_list[2].distr_list.push_back(mom_k);
	m_etah_list[2].distr_list.push_back(to_m_etah);
	aH_list[2].distr_list.push_back( a_distr);

	DH_mass_list[2].distr_list.push_back( DH_mass_interpol);
	etah_mass_list[2].distr_list.push_back( eta_mass_rest/a_distr);

	Ens_list[2].push_back(Ens);
	
      }
      else if(m=="M4") {
	C_beta[3].distr_list.push_back( C);
	if(Ens != "cA211a.14.48") {
	if(Ens == "cD211a.054.96") {
	  Corr.Tmin = (int)(0.70/(a_distr.ave()/fm_to_inv_Gev));
	  Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	}
	else if(Ens != "cE211a.044.114") {
	  Corr.Tmin = (int)(0.72/(a_distr.ave()/fm_to_inv_Gev));
	  Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	}
	else {
	  Corr.Tmin = (int)(0.49/(a_distr.ave()/fm_to_inv_Gev));
	  Corr.Tmax = (int)(0.7/(a_distr.ave()/fm_to_inv_Gev));
	}
	distr_t R3 = Corr.Fit_distr(ratio_corr_RED_interpol);
	R[2].distr_list.push_back(R3);
	aR[2].distr_list.push_back(a_distr);

	mom_k_list[3].distr_list.push_back(mom_k);
	m_etah_list[3].distr_list.push_back(to_m_etah);
	aH_list[3].distr_list.push_back( a_distr);

	DH_mass_list[3].distr_list.push_back( DH_mass_interpol);
	etah_mass_list[3].distr_list.push_back( eta_mass_rest/a_distr);

	Ens_list[3].push_back(Ens);
	}

      }
      else if(m=="M5") {
	C_beta[4].distr_list.push_back( C);
	if(Ens != "cB211b.072.64" && Ens != "cA211a.12.48") {
	  if(Ens=="cC211a.06.80") {
	    Corr.Tmin = (int)(0.5/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(0.7/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else if(Ens=="cD211a.054.96") {
	    Corr.Tmin = (int)(1.2/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.5/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else if(Ens=="cE211a.044.112") {
	    Corr.Tmin = (int)(0.9/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.2/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else crash("Ens: "+Ens+" not found") ;

	  distr_t R4 = Corr.Fit_distr(ratio_corr_RED_interpol);
	  R[3].distr_list.push_back(R4);
	  aR[3].distr_list.push_back(a_distr);

	  mom_k_list[4].distr_list.push_back(mom_k);
	  m_etah_list[4].distr_list.push_back(to_m_etah);
	  aH_list[4].distr_list.push_back( a_distr);

	  DH_mass_list[4].distr_list.push_back( DH_mass_interpol);
	  etah_mass_list[4].distr_list.push_back( eta_mass_rest/a_distr);
	  Ens_list[4].push_back(Ens);
	  

	}
      }
      else if(m=="M6") {
	C_beta[5].distr_list.push_back( C);
	if(Ens != "cB211b.072.64" && Ens != "cC211a.06.80" && Ens != "cA211a.12.48") {
	  if(Ens=="cD211a.054.96") {
	    Corr.Tmin = (int)(0.6/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else if(Ens=="cE211a.044.112") {
	    Corr.Tmin = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.5/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else crash("Ens: "+Ens+" not found");

	  
	  distr_t R5 = Corr.Fit_distr(ratio_corr_RED_interpol);
	  R[4].distr_list.push_back(R5);
	  aR[4].distr_list.push_back(a_distr);

	  mom_k_list[5].distr_list.push_back(mom_k);
	  m_etah_list[5].distr_list.push_back(to_m_etah);
	  aH_list[5].distr_list.push_back( a_distr);

	  DH_mass_list[5].distr_list.push_back( DH_mass_interpol);
	  etah_mass_list[5].distr_list.push_back( eta_mass_rest/a_distr);
	  Ens_list[5].push_back(Ens);
	
	}
	
      }




      

      
      

      

      //####################### PRINT ADDITIONAL INFO ON SCREEN  ########################################
      cout<<m<<": m_eta: "<<(eta_mass_rest/a_distr).ave()<<" +- "<<(eta_mass_rest/a_distr).err()<<" m_h: "<<(H_mass/a_distr).ave()<<" +- "<<(H_mass/a_distr).err()<<" Delta: "<<((H_mass-eta_mass_rest)/(a_distr)).ave()<<" +- "<<((H_mass-eta_mass_rest)/(a_distr)).err()<<endl;
      distr_t k = (H_mass*H_mass - eta_mass_rest*eta_mass_rest)/(2.0*a_distr*H_mass);
      cout<<"k_lat: "<<k.ave()<<" +- "<<k.err()<<endl;
      cout<<"ZV_ratio("<<m<<") : a^2: "<<(a_distr*a_distr).ave() <<" "<<(ZV_had/ZV).ave()<<"  "<<(ZV_had/ZV).err()<<endl;
      //#################################################################################################
      
      

    }
  }
  

  
  //GET AP_splitting
  for(int ir=0;ir<6;ir++) {
    
    AP_splitting[ir] = -0.0*etah_mass_list[ir];
    
    for(int iens=0;iens<(signed)Ens_list[ir].size(); iens++) {
      
      
      for(int is=0;is<=ir;is++) {
	//find ensemble corresponding to Ensemble
	int ens_id=-1;
	for(int s=0;s<(signed)Ens_list[is].size(); s++) {
	  if(Ens_list[is][s] == Ens_list[ir][iens] ) ens_id=s;
	}
	assert(ens_id != -1);

	cout<<"AP splitting("<<Ens_list[ir][iens]<<", ir="<<ir<<") : "<<AP_splitting[ir][iens].ave()<<endl;
	AP_splitting[ir].distr_list[iens] = AP_splitting[ir].distr_list[iens] + DH_mass_list[is].distr_list[ens_id];
	cout<<"DH_mass_list["<<is<<","<<ens_id<<"]: "<<DH_mass_list[is][ens_id].ave()<<endl;
	cout<<"AP splitting("<<Ens_list[ir][iens]<<", ir="<<ir<<") : "<<AP_splitting[ir][iens].ave()<<endl;
      }
    }
  }

  

  for(int ir=0;ir<6;ir++) {

    Print_To_File({}, {aH_list[ir].ave(), (AP_splitting[ir]-Mh_masses[ir+1]).ave(), AP_splitting[ir].err()}, "../data/heavy_radiative/Fits_R/H_mass_r_"+to_string(ir), "", "");

  }

  //now

  vector<distr_t_list> mom_k(6);
    for(int ir=0;ir<6;ir++) {
      mom_k[ir] = AP_splitting[ir] - Mh_masses[ir+1];
    }
  
    mom_k_list = mom_k;
  
  
  //m_etah_list = AP_splitting;

  //print C_beta

  for(int r=0;r<(signed)C_beta.size();r++) {
    Print_To_File({}, {a_distr_list.ave(), C_beta[r].ave(), C_beta[r].err()} ,  "../data/heavy_radiative/Fits_R/C_beta_"+to_string(r+1), "", "");
  }




  //continuum limit extrapolation

  class ipar {
  public:
    ipar() : FF(0.0), FF_err(0.0), FM1(0.0), FM1_err(0.0),  metac(0.0), metac_err(0.0), mhc(0.0), mhc_err(0.0) {}
    double FF, FF_err;
    double FM1, FM1_err;
    double metac, metac_err;
    double mhc, mhc_err;
    double a;
  };
  
  class fpar {
  public:
    fpar() {}
    fpar(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar  class constructor Vfloat par has size != 2");
      R=par[0];
      A=par[1];
      A2=par[2];
    }
    double R,A,A2;
  };

  int Nmeas= a_distr_list.size();

  //fit on mean values to get ch2
  bootstrap_fit<fpar,ipar> bf_mh(Njacks);
  bf_mh.set_warmup_lev(1); //sets warmup
  bf_mh.Set_number_of_measurements(Nmeas);
  bf_mh.Set_verbosity(1);
  bf_mh.Add_par("R", 1.0, 0.1);
  bf_mh.Add_par("A", 1.0, 0.1);
  bf_mh.Add_par("A2", 1.0, 0.1);
  //for mean values
  bootstrap_fit<fpar,ipar> bf_mh_ch2(1);
  bf_mh_ch2.set_warmup_lev(1); //sets warmup
  bf_mh_ch2.Set_number_of_measurements(Nmeas);
  bf_mh_ch2.Set_verbosity(1);
  bf_mh_ch2.Add_par("R", 1.0, 0.1);
  bf_mh_ch2.Add_par("A", 1.0, 0.1);
  bf_mh_ch2.Add_par("A2", 1.0, 0.1);
  bf_mh.Fix_par("A2", 0.0);
  bf_mh_ch2.Fix_par("A2",0.0);
  
  //##############################//
  
 

  string FIT_TAG="Hc";

  //ansatz
  bf_mh.ansatz=  [](const fpar &p, const ipar &ip) {
  		  
    return p.R + p.A*pow(ip.a,2) + p.A2*pow(ip.a,4);

  };

  
  bf_mh.measurement=  [&FIT_TAG ](const fpar &p, const ipar &ip) {
    
    if(FIT_TAG=="Hc") return ip.mhc;
    else if(FIT_TAG=="etac") return ip.metac;
    else if(FIT_TAG=="FM1") return ip.FM1;
    return ip.FF;
       
  };
  bf_mh.error=  [&FIT_TAG ](const fpar &p, const ipar &ip) {

    if(FIT_TAG=="Hc") return ip.mhc_err;
    else if(FIT_TAG=="etac") return ip.metac_err;
    else if(FIT_TAG=="FM1") return ip.FM1_err;
    return ip.FF_err;
    
  };

  bf_mh_ch2.ansatz= bf_mh.ansatz;
  bf_mh_ch2.measurement = bf_mh.measurement;
  bf_mh_ch2.error = bf_mh.error;

  int Nens= a_distr_list.size();

  //fill the data
  vector<vector<ipar>> data(Njacks);
  vector<vector<ipar>> data_ch2(1);
  //allocate space for output result
 
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= Fc.distr_list[iens].distr[ijack];
      data[ijack][iens].FF_err = Fc.err(iens);
      data[ijack][iens].FM1= Fm1.distr_list[iens].distr[ijack];
      data[ijack][iens].FM1_err = Fm1.err(iens);
      data[ijack][iens].metac = metac.distr_list[iens].distr[ijack];
      data[ijack][iens].metac_err= metac.err(iens);
      data[ijack][iens].mhc= mhc.distr_list[iens].distr[ijack];
      data[ijack][iens].mhc_err= mhc.err(iens);
      data[ijack][iens].a = a_distr_list.distr_list[iens].distr[ijack];
      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= Fc.ave(iens);
	data_ch2[ijack][iens].FF_err = Fc.err(iens);
	data_ch2[ijack][iens].FM1= Fm1.ave(iens);
	data_ch2[ijack][iens].FM1_err = Fm1.err(iens);
	data_ch2[ijack][iens].metac = metac.ave(iens);
	data_ch2[ijack][iens].metac_err= metac.err(iens);
	data_ch2[ijack][iens].mhc= mhc.ave(iens);
	data_ch2[ijack][iens].mhc_err= mhc.err(iens);
	data_ch2[ijack][iens].a = a_distr_list.ave(iens);
	
      }
    }
  }


  //append
  bf_mh.Append_to_input_par(data);
  bf_mh_ch2.Append_to_input_par(data_ch2);
  //fit
  boot_fit_data<fpar> Bt_fit_Fc;
  boot_fit_data<fpar> Bt_fit_Fc_ch2;
  FIT_TAG="";
  Bt_fit_Fc= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Fc_ch2= bf_mh_ch2.Perform_bootstrap_fit();

  boot_fit_data<fpar> Bt_fit_Fm1;
  boot_fit_data<fpar> Bt_fit_Fm1_ch2;
  FIT_TAG="FM1";
  Bt_fit_Fm1= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Fm1_ch2= bf_mh_ch2.Perform_bootstrap_fit();
  FIT_TAG="Hc";
  boot_fit_data<fpar> Bt_fit_Hc;
  boot_fit_data<fpar> Bt_fit_Hc_ch2;
  Bt_fit_Hc= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Hc_ch2= bf_mh_ch2.Perform_bootstrap_fit();
  FIT_TAG="etac";
  boot_fit_data<fpar> Bt_fit_etac;
  boot_fit_data<fpar> Bt_fit_etac_ch2;
  Bt_fit_etac= bf_mh.Perform_bootstrap_fit();
  Bt_fit_etac_ch2= bf_mh_ch2.Perform_bootstrap_fit();


  //set to zero the linear term and refit

  bf_mh.Fix_par("A",0.0);
  bf_mh_ch2.Fix_par("A",0.0);

  boot_fit_data<fpar> Bt_fit_Fc_const;
  boot_fit_data<fpar> Bt_fit_Fc_ch2_const;
  FIT_TAG="";
  Bt_fit_Fc_const= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Fc_ch2_const= bf_mh_ch2.Perform_bootstrap_fit();
  boot_fit_data<fpar> Bt_fit_Fm1_const;
  boot_fit_data<fpar> Bt_fit_Fm1_ch2_const;
  FIT_TAG="FM1";
  Bt_fit_Fm1_const= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Fm1_ch2_const= bf_mh_ch2.Perform_bootstrap_fit();
  FIT_TAG="Hc";
  boot_fit_data<fpar> Bt_fit_Hc_const;
  boot_fit_data<fpar> Bt_fit_Hc_ch2_const;
  Bt_fit_Hc_const= bf_mh.Perform_bootstrap_fit();
  Bt_fit_Hc_ch2_const= bf_mh_ch2.Perform_bootstrap_fit();
  FIT_TAG="etac";
  boot_fit_data<fpar> Bt_fit_etac_const;
  boot_fit_data<fpar> Bt_fit_etac_ch2_const;
  Bt_fit_etac_const= bf_mh.Perform_bootstrap_fit();
  Bt_fit_etac_ch2_const= bf_mh_ch2.Perform_bootstrap_fit();

  //Hc pars
  distr_t M_Hc, M_Hc_const, D_Hc;
  distr_t M_etac, M_etac_const, D_etac;
  distr_t FFc, FFc_const, D_FFc;
  distr_t FFm1, FFm1_const, D_FFm1;

  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) {

    //Hc
    M_Hc.distr.push_back( Bt_fit_Hc.par[ijack].R);
    M_Hc_const.distr.push_back( Bt_fit_Hc_const.par[ijack].R);
    D_Hc.distr.push_back( Bt_fit_Hc.par[ijack].A);
    //etac
    M_etac.distr.push_back( Bt_fit_etac.par[ijack].R);
    M_etac_const.distr.push_back( Bt_fit_etac_const.par[ijack].R);
    D_etac.distr.push_back( Bt_fit_etac.par[ijack].A);
    //Fc
    FFc.distr.push_back( Bt_fit_Fc.par[ijack].R);
    FFc_const.distr.push_back( Bt_fit_Fc_const.par[ijack].R);
    D_FFc.distr.push_back( Bt_fit_Fc.par[ijack].A);
    //Fm1
    FFm1.distr.push_back( Bt_fit_Fm1.par[ijack].R);
    FFm1_const.distr.push_back( Bt_fit_Fm1_const.par[ijack].R);
    D_FFm1.distr.push_back( Bt_fit_Fm1.par[ijack].A);
    
  }
  //get ch2
  double ch2_etac_lin = Bt_fit_etac_ch2.get_ch2_ave();
  double ch2_etac_const = Bt_fit_etac_ch2_const.get_ch2_ave();
  double ch2_Hc_lin = Bt_fit_Hc_ch2.get_ch2_ave();
  double ch2_Hc_const = Bt_fit_Hc_ch2_const.get_ch2_ave();
  double ch2_Fc_lin = Bt_fit_Fc.get_ch2_ave();
  double ch2_Fc_const = Bt_fit_Fc_const.get_ch2_ave();
  double ch2_Fm1_lin = Bt_fit_Fm1.get_ch2_ave();
  double ch2_Fm1_const = Bt_fit_Fm1_const.get_ch2_ave();

  //build fitting function

  //Dof(linear_fit)= 2 , Dof(constant_fit)=3

  int dof_lin= Nmeas - 2;
  int dof_const = Nmeas-1;
  

  double w_etac1 = exp(-0.5*(ch2_etac_lin - 2*dof_lin))/( exp(-0.5*(ch2_etac_lin - 2*dof_lin)) + exp(-0.5*(ch2_etac_const - 2*dof_const)));
  double w_etac2 = 1.0 - w_etac1;

  double w_Hc1 = exp(-0.5*(ch2_Hc_lin - 2*dof_lin))/( exp(-0.5*(ch2_Hc_lin - 2*dof_lin)) + exp(-0.5*(ch2_Hc_const - 2*dof_const)));
  double w_Hc2 = 1.0 - w_Hc1;

  double w_Fc1 = exp(-0.5*(ch2_Fc_lin - 2*dof_lin))/( exp(-0.5*(ch2_Fc_lin - 2*dof_lin)) + exp(-0.5*(ch2_Fc_const - 2*dof_const)));
  double w_Fc2 = 1.0 - w_Fc1;

  double w_Fm11 = exp(-0.5*(ch2_Fm1_lin - 2*dof_lin))/( exp(-0.5*(ch2_Fm1_lin - 2*dof_lin)) + exp(-0.5*(ch2_Fm1_const - 2*dof_const)));
  double w_Fm12 = 1.0 - w_Fm11;


  //how do we fit

  w_etac1=1.0; w_etac2=0.0;
  w_Hc1=1.0; w_Hc2=0.0;
  w_Fc1=1.0; w_Fc2=0.0;
  w_Fm11=1.0; w_Fm12=0.0;

  cout<<"reduced ch2 (etac): "<<ch2_etac_lin/dof_lin<<endl;
  cout<<"reduced ch2 (hc): "<<ch2_Hc_lin/dof_lin<<endl;
  cout<<"reduced ch2 (Fc): "<<ch2_Fc_lin/dof_lin<<endl;
  cout<<"reduced ch2 (Fm1): "<<ch2_Fm1_lin/dof_lin<<endl;

  distr_t M_Hc_ave= w_Hc1*M_Hc + w_Hc2*M_Hc_const;
  distr_t M_etac_ave= w_etac1*M_etac + w_etac2*M_etac_const;
  distr_t Fc_ave = w_Fc1*FFc + w_Fc2*FFc_const;
  distr_t Fm1_ave = w_Fm11*FFm1 + w_Fm12*FFm1_const;
  
  distr_t D_Hc_ave = w_Hc1*D_Hc;
  distr_t D_etac_ave = w_etac1*D_etac;
  distr_t D_FFc_ave= w_Fc1*D_FFc;
  distr_t D_FFm1_ave= w_Fm11*D_FFm1;
   

  
  //now combine and print the results


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list Metac_to_print(UseJack);
  distr_t_list MHc_to_print(UseJack);
  distr_t_list Fc_to_print(UseJack);
  distr_t_list Fm1_to_print(UseJack);
  Vfloat  a_to_print;

  //print in step of 0.01 fm;
  double a_max = 0.09;
  
  for(int i=0;i< (int)(a_max/0.001);i++) {
    double x = i*0.001;
    Metac_to_print.distr_list.push_back(  M_etac_ave + D_etac_ave*pow(x*fm_to_inv_Gev,2) );
    MHc_to_print.distr_list.push_back(  M_Hc_ave + D_Hc_ave*pow(x*fm_to_inv_Gev,2) );
    Fc_to_print.distr_list.push_back( Fc_ave + D_FFc_ave*pow(x*fm_to_inv_Gev,2));
    Fm1_to_print.distr_list.push_back( Fm1_ave + D_FFm1_ave*pow(x*fm_to_inv_Gev,2));
    a_to_print.push_back(x);
  }


  boost::filesystem::create_directory("../data/heavy_radiative/Fits_Hc");

  Print_To_File({}, {a_to_print, Metac_to_print.ave(), Metac_to_print.err()} , "../data/heavy_radiative/Fits_Hc/Metac_extr.fit", "", "");
  Print_To_File({}, {a_to_print, MHc_to_print.ave(), MHc_to_print.err()} , "../data/heavy_radiative/Fits_Hc/MHc_extr.fit", "", "");
  Print_To_File({}, {a_to_print, Fc_to_print.ave(), Fc_to_print.err() }, "../data/heavy_radiative/Fits_Hc/Fc.fit", "", "");
  Print_To_File({}, {a_to_print, Fm1_to_print.ave(), Fm1_to_print.err() }, "../data/heavy_radiative/Fits_Hc/Fm1.fit", "", "");

  //print data

  Print_To_File({} , { (a_distr_list/fm_to_inv_Gev).ave(), metac.ave(), metac.err(), mhc.ave(), mhc.err(), Fc.ave(), Fc.err(), Fm1.ave(), Fm1.err() }, "../data/heavy_radiative/Fits_Hc/Fc.data", "", "");


  cout<<"continuum limit extrapolation of k_mom_transf and m_etah"<<endl;

  for(int b=0;b<2;b++) {

    for(int ir=0;ir<(signed)m_etah_list.size();ir++) {

      if(b==0) cout<<"Extrapolation of metah: "<<ir+1<<endl;
      else cout<<"Extrapolation of mom_k: "<<ir+1<<endl;
       
      class ipar_r {
      public:
	ipar_r() : R(0.0), R_err(0.0), a(0.0) {}
	double R, R_err;
	double a;
      };
      
      class fpar_r {
      public:
	fpar_r() {}
	fpar_r(const Vfloat &par) {
	  if((signed)par.size() != 3) crash("In class fpar_r  class constructor Vfloat par has size != 3");
	  R=par[0];
	  A=par[1];
	  A2=par[2];
	}
	double R,A,A2;
      };
      
      int Nmeas_rat= m_etah_list[ir].size();
      
      
      //fit on mean values to get ch2
      bootstrap_fit<fpar_r,ipar_r> bf_r(Njacks);
      bf_r.set_warmup_lev(1); //sets warmup
      bf_r.Set_number_of_measurements(Nmeas_rat);
      bf_r.Set_verbosity(1);
      bf_r.Add_par("R", 1.0, 0.1);
      bf_r.Add_par("A", 1.0, 0.1);
      bf_r.Add_par("A2", 1.0, 0.1);
      //for mean values
      bootstrap_fit<fpar_r,ipar_r> bf_r_ch2(1);
      bf_r_ch2.set_warmup_lev(1); //sets warmup
      bf_r_ch2.Set_number_of_measurements(Nmeas_rat);
      bf_r_ch2.Set_verbosity(1);
      bf_r_ch2.Add_par("R", 1.0, 0.1);
      bf_r_ch2.Add_par("A", 1.0, 0.1);
      bf_r_ch2.Add_par("A2", 1.0, 0.1);

      if(ir > 3) {
	bf_r.Fix_par("A2",0.0);
	bf_r_ch2.Fix_par("A2",0.0);
      }
      
     
      //##############################//
      
      

   
       
      //ansatz
      bf_r.ansatz=  [](const fpar_r &p, const ipar_r &ip) {
	 
	return p.R + p.A*pow(ip.a,2) + p.A2*pow(ip.a,4);
	
      };
      
       
      bf_r.measurement =   [](const fpar_r &p, const ipar_r &ip) {
	 
	return ip.R;
	
      };
      bf_r.error=  [](const fpar_r &p, const ipar_r &ip) {
	 
	
	return ip.R_err;
	 
      };
      
      bf_r_ch2.ansatz= bf_r.ansatz;
      bf_r_ch2.measurement = bf_r.measurement;
      bf_r_ch2.error = bf_r.error;
       
      
      //fill the data
      vector<vector<ipar_r>> data(Njacks);
      vector<vector<ipar_r>> data_ch2(1);
      //allocate space for output result
      

      
      for(auto &data_iboot: data) data_iboot.resize(Nmeas_rat);
      for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas_rat);
      for(int ijack=0;ijack<Njacks;ijack++) {
       for(int iens=0;iens<Nmeas_rat;iens++) {
	 data[ijack][iens].R= (b==0)?(m_etah_list[ir].distr_list[iens].distr[ijack]):(mom_k_list[ir].distr_list[iens].distr[ijack]);
	 data[ijack][iens].R_err = (b==0)?(m_etah_list[ir].err(iens)):(mom_k_list[ir].err(iens));
	 data[ijack][iens].a = aH_list[ir].distr_list[iens].distr[ijack];
	 //mean values
	 if(ijack==0) {
	   data_ch2[ijack][iens].R= (b==0)?(m_etah_list[ir].ave(iens)):(mom_k_list[ir].ave(iens)); 
	   data_ch2[ijack][iens].R_err = (b==0)?(m_etah_list[ir].err(iens)):(mom_k_list[ir].ave(iens));
	   data_ch2[ijack][iens].a = aH_list[ir].ave(iens);
	   
	 }
       }
      }
      
      
      //append
      bf_r.Append_to_input_par(data);
      bf_r_ch2.Append_to_input_par(data_ch2);
      

      //fit
      boot_fit_data<fpar_r> Bt_fit_R;
      boot_fit_data<fpar_r> Bt_fit_R_ch2;
      Bt_fit_R= bf_r.Perform_bootstrap_fit();
      Bt_fit_R_ch2= bf_r_ch2.Perform_bootstrap_fit();

      distr_t r,  D_r, D2_r;
     
       //retrieve parameters
       for(int ijack=0;ijack<Njacks;ijack++) {
         
	 r.distr.push_back( Bt_fit_R.par[ijack].R);
	 D_r.distr.push_back( Bt_fit_R.par[ijack].A);
	 D2_r.distr.push_back( Bt_fit_R.par[ijack].A2);
	 
       }
       //get ch2
       double ch2_R_lin = Bt_fit_R.get_ch2_ave();


       //print fitting functions
       distr_t_list R_to_print(UseJack);
       Vfloat  a_to_print;
       //print in step of 0.01 fm;
       double a_max = 0.1;
       
       for(int i=0;i< (int)(a_max/0.001);i++) {
	 double x = i*0.001;
	 R_to_print.distr_list.push_back( r + D_r*pow(x*fm_to_inv_Gev,2)+ D2_r*pow(x*fm_to_inv_Gev,4) );
	 a_to_print.push_back(x);
       }
       
      
       boost::filesystem::create_directory("../data/heavy_radiative/Fits_R");

       string out_string= (b==0)?("metah_"+to_string(ir)):("mom_k_"+to_string(ir));
       
       Print_To_File({}, {a_to_print, (R_to_print).ave(), (R_to_print).err()} , "../data/heavy_radiative/Fits_R/"+out_string+"_extr.fit", "", "");

       if(b==0) {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (m_etah_list[ir]).ave(), (m_etah_list[ir]).err() }, "../data/heavy_radiative/Fits_R/metah_"+to_string(ir)+".data", "", "");
       }
       else {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (mom_k_list[ir]).ave(), (mom_k_list[ir]).err() }, "../data/heavy_radiative/Fits_R/mom_k_"+to_string(ir)+".data", "", "");
       }


    }
  }


  for(int a=0;a<(signed)m_etah_list_glb.size();a++) {

    Print_To_File({}, { (aH_glb_list[a]/fm_to_inv_Gev).ave(), m_etah_list_glb[a].ave(), m_etah_list_glb[a].err()}, "../data/heavy_radiative/Fits_R/metah_glb_"+to_string(a)+".data", "","");
    
  }

  
   //global mass-cont extrapolation

  //#######################################################################################################


   

  class ipar_r_glb_mass {
  public:
    ipar_r_glb_mass() : R(0.0), R_err(0.0), m_r(0.0), a(0.0) {}
    double R, R_err;
    double m_r;
    double a;
    int n;
  };
  
  class fpar_r_glb_mass {
  public:
    fpar_r_glb_mass() {}
    fpar_r_glb_mass(const Vfloat &par) {
      if((signed)par.size() != 12) crash("In class fpar_r  class constructor Vfloat par has size != 12");
      R1=par[0];
      R2=par[1];
      R3=par[2];
      R4=par[3];
      R5=par[4];
      R6=par[5];
      A=par[6];
      Am=par[7];
      Am2=par[8];
      A4=par[9];
      A4m2=par[10];
      A4m4=par[11];
    }
    double R1,R2,R3,R4,R5,R6,A, Am, Am2, A4, A4m2, A4m4;
  };

  int Nmeas_glb_mass=0;

  for(int r=0;r<(signed)m_etah_list.size();r++) {
    Nmeas_glb_mass += m_etah_list[r].size();
  }
  

  //fit on mean values to get ch2
  bootstrap_fit<fpar_r_glb_mass,ipar_r_glb_mass> bf_glb_mass(Njacks);
  bf_glb_mass.set_warmup_lev(1); //sets warmup
  bf_glb_mass.Set_number_of_measurements(Nmeas_glb_mass);
  bf_glb_mass.Set_verbosity(1);
  bf_glb_mass.Add_par("R1", 3.0, 0.1);
  bf_glb_mass.Add_par("R2", 3.5, 0.1);
  bf_glb_mass.Add_par("R3", 4.0, 0.1);
  bf_glb_mass.Add_par("R4", 4.5, 0.1);
  bf_glb_mass.Add_par("R5", 5.0, 0.1);
  bf_glb_mass.Add_par("R6", 5.5, 0.1);  
  bf_glb_mass.Add_par("A", 1.0, 0.1);
  bf_glb_mass.Add_par("Am", 1.0, 0.1);
  bf_glb_mass.Add_par("Am2", 1.0, 0.1);
  bf_glb_mass.Add_par("A4", 1.0, 0.1);
  bf_glb_mass.Add_par("A4m2", 1.0, 0.1);
  bf_glb_mass.Add_par("A4m4", 1.0, 0.1);
  //for mean values
  bootstrap_fit<fpar_r_glb_mass,ipar_r_glb_mass> bf_glb_mass_ch2(1);
  bf_glb_mass_ch2.set_warmup_lev(1); //sets warmup
  bf_glb_mass_ch2.Set_number_of_measurements(Nmeas_glb_mass);
  bf_glb_mass_ch2.Set_verbosity(1);
  bf_glb_mass_ch2.Add_par("R1", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("R2", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("R3", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("R4", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("R5", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("R6", 1.0, 0.1);  
  bf_glb_mass_ch2.Add_par("A", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("Am", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("Am2", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("A4", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("A4m2", 1.0, 0.1);
  bf_glb_mass_ch2.Add_par("A4m4", 1.0, 0.1);

  //##############################//

  bf_glb_mass.Fix_par("Am",0.0);
  bf_glb_mass_ch2.Fix_par("A4m4");



    //ansatz
    bf_glb_mass.ansatz=  [](const fpar_r_glb_mass &p, const ipar_r_glb_mass &ip) {

      double x=0.0;
      
      if(ip.n==1) { x = p.R1; return x;}
      else if(ip.n==2) x = p.R2;
      else if(ip.n==3) x = p.R3;
      else if(ip.n==4) x = p.R4;
      else if(ip.n==5) x = p.R5;
      else if(ip.n==6) x = p.R6;
      else crash("n not found");
      
      return x*( 1 + p.A*pow(ip.a,2) + p.Am2*pow(ip.a*ip.m_r,2) + p.Am*pow(ip.a,2)*ip.m_r +  p.A4*pow(ip.a,4) + p.A4m2*pow(ip.a*ip.m_r,2)*pow(ip.a,2) + p.A4m4*pow(ip.a*ip.m_r,4) ) ;
      
    };
    
    
    bf_glb_mass.measurement =   [](const fpar_r_glb_mass &p, const ipar_r_glb_mass &ip) {
      
      return ip.R;
      
    };
    bf_glb_mass.error=  [](const fpar_r_glb_mass &p, const ipar_r_glb_mass &ip) {
      
      
      return ip.R_err;
      
    };
    
    bf_glb_mass_ch2.ansatz= bf_glb_mass.ansatz;
    bf_glb_mass_ch2.measurement = bf_glb_mass.measurement;
    bf_glb_mass_ch2.error = bf_glb_mass.error;


    //fill the data
    vector<vector<ipar_r_glb_mass>> data_glb_mass(Njacks);
    vector<vector<ipar_r_glb_mass>> data_glb_mass_ch2(1);
    //allocate space for output result



     for(auto &data_iboot: data_glb_mass) data_iboot.resize(Nmeas_glb_mass);
     for(auto &data_iboot: data_glb_mass_ch2) data_iboot.resize(Nmeas_glb_mass);

      
     for(int ijack=0;ijack<Njacks;ijack++) {
       int it=0;
       for(int rr=0;rr<(signed)m_etah_list.size(); rr++) {
	 for(int iens=0;iens<(signed)m_etah_list[rr].size();iens++) {
	   data_glb_mass[ijack][it].R= m_etah_list[rr].distr_list[iens].distr[ijack]; 
	   data_glb_mass[ijack][it].R_err = m_etah_list[rr].err(iens);
	   data_glb_mass[ijack][it].m_r = pow(lambda, rr+1);
	   data_glb_mass[ijack][it].a = aH_list[rr].distr_list[iens].distr[ijack];
	   data_glb_mass[ijack][it].n = rr+1;
	   //mean values
	   if(ijack==0) {
	     data_glb_mass_ch2[ijack][it].R= m_etah_list[rr].ave(iens); 
	     data_glb_mass_ch2[ijack][it].R_err = m_etah_list[rr].err(iens);
	     data_glb_mass_ch2[ijack][it].m_r = pow(lambda, rr+1);
	     data_glb_mass_ch2[ijack][it].a = aH_list[rr].ave(iens);
	     data_glb_mass_ch2[ijack][it].n = rr+1;
	   }
	   it++;
	 }
       }
     }

  

     //append
     bf_glb_mass.Append_to_input_par(data_glb_mass);
     bf_glb_mass_ch2.Append_to_input_par(data_glb_mass_ch2);
     
     //fit
     boot_fit_data<fpar_r_glb_mass> Bt_fit_glb_mass;
     boot_fit_data<fpar_r_glb_mass> Bt_fit_glb_mass_ch2;
     Bt_fit_glb_mass= bf_glb_mass.Perform_bootstrap_fit();
     Bt_fit_glb_mass_ch2= bf_glb_mass_ch2.Perform_bootstrap_fit();


     distr_t pR1, pR2, pR3, pR4, pR5, pR6, pHA, pHAm,  pHAm2, pHA4, pHA4m2, pHA4m4;
     
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {
              
       pR1.distr.push_back( Bt_fit_glb_mass.par[ijack].R1);
       pR2.distr.push_back( Bt_fit_glb_mass.par[ijack].R2);
       pR3.distr.push_back( Bt_fit_glb_mass.par[ijack].R3);
       pR4.distr.push_back( Bt_fit_glb_mass.par[ijack].R4);
       pR5.distr.push_back( Bt_fit_glb_mass.par[ijack].R5);
       pR6.distr.push_back( Bt_fit_glb_mass.par[ijack].R6);
    
       pHA.distr.push_back( Bt_fit_glb_mass.par[ijack].A);
       pHAm.distr.push_back( Bt_fit_glb_mass.par[ijack].Am);
       pHAm2.distr.push_back( Bt_fit_glb_mass.par[ijack].Am2);
       pHA4.distr.push_back( Bt_fit_glb_mass.par[ijack].A4);
       pHA4m2.distr.push_back( Bt_fit_glb_mass.par[ijack].A4m2);
       pHA4m4.distr.push_back( Bt_fit_glb_mass.par[ijack].A4m4);
           
     }
     //get ch2
     double ch2_glb_mass = Bt_fit_glb_mass_ch2.get_ch2_ave();

     cout<<"ch2/dof [glb_mass fit]: "<<ch2_glb_mass/(Nmeas_glb_mass - bf_glb_mass.Get_number_of_fit_pars())<<endl;



     //print fit function at non-zero lattice spacing for each mass

     for(int rr=0; rr<(signed)m_etah_list.size();rr++) {

       //print fitting functions for MP, FP, MP_ov_FP
       distr_t_list func_to_print(UseJack);
       Vfloat  a_to_print;
       //print in step of 0.01 fm;
       double a_max = 0.09;

       double mr= pow(lambda, rr+1);
     
       for(int i=0;i< (int)(a_max/0.001);i++) {
	 double x = i*0.001;
	 double a= x*fm_to_inv_Gev;
	 distr_t cont;
	 if(rr==0) {cont = pR1;}
	 else if(rr==1) {cont=pR2;}
	 else if(rr==2) {cont=pR3;}
	 else if(rr==3) {cont=pR4;}
	 else if(rr==4) {cont=pR5;}
	 else if(rr==5) {cont=pR6;}
	 else crash("crash cannot find rr");
	 
	 func_to_print.distr_list.push_back( cont*( 1  + pHA*pow(a,2) + pHAm*pow(a,2)*mr+ pHAm2*pow(a*mr,2) + pHA4*pow(a,4) + pHA4m2*pow(a*mr,2)*pow(a,2) +  pHA4m4*pow(a*mr,4) )); 
	 a_to_print.push_back(x);
       }

       Print_To_File({}, {a_to_print, (func_to_print).ave(), (func_to_print).err()} , "../data/heavy_radiative/Fits_R/m_etah_"+to_string(rr)+"_extr_glb_mass.fit", "", "");
     }
        
     
  //#######################################################################################################



  
  //now extrapolate to the continuum limit the ratios

  cout<<"Continuum extrapolation of ratios..."<<endl;


  distr_t_list R_extr(UseJack);


  
  for(int ir=0;ir<(signed)R.size()-1;ir++) {
    cout<<"Extrapolation ratio: "<<ir+1<<endl;

    class ipar_r {
    public:
      ipar_r() : R(0.0), R_err(0.0), a(0.0) {}
      double R, R_err;
      double a;
    };
    
    class fpar_r {
    public:
      fpar_r() {}
      fpar_r(const Vfloat &par) {
	if((signed)par.size() != 3) crash("In class fpar_r  class constructor Vfloat par has size != 3");
	R=par[0];
	A=par[1];
	A2=par[2];
      }
      double R,A,A2;
    };
    
    int Nmeas_rat= R[ir].size();


     //fit on mean values to get ch2
    bootstrap_fit<fpar_r,ipar_r> bf_r(Njacks);
    bf_r.set_warmup_lev(1); //sets warmup
    bf_r.Set_number_of_measurements(Nmeas_rat);
    bf_r.Set_verbosity(1);
    bf_r.Add_par("R", 1.0, 0.1);
    bf_r.Add_par("A", 1.0, 0.1);
    bf_r.Add_par("A2", 1.0, 0.1);
    //for mean values
    bootstrap_fit<fpar_r,ipar_r> bf_r_ch2(1);
    bf_r_ch2.set_warmup_lev(1); //sets warmup
    bf_r_ch2.Set_number_of_measurements(Nmeas_rat);
    bf_r_ch2.Set_verbosity(1);
    bf_r_ch2.Add_par("R", 1.0, 0.1);
    bf_r_ch2.Add_par("A", 1.0, 0.1);
    bf_r_ch2.Add_par("A2", 1.0, 0.1);
    
    bf_r.Fix_par("A2", 0.0);
    bf_r_ch2.Fix_par("A2",0.0);
    //##############################//
    
    

   
    
    //ansatz
    bf_r.ansatz=  [](const fpar_r &p, const ipar_r &ip) {
      
      return p.R + p.A*pow(ip.a,2) + p.A2*pow(ip.a,4);
      
    };
    
    
    bf_r.measurement =   [](const fpar_r &p, const ipar_r &ip) {
      
      return ip.R;
      
    };
    bf_r.error=  [](const fpar_r &p, const ipar_r &ip) {
      
      
      return ip.R_err;
      
    };
    
    bf_r_ch2.ansatz= bf_r.ansatz;
    bf_r_ch2.measurement = bf_r.measurement;
    bf_r_ch2.error = bf_r.error;
    
   
    //fill the data
    vector<vector<ipar_r>> data(Njacks);
    vector<vector<ipar_r>> data_ch2(1);
    //allocate space for output result



     for(auto &data_iboot: data) data_iboot.resize(Nmeas_rat);
     for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas_rat);
     for(int ijack=0;ijack<Njacks;ijack++) {
       for(int iens=0;iens<Nmeas_rat;iens++) {
	 data[ijack][iens].R= R[ir].distr_list[iens].distr[ijack];
	 data[ijack][iens].R_err = R[ir].err(iens);
	 data[ijack][iens].a = aR[ir].distr_list[iens].distr[ijack];
	 //mean values
	 if(ijack==0) {
	   data_ch2[ijack][iens].R= R[ir].ave(iens); 
	   data_ch2[ijack][iens].R_err = R[ir].err(iens);
	   data_ch2[ijack][iens].a = aR[ir].ave(iens);
	   
	 }
       }
     }
     

     //append
     bf_r.Append_to_input_par(data);
     bf_r_ch2.Append_to_input_par(data_ch2);


     //fit
     boot_fit_data<fpar_r> Bt_fit_R;
     boot_fit_data<fpar_r> Bt_fit_R_ch2;
     Bt_fit_R= bf_r.Perform_bootstrap_fit();
     Bt_fit_R_ch2= bf_r_ch2.Perform_bootstrap_fit();
     bf_r.Fix_par("A", 0.0);
     bf_r_ch2.Fix_par("A",0.0);
     boot_fit_data<fpar_r> Bt_fit_R_const;
     boot_fit_data<fpar_r> Bt_fit_R_ch2_const;
     Bt_fit_R_const= bf_r.Perform_bootstrap_fit();
     Bt_fit_R_ch2_const= bf_r_ch2.Perform_bootstrap_fit();


     distr_t r, r_const, D_r;
     
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {
              
       r.distr.push_back( Bt_fit_R.par[ijack].R);
       r_const.distr.push_back( Bt_fit_R_const.par[ijack].R);
       D_r.distr.push_back( Bt_fit_R.par[ijack].A);
       
     }
     //get ch2
     double ch2_R_lin = Bt_fit_R.get_ch2_ave();
     double ch2_R_const = Bt_fit_R_const.get_ch2_ave();

     int dof_l= Nmeas_rat - 2;
     int dof_c = Nmeas_rat-1;
          
     double w_R1 = exp(-0.5*(ch2_R_lin - 2*dof_l))/( exp(-0.5*(ch2_R_lin - 2*dof_l)) + exp(-0.5*(ch2_R_const - 2*dof_c)));
     double w_R2 = 1.0 - w_R1;

     


     //exclude constant fits
     w_R1=1.0;
     w_R2=0.0;
    


     //how do we fit
     
     distr_t R_ave = w_R1*r + w_R2*r_const;
     distr_t D_R_ave = w_R1*D_r;


     //print fitting functions for MP, FP, MP_ov_FP
     distr_t_list R_to_print(UseJack);
     Vfloat  a_to_print;
     //print in step of 0.01 fm;
     double a_max = 0.09;
     
     for(int i=0;i< (int)(a_max/0.001);i++) {
       double x = i*0.001;
       R_to_print.distr_list.push_back( R_ave + D_R_ave*pow(x*fm_to_inv_Gev,2));
       a_to_print.push_back(x);
     }

     double F_R = lambda*Get_F_R(ir+1)/Get_F_R(ir);
     
     

     boost::filesystem::create_directory("../data/heavy_radiative/Fits_R");
     
     Print_To_File({}, {a_to_print, (R_to_print*F_R).ave(), (R_to_print*F_R).err()} , "../data/heavy_radiative/Fits_R/R"+to_string(ir)+"_extr.fit", "", "");
     
     Print_To_File({} , { (aR[ir]/fm_to_inv_Gev).ave(), (R[ir]*F_R).ave(), (R[ir]*F_R).err() }, "../data/heavy_radiative/Fits_R/R"+to_string(ir)+".data", "", "");

     R_extr.distr_list.push_back( R_ave*F_R);

  }

  //global mass-cont extrapolation

  //#######################################################################################################


   

  class ipar_r_glb {
  public:
    ipar_r_glb() : R(0.0), R_err(0.0), m_r(0.0), a(0.0) {}
    double R, R_err;
    double m_r;
    double a;
  };
  
  class fpar_r_glb {
  public:
    fpar_r_glb() {}
    fpar_r_glb(const Vfloat &par) {
      if((signed)par.size() != 8) crash("In class fpar_r  class constructor Vfloat par has size != 7");
      R=par[0];
      Rm=par[1];
      Rm2=par[2];
      Rm3=par[3];
      A=par[4];
      Am2=par[5];
      A4=par[6];
      A4m4=par[7];
    }
    double R,Rm, Rm2, Rm3, A,Am2, A4, A4m4;
  };

  int Nmeas_glb=0;

  for(int r=0;r<(signed)R.size();r++) {
    Nmeas_glb += R[r].size();
  }
  

  //fit on mean values to get ch2
  bootstrap_fit<fpar_r_glb,ipar_r_glb> bf_glb(Njacks);
  bf_glb.set_warmup_lev(1); //sets warmup
  bf_glb.Set_number_of_measurements(Nmeas_glb);
  bf_glb.Set_verbosity(1);
  bf_glb.Add_par("R", 1.0, 0.1);
  bf_glb.Add_par("Rm", 1.0, 0.1);
  bf_glb.Add_par("Rm2", 1.0, 0.1);
  bf_glb.Add_par("Rm3", 1.0, 0.1);
  bf_glb.Add_par("A", 1.0, 0.1);
  bf_glb.Add_par("Am2", 1.0, 0.1);
  bf_glb.Add_par("A4", 1.0, 0.1);
  bf_glb.Add_par("A4m4", 1.0, 0.1);
  //for mean values
  bootstrap_fit<fpar_r_glb,ipar_r_glb> bf_glb_ch2(1);
  bf_glb_ch2.set_warmup_lev(1); //sets warmup
  bf_glb_ch2.Set_number_of_measurements(Nmeas_glb);
  bf_glb_ch2.Set_verbosity(1);
  bf_glb_ch2.Add_par("R", 1.0, 0.1);
  bf_glb_ch2.Add_par("Rm", 1.0, 0.1);
  bf_glb_ch2.Add_par("Rm2", 1.0, 0.1);
  bf_glb_ch2.Add_par("Rm3", 1.0, 0.1);
  bf_glb_ch2.Add_par("A", 1.0, 0.1);
  bf_glb_ch2.Add_par("Am2", 1.0, 0.1);
  bf_glb_ch2.Add_par("A4", 1.0, 0.1);
  bf_glb_ch2.Add_par("A4m4", 1.0, 0.1);


  bf_glb.Fix_par("R", 1.0);
  bf_glb_ch2.Fix_par("R", 1.0);
  bf_glb.Fix_par("Rm2", 0.0);
  bf_glb_ch2.Fix_par("Rm2",0.0);
  bf_glb.Fix_par("Rm3", 0.0);
  bf_glb_ch2.Fix_par("Rm3",0.0);
  bf_glb.Fix_par("A4m4", 0.0);
  bf_glb_ch2.Fix_par("A4m4",0.0);
  bf_glb.Fix_par("A4", 0.0);
  bf_glb_ch2.Fix_par("A4",0.0);
  //##############################//



    //ansatz
    bf_glb.ansatz=  [](const fpar_r_glb &p, const ipar_r_glb &ip) {
      
      return p.R + p.Rm/pow(ip.m_r,1) + p.Rm2/pow(ip.m_r,2) + p.Rm3/pow(ip.m_r,3) + p.A*pow(ip.a,2) + p.Am2*pow(ip.a*ip.m_r,2) + p.A4*pow(ip.a,4) + p.A4m4*pow(ip.a*ip.m_r,4);
      
    };
    
    
    bf_glb.measurement =   [](const fpar_r_glb &p, const ipar_r_glb &ip) {
      
      return ip.R;
      
    };
    bf_glb.error=  [](const fpar_r_glb &p, const ipar_r_glb &ip) {
      
      
      return ip.R_err;
      
    };
    
    bf_glb_ch2.ansatz= bf_glb.ansatz;
    bf_glb_ch2.measurement = bf_glb.measurement;
    bf_glb_ch2.error = bf_glb.error;


    //fill the data
    vector<vector<ipar_r_glb>> data_glb(Njacks);
    vector<vector<ipar_r_glb>> data_glb_ch2(1);
    //allocate space for output result



     for(auto &data_iboot: data_glb) data_iboot.resize(Nmeas_glb);
     for(auto &data_iboot: data_glb_ch2) data_iboot.resize(Nmeas_glb);

     Vfloat FTS;
     for(int rr=0;rr<(signed)R.size(); rr++) FTS.push_back( Get_F_R(rr+1)/Get_F_R(rr));
   
     for(int ijack=0;ijack<Njacks;ijack++) {
       int it=0;
       for(int rr=0;rr<(signed)R.size(); rr++) {
	 for(int iens=0;iens<(signed)R[rr].size();iens++) {
	   data_glb[ijack][it].R= R[rr].distr_list[iens].distr[ijack]*lambda*FTS[rr]; 
	   data_glb[ijack][it].R_err = R[rr].err(iens)*lambda*FTS[rr];
	   data_glb[ijack][it].m_r = pow(lambda, rr+1);
	   data_glb[ijack][it].a = aR[rr].distr_list[iens].distr[ijack];
	   //mean values
	   if(ijack==0) {
	     data_glb_ch2[ijack][it].R= R[rr].ave(iens)*lambda*FTS[rr]; 
	     data_glb_ch2[ijack][it].R_err = R[rr].err(iens)*lambda*FTS[rr];
	     data_glb_ch2[ijack][it].m_r = pow(lambda, rr+1);
	     data_glb_ch2[ijack][it].a = aR[rr].ave(iens);
	   }
	   it++;
	 }
       }
     }

  

     //append
     bf_glb.Append_to_input_par(data_glb);
     bf_glb_ch2.Append_to_input_par(data_glb_ch2);
     
     //fit
     boot_fit_data<fpar_r_glb> Bt_fit_glb;
     boot_fit_data<fpar_r_glb> Bt_fit_glb_ch2;
     Bt_fit_glb= bf_glb.Perform_bootstrap_fit();
     Bt_fit_glb_ch2= bf_glb_ch2.Perform_bootstrap_fit();


     distr_t pR, pRm, pRm2, pRm3, pA, pAm2, pA4, pA4m4;
     
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {
              
       pR.distr.push_back( Bt_fit_glb.par[ijack].R);

       pRm.distr.push_back( Bt_fit_glb.par[ijack].Rm);
       pRm2.distr.push_back( Bt_fit_glb.par[ijack].Rm2);
       pRm3.distr.push_back( Bt_fit_glb.par[ijack].Rm3);

       pA.distr.push_back( Bt_fit_glb.par[ijack].A);
       pAm2.distr.push_back( Bt_fit_glb.par[ijack].Am2);
       pA4.distr.push_back( Bt_fit_glb.par[ijack].A4);
       pA4m4.distr.push_back( Bt_fit_glb.par[ijack].A4m4);
       
     }
     //get ch2
     double ch2_glb = Bt_fit_glb_ch2.get_ch2_ave();

     cout<<"ch2/dof [glb fit]: "<<ch2_glb/(Nmeas_glb - bf_glb.Get_number_of_fit_pars())<<endl;



     //print fit function at non-zero lattice spacing for each mass

     for(int rr=0; rr<(signed)R.size();rr++) {

       //print fitting functions for MP, FP, MP_ov_FP
       distr_t_list R_to_print(UseJack);
       Vfloat  a_to_print;
       //print in step of 0.01 fm;
       double a_max = 0.09;

       double mr= pow(lambda, rr+1);
     
       for(int i=0;i< (int)(a_max/0.001);i++) {
	 double x = i*0.001;
	 double a= x*fm_to_inv_Gev;
	 R_to_print.distr_list.push_back( pR + pRm/mr + pRm2/pow(mr,2) + pRm3/pow(mr,3)  + pA*pow(a,2) + pAm2*pow(a*mr,2) + pA4*pow(a,4) + pA4m4*pow(a*mr,4)); 
	 a_to_print.push_back(x);
       }

       Print_To_File({}, {a_to_print, (R_to_print).ave(), (R_to_print).err()} , "../data/heavy_radiative/Fits_R/R"+to_string(rr)+"_extr_glb.fit", "", "");
     }
        

     //print fit function in the continuum as a function of the mass


     Vfloat  linv_glb_to_print;
     //print in step of 0.001;
     double linv_glb_max = 1.1;
     distr_t_list ratios_glb_to_print(UseJack);
     
     for(int i=0;i< (int)(linv_glb_max/0.001);i++) {
       double x = i*0.001;
       ratios_glb_to_print.distr_list.push_back( 1.0 + pRm*x + pRm2*pow(x,2) + pRm3*pow(x,3) );
       linv_glb_to_print.push_back(x);
     }

     Print_To_File({}, {linv_glb_to_print, (ratios_glb_to_print).ave(), (ratios_glb_to_print).err()} , "../data/heavy_radiative/Fits_B/ratio_extr_glb.fit", "", "");

     
  //#######################################################################################################










  
 
  //perform b-quark mass extrapolation from separate cont. extr.


  int Nmeas_R= R_extr.size();

  


  class ipar_r {
    public:
      ipar_r() : R(0.0), R_err(0.0), m_r(0.0) {}
      double R, R_err;
      double m_r;
    };
    
    class fpar_r {
    public:
      fpar_r() {}
      fpar_r(const Vfloat &par) {
	if((signed)par.size() != 3) crash("In class fpar_r  class constructor Vfloat par has size != 3");
	R=par[0];
	A=par[1];
	A2=par[2];
      }
      double R,A,A2;
    };
    
    //fit on mean values to get ch2
    bootstrap_fit<fpar_r,ipar_r> bf_r(Njacks);
    bf_r.set_warmup_lev(1); //sets warmup
    bf_r.Set_number_of_measurements(Nmeas_R);
    bf_r.Set_verbosity(1);
    bf_r.Add_par("R", 1.0, 0.1);
    bf_r.Add_par("A", 1.0, 0.1);
    bf_r.Add_par("A2", 1.0, 0.1);
    //for mean values
    bootstrap_fit<fpar_r,ipar_r> bf_r_ch2(1);
    bf_r_ch2.set_warmup_lev(1); //sets warmup
    bf_r_ch2.Set_number_of_measurements(Nmeas_R);
    bf_r_ch2.Set_verbosity(1);
    bf_r_ch2.Add_par("R", 1.0, 0.1);
    bf_r_ch2.Add_par("A", 1.0, 0.1);
    bf_r_ch2.Add_par("A2", 1.0, 0.1);

    bf_r.Fix_par("R", 1.0);
    bf_r_ch2.Fix_par("R", 1.0);
    bf_r.Fix_par("A2", 0.0);
    bf_r_ch2.Fix_par("A2",0.0);
    //##############################//


    //ansatz
    bf_r.ansatz=  [](const fpar_r &p, const ipar_r &ip) {
      
      return p.R + p.A/pow(ip.m_r,1) + p.A2/pow(ip.m_r,2);
      
    };
    
    
    bf_r.measurement =   [](const fpar_r &p, const ipar_r &ip) {
      
      return ip.R;
      
    };
    bf_r.error=  [](const fpar_r &p, const ipar_r &ip) {
      
      
      return ip.R_err;
      
    };
    
    bf_r_ch2.ansatz= bf_r.ansatz;
    bf_r_ch2.measurement = bf_r.measurement;
    bf_r_ch2.error = bf_r.error;


    //fill the data
    vector<vector<ipar_r>> data_r(Njacks);
    vector<vector<ipar_r>> data_r_ch2(1);
    //allocate space for output result



     for(auto &data_iboot: data_r) data_iboot.resize(Nmeas_R);
     for(auto &data_iboot: data_r_ch2) data_iboot.resize(Nmeas_R);
     for(int ijack=0;ijack<Njacks;ijack++) {
       for(int iens=0;iens<Nmeas_R;iens++) {
	 data_r[ijack][iens].R= R_extr.distr_list[iens].distr[ijack];
	 data_r[ijack][iens].R_err = R_extr.err(iens);
	 data_r[ijack][iens].m_r = pow(lambda, iens+1);
	 //mean values
	 if(ijack==0) {
	   data_r_ch2[ijack][iens].R= R_extr.ave(iens); 
	   data_r_ch2[ijack][iens].R_err = R_extr.err(iens);
	   data_r_ch2[ijack][iens].m_r = pow(lambda, iens+1);
	   
	 }
       }
     }
     

     //append
     bf_r.Append_to_input_par(data_r);
     bf_r_ch2.Append_to_input_par(data_r_ch2);


     //fit
     boot_fit_data<fpar_r> Bt_fit_R;
     boot_fit_data<fpar_r> Bt_fit_R_ch2;
     Bt_fit_R= bf_r.Perform_bootstrap_fit();
     Bt_fit_R_ch2= bf_r_ch2.Perform_bootstrap_fit();
     bf_r.Release_par("A2");
     bf_r_ch2.Release_par("A2");
     boot_fit_data<fpar_r> Bt_fit_R2;

     boot_fit_data<fpar_r> Bt_fit_R2_ch2;
     Bt_fit_R2= bf_r.Perform_bootstrap_fit();
     Bt_fit_R2_ch2= bf_r_ch2.Perform_bootstrap_fit();


     distr_t K, K_2, J;
     
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {
              
       K.distr.push_back( Bt_fit_R.par[ijack].A);
       K_2.distr.push_back( Bt_fit_R2.par[ijack].A);
       J.distr.push_back( Bt_fit_R2.par[ijack].A2);
       
     }
     //get ch2
     double ch2_R = Bt_fit_R_ch2.get_ch2_ave();
     double ch2_R2 = Bt_fit_R2_ch2.get_ch2_ave();

     int dof_R= Nmeas_R - 1;
     int dof_R2 = Nmeas_R-2;
          
     double w_R = exp(-0.5*(ch2_R - 2*dof_R))/( exp(-0.5*(ch2_R - 2*dof_R)) + exp(-0.5*(ch2_R2 - 2*dof_R2)));
     double w_R2 = 1.0 - w_R;


     //how do we fit
     
     distr_t K_ave = w_R*K + w_R2*K_2;
     distr_t J_ave = w_R2*J;


     
     //print fitting functions for MP, FP, MP_ov_FP
     distr_t_list ratios_to_print(UseJack);
     Vfloat  linv_to_print;
     //print in step of 0.001;
     double linv_max = 1.1;
     
     for(int i=0;i< (int)(linv_max/0.001);i++) {
       double x = i*0.001;
       ratios_to_print.distr_list.push_back( 1.0 + K_ave*x + J_ave*x*x);
       linv_to_print.push_back(x);
     }

     auto FUNC_R = [&K_ave, &J_ave] (int n)  -> distr_t {
       
       return 1.0+ K_ave/pow(lambda,n) + J_ave/pow(lambda,2*n);
     };

     auto FUNC_R_glb = [&pRm, &pRm2, &pRm3] (int n) -> distr_t {

       return 1.0 + pRm/pow(lambda,n) + pRm2/pow(lambda,2*n) + pRm3/pow(lambda,3*n);
     };
     

     boost::filesystem::create_directory("../data/heavy_radiative/Fits_B");

     Vfloat l_list;
     for(int i=0;i<Nmeas_R;i++) l_list.push_back( pow(lambda, -1.0*(i+1)));
     
     Print_To_File({}, {linv_to_print, (ratios_to_print).ave(), (ratios_to_print).err()} , "../data/heavy_radiative/Fits_B/ratio_extr.fit", "", "");
     
     Print_To_File({} , { l_list, R_extr.ave(), R_extr.err() }, "../data/heavy_radiative/Fits_B/ratio.data", "", "");

     //print result
     
    
     //final result for the form factor

     distr_t Fb_ave= Fm1_ave;

    
     Fb_ave = Fb_ave*FUNC_R(1)*FUNC_R(2)*FUNC_R(3)*FUNC_R(4)*FUNC_R(5)*FUNC_R(6)*FUNC_R(7);

     distr_t Fb_glb= Fm1_ave*FUNC_R_glb(1)*FUNC_R_glb(2)*FUNC_R_glb(3)*FUNC_R_glb(4)*FUNC_R_glb(5)*FUNC_R_glb(6)*FUNC_R_glb(7);


     //multiply Fb_ave by trivial factors

     if(TYPE=="WEAK") {

       Fb_ave = Fb_ave/pow(lambda,7);

       Fb_ave = Fb_ave*Get_F_R(0)/Get_F_R(7);

       Fb_glb = Fb_glb/pow(lambda,7);
       Fb_glb = Fb_glb*Get_F_R(0)/Get_F_R(7);

     }

     if(TYPE=="STRONG") {

       Fb_ave = Fb_ave/pow(lambda,7);

       Fb_ave = Fb_ave*sqrt(m_hb/m_hc);

       Fb_glb = Fb_glb/pow(lambda,7);
       Fb_glb = Fb_glb*sqrt(m_hb/m_hc);

     }

     cout<<"Nmeas_R: "<<Nmeas_R<<endl;

     cout<<"F1(b): "<<(0.5*Fb_ave).ave()<<" "<<(0.5*Fb_ave).err()<<endl;
     cout<<"F1_glb(b): "<<(0.5*Fb_glb).ave()<<" "<<(0.5*Fb_glb).err()<<endl;
     cout<<"F1(c): "<<(0.5*Fc_ave).ave()<<" "<<(0.5*Fc_ave).err()<<endl;


     cout<<"Gamma[ hb -> etab gamma ]: "<< (1e6)*(2.0/3.0)*Qb*Qb*alpha*(0.5*0.5*Fb_ave*Fb_ave).ave()*( pow(m_hb,2) - pow(m_etab,2))/(m_hb)<<" "<<(1e6)*(2.0/3.0)*Qb*Qb*alpha*(0.5*0.5*Fb_ave*Fb_ave).err()*( pow(m_hb,2) - pow(m_etab,2))/(m_hb)<<" [KeV]"<<endl;

     cout<<"Gamma[ hb -> etab gamma ]_glb: "<< (1e6)*(2.0/3.0)*Qb*Qb*alpha*(0.5*0.5*Fb_glb*Fb_glb).ave()*( pow(m_hb,2) - pow(m_etab,2))/(m_hb)<<" "<<(1e6)*(2.0/3.0)*Qb*Qb*alpha*(0.5*0.5*Fb_glb*Fb_glb).err()*( pow(m_hb,2) - pow(m_etab,2))/(m_hb)<<" [KeV]"<<endl;

     cout<<"Gamma[ hc -> etac gamma ]: "<< (1e6)*(2.0/3.0)*Qc*Qc*alpha*(0.5*0.5*Fc_ave*Fc_ave).ave()*( pow(m_hc,2) - pow(m_etac,2))/(m_hc)<<" "<<(1e6)*(2.0/3.0)*Qc*Qc*alpha*(0.5*0.5*Fc_ave*Fc_ave).err()*( pow(m_hc,2) - pow(m_etac,2))/(m_hc)<<" [KeV]"<<endl;

     
     

  



  return;

}



void heavy_radiative_test() {



  //test


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

  


  vector<string> M({"M1", "M2", "M3", "M4", "M5", "M6"});


   for(int im=0;im<(signed)M.size();im++) {

    

    string m=M[im];

    boost::filesystem::create_directory("../data/heavy_radiative_test");
    
    boost::filesystem::create_directory("../data/heavy_radiative_test/"+m);
    
    
    data_t P5P5_sm_REST;
    data_t P5P5_sm_BIS_REST;
    
    P5P5_sm_REST.Read("../heavy_radiative_test", "mes_contr_2PT_"+m, "P5P5", Sort_light_confs);
    P5P5_sm_BIS_REST.Read("../heavy_radiative_test", "mes_contr_2PT_"+m+"_BIS", "P5P5", Sort_light_confs);

    //set up statistical analysis
    CorrAnalysis Corr(UseJack,Njacks,100);
    Corr.Nt= P5P5_sm_REST.nrows[0];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;


    distr_t_list eff_mass_P5P5 = Corr.effective_mass_t(P5P5_sm_REST.col(0)[0], "../data/heavy_radiative_test/"+m+"/eta_1");
    distr_t_list eff_mass_BIS_P5P5 = Corr.effective_mass_t(P5P5_sm_BIS_REST.col(0)[0], "../data/heavy_radiative_test/"+m+"/eta_2");

    

   }


   exit(-1);

  return;


}
