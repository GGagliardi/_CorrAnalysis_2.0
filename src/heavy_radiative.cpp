#include "../include/heavy_radiative.h"
#include "numerics.h"
#include "stat.h"
#include <boost/optional/detail/optional_reference_spec.hpp>

const double alpha = 1.0 / 137.035999;
const bool UseJack = true;
const int Njacks =  200; //1000;// 93;
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
const double syst_mc_E =0.002;

using namespace std;


void Get_plateaux_int_2pt_H_loc(string Ens, CorrAnalysis &Corr) {

  if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
  else if(Ens=="cD211a.054.96") {Corr.Tmin = 22; Corr.Tmax = 26;}
  else if(Ens=="cE211a.044.112") {Corr.Tmin = 25; Corr.Tmax = 29;}
  else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
  else if(Ens=="cA211a.12.48") {Corr.Tmin = 16; Corr.Tmax = 19;}
  else crash("Ens: "+Ens+" not found");
  
  return;
}

void Get_plateaux_int_2pt_H_sm(string Ens, CorrAnalysis &Corr) {

    if(Ens=="cB211b.072.64") {Corr.Tmin=16; Corr.Tmax=20;}
    else if(Ens=="cD211a.054.96") {Corr.Tmin = 21; Corr.Tmax = 26;}
    else if(Ens=="cE211a.044.112") { Corr.Tmin=25; Corr.Tmax=30;}  //24-30
    else if(Ens=="cC211a.06.80") {Corr.Tmin = 19; Corr.Tmax = 24;}
    else if(Ens=="cA211a.12.48") {Corr.Tmin = 16; Corr.Tmax = 19;}
    else crash("Ens: "+Ens+" not found");

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
    if(M=="M2") {Corr.Tmin=14; Corr.Tmax=18;}
    else if(M=="M3") {Corr.Tmin=10; Corr.Tmax=13;}
    else if(M=="M4") {Corr.Tmin=16; Corr.Tmax=17;}
    else if(M=="M5")  {Corr.Tmin=16; Corr.Tmax=20;}
    else if(M=="M6")  {Corr.Tmin=12; Corr.Tmax=18;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cC211a.06.80") {
    if(M=="M2") {Corr.Tmin=20; Corr.Tmax=25;}
    else if(M=="M3") {Corr.Tmin=19; Corr.Tmax=25;}
    else if(M=="M4") {Corr.Tmin=21; Corr.Tmax=25;}
    else if(M=="M5")  {Corr.Tmin=19; Corr.Tmax=22;}
    else if(M=="M6")  {Corr.Tmin=20; Corr.Tmax=26;}
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
    if(M=="M2") {Corr.Tmin=15; Corr.Tmax=18;}
    else if(M=="M3") {Corr.Tmin=11; Corr.Tmax=13;}
    else if(M=="M4") {Corr.Tmin=13; Corr.Tmax=16;}
    else if(M=="M5")  {Corr.Tmin=12; Corr.Tmax=20;}
    else if(M=="M6")  {Corr.Tmin=15; Corr.Tmax=20;}
    else crash("ERROR IN PLATEAUX");
  }
  else if(Ens=="cC211a.06.80") {
    if(M=="M2") {Corr.Tmin=16; Corr.Tmax=20;}
    else if(M=="M3") {Corr.Tmin=15; Corr.Tmax=20;}
    else if(M=="M4") {Corr.Tmin=18; Corr.Tmax=25;}
    else if(M=="M5")  {Corr.Tmin=21; Corr.Tmax=23;}
    else if(M=="M6")  {Corr.Tmin=20; Corr.Tmax=26;}
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


  Get_chi_c1_decay();


  auto Get_MRS_mass = [](int ir) {

    double L_QCD=Get_Lambda_MS_bar(4);
    double M= pow(lambda,ir)*mc_3GeV_MS; //mc(3GeV) ~ 0.99 GeV FLAG'24
    //solve v= alphas(mv)
    double MRS= MS_bar_to_MRS_mass(3.0, 4, L_QCD, M,  mc_3GeV_MS);
    return MRS;
    
  };

  auto Get_MSMS_mass = [](int ir) {

    double L_QCD=Get_Lambda_MS_bar(4);
    double M= pow(lambda,ir)*mc_3GeV_MS; //mc(3GeV) ~ 0.99 GeV FLAG'24
    //solve v= alphas(mv)
    double MSMS=  m_MS_bar_m( 3.0, 4, L_QCD, M );
    return MSMS;
    
  };

  auto Get_pole_mass = [](int ir) {

    double L_QCD=Get_Lambda_MS_bar(4);
    double M= pow(lambda,ir)*mc_3GeV_MS; //mc(3GeV) ~ 0.99 GeV FLAG'24
    //solve v= alphas(mv)
    double MSMS=  m_MS_bar_m( 3.0, 4, L_QCD, M );
    return MSMS;
    
  };


  cout<<"MRS(M1): "<<2*Get_MRS_mass(0)<<" "<<2*Get_MSMS_mass(0)<<endl;
  cout<<"MRS(M2): "<<2*Get_MRS_mass(1)<<" "<<2*Get_MSMS_mass(1)<<endl;
  cout<<"MRS(M3): "<<2*Get_MRS_mass(2)<<" "<<2*Get_MSMS_mass(2)<<endl;
  cout<<"MRS(M4): "<<2*Get_MRS_mass(3)<<" "<<2*Get_MSMS_mass(3)<<endl;
  cout<<"MRS(M5): "<<2*Get_MRS_mass(4)<<" "<<2*Get_MSMS_mass(4)<<endl;
  cout<<"MRS(M6): "<<2*Get_MRS_mass(5)<<" "<<2*Get_MSMS_mass(5)<<endl;


  auto Get_alpha_vM = [](double M) {


     double L_QCD=Get_Lambda_MS_bar(4);
     //solve v= alphas(0.5*M*v)
     double Mhalf=0.5*M;
     auto lambda_func= [&Mhalf,&L_QCD](double x) {

       return (2.0/3.0)*Get_4l_alpha_s(Mhalf*x, 4 , L_QCD) - x;

     };
     double v=R_brent( lambda_func, 3.0*Get_4l_alpha_s(Mhalf,4,L_QCD), 0.7*Get_4l_alpha_s(Mhalf,4,L_QCD));
    
     return v;  
  };


  

  


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

 
  vector<vector<string>> Ens_list(6);
  vector<string> Ens_list_Hc;

  vector<distr_t_list> H_mass_list;
  for(int i=0;i<6;i++) H_mass_list.emplace_back(UseJack);

  
  
   

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
  
    

    data_t P5P5_sm_loc_M_minus, A0P5_sm_loc_M_minus;
    data_t P5P5_sm_M_minus, P5P5_sm_REST_M_minus;
    
    data_t B1B1_sm_M_minus, B3B3_sm_M_minus, B2B2_sm_M_minus, B2B2_MOT_sm_M_minus, V1V1_sm_M_minus;
    data_t B1B1_sm_loc_M_minus, B2B2_sm_loc_M_minus, B3B3_sm_loc_M_minus;
    data_t PT3_B2P5_tw1_M_minus, PT3_B2P5_tw2_M_minus;
    
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
    PT3_B2P5_tw2_M_minus.Read("../heavy_radiative_tw2", "mes_contr_3PT_"+m_minus+"_tw2", "B1P5", Sort_light_confs);

    B2B2_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B1B1", Sort_light_confs);
    B1B1_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B2B2", Sort_light_confs);
    B3B3_sm_loc_M_minus.Read("../heavy_radiative", "mes_contr_2PT_LOC_"+m_minus, "B3B3", Sort_light_confs);


   

    //#####################################################################################################

    

    
    int Nens=P5P5_sm_loc.size;

    for(int iens=0;iens<Nens;iens++) {

      string Ens =P5P5_sm_loc.Tag[iens];

      bool compute_2nd_tw= (Ens=="cB211b.072.64");

      boost::filesystem::create_directory("../data/heavy_radiative/"+m+"/"+Ens);

      cout<<"########### ANALYZING ENSEMBLE "<<Ens<<" #############"<<endl;


      //########################################## ENSEMBLE INFO ############################################
      int tw1, tw2;
      double amc, amh;
      distr_t a_distr(UseJack);
      distr_t ZV_had(UseJack);
      double syst_mc;
      if(Ens=="cB211b.072.64") {
	tw1=20; tw2=30;
	a_distr=a_B;
	amc= 0.231567;
	ZV_had= ZV_B;
	syst_mc= syst_mc_B;
      }
      else if(Ens=="cA211a.12.48") {
	tw1=19;
	a_distr=a_A;
	amc= 0.262;
	ZV_had= ZV_A;
	syst_mc= syst_mc_A;
      }
      else if(Ens=="cD211a.054.96") {
	tw1=30; //tw2=30;
	a_distr=a_D;
	amc= 0.164898;
	ZV_had= ZV_D;
	syst_mc= syst_mc_D;
      }
      else if(Ens=="cC211a.06.80") {
	tw1=25; //tw2=30;
	a_distr=a_C;
	amc= 0.1984;
	ZV_had=ZV_C;
	syst_mc= syst_mc_C;
      }
      else if(Ens=="cE211a.044.112") {
	tw1=35; //tw2=30;
	a_distr=a_E;
	amc= 0.141255;
	ZV_had=ZV_E;
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
      distr_t_list PT3_tw2_M_minus(UseJack);
      if(compute_2nd_tw) {
	PT3_tw2=Corr.corr_t(PT3_B2P5_tw2.col(0)[0], "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_tw2");
	PT3_tw2_M_minus=Corr.corr_t(PT3_B2P5_tw2_M_minus.col(0)[0], "");
      }
      Corr.Perform_Nt_t_average=1;



  


      
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


           
      //#################################################################################################

      
      //determine DH mass
      Get_plateaux_int_2pt_dH_loc(Ens,m,Corr);
      distr_t DH_sm_loc = Corr.Fit_distr(DH_sm_loc_distr);
      Get_plateaux_int_2pt_dH_sm(Ens,m,Corr);
      distr_t DH_sm = Corr.Fit_distr(DH_distr);
      //Corr.Tmin=12; Corr.Tmax=15;
                 
      double w1= 0.5; // 1/pow(DH_sm.err(),2);
      double w2= 0.5; // 1/pow(DH_sm_loc.err(),2);
      w1= w1/(w1+w2); w2 = 1- w1;
      distr_t DH_mass = w1*DH_sm + w2*DH_sm_loc;
      double syst = fabs(DH_sm.ave() - DH_sm_loc.ave());  //w1*sqrt(pow(DH_sm.ave() - DH_mass.ave(),2)) + w2*sqrt(pow(DH_sm_loc.ave() - DH_mass.ave(),2));
      //syst = sqrt(syst);
      DH_mass= DH_mass.ave() + (DH_mass-DH_mass.ave())*sqrt( 1 + pow( syst/DH_mass.err(),2));
                 
      Get_plateaux_int_2pt_dH_sm(Ens,m,Corr);
     

      //Print DH_mass
      distr_t_list DH_mass_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks,UseJack) ;
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

            
      //   e^-(t-tw1)*MK  * |Z1|/(2MK)
      // C(t) = |Z1|^2/(2MK) e^-Mk*t    sqrt(C(t)) = 

      //determine Hb mass
      Get_plateaux_int_2pt_H_loc(Ens,Corr);
      distr_t H_mass_sm_loc = Corr.Fit_distr(H_mass_sm_loc_averaged_distr);
      Get_plateaux_int_2pt_H_sm(Ens,Corr);
      distr_t H_mass =  Corr.Fit_distr(H_mass_averaged_distr);
      distr_t H_mass_ave= 0.5*(H_mass+H_mass_sm_loc);
      H_mass= H_mass_ave.ave() + (H_mass_ave - H_mass_ave.ave())*sqrt( 1 + pow( 0.5*(H_mass.ave() - H_mass_sm_loc.ave())/H_mass_ave.err(),2));

      
    
      //Print H_mass_fit
      distr_t_list H_mass_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks, UseJack) ;
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	H_mass_fit.distr_list[t] = H_mass;
      }

      Print_To_File({}, { (H_mass_fit/a_distr).ave(), (H_mass_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_fit_pu", "", "");
      Print_To_File({}, { (H_mass_fit).ave(), (H_mass_fit).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/H_mass_fit", "", "");

      distr_t H_mass_minus= H_mass - DH_mass;
      
      distr_t_list MEM_ratio= VEV_ratio*SQRT_D( H_mass_minus/H_mass)*EXPT_D(-1.0*DH_mass, Corr.Nt)*EXP_D(DH_mass*tw1);
      distr_t_list MEM_ratio_tw2(UseJack);
      if(compute_2nd_tw) MEM_ratio_tw2 = VEV_ratio*SQRT_D( H_mass_minus/H_mass)*EXPT_D(-1.0*DH_mass, Corr.Nt)*EXP_D(DH_mass*tw2);
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
      distr_t eta_mass = Corr.Fit_distr(eta_mass_distr);
      distr_t eta_mass_M_minus = Corr.Fit_distr(eta_mass_distr_M_minus);
      distr_t eta_mass_rest= Corr.Fit_distr(eta_mass_rest_distr);
      distr_t eta_mass_M_minus_rest=Corr.Fit_distr(eta_mass_rest_distr_M_minus);
      Print_To_File({}, { ((eta_mass_rest_distr-eta_mass_rest_distr_M_minus)/a_distr).ave(), ((eta_mass_rest_distr-eta_mass_rest_distr_M_minus)/a_distr).err() },  "../data/heavy_radiative/"+m+"/"+Ens+"/Deta_eff_mass_rest_pu", "", "");
      Print_To_File({}, { ((Jpsi_mass_rest_distr- eta_mass_rest)/(a_distr)).ave(), ((Jpsi_mass_rest_distr - eta_mass_rest)/(a_distr)).err() },  "../data/heavy_radiative/"+m+"/"+Ens+"/hyperfine_splitting_pu", "", "");
            

      cout<<"Ens: "<<Ens<<" a: "<<(a_distr/fm_to_inv_Gev).ave()<<" fm"<<endl;
      distr_t mom_k = (2.0/a_distr)*ASIN_D(SQRT_D(SINH_D(eta_mass/2.0)*SINH_D(eta_mass/2.0) -SINH_D(eta_mass_rest/2.0)*SINH_D(eta_mass_rest/2.0) ));
      cout<<"K: "<<mom_k.ave()<<" "<<mom_k.err()<<endl;
      mom_k = eta_mass/a_distr;
      distr_t k_cont= 0.488*a_distr;
      //distr_t C = (POW_D(SINH_D(eta_mass/2.0),2) - POW_D( SIN_D( k_cont/2),2))/POW_D( SINH_D(eta_mass_rest/2.0),2);
      distr_t C = (POW_D(eta_mass,2) - POW_D(k_cont,2))/POW_D(eta_mass_rest,2);
      cout<<"C: "<<C.ave()<<" "<<C.err()<<endl;


      //Print eta_mass_rest_fit
      distr_t_list eta_mass_rest_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks,UseJack) ;
      distr_t_list Deta_mass_rest_fit = 0.0*Get_id_distr_list(Corr.Nt,Njacks,UseJack);
      for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	eta_mass_rest_fit.distr_list[t] = eta_mass_rest;
	Deta_mass_rest_fit.distr_list[t] = eta_mass_rest - eta_mass_M_minus_rest;
      }

      Print_To_File({}, { (eta_mass_rest_fit/a_distr).ave(), (eta_mass_rest_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/eta_mass_rest_fit", "", "");
      Print_To_File({}, { (Deta_mass_rest_fit/a_distr).ave(), (Deta_mass_rest_fit/a_distr).err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/Deta_mass_rest_fit", "", "");

      
      distr_t_list VV_ETA= Corr.matrix_element_t(P5P5_sm_distr, "");
      distr_t ETA_MEM_tw1 = EXP_D(-1.0*eta_mass*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr, ""))/(2.0*eta_mass);
      distr_t ETA_MEM_tw1_M_minus = EXP_D(-1.0*eta_mass_M_minus*tw1)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr_M_minus, ""))/(2.0*eta_mass_M_minus);


            
      
      distr_t ETA_MEM_tw2(UseJack), ETA_MEM_tw2_M_minus(UseJack);
      if(compute_2nd_tw) {
	ETA_MEM_tw2 = EXP_D(-1.0*eta_mass*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr, ""))/(2.0*eta_mass);
	ETA_MEM_tw2_M_minus = EXP_D(-1.0*eta_mass_M_minus*tw2)*Corr.Fit_distr(Corr.matrix_element_t(P5P5_sm_distr_M_minus, ""))/(2.0*eta_mass_M_minus);
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
      distr_t_list ZV_distr_M_minus(UseJack);
           
      if(mm != "M1") {
	ZV_distr_M_minus= 2*(amh/lambda)*P5P5_sm_loc_distr_M_minus/(distr_t_list::derivative(A0P5_sm_loc_distr_M_minus, 0));
      }
      else {  ZV_distr_M_minus= ZV_distr;}
      Print_To_File({}, {ZV_distr.ave(), ZV_distr.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/ZV", "", "");
      if(Ens=="cB211b.072.64") {Corr.Tmin=30; Corr.Tmax=50; }
      else if(Ens=="cA211a.12.48") {Corr.Tmin=26; Corr.Tmax=38; }
      else if(Ens=="cD211a.054.96") {Corr.Tmin = 35; Corr.Tmax = 60;}
      else if(Ens=="cC211a.06.80") {Corr.Tmin = 40; Corr.Tmax = 55;}
      else if(Ens=="cE211a.044.112") {Corr.Tmin = 45; Corr.Tmax = 66;}
      else crash("Ens "+Ens+" not found");
      distr_t ZV= Corr.Fit_distr(ZV_distr);
      distr_t ZV_M_minus= Corr.Fit_distr(ZV_distr_M_minus);
      if(m=="Hc") ZV_M1.distr_list.push_back(ZV);
      cout<<"eta_fact_ratio: "<<(ETA_MEM_tw1/ETA_MEM_tw1_M_minus).ave()<<" +- "<<(ETA_MEM_tw1/ETA_MEM_tw1_M_minus).err()<<endl;
      cout<<"Z[m]/Z[mh-1]: "<<(ZV/ZV_M_minus).ave()<<" +- "<<(ZV/ZV_M_minus).err()<<endl;
      //################################################################################################

      
      //#####################       EVALUATE FORM FACTOR AND RATIOS  ###################################
      distr_t_list PT3_tw1_ratio = PT3_tw1/PT3_tw1_M_minus;
      distr_t_list PT3_tw2_ratio(UseJack);
      if(compute_2nd_tw) PT3_tw2_ratio = PT3_tw2/PT3_tw2_M_minus;
      Print_To_File({}, { PT3_tw1_ratio.ave(), PT3_tw1_ratio.err()}, "../data/heavy_radiative/"+m+"/"+Ens+"/PT3_ratio_tw1", "", "");
      distr_t_list F0_tw1 = 2.0*ZV*PT3_tw1/(ETA_MEM_tw1*H_MEM_tw1*H_mass);
      distr_t_list F0_tw1_II = 2.0*ZV_had*(ZV/ZV_M1[iens])*PT3_tw1/(ETA_MEM_tw1*H_MEM_tw1*H_mass);
      //distr_t SCALE_FACT = SQRT_D( eta_mass*H_mass/(eta_mass_M_minus*H_mass_minus));
      distr_t SCALE_FACT= SQRT_D(H_mass/H_mass_minus);

           
      if(TYPE=="WEAK") SCALE_FACT= 1.0*Get_id_distr(Njacks,UseJack);  // SCALE_FACT= 1.0/(eta_mass_rest/eta_mass_M_minus_rest);  // SCALE_FACT=1.0*Get_id_jack_distr(Njacks);
      SCALE_FACT= 1.0*Get_id_distr(Njacks,UseJack);
      distr_t_list ratio_corr= ((ZV/ZV_M_minus)*PT3_tw1_ratio/( (ETA_MEM_tw1/ETA_MEM_tw1_M_minus)*(H_mass/H_mass_minus)*MEM_ratio))/SCALE_FACT;
        
      distr_t_list F0_tw2(UseJack);
      distr_t_list F0_tw2_II(UseJack);
      distr_t_list ratio_corr_tw2(UseJack);
           
      if(compute_2nd_tw) {
	F0_tw2 = 2.0*ZV*PT3_tw2/(ETA_MEM_tw2*H_MEM_tw2*H_mass);
	F0_tw2_II = 2.0*ZV_had*(ZV/ZV_M1[iens])*PT3_tw2/(ETA_MEM_tw2*H_MEM_tw2*H_mass);
	ratio_corr_tw2=  ((ZV/ZV_M_minus)*PT3_tw2_ratio/( (ETA_MEM_tw2/ETA_MEM_tw2_M_minus)*(H_mass/H_mass_minus)*MEM_ratio_tw2))/SCALE_FACT;

      }
      //################################################################################################


      

      //################################# PRINT FORM FACTORS ###########################################
      //reduce
      distr_t_list F0_tw1_RED(UseJack);
      distr_t_list F0_tw2_RED(UseJack);
      distr_t_list F0_tw1_II_RED(UseJack);
      distr_t_list F0_tw2_II_RED(UseJack);
      distr_t_list ratio_corr_red(UseJack);
      distr_t_list ratio_corr_red_tw2(UseJack);
 
      
      for(int t=tw1;t<Corr.Nt;t++)  {	F0_tw1_RED.distr_list.push_back( F0_tw1[t]);  F0_tw1_II_RED.distr_list.push_back( F0_tw1_II[t]);  ratio_corr_red.distr_list.push_back( ratio_corr[t]);}
      if(compute_2nd_tw) {
	for(int t=tw2;t<Corr.Nt;t++) { F0_tw2_RED.distr_list.push_back( F0_tw2[t]);  F0_tw2_II_RED.distr_list.push_back( F0_tw2_II[t]); ratio_corr_red_tw2.distr_list.push_back( ratio_corr_tw2[t]); }
      }
      Print_To_File({}, {ratio_corr_red.ave(), ratio_corr_red.err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/F_ratio", "", "");
      if(compute_2nd_tw)    Print_To_File({}, {ratio_corr_red_tw2.ave(), ratio_corr_red_tw2.err()},  "../data/heavy_radiative/"+m+"/"+Ens+"/F_ratio_tw2", "", ""); 
          
      Print_To_File({}, {F0_tw1_RED.ave(), F0_tw1_RED.err(), F0_tw1_II_RED.ave(), F0_tw1_II_RED.err() }, "../data/heavy_radiative/"+m+"/"+Ens+"/F0_tw1", "", "");
      if(compute_2nd_tw) Print_To_File({}, {F0_tw2_RED.ave(), F0_tw2_RED.err(), F0_tw2_II_RED.ave(), F0_tw2_II_RED.err() }, "../data/heavy_radiative/"+m+"/"+Ens+"/F0_tw2", "", "");
      //################################################################################################


      //fit Fc

      Corr.Tmin= (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
      Corr.Tmax= (int)(1.4/(a_distr.ave()/fm_to_inv_Gev));

      distr_t F1=Corr.Fit_distr(F0_tw1_II_RED);

      
      //push back results for Fc
      if(m=="Hc") {
	Fc.distr_list.push_back(F1);
	metac.distr_list.push_back( eta_mass_rest/a_distr);
	mhc.distr_list.push_back( H_mass/a_distr);
	a_distr_list.distr_list.push_back(a_distr);
	Ens_list_Hc.push_back(Ens);
      }

      //push back result for Fm1
      if(m=="M1") {
	Fm1.distr_list.push_back(F1);
      }

    
      distr_t to_m_etah= eta_mass_rest/a_distr;

      //to_m_etah = to_m_etah.ave() + (to_m_etah- to_m_etah.ave())*sqrt( 1 + pow(to_m_etah.ave()*syst_mc/to_m_etah.err(),2));
      
      if(m != "Hc") {

	m_etah_list_glb[im-1].distr_list.push_back(to_m_etah) ;
	aH_glb_list[im-1].distr_list.push_back( a_distr);
      }

      if(m=="M1") {
	mom_k_list[0].distr_list.push_back(mom_k);
	m_etah_list[0].distr_list.push_back(to_m_etah);
	aH_list[0].distr_list.push_back( a_distr);
	C_beta[0].distr_list.push_back( C);
	DH_mass_list[0].distr_list.push_back( DH_mass/a_distr);
	
	Ens_list[0].push_back(Ens);
	
      }

      if(m=="M2") {
	C_beta[1].distr_list.push_back( C);
	Corr.Tmin = (int)(0.6/(a_distr.ave()/fm_to_inv_Gev));
	Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	distr_t R1= Corr.Fit_distr(ratio_corr_red);
	R[0].distr_list.push_back(R1);
	aR[0].distr_list.push_back(a_distr);

	mom_k_list[1].distr_list.push_back(mom_k);
	m_etah_list[1].distr_list.push_back(to_m_etah);
	aH_list[1].distr_list.push_back( a_distr);

	DH_mass_list[1].distr_list.push_back( DH_mass/a_distr);
	
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
	distr_t R2 = Corr.Fit_distr(ratio_corr_red);
	R[1].distr_list.push_back(R2);
	aR[1].distr_list.push_back(a_distr);

	mom_k_list[2].distr_list.push_back(mom_k);
	m_etah_list[2].distr_list.push_back(to_m_etah);
	aH_list[2].distr_list.push_back( a_distr);

	DH_mass_list[2].distr_list.push_back( DH_mass/a_distr);

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
	distr_t R3 = Corr.Fit_distr(ratio_corr_red);
	R[2].distr_list.push_back(R3);
	aR[2].distr_list.push_back(a_distr);

	mom_k_list[3].distr_list.push_back(mom_k);
	m_etah_list[3].distr_list.push_back(to_m_etah);
	aH_list[3].distr_list.push_back( a_distr);

	DH_mass_list[3].distr_list.push_back( DH_mass/a_distr);

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

	  distr_t R4 = Corr.Fit_distr(ratio_corr_red);
	  R[3].distr_list.push_back(R4);
	  aR[3].distr_list.push_back(a_distr);

	  mom_k_list[4].distr_list.push_back(mom_k);
	  m_etah_list[4].distr_list.push_back(to_m_etah);
	  aH_list[4].distr_list.push_back( a_distr);

	  DH_mass_list[4].distr_list.push_back( DH_mass/a_distr);
	  Ens_list[4].push_back(Ens);
	  

	}
      }
      else if(m=="M6") {
	C_beta[5].distr_list.push_back( C);
	//if(Ens != "cB211b.072.64" && Ens != "cC211a.06.80" && Ens != "cA211a.12.48") {
	if(Ens != "cB211b.072.64" &&  Ens != "cA211a.12.48") {
	  if(Ens=="cD211a.054.96") {
	    Corr.Tmin = (int)(0.6/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else if(Ens=="cE211a.044.112") {
	    Corr.Tmin = (int)(0.65/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.0/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else if(Ens =="cC211a.06.80") {
	    Corr.Tmin = (int)(0.9/(a_distr.ave()/fm_to_inv_Gev));
	    Corr.Tmax = (int)(1.5/(a_distr.ave()/fm_to_inv_Gev));
	  }
	  else crash("Ens: "+Ens+" not found");

	  
	  distr_t R5 = Corr.Fit_distr(ratio_corr_red);
	  R[4].distr_list.push_back(R5);
	  aR[4].distr_list.push_back(a_distr);

	  mom_k_list[5].distr_list.push_back(mom_k);
	  m_etah_list[5].distr_list.push_back(to_m_etah);
	  aH_list[5].distr_list.push_back( a_distr);

	  DH_mass_list[5].distr_list.push_back( DH_mass/a_distr);
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
  

  
  //GET H_mass_list
  for(int ir=0;ir<6;ir++) {
    
    H_mass_list[ir] = 0.0*m_etah_list[ir];
    
    for(int iens=0;iens<(signed)Ens_list[ir].size(); iens++) {

      //find mHc
      int ens_id_hc=-1;
      for(int s=0;s<(signed)Ens_list_Hc.size(); s++) {
	if(Ens_list_Hc[s] == Ens_list[ir][iens]) ens_id_hc=s;
      }
      assert(ens_id_hc != -1);
      H_mass_list[ir].distr_list[iens] = mhc.distr_list[ens_id_hc];
      
      
      for(int is=1;is<=ir;is++) {
	//find ensemble corresponding to Ensemble
	int ens_id=-1;
	for(int is_ens=0;is_ens<(signed)Ens_list[is].size(); is_ens++) {
	  if(Ens_list[is][is_ens] == Ens_list[ir][iens] ) ens_id=is_ens;
	}
	assert(ens_id != -1);

	H_mass_list[ir].distr_list[iens] = H_mass_list[ir].distr_list[iens] + DH_mass_list[is].distr_list[ens_id];
      }
    }
  }

  

  for(int ir=0;ir<6;ir++) {

    Print_To_File({}, {aH_list[ir].ave(), (H_mass_list[ir]-m_etah_list[ir]).ave(), (H_mass_list[ir]-m_etah_list[ir]).err()}, "../data/heavy_radiative/Fits_R/mass_splitting_r_"+to_string(ir), "", "");

  }

  

  vector<distr_t_list> P_AT_splitting(6);
    for(int ir=0;ir<6;ir++) {
      P_AT_splitting[ir] = H_mass_list[ir] - m_etah_list[ir];
    }
  
      
  
 
  //print C_beta

  for(int r=0;r<(signed)C_beta.size();r++) {
    Print_To_File({}, {a_distr_list.ave(), C_beta[r].ave(), C_beta[r].err()} ,  "../data/heavy_radiative/Fits_R/C_beta_"+to_string(r+1), "", "");
  }




  //continuum limit extrapolation of F1_c F1_M1 mc and mhc

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


  bf_mh.Fix_par("A", 0.0);
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

  bf_mh.Release_par("A2");
  bf_mh_ch2.Release_par("A2");
  bf_mh.Release_par("A");
  bf_mh_ch2.Release_par("A");
  
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
  distr_t M_Hc(UseJack), D_Hc(UseJack), M_Hc_4(UseJack), D_Hc_4(UseJack), D2_Hc_4(UseJack);
  distr_t M_etac(UseJack),  D_etac(UseJack), M_etac_4(UseJack), D_etac_4(UseJack), D2_etac_4(UseJack);
  distr_t FFc(UseJack), FFc_const(UseJack), D_FFc(UseJack);
  distr_t FFm1(UseJack), FFm1_const(UseJack), D_FFm1(UseJack);

  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) {

    //Hc
    M_Hc.distr.push_back( Bt_fit_Hc.par[ijack].R);
    M_Hc_4.distr.push_back( Bt_fit_Hc_const.par[ijack].R);
    D_Hc.distr.push_back( Bt_fit_Hc.par[ijack].A);
    D_Hc_4.distr.push_back( Bt_fit_Hc_const.par[ijack].A);
    D2_Hc_4.distr.push_back( Bt_fit_Hc_const.par[ijack].A2);
    //etac
    M_etac.distr.push_back( Bt_fit_etac.par[ijack].R);
    M_etac_4.distr.push_back( Bt_fit_etac_const.par[ijack].R);
    D_etac.distr.push_back( Bt_fit_etac.par[ijack].A);
    D_etac_4.distr.push_back( Bt_fit_etac_const.par[ijack].A);
    D2_etac_4.distr.push_back( Bt_fit_etac_const.par[ijack].A2);
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
  int dof_const = Nmeas-3;
  

  double w_etac1 = exp(-0.5*(ch2_etac_lin - 2*dof_lin))/( exp(-0.5*(ch2_etac_lin - 2*dof_lin)) + exp(-0.5*(ch2_etac_const - 2*dof_const)));
  double w_etac2 = 1.0 - w_etac1;

  double w_Hc1 = exp(-0.5*(ch2_Hc_lin - 2*dof_lin))/( exp(-0.5*(ch2_Hc_lin - 2*dof_lin)) + exp(-0.5*(ch2_Hc_const - 2*dof_const)));
  double w_Hc2 = 1.0 - w_Hc1;

  double w_Fc1 = exp(-0.5*(ch2_Fc_lin - 2*dof_lin))/( exp(-0.5*(ch2_Fc_lin - 2*dof_lin)) + exp(-0.5*(ch2_Fc_const - 2*dof_const)));
  double w_Fc2 = 1.0 - w_Fc1;

  double w_Fm11 = exp(-0.5*(ch2_Fm1_lin - 2*dof_lin))/( exp(-0.5*(ch2_Fm1_lin - 2*dof_lin)) + exp(-0.5*(ch2_Fm1_const - 2*dof_const)));
  double w_Fm12 = 1.0 - w_Fm11;


  //how do we fit

  //w_etac1=1.0; w_etac2=0.0;
  //w_Hc1=1.0; w_Hc2=0.0;
  w_Fc1=1.0; w_Fc2=0.0;
  w_Fm11=1.0; w_Fm12=0.0;

  cout<<"reduced ch2 (etac): "<<ch2_etac_lin/dof_lin<<endl;
  cout<<"reduced ch2 (hc): "<<ch2_Hc_lin/dof_lin<<endl;
  cout<<"reduced ch2 (Fc): "<<ch2_Fc_lin/dof_lin<<endl;
  cout<<"reduced ch2 (Fm1): "<<ch2_Fm1_lin/dof_lin<<endl;

  distr_t M_Hc_ave= w_Hc1*M_Hc + w_Hc2*M_Hc_4;
  distr_t M_etac_ave= w_etac1*M_etac + w_etac2*M_etac_4;
  distr_t Fc_ave = w_Fc1*FFc + w_Fc2*FFc_const;
  distr_t Fm1_ave = w_Fm11*FFm1 + w_Fm12*FFm1_const;
  
  distr_t D_Hc_ave = w_Hc1*D_Hc + w_Hc2*D_Hc_4;
  distr_t D_etac_ave = w_etac1*D_etac + w_etac2*D_etac_4;
  distr_t D_FFc_ave= w_Fc1*D_FFc;
  distr_t D_FFm1_ave= w_Fm11*D_FFm1;

  distr_t D2_Hc_ave = w_Hc2*D2_Hc_4;   
  distr_t D2_etac_ave = w_etac2*D2_etac_4;
  
  //now combine and print the results


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list Metac_to_print(UseJack);
  distr_t_list MHc_to_print(UseJack);
  distr_t_list Fc_to_print(UseJack);
  distr_t_list Fm1_to_print(UseJack);
  Vfloat  a_to_print;

  //print in step of 0.01 fm;
  double a_max = 0.1;
  
  for(int i=0;i< (int)(a_max/0.001);i++) {
    double x = i*0.001;
    Metac_to_print.distr_list.push_back(  M_etac_ave + D_etac_ave*pow(x*fm_to_inv_Gev,2) + D2_etac_ave*pow(x*fm_to_inv_Gev,4) );
    MHc_to_print.distr_list.push_back(  M_Hc_ave + D_Hc_ave*pow(x*fm_to_inv_Gev,2) + D2_Hc_ave*pow(x*fm_to_inv_Gev,4) );
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


  cout<<"continuum limit extrapolation of P_AT_splitting and m_etah"<<endl;

  distr_t_list  extrapolated_m_etah_SF(UseJack), extrapolated_m_Hh_SF(UseJack), extrapolated_ratio_metah(UseJack), extrapolated_ratio_mHh(UseJack);


  vector<distr_t_list> ratio_m_etah_list = m_etah_list;
  vector<distr_t_list> ratio_m_Hh_list  =  H_mass_list;

  for(int ir=0;ir<(signed)m_etah_list.size();ir++) {

    if(ir==0)  { ratio_m_etah_list[ir] = ratio_m_etah_list[ir]/ratio_m_etah_list[ir];  ratio_m_Hh_list[ir] = ratio_m_Hh_list[ir]/ratio_m_Hh_list[ir]; }
    else {
      for(int ins=0;ins<(signed)m_etah_list[ir].size() ; ins++) {
	string Ens_A = Ens_list[ir][ins];
	bool fnd=false;
	for(int jns=0;jns<(signed)m_etah_list[ir-1].size(); jns++) { 
	  if(Ens_list[ir-1][jns] == Ens_A) { fnd=true;
	    ratio_m_etah_list[ir].distr_list[ins] = ratio_m_etah_list[ir].distr_list[ins]/m_etah_list[ir-1][jns] ;
	    ratio_m_Hh_list[ir].distr_list[ins] = ratio_m_Hh_list[ir].distr_list[ins]/H_mass_list[ir-1][jns] ;
	  }
	}
	assert(fnd==true);
      }

    }
  }


  //assign fiducial error of 0.02% to the ratios

  for(int ir=0;ir<(signed)m_etah_list.size(); ir++) {
    for(int ins=0;ins<(signed)ratio_m_etah_list[ir].size();ins++) {
      distr_t A= ratio_m_etah_list[ir].distr_list[ins];
      distr_t B= ratio_m_Hh_list[ir].distr_list[ins];
      A= A.ave() + (A-A.ave())*sqrt( 1.0 + pow(0.01*0.015*A.ave()/A.err(),2));
      B= B.ave() + (B-B.ave())*sqrt( 1.0 + pow(0.01*0.015*B.ave()/B.err(),2));
      ratio_m_etah_list[ir].distr_list[ins] = A;
      ratio_m_Hh_list[ir].distr_list[ins]=B;
    }
  }

  
  
    
     

  for(int b=0;b<5;b++) {

    for(int ir=0;ir<(signed)m_etah_list.size();ir++) {

      if(b==0) cout<<"Extrapolation of metah: "<<ir+1<<endl;
      else if(b==1) cout<<"Extrapolation of P_AT_splitting: "<<ir+1<<endl;
      else if(b==2)  cout<<"Extrapolation of mhh: "<<ir+1<<endl;
      else if(b==3) cout<<"Extrapolation of meta(n)/metah(n-1): "<<ir+1<<endl;
      else  cout<<"Extrapolation of mH(n)/mH(n-1): "<<ir+1<<endl;
      
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

      bf_r.Fix_par("A2",0.0);
      bf_r_ch2.Fix_par("A2",0.0);

     
       
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

	 distr_t_list Dat(UseJack);
	 if(b==0) Dat= m_etah_list[ir];
	 else if(b==1) Dat= P_AT_splitting[ir];
	 else if(b==2)  Dat= H_mass_list[ir];
	 else if(b==3)  Dat= ratio_m_etah_list[ir];
	 else Dat= ratio_m_Hh_list[ir];
	 
	 data[ijack][iens].R= Dat.distr_list[iens].distr[ijack]; 
	 data[ijack][iens].R_err = Dat.err(iens); 
	 data[ijack][iens].a = aH_list[ir].distr_list[iens].distr[ijack];
	 //mean values
	 if(ijack==0) {
	   data_ch2[ijack][iens].R= Dat.ave(iens) ;
	   data_ch2[ijack][iens].R_err = Dat.err(iens); 
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

      boot_fit_data<fpar_r> Bt_fit_R_a4;
      boot_fit_data<fpar_r> Bt_fit_R_ch2_a4;

      if(Nmeas_rat >= 4 ) {
	bf_r.Release_par("A2");
	bf_r_ch2.Release_par("A2");
      }
      
      Bt_fit_R_a4= bf_r.Perform_bootstrap_fit();
      Bt_fit_R_ch2_a4= bf_r_ch2.Perform_bootstrap_fit();

      distr_t r_a2(UseJack),  D_r_a2(UseJack), D2_r_a2(UseJack);

      distr_t r_a4(UseJack), D_r_a4(UseJack), D2_r_a4(UseJack);  
     
       //retrieve parameters
       for(int ijack=0;ijack<Njacks;ijack++) {
         
	 r_a2.distr.push_back( Bt_fit_R.par[ijack].R);
	 D_r_a2.distr.push_back( Bt_fit_R.par[ijack].A);
	 D2_r_a2.distr.push_back( Bt_fit_R.par[ijack].A2);

	 r_a4.distr.push_back( Bt_fit_R_a4.par[ijack].R);
	 D_r_a4.distr.push_back( Bt_fit_R_a4.par[ijack].A);
	 D2_r_a4.distr.push_back( Bt_fit_R_a4.par[ijack].A2);
	 
	 
       }
       //get ch2
       double ch2_R_lin = Bt_fit_R_ch2.get_ch2_ave();
       double ch2_R_quad = Bt_fit_R_ch2_a4.get_ch2_ave();

       int Npars_lin= 2;
       int Npars_quad= (Nmeas_rat >= 4)?3:2;

       double w_lin= exp(-0.5*(ch2_R_lin +2*Npars_lin))/(  exp(-0.5*(ch2_R_lin + 2*Npars_lin)) + exp(-0.5*(ch2_R_quad + 2*Npars_quad)));
       double w_quad = 1.0 -w_lin;


       cout<<"Printing ch2 for b: "<<b<<" ir: "<<ir<<endl;
       cout<<"ch2(lin): "<<ch2_R_lin/(Nmeas_rat - Npars_lin)<<endl;
       cout<<"ch2(quad): "<<ch2_R_quad/(Nmeas_rat - Npars_quad)<<endl;


       distr_t r(UseJack), D_r(UseJack), D2_r(UseJack);

       r = w_lin*r_a2 + w_quad*r_a4;
       D_r = w_lin*D_r_a2 + w_quad*D_r_a4;
       D2_r = w_lin*D2_r_a2 + w_quad*D2_r_a4;

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

       string out_string;
       if(b==0) out_string = "metah_"+to_string(ir);
       else if(b==1) out_string = "P_AT_splitting_"+to_string(ir);
       else if(b==2)  out_string = "mHh_"+to_string(ir);
       else if(b==3) out_string = "ratio_metah_"+to_string(ir);
       else out_string = "ratio_Hh_"+to_string(ir);

            
       Print_To_File({}, {a_to_print, (R_to_print).ave(), (R_to_print).err()} , "../data/heavy_radiative/Fits_R/"+out_string+"_extr.fit", "", "");

       if(b==0) {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (m_etah_list[ir]).ave(), (m_etah_list[ir]).err() }, "../data/heavy_radiative/Fits_R/"+out_string+".data", "", "");
	 extrapolated_m_etah_SF.distr_list.push_back(r);
       }
       else if(b==1) {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (P_AT_splitting[ir]).ave(), (P_AT_splitting[ir]).err() }, "../data/heavy_radiative/Fits_R/"+out_string+".data", "", "");
       }
       else if(b==2) {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (H_mass_list[ir]).ave(), (H_mass_list[ir]).err() }, "../data/heavy_radiative/Fits_R/"+out_string+".data", "", "");
	 extrapolated_m_Hh_SF.distr_list.push_back(r);
       }
       else if(b==3) {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (ratio_m_etah_list[ir]).ave(), (ratio_m_etah_list[ir]).err() }, "../data/heavy_radiative/Fits_R/"+out_string+".data", "", "");
	 extrapolated_ratio_metah.distr_list.push_back(r);
       }
       else  {
	 Print_To_File({} , { (aH_list[ir]/fm_to_inv_Gev).ave(), (ratio_m_Hh_list[ir]).ave(), (ratio_m_Hh_list[ir]).err() }, "../data/heavy_radiative/Fits_R/"+out_string+".data", "", "");
	 extrapolated_ratio_mHh.distr_list.push_back(r);
       }

    }
  }

  //reconstruct m_etah from ratios


  distr_t_list extrapolated_m_etah_FR=extrapolated_m_etah_SF;
  distr_t_list extrapolated_m_Hh_FR= extrapolated_m_Hh_SF;
  for(int ir=1; ir<(signed)extrapolated_m_etah_FR.size();ir++) {
    extrapolated_m_etah_FR.distr_list[ir] = extrapolated_m_etah_FR[ir-1]*extrapolated_ratio_metah.distr_list[ir];
    extrapolated_m_Hh_FR.distr_list[ir] = extrapolated_m_Hh_FR[ir-1]*extrapolated_ratio_mHh.distr_list[ir];
  }


  for(int a=0;a<(signed)m_etah_list_glb.size();a++) {
    Print_To_File({}, { (aH_glb_list[a]/fm_to_inv_Gev).ave(), m_etah_list_glb[a].ave(), m_etah_list_glb[a].err()}, "../data/heavy_radiative/Fits_R/metah_glb_"+to_string(a)+".data", "","");
  }

  
   //global mass-cont extrapolation for metah and Mh

  //#######################################################################################################


  distr_t_list extrapolated_m_etah(UseJack), extrapolated_m_Hh(UseJack);
  distr_t_list HQV_m_etah(UseJack), HQV_m_Hh(UseJack);


  for(int b=0;b<2;b++) {

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
    bf_glb_mass_ch2.Add_par("R1", 3.0, 0.1);
    bf_glb_mass_ch2.Add_par("R2", 3.5, 0.1);
    bf_glb_mass_ch2.Add_par("R3", 4.0, 0.1);
    bf_glb_mass_ch2.Add_par("R4", 4.5, 0.1);
    bf_glb_mass_ch2.Add_par("R5", 5.0, 0.1);
    bf_glb_mass_ch2.Add_par("R6", 5.5, 0.1);  
    bf_glb_mass_ch2.Add_par("A", 1.0, 0.1);
    bf_glb_mass_ch2.Add_par("Am", 1.0, 0.1);
    bf_glb_mass_ch2.Add_par("Am2", 1.0, 0.1);
    bf_glb_mass_ch2.Add_par("A4", 1.0, 0.1);
    bf_glb_mass_ch2.Add_par("A4m2", 1.0, 0.1);
    bf_glb_mass_ch2.Add_par("A4m4", 1.0, 0.1);

    //##############################//

    bf_glb_mass.Fix_par("Am",0.0);
    bf_glb_mass_ch2.Fix_par("Am",0.0);
    bf_glb_mass.Fix_par("A4m4",0.0);
    bf_glb_mass_ch2.Fix_par("A4m4",0.0);



    //ansatz
    bf_glb_mass.ansatz=  [](const fpar_r_glb_mass &p, const ipar_r_glb_mass &ip) {

      double x=0.0;
      
      if(ip.n==1) { x = p.R1;}
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

	distr_t_list Dat(UseJack);
	if(b==0) Dat = m_etah_list[rr];
	else Dat = H_mass_list[rr];
	
	for(int iens=0;iens<(signed)Dat.size();iens++) {
	  
	  data_glb_mass[ijack][it].R= Dat.distr_list[iens].distr[ijack]; 
	  data_glb_mass[ijack][it].R_err = Dat.err(iens);
	  data_glb_mass[ijack][it].m_r = pow(lambda, rr);
	  data_glb_mass[ijack][it].a = aH_list[rr].distr_list[iens].distr[ijack];
	  data_glb_mass[ijack][it].n = rr+1;
	  //mean values
	  if(ijack==0) {
	    data_glb_mass_ch2[ijack][it].R= Dat.ave(iens); 
	    data_glb_mass_ch2[ijack][it].R_err = Dat.err(iens);
	    data_glb_mass_ch2[ijack][it].m_r = pow(lambda, rr);
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


    distr_t pR1(UseJack), pR2(UseJack), pR3(UseJack), pR4(UseJack), pR5(UseJack), pR6(UseJack), pHA(UseJack), pHAm(UseJack),  pHAm2(UseJack), pHA4(UseJack), pHA4m2(UseJack), pHA4m4(UseJack);

    
     
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
    vector<distr_t> pR({pR1,pR2,pR3,pR4,pR5,pR6});

    cout<<"ch2/dof [glb_mass fit, b="<<b<<"]: "<<ch2_glb_mass/(Nmeas_glb_mass - bf_glb_mass.Get_number_of_fit_pars())<<endl;



    //print fit function at non-zero lattice spacing for each mass

    for(int rr=0; rr<(signed)m_etah_list.size();rr++) {

      //print fitting functions for MP, FP, MP_ov_FP
      distr_t_list func_to_print(UseJack);
      Vfloat  a_to_print;
      //print in step of 0.01 fm;
      double a_max = 0.1;

      double mr= pow(lambda, rr);
     
      for(int i=0;i< (int)(a_max/0.001);i++) {
	double x = i*0.001;
	double a= x*fm_to_inv_Gev;
	distr_t cont(UseJack);
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

      string out_string;
      if(b==0) out_string= "m_etah_"+to_string(rr);
      else out_string = "m_Hh_"+to_string(rr);
      Print_To_File({}, {a_to_print, (func_to_print).ave(), (func_to_print).err()} , "../data/heavy_radiative/Fits_R/"+out_string+"_extr_glb_mass.fit", "", "");
    }

    if(b==0) {
      for(int i=0;i<6;i++) {
	extrapolated_m_etah.distr_list.push_back( pR[i]);
	HQV_m_etah.distr_list.push_back( sqrt(2.0)*SQRT_D(pR[i]/Get_MRS_mass(i) -2));   
      }
    }
    else {
      for(int i=0;i<6;i++) {
	extrapolated_m_Hh.distr_list.push_back( pR[i]);
	HQV_m_Hh.distr_list.push_back( sqrt(2.0)*SQRT_D( pR[i]/Get_MRS_mass(i) -2));   
      }
    }
    
    
  }

  //print extrapolated masses
  Print_To_File({}, { extrapolated_m_etah.ave(), extrapolated_m_etah.err() , HQV_m_etah.ave(), HQV_m_etah.err() }, "../data/heavy_radiative/Fits_R/m_etah_extr_glb_mass.res","","");
  Print_To_File({}, { extrapolated_m_Hh.ave(), extrapolated_m_Hh.err() , HQV_m_Hh.ave(), HQV_m_Hh.err() }, "../data/heavy_radiative/Fits_R/m_Hh_extr_glb_mass.res","","");
     
  //#######################################################################################################


  //combine extrapolated_* and extrapolate_SF_*


  distr_t_list extrapolated_m_etah_AVE(UseJack), extrapolated_m_Hh_AVE(UseJack);

  for(int i=0;i<(signed)extrapolated_m_etah.size();i++) {

    extrapolated_m_etah_AVE.distr_list.push_back( 0.5*(extrapolated_m_etah.distr_list[i] + extrapolated_m_etah_SF.distr_list[i]));
    double syst_eta= 0.5*pow( extrapolated_m_etah_AVE.ave(i) - extrapolated_m_etah.ave(i),2) +  0.5*pow( extrapolated_m_etah_AVE.ave(i) - extrapolated_m_etah_SF.ave(i),2);
    syst_eta=sqrt(syst_eta);
    extrapolated_m_etah_AVE.distr_list[i] = extrapolated_m_etah_AVE.distr_list[i].ave() + ( extrapolated_m_etah_AVE.distr_list[i] -  extrapolated_m_etah_AVE.distr_list[i].ave())*sqrt( 1.0 + pow(syst_eta/extrapolated_m_etah_AVE.err(i),2));


    extrapolated_m_Hh_AVE.distr_list.push_back( 0.5*(extrapolated_m_Hh.distr_list[i] + extrapolated_m_Hh_SF.distr_list[i]));
    double syst_H= 0.5*pow( extrapolated_m_Hh_AVE.ave(i) - extrapolated_m_Hh.ave(i),2) +  0.5*pow( extrapolated_m_Hh_AVE.ave(i) - extrapolated_m_Hh_SF.ave(i),2);
    syst_H=sqrt(syst_H);
    extrapolated_m_Hh_AVE.distr_list[i] = extrapolated_m_Hh_AVE.distr_list[i].ave() + ( extrapolated_m_Hh_AVE.distr_list[i] -  extrapolated_m_Hh_AVE.distr_list[i].ave())*sqrt( 1.0 + pow(syst_H/extrapolated_m_Hh_AVE.err(i),2));
    
    

  }
  

  Print_To_File({}, { extrapolated_m_etah_AVE.ave(), extrapolated_m_etah_AVE.err(), extrapolated_m_etah_SF.ave(), extrapolated_m_etah_SF.err(), extrapolated_m_etah.ave(), extrapolated_m_etah.err(), extrapolated_m_etah_FR.ave(), extrapolated_m_etah_FR.err()}, "../data/heavy_radiative/Fits_R/masses_etah.fin","","");
  Print_To_File({}, { extrapolated_m_Hh_AVE.ave(), extrapolated_m_Hh_AVE.err(), extrapolated_m_Hh_SF.ave(), extrapolated_m_Hh_SF.err(), extrapolated_m_Hh.ave(), extrapolated_m_Hh.err(), extrapolated_m_Hh_FR.ave(), extrapolated_m_Hh_FR.err()}, "../data/heavy_radiative/Fits_R/masses_Hh.fin","","");



  
  //now extrapolate to the continuum limit the ratios

  cout<<"Continuum extrapolation of ratios..."<<endl;


  distr_t_list R_extr(UseJack);


  
  for(int ir=0;ir<(signed)R.size();ir++) {
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


     distr_t r(UseJack), r_const(UseJack), D_r(UseJack);
     
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
     double a_max = 0.1;
     
     for(int i=0;i< (int)(a_max/0.001);i++) {
       double x = i*0.001;
       R_to_print.distr_list.push_back( R_ave + D_R_ave*pow(x*fm_to_inv_Gev,2));
       a_to_print.push_back(x);
     }

         

     boost::filesystem::create_directory("../data/heavy_radiative/Fits_R");
     
     Print_To_File({}, {a_to_print, (R_to_print).ave(), (R_to_print).err()} , "../data/heavy_radiative/Fits_R/R"+to_string(ir)+"_extr.fit", "", "");
     
     Print_To_File({} , { (aR[ir]/fm_to_inv_Gev).ave(), (R[ir]).ave(), (R[ir]).err() }, "../data/heavy_radiative/Fits_R/R"+to_string(ir)+".data", "", "");

     R_extr.distr_list.push_back( R_ave);

  }

  //global mass-cont extrapolation for ratios


  //Print FFs

  boost::filesystem::create_directory("../data/heavy_radiative/Fits_B");


  distr_t_list m_etah_to_use= extrapolated_m_etah_AVE;
  for(int i=0;i<4;i++) m_etah_to_use.distr_list[i]= extrapolated_m_etah_SF.distr_list[i];
  

  //get perturbative velocity v= alpha( v*Meta/2);
  Vfloat v_heavy_pert;
  for(int i=0;i<(signed)m_etah_to_use.size();i++) {
    v_heavy_pert.push_back( Get_alpha_vM( m_etah_to_use.ave(i)));
  }

  distr_t_list F_H(UseJack);
  F_H.distr_list.push_back( Fm1_ave);
  for(int i=0;i<(signed)R_extr.size();i++) {
    F_H.distr_list.push_back( F_H.distr_list[i]*R_extr.distr_list[i]);
  }

  Print_To_File({}, { m_etah_to_use.ave(), m_etah_to_use.err(), extrapolated_m_Hh.ave(), extrapolated_m_Hh.err(),  F_H.ave(), F_H.err(), v_heavy_pert }, "../data/heavy_radiative/Fits_B/F_H.dat", "", "");


  //continuum limit extrapolation of F_H


  int Nfits=10; //13
  int Nfits_glb=Nfits*3;
  //Vint Npars({1,2,2,3});
  //Vint Npars({1,2,2,3,1,2,2,3,2,3,2,2,3}); //npars for each fit
  Vint Npars({1,2,2,3,2,2,3,2,2,3});

  Vfloat  metah_to_print;
  int Nsteps=300;
  double metah_min=2.0;
  double metah_max=15.0;
  for(int istep=0;istep<Nsteps;istep++) { metah_to_print.push_back( metah_min + istep*(metah_max-metah_min)/Nsteps ) ; }
  Vfloat vhpert_to_print;
  for(int istep=0;istep<Nsteps;istep++) { vhpert_to_print.push_back( Get_alpha_vM(metah_to_print[istep])); }

  distr_t_list Fb_vals(UseJack);
  distr_t_list Rpars(UseJack);
  Vfloat ch2_vals, npars_vals, ndof_vals;

  for(int ifit=0;ifit<Nfits_glb;ifit++) {

    
  class ipar_Fb {
  public:
    ipar_Fb() : FH(0.0), FH_err(0.0), m_r(0.0), vpert_r(0.0) {}
    double FH, FH_err;
    double m_r, vpert_r;
  };
  
  class fpar_Fb {
  public:
    fpar_Fb() {}
    fpar_Fb(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_Fb  class constructor Vfloat par has size != 3");
      R=par[0];
      A=par[1];
      A2=par[2];
    }
    double R,A,A2;
  };

  int offset=0;
  int Nmeas_FH= F_H.size();
  if(ifit/Nfits == 1) { Nmeas_FH--; offset=1;}
  if(ifit/Nfits == 2) { Nmeas_FH -= 2; offset=2;}
  if(ifit/Nfits == 3) { Nmeas_FH -=3; offset=3;}
    
    //fit on mean values to get ch2
    bootstrap_fit<fpar_Fb,ipar_Fb> bf_Fb(Njacks);
    bf_Fb.set_warmup_lev(1); //sets warmup
    bf_Fb.Set_number_of_measurements(Nmeas_FH);
    bf_Fb.Set_verbosity(1);
    bf_Fb.Add_par("R", 1.0, 0.1);
    bf_Fb.Add_par("A", 1.0, 0.1);
    bf_Fb.Add_par("A2", 1.0, 0.1);
    //for mean values
    bootstrap_fit<fpar_Fb,ipar_Fb> bf_Fb_ch2(1);
    bf_Fb_ch2.set_warmup_lev(1); //sets warmup
    bf_Fb_ch2.Set_number_of_measurements(Nmeas_FH);
    bf_Fb_ch2.Set_verbosity(1);

       
    bf_Fb_ch2.Add_par("R", 1.0, 0.1);
    bf_Fb_ch2.Add_par("A", 1.0, 0.1);
    bf_Fb_ch2.Add_par("A2", 1.0, 0.1);
    //##############################//

    //FIX FIT PARAMETERS
    if(Npars[(ifit+Nfits)%Nfits] == 1) { bf_Fb.Fix_par("A2",0.0); bf_Fb_ch2.Fix_par("A2",0.0);   bf_Fb.Fix_par("A",0.0); bf_Fb_ch2.Fix_par("A",0.0) ; }
    if(Npars[(ifit+Nfits)%Nfits] == 2) { bf_Fb.Fix_par("A2",0.0); bf_Fb_ch2.Fix_par("A2",0.0);}


 

    //ansatz
    bf_Fb.ansatz=  [&ifit, &Nfits](const fpar_Fb &p, const ipar_Fb &ip) {

      if((ifit+Nfits)%Nfits==0) { return p.R/pow(ip.m_r,0.5)  ; }   //1par
      else if((ifit+Nfits)%Nfits==1) { return (p.R/pow(ip.m_r,0.5))*( 1 + p.A/ip.m_r ) ; }  //2par
      else if((ifit+Nfits)%Nfits==2) { return (p.R/pow(ip.m_r,0.5))*( 1 + p.A/pow(ip.m_r,2));} //2par
      else if((ifit+Nfits)%Nfits==3) { return p.R/pow(ip.m_r,0.5)*( 1 + p.A/ip.m_r + p.A2/pow(ip.m_r,2)) ; }  //3par
      //else if((ifit+Nfits)%Nfits==4) { return p.R/(ip.m_r*ip.vpert_r); } //1par
      //else if((ifit+Nfits)%Nfits==5) { return (p.R/(ip.m_r*ip.vpert_r))*( 1 + p.A/ip.m_r) ; } //2par
      //else if((ifit+Nfits)%Nfits==6) { return (p.R/(ip.m_r*ip.vpert_r))*( 1 + p.A*pow(ip.vpert_r,2));} //2par
      //else if((ifit+Nfits)%Nfits==7) { return (p.R/(ip.m_r*ip.vpert_r))*( 1 + p.A/ip.m_r + p.A2*pow(ip.vpert_r,2));} //3par
      else if((ifit+Nfits)%Nfits==4) { return p.R + p.A/ip.m_r     ; } //2par
      else if((ifit+Nfits)%Nfits==5) { return p.R + p.A/pow(ip.m_r,2); } //2par
      else if((ifit+Nfits)%Nfits==6) { return p.R + p.A/ip.m_r + p.A2/pow(ip.m_r,2) ; }  //3par
      else if((ifit+Nfits)%Nfits==7) { return (p.R/(ip.m_r/pow(ip.m_r,2.0/3.0)))*(1.0 + p.A/pow(ip.m_r,2.0)) ; } //2par
      else if((ifit+Nfits)%Nfits==8) { return (p.R/(ip.m_r/pow(ip.m_r,2.0/3.0)))*(1.0 + p.A/pow(ip.m_r,1.0));} //2par
      else if((ifit+Nfits)%Nfits==9) { return  (p.R/(ip.m_r/pow(ip.m_r,2.0/3.0)))*(1.0 + p.A/pow(ip.m_r,2.0)  + p.A2/pow(ip.m_r,1.0)) ; } //3par
      else crash("ifit: "+to_string(ifit)+" not found") ;

    

    
      
      return 0.0;
           
    };
    
    
    bf_Fb.measurement =   [](const fpar_Fb &p, const ipar_Fb &ip) {
      return ip.FH;
    };
    bf_Fb.error=  [](const fpar_Fb &p, const ipar_Fb &ip) {
      return ip.FH_err;
    };
    
    bf_Fb_ch2.ansatz= bf_Fb.ansatz;
    bf_Fb_ch2.measurement = bf_Fb.measurement;
    bf_Fb_ch2.error = bf_Fb.error;


    //fill the data
    vector<vector<ipar_Fb>> data_Fb(Njacks);
    vector<vector<ipar_Fb>> data_Fb_ch2(1);
    //allocate space for output result



    //insert covariance matrix


    Eigen::MatrixXd Cov_Matrix(Nmeas_FH,Nmeas_FH);
    Eigen::MatrixXd Corr_Matrix(Nmeas_FH,Nmeas_FH);


    ifstream Load_Corr;

    if(UseJack) Load_Corr.open("../data/heavy_radiative/extr_Fb/corr_offset_"+to_string(offset)+".boot");


    for(int i=0;i<Nmeas_FH;i++) for(int j=0;j<Nmeas_FH;j++) {Cov_Matrix(i,j)=0; Corr_Matrix(i,j)=0;}
     
     

     for(int i=0; i<Nmeas_FH;i++) {
       for(int j=0; j<Nmeas_FH;j++) {

	 if(!UseJack) {
	   Cov_Matrix(i,j) = 0.25*(F_H.distr_list[i+offset]%F_H.distr_list[j+offset]);
	   Corr_Matrix(i,j) = Cov_Matrix(i,j)/(0.25*F_H.err(i+offset)*F_H.err(j+offset));
	 }
	 else {
	   double c;
	   Load_Corr >> c;
	   Corr_Matrix(i,j) = c;
	   Cov_Matrix(i,j) = Corr_Matrix(i,j)*0.25*(F_H.err(i+offset)*F_H.err(j+offset));
	   
	 }

       }
     }

     if(!UseJack) {
       ofstream Print_Corr("../data/heavy_radiative/extr_Fb/corr_ifit_"+to_string(ifit)+".boot");
       Print_Corr<<Corr_Matrix<<endl;
       Print_Corr.close();
     }
     else Load_Corr.close();
     

     bf_Fb.Add_covariance_matrix(Cov_Matrix);
     bf_Fb_ch2.Add_covariance_matrix(Cov_Matrix);

       

  


     for(auto &data_iboot: data_Fb) data_iboot.resize(Nmeas_FH);
     for(auto &data_iboot: data_Fb_ch2) data_iboot.resize(Nmeas_FH);
     for(int ijack=0;ijack<Njacks;ijack++) {
       for(int iens=0;iens<Nmeas_FH;iens++) {
	 data_Fb[ijack][iens].FH= 0.5*F_H.distr_list[iens+offset].distr[ijack];
	 data_Fb[ijack][iens].FH_err = 0.5*F_H.err(iens+offset);
	 data_Fb[ijack][iens].m_r = m_etah_to_use.distr_list[iens+offset].distr[ijack];
	 data_Fb[ijack][iens].vpert_r =  v_heavy_pert[iens+offset];
	 //mean values
	 if(ijack==0) {
	   data_Fb_ch2[ijack][iens].FH= 0.5*F_H.ave(iens+offset); 
	   data_Fb_ch2[ijack][iens].FH_err = 0.5*F_H.err(iens+offset);
	   data_Fb_ch2[ijack][iens].m_r = m_etah_to_use.ave(iens+offset);
	   data_Fb_ch2[ijack][iens].vpert_r =  v_heavy_pert[iens+offset];
	 }
       }
     }
     

     //append
     bf_Fb.Append_to_input_par(data_Fb);
     bf_Fb_ch2.Append_to_input_par(data_Fb_ch2);



     //fit
     boot_fit_data<fpar_Fb> Bt_fit_R;
     boot_fit_data<fpar_Fb> Bt_fit_R_ch2;
     Bt_fit_R= bf_Fb.Perform_bootstrap_fit();
     Bt_fit_R_ch2= bf_Fb_ch2.Perform_bootstrap_fit();
     
     distr_t Rp(UseJack), Ap(UseJack), A2p(UseJack);
     
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {
       Rp.distr.push_back( Bt_fit_R.par[ijack].R);
       Ap.distr.push_back( Bt_fit_R.par[ijack].A);
       A2p.distr.push_back( Bt_fit_R.par[ijack].A2);
     }

     Rpars.distr_list.push_back(Rp);
     
     //get ch2
     double ch2_R = Bt_fit_R_ch2.get_ch2_ave();
     double Npars= bf_Fb.Get_number_of_fit_pars();
     double Nmeas= bf_Fb.Get_number_of_measurements();
     double red_ch2_R= ch2_R/(Nmeas-Npars);
     ch2_vals.push_back(red_ch2_R);
     npars_vals.push_back( 1.0*Npars);
     ndof_vals.push_back( 1.0*(Nmeas-Npars));
     
     
     //print fitting function
    
    

     distr_t_list Fh_to_print(UseJack, Nsteps);
     distr_t Fb(UseJack);

     double mphys=m_etab;
     double vphys=Get_alpha_vM(mphys);
     
     for(int i=0;i< Nsteps;i++) {
       for(int ij=0;ij<Njacks;ij++) {

	 double mh=metah_to_print[i];
	 //double vh=  Get_alpha_vM(mh);

	 ipar_Fb iFb;
	 iFb.m_r = mh;
	 iFb.vpert_r = vhpert_to_print[i];

	 ipar_Fb iFb_phys;
	 iFb_phys.m_r = mphys;
	 iFb_phys.vpert_r= vphys;

	 fpar_Fb fFb;
	 fFb.R = Rp.distr[ij];
	 fFb.A = Ap.distr[ij];
	 fFb.A2 = A2p.distr[ij];

	 Fh_to_print.distr_list[i].distr.push_back( bf_Fb.ansatz(fFb,iFb) );

	 if(i==0) Fb.distr.push_back( bf_Fb.ansatz(fFb,iFb_phys));

       }
     }


     Fb_vals.distr_list.push_back(Fb);

     //print ansatz

     boost::filesystem::create_directory("../data/heavy_radiative/extr_Fb");

     Print_To_File({}, {metah_to_print, Fh_to_print.ave(), Fh_to_print.err() }, "../data/heavy_radiative/extr_Fb/ifit_"+to_string(ifit)+".extr", "", "");


     //print covariance matrix

     cout<<"Printing correlation matrix: "<<endl;
     cout<<Corr_Matrix<<endl;

   
     
  }


  double alpha_mc_ave = 1.0/133.5;
  double alpha_mc_err = 0.0000280548;
  double alpha_mb_ave = 1.0/132.2;
  double alpha_mb_err = 0.00001716557;

  distr_t alpha_mc(UseJack), alpha_mb(UseJack);

  for(int n=0;n<Njacks;n++) {

    alpha_mc.distr.push_back( alpha_mc_ave + alpha_mc_err*GM()/sqrt(Njacks-1.0));
    alpha_mb.distr.push_back( alpha_mb_ave + alpha_mb_err*GM()/sqrt(Njacks-1.0));

  }

  //lambda function to compute decay width

  auto Gamma = [](double Q, const distr_t &F2, double mh, double meta) -> distr_t {

    return (1e6)*(2.0/3.0)*Q*Q*alpha*F2*( pow(mh,2) - pow(meta,2))/(mh) ;

  };

  auto Gamma_b = [&alpha_mb](double Q, const distr_t &F2, double mh, double meta) -> distr_t {

    return (1e6)*(2.0/3.0)*Q*Q*alpha_mb*F2*( pow(mh,2) - pow(meta,2))/(mh) ;

  };

  auto Gamma_c = [&alpha_mc](double Q, const distr_t &F2, double mh, double meta) -> distr_t {

    return (1e6)*(2.0/3.0)*Q*Q*alpha_mc*F2*( pow(mh,2) - pow(meta,2))/(mh) ;

  };

  //print form factors from all the Nfits fits

  for(int ifit=0;ifit<Nfits_glb;ifit++) {

    if(ndof_vals[ifit] > 1) {

    cout<<"#### ifit: "<<ifit<<"####"<<endl;
    cout<<"Fb: "<<Fb_vals.ave(ifit)<<" +- "<<Fb_vals.err(ifit)<<endl;
    distr_t F = Fb_vals.distr_list[ifit];
    distr_t G= Gamma(Qb, F*F, m_hb, m_etab);
    cout<<"Gamma[hb->etab+gamma]: "<<G.ave()<<" +- "<<G.err()<<" [keV]"<<endl;
    cout<<"ch2/dof: "<<ch2_vals[ifit]<<endl;
    cout<<"R: "<<Rpars.ave(ifit)<<" +- "<<Rpars.err(ifit)<<endl;

    }
    
  }


  //grand average

  distr_t Fb_ave = 0.0*Get_id_distr(Njacks,UseJack);
  distr_t Fb_AIC = 0.0*Get_id_distr(Njacks,UseJack);
  distr_t Fb_A_AIC = 0.0*Get_id_distr(Njacks,UseJack);
  distr_t Fb_B_AIC = 0.0*Get_id_distr(Njacks,UseJack);
  distr_t Fb_C_AIC = 0.0*Get_id_distr(Njacks,UseJack);
  double sum_AIC=0.0;
  double sum_AIC_A=0.0;
  double sum_AIC_B=0.0;
  double sum_AIC_C=0.0;
  Vfloat ws;
  Vfloat ws_A;
  Vfloat ws_B;
  Vfloat ws_C;

  int Nfits_eff=0;

  for(int i=0;i<Nfits_glb;i++) {

    if(ndof_vals[i] > 1) {

      if( (ch2_vals[i] < (1.0 + sqrt(2.0/ndof_vals[i])))) {
	Fb_ave = Fb_ave + Fb_vals[i]; Nfits_eff++;
      }
      
      double w_AIC = exp(-0.5*(ch2_vals[i]*ndof_vals[i] + 2*npars_vals[i] -2*(npars_vals[i]+ndof_vals[i])  ));
      
      ws.push_back(w_AIC);

      if( (i+Nfits)%Nfits < 4) { ws_A.push_back( w_AIC ); ws_B.push_back( 0.0); ws_C.push_back( 0.0) ; }
      else if( (i+Nfits)%Nfits < 7) { ws_A.push_back( 0.0); ws_B.push_back( w_AIC ); ws_C.push_back( 0.0 ); }
      else if( (i+Nfits)%Nfits < 10) { ws_A.push_back( 0.0); ws_B.push_back( 0.0) ; ws_C.push_back( w_AIC) ; }
      else crash("FIT NOT FOUND");
      
      Fb_AIC = Fb_AIC + w_AIC*Fb_vals[i];
      Fb_A_AIC = Fb_A_AIC + ws_A[i]*Fb_vals[i];
      Fb_B_AIC = Fb_B_AIC + ws_B[i]*Fb_vals[i];
      Fb_C_AIC = Fb_C_AIC + ws_C[i]*Fb_vals[i];
      
      sum_AIC += w_AIC;

      sum_AIC_A += ws_A[i];
      sum_AIC_B += ws_B[i];
      sum_AIC_C += ws_C[i];

    }
    else { ws.push_back(0.0); ws_A.push_back(0.0) ; ws_B.push_back(0.0); ws_C.push_back(0.0) ; }
    
  }

  Fb_ave = Fb_ave/Nfits_eff;

  Fb_AIC = Fb_AIC/sum_AIC;
  Fb_A_AIC = Fb_A_AIC/sum_AIC_A;
  Fb_B_AIC = Fb_B_AIC/sum_AIC_B;
  Fb_C_AIC = Fb_C_AIC/sum_AIC_C;
  

  double syst=0, syst_AIC=0, syst_AIC_A=0, syst_AIC_B=0, syst_AIC_C=0;
  for(int i=0;i<Nfits_glb;i++) {

    if(ndof_vals[i] > 1) {
    
      if(ch2_vals[i] < (1.0 + sqrt(2.0/ndof_vals[i]))) syst += (1.0/Nfits_eff)*pow( Fb_vals[i].ave() - Fb_ave.ave(),2);
    
      syst_AIC += (ws[i]/sum_AIC)*pow( Fb_vals[i].ave() - Fb_AIC.ave(),2);
      syst_AIC_A += (ws_A[i]/sum_AIC_A)*pow( Fb_vals[i].ave() - Fb_A_AIC.ave(),2);
      syst_AIC_B += (ws_B[i]/sum_AIC_B)*pow( Fb_vals[i].ave() - Fb_B_AIC.ave(),2);
      syst_AIC_C += (ws_C[i]/sum_AIC_C)*pow( Fb_vals[i].ave() - Fb_C_AIC.ave(),2);
    }
  }
  syst = sqrt(syst);
  syst_AIC = sqrt(syst_AIC);
  syst_AIC_A = sqrt(syst_AIC_A);
  syst_AIC_B = sqrt(syst_AIC_B);
  syst_AIC_C = sqrt(syst_AIC_C);

  //add systematic error
  Fb_ave = Fb_ave.ave() + (Fb_ave - Fb_ave.ave())*sqrt( 1.0 + pow( syst/Fb_ave.err(),2));
  Fb_AIC = Fb_AIC.ave() + (Fb_AIC - Fb_AIC.ave())*sqrt( 1.0 + pow( syst_AIC/Fb_AIC.err(),2));
  Fb_A_AIC = Fb_A_AIC.ave() + (Fb_A_AIC - Fb_A_AIC.ave())*sqrt( 1.0 + pow( syst_AIC_A/Fb_A_AIC.err(),2));
  Fb_B_AIC = Fb_B_AIC.ave() + (Fb_B_AIC - Fb_B_AIC.ave())*sqrt( 1.0 + pow( syst_AIC_B/Fb_B_AIC.err(),2));
  Fb_C_AIC = Fb_C_AIC.ave() + (Fb_C_AIC - Fb_C_AIC.ave())*sqrt( 1.0 + pow( syst_AIC_C/Fb_C_AIC.err(),2));
  
  distr_t Fb_final= 0.5*(Fb_ave+Fb_AIC);
  Fb_final = Fb_final.ave() + (Fb_final - Fb_final.ave())*sqrt( 1.0 + pow( (Fb_AIC.ave()-Fb_ave.ave())/Fb_final.err(),2));
  
  
  cout<<"### AVERAGE ###"<<endl;
  cout<<"Fb(uni-ave): "<<Fb_ave.ave()<<" +- "<<Fb_ave.err()<<endl;
  cout<<"Fb(AIC): "<<Fb_AIC.ave()<<" +- "<<Fb_AIC.err()<<endl;
  cout<<"Fb(AIC-A): "<<Fb_A_AIC.ave()<<" +- "<<Fb_A_AIC.err()<<endl;
  cout<<"Fb(AIC-B): "<<Fb_B_AIC.ave()<<" +- "<<Fb_B_AIC.err()<<endl;
  cout<<"Fb(AIC-C): "<<Fb_C_AIC.ave()<<" +- "<<Fb_C_AIC.err()<<endl;
     
  cout<<"Fb(combined): "<<Fb_final.ave()<<" +- "<<Fb_final.err()<<endl;
  cout<<"Gamma[hb->etab+gamma]: "<<Gamma(Qb, Fb_AIC*Fb_AIC, m_hb, m_etab).ave()<<" +- "<<Gamma(Qb, Fb_AIC*Fb_AIC, m_hb, m_etab).err()<<" [KeV] "<<endl;
  cout<<"Gamma[hb->etab+gamma](alpha(mb)): "<<Gamma_b(Qb, Fb_AIC*Fb_AIC, m_hb, m_etab).ave()<<" +- "<<Gamma_b(Qb, Fb_AIC*Fb_AIC, m_hb, m_etab).err()<<" [KeV] "<<endl;
  cout<<"alpha(mb)^-1 : "<<(1.0/alpha_mb).ave()<<" +- "<<(1.0/alpha_mb).err()<<endl;

  //distr_t alpha_mc(UseJack), alpha_mb(UseJack);
  
  
  //print form factor and decay width at the charm


  //print grand average

  ofstream Print("../data/heavy_radiative/Fits_B/Fb_fin.ave");

  Print<<m_etab<<"  "<<Fb_ave.ave()<<" +- "<<Fb_ave.err()<<"  "<<Fb_AIC.ave()<<"  "<<Fb_AIC.err()<<endl;

  Print.close();


  

  cout<<"#####  CHARM #######"<<endl;
  cout<<"Fc: "<<(0.5*Fc_ave).ave()<<" "<<(0.5*Fc_ave).err()<<endl;
  distr_t Gc = Gamma(Qc, 0.25*Fc_ave*Fc_ave, m_hc, m_etac);
  distr_t Gc_c = Gamma_c(Qc, 0.25*Fc_ave*Fc_ave, m_hc, m_etac);
  cout<<"Gamma[hc->etac+gamma]: "<<Gc.ave()<<" +- "<<Gc.err()<<" [keV]"<<endl;
  cout<<"Gamma[hc->etac+gamma](alpha(mc)): "<<Gc_c.ave()<<" +- "<<Gc_c.err()<<" [keV]"<<endl;
  cout<<"alpha(mc)^-1 : "<<(1.0/alpha_mc).ave()<<" +- "<<(1.0/alpha_mc).err()<<endl;
  
  
  

        

  



  return;

}


void Get_chi_c1_decay() {


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

    


  


   data_t A1A1, V1V1, V1V1_mot, V2V2_mot, V3V3_mot, P5P5_mot;
   data_t A2A2, A3A3, V2V2, V3V3;
   data_t A1A1_loc, A2A2_loc, A3A3_loc;
   data_t A1A1_TM;

   data_t A2V3_PT3_tw1, A3V2_PT3_tw1;
   data_t A2V3_PT3_tw2, A3V2_PT3_tw2;
   data_t A2V3_PT3_rev_tw1, A3V2_PT3_rev_tw1;

   

     A1A1.Read("../chi_c/tw1", "mes_contr_2PT", "A1A1", Sort_light_confs);
     A2A2.Read("../chi_c/tw1", "mes_contr_2PT", "A2A2", Sort_light_confs);
     A3A3.Read("../chi_c/tw1", "mes_contr_2PT", "A3A3", Sort_light_confs);
     A1A1_loc.Read("../chi_c/tw1", "mes_contr_2PT_LOC", "A1A1", Sort_light_confs);
     A2A2_loc.Read("../chi_c/tw1", "mes_contr_2PT_LOC", "A2A2", Sort_light_confs);
     A3A3_loc.Read("../chi_c/tw1", "mes_contr_2PT_LOC", "A3A3", Sort_light_confs);
     V1V1.Read("../chi_c/tw1", "mes_contr_2PT", "V1V1", Sort_light_confs);
     V2V2.Read("../chi_c/tw1", "mes_contr_2PT", "V2V2", Sort_light_confs);
     V3V3.Read("../chi_c/tw1", "mes_contr_2PT", "V3V3", Sort_light_confs);
     V1V1_mot.Read("../chi_c/tw1", "mes_contr_2PT_MOT", "V1V1", Sort_light_confs);
     V2V2_mot.Read("../chi_c/tw1", "mes_contr_2PT_MOT", "V2V2", Sort_light_confs);
     V3V3_mot.Read("../chi_c/tw1", "mes_contr_2PT_MOT", "V3V3", Sort_light_confs);
     P5P5_mot.Read("../chi_c/tw1", "mes_contr_2PT_MOT", "P5P5", Sort_light_confs);

     A1A1_TM.Read("../heavy_radiative", "mes_contr_2PT_M1", "A1A1", Sort_light_confs);

     A2V3_PT3_tw1.Read("../chi_c/tw1", "mes_contr_3PT_J1_tw1", "A2V3", Sort_light_confs);
     A3V2_PT3_tw1.Read("../chi_c/tw1", "mes_contr_3PT_J1_tw1", "A3V2", Sort_light_confs);

     A2V3_PT3_tw1.Read("../chi_c/tw1", "mes_contr_3PT_J1_tw1", "A2V3", Sort_light_confs);
     A3V2_PT3_tw1.Read("../chi_c/tw1", "mes_contr_3PT_J1_tw1", "A3V2", Sort_light_confs);

     A2V3_PT3_tw2.Read("../chi_c/tw2", "mes_contr_3PT_J1_tw2", "A2V3", Sort_light_confs);
     A3V2_PT3_tw2.Read("../chi_c/tw2", "mes_contr_3PT_J1_tw2", "A3V2", Sort_light_confs);

     A2V3_PT3_rev_tw1.Read("../chi_c_rev/tw1", "mes_contr_3PT_J2_tw1", "A1V3", Sort_light_confs);
     A3V2_PT3_rev_tw1.Read("../chi_c_rev/tw1", "mes_contr_3PT_J2_tw1", "A3V1", Sort_light_confs);
     

     //data for cont.extrapolation
     distr_t_list a_list(UseJack);
     distr_t_list MJPSI_list(UseJack);
     distr_t_list MCHI_list(UseJack) ;
     distr_t_list E1_list(UseJack) ;
     distr_t_list M2_ov_E1_list(UseJack) ;
     distr_t_list M2_list(UseJack);
     
     double syst_M1=200;
     double syst_E1=200;

     
     int Nens= A1A1.size;

     boost::filesystem::create_directory("../data/chi_c1");

     for(int iens=0;iens<Nens;iens++) {

       string Ens= A1A1.Tag[iens];
       
       boost::filesystem::create_directory("../data/chi_c1/"+Ens);


       bool compute_2nd_tw= (Ens=="cB211b.072.64");

       bool compute_rev =  (Ens=="cD211a.054.96" || Ens=="cA211a.12.48");

       
       cout<<"########### ANALYZING ENSEMBLE "<<Ens<<" #############"<<endl;


       //########################################## ENSEMBLE INFO ############################################
       int tw1, tw2;
       double amc;
       distr_t a_distr(UseJack);
       distr_t ZV_had(UseJack);
       double syst_mc;
       if(Ens=="cB211b.072.64") {
	 tw1=20;  tw2=30;
	 a_distr=a_B;
	 amc= 0.231567;
	 ZV_had= ZV_B;
	 syst_mc= syst_mc_B;
       }
       else if(Ens=="cA211a.12.48") {
	 tw1=19;
	 a_distr=a_A;
	 amc= 0.262;
	 ZV_had= ZV_A;
	 syst_mc= syst_mc_A;
       }
       else if(Ens=="cD211a.054.96") {
	 tw1=30; //tw2=30;
	 a_distr=a_D;
	 amc= 0.164898;
	 ZV_had= ZV_D;
	 syst_mc= syst_mc_D;
       }
       else if(Ens=="cC211a.06.80") {
	 tw1=25; //tw2=30;
	 a_distr=a_C;
	 amc= 0.1984;
	 ZV_had=ZV_C;
	 syst_mc= syst_mc_C;
       }
       else if(Ens=="cE211a.044.112") {
	 tw1=35; //tw2=30;
	 a_distr=a_E;
	 amc= 0.141255;
	 ZV_had=ZV_E;
	 syst_mc= syst_mc_E;
       }
       else crash("Error");
       //#################################################################################################

       distr_t ZV= ZV_had;
       
       
       CorrAnalysis Corr(UseJack,Njacks,Njacks);
       Corr.Nt= A1A1.nrows[iens];
       Corr.Reflection_sign=1;
       Corr.Perform_Nt_t_average=1;

       if(Ens=="cB211b.072.64") {
	 int hens=-1;
	 for(int b=0;b<(signed)A1A1_TM.size; b++ ){
	   if(A1A1_TM.Tag[b] == Ens)  hens=b;
	 }
	 assert(hens != -1);
	 distr_t_list eff_A1A1= Corr.effective_mass_t(A1A1_TM.col(0)[hens], "../data/chi_c1/"+Ens+"/A1A1_TM");
       }
       
       
       distr_t_list eff_A1A1= Corr.effective_mass_t(A1A1.col(0)[iens], "../data/chi_c1/"+Ens+"/A1A1");
       distr_t_list eff_V1V1= Corr.effective_mass_t(V1V1.col(0)[iens], "../data/chi_c1/"+Ens+"/V1V1");
       distr_t_list eff_V1V1_mot= Corr.effective_mass_t(V1V1_mot.col(0)[iens], "../data/chi_c1/"+Ens+"/V1V1_mot");
       distr_t_list eff_V3V3_mot= Corr.effective_mass_t(V3V3_mot.col(0)[iens], "../data/chi_c1/"+Ens+"/V3V3_mot");
       distr_t_list eff_P5P5_mot= Corr.effective_mass_t(P5P5_mot.col(0)[iens], "../data/chi_c1/"+Ens+"/P5P5_mot");

       distr_t_list C2_chic1= Corr.corr_t( Multiply_Vvector_by_scalar( summ_master( A1A1.col(0)[iens], A2A2.col(0)[iens], A3A3.col(0)[iens] ), 1.0/3.0), "");
       distr_t_list C2_chic1_LOC= Corr.corr_t( Multiply_Vvector_by_scalar( summ_master( A1A1_loc.col(0)[iens], A2A2_loc.col(0)[iens], A3A3_loc.col(0)[iens] ), 1.0/3.0), "");
       distr_t_list C2_JPSI_T= Corr.corr_t( Multiply_Vvector_by_scalar( summ_master( V1V1_mot.col(0)[iens], V2V2_mot.col(0)[iens] ), 1.0/2.0), "");
       distr_t_list C2_JPSI_L= Corr.corr_t( V3V3_mot.col(0)[iens], "");
       distr_t_list C2_JPSI_REST= Corr.corr_t( Multiply_Vvector_by_scalar( summ_master( V1V1.col(0)[iens], V2V2.col(0)[iens], V3V3.col(0)[iens] ), 1.0/3.0), "");

       distr_t_list eff_chic1 = Corr.effective_mass_t(C2_chic1, "../data/chi_c1/"+Ens+"/eff_mass_chi_c1");
       distr_t_list eff_chic1_LOC = Corr.effective_mass_t(C2_chic1_LOC, "../data/chi_c1/"+Ens+"/eff_mass_chi_c1_loc");
       distr_t_list eff_JPSI_T = Corr.effective_mass_t( C2_JPSI_T, "../data/chi_c1/"+Ens+"/eff_mass_JPSI_T");
       distr_t_list eff_JPSI_L = Corr.effective_mass_t( C2_JPSI_L, "../data/chi_c1/"+Ens+"/eff_mass_JPSI_L");
       distr_t_list eff_JPSI_REST = Corr.effective_mass_t( C2_JPSI_REST, "");

       Print_To_File({}, { (eff_chic1/a_distr).ave(), (eff_chic1/a_distr).err() },  "../data/chi_c1/"+Ens+"/eff_mass_chi_c1_pu", "", "");
       Print_To_File({}, { (eff_JPSI_REST/a_distr).ave(), (eff_JPSI_REST/a_distr).err() },  "../data/chi_c1/"+Ens+"/eff_mass_JPSI_pu", "", "");

       if(Ens=="cB211b.072.64") {Corr.Tmin=22; Corr.Tmax=27;}
       else if(Ens=="cA211a.12.48") {Corr.Tmin=18; Corr.Tmax=22;}
       else if(Ens=="cD211a.054.96") {Corr.Tmin = 28; Corr.Tmax = 34;}
       else if(Ens=="cC211a.06.80") {Corr.Tmin = 22; Corr.Tmax = 28;}
       else if(Ens=="cE211a.044.112") {Corr.Tmin = 24; Corr.Tmax = 31;}
       else crash("Ens: "+Ens+" not found");

       

       distr_t_list Z_chic1 = Corr.matrix_element_t(C2_chic1, "../data/chi_c1/"+Ens+"/Z_chi_c1");
       distr_t M_chic1 = Corr.Fit_distr(eff_chic1);
       distr_t M_chic1_loc= Corr.Fit_distr(eff_chic1_LOC);
       //add syst uncertainty from sm-loc
       double syst= 0.0; // fabs( M_chic1.ave() - M_chic1_loc.ave());
       M_chic1 = M_chic1.ave() + (M_chic1 -M_chic1.ave())*sqrt( 1.0 + pow( syst/M_chic1.err(),2));
       distr_t_list CHI_MEM= EXPT_D(-1.0*M_chic1, Corr.Nt)*Corr.Fit_distr(Z_chic1)/(2*M_chic1);
       distr_t_list CHI_MEM_tw1= EXP_D(1.0*M_chic1*tw1)*CHI_MEM;
       distr_t_list CHI_MEM_tw2= EXP_D(1.0*M_chic1*tw2)*CHI_MEM;

       
       distr_t_list chic1_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks,UseJack) ;
       for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	 chic1_fit.distr_list[t] = M_chic1;
       }

       

       Corr.Tmin= (int)( 2.2*fm_to_inv_Gev/a_distr.ave() );
       Corr.Tmax= (int)( 3.2*fm_to_inv_Gev/a_distr.ave() );
       
       distr_t_list Z_JPSI_T = Corr.matrix_element_t(C2_JPSI_T, "");
       distr_t_list Z_JPSI_L = Corr.matrix_element_t(C2_JPSI_L, "");
       distr_t E_JPSI_T = Corr.Fit_distr( eff_JPSI_T);
       distr_t E_JPSI_L = Corr.Fit_distr( eff_JPSI_L);
       distr_t M_JPSI= Corr.Fit_distr( eff_JPSI_REST);

       distr_t_list MJPSI_fit = 0.0*Get_id_distr_list(Corr.Nt, Njacks,UseJack) ;
       for(int t=Corr.Tmin;t<= Corr.Tmax;t++) {
	 MJPSI_fit.distr_list[t] = M_JPSI;
       }

     
       Print_To_File({}, { (chic1_fit/a_distr).ave(), (chic1_fit/a_distr).err() },  "../data/chi_c1/"+Ens+"/chi_c1_fit_pu", "", "");
       Print_To_File({}, { (MJPSI_fit/a_distr).ave(), (MJPSI_fit/a_distr).err() },  "../data/chi_c1/"+Ens+"/JPSI_fit_pu", "", "");
       

       distr_t k = (M_chic1-E_JPSI_T)/a_distr;
       distr_t k_bis = (M_chic1*M_chic1 - M_JPSI*M_JPSI)/(2.0*M_chic1)/a_distr ;
           
      
    
       
       distr_t JPSI_L_MEM_tw1 =  EXP_D(-1.0*E_JPSI_L*tw1)*Corr.Fit_distr(Z_JPSI_L)/(2*E_JPSI_L);
       distr_t JPSI_T_MEM_tw1 =  EXP_D(-1.0*E_JPSI_T*tw1)*Corr.Fit_distr(Z_JPSI_T)/(2*E_JPSI_T);
       distr_t JPSI_L_MEM_tw2 =  EXP_D(-1.0*E_JPSI_L*tw2)*Corr.Fit_distr(Z_JPSI_L)/(2*E_JPSI_L);
       distr_t JPSI_T_MEM_tw2 =  EXP_D(-1.0*E_JPSI_T*tw2)*Corr.Fit_distr(Z_JPSI_T)/(2*E_JPSI_T);
       

       Corr.Perform_Nt_t_average=0;


       distr_t_list PT3_LONG_distr_RE_tw1 = 2.0*ZV*Corr.corr_t(A2V3_PT3_tw1.col(0)[iens], "../data/chi_c1/"+Ens+"/LONG_PT3_tw1_RE")/(CHI_MEM_tw1*JPSI_L_MEM_tw1);
       //distr_t_list PT3_LONG_distr_IM_tw1 = Corr.corr_t(A2V3_PT3_tw1.col(1)[iens], "../data/chi_c1/"+Ens+"/LONG_PT3_tw1_IM");

       distr_t_list PT3_TR_distr_RE_tw1 = 2.0*ZV*Corr.corr_t(A3V2_PT3_tw1.col(0)[iens], "../data/chi_c1/"+Ens+"/TR_PT3_tw1_RE")/(CHI_MEM_tw1*JPSI_T_MEM_tw1);
       //distr_t_list PT3_TR_distr_IM_tw1 = Corr.corr_t(A3V2_PT3_tw1.col(1)[iens], "../data/chi_c1/"+Ens+"/TR_PT3_tw1_IM");


       distr_t_list PT3_LONG_distr_RE_rev_tw1(UseJack), PT3_TR_distr_RE_rev_tw1(UseJack);

       if(compute_rev) {

	 //find ens id
	 int j_ens=-1;
	 for(int b=0;b<(signed)A2V3_PT3_rev_tw1.size; b++) {
	   if( A2V3_PT3_rev_tw1.Tag[b] == Ens) j_ens = b ;
	 }
	 assert(j_ens != -1 );
	 
	 
	 PT3_LONG_distr_RE_rev_tw1 = -2.0*ZV*Corr.corr_t(A2V3_PT3_rev_tw1.col(0)[j_ens], "../data/chi_c1/"+Ens+"/LONG_PT3_tw1_rev_RE")/(CHI_MEM_tw1*JPSI_L_MEM_tw1);
	 PT3_TR_distr_RE_rev_tw1 = -2.0*ZV*Corr.corr_t(A3V2_PT3_rev_tw1.col(0)[j_ens], "../data/chi_c1/"+Ens+"/TR_PT3_tw1_rev_RE")/(CHI_MEM_tw1*JPSI_T_MEM_tw1);
       }
       

       distr_t_list PT3_LONG_distr_RE_tw2(UseJack), PT3_LONG_distr_IM_tw2(UseJack);
       distr_t_list PT3_TR_distr_RE_tw2(UseJack), PT3_TR_distr_IM_tw2(UseJack);

  
       if(compute_2nd_tw) {
	 PT3_LONG_distr_RE_tw2 = 2.0*ZV*Corr.corr_t(A2V3_PT3_tw2.col(0)[0], "../data/chi_c1/"+Ens+"/LONG_PT3_tw2_RE")/(CHI_MEM_tw2*JPSI_L_MEM_tw2);
	 PT3_TR_distr_RE_tw2 = 2.0*ZV*Corr.corr_t(A3V2_PT3_tw2.col(0)[0], "../data/chi_c1/"+Ens+"/TR_PT3_tw2_RE")/(CHI_MEM_tw2*JPSI_T_MEM_tw2);
       }


       //averaged J1 and J2
       if(compute_rev) {
	 PT3_LONG_distr_RE_tw1 = 0.5*( PT3_LONG_distr_RE_tw1 + PT3_LONG_distr_RE_rev_tw1 );
	 PT3_TR_distr_RE_tw1 = 0.5* ( PT3_TR_distr_RE_tw1 + PT3_TR_distr_RE_rev_tw1 );
       }
       

       //summ and diff
  

       distr_t_list PT3_M1_distr_tw1 = (PT3_TR_distr_RE_tw1+PT3_LONG_distr_RE_tw1)/(sqrt(2.0)*M_chic1);
       distr_t_list PT3_E1_distr_tw1 = (PT3_TR_distr_RE_tw1-PT3_LONG_distr_RE_tw1)/(sqrt(2.0)*M_chic1);
       distr_t_list M1_E1_ratio_distr_tw1 = PT3_M1_distr_tw1/PT3_E1_distr_tw1;
       
      
       distr_t_list PT3_E1_distr_tw2(UseJack), PT3_M1_distr_tw2(UseJack), M1_E1_ratio_distr_tw2(UseJack);

       if(compute_2nd_tw) {
	 PT3_M1_distr_tw2 = (PT3_TR_distr_RE_tw2+PT3_LONG_distr_RE_tw2)/(sqrt(2.0)*M_chic1);
	 PT3_E1_distr_tw2 = (PT3_TR_distr_RE_tw2-PT3_LONG_distr_RE_tw2)/(sqrt(2.0)*M_chic1);
	 M1_E1_ratio_distr_tw2 = PT3_M1_distr_tw2/PT3_E1_distr_tw2;
       }



       distr_t_list PT3_E1_distr_tw1_red(UseJack), PT3_M1_distr_tw1_red(UseJack), M1_E1_ratio_distr_tw1_red(UseJack);

       distr_t_list PT3_E1_distr_tw2_red(UseJack), PT3_M1_distr_tw2_red(UseJack), M1_E1_ratio_distr_tw2_red(UseJack);

     

       for(int t=tw1;t<Corr.Nt;t++) {
	 PT3_E1_distr_tw1_red.distr_list.push_back( PT3_E1_distr_tw1.distr_list[t] );
	 PT3_M1_distr_tw1_red.distr_list.push_back( PT3_M1_distr_tw1.distr_list[t] );
	 M1_E1_ratio_distr_tw1_red.distr_list.push_back( M1_E1_ratio_distr_tw1.distr_list[t]) ;
       }

       if(compute_2nd_tw) {
	   for(int t=tw2;t<Corr.Nt;t++) {
	     PT3_E1_distr_tw2_red.distr_list.push_back( PT3_E1_distr_tw2.distr_list[t] );
	     PT3_M1_distr_tw2_red.distr_list.push_back( PT3_M1_distr_tw2.distr_list[t] );
	     M1_E1_ratio_distr_tw2_red.distr_list.push_back( M1_E1_ratio_distr_tw2.distr_list[t]);
	   }
	   
       }


       Print_To_File({}, {PT3_E1_distr_tw1_red.ave(), PT3_E1_distr_tw1_red.err(), PT3_M1_distr_tw1_red.ave(), PT3_M1_distr_tw1_red.err(), M1_E1_ratio_distr_tw1_red.ave(), M1_E1_ratio_distr_tw1_red.err() },  "../data/chi_c1/"+Ens+"/FF_tw1", "", "");

      
       if(compute_2nd_tw) {
	 Print_To_File({}, {PT3_E1_distr_tw2_red.ave(), PT3_E1_distr_tw2_red.err(), PT3_M1_distr_tw2_red.ave(), PT3_M1_distr_tw2_red.err(), M1_E1_ratio_distr_tw2_red.ave(), M1_E1_ratio_distr_tw2_red.err() },  "../data/chi_c1/"+Ens+"/FF_tw2", "", "");
       }

       
       int t1= (int)( 1.1*fm_to_inv_Gev/(a_distr.ave()));
       int t2= (int)( 1.5*fm_to_inv_Gev/(a_distr.ave()));

       
       
       Corr.Tmin=t1; Corr.Tmax=t2;
       
       distr_t E1_tw1= Corr.Fit_distr(PT3_E1_distr_tw1_red);
       distr_t E1_tw2(UseJack);

       
       
       if(compute_2nd_tw) {
	   E1_tw2= Corr.Fit_distr(PT3_E1_distr_tw2_red);
	   double dx= fabs( E1_tw1.ave() - E1_tw2.ave());
	   syst_E1 = dx*erf( dx/(sqrt(2.0)*sqrt( pow(E1_tw1.err(),2) + pow(E1_tw2.err(),2))))/E1_tw1.ave();
	   
       }
       
     

       if(Ens=="cB211b.072.64") {Corr.Tmin=8; Corr.Tmax=12;}
       else if(Ens=="cA211a.12.48") {Corr.Tmin=8; Corr.Tmax=11;}
       else if(Ens=="cD211a.054.96") {Corr.Tmin = 13; Corr.Tmax = 20;}
       else if(Ens=="cC211a.06.80") {Corr.Tmin = 10; Corr.Tmax = 15;}
       else if(Ens=="cE211a.044.112") {Corr.Tmin = 24; Corr.Tmax = 31;}
       else crash("Ens: "+Ens+" not found");
       
       distr_t M1_tw1= Corr.Fit_distr(PT3_M1_distr_tw1_red);
       distr_t a2_tw1 =  M1_tw1/SQRT_D( M1_tw1*M1_tw1 + E1_tw1*E1_tw1);
       distr_t a2_tw1_impr = Corr.Fit_distr(M1_E1_ratio_distr_tw1_red);

       distr_t M1_tw2(UseJack), a2_tw2(UseJack), a2_tw2_impr(UseJack);

       if(compute_2nd_tw) {
	
	  M1_tw2= Corr.Fit_distr(PT3_M1_distr_tw2_red);
	  a2_tw2 =  M1_tw2/SQRT_D( M1_tw2*M1_tw2 + E1_tw2*E1_tw2);
	  a2_tw2_impr = Corr.Fit_distr( M1_E1_ratio_distr_tw2_red);
	  double dx= fabs(M1_tw1.ave() - M1_tw2.ave());
	  syst_M1 = dx*erf( dx/(sqrt(2.0)*sqrt( pow(M1_tw1.err(),2) + pow(M1_tw2.err(),2))))/M1_tw1.ave();
	  
	  
       }

       //add systematic error

       cout<<"Ens: "<<Ens<<endl;
       cout<<"before adding syst: "<<endl;
       cout<<"E1: "<<E1_tw1.ave()<<E1_tw1.err()<<endl;
       cout<<"M2: "<<M1_tw1.ave()<<M1_tw1.err()<<endl;

       E1_tw1 = E1_tw1.ave() + (E1_tw1-E1_tw1.ave())*sqrt( 1.0 + pow(syst_E1*E1_tw1.ave()/E1_tw1.err(),2));
       M1_tw1 = M1_tw1.ave() + (M1_tw1-M1_tw1.ave())*sqrt( 1.0 + pow(syst_M1*M1_tw1.ave()/M1_tw1.err(),2));

       cout<<"after adding syst: "<<endl;
       cout<<"E1: "<<E1_tw1.ave()<<E1_tw1.err()<<endl;
       cout<<"M2: "<<M1_tw1.ave()<<M1_tw1.err()<<endl;
       
       cout<<"syst E1: "<<syst_E1<<endl;
       cout<<"syst M2: "<<syst_M1<<endl;

       
       

       //push-back
       MJPSI_list.distr_list.push_back(M_JPSI/a_distr);
       MCHI_list.distr_list.push_back(M_chic1/a_distr);
       E1_list.distr_list.push_back(E1_tw1);
       M2_ov_E1_list.distr_list.push_back( (M1_tw1/E1_tw1) );
       M2_list.distr_list.push_back( M1_tw1);
       a_list.distr_list.push_back( a_distr);

       
       distr_t rate= 389.4*alpha*(4/27.0)*( E1_tw1*E1_tw1 + M1_tw1*M1_tw1); //MeV


       
       
       cout<<"###################################"<<endl; 
       cout<<"Printing form factors for Ensemble: "<<Ens<<endl;
       cout<<"k: "<< k.ave()<<" +- "<<k.err()<<endl;
       cout<<"k: "<< k_bis.ave()<<" +- "<<k_bis.err()<<endl;
       cout<<"E1(tw1): "<<E1_tw1.ave()<<" +- "<<E1_tw1.err()<<endl;
       cout<<"M1(tw1): "<<M1_tw1.ave()<<" +- "<<M1_tw1.err()<<endl;
       cout<<"a2(tw1): "<<a2_tw1.ave()<<" +- "<<a2_tw1.err()<<endl;
       cout<<"a2_impr(tw2): "<<a2_tw1_impr.ave()<<" +- "<<a2_tw1_impr.err()<<endl;
       cout<<"rate: "<<rate.ave()<<" +- "<<rate.err()<<" [MeV]"<<endl;
       if(compute_2nd_tw) {
	 cout<<"E1(tw2): "<<E1_tw2.ave()<<" +- "<<E1_tw2.err()<<endl;
	 cout<<"M1(tw2): "<<M1_tw2.ave()<<" +- "<<M1_tw2.err()<<endl;
	 cout<<"a2(tw2): "<<a2_tw2.ave()<<" +- "<<a2_tw2.err()<<endl;
	 cout<<"a2_impr(tw2): "<<a2_tw2_impr.ave()<<" +- "<<a2_tw2_impr.err()<<endl;
       }
       cout<<"###################################"<<endl; 
          
  
       
       
       
       
       
     }
  
     cout<<"Eff mass computed for chi_c1 -> J/Psi + gamma"<<endl;
 


     //continuum-limit extrapolation
     //###################################################################################################



     
     class ipar {
     public:
       ipar() : E1(0.0), E1_err(0.0), M2E1(0.0), M2E1_err(0.0),  mJ(0.0), mJ_err(0.0), mChi(0.0), mChi_err(0.0), a(0.0) {}
       double E1, E1_err;
       double M2E1, M2E1_err;
       double mJ, mJ_err;
       double mChi, mChi_err;
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


     int Nmeas= a_list.size();
     
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

          
     
     string FIT_TAG="";

     bool exclude_coarsest=false;
     
     //ansatz
     bf_mh.ansatz=  [&exclude_coarsest](const fpar &p, const ipar &ip) {

       if(exclude_coarsest && ip.a > 0.084*fm_to_inv_Gev ) return 0.0;
       
       return p.R + p.A*pow(ip.a,2) + p.A2*pow(ip.a,4);
       
     };

  
     bf_mh.measurement=  [&FIT_TAG, &exclude_coarsest ](const fpar &p, const ipar &ip) {

       if(exclude_coarsest && ip.a > 0.084*fm_to_inv_Gev ) return 0.0;
       
       if(FIT_TAG=="Chi") return ip.mChi;
       else if(FIT_TAG=="J") return ip.mJ;
       else if(FIT_TAG=="E1") return ip.E1;
       return ip.M2E1;
       
     };
     bf_mh.error=  [&FIT_TAG ](const fpar &p, const ipar &ip) {
       
       if(FIT_TAG=="Chi") return ip.mChi_err;
       else if(FIT_TAG=="J") return ip.mJ_err;
       else if(FIT_TAG=="E1") return ip.E1_err;
       return ip.M2E1_err;
       
     };

     bf_mh_ch2.ansatz= bf_mh.ansatz;
     bf_mh_ch2.measurement = bf_mh.measurement;
     bf_mh_ch2.error = bf_mh.error;

         
  
     //fill the data
     vector<vector<ipar>> data(Njacks);
     vector<vector<ipar>> data_ch2(1);
     //allocate space for output result

     cout<<"Nens: "<<Nens<<endl;
     
     for(auto &data_iboot: data) data_iboot.resize(Nens);
     for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
     for(int ijack=0;ijack<Njacks;ijack++) {
       for(int iens=0;iens<Nens;iens++) {

	 cout<<"iens: "<<iens<<endl;
	 
	 data[ijack][iens].E1= E1_list.distr_list[iens].distr[ijack];
	 data[ijack][iens].E1_err = E1_list.err(iens);
	 data[ijack][iens].M2E1= M2_ov_E1_list.distr_list[iens].distr[ijack];
	 data[ijack][iens].M2E1_err = M2_ov_E1_list.err(iens);
	 data[ijack][iens].mJ = MJPSI_list.distr_list[iens].distr[ijack];
	 data[ijack][iens].mJ_err= MJPSI_list.err(iens);
	 data[ijack][iens].mChi= MCHI_list.distr_list[iens].distr[ijack];
	 data[ijack][iens].mChi_err= MCHI_list.err(iens);
	 data[ijack][iens].a = a_list.distr_list[iens].distr[ijack];
	 //mean values
	 if(ijack==0) {
	   data_ch2[ijack][iens].E1= E1_list.ave(iens);
	   data_ch2[ijack][iens].E1_err = E1_list.err(iens);
	   data_ch2[ijack][iens].M2E1= M2_ov_E1_list.ave(iens);
	   data_ch2[ijack][iens].M2E1_err = M2_ov_E1_list.err(iens);
	   data_ch2[ijack][iens].mJ = MJPSI_list.ave(iens);
	   data_ch2[ijack][iens].mJ_err= MJPSI_list.err(iens);
	   data_ch2[ijack][iens].mChi=  MCHI_list.ave(iens);
	   data_ch2[ijack][iens].mChi_err= MCHI_list.err(iens);
	   data_ch2[ijack][iens].a = a_list.ave(iens);
	   
	 }
       }
     }
     
         
     //append
     bf_mh.Append_to_input_par(data);
     bf_mh_ch2.Append_to_input_par(data_ch2);
     //fit
     boot_fit_data<fpar> Bt_fit_M2E1;
     boot_fit_data<fpar> Bt_fit_M2E1_ch2;
     FIT_TAG="";
     Bt_fit_M2E1= bf_mh.Perform_bootstrap_fit();
     Bt_fit_M2E1_ch2= bf_mh_ch2.Perform_bootstrap_fit();     
     boot_fit_data<fpar> Bt_fit_E1;
     boot_fit_data<fpar> Bt_fit_E1_ch2;
     FIT_TAG="E1";
     Bt_fit_E1= bf_mh.Perform_bootstrap_fit();
     Bt_fit_E1_ch2= bf_mh_ch2.Perform_bootstrap_fit();
     FIT_TAG="Chi";
     boot_fit_data<fpar> Bt_fit_Chi;
     boot_fit_data<fpar> Bt_fit_Chi_ch2;
     Bt_fit_Chi= bf_mh.Perform_bootstrap_fit();
     Bt_fit_Chi_ch2= bf_mh_ch2.Perform_bootstrap_fit();
     FIT_TAG="J";
     boot_fit_data<fpar> Bt_fit_J;
     boot_fit_data<fpar> Bt_fit_J_ch2;
     Bt_fit_J= bf_mh.Perform_bootstrap_fit();
     Bt_fit_J_ch2= bf_mh_ch2.Perform_bootstrap_fit();
     
     
     //release a4 coefficient
     
  

     //here I have to modify the fit ansatz by doing a constant fit in the case of E1 and by excluding the coarsest data point in the case of a2


     
     
     //fit

     exclude_coarsest=true;
     
     boot_fit_data<fpar> Bt_fit_M2E1_Q;
     boot_fit_data<fpar> Bt_fit_M2E1_ch2_Q;
     FIT_TAG="";
     Bt_fit_M2E1_Q= bf_mh.Perform_bootstrap_fit();
     Bt_fit_M2E1_ch2_Q = bf_mh_ch2.Perform_bootstrap_fit();     
     boot_fit_data<fpar> Bt_fit_E1_Q;
     boot_fit_data<fpar> Bt_fit_E1_ch2_Q;
     FIT_TAG="E1";

     exclude_coarsest=true;

     //bf_mh.Fix_par("A",0.0);
     //bf_mh_ch2.Fix_par("A",0.0);
     
     Bt_fit_E1_Q= bf_mh.Perform_bootstrap_fit();
     Bt_fit_E1_ch2_Q= bf_mh_ch2.Perform_bootstrap_fit();

     
     //bf_mh.Release_par("A");
     //bf_mh_ch2.Release_par("A");
     
     
     bf_mh.Release_par("A2");
     bf_mh_ch2.Release_par("A2");

     
     FIT_TAG="Chi";
     boot_fit_data<fpar> Bt_fit_Chi_Q;
     boot_fit_data<fpar> Bt_fit_Chi_ch2_Q;
     Bt_fit_Chi_Q= bf_mh.Perform_bootstrap_fit();
     Bt_fit_Chi_ch2_Q= bf_mh_ch2.Perform_bootstrap_fit();
     FIT_TAG="J";
     boot_fit_data<fpar> Bt_fit_J_Q;
     boot_fit_data<fpar> Bt_fit_J_ch2_Q;
     Bt_fit_J_Q= bf_mh.Perform_bootstrap_fit();
     Bt_fit_J_ch2_Q= bf_mh_ch2.Perform_bootstrap_fit();

     //pars 
     distr_t M_J(UseJack), M_J_Q(UseJack),  D_J(UseJack), D_J_Q(UseJack), D2_J(UseJack), D2_J_Q(UseJack);
     distr_t M_Chi(UseJack), M_Chi_Q(UseJack),  D_Chi(UseJack), D_Chi_Q(UseJack), D2_Chi(UseJack), D2_Chi_Q(UseJack);
     distr_t E1(UseJack), E1_Q(UseJack), D_E1(UseJack), D_E1_Q(UseJack), D2_E1(UseJack), D2_E1_Q(UseJack);
     distr_t M2E1(UseJack), M2E1_Q(UseJack),  D_M2E1(UseJack), D_M2E1_Q(UseJack), D2_M2E1(UseJack), D2_M2E1_Q(UseJack);

  
     //retrieve parameters
     for(int ijack=0;ijack<Njacks;ijack++) {


       //J/Psi
       M_J.distr.push_back( Bt_fit_J.par[ijack].R);
       M_J_Q.distr.push_back( Bt_fit_J_Q.par[ijack].R);
       D_J.distr.push_back( Bt_fit_J.par[ijack].A);
       D_J_Q.distr.push_back( Bt_fit_J_Q.par[ijack].A);
       D2_J.distr.push_back( Bt_fit_J.par[ijack].A2);
       D2_J_Q.distr.push_back( Bt_fit_J_Q.par[ijack].A2);

       //Chi
       M_Chi.distr.push_back( Bt_fit_Chi.par[ijack].R);
       M_Chi_Q.distr.push_back( Bt_fit_Chi_Q.par[ijack].R);
       D_Chi.distr.push_back( Bt_fit_Chi.par[ijack].A);
       D_Chi_Q.distr.push_back( Bt_fit_Chi_Q.par[ijack].A);
       D2_Chi.distr.push_back( Bt_fit_Chi.par[ijack].A2);
       D2_Chi_Q.distr.push_back( Bt_fit_Chi_Q.par[ijack].A2);

       //E1
       E1.distr.push_back( Bt_fit_E1.par[ijack].R);
       E1_Q.distr.push_back( Bt_fit_E1_Q.par[ijack].R);
       D_E1.distr.push_back( Bt_fit_E1.par[ijack].A);
       D_E1_Q.distr.push_back( Bt_fit_E1_Q.par[ijack].A);
       D2_E1.distr.push_back( Bt_fit_E1.par[ijack].A2);
       D2_E1_Q.distr.push_back( Bt_fit_E1_Q.par[ijack].A2);

       //M2/E1
       M2E1.distr.push_back( Bt_fit_M2E1.par[ijack].R);
       M2E1_Q.distr.push_back( Bt_fit_M2E1_Q.par[ijack].R);
       D_M2E1.distr.push_back( Bt_fit_M2E1.par[ijack].A);
       D_M2E1_Q.distr.push_back( Bt_fit_M2E1_Q.par[ijack].A);
       D2_M2E1.distr.push_back( Bt_fit_M2E1.par[ijack].A2);
       D2_M2E1_Q.distr.push_back( Bt_fit_M2E1_Q.par[ijack].A2);
        
     }


  
  
     //get ch2
     double ch2_J_lin= Bt_fit_J_ch2.get_ch2_ave();
     double ch2_J_Q = Bt_fit_J_ch2_Q.get_ch2_ave();

     double ch2_Chi_lin = Bt_fit_Chi_ch2.get_ch2_ave();
     double ch2_Chi_Q = Bt_fit_Chi_ch2_Q.get_ch2_ave();

     double ch2_E1_lin = Bt_fit_E1_ch2.get_ch2_ave();
     double ch2_E1_Q = Bt_fit_E1_ch2_Q.get_ch2_ave();

     double ch2_M2E1_lin = Bt_fit_M2E1_ch2.get_ch2_ave();
     double ch2_M2E1_Q = Bt_fit_M2E1_ch2_Q.get_ch2_ave();

     //build fitting function

     //Dof(linear_fit)= 2 , Dof(constant_fit)=3

     int dof_lin= Nmeas - 2;
     int dof_Q = Nmeas-3;

     cout<<"dof lin: "<<dof_lin<<endl;
     cout<<"dof Q: "<<dof_Q<<endl;
  

     double w_J_1 = exp(-0.5*(ch2_J_lin - 2*dof_lin))/( exp(-0.5*(ch2_J_lin - 2*dof_lin)) + exp(-0.5*(ch2_J_Q - 2*dof_Q)));
     double w_J_2 = 1.0 - w_J_1;

     double w_Chi_1 = exp(-0.5*(ch2_Chi_lin - 2*dof_lin))/( exp(-0.5*(ch2_Chi_lin - 2*dof_lin)) + exp(-0.5*(ch2_Chi_Q - 2*dof_Q)));
     double w_Chi_2 = 1.0 - w_Chi_1;

     double w_E1_1 = exp(-0.5*(ch2_E1_lin - 2*dof_lin))/( exp(-0.5*(ch2_E1_lin - 2*dof_lin)) + exp(-0.5*(ch2_E1_Q - 2*dof_Q)));
     double w_E1_2 = 1.0 - w_E1_1;

     double w_M2E1_1 = exp(-0.5*(ch2_M2E1_lin - 2*dof_lin))/( exp(-0.5*(ch2_M2E1_lin - 2*dof_lin)) + exp(-0.5*(ch2_M2E1_Q - 2*dof_Q)));
     double w_M2E1_2 = 1.0 - w_M2E1_1;


     //how do we fit
     w_J_1=1.0; w_J_2=0.0;
     w_Chi_1=1.0; w_Chi_2=0.0;
     //w_E1_1 = 1.0; w_E1_2 = 0.0;
     //w_M2E1_1 = 1.0; w_M2E1_2 = 0.0;

     cout<<"reduced ch2 (J): "<<ch2_J_lin/dof_lin<<endl;
     cout<<"reduced ch2 (Chi): "<<ch2_Chi_lin/dof_lin<<endl;
     cout<<"reduced ch2 (E1): "<<ch2_E1_lin/dof_lin<<" w: "<<w_E1_1<<endl;
     cout<<"reduced ch2 (E1-Q): "<<ch2_E1_Q/dof_Q<<" w: "<<w_E1_2<<endl;
     cout<<"reduced ch2 (M2/E1): "<<ch2_M2E1_lin/dof_Q<<" w: "<<w_M2E1_1<<endl;
     cout<<"reduced ch2 (M2/E1-Q): "<<ch2_M2E1_Q/dof_Q<<" w: "<<w_M2E1_2<<endl;

     distr_t M_Chi_ave= w_Chi_1*M_Chi + w_Chi_2*M_Chi_Q;
     distr_t M_J_ave= w_J_1*M_J + w_J_2*M_J_Q;
     distr_t E1_ave = w_E1_1*E1 + w_E1_2*E1_Q;
     distr_t M2E1_ave = w_M2E1_1*M2E1 + w_M2E1_2*M2E1_Q;
  
     distr_t D_Chi_ave = w_Chi_1*D_Chi + w_Chi_2*D_Chi_Q;
     distr_t D_J_ave = w_J_1*D_J + w_J_2*D_J_Q;
     distr_t D_E1_ave= w_E1_1*D_E1 + w_E1_2*D_E1_Q;
     distr_t D_M2E1_ave= w_M2E1_1*D_M2E1 + w_M2E1_2*D_M2E1_Q;

     distr_t D2_Chi_ave = w_Chi_1*D2_Chi + w_Chi_2*D2_Chi_Q;
     distr_t D2_J_ave = w_J_1*D2_J + w_J_2*D2_J_Q;
     distr_t D2_E1_ave= w_E1_1*D2_E1 + w_E1_2*D2_E1_Q;
     distr_t D2_M2E1_ave= w_M2E1_1*D2_M2E1 + w_M2E1_2*D2_M2E1_Q;


     distr_t syst_M_E1= E1_Q;
     distr_t syst_D_E1= D_E1_Q;
     distr_t syst_M_a2= M2E1_Q;
     distr_t syst_D_a2= D_M2E1_Q;

  
     //now combine and print the results


     //print fitting functions for MJ, MChi, E1, M2/E1
     distr_t_list MJ_to_print(UseJack);
     distr_t_list MChi_to_print(UseJack);
     distr_t_list E1_to_print(UseJack);
     distr_t_list M2E1_to_print(UseJack);
     Vfloat  a_to_print;

     distr_t_list E1_2nd_fit_to_print(UseJack);
     distr_t_list M2E1_2nd_fit_to_print(UseJack);
     
     //print in step of 0.01 fm;
     double a_max = 0.1;
  
     for(int i=0;i< (int)(a_max/0.001);i++) {
       double x = i*0.001;
       MJ_to_print.distr_list.push_back(  M_J_ave + D_J_ave*pow(x*fm_to_inv_Gev,2) + D2_J_ave*pow(x*fm_to_inv_Gev,4) );
       MChi_to_print.distr_list.push_back(  M_Chi_ave + D_Chi_ave*pow(x*fm_to_inv_Gev,2) + D2_Chi_ave*pow(x*fm_to_inv_Gev,4) );
       E1_to_print.distr_list.push_back( E1 + D_E1*pow(x*fm_to_inv_Gev,2) +   D2_E1*pow(x*fm_to_inv_Gev,4) );
       M2E1_to_print.distr_list.push_back( M2E1 + D_M2E1*pow(x*fm_to_inv_Gev,2) +  D2_M2E1*pow(x*fm_to_inv_Gev,4) );

       E1_2nd_fit_to_print.distr_list.push_back( syst_M_E1 + syst_D_E1*pow(x*fm_to_inv_Gev,2));
       M2E1_2nd_fit_to_print.distr_list.push_back( syst_M_a2 + syst_D_a2*pow(x*fm_to_inv_Gev,2));
       
       a_to_print.push_back(x);
     }


     boost::filesystem::create_directory("../data/chi_c1/Fits");

     Print_To_File({}, {a_to_print, MJ_to_print.ave(), MJ_to_print.err()} , "../data/chi_c1/Fits/MJ.fit", "", "");
     Print_To_File({}, {a_to_print, MChi_to_print.ave(), MChi_to_print.err()} , "../data/chi_c1/Fits/MChi.fit", "", "");
     Print_To_File({}, {a_to_print, E1_to_print.ave(), E1_to_print.err(), E1_2nd_fit_to_print.ave(), E1_2nd_fit_to_print.err() }, "../data/chi_c1/Fits/E1.fit", "", "");
     Print_To_File({}, {a_to_print, M2E1_to_print.ave(), M2E1_to_print.err(), M2E1_2nd_fit_to_print.ave(), M2E1_2nd_fit_to_print.err() }, "../data/chi_c1/Fits/M2E1.fit", "", "");

     //print data

     Print_To_File({} , { (a_list/fm_to_inv_Gev).ave(), MJPSI_list.ave(), MJPSI_list.err(), MCHI_list.ave(), MCHI_list.err(), E1_list.ave(), E1_list.err(), M2_ov_E1_list.ave(), M2_ov_E1_list.err(), M2_list.ave(), M2_list.err() }, "../data/chi_c1/Fits/fit.data", "", "");

     cout<<"Printing final rate: "<<endl;
     distr_t m2ovE1= M2E1_ave; // /SQRT_D( 1.0 -M2E1_ave*M2E1_ave);

     cout<<"E1_ave: "<<E1_ave.ave()<<" +- "<<E1_ave.err()<<endl;
     cout<<"M2/E1_ave: "<<M2E1_ave.ave()<<" +- "<<M2E1_ave.err()<<endl;
     
     distr_t Gamma = 389.4*alpha*(4/27.0)*( E1_ave*E1_ave*( 1.0 + m2ovE1*m2ovE1)); //MeV
     cout<<"G: "<<Gamma.ave()<<" +- "<<Gamma.err()<<endl;

  
     //##################################################################################################


     cout<<"Continuum-extrapolation performed!"<<endl;
     exit(-1);

     return;

     }
