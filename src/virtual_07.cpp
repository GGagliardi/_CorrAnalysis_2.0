#include "../include/virtual_07.h"
#include "Corr_analysis.h"
#include "Spectral.h"
#include "numerics.h"
#include "stat.h"
using namespace std;


bool verbose_lev_07=1;
Vfloat sigmas_07({ 1.0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.5});
//Vfloat sigmas_07({ 2, 3.5});
Vfloat sigmas_07_w0({0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,  0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0,2.25,2.5,2.75,3.0});
int prec_07=128;
const string MODE_FF="TANT";
const bool Skip_spectral_reconstruction_07 = false;
const bool virtuality_scan = false;
const bool Use_preconditioning = true;
const string  preco_tag= (Use_preconditioning)?"prec_":"";
const double Mjpsi= 3.0969; //GeV
const double Mphi= 1.019461; //GeV
const double MDs_phys = 1.96847; // GeV
const double MBs = 5.36692;
const double E0_fact = 0.90;
const bool CONS_EM_CURRENT = true;
const string SM_TYPE = "FF_Exp";
const double QU = -1.0/3.0;
const double QD = -1.0 / 3.0;
const double mh0 = 1.0 / 0.5074431945705498;
const double mh1 = 1.0 / 0.4007281135883006;
const double mh3 = 1.0 / 0.2845905478176752;
const Vfloat masses({mh0, mh1, mh3});


rt_07_Bs Get_virtual_tensor_FF(int n_xg, bool UseJack, int Njacks, string MESON,  string Corr_path,string path_out) {

  Njacks=25;

  Vfloat virtualities;
  for(int i=0;i<75;i++) { virtualities.push_back( 3*i/74.0) ; }
  
  rt_07_Bs return_class;


  Vfloat y_eff({0.67, 0.65, 0.60, 0.55});
  return_class.y_eff= y_eff;



  PrecFloat::setDefaultPrecision(prec_07);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;
  cout<<"Number of xg to analyze: "<<n_xg<<endl;

  int t_07_s, t_07_s_HLT, t_07_c;
 
  string TAG_CURR="";
  if(CONS_EM_CURRENT==false) TAG_CURR="LOC_";

 
  int size_mu_nu= 4;

  double sign_kz=-1.0; //correct one is -1

  //BK
  vector<vector<vector<data_t>>> C_B_d_data(size_mu_nu);
  //T
  vector<vector<vector<data_t>>>  C_T_d_data(size_mu_nu);


  vector<vector<vector<data_t>>> C_B_u_data_std(size_mu_nu), C_B_d_data_std(size_mu_nu);
  vector<vector<vector<data_t>>> C_T_u_data_std(size_mu_nu), C_T_d_data_std(size_mu_nu);



  data_t data_2pts_SM, data_2pts_SMSM;

 
  


  
  for(int mu=0;mu<size_mu_nu;mu++) {

    C_B_d_data[mu].resize(size_mu_nu);
    C_T_d_data[mu].resize(size_mu_nu);

    C_B_u_data_std[mu].resize(size_mu_nu);
    C_B_d_data_std[mu].resize(size_mu_nu);
    C_T_u_data_std[mu].resize(size_mu_nu);
    C_T_d_data_std[mu].resize(size_mu_nu);

   
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_B_d_data[mu][nu].resize(n_xg);
      C_T_d_data[mu][nu].resize(n_xg);

      C_B_u_data_std[mu][nu].resize(n_xg);
      C_B_d_data_std[mu][nu].resize(n_xg);
      C_T_u_data_std[mu][nu].resize(n_xg);
      C_T_d_data_std[mu][nu].resize(n_xg);

    
    }
  }

  //custom sorting of gauge confs
  auto Sort_confs = [](string A, string B) {

			   

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
  
  //read data

  data_2pts_SM.Read(Corr_path+"/spectre", "mes_contr_2pts_SM_3", "P5P5", Sort_confs);
  data_2pts_SMSM.Read(Corr_path+"/spectre", "mes_contr_2pts_SMSM_3", "P5P5", Sort_confs);

  
 
  //loop over mu and nu axial
  vector<pair<int,int>> mu_nu_pair_B({make_pair(1,1),make_pair(2,2)});
  vector<pair<int,int>> mu_nu_pair_T({make_pair(1,2),make_pair(2,1)});
   
  
   
  for(int ixg=0;ixg<n_xg;ixg++) {


   
    //B
    for(auto &pair_B : mu_nu_pair_B) {
      int mu=pair_B.first;
      int nu=pair_B.second;
  
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //B
      //d
      C_B_d_data[mu][nu][ixg].Read(Corr_path+"/spectre", TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      //B
      //u
      C_B_u_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_u_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_B_d_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }

    //T
    for(auto &pair_T : mu_nu_pair_T) {
      int mu=pair_T.first;
      int nu=pair_T.second;
    
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //T
      //d
      C_T_d_data[mu][nu][ixg].Read(Corr_path+"/spectre", TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      //u
      C_T_u_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_u_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_T_d_data_std[mu][nu][ixg].Read(Corr_path+"/standard", TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }
  }

  
   
  GaussianMersenne GM_07(4455);

  //resample RCs
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t ZT_A_boot(0), ZT_B_boot(0), ZT_C_boot(0), ZT_D_boot(0);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");

   

  for(int ijack=0; ijack<Njacks;ijack++) {

   

    
    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM_07()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM_07()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM_07()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM_07()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM_07()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM_07()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM_07()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM_07()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

   

       
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM_07()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM_07()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM_07()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM_07()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    
    
  }

  for(int iboot=0;iboot<1000;iboot++) {

    ZT_A_boot.distr.push_back( L_info_A.ZT_RI2 + GM_07()*L_info_A.ZT_RI2_err);
    ZT_B_boot.distr.push_back( L_info_B.ZT_RI2 + GM_07()*L_info_B.ZT_RI2_err);
    ZT_C_boot.distr.push_back( L_info_C.ZT_RI2 + GM_07()*L_info_C.ZT_RI2_err);
    ZT_D_boot.distr.push_back( L_info_D.ZT_RI2 + GM_07()*L_info_D.ZT_RI2_err);

  }

  cout<<"a_A: "<<a_A.ave()/fmTGeV << " +- "<<a_A.err()/fmTGeV<<endl;
  cout<<"a_B: "<<a_B.ave()/fmTGeV << " +- "<<a_B.err()/fmTGeV<<endl;
  cout<<"a_C: "<<a_C.ave()/fmTGeV << " +- "<<a_C.err()/fmTGeV<<endl;
  cout<<"a_D: "<<a_D.ave()/fmTGeV << " +- "<<a_D.err()/fmTGeV<<endl;

  cout<<"ZT_A : "<<ZT_A.ave()<< " +- "<<ZT_A.err()<<endl;
  cout<<"ZT_B : "<<ZT_B.ave()<< " +- "<<ZT_B.err()<<endl;
  cout<<"ZT_C : "<<ZT_C.ave()<< " +- "<<ZT_C.err()<<endl;
  cout<<"ZT_D : "<<ZT_D.ave()<< " +- "<<ZT_D.err()<<endl;

  cout<<"ZV_A : "<<ZV_A.ave()<< " +- "<<ZV_A.err()<<endl;
  cout<<"ZV_B : "<<ZV_B.ave()<< " +- "<<ZV_B.err()<<endl;
  cout<<"ZV_C : "<<ZV_C.ave()<< " +- "<<ZV_C.err()<<endl;
  cout<<"ZV_D : "<<ZV_D.ave()<< " +- "<<ZV_D.err()<<endl;


  
  
  
 
  int Nens= C_B_d_data[1][1][0].size;
  vector<string> Ens_tags= C_B_d_data[1][1][0].Tag;
  Vint Nts=C_B_d_data[1][1][0].nrows;

  distr_t_list xg_list(UseJack);
  for(int ixg=0;ixg<4;ixg++) xg_list.distr_list.push_back( Get_id_distr(Njacks, UseJack)*(ixg+1)*0.1);
  vector<distr_t_list> F_T_u_list;
  vector<distr_t_list> F_T_d_I_list;
  vector<distr_t_list> F_T_d_sp_I_list;
  vector<distr_t_list> FV_T_u_real_list;
  vector<distr_t_list> FV_T_d_real_list;
  vector<distr_t_list> FA_T_u_real_list;
  vector<distr_t_list> FA_T_d_real_list;

  vector<distr_t_list> FA_T_d_I_real_list;
  vector<distr_t_list> FA_T_d_II_real_list;

  vector<distr_t_list> FV_T_d_I_real_list;
  vector<distr_t_list> FV_T_d_II_real_list;

  vector<distr_t_list> FV_T_d_sp_real_list;
  vector<distr_t_list> FA_T_d_sp_real_list;
  

  vector<vector<distr_t_list>> F_T_d_MB_RE_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_sm_list(n_xg);

  

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_sm_list(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_II_state_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_II_state_sm_list(n_xg);

  vector<vector<distr_t_list>> F_T_d_MB_RE_VMD_III_state_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_MB_IM_VMD_III_state_sm_list(n_xg);

  vector<distr_t_list> F_T_u_VMD_list;
  vector<distr_t_list> F_T_u_VMD_spectre_list;

  vector<vector<double>> sigma_simulated(n_xg);

  for(int ixg=0; ixg<n_xg;ixg++) {

    sigma_simulated[ixg].resize(sigmas_07.size());

    F_T_u_list.emplace_back(UseJack);
    F_T_d_I_list.emplace_back(UseJack);
    F_T_d_sp_I_list.emplace_back(UseJack);
    FV_T_u_real_list.emplace_back(UseJack);
    FV_T_d_real_list.emplace_back(UseJack);
    FA_T_u_real_list.emplace_back(UseJack);
    FA_T_d_real_list.emplace_back(UseJack);

    FA_T_d_I_real_list.emplace_back(UseJack);
    FA_T_d_II_real_list.emplace_back(UseJack);

    FV_T_d_I_real_list.emplace_back(UseJack);
    FV_T_d_II_real_list.emplace_back(UseJack);

    F_T_u_VMD_list.emplace_back(UseJack);
    F_T_u_VMD_spectre_list.emplace_back(UseJack);

    FV_T_d_sp_real_list.emplace_back(UseJack);
    FA_T_d_sp_real_list.emplace_back(UseJack);

    for(int iss=0; iss<(signed)sigmas_07_w0.size(); iss ++) {

   
    F_T_d_MB_RE_VMD_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_sm_list[ixg].emplace_back(UseJack);

    F_T_d_MB_RE_VMD_II_state_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_II_state_sm_list[ixg].emplace_back(UseJack);

    F_T_d_MB_RE_VMD_III_state_sm_list[ixg].emplace_back(UseJack);
    F_T_d_MB_IM_VMD_III_state_sm_list[ixg].emplace_back(UseJack);


    }
    
    for(int is=0; is<(signed)sigmas_07.size(); is++) {
   
      F_T_d_MB_RE_sm_list[ixg].emplace_back(UseJack);
      F_T_d_MB_IM_sm_list[ixg].emplace_back(UseJack);

      
    }
  }




   auto K_RE_distr= [](const distr_t &E, const distr_t &m, double s) -> distr_t {

     distr_t ret(1);
     if(SM_TYPE=="FF_Exp")  {
       for(int ijack=0;ijack<E.size();ijack++) {
	 double x= (E.distr[ijack]-m.distr[ijack]);
	 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)));
	 //ret.distr.push_back(  ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
       }
     }

     
     else if(SM_TYPE=="FF_Gauss") {
       for(int ijack=0;ijack<E.size();ijack++) {
	 double x= (E.distr[ijack]-m.distr[ijack]);
	 ret.distr.push_back(  ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
	 //ret.distr.push_back(sqrt(2)*DawsonF(PrecFloat(x/(sqrt(2.0)*s))).get()/s);
       }
     }
     else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
     
     return ret;
  
  };


   auto K_RE_distr_sub= [](const distr_t &E, const distr_t &m, double s) -> distr_t {

     distr_t ret(1);
     if(SM_TYPE=="FF_Exp")  {
       for(int ijack=0;ijack<E.size();ijack++) {
	 double x= (E.distr[ijack]-m.distr[ijack]);
	 double y= E.distr[ijack];
	 ret.distr.push_back( 2*sinh(x)*cos(s)/(cosh(2*x) - cos(2*s)) - 2*sinh(y)*cos(0)/(cosh(2*y) - cos(2*0))  );
	 //ret.distr.push_back(  ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
       }
     }

     
     else if(SM_TYPE=="FF_Gauss") {
       for(int ijack=0;ijack<E.size();ijack++) {
	 double x= (E.distr[ijack]-m.distr[ijack]);
	 ret.distr.push_back(  ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
	 //ret.distr.push_back(sqrt(2)*DawsonF(PrecFloat(x/(sqrt(2.0)*s))).get()/s);
       }
     }
     else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");
     
     return ret;
  
  };

   auto K_IM_distr= [](const distr_t &E, const distr_t &m, double s) -> distr_t {


       distr_t ret(1);
       if(SM_TYPE=="FF_Exp") {
	 for(int ijack=0;ijack<E.size();ijack++) {
	   double x= (E.distr[ijack]-m.distr[ijack]);
	   ret.distr.push_back ( 2*cosh(x)*sin(s)/(cosh(2*x) - cos(2*s)) );
	   //ret.distr.push_back(exp(-x)*sin(s)/( 1 + exp(-2*x) -2*cos(s)*exp(-x)));
	 }
       }
       else if(SM_TYPE=="FF_Gauss") {
	 for(int ijack=0;ijack<E.size();ijack++) {
	   ret.distr.push_back( M_PI*Get_exact_gauss( E.distr[ijack], m.distr[ijack],  s, 0.0 ));
	 }
       }
       else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");

        
       return ret;
       
  };

 

  auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {


    if(SM_TYPE=="FF_Gauss") {
      PrecFloat x = (E-m)/(sqrt(PrecFloat(2))*s);
      return  cos(s/2)/(cosh(x)/sinh(x/2) - cos(s)/sinh(x/2));
      //return ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
      //return sqrt(2)*DawsonF(x)/s;
    }


    if(SM_TYPE=="FF_Gauss_Schwartz") {
      
      PrecFloat x = (E-m);
      if( abs(x)/s > 0.1) return ( 1 - exp( -x*x/(2*s*s)))/x;
      else {
        int n=2;
        bool converged=false;
        PrecFloat fact = -pow(x/(sqrt(PrecFloat(2))*s),2);
        PrecFloat sum= x/(2*s*s);
        PrecFloat M= sum;
        PrecFloat prec_sum;
        while(!converged) {
          M *= fact/n;
          prec_sum=sum;
          sum += M;
          n++;
          if(prec_sum==sum) converged=true;
        }
        return sum;
      }
      
    }
    
    
    if(SM_TYPE=="FF_Cauchy") {
      PrecFloat t = (E-m);
      return t/( t*t + s*s);
    }

    if(SM_TYPE=="FF_Exp") {
      
      PrecFloat x= (E-m);

      PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 

      return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
      
      //return ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
    }
    
    if(SM_TYPE=="FF_Sinh_half") {

      PrecFloat x = (E-m);

      return exp(-x/2)*2*sinh(x/2)/( pow(2*sinh(x/2),2) + s*s);

    }

    PrecFloat norm;
    if( s > 1) norm= PrecFloat(2)*log( s + sqrt( s*s -1))/sqrt(s*s -1);
    else if(s==1) norm=PrecFloat(2);
    else {
      PrecFloat phi= abs(atan( sqrt( 1 - s*s)/s));
      norm= PrecFloat(2)*phi/sqrt( 1 - s*s);
    }
    norm /= precPi();
    norm = 1/norm;
    
       
    PrecFloat t = (E-m);
    PrecFloat x=sinh(t);
    PrecFloat res= (x + (s*s/x));
    res=1/res;

    if( abs(t) >= 1) return norm*res;
    else return norm*x/( s*s + x*x);
      
    exit(-1);
    return 0;
  };

  //smeared kernel of the immaginary part
  auto K_IM = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

    if(SM_TYPE=="FF_Gauss") return precPi()*Get_exact_gauss(E, m, s, E0);
    //return precPi()*Get_exact_gauss(E, m, s, E0);
    else if(SM_TYPE=="FF_Gauss_Schwartz") return precPi()*Get_exact_gauss(E, m, s, E0);
    else if(SM_TYPE=="FF_Cauchy") {
      PrecFloat t= (E-m);
      return s/( t*t + s*s);
    }

    else if(SM_TYPE=="FF_Exp") {

      PrecFloat x= (E-m);

      //PrecFloat phi= abs(atan( sin(s)/(1-cos(s))));
      
      //PrecFloat norm= PrecFloat(2)*phi/precPi() ;

      PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 

      return 2*sin(s)/(cosh_ov_cosh_half - cos(2*s)/cosh(x));

      //return (exp(-x)*sin(s))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
      
    }
    
    else if(SM_TYPE=="FF_Sinh_half") {
      
      PrecFloat x = (E-m);
      
      return exp(-x/2)*s/( pow(2*sinh(x/2),2) + s*s);

    }
    

    return PrecFloat(0);
    
  };
  
  
  boost::filesystem::create_directory( path_out+"/corr_2pts");
  boost::filesystem::create_directory( path_out+"/corr_3pts");
  boost::filesystem::create_directory( path_out+"/mass");
  boost::filesystem::create_directory( path_out+"/covariance");
  boost::filesystem::create_directory( path_out+"/FF_d_I");
  boost::filesystem::create_directory( path_out+"/FF_d_II");
  boost::filesystem::create_directory( path_out+"/FF_u");
  boost::filesystem::create_directory( path_out+"/FF_d");
  boost::filesystem::create_directory( path_out+"/FF");


  
  
  //loop over ensembles
  for(int iens=0; iens<Nens;iens++) {

    double Mp_cont;
    
    if(MESON=="B0s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=14; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=22; t_07_c=35;    } //t_07_s_HLT_std=22
      else crash("B0s-crash-virtual");
      Mp_cont = masses[0];
    }
    else if(MESON=="B1s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=14; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=22; t_07_c=35;    }
      else crash("B1s-crash-virtual");
      Mp_cont = masses[1];
    }
    else if(MESON=="B3s") {
      if(Ens_tags[iens] == "cB211b.072.64") {  t_07_s=25; t_07_s_HLT=10; t_07_c=25;    }
      else if(Ens_tags[iens] == "cD211a.054.96") {  t_07_s=35; t_07_s_HLT=21; t_07_c=35;    }
      else crash("B3s-crash-virtual");
      Mp_cont = masses[2];
    }
    else crash("Meson: "+MESON+" not yet simulated");

    boost::filesystem::create_directory( path_out+"/corr");
    boost::filesystem::create_directory( path_out+"/mass/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_I/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_II/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d_II/"+data_2pts_SM.Tag[iens]+"/VMD_virt_scan");
    boost::filesystem::create_directory( path_out+"/FF_u/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF_d/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/FF/"+data_2pts_SM.Tag[iens]);
    
    CorrAnalysis Corr(UseJack,Njacks,1000);
    Corr.Nt = Nts[iens];
    CorrAnalysis Corr_boot(false,Njacks, 1000);
    Corr_boot.Nt= Nts[iens];
   

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts_SM.Tag[iens]);
   

    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d;

    thetas= Read_From_File(Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 1 , 5);
    masses_u= Read_From_File(Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File( Corr_path+"/spectre/"+Ens_tags[iens]+"/pars_list.dat", 4 , 5);

    double mu= masses_u[0];
    double md= masses_d[0];


    //RCs
    
    distr_t ZT, a_distr, ZV;
    distr_t ZT_boot;
    if(data_2pts_SM.Tag[iens].substr(1,1)=="A") { ZT=ZT_A; ZV=ZV_A; ZT_boot = ZT_A_boot; a_distr=a_A;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="B") { ZT= ZT_B; ZV=ZV_B; ZT_boot = ZT_B_boot; a_distr=a_B;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="C") { ZT= ZT_C; ZV=ZV_C; ZT_boot = ZT_C_boot; a_distr=a_C;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="D") { ZT= ZT_D; ZV=ZV_D; ZT_boot = ZT_D_boot; a_distr=a_D;}
    else crash("Ensemble: "+data_2pts_SM.Tag[iens]+" not recognised");


    if(CONS_EM_CURRENT==false) ZT = ZT*ZV;


    //2pts plateaux
    if(data_2pts_SM.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=14; Corr.Tmax= 24;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 24;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=20;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=19;}
      else { Corr.Tmin=14; Corr.Tmax=27;}
      
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=24; Corr.Tmax=38;}
      else if(MESON == "B2s") {Corr.Tmin=22; Corr.Tmax=34;}
      else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=23;}
      else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=22;}
      else {Corr.Tmin=23; Corr.Tmax=33;}
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON == "B2s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=31;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=30;}
      else {Corr.Tmin=29; Corr.Tmax=42;}
    }
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=32; Corr.Tmax=53;}
      else if(MESON == "B2s") { Corr.Tmin=28; Corr.Tmax= 37;}
      else if(MESON=="B3s") {Corr.Tmin=28; Corr.Tmax=37;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=37;}
      else {Corr.Tmin=38; Corr.Tmax=52;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts_SM.Tag[iens]+" not recognized");

    Corr_boot.Tmin= Corr.Tmin; Corr_boot.Tmax=Corr.Tmax;


    //mass and decay constants
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, path_out+"/masses/"+data_2pts_SM.Tag[iens]+"/eff_mass_SM.dat");
    distr_t_list pt2_distr_SMSM= Corr.corr_t(data_2pts_SMSM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SMSM.dat");
    distr_t_list eff_mass_SMSM = Corr.effective_mass_t(pt2_distr_SMSM, path_out+"/masses/"+data_2pts_SM.Tag[iens]+"/eff_mass_SMSM.dat");
    distr_t mel_SMSM= Corr.Fit_distr(Corr.mel_ov_mass_t(pt2_distr_SMSM, ""))/2.0;

    distr_t_list pt2_distr_SM_boot= Corr_boot.corr_t(data_2pts_SM.col(0)[iens], "");
    distr_t_list eff_mass_SM_boot = Corr_boot.effective_mass_t(pt2_distr_SM_boot, "");
    distr_t_list pt2_distr_SMSM_boot= Corr_boot.corr_t(data_2pts_SMSM.col(0)[iens], "");
    distr_t_list eff_mass_SMSM_boot = Corr_boot.effective_mass_t(pt2_distr_SMSM_boot,"");


    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    distr_t M_P_SMSM = Corr.Fit_distr(eff_mass_SMSM);
    distr_t M_P= 0.5*(M_P_SM+M_P_SMSM);
    auto SINH= [](double x) { return sinh(x);};
    distr_t_list FP_SM_distr_list = (mu + md )*Corr.residue_t( pt2_distr_SM, "")/(M_P*distr_t::f_of_distr(SINH, M_P)*Corr.matrix_element_t(pt2_distr_SMSM, ""));
    Print_To_File({}, {FP_SM_distr_list.ave(), FP_SM_distr_list.err()}, path_out+"/decay_const/"+data_2pts_SM.Tag[iens]+"/decay_const_SM.dat.t", "", "");
    distr_t FP_SM= Corr.Fit_distr( FP_SM_distr_list  );

    distr_t M_P_boot = 0.5*( Corr_boot.Fit_distr( eff_mass_SMSM_boot) + Corr_boot.Fit_distr(eff_mass_SM_boot));

  
    distr_t mel_SMSM_boot= Corr_boot.Fit_distr(Corr_boot.mel_ov_mass_t(pt2_distr_SMSM_boot, ""))/2.0;
  


    cout<<"MP: "<<(M_P/a_distr).ave()<<" +- "<<(M_P/a_distr).err()<<endl;
    cout<<"FP_SM: "<<(FP_SM/a_distr).ave()<<" +- "<<(FP_SM/a_distr).err()<<endl;

    
    distr_t F_P= FP_SM;

    //loop over photon momentum
    for(int ixg=0;ixg<n_xg;ixg++) {

      //get xg, Eg, kz from thetas
      double theta=thetas[ixg];
      pt3_momenta pt3_mom_07(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], 0.0, L_info.L, L_info.T);
      double Eg= pt3_mom_07.Egamma();
      distr_t Eg_off = M_P - Eg;
      double kz = pt3_mom_07.k()[2];
      distr_t xg= pt3_mom_07.x_gamma(M_P);
      distr_t xg_boot= pt3_mom_07.x_gamma(M_P_boot);
      distr_t Eg_MB_off = MBs*a_distr*(1.0 -xg/2.0);
      //double y_eff= y*(1-xg.ave()/2.0)/(1-0.1/2.0);
     
      distr_t Eg_MB_off_new= (y_eff[ixg]*Mp_cont + (1-y_eff[ixg])*MBs)*(1.0 -0.5*xg);
      Eg_MB_off = y_eff[ixg]*Eg_off + (1-y_eff[ixg])*Eg_MB_off;
      //xg_list.distr_list.push_back(xg);

      
    
      cout<<"##### Considering kinematic with..."<<endl;
      cout<<"Eg (on-shell): "<<Eg<<endl;
      cout<<"Eg (virt): "<<Eg_off.ave()<<" +- "<<Eg_off.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      int Im_Re;

      //tensor- electric part
      Corr.Reflection_sign= 1;
      Im_Re=1;
      Corr.Perform_Nt_t_average = 0;
      Corr_boot.Reflection_sign=1;
      Corr_boot.Perform_Nt_t_average=0;
	 
	 
      distr_t_list T_u_std = 0.5*QU*Corr.corr_t(summ_master(C_T_u_data_std[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_u_data_std[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_u_xg_"+to_string(ixg));
      	 
      distr_t_list T_d_std = 0.5*QD*Corr.corr_t(summ_master(C_T_d_data_std[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data_std[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_d_xg_"+to_string(ixg));
      
      distr_t_list T_d = 0.5*QD*Corr.corr_t(summ_master(C_T_d_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_T_d_sp_xg_"+to_string(ixg));
      distr_t_list T_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_T_d_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)),"");

      Im_Re=0;

      distr_t_list B_u_std = 0.5*QU*Corr.corr_t(summ_master(C_B_u_data_std[1][1][ixg].col(Im_Re)[iens], C_B_u_data_std[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_u_xg_"+to_string(ixg));
      distr_t_list B_d_std = 0.5*QD*Corr.corr_t(summ_master(C_B_d_data_std[1][1][ixg].col(Im_Re)[iens], C_B_d_data_std[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_d_xg_"+to_string(ixg));
      distr_t_list B_d = 0.5*QD*Corr.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+TAG_CURR+""+data_2pts_SM.Tag[iens]+"_B_d_sp_xg_"+to_string(ixg));
      distr_t_list B_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]),"");


     
      //determine the form factors that do not need spectral-reconstruction techniques, i.e. up-type quark (heavy) contribution and down-type quark (s-quark) contribution in I-TO


      //#############     VMD PRED  FT_u  #####################
      //define corr for Fu in 2nd TO for VMD pred
      distr_t_list Corr_Tu_2TO(UseJack);
      for(int t=t_07_c; t <=Corr.Nt/2;t++) {
	Corr_Tu_2TO.distr_list.push_back( (xg/2.0)*(sign_kz*T_u_std.distr_list[t] + B_u_std.distr_list[t]));
      }
      CorrAnalysis Corr_VMD_anal(UseJack, Njacks,1000);
      if(MESON=="B0s") {
	
	if(Ens_tags[iens] == "cB211b.072.64") {
	  Corr_VMD_anal.Tmin= 17;
	  Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=25;
	  Corr_VMD_anal.Tmax=46;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else if(MESON=="B1s") {
	if(Ens_tags[iens] == "cB211b.072.64") {
	Corr_VMD_anal.Tmin= 17;
	Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=27;
	  Corr_VMD_anal.Tmax=50;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else if(MESON=="B3s") {
	if(Ens_tags[iens] == "cB211b.072.64") {
	  Corr_VMD_anal.Tmin= 17;
	  Corr_VMD_anal.Tmax= 30;
	}
	else if(Ens_tags[iens] == "cD211a.054.96") {
	  Corr_VMD_anal.Tmin=25;
	  Corr_VMD_anal.Tmax=50;
	}
	else crash("Ensemble: "+Ens_tags[iens]+" not yet implemented");
      }
      else crash("Meson: "+MESON+" not yet implemented");
      
      Corr_VMD_anal.Nt=2*(Corr_Tu_2TO.size()-1);
      distr_t_list Corr_Tu_2TO_symm= Corr_Tu_2TO;
      for(int t=Corr_VMD_anal.Nt/2 +1;t<Corr_VMD_anal.Nt;t++) Corr_Tu_2TO_symm.distr_list.push_back( Corr_Tu_2TO.distr_list[Corr_VMD_anal.Nt -t]);
      distr_t_list eff_M_Tu_distr= Corr_VMD_anal.effective_mass_t(Corr_Tu_2TO_symm, "" );
      distr_t eff_M_Tu= Corr_VMD_anal.Fit_distr(eff_M_Tu_distr)/a_distr;
      Print_To_File({}, {eff_M_Tu_distr.ave(), eff_M_Tu_distr.err(), (eff_M_Tu_distr-eff_M_Tu_distr[30]).ave(),  (eff_M_Tu_distr-eff_M_Tu_distr[30]).err() }, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vb_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      distr_t Ups_MT = Corr_VMD_anal.Fit_distr( Corr_Tu_2TO_symm*(EXPT_DL(eff_M_Tu_distr)));
      F_T_u_VMD_spectre_list[ixg].distr_list.push_back( Ups_MT*ZT*(1.0/(mel_SMSM*Eg))*EXP_D( M_P*t_07_c)*K_RE_distr_sub( eff_M_Tu*a_distr, Eg_off, 1e-5));
      distr_t_list Corr_F_T_u_VMD = Ups_MT*EXPT_D(-1.0*eff_M_Tu*a_distr, Corr_Tu_2TO.size());
      distr_t F_T_u_VMD_distr(UseJack, Njacks);
      for(int ty=1;ty<Corr_Tu_2TO.size(); ty++) {
	F_T_u_VMD_distr = F_T_u_VMD_distr + Corr_F_T_u_VMD[ty]*EXP_D(ty*Eg_off);
      }
      F_T_u_VMD_list[ixg].distr_list.push_back( F_T_u_VMD_distr*ZT*(1.0/(mel_SMSM*Eg))*EXP_D( M_P*t_07_c));


      
      //#######################################################
      
      
      distr_t F_T_u(UseJack, Njacks);
      distr_t F_T_d_I(UseJack, Njacks);
      distr_t F_T_d_I_sp(UseJack, Njacks);
      distr_t FV_T_u_real(UseJack, Njacks);
      distr_t FV_T_d_real(UseJack, Njacks);
      distr_t FA_T_u_real(UseJack, Njacks);
      distr_t FA_T_d_real(UseJack, Njacks);
      distr_t FV_T_d_sp_real(UseJack,Njacks);
      distr_t FA_T_d_sp_real(UseJack,Njacks);

      distr_t FA_T_d_I_real(UseJack, Njacks);
      distr_t FA_T_d_II_real(UseJack, Njacks);

      distr_t FV_T_d_I_real(UseJack, Njacks);
      distr_t FV_T_d_II_real(UseJack, Njacks);

      distr_t_list F_T_d_ty(UseJack), F_T_d_sp_ty(UseJack), F_T_u_ty(UseJack);
      distr_t_list F_T_d_ty_psum(UseJack), F_T_d_sp_ty_psum(UseJack), F_T_u_ty_psum(UseJack);

      distr_t_list FA_T_d_real_ty(UseJack), FA_T_d_sp_real_ty(UseJack), FA_T_u_real_ty(UseJack);
      distr_t_list FA_T_d_real_ty_psum(UseJack), FA_T_d_sp_real_ty_psum(UseJack), FA_T_u_real_ty_psum(UseJack);

      distr_t_list FV_T_d_real_ty(UseJack), FV_T_d_sp_real_ty(UseJack), FV_T_u_real_ty(UseJack);
      distr_t_list FV_T_d_real_ty_psum(UseJack), FV_T_d_sp_real_ty_psum(UseJack), FV_T_u_real_ty_psum(UseJack);

      //reminder of

      distr_t reminder_FV_T_d_I_real(UseJack,Njacks), reminder_FA_T_d_I_real(UseJack,Njacks) , reminder_F_T_d_I(UseJack,Njacks), reminder_F_T_d_I_sp(UseJack,Njacks); 

      
      auto HeavyTheta=[](const int x) {  return ((x>=0)+(x>0))/2.0;   };
      auto Exp= [&Njacks, &UseJack](const distr_t &A) -> distr_t { distr_t ret(UseJack); for(int ijack=0;ijack<Njacks;ijack++) ret.distr.push_back( exp(A.distr[ijack])); return ret;};
      double T=Corr.Nt;
      
      for(int ty=0; ty <= Corr.Nt/2; ty++) {
	
	const distr_t f1=  Exp(-(T/2-ty)*Eg_off);

	const distr_t f2= Exp(-((3*T/2)-ty)*Eg_off);

	const double f1_real= exp(-(T/2-ty)*Eg);
	const double f2_real= exp(-((3*T/2)-ty)*Eg);

	const double h1=HeavyTheta((T/2)-ty);
	
	const double h2=HeavyTheta(ty-(T/2));

	const distr_t f1_sub= f1 - Exp( -1.0*Eg_off*abs(T/2 - t_07_s));
	const distr_t f2_sub= f2 - Exp( -1.0*Eg_off*abs(T/2 - t_07_s));

	const distr_t f1_sub_sp= f1 - Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));
	const distr_t f2_sub_sp= f2 - Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));

	const double f1_real_sub= f1_real - exp( -1.0*Eg*abs(T/2 - t_07_s));
	const double f2_real_sub= f2_real - exp( -1.0*Eg*abs(T/2 - t_07_s));



	F_T_d_ty.distr_list.push_back( (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2)            );
	F_T_d_sp_ty.distr_list.push_back( (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1*f1+h2*f2)   );
	F_T_u_ty.distr_list.push_back(  (sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2));


	F_T_d_ty_psum.distr_list.push_back( ((ty==0)?F_T_d_ty[ty]:(F_T_d_ty_psum[ty-1] + F_T_d_ty[ty])));
	F_T_d_sp_ty_psum.distr_list.push_back( ((ty==0)?F_T_d_sp_ty[ty]:(F_T_d_sp_ty_psum[ty-1] + F_T_d_sp_ty[ty])));
	F_T_u_ty_psum.distr_list.push_back( ((ty==0)?F_T_u_ty[ty]:(F_T_u_ty_psum[ty-1] + F_T_u_ty[ty])));



	FA_T_d_real_ty.distr_list.push_back( ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real)    )      ;

        FA_T_d_sp_real_ty.distr_list.push_back( ( B_d.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FA_T_u_real_ty.distr_list.push_back(  ( B_u_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));



        FA_T_d_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_d_real_ty[ty]:(FA_T_d_real_ty_psum[ty-1] + FA_T_d_real_ty[ty])));
	FA_T_d_sp_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_d_sp_real_ty[ty]:(FA_T_d_sp_real_ty_psum[ty-1] + FA_T_d_sp_real_ty[ty])));
	FA_T_u_real_ty_psum.distr_list.push_back( ((ty==0)?FA_T_u_real_ty[ty]:(FA_T_u_real_ty_psum[ty-1] + FA_T_u_real_ty[ty])));



	FV_T_d_real_ty.distr_list.push_back( (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FV_T_d_sp_real_ty.distr_list.push_back( (sign_kz*T_d.distr_list[ty]*(1.0- xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));
	FV_T_u_real_ty.distr_list.push_back(  (sign_kz*T_u_std.distr_list[ty]*(1.0- xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real));


	FV_T_d_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_d_real_ty[ty]:(FV_T_d_real_ty_psum[ty-1] + FV_T_d_real_ty[ty])));
	FV_T_d_sp_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_d_sp_real_ty[ty]:(FV_T_d_sp_real_ty_psum[ty-1] + FV_T_d_sp_real_ty[ty])));
	FV_T_u_real_ty_psum.distr_list.push_back( ((ty==0)?FV_T_u_real_ty[ty]:(FV_T_u_real_ty_psum[ty-1] + FV_T_u_real_ty[ty])));

	
	F_T_u = F_T_u +  (sign_kz*T_u_std.distr_list[ty]*(xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	FV_T_u_real = FV_T_u_real + (sign_kz*T_u_std.distr_list[ty]*(1.0- xg/2.0)   + B_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FV_T_d_real = FV_T_d_real + (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_u_real = FA_T_u_real + ( B_u_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_u_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_d_real = FA_T_d_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_d_sp_real = FA_T_d_sp_real + ( B_d.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FV_T_d_sp_real = FV_T_d_sp_real + (sign_kz*T_d.distr_list[ty]*(1.0- xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);

	if(ty<=t_07_s) {
	FA_T_d_I_real = FA_T_d_I_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	FV_T_d_I_real = FV_T_d_I_real +  (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	}
	else {
	  FA_T_d_II_real = FA_T_d_II_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	  FV_T_d_II_real = FV_T_d_II_real +  (sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_real_sub+ h2*f2_real_sub);
	}
	
	
	if(ty<= t_07_s) {  //SHOULD BE <= t_07_s
	  F_T_d_I = F_T_d_I + (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1*f1_sub+ h2*f2_sub);
	}
	if(ty <= t_07_s_HLT) {
	  F_T_d_I_sp = F_T_d_I_sp + (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1*f1_sub_sp+h2*f2_sub_sp);
	}

	//collect reminders
	reminder_FA_T_d_I_real = reminder_FA_T_d_I_real + ( B_d_std.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d_std.distr_list[ty]*xg/2.0 )*(h1 + h2)*exp( -Eg*abs(T/2 - t_07_s));
	reminder_FV_T_d_I_real = reminder_FV_T_d_I_real + ( sign_kz*T_d_std.distr_list[ty]*(1.0- xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1 + h2)*exp( -Eg*abs(T/2 - t_07_s));

	reminder_F_T_d_I= reminder_F_T_d_I +  (sign_kz*T_d_std.distr_list[ty]*(xg/2.0)   + B_d_std.distr_list[ty]*xg/2.0 )*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s));
	reminder_F_T_d_I_sp= reminder_F_T_d_I_sp +   (sign_kz*T_d.distr_list[ty]*(xg/2.0) + B_d.distr_list[ty]*xg/2.0)*(h1+h2)*Exp( -1.0*Eg_off*abs(T/2 - t_07_s_HLT));
	
      }

      //sum reminders 
      FA_T_d_I_real = FA_T_d_I_real + reminder_FA_T_d_I_real;
      FV_T_d_I_real = FV_T_d_I_real + reminder_FV_T_d_I_real;
      F_T_d_I = F_T_d_I + reminder_F_T_d_I;
      F_T_d_I_sp = F_T_d_I_sp + reminder_F_T_d_I_sp;

     
      
         
      //normalize
      F_T_u = F_T_u*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      F_T_d_I = F_T_d_I*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FV_T_u_real= FV_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FV_T_d_real= FV_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_u_real= FA_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FA_T_d_real= FA_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_d_I_real= FA_T_d_I_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_d_II_real= FA_T_d_II_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;

      FV_T_d_I_real = FV_T_d_I_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FV_T_d_II_real= FV_T_d_II_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;

      FV_T_d_sp_real= FV_T_d_sp_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT) ;
      FA_T_d_sp_real= FA_T_d_sp_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT) ;
      F_T_d_I_sp = F_T_d_I_sp*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);


      F_T_d_ty = F_T_d_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      F_T_d_sp_ty = F_T_d_sp_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      F_T_u_ty= F_T_u_ty*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      F_T_d_ty_psum = F_T_d_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      F_T_d_sp_ty_psum = F_T_d_sp_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      F_T_u_ty_psum= F_T_u_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);


      FA_T_d_real_ty = FA_T_d_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FA_T_d_sp_real_ty = FA_T_d_sp_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FA_T_u_real_ty= FA_T_u_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FA_T_d_real_ty_psum = FA_T_d_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FA_T_d_sp_real_ty_psum = FA_T_d_sp_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FA_T_u_real_ty_psum= FA_T_u_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FV_T_d_real_ty = FV_T_d_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FV_T_d_sp_real_ty = FV_T_d_sp_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FV_T_u_real_ty= FV_T_u_real_ty*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

      FV_T_d_real_ty_psum = FV_T_d_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s);
      FV_T_d_sp_real_ty_psum = FV_T_d_sp_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s_HLT))*Exp( M_P*t_07_s_HLT);
      FV_T_u_real_ty_psum= FV_T_u_real_ty_psum*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c);

           
      //push_back
      F_T_u_list[ixg].distr_list.push_back( F_T_u);
      F_T_d_I_list[ixg].distr_list.push_back( F_T_d_I);
      F_T_d_sp_I_list[ixg].distr_list.push_back( F_T_d_I_sp);
     
      FV_T_u_real_list[ixg].distr_list.push_back( FV_T_u_real);
      FV_T_d_real_list[ixg].distr_list.push_back( FV_T_d_real);
      FA_T_u_real_list[ixg].distr_list.push_back( FA_T_u_real);
      FA_T_d_real_list[ixg].distr_list.push_back( FA_T_d_real);
      FA_T_d_I_real_list[ixg].distr_list.push_back( FA_T_d_I_real);
      FA_T_d_II_real_list[ixg].distr_list.push_back( FA_T_d_II_real);
      FV_T_d_I_real_list[ixg].distr_list.push_back( FV_T_d_I_real);
      FV_T_d_II_real_list[ixg].distr_list.push_back( FV_T_d_II_real);

      FV_T_d_sp_real_list[ixg].distr_list.push_back( FV_T_d_sp_real);
      FA_T_d_sp_real_list[ixg].distr_list.push_back( FA_T_d_sp_real);

      //print to File ty analysis
      Print_To_File( { }, { F_T_d_ty.ave(), F_T_d_ty.err(), F_T_d_ty_psum.ave(), F_T_d_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { F_T_d_sp_ty.ave(), F_T_d_sp_ty.err(), F_T_d_sp_ty_psum.ave(), F_T_d_sp_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { F_T_u_ty.ave(), F_T_u_ty.err(),  F_T_u_ty_psum.ave(), F_T_u_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      Print_To_File( { }, { FA_T_d_real_ty.ave(), FA_T_d_real_ty.err(), FA_T_d_real_ty_psum.ave(), FA_T_d_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FA_T_d_sp_real_ty.ave(), FA_T_d_sp_real_ty.err(), FA_T_d_sp_real_ty_psum.ave(), FA_T_d_sp_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FA_T_u_real_ty.ave(), FA_T_u_real_ty.err(),  FA_T_u_real_ty_psum.ave(), FA_T_u_real_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      Print_To_File( { }, { FV_T_d_real_ty.ave(), FV_T_d_real_ty.err(), FV_T_d_real_ty_psum.ave(), FV_T_d_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FV_T_d_sp_real_ty.ave(), FV_T_d_sp_real_ty.err(), FV_T_d_sp_real_ty_psum.ave(), FV_T_d_sp_real_ty_psum.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_sp_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { FV_T_u_real_ty.ave(), FV_T_u_real_ty.err(),  FV_T_u_real_ty_psum.ave(), FV_T_u_real_ty_psum.err() }, path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");


      Print_To_File( { }, { (FA_T_d_real_ty+FA_T_u_real_ty).ave(), (FA_T_d_real_ty+FA_T_u_real_ty).err(), (FA_T_d_real_ty_psum+FA_T_u_real_ty_psum).ave(), (FA_T_d_real_ty_psum+FA_T_u_real_ty).err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"TA_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
      Print_To_File( { }, { (FV_T_d_real_ty+FV_T_u_real_ty).ave(), (FV_T_d_real_ty+FV_T_u_real_ty).err(), (FV_T_d_real_ty_psum+FV_T_u_real_ty_psum).ave(), (FV_T_d_real_ty_psum+FV_T_u_real_ty).err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"TV_real_ty_analysis_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
        
    
      
   

   
            

      //define corrs to be used in spec dens reconstruction
      distr_t_list Corr_T(UseJack);
      distr_t_list Corr_T_boot(0);

      distr_t_list Corr_TV_real(UseJack), Corr_TV_real_boot(0);

      for(int t=t_07_s_HLT; t <= Corr.Nt/2; t++) {
	
	Corr_T.distr_list.push_back( sign_kz*T_d.distr_list[t] + B_d.distr_list[t]);
	Corr_T_boot.distr_list.push_back( sign_kz*T_d_boot.distr_list[t] + B_d_boot.distr_list[t]);

	Corr_TV_real.distr_list.push_back(   sign_kz*T_d.distr_list[t]*(1 - xg/2.0) + B_d.distr_list[t]*(xg/2.0));
	Corr_TV_real_boot.distr_list.push_back(   sign_kz*T_d_boot.distr_list[t]*(1 - xg_boot/2.0) + B_d_boot.distr_list[t]*(xg_boot/2.0));
      }

  

      cout<<"ZT: "<<ZT_boot.ave()<<" +- "<<ZT_boot.err()<<" [boot], "<<ZT.ave()<<" +- "<<ZT.err()<<" [jack]"<<endl;
      cout<<"mel: "<<mel_SMSM_boot.ave()<<" +- "<<mel_SMSM_boot.err()<<" [boot], "<<mel_SMSM.ave()<<" +- "<<mel_SMSM.err()<<" [jack]"<<endl;
      cout<<"M_P: "<<M_P_boot.ave()<<" +- "<<M_P_boot.err()<<" [boot], "<<M_P.ave()<<" +- "<<M_P.err()<<" [jack]"<<endl;
      cout<<"xg: "<<xg_boot.ave()<<" +- "<<xg_boot.err()<<" [boot], "<<xg.ave()<<" +- "<<xg.err()<<" [jack]"<<endl;
      cout<<"Exp: "<<EXP_D(M_P*t_07_s_HLT).ave()<<" +- "<<EXP_D(M_P*t_07_s_HLT).err()<<" [boot], "<<Exp(M_P*t_07_s_HLT).ave()<<" +- "<<Exp(M_P*t_07_s_HLT).err()<<" [jack]"<<endl;

      
      
      cout<<"Exp(m*t)_boot: "<<EXP_D(M_P_boot*t_07_s_HLT).size()<<endl;

      //rescale error of corr_boot

    
      distr_t FACT_boot= ZT_boot*(1.0/(mel_SMSM_boot*Eg))*EXP_D( M_P_boot*t_07_s_HLT)*xg_boot/2.0;
      distr_t_list Corr_T_boot_to_print=FACT_boot*Corr_T_boot;
      distr_t FACT= ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT)*xg/2.0;
      distr_t FACT_real = ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT);

      for(int t=t_07_s_HLT; t<= Corr.Nt/2; t++) {

	Corr_T_boot_to_print.distr_list[t-t_07_s_HLT] = Corr_T_boot_to_print.ave(t-t_07_s_HLT) + (Corr_T_boot_to_print.distr_list[t-t_07_s_HLT] - Corr_T_boot_to_print.ave(t-t_07_s_HLT))*(FACT*Corr_T).err(t-t_07_s_HLT)/Corr_T_boot_to_print.err(t-t_07_s_HLT);

	Corr_T_boot.distr_list[t-t_07_s_HLT] =  Corr_T_boot.ave(t-t_07_s_HLT) + (Corr_T_boot.distr_list[t-t_07_s_HLT] - Corr_T_boot.ave(t-t_07_s_HLT))*(Corr_T).err(t-t_07_s_HLT)/Corr_T_boot.err(t-t_07_s_HLT);
	
      }
      
      
      

   

      //print correlator to file
      boost::filesystem::create_directory("../Nazario_data_C_HLT");
      ofstream Print_boot("../Nazario_data_C_HLT/C_mh_"+to_string_with_precision(M_P.ave(), 6)+"_"+to_string_with_precision(xg_list.ave(ixg),2)+"_"+Ens_tags[iens]+".boot");
      Print_boot.precision(10);
      Print_boot<<"#Nconfs=1000     T="<<Corr.Nt/2 - t_07_s_HLT<<endl;
      for(int iboot=0;iboot<1000;iboot++) {
	for(int t=t_07_s_HLT;t<=Corr.Nt/2; t++) Print_boot<<Corr_T_boot_to_print.distr_list[t-t_07_s_HLT].distr[iboot]<<endl ;
      }
      Print_boot.close();

      int tmax= Corr_T.size();

     
      Print_To_File({}, { (Corr_T*FACT).ave(),  (Corr_T*FACT).err(),  (Corr_T_boot_to_print).ave(),  (Corr_T_boot_to_print).err()},  path_out+"/corr_2pts/HLT_boot_jack_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "",""); 

      
      CorrAnalysis Corr_HLT(UseJack, Njacks,1000);
      CorrAnalysis Corr_HLT_boot(0, Njacks, 1000);
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 14;
	    Corr_HLT.Tmax= 25;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 14;
	    Corr_HLT.Tmax= 25;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 13;
	    Corr_HLT.Tmax= 19;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 22; //27; //22;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 13;
	    Corr_HLT.Tmax= 18;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 23; //27; //23;
	    Corr_HLT.Tmax= 26; //34; //26;
	  }
	  
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
	
      }
      else if(MESON=="B1s") {
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 17;
	     Corr_HLT.Tmax= 21;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin= 19;
	    Corr_HLT.Tmax= 26;
	   }
	 }
	 else if(ixg==1) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 17;
	      Corr_HLT.Tmax= 22;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else if(ixg==2) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 16;
	      Corr_HLT.Tmax= 21;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else if(ixg==3) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 16;
	      Corr_HLT.Tmax= 22;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 20;
	      Corr_HLT.Tmax= 26;
	    }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 14;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 19;
	     Corr_HLT.Tmax= 27;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 14;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 15;
	     Corr_HLT.Tmax= 24;
	   }
	 }
	 else if(ixg==2) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 14;
	      Corr_HLT.Tmax= 18;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin= 15;
	      Corr_HLT.Tmax= 21;
	    }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 14;
	     Corr_HLT.Tmax= 18;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin= 15;
	     Corr_HLT.Tmax= 21;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      Corr_HLT.Nt=2*(Corr_T.size()-1);
      Corr_HLT_boot.Nt= Corr_HLT.Nt;
      distr_t_list Corr_T_symm= Corr_T;
      distr_t_list Corr_T_symm_boot= Corr_T_boot;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_symm.distr_list.push_back( Corr_T.distr_list[Corr_HLT.Nt -t]);
	Corr_T_symm_boot.distr_list.push_back( Corr_T_boot.distr_list[Corr_HLT_boot.Nt -t]);
      }
      distr_t_list eff_M_Td_distr= Corr_HLT.effective_mass_t(Corr_T_symm, "" );
      distr_t_list eff_M_Td_boot_distr= Corr_HLT_boot.effective_mass_t( Corr_T_symm_boot, "");
      distr_t eff_M_Td= Corr_HLT.Fit_distr(eff_M_Td_distr)/a_distr;
      distr_t eff_M_Td_boot= Corr_HLT.Fit_distr(eff_M_Td_boot_distr);
      double Mphi_motion= sqrt( Mphi*Mphi + pow(kz/a_distr.ave(),2));
     
      Print_To_File({}, {eff_M_Td_distr.ave(), eff_M_Td_distr.err()}, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
     
      Print_To_File({}, { (Corr_T*FACT).ave(), (Corr_T*FACT).err() }, path_out+"/corr_2pts/"+TAG_CURR+"HLT_"+Ens_tags[iens]+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      distr_t M_phi_eff= SQRT_D( POW_D(eff_M_Td*a_distr,2) - pow(kz,2))/a_distr;
      distr_t Mphi_at_mb_half= SQRT_D(M_phi_eff*M_phi_eff + pow(0.5*5.367,2));

      cout<<"eff_M_T: "<<eff_M_Td.ave()<<" +- "<<eff_M_Td.err()<<" expected: "<<Mphi_motion<<endl;
      cout<<"M_phi: "<<M_phi_eff.ave()<<" "<<M_phi_eff.err()<<endl;
      distr_t_list phi_MT_distr= Corr_T_symm*EXPT_DL(eff_M_Td_distr);
      distr_t phi_MT = Corr_HLT.Fit_distr( Corr_T_symm*(EXPT_DL(eff_M_Td_distr)));
      distr_t phi_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_symm_boot*(EXPT_DL(eff_M_Td_boot_distr)));


      distr_t_list fake_distr_with_hole= Get_id_distr_list(Corr_T.size(), Njacks, UseJack);
      fake_distr_with_hole.distr_list[0] = 0.0*fake_distr_with_hole.distr_list[0];

      distr_t_list fake_boot_distr_with_hole= Get_id_distr_list(Corr_T.size(), 1000, 0);
      fake_boot_distr_with_hole.distr_list[0] = 0.0*fake_boot_distr_with_hole.distr_list[0];
      
      distr_t_list Corr_T_VMD = fake_distr_with_hole*phi_MT*EXPT_D(-1.0*eff_M_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_VMD_boot = fake_boot_distr_with_hole*phi_MT_boot*EXPT_D(-1.0*eff_M_Td_boot, Corr_T_boot.size());

     
      for(int t=1; t < Corr_T_VMD_boot.size(); t++) Corr_T_VMD_boot.distr_list[t] = Corr_T_VMD_boot.ave(t) + (Corr_T_VMD_boot.distr_list[t] - Corr_T_VMD_boot.ave(t))*Corr_T_VMD.err(t)/Corr_T_VMD_boot.err(t);
     
           
      distr_t_list Corr_T_sub = Corr_T - Corr_T_VMD;
      distr_t_list Corr_T_boot_sub = Corr_T_boot - Corr_T_VMD_boot;


      cout<<"PRINTING VMD PREDICTION FOR THE COUPLING ofr meson: "<<MESON+" xg: "<<xg_list.ave(ixg)<<" Ensemble: "<<Ens_tags[iens]<<endl;
      cout<<"g+ * f_V : "<<(FACT*phi_MT*2*eff_M_Td/(a_distr*M_phi_eff)).ave()<<" +- " <<(FACT*phi_MT*2*eff_M_Td/(a_distr*M_phi_eff)).err()<<endl;
      cout<<"Form factor at q^2=q'^2 =0: "<<(FACT*phi_MT*(eff_M_Td/Mphi_at_mb_half)*1.0/(  Mphi_at_mb_half*a_distr -  0.5*5.367*a_distr)).ave()<<endl;
      cout<<"##############################################################################"<<endl;


      //#########################################################################################################



      //Determine second-state 
      
      //ari-symmetrize
      distr_t_list Corr_T_ary_symm= Corr_T_sub;
      distr_t_list Corr_T_ary_symm_boot= Corr_T_boot_sub;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_ary_symm.distr_list.push_back( Corr_T_sub.distr_list[Corr_HLT.Nt -t]);
	Corr_T_ary_symm_boot.distr_list.push_back( Corr_T_boot_sub.distr_list[Corr_HLT.Nt -t]);
      }
      distr_t_list eff_M_prime_Td_distr= Corr_HLT.effective_mass_t( Corr_T_ary_symm, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_exc_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2));
      distr_t_list eff_M_prime_Td_boot_distr= Corr_HLT.effective_mass_t( Corr_T_ary_symm_boot, "");
      int Tmin_old= Corr_HLT.Tmin; int Tmax_old= Corr_HLT.Tmax;

      
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 5;
	    Corr_HLT.Tmax= 9;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=9;
	    Corr_HLT.Tmax=14;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 5;
	    Corr_HLT.Tmax= 9;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 4;
	    Corr_HLT.Tmax= 6;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }

	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 4;
	    Corr_HLT.Tmax= 6;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=10;
	    Corr_HLT.Tmax=14;
	  }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");

      }
      else if(MESON=="B1s") {
	 
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 6;
	     Corr_HLT.Tmax= 9;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=11;
	     Corr_HLT.Tmax=19;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 6;
	     Corr_HLT.Tmax= 9;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=11;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 8;
	     Corr_HLT.Tmax= 11;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 8;
	     Corr_HLT.Tmax= 11;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=17;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 8;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=10;
	     Corr_HLT.Tmax=16;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=7;
	     Corr_HLT.Tmax=11;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 5;
	     Corr_HLT.Tmax= 7;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      
      distr_t eff_M_prime_Td= Corr_HLT.Fit_distr(eff_M_prime_Td_distr)/a_distr;
      distr_t eff_M_prime_Td_boot= Corr_HLT_boot.Fit_distr(eff_M_prime_Td_boot_distr);
      distr_t_list phi_prime_MT_distr=  Corr_T_ary_symm*(EXPT_DL(eff_M_prime_Td_distr));
      distr_t_list phi_prime_MT_boot_distr=  Corr_T_ary_symm_boot*(EXPT_DL(eff_M_prime_Td_boot_distr));
      distr_t phi_prime_MT= Corr_HLT.Fit_distr( Corr_T_ary_symm*(EXPT_DL(eff_M_prime_Td_distr)));
      distr_t phi_prime_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_ary_symm_boot*(EXPT_DL(eff_M_prime_Td_boot_distr)));
      
      Corr_HLT.Tmin = Tmin_old; Corr_HLT.Tmax= Tmax_old;
      Corr_HLT_boot.Tmin = Tmin_old; Corr_HLT_boot.Tmax= Tmax_old;
      
      distr_t_list Corr_T_VMD_II = fake_distr_with_hole*phi_prime_MT*EXPT_D(-1.0*eff_M_prime_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_sub_II = Corr_T_sub - Corr_T_VMD_II;
      distr_t_list Corr_T_boot_VMD_II = fake_boot_distr_with_hole*phi_prime_MT_boot*EXPT_D(-1.0*eff_M_prime_Td_boot, Corr_T.size());
    
      for(int t=1; t < Corr_T_boot_VMD_II.size(); t++) Corr_T_boot_VMD_II.distr_list[t] = Corr_T_boot_VMD_II.ave(t) + (Corr_T_boot_VMD_II.distr_list[t] - Corr_T_boot_VMD_II.ave(t))*Corr_T_VMD_II.err(t)/Corr_T_boot_VMD_II.err(t);
    
      
      distr_t_list Corr_T_boot_sub_II = Corr_T_boot_sub - Corr_T_boot_VMD_II;



      //#########################################################################################################



      //Determine third-state 
      
      //ari-symmetrize
      distr_t_list Corr_T_aryary_symm= Corr_T_sub_II;
      distr_t_list Corr_T_aryary_symm_boot= Corr_T_boot_sub_II;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) {
	Corr_T_aryary_symm.distr_list.push_back( Corr_T_sub_II.distr_list[Corr_HLT.Nt -t]);
	Corr_T_aryary_symm_boot.distr_list.push_back( Corr_T_boot_sub_II.distr_list[Corr_HLT.Nt -t]);
      }
      distr_t_list eff_M_second_Td_distr= Corr_HLT.effective_mass_t( Corr_T_aryary_symm, path_out+"/mass/"+Ens_tags[iens]+"/"+TAG_CURR+"Vs_exc2_eff_xg_"+to_string_with_precision(xg_list.ave(ixg),2));
      distr_t_list eff_M_second_Td_boot_distr= Corr_HLT.effective_mass_t( Corr_T_aryary_symm_boot, "");
      Tmin_old= Corr_HLT.Tmin; Tmax_old= Corr_HLT.Tmax;
      if(MESON=="B0s") {
	if(ixg==0) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==1) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==2) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else if(ixg==3) {
	  if(Ens_tags[iens]=="cB211b.072.64") {
	    Corr_HLT.Tmin= 1;
	    Corr_HLT.Tmax= 3;
	  }
	  else if(Ens_tags[iens]=="cD211a.054.96") {
	    Corr_HLT.Tmin=2;
	    Corr_HLT.Tmax=4;
	  }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
	
      }
      else if(MESON=="B1s") {
	 
	 if(ixg==0) {
	    if(Ens_tags[iens]=="cB211b.072.64") {
	      Corr_HLT.Tmin= 1;
	      Corr_HLT.Tmax= 3;
	    }
	    else if(Ens_tags[iens]=="cD211a.054.96") {
	      Corr_HLT.Tmin=6;
	      Corr_HLT.Tmax=9;
	    }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=8;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=4;
	     Corr_HLT.Tmax=7;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else if(MESON=="B3s") {
	
	 if(ixg==0) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=2;
	     Corr_HLT.Tmax=6;
	   }
	 }
	 else if(ixg==1) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=4;
	   }
	 }
	 else if(ixg==2) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=3;
	   }
	 }
	 else if(ixg==3) {
	   if(Ens_tags[iens]=="cB211b.072.64") {
	     Corr_HLT.Tmin= 1;
	     Corr_HLT.Tmax= 3;
	   }
	   else if(Ens_tags[iens]=="cD211a.054.96") {
	     Corr_HLT.Tmin=1;
	     Corr_HLT.Tmax=3;
	   }
	 }
	 else crash("ixg: "+to_string(ixg)+" not yet implemented");
	 
      }
      else crash("Meson: "+MESON+" not yet implemented");

      Corr_HLT_boot.Tmin= Corr_HLT.Tmin;
      Corr_HLT_boot.Tmax= Corr_HLT.Tmax;

      
      
      distr_t eff_M_second_Td= Corr_HLT.Fit_distr(eff_M_second_Td_distr)/a_distr;
      distr_t eff_M_second_Td_boot= Corr_HLT_boot.Fit_distr(eff_M_second_Td_boot_distr);
      distr_t_list phi_second_MT_distr=  Corr_T_aryary_symm*(EXPT_DL(eff_M_second_Td_distr));
      distr_t_list phi_second_MT_boot_distr=  Corr_T_aryary_symm_boot*(EXPT_DL(eff_M_second_Td_boot_distr));
      distr_t phi_second_MT= Corr_HLT.Fit_distr( Corr_T_aryary_symm*(EXPT_DL(eff_M_second_Td_distr)));
      distr_t phi_second_MT_boot= Corr_HLT_boot.Fit_distr( Corr_T_aryary_symm_boot*(EXPT_DL(eff_M_second_Td_boot_distr)));
      
      Corr_HLT.Tmin = Tmin_old; Corr_HLT.Tmax= Tmax_old;
      Corr_HLT_boot.Tmin = Tmin_old; Corr_HLT_boot.Tmax= Tmax_old;

    
      
      distr_t_list Corr_T_VMD_III = fake_distr_with_hole*phi_second_MT*EXPT_D(-1.0*eff_M_second_Td*a_distr, Corr_T.size());
      distr_t_list Corr_T_sub_III = Corr_T_sub_II - Corr_T_VMD_III;
      distr_t_list Corr_T_boot_VMD_III = fake_boot_distr_with_hole*phi_second_MT_boot*EXPT_D(-1.0*eff_M_second_Td_boot, Corr_T.size());
  
      for(int t=1; t < Corr_T_boot_VMD_III.size(); t++) Corr_T_boot_VMD_III.distr_list[t] = Corr_T_boot_VMD_III.ave(t) + (Corr_T_boot_VMD_III.distr_list[t] - Corr_T_boot_VMD_III.ave(t))*Corr_T_VMD_III.err(t)/Corr_T_boot_VMD_III.err(t);
  
      distr_t_list Corr_T_boot_sub_III = Corr_T_boot_sub_II - Corr_T_boot_VMD_III;





      //#########################################################################################################

      
      Print_To_File({}, { (Corr_T_sub*FACT).ave(), (Corr_T_sub*FACT).err(), (Corr_T_sub_II*FACT).ave(), (Corr_T_sub_II*FACT).err(), (Corr_T_sub_III*FACT).ave(), (Corr_T_sub_II*FACT).err() }, path_out+"/corr_2pts/"+TAG_CURR+"sub_HLT_"+Ens_tags[iens]+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

      distr_t_list Corr_T_ori= Corr_T;

      if(Use_preconditioning) {
	Corr_T= Corr_T_sub;
	Corr_T_boot= Corr_T_boot_sub;
      }

      Print_To_File({}, { (FACT*phi_MT_distr/a_distr).ave(), (FACT*phi_MT_distr/a_distr).err(), (FACT*phi_prime_MT_distr/a_distr).ave(), (FACT*phi_prime_MT_distr/a_distr).err(), (FACT*phi_second_MT_distr/a_distr).ave(), (FACT*phi_second_MT_distr/a_distr).err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_MT_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
 
      //Generate Covariance matrix
      Vfloat cov_T, corr_T;
      Vfloat TT, RR;
      for(int tt=0;tt< tmax;tt++)
	for(int rr=0;rr< tmax;rr++) {
	  TT.push_back(tt);
	  RR.push_back(rr);
	  
	  corr_T.push_back( ( Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr])/( Corr_T_boot.err(tt)*Corr_T_boot.err(rr)));
	  cov_T.push_back( (Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr])*(Corr_T.err(tt)*Corr_T.err(rr)/(Corr_T_boot.err(tt)*Corr_T_boot.err(rr))));
	  
	  
	}

      int tmax_new=0;
      //tmax_new = tmax-1;
      for(int t=1;t<tmax;t++) if( Corr_T_ori.err(t)/fabs(Corr_T_ori.ave(t)) < 0.4) tmax_new++;
      int tmax_new_im=0;
      for(int t=1;t<tmax;t++) if( Corr_T_ori.err(t)/fabs(Corr_T_ori.ave(t)) < 0.4) tmax_new_im++;
   
   

      //print covariance matrix
      Print_To_File({},{TT,RR, cov_T, corr_T}, path_out+"/covariance/"+Ens_tags[iens]+"/"+preco_tag+"cov_T_xg_"+to_string_with_precision(xg.ave(),2)+".cov", "" , "");
  

      cout<<"Starting spectral reconstruction:..."<<endl;

      for(int iss=0; iss <(signed)sigmas_07_w0.size(); iss++) {	


	F_T_d_MB_RE_VMD_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()));
	F_T_d_MB_IM_VMD_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave()));

	F_T_d_MB_RE_VMD_II_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   );
	F_T_d_MB_IM_VMD_II_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())    );

	F_T_d_MB_RE_VMD_III_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())   + FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()) +  FACT*phi_second_MT*K_RE_distr_sub( eff_M_second_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())  );
	F_T_d_MB_IM_VMD_III_state_sm_list[ixg][iss].distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off,  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave()) +  FACT*phi_second_MT*K_IM_distr( eff_M_second_Td*a_distr, Eg_MB_off, sigmas_07_w0[iss]*a_distr.ave())    );

	distr_t_list VMD_VIRT_SCAN_RE(UseJack), VMD_VIRT_SCAN_IM(UseJack);
	distr_t_list VMD_II_VIRT_SCAN_RE(UseJack), VMD_II_VIRT_SCAN_IM(UseJack);
	distr_t_list VMD_III_VIRT_SCAN_RE(UseJack), VMD_III_VIRT_SCAN_IM(UseJack);

	for(int vir=0;vir<(signed)virtualities.size(); vir++) {

	  VMD_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave()));
	  VMD_VIRT_SCAN_IM.distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave()) );

	  VMD_II_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   );
	  VMD_II_VIRT_SCAN_IM.distr_list.push_back(  FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr , M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())    );

	  VMD_III_VIRT_SCAN_RE.distr_list.push_back(  FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr , M_P*virtualities[vir] , sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   +  FACT*phi_second_MT*K_RE_distr_sub( eff_M_second_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())   );
	  VMD_III_VIRT_SCAN_IM.distr_list.push_back( FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir],  sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())  + FACT*phi_second_MT*K_IM_distr( eff_M_second_Td*a_distr, M_P*virtualities[vir], sigmas_07_w0[iss]*a_distr.ave())    );

	}

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_VIRT_SCAN_RE.ave(), VMD_VIRT_SCAN_RE.err(), VMD_VIRT_SCAN_IM.ave(), VMD_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_II_VIRT_SCAN_RE.ave(), VMD_II_VIRT_SCAN_RE.err(), VMD_II_VIRT_SCAN_IM.ave(), VMD_II_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/II_xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");

	Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(), Njacks)*M_P/a_distr).ave(), VMD_III_VIRT_SCAN_RE.ave(), VMD_III_VIRT_SCAN_RE.err(), VMD_III_VIRT_SCAN_IM.ave(), VMD_III_VIRT_SCAN_IM.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VMD_virt_scan/III_xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_sm_"+to_string_with_precision(sigmas_07_w0[iss],3)+"_"+Ens_tags[iens], "", "");
      }

      

      if(!Skip_spectral_reconstruction_07) {

	vector<distr_t_list> F_T_d_RE_virt_scan;
	vector<distr_t_list> F_T_d_IM_virt_scan;

	for(int isg=0;isg<(signed)sigmas_07.size();isg++) {
	  F_T_d_RE_virt_scan.emplace_back(UseJack, virtualities.size(), Njacks);
	  F_T_d_IM_virt_scan.emplace_back(UseJack, virtualities.size(), Njacks);
	}

	#pragma omp parallel for schedule(dynamic)
	//spectral reconstruction for second time ordering
	for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

	  double s= sigmas_07[isg]*a_distr.ave();

	  double syst_T;
	  double mult_T_IM= 0.5;

	  if(SM_TYPE=="FF_Gauss") mult_T_IM=0.3;
	  double mult_T_RE=0.05;
	  double Ag_target=1e-2;
	  if(sigmas_07[isg] < 0.5) Ag_target=1e-2;
	  else if(sigmas_07[isg] < 1.5) Ag_target=1e-2;
	  double th= E0_fact*Mphi_motion;
	  double l_re_T;
	  mult_T_RE=1;
	  
	  
	
	  if(SM_TYPE=="FF_Exp") {
	  //select multiplicities
	    if(MESON=="B0s" || MESON=="B1s") mult_T_IM=0.1;
	    if(MESON=="B3s") mult_T_IM=0.3;
	    //if(MESON=="B3s" && ixg == 0) mult_T_RE=0.005;
	    //if(MESON=="B3s" && ixg == 1) mult_T_RE=0.01;
	    if(MESON=="B3s" && ixg >= 2 && Ens_tags[iens] == "cB211b.072.64") mult_T_RE=0.05;
	    if(MESON=="B0s" && ixg==3 && Ens_tags[iens] == "cD211a.054.96") mult_T_RE=10;


	    if( isg < sigmas_07.size()/2 && MESON == "B1s"  ) mult_T_IM=2.5;
	    if( MESON == "B3s") mult_T_IM=4;

	   
	    

	  }
	  
	 
    

	  distr_t RE_sm= Get_id_jack_distr(Njacks);
	  distr_t IM_sm= Get_id_jack_distr(Njacks);
	  

	  cout<<"#########################  E= MB - Egamma ######################### "<<endl<<flush;

	  //modify s in order to make it scale with the energy oE

	  s= s*Eg_MB_off.ave()/(MBs*a_distr.ave());
	  if(iens==0) sigma_simulated[ixg][isg] = sigmas_07[isg]*Eg_MB_off_new.ave()/MBs;

	  distr_t RE_MB_sm = Get_Laplace_transfo_tmin( Eg_MB_off.ave(),  s, th*a_distr.ave(),  Nts[iens], 3, tmax_new, prec_07, SM_TYPE+"_RE_ixg_"+to_string(ixg+1),K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1);

	  distr_t REMINDER_MB(UseJack,Njacks);

	  for(int t=1;t<Corr_T.size();t++) REMINDER_MB = REMINDER_MB + FACT*Corr_T[t]; 
	  
	  
	  RE_MB_sm = RE_MB_sm.ave() +  (RE_MB_sm - RE_MB_sm.ave())*(sqrt( pow(syst_T,2) + pow(RE_MB_sm.err(),2)))/RE_MB_sm.err() -REMINDER_MB;


	  distr_t IM_MB_sm = Get_Laplace_transfo_tmin( Eg_MB_off.ave(),  s, th*a_distr.ave(),  Nts[iens], 3, tmax_new_im, prec_07, SM_TYPE+"_IM_ixg_"+to_string(ixg+1),K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"MB_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , preco_tag+MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1);
	  
	  IM_MB_sm = IM_MB_sm.ave() +  (IM_MB_sm - IM_MB_sm.ave())*(sqrt( pow(syst_T,2) + pow(IM_MB_sm.err(),2)))/IM_MB_sm.err();
	  
	  F_T_d_MB_RE_sm_list[ixg][isg].distr_list.push_back( RE_MB_sm + ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, Eg_MB_off, s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
	  
	  F_T_d_MB_IM_sm_list[ixg][isg].distr_list.push_back( IM_MB_sm +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, Eg_MB_off, s)+ 0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, Eg_MB_off, s) ));
	  

	  if(virtuality_scan) {
	    
	    for(int vir=0;vir<(signed)virtualities.size(); vir++) {

	      F_T_d_RE_virt_scan[isg].distr_list[vir] = Get_Laplace_transfo ( M_P.ave()*virtualities[vir],  s, th*a_distr.ave(),  Nts[iens], tmax_new, prec_07, SM_TYPE+"_RE",K_RE, Corr_T, syst_T, mult_T_RE ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"VIRT_SCAN_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1) +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_RE_distr_sub( eff_M_Td*a_distr, M_P*virtualities[vir], s) + 0.0*FACT*phi_prime_MT*K_RE_distr_sub( eff_M_prime_Td*a_distr, M_P*virtualities[vir], s) - REMINDER_MB);
	      
	      F_T_d_IM_virt_scan[isg].distr_list[vir] = Get_Laplace_transfo ( M_P.ave()*virtualities[vir],  s, th*a_distr.ave(),  Nts[iens], tmax_new_im, prec_07, SM_TYPE+"_IM",K_IM, Corr_T, syst_T, mult_T_IM ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"VIRT_SCAN_T_"+Ens_tags[iens], Ag_target,0, FACT, 0.0 , MESON+"_07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.0,1)  +  ((Use_preconditioning==false)?0.0:1.0)*(FACT*phi_MT*K_IM_distr( eff_M_Td*a_distr, M_P*virtualities[vir], s)+0.0*FACT*phi_prime_MT*K_IM_distr( eff_M_prime_Td*a_distr, M_P*virtualities[vir], s));

	    }

	  }
	 	  
	  
	}
	//if virtuality scan, print results
	if(virtuality_scan) {
	  for(int isg=0;isg<(signed)sigmas_07.size();isg++) {
	    Print_To_File({}, {(virtualities*Get_id_jack_distr_list(virtualities.size(),Njacks)*M_P/a_distr).ave(), F_T_d_RE_virt_scan[isg].ave(), F_T_d_RE_virt_scan[isg].err(), F_T_d_IM_virt_scan[isg].ave(), F_T_d_IM_virt_scan[isg].err() }, path_out+"/FF_d_II/"+Ens_tags[iens]+"/VIRT_SCAN_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(0.1+0.1*ixg,2), "", "");
	  }
	}
      }  
      
    }

  }

  

  //print results
  //per kinematic
  for(int ixg=0;ixg<n_xg;ixg++) {

    Print_To_File( Ens_tags, { F_T_u_list[ixg].ave(), F_T_u_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { F_T_d_I_list[ixg].ave(), F_T_d_I_list[ixg].err() , F_T_d_sp_I_list[ixg].ave(), F_T_d_sp_I_list[ixg].err()    } , path_out+"/FF_d_I/"+TAG_CURR+"F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_u_real_list[ixg].ave(), FV_T_u_real_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_d_real_list[ixg].ave(), FV_T_d_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_u_real_list[ixg].ave(), FA_T_u_real_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_d_real_list[ixg].ave(), FA_T_d_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    Print_To_File( Ens_tags, { FV_T_d_sp_real_list[ixg].ave(), FV_T_d_sp_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FV_T_sp_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_d_sp_real_list[ixg].ave(), FA_T_d_sp_real_list[ixg].err()     } , path_out+"/FF_d/"+TAG_CURR+"FA_T_sp_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    Print_To_File( Ens_tags, { F_T_u_VMD_list[ixg].ave(), F_T_u_VMD_list[ixg].err(), F_T_u_VMD_spectre_list[ixg].ave(), F_T_u_VMD_spectre_list[ixg].err()     } , path_out+"/FF_u/"+TAG_CURR+"VMD_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");

    
    if(!Skip_spectral_reconstruction_07) {
    for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

   
      Print_To_File( Ens_tags, { F_T_d_MB_RE_sm_list[ixg][isg].ave(), F_T_d_MB_RE_sm_list[ixg][isg].err(), F_T_d_MB_IM_sm_list[ixg][isg].ave(), F_T_d_MB_IM_sm_list[ixg][isg].err()     } , path_out+"/FF_d_II/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");


    }
    }
  }
  //per ensemble
  for(int iens=0; iens<Nens;iens++) {
    distr_t_list F_T_u_per_ens(UseJack),  F_T_d_I_per_ens(UseJack);
    distr_t_list FV_T_u_real_per_ens(UseJack),  FV_T_d_real_per_ens(UseJack);
    distr_t_list FA_T_u_real_per_ens(UseJack),  FA_T_d_real_per_ens(UseJack);
    distr_t_list FA_T_d_I_real_per_ens(UseJack), FA_T_d_II_real_per_ens(UseJack);
    distr_t_list FV_T_d_I_real_per_ens(UseJack), FV_T_d_II_real_per_ens(UseJack);
    distr_t_list FA_T_d_sp_real_per_ens(UseJack), FV_T_d_sp_real_per_ens(UseJack);
    distr_t_list F_T_u_VMD_per_ens(UseJack), F_T_u_VMD_spectre_per_ens(UseJack);
    distr_t_list F_T_d_sp_I_per_ens(UseJack);

    for(int ixg=0;ixg<n_xg;ixg++) {
      F_T_u_per_ens.distr_list.push_back( F_T_u_list[ixg].distr_list[iens]);
      F_T_d_I_per_ens.distr_list.push_back( F_T_d_I_list[ixg].distr_list[iens]);
      F_T_d_sp_I_per_ens.distr_list.push_back( F_T_d_sp_I_list[ixg].distr_list[iens]);
      FV_T_u_real_per_ens.distr_list.push_back( FV_T_u_real_list[ixg].distr_list[iens]);
      FV_T_d_real_per_ens.distr_list.push_back( FV_T_d_real_list[ixg].distr_list[iens]);
      FA_T_u_real_per_ens.distr_list.push_back( FA_T_u_real_list[ixg].distr_list[iens]);
      FA_T_d_real_per_ens.distr_list.push_back( FA_T_d_real_list[ixg].distr_list[iens]);

      FA_T_d_I_real_per_ens.distr_list.push_back( FA_T_d_I_real_list[ixg].distr_list[iens]);
      FA_T_d_II_real_per_ens.distr_list.push_back( FA_T_d_II_real_list[ixg].distr_list[iens]);

      FV_T_d_I_real_per_ens.distr_list.push_back( FV_T_d_I_real_list[ixg].distr_list[iens]);
      FV_T_d_II_real_per_ens.distr_list.push_back( FV_T_d_II_real_list[ixg].distr_list[iens]);

      FV_T_d_sp_real_per_ens.distr_list.push_back( FV_T_d_sp_real_list[ixg].distr_list[iens]);
      FA_T_d_sp_real_per_ens.distr_list.push_back( FA_T_d_sp_real_list[ixg].distr_list[iens]);

      F_T_u_VMD_per_ens.distr_list.push_back( F_T_u_VMD_list[ixg].distr_list[iens]);
      F_T_u_VMD_spectre_per_ens.distr_list.push_back( F_T_u_VMD_spectre_list[ixg].distr_list[iens]);

      
      distr_t_list F_T_d_MB_RE_sm_per_ens_per_kin, F_T_d_MB_IM_sm_per_ens_per_kin;
      
      
      distr_t_list F_T_d_MB_RE_VMD_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_sm_per_ens_per_kin;
      
      
      distr_t_list F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin;
      
      
      distr_t_list F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin, F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin;
       
      
      if(!Skip_spectral_reconstruction_07) {
	for(int iss=0; iss<(signed)sigmas_07_w0.size(); iss++) {
	
	  F_T_d_MB_RE_VMD_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_sm_list[ixg][iss].distr_list[iens]);

	  F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_II_state_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_II_state_sm_list[ixg][iss].distr_list[iens]);

	 
	  F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_VMD_III_state_sm_list[ixg][iss].distr_list[iens]);
	  F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_VMD_III_state_sm_list[ixg][iss].distr_list[iens]);
	}
	for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {
	

	  F_T_d_MB_RE_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens]);
	  F_T_d_MB_IM_sm_per_ens_per_kin.distr_list.push_back( F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens]);
	}

	
      
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], F_T_d_MB_RE_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_sm_per_ens_per_kin.err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_II_state_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_II_state_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_MB_II_state_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07_w0, F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.ave(), F_T_d_MB_RE_VMD_III_state_sm_per_ens_per_kin.err(), F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_VMD_III_state_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_MB_III_state_F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_u_per_ens[ixg]+ F_T_d_I_per_ens[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_u_per_ens[ixg]+F_T_d_I_per_ens[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	Print_To_File( {}, {sigmas_07, sigma_simulated[ixg], (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_d_I_per_ens[ixg]).ave(), (F_T_d_MB_RE_sm_per_ens_per_kin+F_T_d_I_per_ens[ixg]).err(), F_T_d_MB_IM_sm_per_ens_per_kin.ave(), F_T_d_MB_IM_sm_per_ens_per_kin.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
	
		
      }
    }

    Print_To_File( {}, {xg_list.ave(),  F_T_u_per_ens.ave(), F_T_u_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"F_T", "", "");
    Print_To_File( {}, {xg_list.ave(), F_T_d_I_per_ens.ave(), F_T_d_I_per_ens.err(), F_T_d_sp_I_per_ens.ave(), F_T_d_sp_I_per_ens.err()     } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"F_T", "", "");
    Print_To_File( {}, {xg_list.ave(),  FV_T_u_real_per_ens.ave(), FV_T_u_real_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_real_per_ens.ave(), FV_T_d_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(),  FA_T_u_real_per_ens.ave(), FA_T_u_real_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FA_T_d_real_per_ens.ave(), FA_T_d_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(), FA_T_d_I_real_per_ens.ave(), FA_T_d_I_real_per_ens.err()     } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FA_T_d_II_real_per_ens.ave(), FA_T_d_II_real_per_ens.err()     } , path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(), FV_T_d_I_real_per_ens.ave(), FV_T_d_I_real_per_ens.err()     } , path_out+"/FF_d_I/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_II_real_per_ens.ave(), FV_T_d_II_real_per_ens.err()     } , path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_real", "", "");

    Print_To_File( {}, {xg_list.ave(),  FA_T_d_sp_real_per_ens.ave(), FA_T_d_sp_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FA_T_sp_real", "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_sp_real_per_ens.ave(), FV_T_d_sp_real_per_ens.err()     } , path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+"FV_T_sp_real", "", "");
    
    Print_To_File( {}, {xg_list.ave(),  F_T_u_VMD_per_ens.ave(), F_T_u_VMD_per_ens.err(), F_T_u_VMD_spectre_per_ens.ave(), F_T_u_VMD_spectre_per_ens.err()     } , path_out+"/FF_u/"+Ens_tags[iens]+"/"+TAG_CURR+"VMD_F_T", "", "");
    
    
    if(!Skip_spectral_reconstruction_07) {
    for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {

    
      distr_t_list F_T_d_MB_RE_sm_per_ens_per_sigma, F_T_d_MB_IM_sm_per_ens_per_sigma;
      
       for(int ixg=0;ixg<n_xg;ixg++) {

	 F_T_d_MB_RE_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_MB_RE_sm_list[ixg][isg].distr_list[iens]);
	 F_T_d_MB_IM_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_MB_IM_sm_list[ixg][isg].distr_list[iens]);
       }


       Print_To_File( {}, {xg_list.ave(), F_T_d_MB_RE_sm_per_ens_per_sigma.ave(), F_T_d_MB_RE_sm_per_ens_per_sigma.err(), F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF_d_II/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       Print_To_File( {}, {xg_list.ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).err(),  F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       Print_To_File( {}, {xg_list.ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_d_I_per_ens).ave(), (F_T_d_MB_RE_sm_per_ens_per_sigma + F_T_d_I_per_ens).err(),  F_T_d_MB_IM_sm_per_ens_per_sigma.ave(), F_T_d_MB_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF_d/"+Ens_tags[iens]+"/"+TAG_CURR+preco_tag+"F_MB_T_sm_"+to_string_with_precision(sigmas_07[isg],3), "", "");
       
    }
    }
  }


  //extrapolate the results for each sigma and xg to the continuum limit, employing either a constant or linear fit in a^2 (to be combined with BMA with uniform weight 1/2 )

  D(1);

  for(int ixg=0; ixg<n_xg;ixg++) {

    //################################################################################
    // b-quark contribution
    distr_t b_B64(UseJack), b_D96(UseJack) ;
    for(int iens=0;iens<Ens_tags.size();iens++) {
      if(Ens_tags[iens] == "cB211b.072.64") {
	b_B64= F_T_u_list[ixg][iens];
      }
      else b_D96= F_T_u_list[ixg][iens];
    }

         
    distr_t FT_b_const= (1/pow(b_B64.err(),2))*b_B64 + (1/pow(b_D96.err(),2))*b_D96;
    FT_b_const = FT_b_const/( (1/pow(b_B64.err(),2)) +  (1/pow(b_D96.err(),2)));
    double ch2_const = (FT_b_const.ave() - b_B64.ave())/pow(b_B64.err(),2) +  (FT_b_const.ave() - b_D96.ave())/pow(b_D96.err(),2) ;
       
    distr_t FT_b_a2= b_D96*a_B*a_B - b_B64*a_D*a_D;
    FT_b_a2 = FT_b_a2/( a_B*a_B - a_D*a_D);


        
    //combine using Eq. 29
    distr_t FT_b(UseJack);
    if(ch2_const < 2) {
      FT_b  = BMA_Eq_29( {FT_b_const, FT_b_a2});
    }
    else FT_b= FT_b_a2;
    //################################################################################

    D(2);
    

   
    // s-quark contribution

    distr_t_list FT_s_RE_sigma(UseJack), FT_s_IM_sigma(UseJack);
    

    for(int is=0;is<(signed)sigma_simulated[ixg].size();is++) {

     
      
	distr_t s_I_B64(UseJack), s_I_D96(UseJack) ;
	distr_t s_II_B64(UseJack), s_II_D96(UseJack);
	distr_t s_B64_IM(UseJack), s_D96_IM(UseJack);
	for(int iens=0;iens<Ens_tags.size();iens++) {
	  if(Ens_tags[iens] == "cB211b.072.64") {
	    s_I_B64= F_T_d_I_list[ixg][iens];

            s_II_B64= F_T_d_MB_RE_sm_list[ixg][is][iens];
	    s_B64_IM= F_T_d_MB_IM_sm_list[ixg][is][iens];
	  }
	  else {
	    s_I_D96= F_T_d_I_list[ixg][iens];
	    s_II_D96= F_T_d_MB_RE_sm_list[ixg][is][iens];
	    s_D96_IM= F_T_d_MB_IM_sm_list[ixg][is][iens];
	  }
	}

	//add sigma contribution
	distr_t s_B64_RE = s_I_B64 + s_II_B64;
	distr_t s_D96_RE = s_I_D96 + s_II_D96;

	//extrapolate


	distr_t FT_s_RE_const =  (1/pow(s_B64_RE.err(),2))*s_B64_RE + (1/pow(s_D96_RE.err(),2))*s_D96_RE;
	distr_t FT_s_RE_a2= s_D96_RE*a_B*a_B - s_B64_RE*a_D*a_D;
	FT_s_RE_a2 = FT_s_RE_a2/( a_B*a_B - a_D*a_D);
	//combine using Eq. 29
	distr_t FT_s_RE(UseJack);
	if(ch2_const < 2) {
	  FT_s_RE  = BMA_Eq_29( {FT_s_RE_const, FT_s_RE_a2});
	}
	else FT_s_RE= FT_s_RE_a2;


	distr_t FT_s_IM_const =  (1/pow(s_B64_IM.err(),2))*s_B64_IM + (1/pow(s_D96_IM.err(),2))*s_D96_IM;
	distr_t FT_s_IM_a2= s_D96_IM*a_B*a_B - s_B64_IM*a_D*a_D;
	FT_s_IM_a2 = FT_s_IM_a2/( a_B*a_B - a_D*a_D);
	//combine using Eq. 29
	distr_t FT_s_IM(UseJack);
	if(ch2_const < 2) {
	  FT_s_IM  = BMA_Eq_29( {FT_s_IM_const, FT_s_IM_a2});
	}
	else FT_s_IM= FT_s_IM_a2;

	FT_s_RE_sigma.distr_list.push_back(FT_s_RE);
	FT_s_IM_sigma.distr_list.push_back(FT_s_IM);

    }

  


    //extrapolate to vanishing sigma
    class ipar_07 {
    public:
      ipar_07() : FF_RE(0.0), FF_RE_err(0.0), FF_IM(0.0), FF_IM_err(0.0) {}
      double FF_RE, FF_RE_err, FF_IM, FF_IM_err, sigma;
    };
  
    class fpar_07 {
    public:
      fpar_07() {}
      fpar_07(const Vfloat &par) {
	if((signed)par.size() != 3) crash("In class fpar_07  class constructor Vfloat par has size != 3");
	R=par[0];
	A=par[1];
	B=par[2];
	
      }
      double R,A,B;
    };
    
    int sigma_to_fit=4;
    //init bootstrap fit
    bootstrap_fit<fpar_07,ipar_07> bf_07(Njacks);
    bf_07.set_warmup_lev(1); //sets warmup
    bf_07.Set_number_of_measurements(sigma_to_fit);
    bf_07.Set_verbosity(1);

    bf_07.Add_par("R", -0.05, 0.1);
    bf_07.Add_par("A", 1.0, 0.1);
    bf_07.Add_par("B", 1.0 , 0.1);
   
  
    //fit on mean values to get ch2
    bootstrap_fit<fpar_07,ipar_07> bf_07_ch2(1);
    bf_07_ch2.set_warmup_lev(1); //sets warmup
    bf_07_ch2.Set_number_of_measurements(sigma_to_fit);
    bf_07_ch2.Set_verbosity(1);
    bf_07_ch2.Add_par("R", -0.05, 0.1);
    bf_07_ch2.Add_par("A", 1.0, 0.1);
    bf_07_ch2.Add_par("B", 1.0, 0.1);


    bool Is_RE=true;



    //ansatz
    bf_07.ansatz=  [](const fpar_07 &p, const ipar_07 &ip) {

      return p.R + p.A*(ip.sigma) + p.B*pow(ip.sigma,2);
  
    	 	  
    };

  
    bf_07.measurement=  [&Is_RE ](const fpar_07 &p, const ipar_07 &ip) {

      if(Is_RE) return ip.FF_RE;

      return ip.FF_IM;
      
    };
    
    bf_07.error=  [&Is_RE ](const fpar_07 &p, const ipar_07 &ip) {

      if(Is_RE) return ip.FF_RE_err;

      return ip.FF_IM_err;
    
  };
	
  bf_07_ch2.ansatz= bf_07.ansatz;
  bf_07_ch2.measurement = bf_07.measurement;
  bf_07_ch2.error = bf_07.error;

  D(2);

  //start fitting
  //fill the data
  vector<vector<ipar_07>> data_07(Njacks);
  vector<vector<ipar_07>> data_07_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_07> Bt_fit_07_RE;
  boot_fit_data<fpar_07> Bt_fit_07_RE_ch2;
  boot_fit_data<fpar_07> Bt_fit_07_IM;
  boot_fit_data<fpar_07> Bt_fit_07_IM_ch2;
  for(auto &data_iboot: data_07) data_iboot.resize(sigma_to_fit);
  for(auto &data_iboot: data_07_ch2) data_iboot.resize(sigma_to_fit);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int is=0;is<sigma_to_fit;is++) {
      data_07[ijack][is].FF_RE = (FT_s_RE_sigma[is]).distr[ijack];
      data_07[ijack][is].FF_RE_err= (FT_s_RE_sigma[is]).err();
      data_07[ijack][is].FF_IM = (FT_s_IM_sigma[is]).distr[ijack];
      data_07[ijack][is].FF_IM_err= (FT_s_IM_sigma[is]).err();
      data_07[ijack][is].sigma = sigma_simulated[ixg][is];
      if(ijack==0) {
	data_07_ch2[ijack][is].FF_RE = (FT_s_RE_sigma[is]).ave();
	data_07_ch2[ijack][is].FF_RE_err= (FT_s_RE_sigma[is]).err();
	data_07_ch2[ijack][is].FF_IM = (FT_s_IM_sigma[is]).ave();
	data_07_ch2[ijack][is].FF_IM_err= (FT_s_IM_sigma[is]).err();
	data_07_ch2[ijack][is].sigma = sigma_simulated[ixg][is];
      }
    }
  }
  
  //append
  bf_07.Append_to_input_par(data_07);
  bf_07_ch2.Append_to_input_par(data_07_ch2);
  //fit
  cout<<"Fitting real part of the form factor for ixg: "<<ixg<<endl;

  D(3);

  Bt_fit_07_RE= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_RE_ch2= bf_07_ch2.Perform_bootstrap_fit();    
  double ch2_red_RE= Bt_fit_07_RE_ch2.get_ch2_ave()/( sigma_to_fit - bf_07.Get_number_of_fit_pars());
  cout<<"Reduced ch2: "<<ch2_red_RE<<endl;
  Is_RE=false;
  cout<<"Fitting imag part of the form factor for ixg: "<<ixg<<endl;
  Bt_fit_07_IM= bf_07.Perform_bootstrap_fit();
  Bt_fit_07_IM_ch2= bf_07_ch2.Perform_bootstrap_fit();    
  double ch2_red_IM= Bt_fit_07_IM_ch2.get_ch2_ave()/( sigma_to_fit - bf_07.Get_number_of_fit_pars());
  cout<<"Reduced ch2: "<<ch2_red_IM<<endl;

  D(4);

  //retrieve parameters

  distr_t R_RE(UseJack), R_IM(UseJack), A_RE(UseJack), A_IM(UseJack), B_RE(UseJack), B_IM(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++) {
    R_RE.distr.push_back( Bt_fit_07_RE.par[ijack].R);
    A_RE.distr.push_back( Bt_fit_07_RE.par[ijack].A);
    B_RE.distr.push_back( Bt_fit_07_RE.par[ijack].B);

    R_IM.distr.push_back( Bt_fit_07_IM.par[ijack].R);
    A_IM.distr.push_back( Bt_fit_07_IM.par[ijack].A);
    B_IM.distr.push_back( Bt_fit_07_IM.par[ijack].B);
    
  }

  //get fit function

  Vfloat sigmas_to_print;
  int Nsigmas_to_print=200;
  for(int n=0;n<Nsigmas_to_print;n++) sigmas_to_print.push_back( 3.0*n/(Nsigmas_to_print-1));

  distr_t_list F_print_RE(UseJack,Nsigmas_to_print,Njacks), F_print_IM(UseJack,Nsigmas_to_print, Njacks);
  
  for(int n=0;n<Nsigmas_to_print;n++) {

    ipar_07 fS;
    fS.sigma= sigmas_to_print[n];
    for(int ijack=0;ijack<Njacks;ijack++) {
      F_print_RE[n].distr[ijack] = bf_07.ansatz( Bt_fit_07_RE.par[ijack], fS);
      F_print_IM[n].distr[ijack] = bf_07.ansatz( Bt_fit_07_IM.par[ijack], fS);

    }
    
  }


  //print data and fit parameters
  boost::filesystem::create_directory(path_out+"/FF_d_extr");
  string out_path_data = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+".data" ;
  string out_path_fit_func = path_out+"/FF_d_extr/FT_ixg_"+to_string(ixg)+".fit_func";

  Print_To_File({}, {sigma_simulated[ixg], FT_s_RE_sigma.ave(), FT_s_RE_sigma.err(), FT_s_IM_sigma.ave(), FT_s_IM_sigma.err()}, out_path_data, "", "");
  Print_To_File({}, {sigmas_to_print, F_print_RE.ave(), F_print_RE.err(), F_print_IM.ave(), F_print_IM.err()}, out_path_fit_func, "", "");


  
  
      
      
	

  return_class.F_T_u.distr_list.push_back( FT_b);
  return_class.F_T_d_RE.distr_list.push_back( R_RE);
  return_class.F_T_d_IM.distr_list.push_back( R_IM);
    
    

  

    
  }

  
  return_class.num_xg= n_xg;
  return_class.MESON= MESON;
  

    
 




  return return_class;

}
