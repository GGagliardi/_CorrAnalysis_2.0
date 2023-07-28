#include "../include/virtual_07.h"
#include "Corr_analysis.h"
using namespace std;

bool verbose_lev_07=1;
Vfloat sigmas_07({0.6,0.5,0.4,0.35,0.3,0.25,0.2,0.15, 0.1}); //sigma in GeV  
int prec_07=128;
const string MODE_FF="TANT";
const bool Skip_spectral_reconstruction=true;
const bool Reconstruct_axial_part=true;
const bool Reconstruct_vector_part=true;
const bool Perform_FF_and_Br_reconstruction=false;
const bool Perform_eps_extrapolation=true;
const double Mjpsi= 3.0969; //GeV
const double Mphi= 1.019461; //GeV
const double MDs_phys = 1.96847; // GeV
const double E0_fact = 0.9;
const bool CONS_EM_CURRENT = true;
const string SM_TYPE= "FF_Exp";



rt_07_Bs Get_virtual_tensor_FF(int n_xg, bool UseJack, int Njacks, string Meson,  string Corr_path, string path_out) {


  rt_07_Bs return_class;

  PrecFloat::setDefaultPrecision(prec_07);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;

  int t_07=18;

  string TAG_CURR="";
  if(CONS_EM_CURRENT==false) TAG_CURR="LOC_";

 
  int size_mu_nu= 4;

  //BK
  vector<vector<vector<data_t>>> C_B_u_data(size_mu_nu), C_B_d_data(size_mu_nu);
  //T
  vector<vector<vector<data_t>>> C_T_u_data(size_mu_nu), C_T_d_data(size_mu_nu);


  
  for(int mu=0;mu<size_mu_nu;mu++) {

    C_B_u_data[mu].resize(size_mu_nu);
    C_B_d_data[mu].resize(size_mu_nu);
    C_T_u_data[mu].resize(size_mu_nu);
    C_T_d_data[mu].resize(size_mu_nu);

   
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_B_u_data[mu][nu].resize(n_xg);
      C_B_d_data[mu][nu].resize(n_xg);
      C_T_u_data[mu][nu].resize(n_xg);
      C_T_d_data[mu][nu].resize(n_xg);

    
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

 
  //loop over mu and nu axial
  vector<pair<int,int>> mu_nu_pair_A({make_pair(1,1),make_pair(2,2)});
  vector<pair<int,int>> mu_nu_pair_V({make_pair(1,2),make_pair(2,1)});
  vector<pair<int,int>> red_mu_nu_pair_V({make_pair(1,2)});
  vector<pair<int,int>> red_mu_nu_pair_A({make_pair(1,1) });
  
  
  mu_nu_pair_A = { make_pair(1,1), make_pair(2,2) };
  
  for(int ixg=0;ixg<n_xg;ixg++) {


   
    //B
    for(auto &pair_A : mu_nu_pair_A) {
      int mu=pair_A.first;
      int nu=pair_A.second;
  
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //B
      //u
      C_B_u_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_u_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_B_d_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }

    //T
    for(auto &pair_V : mu_nu_pair_V) {
      int mu=pair_V.first;
      int nu=pair_V.second;
    
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //T
      //u
      C_T_u_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_u_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_T_d_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }
  }



   
  GaussianMersenne GM(652205123);

  //resample RCs
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

        
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

  }

  int Nens= C_B_u_data[1][1][0].size;
  vector<string> Ens_tags= C_B_u_data[1][1][0].Tag;
  Vint Nts=C_B_u_data[1][1][0].nrows;


  vector<distr_t_list> FV_T_u, FA_T_u;
  vector<distr_t_list> FV_T_d_I, FA_T_d_I;
  vector<vector<distr_t_list>> FV_T_d_sm(n_xg), FA_T_d_sm(n_xg);

  for(int ixg=0; ixg<n_xg;ixg++) {

    FV_T_u.emplace_back(UseJack);
    FA_T_u.emplace_back(UseJack);
    FV_T_d_I.emplace_back(UseJack);
    FA_T_d_I.emplace_back(UseJack);
    for(int is=0; is<(signed)sigmas_07.size(); is++) {
      FV_T_d_sm[ixg].emplace_back(UseJack);
      FA_T_d_sm[ixg].emplace_back(UseJack);
    }
  }


    auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {


    if(SM_TYPE=="FF_Gauss") {
      PrecFloat x = (E-m)/(sqrt(PrecFloat(2))*s);
      return sqrt(2)*DawsonF(x)/s;
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

     
      return ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
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
    if(SM_TYPE=="FF_Gauss_Schwartz") return precPi()*Get_exact_gauss(E, m, s, E0);
    if(SM_TYPE=="FF_Cauchy") {
      PrecFloat t= (E-m);
      return s/( t*t + s*s);
    }

    if(SM_TYPE=="FF_Exp") {

      PrecFloat x= (E-m);

      //PrecFloat phi= abs(atan( sin(s)/(1-cos(s))));
      
      //PrecFloat norm= PrecFloat(2)*phi/precPi() ;

      return (exp(-x)*sin(s))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));

    }

    if(SM_TYPE=="FF_Sinh_half") {

      PrecFloat x = (E-m);

      return exp(-x/2)*s/( pow(2*sinh(x/2),2) + s*s);

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

    PrecFloat t = E-m;
    PrecFloat x = sinh(t);
    return norm*s/( pow(x,2) + pow(s,2));

  };



  


  //loop over ensembles
  for(int iens=0; iens<Nens;iens++) {

    CorrAnalysis Corr;
    Corr.Nt = Nts[iens];
    

    
      

  }


  

 





  return return_class;

}
