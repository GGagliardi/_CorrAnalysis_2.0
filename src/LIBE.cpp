#include "../include/LIBE.h"
#include "numerics.h"
using namespace std;


const int Nboots = 800;
const bool UseJack=true;
const int Njacks = 25;
const vector<string> muval({"+7.2000e-04", "+2.1600e-03", "+3.6000e-03", "+5.0400e-03",  "+6.8400e-03","+1.8250e-02"});

const vector<string> muval_OS({"7.2000e-04", "2.1600e-03", "3.6000e-03", "5.0400e-03",  "6.8400e-03","1.8250e-02"});
const string musea = "+1.8250e-02";
const int mus = 6;
const double em2 = 1e6;
const double alpha = 1.0 / 137.035999;
const double qu = 2.0 / 3;
const double qd = -qu / 2;
const double MK_plus = 0.493677;
const double MK_zero = 0.497611;
const double DMK = -0.003934;
const double kappa = 2.837297;
const vector<double> mu_n({7.2000e-04, 2.1600e-03, 3.6000e-03, 5.0400e-03, 6.8400e-03,1.8250e-02});


void Compute_LIBE() {



  
  auto Sort_light_confs = [](string A, string B) {
    
    
    //return A<B;
    
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
    
  
    if(A_bis.length() <= 4) return A_bis < B_bis;
    string rA = A_bis.substr(A_bis.length()-2);
    string rB = B_bis.substr(B_bis.length()-2);
    
        
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(rA.substr(1,1));
      int n2 = stoi(rB.substr(1,1));
      if(rA == rB) {
	if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
	else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      
      else return n1<n2;
    }
    return A_bis<B_bis;
  };


  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");

  GaussianMersenne GM(78821);
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

   
    
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

  }
  
  cout<<"RC generated!"<<endl;


  vector<data_t> P5P5_ll_self_up(mus), P5P5_ll_self_down(mus), P5P5_ll_exchange(mus), P5P5_ll_dm_up(mus), P5P5_ll_dm_down(mus), P5P5_ll_dmcr_up(mus), P5P5_ll_dmcr_down(mus);
  vector<data_t> VKVK_ll_self_up(mus), VKVK_ll_self_down(mus), VKVK_ll_exchange(mus), VKVK_ll_dm_up(mus), VKVK_ll_dm_down(mus), VKVK_ll_dmcr_up(mus), VKVK_ll_dmcr_down(mus);
  vector<data_t> V0P5_ll_self_up(mus), V0P5_ll_self_down(mus), V0P5_ll_exchange(mus), V0P5_ll_dm_up(mus), V0P5_ll_dm_down(mus), V0P5_ll_dmcr_up(mus), V0P5_ll_dmcr_down(mus);
  vector<data_t> S0P5_ll_self_up(mus), S0P5_ll_self_down(mus), S0P5_ll_exchange(mus), S0P5_ll_dm_up(mus), S0P5_ll_dm_down(mus), S0P5_ll_dmcr_up(mus), S0P5_ll_dmcr_down(mus);
  vector<data_t> A0P5_ll_self_up(mus), A0P5_ll_self_down(mus), A0P5_ll_exchange(mus), A0P5_ll_dm_up(mus), A0P5_ll_dm_down(mus), A0P5_ll_dmcr_up(mus), A0P5_ll_dmcr_down(mus);
  

  vector<data_t> P5P5_ll_iso(mus), VKVK_ll_iso(mus), V0P5_ll_iso(mus), S0P5_ll_iso(mus), A0P5_ll_iso(mus);


  //OS correlators
  vector<data_t> OS_P5P5_ll_self_up(mus), OS_P5P5_ll_self_down(mus), OS_P5P5_ll_exchange(mus), OS_P5P5_ll_dm_up(mus), OS_P5P5_ll_dm_down(mus), OS_P5P5_ll_dmcr_up(mus), OS_P5P5_ll_dmcr_down(mus);
  vector<data_t> OS_A0P5_ll_self_up(mus), OS_A0P5_ll_self_down(mus), OS_A0P5_ll_exchange(mus), OS_A0P5_ll_dm_up(mus), OS_A0P5_ll_dm_down(mus), OS_A0P5_ll_dmcr_up(mus), OS_A0P5_ll_dmcr_down(mus);
  vector<data_t> OS_P5P5_ll_iso(mus), OS_A0P5_ll_iso(mus);
  
  vector<data_t> P5P5_ls_self_up(mus), P5P5_ls_self_down(mus), P5P5_ls_exchange(mus), P5P5_ls_dm_up(mus), P5P5_ls_dm_down(mus), P5P5_ls_dmcr_up(mus), P5P5_ls_dmcr_down(mus);
  vector<data_t> A0A0_ls_self_up(mus), A0A0_ls_self_down(mus), A0A0_ls_exchange(mus), A0A0_ls_dm_up(mus), A0A0_ls_dm_down(mus), A0A0_ls_dmcr_up(mus), A0A0_ls_dmcr_down(mus);
  vector<data_t> V0V0_ls_self_up(mus), V0V0_ls_self_down(mus), V0V0_ls_exchange(mus), V0V0_ls_dm_up(mus), V0V0_ls_dm_down(mus), V0V0_ls_dmcr_up(mus), V0V0_ls_dmcr_down(mus);
  vector<data_t> VKVK_ls_self_up(mus), VKVK_ls_self_down(mus), VKVK_ls_exchange(mus), VKVK_ls_dm_up(mus), VKVK_ls_dm_down(mus), VKVK_ls_dmcr_up(mus), VKVK_ls_dmcr_down(mus);
  vector<data_t> AKAK_ls_self_up(mus), AKAK_ls_self_down(mus), AKAK_ls_exchange(mus), AKAK_ls_dm_up(mus), AKAK_ls_dm_down(mus), AKAK_ls_dmcr_up(mus), AKAK_ls_dmcr_down(mus);
  vector<data_t> V0P5_ls_self_up(mus), V0P5_ls_self_down(mus), V0P5_ls_exchange(mus), V0P5_ls_dm_up(mus), V0P5_ls_dm_down(mus), V0P5_ls_dmcr_up(mus), V0P5_ls_dmcr_down(mus);

  vector<data_t> OS_P5P5_ls_self_up(mus), OS_P5P5_ls_self_down(mus), OS_P5P5_ls_exchange(mus), OS_P5P5_ls_dm_up(mus), OS_P5P5_ls_dm_down(mus), OS_P5P5_ls_dmcr_up(mus), OS_P5P5_ls_dmcr_down(mus);
  vector<data_t> OS_A0A0_ls_self_up(mus), OS_A0A0_ls_self_down(mus), OS_A0A0_ls_exchange(mus), OS_A0A0_ls_dm_up(mus), OS_A0A0_ls_dm_down(mus), OS_A0A0_ls_dmcr_up(mus), OS_A0A0_ls_dmcr_down(mus);
  vector<data_t> OS_V0V0_ls_self_up(mus), OS_V0V0_ls_self_down(mus), OS_V0V0_ls_exchange(mus), OS_V0V0_ls_dm_up(mus), OS_V0V0_ls_dm_down(mus), OS_V0V0_ls_dmcr_up(mus), OS_V0V0_ls_dmcr_down(mus);
  vector<data_t> OS_VKVK_ls_self_up(mus), OS_VKVK_ls_self_down(mus), OS_VKVK_ls_exchange(mus), OS_VKVK_ls_dm_up(mus), OS_VKVK_ls_dm_down(mus), OS_VKVK_ls_dmcr_up(mus), OS_VKVK_ls_dmcr_down(mus);
  vector<data_t> OS_AKAK_ls_self_up(mus), OS_AKAK_ls_self_down(mus), OS_AKAK_ls_exchange(mus), OS_AKAK_ls_dm_up(mus), OS_AKAK_ls_dm_down(mus), OS_AKAK_ls_dmcr_up(mus), OS_AKAK_ls_dmcr_down(mus);
  vector<data_t> OS_V0P5_ls_self_up(mus), OS_V0P5_ls_self_down(mus), OS_V0P5_ls_exchange(mus), OS_V0P5_ls_dm_up(mus), OS_V0P5_ls_dm_down(mus), OS_V0P5_ls_dmcr_up(mus), OS_V0P5_ls_dmcr_down(mus);

  vector<data_t> P5P5_ls_iso(mus), VKVK_ls_iso(mus), V0P5_ls_iso(mus), A0A0_ls_iso(mus), AKAK_ls_iso(mus);
  


  for(int i=0; i < (signed)muval.size(); i++) {

    //read files
    //ll
    P5P5_ll_self_up[i].Read("../LIBE_data", "self_plus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_self_down[i].Read("../LIBE_data", "self_minus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_exchange[i].Read("../LIBE_data", "exchange_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_dm_up[i].Read("../LIBE_data", "dm_plus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_dm_down[i].Read("../LIBE_data", "dm_minus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ll_iso[i].Read("../LIBE_data", "iso_ll_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);

    VKVK_ll_self_up[i].Read("../LIBE_data", "self_plus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_self_down[i].Read("../LIBE_data", "self_minus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_exchange[i].Read("../LIBE_data", "exchange_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_dm_up[i].Read("../LIBE_data", "dm_plus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_dm_down[i].Read("../LIBE_data", "dm_minus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ll_iso[i].Read("../LIBE_data", "iso_ll_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);

    V0P5_ll_self_up[i].Read("../LIBE_data", "self_plus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_self_down[i].Read("../LIBE_data", "self_minus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_exchange[i].Read("../LIBE_data", "exchange_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_dm_up[i].Read("../LIBE_data", "dm_plus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_dm_down[i].Read("../LIBE_data", "dm_minus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ll_iso[i].Read("../LIBE_data", "iso_ll_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);

    A0P5_ll_self_up[i].Read("../LIBE_data", "self_plus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_self_down[i].Read("../LIBE_data", "self_minus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_exchange[i].Read("../LIBE_data", "exchange_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_dm_up[i].Read("../LIBE_data", "dm_plus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_dm_down[i].Read("../LIBE_data", "dm_minus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);
    A0P5_ll_iso[i].Read("../LIBE_data", "iso_ll_muval_"+muval[i]+".dat",  "A0P5", Sort_light_confs);


    S0P5_ll_self_up[i].Read("../LIBE_data", "self_plus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_self_down[i].Read("../LIBE_data", "self_minus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_exchange[i].Read("../LIBE_data", "exchange_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_dm_up[i].Read("../LIBE_data", "dm_plus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_dm_down[i].Read("../LIBE_data", "dm_minus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);
    S0P5_ll_iso[i].Read("../LIBE_data", "iso_ll_muval_"+muval[i]+".dat",  "S0P5", Sort_light_confs);



    //ll OS
    OS_P5P5_ll_self_up[i].Read("../LIBE_data", "OS_self_plus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_self_down[i].Read("../LIBE_data", "OS_self_minus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_exchange[i].Read("../LIBE_data", "OS_exchange_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ll_iso[i].Read("../LIBE_data", "OS_iso_ll_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);

    
    OS_A0P5_ll_self_up[i].Read("../LIBE_data", "OS_self_plus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_self_down[i].Read("../LIBE_data", "OS_self_minus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_exchange[i].Read("../LIBE_data", "OS_exchange_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);
    OS_A0P5_ll_iso[i].Read("../LIBE_data", "OS_iso_ll_muval_"+muval_OS[i]+".dat",  "A0P5", Sort_light_confs);

    

    //ls
    P5P5_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);
    P5P5_ls_iso[i].Read("../LIBE_data", "iso_ls_muval_"+muval[i]+".dat",  "P5P5", Sort_light_confs);

    
    A0A0_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);
    A0A0_ls_iso[i].Read("../LIBE_data", "iso_ls_muval_"+muval[i]+".dat",  "A0A0", Sort_light_confs);

    V0V0_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
    V0V0_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "V0V0", Sort_light_confs);
   

    VKVK_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);
    VKVK_ls_iso[i].Read("../LIBE_data", "iso_ls_muval_"+muval[i]+".dat",  "VKVK", Sort_light_confs);


    AKAK_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    AKAK_ls_iso[i].Read("../LIBE_data", "iso_ls_muval_"+muval[i]+".dat",  "AKAK", Sort_light_confs);
    

    V0P5_ls_self_up[i].Read("../LIBE_data", "self_plus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_self_down[i].Read("../LIBE_data", "self_minus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_exchange[i].Read("../LIBE_data", "exchange_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_dm_up[i].Read("../LIBE_data", "dm_plus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_dm_down[i].Read("../LIBE_data", "dm_minus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_dmcr_up[i].Read("../LIBE_data", "dmcr_plus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_dmcr_down[i].Read("../LIBE_data", "dmcr_minus_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);
    V0P5_ls_iso[i].Read("../LIBE_data", "iso_ls_muval_"+muval[i]+".dat",  "V0P5", Sort_light_confs);


    //ls OS

    OS_P5P5_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    OS_P5P5_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "P5P5", Sort_light_confs);
    
    
    OS_A0A0_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
    OS_A0A0_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "A0A0", Sort_light_confs);
   
    OS_V0V0_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
    OS_V0V0_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "V0V0", Sort_light_confs);
   

    OS_VKVK_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
    OS_VKVK_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "VKVK", Sort_light_confs);
   

    OS_AKAK_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
    OS_AKAK_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "AKAK", Sort_light_confs);
   
    OS_V0P5_ls_self_up[i].Read("../LIBE_data", "OS_self_plus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_self_down[i].Read("../LIBE_data", "OS_self_minus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_exchange[i].Read("../LIBE_data", "OS_exchange_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_dm_up[i].Read("../LIBE_data", "OS_dm_plus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_dm_down[i].Read("../LIBE_data", "OS_dm_minus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_dmcr_up[i].Read("../LIBE_data", "OS_dmcr_plus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    OS_V0P5_ls_dmcr_down[i].Read("../LIBE_data", "OS_dmcr_minus_ls_muval_"+muval_OS[i]+".dat",  "V0P5", Sort_light_confs);
    

  }
  

  
  int Nens= P5P5_ll_iso[0].Tag.size();

  boost::filesystem::create_directory("../data/LIBE");
  

  //loop over ensembles
  for(int iens=0;iens<Nens;iens++) {

    string Ens_tag= P5P5_ll_iso[0].Tag[iens];
    cout<<"Analyzing ensemble: "<<Ens_tag<<endl;
    boost::filesystem::create_directory("../data/LIBE/"+Ens_tag);
    for(int imu=0;imu<mus;imu++) {
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]);
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr");
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope");
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/mass");
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN");
      boost::filesystem::create_directory("../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/amu");
    }
      
      
    
    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(P5P5_ll_iso[0].Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = L_info.T;
    double L = (double)L_info.L;
    double V= pow(L_info.L,3);
    


    //find ensemble info
    distr_t Za, Zv, a_distr;
    if(Ens_tag.substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A; a_distr=a_A;}
    else if(Ens_tag.substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; a_distr=a_B;}
    else if(Ens_tag.substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(Ens_tag.substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; a_distr=a_D;}
    else crash("Ensemble: "+Ens_tag+" not recognised");


    distr_t dmud_ph(UseJack), DM_ph(UseJack), dmcr_ph(UseJack), dmus_ph(UseJack);
    
    //loop over muval
    for(int imu=0;imu<mus;imu++) {

      double muiso= mu_n[imu];
      
      //iso
      distr_t_list P5P5_ll_iso_distr = Corr.corr_t( P5P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ll");
      distr_t_list OS_P5P5_ll_iso_distr = Corr.corr_t( OS_P5P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_P5P5_ll");
      distr_t_list VKVK_ll_iso_distr = Corr.corr_t( VKVK_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ll");      
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ll_iso_distr = -1.0*Corr.corr_t( V0P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ll");
      distr_t_list A0P5_ll_iso_distr = -1.0*Corr.corr_t( A0P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll");
      distr_t_list OS_A0P5_ll_iso_distr = -1.0*Corr.corr_t( OS_A0P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll_OS");
      Corr.Reflection_sign=1;
      //exchange
      distr_t_list P5P5_ll_exch_distr = Corr.corr_t( P5P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ll_exchange");
      distr_t_list VKVK_ll_exch_distr = Corr.corr_t( VKVK_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ll_exchange");
      distr_t_list OS_P5P5_ll_exch_distr = Corr.corr_t( OS_P5P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_P5P5_ll_exchange");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ll_exch_distr =  -1.0*Corr.corr_t( V0P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ll_exchange");
      distr_t_list A0P5_ll_exch_distr =  -1.0*Corr.corr_t( A0P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll_exchange");
      distr_t_list OS_A0P5_ll_exch_distr =  -1.0*Corr.corr_t( OS_A0P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0P5_ll_exchange");
      Corr.Reflection_sign=1;
      //self
      distr_t_list P5P5_ll_self_distr = Corr.corr_t( P5P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ll_self");
      distr_t_list VKVK_ll_self_distr = Corr.corr_t( VKVK_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ll_self");
      distr_t_list OS_P5P5_ll_self_distr = Corr.corr_t( OS_P5P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_P5P5_ll_self");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ll_self_distr =  -1.0*Corr.corr_t( V0P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ll_self");
      distr_t_list A0P5_ll_self_distr =  -1.0*Corr.corr_t( A0P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll_self");
      distr_t_list OS_A0P5_ll_self_distr =  -1.0*Corr.corr_t( OS_A0P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0P5_ll_self");
      Corr.Reflection_sign=1;
      //dm
      distr_t_list P5P5_ll_dm_distr = Corr.corr_t( P5P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ll_dm");
      distr_t_list VKVK_ll_dm_distr = Corr.corr_t( VKVK_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ll_dm");
      distr_t_list OS_P5P5_ll_dm_distr = Corr.corr_t( P5P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_P5P5_ll_dm");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ll_dm_distr =  -1.0*Corr.corr_t( V0P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ll_dm");
      distr_t_list A0P5_ll_dm_distr =  -1.0*Corr.corr_t( A0P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll_dm");
      distr_t_list OS_A0P5_ll_dm_distr =  -1.0*Corr.corr_t( OS_A0P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0P5_ll_dm");
      Corr.Reflection_sign=1;
      //dmcr
      distr_t_list P5P5_ll_dmcr_distr = Corr.corr_t( P5P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ll_dmcr");
      distr_t_list VKVK_ll_dmcr_distr = Corr.corr_t( VKVK_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ll_dmcr");
      distr_t_list OS_P5P5_ll_dmcr_distr = Corr.corr_t( OS_P5P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_P5P5_ll_dmcr");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ll_dmcr_distr =  -1.0*Corr.corr_t( V0P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ll_dmcr");
      distr_t_list A0P5_ll_dmcr_distr =  -1.0*Corr.corr_t( A0P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0P5_ll_dmcr");
      distr_t_list OS_A0P5_ll_dmcr_distr =  -1.0*Corr.corr_t( OS_A0P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0P5_ll_dmcr");
      Corr.Reflection_sign=1;


      //S0P5 #################
      Corr.Perform_Nt_t_average=0;
      distr_t_list S0P5_ll_iso_distr = Corr.corr_t( S0P5_ll_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/S0P5_ll");
      distr_t_list S0P5_ll_exch_distr = Corr.corr_t( S0P5_ll_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/S0P5_ll_exchange");
      distr_t_list S0P5_ll_self_distr = Corr.corr_t( S0P5_ll_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/S0P5_ll_self");
      distr_t_list S0P5_ll_dm_distr = Corr.corr_t( S0P5_ll_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/S0P5_ll_dm");
      distr_t_list S0P5_ll_dmcr_distr = Corr.corr_t( S0P5_ll_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/S0P5_ll_dmcr");
      Corr.Perform_Nt_t_average=1;
      //######################
      

      //ls
      
      //iso
      distr_t_list P5P5_ls_iso_distr = Corr.corr_t( P5P5_ls_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls");
      distr_t_list VKVK_ls_iso_distr = Corr.corr_t( VKVK_ls_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls");
      distr_t_list A0A0_ls_iso_distr = Corr.corr_t( A0A0_ls_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls");
      distr_t_list AKAK_ls_iso_distr = Corr.corr_t( AKAK_ls_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_iso_distr =  -1.0*Corr.corr_t( V0P5_ls_iso[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls");
      Corr.Reflection_sign=1;
      //exchange
      distr_t_list P5P5_ls_exch_distr = Corr.corr_t( P5P5_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_exchange");
      distr_t_list VKVK_ls_exch_distr = Corr.corr_t( VKVK_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_exchange");
      distr_t_list A0A0_ls_exch_distr = Corr.corr_t( A0A0_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_exchange");
      distr_t_list V0V0_ls_exch_distr = Corr.corr_t( V0V0_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_exchange");
      distr_t_list AKAK_ls_exch_distr = Corr.corr_t( AKAK_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_exchange");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_exch_distr =  -1.0*Corr.corr_t( V0P5_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_exchange");
      Corr.Reflection_sign=1;
      //self up
      distr_t_list P5P5_ls_self_up_distr = Corr.corr_t( P5P5_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_self_up");
      distr_t_list VKVK_ls_self_up_distr = Corr.corr_t( VKVK_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_self_up");
      distr_t_list A0A0_ls_self_up_distr = Corr.corr_t( A0A0_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_self_up");
      distr_t_list V0V0_ls_self_up_distr = Corr.corr_t( V0V0_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_self_up");
      distr_t_list AKAK_ls_self_up_distr = Corr.corr_t( AKAK_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_self_up");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_self_up_distr =  -1.0*Corr.corr_t( V0P5_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_self_up");
      Corr.Reflection_sign=1;
      //self down
      distr_t_list P5P5_ls_self_down_distr = Corr.corr_t( P5P5_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_self_down");
      distr_t_list VKVK_ls_self_down_distr = Corr.corr_t( VKVK_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_self_down");
      distr_t_list A0A0_ls_self_down_distr = Corr.corr_t( A0A0_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_self_down");
      distr_t_list V0V0_ls_self_down_distr = Corr.corr_t( V0V0_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_self_down");
      distr_t_list AKAK_ls_self_down_distr = Corr.corr_t( AKAK_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_self_down");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_self_down_distr =  -1.0*Corr.corr_t( V0P5_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_self_down");
      Corr.Reflection_sign=1;
      
      //dm up
      distr_t_list P5P5_ls_dm_up_distr = Corr.corr_t( P5P5_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_dm_up");
      distr_t_list VKVK_ls_dm_up_distr = Corr.corr_t( VKVK_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_dm_up");

      distr_t_list A0A0_ls_dm_up_distr = Corr.corr_t( A0A0_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_dm_up");
      distr_t_list V0V0_ls_dm_up_distr = Corr.corr_t( V0V0_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_dm_up");
      distr_t_list AKAK_ls_dm_up_distr = Corr.corr_t( AKAK_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_dm_up");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_dm_up_distr =  -1.0*Corr.corr_t( V0P5_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_dm_up");
      Corr.Reflection_sign=1;
      //dm down
      distr_t_list P5P5_ls_dm_down_distr = Corr.corr_t( P5P5_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_dm_down");
      distr_t_list VKVK_ls_dm_down_distr = Corr.corr_t( VKVK_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_dm_down");
      distr_t_list V0V0_ls_dm_down_distr = Corr.corr_t( V0V0_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_dm_down");
      distr_t_list A0A0_ls_dm_down_distr = Corr.corr_t( A0A0_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_dm_down");
      distr_t_list AKAK_ls_dm_down_distr = Corr.corr_t( AKAK_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_dm_down");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_dm_down_distr =  -1.0*Corr.corr_t( V0P5_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_dm_down");
      Corr.Reflection_sign=1;
      
      //dmcr up
      distr_t_list P5P5_ls_dmcr_up_distr = Corr.corr_t( P5P5_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_dmcr_up");
      distr_t_list VKVK_ls_dmcr_up_distr = Corr.corr_t( VKVK_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_dmcr_up");
      distr_t_list V0V0_ls_dmcr_up_distr = Corr.corr_t( V0V0_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_dmcr_up");
      distr_t_list A0A0_ls_dmcr_up_distr = Corr.corr_t( A0A0_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_dmcr_up");
      distr_t_list AKAK_ls_dmcr_up_distr = Corr.corr_t( AKAK_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_dmcr_up");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_dmcr_up_distr =  -1.0*Corr.corr_t( V0P5_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_dmcr_up");
      Corr.Reflection_sign=1;

      //dmcr down
      distr_t_list P5P5_ls_dmcr_down_distr = Corr.corr_t( P5P5_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/P5P5_ls_dmcr_down");
      distr_t_list VKVK_ls_dmcr_down_distr = Corr.corr_t( VKVK_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/VKVK_ls_dmcr_down");
      distr_t_list V0V0_ls_dmcr_down_distr = Corr.corr_t( V0V0_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0V0_ls_dmcr_down");
      distr_t_list A0A0_ls_dmcr_down_distr = Corr.corr_t( A0A0_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/A0A0_ls_dmcr_down");
      distr_t_list AKAK_ls_dmcr_down_distr = Corr.corr_t( AKAK_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/AKAK_ls_dmcr_down");
      Corr.Reflection_sign=-1;
      distr_t_list V0P5_ls_dmcr_down_distr =  -1.0*Corr.corr_t( V0P5_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/V0P5_ls_dmcr_down");
      Corr.Reflection_sign=1;

      //OS
      //#############################################################################################################


     
      //exchange
      distr_t_list OS_VKVK_ls_exch_distr = Corr.corr_t( OS_VKVK_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_exchange");
      distr_t_list OS_A0A0_ls_exch_distr = Corr.corr_t( OS_A0A0_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_exchange");
      distr_t_list OS_V0V0_ls_exch_distr = Corr.corr_t( OS_V0V0_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_exchange");
      distr_t_list OS_AKAK_ls_exch_distr = Corr.corr_t( OS_AKAK_ls_exchange[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_exchange");
      //self up
      distr_t_list OS_VKVK_ls_self_up_distr = Corr.corr_t( OS_VKVK_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_self_up");
      distr_t_list OS_A0A0_ls_self_up_distr = Corr.corr_t( OS_A0A0_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_self_up");
      distr_t_list OS_V0V0_ls_self_up_distr = Corr.corr_t( OS_V0V0_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_self_up");
      distr_t_list OS_AKAK_ls_self_up_distr = Corr.corr_t( OS_AKAK_ls_self_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_self_up");
      //self down
      distr_t_list OS_VKVK_ls_self_down_distr = Corr.corr_t( OS_VKVK_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_self_down");
      distr_t_list OS_A0A0_ls_self_down_distr = Corr.corr_t( OS_A0A0_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_self_down");
      distr_t_list OS_V0V0_ls_self_down_distr = Corr.corr_t( OS_V0V0_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_self_down");
      distr_t_list OS_AKAK_ls_self_down_distr = Corr.corr_t( OS_AKAK_ls_self_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_self_down");
      //dm up
      distr_t_list OS_VKVK_ls_dm_up_distr = Corr.corr_t( OS_VKVK_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_dm_up");
      distr_t_list OS_A0A0_ls_dm_up_distr = Corr.corr_t( OS_A0A0_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_dm_up");
      distr_t_list OS_V0V0_ls_dm_up_distr = Corr.corr_t( OS_V0V0_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_dm_up");
      distr_t_list OS_AKAK_ls_dm_up_distr = Corr.corr_t( OS_AKAK_ls_dm_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_dm_up");
      //dm down
      distr_t_list OS_VKVK_ls_dm_down_distr = Corr.corr_t( OS_VKVK_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_dm_down");
      distr_t_list OS_V0V0_ls_dm_down_distr = Corr.corr_t( OS_V0V0_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_dm_down");
      distr_t_list OS_A0A0_ls_dm_down_distr = Corr.corr_t( OS_A0A0_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_dm_down");
      distr_t_list OS_AKAK_ls_dm_down_distr = Corr.corr_t( OS_AKAK_ls_dm_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_dm_down");
      
      
      //dmcr up
      distr_t_list OS_VKVK_ls_dmcr_up_distr = Corr.corr_t( OS_VKVK_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_dmcr_up");
      distr_t_list OS_V0V0_ls_dmcr_up_distr = Corr.corr_t( OS_V0V0_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_dmcr_up");
      distr_t_list OS_A0A0_ls_dmcr_up_distr = Corr.corr_t( OS_A0A0_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_dmcr_up");
      distr_t_list OS_AKAK_ls_dmcr_up_distr = Corr.corr_t( OS_AKAK_ls_dmcr_down[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_dmcr_up");
      
      //dmcr down
      distr_t_list OS_VKVK_ls_dmcr_down_distr = Corr.corr_t( OS_VKVK_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_VKVK_ls_dmcr_down");
      distr_t_list OS_V0V0_ls_dmcr_down_distr = Corr.corr_t( OS_V0V0_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_V0V0_ls_dmcr_down");
      distr_t_list OS_A0A0_ls_dmcr_down_distr = Corr.corr_t( OS_A0A0_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_A0A0_ls_dmcr_down");
      distr_t_list OS_AKAK_ls_dmcr_down_distr = Corr.corr_t( OS_AKAK_ls_dmcr_up[imu].col(0)[iens], "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/OS_AKAK_ls_dmcr_down");
      



      //###############################################################################################################


      


      distr_t_list pion_distr= Corr.effective_mass_t( P5P5_ll_iso_distr, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/mass/pion");
      distr_t_list kaon_distr= Corr.effective_mass_t( P5P5_ls_iso_distr, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/mass/kaon");

      if(Ens_tag=="cB211b.072.48") {
	Corr.Tmin=20;
	Corr.Tmax=35;
      }
      else {
	Corr.Tmin=22;
	Corr.Tmax=45;
      }

      distr_t Mpi=Corr.Fit_distr(pion_distr);
      distr_t M_K= Corr.Fit_distr(kaon_distr);
      distr_t Delta_Mpi_iso= 0.13957039*a_distr - 0.135*a_distr ;
      distr_t Delta_MK_iso=  (MK_plus*a_distr - 0.4946*a_distr);
      cout<<"Delta_Mpi_iso: "<<(Delta_Mpi_iso/a_distr).ave()<<" +- "<<(Delta_Mpi_iso/a_distr).err()<<endl;
      cout<<"Mpi_iso: "<<(Mpi/a_distr).ave()<<" +- "<<(Mpi/a_distr).err()<<endl;
      cout<<"MK_iso: "<<(M_K/a_distr).ave()<<" +- "<<(M_K/a_distr).err()<<endl;
      
      //compute pion mass splitting
      distr_t_list Delta_pion_distr= Corr.effective_slope_t( (4*M_PI)*0.5*alpha*em2*P5P5_ll_exch_distr, P5P5_ll_iso_distr , "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/pion_splitting" );

      distr_t_list Delta_pion_subL_distr= Delta_pion_distr + alpha*(2*kappa/(2.0*Mpi*pow(L,2)))*( 1.0 + 0.5*L*Mpi);

      distr_t FSE_charged_pion= alpha*(2*kappa/(2.0*Mpi*pow(L,2)))*( 1.0 + 0.5*L*Mpi);  

      Print_To_File({}, {Delta_pion_distr.ave(), Delta_pion_distr.err(), Delta_pion_subL_distr.ave(), Delta_pion_subL_distr.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/pion_splitting_subL", "", "");
      
      //compute dmcrit
      //distr_t_list dmcr_distr= em2*4.0*M_PI*alpha*( (qu*qu + qd*qd)*V0P5_ll_self_distr )/(V0P5_ll_dmcr_distr) ;

      //compute dmcrit
      distr_t_list FACT_DM_QED= em2*4.0*M_PI*alpha*( distr_t_list::derivative((qu*qu+qd*qd)*V0P5_ll_self_distr + qu*qd*V0P5_ll_exch_distr,0)/distr_t_list::derivative(V0P5_ll_iso_distr,0) -  ((qu*qu+qd*qd)*P5P5_ll_self_distr + qu*qd*P5P5_ll_exch_distr)/P5P5_ll_iso_distr);
      distr_t_list FACT_DM_P = distr_t_list::derivative(V0P5_ll_dmcr_distr,0)/distr_t_list::derivative(V0P5_ll_iso_distr,0) - P5P5_ll_dmcr_distr/P5P5_ll_iso_distr;
      distr_t_list dmcr_distr = FACT_DM_QED/FACT_DM_P;
      
      //cout<<"conv fact: "<<(1/(a_distr.ave()*1e-3))<<endl;

      //dmcr_distr= em2*4.0*M_PI*alpha*( (qu*qu + qd*qd)*V0P5_ll_self_distr )/(V0P5_ll_dmcr_distr) ;
      
      Print_To_File( {}, {dmcr_distr.ave(), dmcr_distr.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/dmcr.eff", "", "");

      if(Ens_tag=="cB211b.072.48") {
	Corr.Tmin = 12;
	Corr.Tmax = 25;
      }
      else {
	Corr.Tmin = 14;
	Corr.Tmax = 30;
      }

      distr_t dmcr = Corr.Fit_distr( dmcr_distr);

      //compute DM from Mpi^+
      distr_t_list DM_distr = 0.5*(em2*4.0*M_PI*alpha*( (qu*qu+qd*qd)*Corr.effective_slope_t(P5P5_ll_self_distr,P5P5_ll_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/self_pi_plus")     +qu*qd*Corr.effective_slope_t(P5P5_ll_exch_distr, P5P5_ll_iso_distr , "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/exchange_pi_plus")  ) -dmcr*Corr.effective_slope_t(P5P5_ll_dmcr_distr, P5P5_ll_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dmcr_pi_plus" ) -FSE_charged_pion   -Delta_Mpi_iso)/Corr.effective_slope_t( P5P5_ll_dm_distr, P5P5_ll_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dm_pi_plus");

      distr_t_list DM_distr_fake_QED = 0.5*(em2*4.0*M_PI*alpha*( (qu*qu+qu*qu)*Corr.effective_slope_t(P5P5_ll_self_distr,P5P5_ll_iso_distr,  "")     +qu*qu*Corr.effective_slope_t(P5P5_ll_exch_distr, P5P5_ll_iso_distr , "")  ) -2.0*dmcr*(4.0/5.0)*Corr.effective_slope_t(P5P5_ll_dmcr_distr, P5P5_ll_iso_distr,  "" ))/Corr.effective_slope_t( P5P5_ll_dm_distr, P5P5_ll_iso_distr, "");

      Print_To_File( {}, {DM_distr.ave(), DM_distr.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/DM.eff", "", "");

      Print_To_File( {}, {DM_distr_fake_QED.ave(), DM_distr_fake_QED.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/DM_fQED.eff", "", "");


      //compute dmud from mK^+ - mK^0
      distr_t_list dmud_distr= 0.5*(em2*4.0*M_PI*alpha*( (qu*qu-qd*qd)*( Corr.effective_slope_t(P5P5_ls_self_up_distr,P5P5_ls_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/self_up_K_plus")) + (qu*qd-qd*qd)*(Corr.effective_slope_t(P5P5_ls_exch_distr,P5P5_ls_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/exchange_K_plus"))) -(dmcr*3.0/5.0)*Corr.effective_slope_t(P5P5_ls_dmcr_up_distr, P5P5_ls_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dmcr_up_K_plus" )  +DMK*a_distr)/Corr.effective_slope_t( P5P5_ls_dm_up_distr, P5P5_ls_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dm_up_K_plus");

      
      Print_To_File( {}, {dmud_distr.ave(), dmud_distr.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/dmud.eff", "", "");

      if(Ens_tag=="cB211b.072.48") {
	Corr.Tmin=20;
	Corr.Tmax=35;
      }
      else {
	Corr.Tmin=22;
	Corr.Tmax=45;
      }

      distr_t DM= Corr.Fit_distr(DM_distr);
      distr_t DM_fQED= Corr.Fit_distr(DM_distr_fake_QED);
      distr_t dmud = Corr.Fit_distr(dmud_distr);


      //compute dmus

      distr_t_list dms_distr = -1.0*(  em2*4.0*M_PI*alpha*( qu*qu*Corr.effective_slope_t( P5P5_ls_self_up_distr, P5P5_ls_iso_distr, "") + qd*qd*Corr.effective_slope_t( P5P5_ls_self_down_distr, P5P5_ls_iso_distr,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/self_down_K_plus") + qu*qd*Corr.effective_slope_t(P5P5_ls_exch_distr,P5P5_ls_iso_distr,""))  -dmcr*(4.0/5.0)*Corr.effective_slope_t(P5P5_ls_dmcr_up_distr, P5P5_ls_iso_distr,  "" ) -dmcr*(1.0/5.0)*Corr.effective_slope_t(P5P5_ls_dmcr_down_distr, P5P5_ls_iso_distr, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dmcr_down_K_plus") - (DM -dmud)*Corr.effective_slope_t( P5P5_ls_dm_up_distr, P5P5_ls_iso_distr,  "")	+ Delta_MK_iso)/Corr.effective_slope_t( P5P5_ls_dm_down_distr, P5P5_ls_iso_distr,   "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/eff_slope/dm_down_K_plus")			       ;		


      Print_To_File( {}, {dms_distr.ave(), dms_distr.err()}, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/dms.eff", "", "");
     
      distr_t dmus = Corr.Fit_distr(dms_distr);

      if( imu == 0) { dmud_ph = dmud; DM_ph=DM; dmcr_ph = dmcr; dmus_ph=dmus;}

         
   
   
      //VKVK ll correlator

      distr_t_list dVK_up= -1.0*(Za*Za/V)*pow(qu,2)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*VKVK_ll_self_distr + pow(qu,2)*VKVK_ll_exch_distr ) +2*(DM_ph - dmud_ph)*VKVK_ll_dm_distr -2*(dmcr_ph*4.0/5.0)*VKVK_ll_dmcr_distr ) ;
      distr_t_list dVK_down = -1.0*(Za*Za/V)*pow(qd,2)*( em2*4.0*M_PI*alpha*( 2*pow(qd,2)*VKVK_ll_self_distr + pow(qd,2)*VKVK_ll_exch_distr ) +2*(DM_ph + dmud_ph)*VKVK_ll_dm_distr -2*(dmcr_ph/5.0)*VKVK_ll_dmcr_distr );
      distr_t_list dVK = dVK_up + dVK_down;


      //distr_t_list OS_dVK_up= -1.0*(Zv*Zv/V)*pow(qu,2)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*VKVK_ll_self_distr + pow(qu,2)*VKVK_ll_exch_distr ) +2*(DM_ph - dmud_ph)*VKVK_ll_dm_distr -2*(dmcr_ph*4.0/5.0)*VKVK_ll_dmcr_distr ) ;
      //distr_t_list OS_dVK_down = -1.0*(Zv*Zv/V)*pow(qd,2)*( em2*4.0*M_PI*alpha*( 2*pow(qd,2)*VKVK_ll_self_distr + pow(qd,2)*VKVK_ll_exch_distr ) +2*(DM_ph + dmud_ph)*VKVK_ll_dm_distr -2*(dmcr_ph/5.0)*VKVK_ll_dmcr_distr );
      //distr_t_list OS_dVK = OS_dVK_up + OS_dVK_down;

      distr_t_list dVK_up_QED = -1.0*(Za*Za/V)*pow(qu,2)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*VKVK_ll_self_distr + pow(qu,2)*VKVK_ll_exch_distr ) -2*(dmcr_ph*4.0/5.0)*VKVK_ll_dmcr_distr);
      distr_t_list dVK_up_SIB = -1.0*(Za*Za/V)*pow(qu,2)*( -2*(DM_ph - dmud_ph)*VKVK_ll_dm_distr  );


      //for tau decay VKVK, AKAK , A0A0 ls correlators

      distr_t_list dVK_ls= -1.0*(Za*Za/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*VKVK_ls_self_up_distr + pow(qd,2)*VKVK_ls_self_down_distr  + qu*qd*VKVK_ls_exch_distr )  +(DM_ph - dmud_ph)*VKVK_ls_dm_up_distr -dmus_ph*VKVK_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*VKVK_ls_dmcr_up_distr -(dmcr_ph/5.0)*VKVK_ls_dmcr_down_distr ) ;
      
      distr_t_list dAK_ls= -1.0*(Zv*Zv/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*AKAK_ls_self_up_distr + pow(qd,2)*AKAK_ls_self_down_distr  + qu*qd*AKAK_ls_exch_distr ) +(DM_ph - dmud_ph)*AKAK_ls_dm_up_distr -dmus_ph*AKAK_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*AKAK_ls_dmcr_up_distr -(dmcr_ph/5.0)*AKAK_ls_dmcr_down_distr ) ;
      
      distr_t_list dA0_ls= 1.0*(Zv*Zv/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*A0A0_ls_self_up_distr + pow(qd,2)*A0A0_ls_self_down_distr  + qu*qd*A0A0_ls_exch_distr ) +(DM_ph - dmud_ph)*A0A0_ls_dm_up_distr -dmus_ph*A0A0_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*A0A0_ls_dmcr_up_distr -(dmcr_ph/5.0)*A0A0_ls_dmcr_down_distr ) ;

      distr_t_list dV0_ls= 1.0*(Za*Za/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*V0V0_ls_self_up_distr + pow(qd,2)*V0V0_ls_self_down_distr  + qu*qd*V0V0_ls_exch_distr ) +(DM_ph - dmud_ph)*V0V0_ls_dm_up_distr -dmus_ph*V0V0_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*V0V0_ls_dmcr_up_distr -(dmcr_ph/5.0)*V0V0_ls_dmcr_down_distr ) ;

      
      Print_To_File( {}, { dVK_ls.ave(), dVK_ls.err(), dAK_ls.ave(), dAK_ls.err(), dA0_ls.ave(), dA0_ls.err(), dV0_ls.ave(), dV0_ls.err()} , "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/tau_Vus_QED", "", "" );


      distr_t_list OS_dVK_ls= -1.0*(Zv*Zv/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*OS_VKVK_ls_self_up_distr + pow(qd,2)*OS_VKVK_ls_self_down_distr  + qu*qd*OS_VKVK_ls_exch_distr )  +(DM_ph - dmud_ph)*OS_VKVK_ls_dm_up_distr -dmus_ph*OS_VKVK_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*OS_VKVK_ls_dmcr_up_distr -(dmcr_ph/5.0)*OS_VKVK_ls_dmcr_down_distr ) ;
      
      distr_t_list OS_dAK_ls= -1.0*(Za*Za/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*OS_AKAK_ls_self_up_distr + pow(qd,2)*OS_AKAK_ls_self_down_distr  + qu*qd*OS_AKAK_ls_exch_distr ) +(DM_ph - dmud_ph)*OS_AKAK_ls_dm_up_distr -dmus_ph*OS_AKAK_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*OS_AKAK_ls_dmcr_up_distr -(dmcr_ph/5.0)*OS_AKAK_ls_dmcr_down_distr ) ;
      
      distr_t_list OS_dA0_ls= 1.0*(Za*Za/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*OS_A0A0_ls_self_up_distr + pow(qd,2)*OS_A0A0_ls_self_down_distr  + qu*qd*OS_A0A0_ls_exch_distr ) +(DM_ph - dmud_ph)*OS_A0A0_ls_dm_up_distr -dmus_ph*OS_A0A0_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*OS_A0A0_ls_dmcr_up_distr -(dmcr_ph/5.0)*OS_A0A0_ls_dmcr_down_distr ) ;

      distr_t_list OS_dV0_ls= 1.0*(Zv*Zv/V)*( em2*4.0*M_PI*alpha*( pow(qu,2)*OS_V0V0_ls_self_up_distr + pow(qd,2)*OS_V0V0_ls_self_down_distr  + qu*qd*OS_V0V0_ls_exch_distr ) +(DM_ph - dmud_ph)*OS_V0V0_ls_dm_up_distr -dmus_ph*OS_V0V0_ls_dm_down_distr  -(dmcr_ph*4.0/5.0)*OS_V0V0_ls_dmcr_up_distr -(dmcr_ph/5.0)*OS_V0V0_ls_dmcr_down_distr ) ;

      
      Print_To_File( {}, { OS_dVK_ls.ave(), OS_dVK_ls.err(), OS_dAK_ls.ave(), OS_dAK_ls.err(), OS_dA0_ls.ave(), OS_dA0_ls.err(), OS_dV0_ls.ave(), OS_dV0_ls.err()} , "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/corr/tau_Vus_QED_OS", "", "" );
      


      // A0P5 and P5P5 correlators to estract Zv

      distr_t_list dA0P5_ll =  (-1.0/V)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*A0P5_ll_self_distr + pow(qu,2)*A0P5_ll_exch_distr ) -2*(DM-dmud)*A0P5_ll_dm_distr -2*(dmcr*4.0/5.0)*A0P5_ll_dmcr_distr ) ;

      distr_t_list dP5P5_ll =  (-1.0/V)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*P5P5_ll_self_distr + pow(qu,2)*P5P5_ll_exch_distr ) -2*(DM-dmud)*P5P5_ll_dm_distr -2*(dmcr*4.0/5.0)*P5P5_ll_dmcr_distr ) ;

      distr_t_list dt_A0P5_ll = distr_t_list::derivative(-1.0*A0P5_ll_iso_distr/V, 0);
      distr_t_list P5P5_ll = -1.0*P5P5_ll_iso_distr/V;
      distr_t_list dt_dA0P5_ll = distr_t_list::derivative( dA0P5_ll, 0);


      //OS
      distr_t_list OS_dA0P5_ll =  (-1.0/V)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*OS_A0P5_ll_self_distr + pow(qu,2)*OS_A0P5_ll_exch_distr ) -2*(DM-dmud)*OS_A0P5_ll_dm_distr -2*(dmcr*4.0/5.0)*OS_A0P5_ll_dmcr_distr ) ;

      distr_t_list OS_dP5P5_ll =  (-1.0/V)*( em2*4.0*M_PI*alpha*( 2*pow(qu,2)*OS_P5P5_ll_self_distr + pow(qu,2)*OS_P5P5_ll_exch_distr ) -2*(DM-dmud)*OS_P5P5_ll_dm_distr -2*(dmcr*4.0/5.0)*OS_P5P5_ll_dmcr_distr ) ;

      distr_t_list OS_dt_A0P5_ll = distr_t_list::derivative(-1.0*OS_A0P5_ll_iso_distr/V, 0);
      distr_t_list OS_P5P5_ll = -1.0*OS_P5P5_ll_iso_distr/V;
      distr_t_list OS_dt_dA0P5_ll = distr_t_list::derivative( OS_dA0P5_ll, 0);

   


      distr_t_list RV_QCED = 2*muiso*dP5P5_ll/dt_A0P5_ll - 2*muiso*P5P5_ll*dt_dA0P5_ll/(dt_A0P5_ll*dt_A0P5_ll) + 2*(DM-dmud)*P5P5_ll/dt_A0P5_ll ;











      
      
      distr_t_list A0P5_QCED_ll = -1.0*A0P5_ll_iso_distr/V + dA0P5_ll;
      distr_t_list P5P5_QCED_ll = -1.0*P5P5_ll_iso_distr/V + dP5P5_ll;

      //distr_t_list RV_QCED = (1.0 + (DM-dmud)/muiso)*P5P5_QCED_ll/(distr_t_list::derivative(A0P5_QCED_ll, 0))/( P5P5_ll_iso_distr/distr_t_list::derivative(A0P5_ll_iso_distr, 0));

      distr_t_list RV_isoQCD=  2.0*muiso*P5P5_ll_iso_distr/distr_t_list::derivative(A0P5_ll_iso_distr, 0);

      RV_QCED= RV_QCED/RV_isoQCD;

      Print_To_File({}, {RV_QCED.ave(), RV_QCED.err(), RV_isoQCD.ave(), RV_isoQCD.err() }, "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/REN/RV", "", "");
      
						      

         
      //compute HVP
      auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
      distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);

   
      distr_t_list Int_dVK(UseJack);
      for(int t=0;t<Corr.Nt/2;t++) Int_dVK.distr_list.push_back( dVK.distr_list[t]*4.0*pow(alpha,2)*w(t,1)*Ker[t]*1e10);

      
      Print_To_File({}, {Int_dVK.ave(), Int_dVK.err() } , "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/amu/Int_TM", "" , "");
      Print_To_File({}, {dVK.ave(), dVK.err(), dVK_up_QED.ave(), dVK_up_QED.err(), dVK_up_SIB.ave(), dVK_up_SIB.err() } ,  "../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/amu/dVK_TM", "" , "");



      //

      

      
      //get amu^HVP
      
      //check dm corrections
      int tmax_15 = 1.5/(a_distr.ave()/fmTGeV);
      int tmax_20 = 2.0/(a_distr.ave()/fmTGeV);
      int tmax_25 = 2.5/(a_distr.ave()/fmTGeV);
      int tmax_30 = 3.0/(a_distr.ave()/fmTGeV);
      int tmax_35 = 3.5/(a_distr.ave()/fmTGeV);
         
      distr_t HVP_LIBE_tm_15(UseJack, UseJack?Njacks:800) ; distr_t HVP_LIBE_tm_20(UseJack, UseJack?Njacks:800); distr_t HVP_LIBE_tm_25(UseJack, UseJack?Njacks:800);   distr_t HVP_LIBE_tm_30(UseJack, UseJack?Njacks:800); distr_t HVP_LIBE_tm_35(UseJack, UseJack?Njacks:800);

      for(int t=1; t<tmax_15;t++) {
	HVP_LIBE_tm_15 = HVP_LIBE_tm_15 + Int_dVK.distr_list[t];
      }

      for(int t=1; t<tmax_20;t++) {
	HVP_LIBE_tm_20 = HVP_LIBE_tm_20 + Int_dVK.distr_list[t];
      }
           
      for(int t=1; t<tmax_25;t++) {
	HVP_LIBE_tm_25 = HVP_LIBE_tm_25 + Int_dVK.distr_list[t];
      }
      
      for(int t=1; t<tmax_30;t++) {
	HVP_LIBE_tm_30 = HVP_LIBE_tm_30 +  Int_dVK.distr_list[t];
      }
      
      for(int t=1; t<tmax_35;t++) {
	HVP_LIBE_tm_35 = HVP_LIBE_tm_35 +  Int_dVK.distr_list[t];
      }


      //bounding attempt
      distr_t p2_mot= 2*SQRT_D( Mpi*Mpi + pow( 2*M_PI/L,2));
      distr_t_list VKVK_iso = -(qu*qu+ qd*qd)*1e10*Za*Za*VKVK_ll_iso_distr/V;
      distr_t_list VKVK_QCED = VKVK_iso + 1e10*dVK;

      distr_t amu_HVP_iso;
      distr_t amu_HVP_QCED;

      int Tcut_opt_iso, Tcut_opt_QCED;
      
      //Bounding_HVP(amu_HVP_iso, Tcut_opt_iso,  VKVK_iso, a_distr,"../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/amu/bounding_iso" , p2_mot);
      //Bounding_HVP(amu_HVP_QCED, Tcut_opt_QCED,  VKVK_QCED, a_distr,"../data/LIBE/"+Ens_tag+"/muval_"+muval[imu]+"/amu/bounding_QCED" , Mpi);

      distr_t amu_HVP_dQCED = amu_HVP_QCED - amu_HVP_iso;

      
      
      
      
      
      cout<<"####  EVALUATING LIBE CORRECTION FOR HVP ENSEMBLE: "<<Ens_tag<<" #######"<<endl;
      cout<<"tcut 1.5fm: "<<HVP_LIBE_tm_15.ave()<<" +- "<<HVP_LIBE_tm_15.err()<<" [TM]"<<endl;
      cout<<"tcut 2.0fm: "<<HVP_LIBE_tm_20.ave()<<" +- "<<HVP_LIBE_tm_20.err()<<" [TM]"<<endl;
      cout<<"tcut 2.5fm: "<<HVP_LIBE_tm_25.ave()<<" +- "<<HVP_LIBE_tm_25.err()<<" [TM]"<<endl;
      cout<<"tcut 3.0fm: "<<HVP_LIBE_tm_30.ave()<<" +- "<<HVP_LIBE_tm_30.err()<<" [TM]"<<endl;
      cout<<"tcut 3.5fm: "<<HVP_LIBE_tm_35.ave()<<" +- "<<HVP_LIBE_tm_35.err()<<" [TM]"<<endl;
      //cout<<"Bounding: "<<amu_HVP_dQCED.ave()<<" +- "<<amu_HVP_dQCED.err()<<" [TM]"<<endl;
      cout<<"####    DONE! ####"<<endl;
      
      
      

    }
  }

  return; 

}
