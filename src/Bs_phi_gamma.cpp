#include "../include/Bs_phi_gamma.h"
#include "numerics.h"
using namespace std;


const int Nboots = 800;
const bool UseJack=true;
const int Njacks=20;
const double Q = -1.0/3.0; //electric charge of d-type quark
const double Lambda_QCD= 0.3; //300 MeV
const int Nmasses = 5;
const int Nmasses_triv=3;
const int Nsou=3;
const Vfloat ratio_mh({1, 1.0 / 1.4835, 1.0 / 2.02827, 1.0 / 2.53531, 1.0 / 3.04236});
// const Vfloat ratio_mh({ 1.0/1.4835, 1.0/2.02827, 1.0/2.53531, 1.0/3.04236});
const double sign_kz=-1;



void Compute_Bs_phi_gamma() {


  
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

  //resample Jpsi tensor decay constant
  distr_t fT_Jpsi_distr(UseJack);

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
  
  vector<vector<data_t>> PT3_T2_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT3_T3_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT3_B1_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT3_B2_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT3_B3_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT2_ss_data_triv(Nmasses_triv);
  vector<vector<data_t>> PT2_bs_data_triv(Nmasses_triv);
  for(int i=0;i<Nmasses_triv;i++) {
    PT3_T2_data_triv[i].resize(Nsou);
    PT3_T3_data_triv[i].resize(Nsou);
    PT3_B1_data_triv[i].resize(Nsou);
    PT3_B2_data_triv[i].resize(Nsou);
    PT3_B3_data_triv[i].resize(Nsou);
    PT2_ss_data_triv[i].resize(Nsou);
    PT2_bs_data_triv[i].resize(Nsou);
  }
 

  for(int i=0;i<Nmasses_triv;i++) {
    for(int j=0;j<Nsou;j++) {
      PT3_T2_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T2P5", Sort_light_confs);
      PT3_T3_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T3P5", Sort_light_confs);
      PT3_B1_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B1P5", Sort_light_confs);
      PT3_B2_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B2P5", Sort_light_confs);
      PT3_B3_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B3P5", Sort_light_confs);
      PT2_ss_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT2_SS_"+to_string(i)+"_SOU_"+to_string(j), "V1V1", Sort_light_confs);
      PT2_bs_data_triv[i][j].Read("../Bs_phi_gamma_data/trivial", "mes_contr_PT2_BS_"+to_string(i)+"_SOU_"+to_string(j), "P5P5", Sort_light_confs);
    }

  }


  vector<vector<data_t>> PT3_T2_data_traj_a(Nmasses);
  vector<vector<data_t>> PT3_T3_data_traj_a(Nmasses);
  vector<vector<data_t>> PT3_B1_data_traj_a(Nmasses);
  vector<vector<data_t>> PT3_B2_data_traj_a(Nmasses);
  vector<vector<data_t>> PT3_B3_data_traj_a(Nmasses);
  vector<vector<data_t>> PT2_ss_data_traj_a(Nmasses);
  vector<vector<data_t>> PT2_bs_data_traj_a(Nmasses);
  for(int i=0;i<Nmasses;i++) {
    PT3_T2_data_traj_a[i].resize(Nsou);
    PT3_T3_data_traj_a[i].resize(Nsou);
    PT3_B1_data_traj_a[i].resize(Nsou);
    PT3_B2_data_traj_a[i].resize(Nsou);
    PT3_B3_data_traj_a[i].resize(Nsou);
    PT2_ss_data_traj_a[i].resize(Nsou);
    PT2_bs_data_traj_a[i].resize(Nsou);
  }
 

  for(int i=0;i<Nmasses;i++) {
    for(int j=0;j<Nsou;j++) {
      PT3_T2_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T2P5", Sort_light_confs);
      PT3_T3_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T3P5", Sort_light_confs);
      PT3_B1_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B1P5", Sort_light_confs);
      PT3_B2_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B2P5", Sort_light_confs);
      PT3_B3_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B3P5", Sort_light_confs);
      PT2_ss_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT2_SS_"+to_string(i)+"_SOU_"+to_string(j), "V1V1", Sort_light_confs);
      PT2_bs_data_traj_a[i][j].Read("../Bs_phi_gamma_data/traj_a", "mes_contr_PT2_BS_"+to_string(i)+"_SOU_"+to_string(j), "P5P5", Sort_light_confs);
    }

  }


  vector<vector<data_t>> PT3_T2_data_traj_b(Nmasses);
  vector<vector<data_t>> PT3_T3_data_traj_b(Nmasses);
  vector<vector<data_t>> PT3_B1_data_traj_b(Nmasses);
  vector<vector<data_t>> PT3_B2_data_traj_b(Nmasses);
  vector<vector<data_t>> PT3_B3_data_traj_b(Nmasses);
  vector<vector<data_t>> PT2_ss_data_traj_b(Nmasses);
  vector<vector<data_t>> PT2_bs_data_traj_b(Nmasses);
  for(int i=0;i<Nmasses;i++) {
    PT3_T2_data_traj_b[i].resize(Nsou);
    PT3_T3_data_traj_b[i].resize(Nsou);
    PT3_B1_data_traj_b[i].resize(Nsou);
    PT3_B2_data_traj_b[i].resize(Nsou);
    PT3_B3_data_traj_b[i].resize(Nsou);
    PT2_ss_data_traj_b[i].resize(Nsou);
    PT2_bs_data_traj_b[i].resize(Nsou);
  }
 

  for(int i=0;i<Nmasses;i++) {
    for(int j=0;j<Nsou;j++) {
      PT3_T2_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T2P5", Sort_light_confs);
      PT3_T3_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T3P5", Sort_light_confs);
      PT3_B1_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B1P5", Sort_light_confs);
      PT3_B2_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B2P5", Sort_light_confs);
      PT3_B3_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B3P5", Sort_light_confs);
      PT2_ss_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT2_SS_"+to_string(i)+"_SOU_"+to_string(j), "V1V1", Sort_light_confs);
      PT2_bs_data_traj_b[i][j].Read("../Bs_phi_gamma_data/traj_b", "mes_contr_PT2_BS_"+to_string(i)+"_SOU_"+to_string(j), "P5P5", Sort_light_confs);
    }

  }


  vector<vector<data_t>> PT3_T2_data_rest(Nmasses);
  vector<vector<data_t>> PT3_T3_data_rest(Nmasses);
  vector<vector<data_t>> PT3_B1_data_rest(Nmasses);
  vector<vector<data_t>> PT3_B2_data_rest(Nmasses);
  vector<vector<data_t>> PT3_B3_data_rest(Nmasses);
  vector<vector<data_t>> PT2_ss_data_rest(Nmasses);
  vector<vector<data_t>> PT2_bs_data_rest(Nmasses);
  for(int i=0;i<Nmasses;i++) {
    PT3_T2_data_rest[i].resize(Nsou);
    PT3_T3_data_rest[i].resize(Nsou);
    PT3_B1_data_rest[i].resize(Nsou);
    PT3_B2_data_rest[i].resize(Nsou);
    PT3_B3_data_rest[i].resize(Nsou);
    PT2_ss_data_rest[i].resize(Nsou);
    PT2_bs_data_rest[i].resize(Nsou);
  }
 

  for(int i=0;i<Nmasses;i++) {
    for(int j=0;j<Nsou;j++) {
      PT3_T2_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T2P5", Sort_light_confs);
      PT3_T3_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "T3P5", Sort_light_confs);
      PT3_B1_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B1P5", Sort_light_confs);
      PT3_B2_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B2P5", Sort_light_confs);
      PT3_B3_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT3_MB"+to_string(i)+"_SOU_"+to_string(j), "B3P5", Sort_light_confs);
      PT2_ss_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT2_SS_0_SOU_"+to_string(j), "V1V1", Sort_light_confs);
      PT2_bs_data_rest[i][j].Read("../Bs_phi_gamma_data/rest", "mes_contr_PT2_BS_"+to_string(i)+"_SOU_"+to_string(j), "P5P5", Sort_light_confs);
    }

  }


  

  cout<<"Input read!"<<endl;
  
  int Nens= PT2_ss_data_triv[0][0].Tag.size();

  for(int iens=0; iens<Nens;iens++) {


    distr_t_list FF_T1_trivial_list(UseJack), FF_T1_traj_a_list(UseJack), FF_T1_traj_b_list(UseJack), FF_T1_rest_list(UseJack);
    distr_t_list FF_T2_trivial_list(UseJack), FF_T2_traj_a_list(UseJack), FF_T2_traj_b_list(UseJack), FF_T2_rest_list(UseJack);
    distr_t_list MBs_traj(UseJack), MBs_trivial(UseJack), MBs_rest(UseJack);

    //RCs
    distr_t Za, Zv, Z_T, a_distr;
    if(PT2_ss_data_triv[0][0].Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A; Z_T=ZT_A; a_distr=a_A;}
    else if(PT2_ss_data_triv[0][0].Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; Z_T=ZT_B; a_distr=a_B;}
    else if(PT2_ss_data_triv[0][0].Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; Z_T=ZT_C; a_distr=a_C;}
    else if(PT2_ss_data_triv[0][0].Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; Z_T=ZT_D; a_distr=a_D;}
    else crash("Ensemble: "+PT2_ss_data_triv[0][0].Tag[iens]+" not recognised");


     distr_t M_phi_rest(UseJack);
    //if(PT2_ss_data_traj_a[0][0].Tag[iens] == "cB211b.072.64") {
    // for(int ijack=0;ijack<Njacks;ijack++)  M_phi_rest.distr.push_back( 0.418 + 0.002*GM()/sqrt(Njacks-1.0)) ; 
    //}
    //else crash("Ensemble: "+PT2_ss_data_traj_a[0][0].Tag[iens]+" not yet implemented");
    

    Vfloat Thetas_triv, Thetas_traj_a, Thetas_traj_b, masses_b;
    vector<int> Tins_triv, Tins_traj_a, Tins_traj_b, Tins_rest;
    double mass_s;
    if(PT2_ss_data_triv[0][0].Tag[iens]=="cB211b.072.64") {
      Thetas_triv = { 1.734914, 2.464726 , 3.167159 };
      Thetas_traj_a = { 0.640817, 1.12876, 1.72705, 2.33141, 2.96339};
      Thetas_traj_b = { 1.05335, 1.66528, 2.33419, 2.95622, 3.56799};
      mass_s=0.0184;
      //Tins_traj_b= {30,30,30,30,30};
      //Tins_traj_a= {38,38,38,37,36};
      //Tins_triv= {38,37,36};
      //Tins_traj_b = {32,32,30,28,26};
      //Tins_traj_a = {32,32,30,28,26};
      //Tins_triv= {32,32,30};
      //Tins_rest= {32,32,32,32,32};
      //Tins_traj_b = {32,32,28,25,22};
      //Tins_traj_a = {32,32,28,25,22};
      //Tins_triv= {32,32,28};
      Tins_traj_b = {32,32,28,23,20};
      Tins_traj_a = {32,32,28,23,20};
      Tins_triv= {32,31,25};
      Tins_rest= {25,25,25,25,25};
      masses_b = {0.237, 0.35160, 0.48070, 0.60087, 0.72104};
    }
    else crash("Ensemble: "+PT2_ss_data_triv[0][0].Tag[iens]+" not found");

  
    

    boost::filesystem::create_directory("../data/Bs_phi_gamma");
    
    cout<<"Analyzing ensemble: "<<PT2_ss_data_triv[0][0].Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(PT2_ss_data_triv[0][0].Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = PT2_ss_data_triv[0][0].nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

   

    for(int im=0;im<Nmasses;im++) {
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im));
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]);
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_traj_a");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/masses_traj_a");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/FF_traj_a");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_traj_b");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/masses_traj_b");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/FF_traj_b");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_rest");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/masses_rest");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/FF_rest");

      
      Vfloat TT, TT_ins_a, TT_ins_b, TT_ins_rest;
      for(int t=0;t<Corr.Nt;t++) {TT_ins_a.push_back( Tins_traj_a[im]-t); };
      for(int t=0;t<Corr.Nt;t++) {TT_ins_b.push_back( Tins_traj_b[im]-t); };
      for(int t=0;t<Corr.Nt;t++) {TT_ins_rest.push_back( Tins_rest[im]-t); };

    
      VVfloat PT3_T2_a, PT3_T3_a,  PT3_B1_a, PT3_B2_a, PT3_B3_a, PT2_SS_a, PT2_BS_a;
      VVfloat PT3_T2_b, PT3_T3_b,  PT3_B1_b, PT3_B2_b, PT3_B3_b, PT2_SS_b, PT2_BS_b;
      VVfloat PT3_T2_rest, PT3_T3_rest,  PT3_B1_rest, PT3_B2_rest, PT3_B3_rest, PT2_SS_rest, PT2_BS_rest;
      for(int is=0;is<Nsou;is++) {

	if(is==0) {
	  PT3_T2_a = PT3_T2_data_traj_a[im][is].col(1)[iens];
	  PT3_T3_a = PT3_T3_data_traj_a[im][is].col(1)[iens];
	  PT3_B1_a= PT3_B1_data_traj_a[im][is].col(0)[iens];
	  PT3_B2_a= PT3_B2_data_traj_a[im][is].col(0)[iens];
	  PT3_B3_a= PT3_B3_data_traj_a[im][is].col(0)[iens];
	  PT2_SS_a=PT2_ss_data_traj_a[im][0].col(0)[iens];
	  PT2_BS_a=PT2_bs_data_traj_a[im][0].col(0)[iens];

	  PT3_T2_b = PT3_T2_data_traj_b[im][is].col(1)[iens];
	  PT3_T3_b = PT3_T3_data_traj_b[im][is].col(1)[iens];
	  PT3_B1_b= PT3_B1_data_traj_b[im][is].col(0)[iens];
	  PT3_B2_b= PT3_B2_data_traj_b[im][is].col(0)[iens];
	  PT3_B3_b= PT3_B3_data_traj_b[im][is].col(0)[iens];
	  PT2_SS_b=PT2_ss_data_traj_b[im][0].col(0)[iens];
	  PT2_BS_b=PT2_bs_data_traj_b[im][0].col(0)[iens];

	  PT3_T2_rest= PT3_T2_data_rest[im][is].col(1)[iens];
	  PT3_T3_rest= PT3_T3_data_rest[im][is].col(1)[iens];
	  PT3_B1_rest= PT3_B1_data_rest[im][is].col(0)[iens];
	  PT3_B2_rest= PT3_B2_data_rest[im][is].col(0)[iens];
	  PT3_B3_rest= PT3_B3_data_rest[im][is].col(0)[iens];
	  PT2_SS_rest= PT2_ss_data_rest[im][0].col(0)[iens];
	  PT2_BS_rest= PT2_bs_data_rest[im][0].col(0)[iens];


	}
	else {
	  PT3_T2_a = summ_master( PT3_T2_a, PT3_T2_data_traj_a[im][is].col(1)[iens]);
	  PT3_T3_a = summ_master( PT3_T3_a, PT3_T3_data_traj_a[im][is].col(1)[iens]);
	  PT3_B1_a = summ_master( PT3_B1_a, PT3_B1_data_traj_a[im][is].col(0)[iens]);
	  PT3_B2_a = summ_master( PT3_B2_a, PT3_B2_data_traj_a[im][is].col(0)[iens]);
	  PT3_B3_a = summ_master( PT3_B3_a, PT3_B3_data_traj_a[im][is].col(0)[iens]);
	  PT2_SS_a = summ_master( PT2_SS_a, PT2_ss_data_traj_a[im][is].col(0)[iens]);
	  PT2_BS_a = summ_master( PT2_BS_a, PT2_bs_data_traj_a[im][is].col(0)[iens]);

	  PT3_T2_b = summ_master( PT3_T2_b, PT3_T2_data_traj_b[im][is].col(1)[iens]);
	  PT3_T3_b = summ_master( PT3_T3_b, PT3_T3_data_traj_b[im][is].col(1)[iens]);
	  PT3_B1_b = summ_master( PT3_B1_b, PT3_B1_data_traj_b[im][is].col(0)[iens]);
	  PT3_B2_b = summ_master( PT3_B2_b, PT3_B2_data_traj_b[im][is].col(0)[iens]);
	  PT3_B3_b = summ_master( PT3_B3_b, PT3_B3_data_traj_b[im][is].col(0)[iens]);
	  PT2_SS_b = summ_master( PT2_SS_b, PT2_ss_data_traj_b[im][is].col(0)[iens]);
	  PT2_BS_b = summ_master( PT2_BS_b, PT2_bs_data_traj_b[im][is].col(0)[iens]);

	  PT3_T2_rest = summ_master( PT3_T2_rest, PT3_T2_data_rest[im][is].col(1)[iens]);
	  PT3_T3_rest = summ_master( PT3_T3_rest, PT3_T3_data_rest[im][is].col(1)[iens]);
	  PT3_B1_rest = summ_master( PT3_B1_rest, PT3_B1_data_rest[im][is].col(0)[iens]);
	  PT3_B2_rest = summ_master( PT3_B2_rest, PT3_B2_data_rest[im][is].col(0)[iens]);
	  PT3_B3_rest = summ_master( PT3_B3_rest, PT3_B3_data_rest[im][is].col(0)[iens]);
	  PT2_SS_rest = summ_master( PT2_SS_rest, PT2_ss_data_rest[im][is].col(0)[iens]);
	  PT2_BS_rest = summ_master( PT2_BS_rest, PT2_bs_data_rest[im][is].col(0)[iens]);

	  
	}
	
      }
      PT3_T2_a= Multiply_Vvector_by_scalar(PT3_T2_a, 1.0/Nsou);
      PT3_T3_a= Multiply_Vvector_by_scalar(PT3_T3_a, 1.0/Nsou);
      PT3_B1_a= Multiply_Vvector_by_scalar(PT3_B1_a, 1.0/Nsou);
      PT3_B2_a= Multiply_Vvector_by_scalar(PT3_B2_a, 1.0/Nsou);
      PT3_B3_a= Multiply_Vvector_by_scalar(PT3_B3_a, 1.0/Nsou);
      PT2_SS_a= Multiply_Vvector_by_scalar(PT2_SS_a, 1.0/Nsou);
      PT2_BS_a= Multiply_Vvector_by_scalar(PT2_BS_a, 1.0/Nsou);

      PT3_T2_b= Multiply_Vvector_by_scalar(PT3_T2_b, 1.0/Nsou);
      PT3_T3_b= Multiply_Vvector_by_scalar(PT3_T3_b, 1.0/Nsou);
      PT3_B1_b= Multiply_Vvector_by_scalar(PT3_B1_b, 1.0/Nsou);
      PT3_B2_b= Multiply_Vvector_by_scalar(PT3_B2_b, 1.0/Nsou);
      PT3_B3_b= Multiply_Vvector_by_scalar(PT3_B3_b, 1.0/Nsou);
      PT2_SS_b= Multiply_Vvector_by_scalar(PT2_SS_b, 1.0/Nsou);
      PT2_BS_b= Multiply_Vvector_by_scalar(PT2_BS_b, 1.0/Nsou);

      PT3_T2_rest= Multiply_Vvector_by_scalar(PT3_T2_rest, 1.0/Nsou);
      PT3_T3_rest= Multiply_Vvector_by_scalar(PT3_T3_rest, 1.0/Nsou);
      PT3_B1_rest= Multiply_Vvector_by_scalar(PT3_B1_rest, 1.0/Nsou);
      PT3_B2_rest= Multiply_Vvector_by_scalar(PT3_B2_rest, 1.0/Nsou);
      PT3_B3_rest= Multiply_Vvector_by_scalar(PT3_B3_rest, 1.0/Nsou);
      PT2_SS_rest= Multiply_Vvector_by_scalar(PT2_SS_rest, 1.0/Nsou);
      PT2_BS_rest= Multiply_Vvector_by_scalar(PT2_BS_rest, 1.0/Nsou);


      //analyze correlators
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_T2_a_distr = Corr.corr_t(PT3_T2_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_T2");
      distr_t_list PT3_T3_a_distr = Corr.corr_t(PT3_T3_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_T3");
      distr_t_list PT3_B1_a_distr = Corr.corr_t(PT3_B1_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_B1");
      distr_t_list PT3_B2_a_distr = Corr.corr_t(PT3_B2_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_B2");
      distr_t_list PT3_B3_a_distr = Corr.corr_t(PT3_B3_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_B3");
      Corr.Perform_Nt_t_average=1;
      distr_t_list PT2_SS_a_distr = Corr.corr_t(PT2_SS_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT2_SS");
      distr_t_list PT2_BS_a_distr = Corr.corr_t(PT2_BS_a, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT2_BS");

      Print_To_File( {}, { (PT3_B1_a_distr  - 0.5*(PT3_B2_a_distr + PT3_B3_a_distr)).ave(), (PT3_B1_a_distr  - 0.5*(PT3_B2_a_distr + PT3_B3_a_distr)).err() } ,  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_a/PT3_BB.t", "", "");
      
      

      distr_t_list Bs_eff_mass_a= Corr.effective_mass_t( PT2_BS_a_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_a[0][0].Tag[iens]+"/masses_traj_a/Bs_mass");
      distr_t_list phi_eff_mass_a= Corr.effective_mass_t( PT2_SS_a_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/masses_traj_a/phi_mass");

      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_T2_b_distr = Corr.corr_t(PT3_T2_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT3_T2");
      distr_t_list PT3_T3_b_distr = Corr.corr_t(PT3_T3_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT3_T3");
      distr_t_list PT3_B1_b_distr = Corr.corr_t(PT3_B1_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT3_B1");
      distr_t_list PT3_B2_b_distr = Corr.corr_t(PT3_B2_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT3_B2");
      distr_t_list PT3_B3_b_distr = Corr.corr_t(PT3_B3_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT3_B3");
      Corr.Perform_Nt_t_average=1;
      distr_t_list PT2_SS_b_distr = Corr.corr_t(PT2_SS_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT2_SS");
      distr_t_list PT2_BS_b_distr = Corr.corr_t(PT2_BS_b, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/corr_traj_b/PT2_BS");

      Print_To_File( {}, { (PT3_B1_b_distr  - 0.5*(PT3_B2_b_distr + PT3_B3_b_distr)).ave(), (PT3_B1_b_distr  - 0.5*(PT3_B2_b_distr + PT3_B3_b_distr)).err() } ,  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_traj_b/PT3_BB.t", "", "");

      distr_t_list Bs_eff_mass_b= Corr.effective_mass_t( PT2_BS_b_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_b[0][0].Tag[iens]+"/masses_traj_b/Bs_mass");
      distr_t_list phi_eff_mass_b= Corr.effective_mass_t( PT2_SS_b_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_b[0][0].Tag[iens]+"/masses_traj_b/phi_mass");

      //######
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_T2_rest_distr = Corr.corr_t(PT3_T2_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_T2");
      distr_t_list PT3_T3_rest_distr = Corr.corr_t(PT3_T3_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_T3");
      distr_t_list PT3_B1_rest_distr = Corr.corr_t(PT3_B1_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_B1");
      distr_t_list PT3_B2_rest_distr = Corr.corr_t(PT3_B2_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_B2");
      distr_t_list PT3_B3_rest_distr = Corr.corr_t(PT3_B3_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_B3");
      Corr.Perform_Nt_t_average=1;
      distr_t_list PT2_SS_rest_distr = Corr.corr_t(PT2_SS_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT2_SS");
      distr_t_list PT2_BS_rest_distr = Corr.corr_t(PT2_BS_rest, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT2_BS");

      Print_To_File( {}, { (PT3_B1_rest_distr  - 0.5*(PT3_B2_rest_distr + PT3_B3_rest_distr)).ave(), (PT3_B1_rest_distr  - 0.5*(PT3_B2_rest_distr + PT3_B3_rest_distr)).err() } ,  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/corr_rest/PT3_BB.t", "", "");
      
      

      distr_t_list Bs_eff_mass_rest= Corr.effective_mass_t( PT2_BS_rest_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_a[0][0].Tag[iens]+"/masses_rest/Bs_mass");
      distr_t_list phi_eff_mass_rest= Corr.effective_mass_t( PT2_SS_rest_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_traj_a[0][0].Tag[iens]+"/masses_rest/phi_mass");


      //######
      

      int Tmin_Bs=0, Tmax_Bs=0;
      int Tmin_phi=0, Tmax_phi=0;
      int Tmin_Bs_rest=0, Tmax_Bs_rest=0;
      
      if(PT2_ss_data_traj_a[0][0].Tag[iens] =="cB211b.072.64") {
	if(im == 0) {
	  Tmin_Bs=19;  Tmax_Bs=30;
	  Tmin_Bs_rest= 16; Tmax_Bs_rest= 27;
	  Tmin_phi=14; Tmax_phi=20;
	}
	else if(im == 1) {
	  Tmin_Bs=17;  Tmax_Bs=24;
	  Tmin_Bs_rest= 15 ; Tmax_Bs_rest= 22 ;
	  Tmin_phi=14; Tmax_phi=18;
	}
	else if(im==2) {
	  Tmin_Bs=17;  Tmax_Bs=23;
	  Tmin_Bs_rest= 14 ; Tmax_Bs_rest= 20 ;
	  Tmin_phi=14; Tmax_phi=17;
	}
	else if(im==3) {
	  Tmin_Bs=17;  Tmax_Bs=23;
	  Tmin_Bs_rest= 14 ; Tmax_Bs_rest= 20 ;
	  Tmin_phi=13; Tmax_phi=16;
	}
	else if(im==4) {
	  Tmin_Bs=12;  Tmax_Bs=16;
	  Tmin_Bs_rest= 9 ; Tmax_Bs_rest= 15 ;
	  Tmin_phi=11; Tmax_phi=13;
	}
      }
      else crash("Ensemble: "+PT2_ss_data_traj_a[0][0].Tag[iens]+" not found");
      
      Corr.Tmin=Tmin_Bs; Corr.Tmax=Tmax_Bs;
      distr_t M_Bs = Corr.Fit_distr(Bs_eff_mass_a);
      Corr.Tmin=Tmin_Bs_rest; Corr.Tmax= Tmax_Bs_rest;
      distr_t M_Bs_rest= Corr.Fit_distr(Bs_eff_mass_rest);
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t M_phi_b = Corr.Fit_distr(phi_eff_mass_b);
      distr_t M_phi_a = Corr.Fit_distr(phi_eff_mass_a);
      M_phi_rest = Corr.Fit_distr(phi_eff_mass_rest); 
      MBs_traj.distr_list.push_back( M_Bs/a_distr);
      MBs_rest.distr_list.push_back( M_Bs_rest/a_distr);

      

      //get photon momentum [traj_a]
      pt3_momenta pt3_mom_a(Vfloat({0.0,0.0,0.0}), Vfloat({0.0, 0.0, 0.0}), Vfloat({Thetas_traj_a[im], Thetas_traj_a[im], Thetas_traj_a[im]}), masses_b[im], mass_s, 0.0, L_info.L, L_info.T);
      double Eg_a= pt3_mom_a.Egamma();
      distr_t Eg_a_virt= M_Bs - M_phi_a;
      distr_t xg_a= pt3_mom_a.x_gamma(M_Bs);
      double kz_a = pt3_mom_a.k()[2];
      cout<<"im: "<<im<<" xg[traj_a]: "<<xg_a.ave()<<endl;
      cout<<"kz[traj_a]: "<<kz_a<<endl;
      //get photon momentum [traj_b]
      pt3_momenta pt3_mom_b(Vfloat({0.0,0.0,0.0}), Vfloat({0.0, 0.0, 0.0}), Vfloat({Thetas_traj_b[im], Thetas_traj_b[im], Thetas_traj_b[im]}), masses_b[im], mass_s, 0.0, L_info.L, L_info.T);
      double Eg_b= pt3_mom_b.Egamma();
      distr_t Eg_b_virt = M_Bs - M_phi_b;
      distr_t xg_b= pt3_mom_b.x_gamma(M_Bs);
      double kz_b = pt3_mom_b.k()[2];
      cout<<"im: "<<im<<" xg[traj_b]: "<<xg_b.ave()<<endl;
      cout<<"kz[traj_b]: "<<kz_b<<endl;
     
      Corr.Tmin=Tmin_Bs; Corr.Tmax=Tmax_Bs;
      distr_t F_B_a= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_BS_a_distr, ""))/2.0;
      distr_t F_B_b= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_BS_b_distr, ""))/2.0;
      distr_t F_B_rest= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_BS_rest_distr, ""))/2.0;
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t F_S_a= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_SS_a_distr, ""))/2.0;
      distr_t F_S_b= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_SS_b_distr, ""))/2.0;
      distr_t F_S_rest= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_SS_rest_distr, ""))/2.0;

     
      distr_t_list Z_FACTOR_a= (Z_T/(F_B_a*F_S_a))*EXPT_D(M_Bs, Corr.Nt)*EXPT_D(-1.0*M_phi_a, Corr.Nt)*EXP_D(M_phi_a*Tins_traj_a[im]);
      distr_t_list Z_FACTOR_b= (Z_T/(F_B_b*F_S_b))*EXPT_D(M_Bs, Corr.Nt)*EXPT_D(-1.0*M_phi_b, Corr.Nt)*EXP_D(M_phi_b*Tins_traj_b[im]);
      distr_t_list Z_FACTOR_rest= (Z_T/(F_B_rest*F_S_rest))*EXPT_D(M_Bs_rest, Corr.Nt)*EXPT_D(-1.0*M_phi_rest, Corr.Nt)*EXP_D(M_phi_rest*Tins_rest[im]);

     
      distr_t GAMMA_a= SQRT_D( 1.0 + kz_a*kz_a/(M_phi_rest*M_phi_rest));
      distr_t GAMMA_b= SQRT_D( 1.0 + kz_b*kz_b/(M_phi_rest*M_phi_rest));

           
      
      distr_t_list T1_a(UseJack), T2_a(UseJack);
      
      //for(int t=0;t<Tins_traj_a[im];t++) FF_a.distr_list.push_back( Z_FACTOR_a.distr_list[t]*(PT3_T2_a_distr.distr_list[t] + PT3_B1_a_distr.distr_list[t])*xg_a/(4.0*kz_a));
      for(int t=0;t<Tins_traj_a[im];t++) {  T1_a.distr_list.push_back(  GAMMA_a*Z_FACTOR_a.distr_list[t]*(-1.0*( sign_kz*(Eg_a_virt/(4.0*kz_a*M_Bs))*(PT3_T3_a_distr.distr_list[t] - PT3_T2_a_distr.distr_list[t]) + (1.0/(4.0*M_Bs))*( 2.0*PT3_B1_a_distr.distr_list[t] - PT3_B2_a_distr.distr_list[t] - PT3_B3_a_distr.distr_list[t])))); }  //a factor of 2 in the denominator comes from the presence of 2T1(q^2) in the FF decomposition of <V | J_T | B> matrix element

      for(int t=0;t<Tins_traj_a[im];t++) {  T2_a.distr_list.push_back(  GAMMA_a*Z_FACTOR_a.distr_list[t]*(-1.0*( sign_kz*(3.0*kz_a/(2.0*( M_Bs*M_Bs- M_phi_rest*M_phi_rest)  ))*(PT3_T3_a_distr.distr_list[t] - PT3_T2_a_distr.distr_list[t]) + (Eg_a_virt/(2.0*(M_Bs*M_Bs-M_phi_rest*M_phi_rest)))*( 2.0*PT3_B1_a_distr.distr_list[t] - PT3_B2_a_distr.distr_list[t] - PT3_B3_a_distr.distr_list[t])))); } 
     
      Print_To_File({}, { T1_a.ave(), T1_a.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_a[0][0].Tag[iens]+"/FF_traj_a/T1", "", "");
      Print_To_File({}, { T2_a.ave(), T2_a.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_a[0][0].Tag[iens]+"/FF_traj_a/T2", "", "");

     
      
      distr_t_list T1_b(UseJack), T2_b(UseJack);
      
      for(int t=0;t<Tins_traj_b[im];t++) {  T1_b.distr_list.push_back(  GAMMA_b*Z_FACTOR_b.distr_list[t]*(-1.0*( sign_kz*(Eg_b_virt/(4.0*kz_b*M_Bs))*(PT3_T3_b_distr.distr_list[t] - PT3_T2_b_distr.distr_list[t]) + (1.0/(4.0*M_Bs))*( 2.0*PT3_B1_b_distr.distr_list[t] - PT3_B2_b_distr.distr_list[t] - PT3_B3_b_distr.distr_list[t])))); }
      
      for(int t=0;t<Tins_traj_b[im];t++) {  T2_b.distr_list.push_back(  GAMMA_b*Z_FACTOR_b.distr_list[t]*(-1.0*( sign_kz*(3.0*kz_b/(2.0*( M_Bs*M_Bs- M_phi_rest*M_phi_rest)  ))*(PT3_T3_b_distr.distr_list[t] - PT3_T2_b_distr.distr_list[t]) + (Eg_b_virt/(2.0*(M_Bs*M_Bs-M_phi_rest*M_phi_rest)))*( 2.0*PT3_B1_b_distr.distr_list[t] - PT3_B2_b_distr.distr_list[t] - PT3_B3_b_distr.distr_list[t])))); } 
     
      Print_To_File({}, { T1_b.ave(), T1_b.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_b[0][0].Tag[iens]+"/FF_traj_b/T1", "", "");
      Print_To_File({}, { T2_b.ave(), T2_b.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_traj_b[0][0].Tag[iens]+"/FF_traj_b/T2", "", "");

     

     
      distr_t_list T2_rest(UseJack);
      for(int t=0;t<Tins_rest[im];t++) {  T2_rest.distr_list.push_back(  Z_FACTOR_rest.distr_list[t]*(-1.0*((1.0/((M_Bs_rest + M_phi_rest)))*(PT3_B1_rest_distr.distr_list[t] ))));}

      Print_To_File({}, { T2_rest.ave(), T2_rest.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_rest[0][0].Tag[iens]+"/FF_rest/T2", "", "");

           
     

      //traj a T1
      if(im==0) { Corr.Tmin=20; Corr.Tmax=25 ; }
      else if(im==1) { Corr.Tmin=17; Corr.Tmax=24;}
      else if(im==2) { Corr.Tmin=15; Corr.Tmax=23;}
      else if(im==3) { Corr.Tmin=11; Corr.Tmax=16;}
      else if(im==4) { Corr.Tmin=12; Corr.Tmax=17;}
      else crash("im="+to_string(im)+" is invalid");

      FF_T1_traj_a_list.distr_list.push_back( Corr.Fit_distr( T1_a));

      //traj a T2
      if(im==0) { Corr.Tmin=20; Corr.Tmax=25 ; }
      else if(im==1) { Corr.Tmin=18; Corr.Tmax=25;}
      else if(im==2) { Corr.Tmin=14; Corr.Tmax=22;}
      else if(im==3) { Corr.Tmin=11; Corr.Tmax=16;}
      else if(im==4) { Corr.Tmin=13; Corr.Tmax=17;}
      else crash("im="+to_string(im)+" is invalid");

      
      FF_T2_traj_a_list.distr_list.push_back( Corr.Fit_distr( T2_a));

      //traj b T1
      if(im==0) { Corr.Tmin=20; Corr.Tmax=25 ; }
      else if(im==1) { Corr.Tmin=18; Corr.Tmax=26;}
      else if(im==2) { Corr.Tmin=19; Corr.Tmax=23;}
      else if(im==3) { Corr.Tmin=14; Corr.Tmax=20;}
      else if(im==4) { Corr.Tmin=13; Corr.Tmax=16;}
      else crash("im="+to_string(im)+" is invalid");
   
      FF_T1_traj_b_list.distr_list.push_back( Corr.Fit_distr( T1_b));

      //traj b T2
      if(im==0) { Corr.Tmin=19; Corr.Tmax=24 ; }
      else if(im==1) { Corr.Tmin=18; Corr.Tmax=23;}
      else if(im==2) { Corr.Tmin=19; Corr.Tmax=23;}
      else if(im==3) { Corr.Tmin=14; Corr.Tmax=20;}
      else if(im==4) { Corr.Tmin=14; Corr.Tmax=19;}
      else crash("im="+to_string(im)+" is invalid");
     
      
      FF_T2_traj_b_list.distr_list.push_back( Corr.Fit_distr( T2_b));


      //rest T2
      if(im==0) { Corr.Tmin=20; Corr.Tmax=23 ; }
      else if(im==1) { Corr.Tmin=19; Corr.Tmax=23;}
      else if(im==2) { Corr.Tmin=19; Corr.Tmax=23;}
      else if(im==3) { Corr.Tmin=17; Corr.Tmax=21;}
      else if(im==4) { Corr.Tmin=16; Corr.Tmax=19;}
      else crash("im="+to_string(im)+" is invalid");
      
            
      FF_T2_rest_list.distr_list.push_back( Corr.Fit_distr( T2_rest));


     
    }


    //######################################## trivial traj #########################################à

     for(int im=0;im<Nmasses_triv;im++) {
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im));
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]);
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/masses_triv");
      boost::filesystem::create_directory("../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/FF_triv");
    
      Vfloat TT, TT_ins;
      for(int t=0;t<Corr.Nt;t++) {TT_ins.push_back( Tins_triv[im]-t); };
          
      VVfloat PT3_T2, PT3_T3,  PT3_B1, PT3_B2, PT3_B3, PT2_SS, PT2_BS;
      
      for(int is=0;is<Nsou;is++) {

	if(is==0) {
	  PT3_T2 = PT3_T2_data_triv[im][is].col(1)[iens];
	  PT3_T3 = PT3_T3_data_triv[im][is].col(1)[iens];
	  PT3_B1= PT3_B1_data_triv[im][is].col(0)[iens];
	  PT3_B2= PT3_B2_data_triv[im][is].col(0)[iens];
	  PT3_B3= PT3_B3_data_triv[im][is].col(0)[iens];
	  PT2_SS=PT2_ss_data_triv[im][0].col(0)[iens];
	  PT2_BS=PT2_bs_data_triv[im][0].col(0)[iens];



	}
	else {
	  PT3_T2 = summ_master( PT3_T2, PT3_T2_data_triv[im][is].col(1)[iens]);
	  PT3_T3 = summ_master( PT3_T3, PT3_T3_data_triv[im][is].col(1)[iens]);
	  PT3_B1 = summ_master( PT3_B1, PT3_B1_data_triv[im][is].col(0)[iens]);
	  PT3_B2 = summ_master( PT3_B2, PT3_B2_data_triv[im][is].col(0)[iens]);
	  PT3_B3 = summ_master( PT3_B3, PT3_B3_data_triv[im][is].col(0)[iens]);
	  PT2_SS = summ_master( PT2_SS, PT2_ss_data_triv[im][is].col(0)[iens]);
	  PT2_BS = summ_master( PT2_BS, PT2_bs_data_triv[im][is].col(0)[iens]);

	 
	}
	
      }
      
      PT3_T2= Multiply_Vvector_by_scalar(PT3_T2, 1.0/Nsou);
      PT3_T3= Multiply_Vvector_by_scalar(PT3_T3, 1.0/Nsou);
      PT3_B1= Multiply_Vvector_by_scalar(PT3_B1, 1.0/Nsou);
      PT3_B2= Multiply_Vvector_by_scalar(PT3_B2, 1.0/Nsou);
      PT3_B3= Multiply_Vvector_by_scalar(PT3_B3, 1.0/Nsou);
      PT2_SS= Multiply_Vvector_by_scalar(PT2_SS, 1.0/Nsou);
      PT2_BS= Multiply_Vvector_by_scalar(PT2_BS, 1.0/Nsou);

      

      //analyze correlators
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_T2_distr = Corr.corr_t(PT3_T2, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT3_T2");
      distr_t_list PT3_T3_distr = Corr.corr_t(PT3_T3, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT3_T3");
      distr_t_list PT3_B1_distr = Corr.corr_t(PT3_B1, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT3_B1");
      distr_t_list PT3_B2_distr = Corr.corr_t(PT3_B2, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT3_B2");
      distr_t_list PT3_B3_distr = Corr.corr_t(PT3_B3, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT3_B3");
      Corr.Perform_Nt_t_average=1;
      distr_t_list PT2_SS_distr = Corr.corr_t(PT2_SS, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT2_SS");
      distr_t_list PT2_BS_distr = Corr.corr_t(PT2_BS, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/corr_triv/PT2_BS");

      distr_t_list Bs_eff_mass= Corr.effective_mass_t( PT2_BS_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_triv[0][0].Tag[iens]+"/masses_triv/Bs_mass");
      distr_t_list phi_eff_mass= Corr.effective_mass_t( PT2_SS_distr, "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_ss_data_triv[0][0].Tag[iens]+"/masses_triv/phi_mass");

      
      int Tmin_Bs=0, Tmax_Bs=0;
      int Tmin_phi=0, Tmax_phi=0;
      
      if(PT2_ss_data_traj_a[0][0].Tag[iens] =="cB211b.072.64") {
	if(im == 0) {
	  Tmin_Bs=19;  Tmax_Bs=30;
	  Tmin_phi=14; Tmax_phi=20;
	}
	else if(im == 1) {
	  Tmin_Bs=17;  Tmax_Bs=24;
	  Tmin_phi=14; Tmax_phi=18;
	}
	else if(im==2) {
	  Tmin_Bs=17;  Tmax_Bs=23;
	  Tmin_phi=14; Tmax_phi=17;
	}
	else if(im==3) {
	  Tmin_Bs=17;  Tmax_Bs=23;
	  Tmin_phi=13; Tmax_phi=16;
	}
	else if(im==4) {
	  Tmin_Bs=12;  Tmax_Bs=16;
	  Tmin_phi=10; Tmax_phi=11;
	}
      }
      else crash("Ensemble: "+PT2_ss_data_traj_a[0][0].Tag[iens]+" not found");
      
      Corr.Tmin=Tmin_Bs; Corr.Tmax=Tmax_Bs;
      distr_t M_Bs = Corr.Fit_distr(Bs_eff_mass);
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t M_phi = Corr.Fit_distr(phi_eff_mass);
      


      MBs_trivial.distr_list.push_back( M_Bs/a_distr);

      //get photon momentum [traj_a]
      pt3_momenta pt3_mom( Vfloat({0.0,0.0,0.0}), Vfloat({0.0, 0.0, 0.0}), Vfloat({Thetas_triv[im], Thetas_triv[im], Thetas_triv[im]}) , masses_b[im], mass_s, 0.0, L_info.L, L_info.T);
      pt3_momenta pt3_mom_test( Vfloat({0.0,0.0,0.0}), Vfloat({0.0, 0.0, 0.0}), Vfloat({0.0, 0.0, sqrt(3.0)*Thetas_triv[im]}) , masses_b[im], mass_s, 0.0, L_info.L, L_info.T);
      double Eg= pt3_mom.Egamma();
      distr_t Eg_virt = M_Bs - M_phi;
      distr_t xg= pt3_mom.x_gamma(M_Bs);
      double kz = pt3_mom.k()[2];
      double kz_test= pt3_mom_test.k()[2];
      cout<<"im: "<<im<<" xg[triv]: "<<xg.ave()<<endl;
      cout<<"kz[triv]: "<<kz<<", k_test: "<<kz_test/sqrt(3.0)<<" expected: "<<((M_Bs*M_Bs-M_phi_rest*M_phi_rest)/(2.0*M_Bs*sqrt(3.0))).ave()<<endl;           
      Corr.Tmin=Tmin_Bs; Corr.Tmax=Tmax_Bs;
      distr_t F_B= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_BS_distr, ""))/2.0;
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t F_S= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_SS_distr, ""))/2.0;

      distr_t_list Z_FACTOR= (Z_T/(F_B*F_S))*EXPT_D(M_Bs, Corr.Nt)*EXPT_D(-1.0*M_phi, Corr.Nt)*EXP_D(M_phi*Tins_triv[im]);




      distr_t GAMMA= SQRT_D( 1.0 + kz*kz/(M_phi_rest*M_phi_rest));







      
      
      
      
      distr_t_list T1(UseJack), T2(UseJack);

      for(int t=0;t<Tins_triv[im];t++) {  T1.distr_list.push_back(  GAMMA*Z_FACTOR.distr_list[t]*(-1.0*( sign_kz*(Eg/(4.0*kz*M_Bs))*(PT3_T3_distr.distr_list[t] - PT3_T2_distr.distr_list[t]) + (1.0/(4.0*M_Bs))*( 2.0*PT3_B1_distr.distr_list[t] - PT3_B2_distr.distr_list[t] - PT3_B3_distr.distr_list[t])))); }

      for(int t=0;t<Tins_triv[im];t++) {  T2.distr_list.push_back(  GAMMA*Z_FACTOR.distr_list[t]*(-1.0*( sign_kz*(3.0*kz/(2.0*( M_Bs*M_Bs- M_phi_rest*M_phi_rest)  ))*(PT3_T3_distr.distr_list[t] - PT3_T2_distr.distr_list[t]) + (Eg/(2.0*(M_Bs*M_Bs-M_phi_rest*M_phi_rest)))*( 2.0*PT3_B1_distr.distr_list[t] - PT3_B2_distr.distr_list[t] - PT3_B3_distr.distr_list[t])))); } 
      
     
      Print_To_File({}, { T1.ave(), T1.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_triv[0][0].Tag[iens]+"/FF_triv/T1", "", "");
      Print_To_File({}, { T2.ave(), T2.err()},  "../data/Bs_phi_gamma/MB_"+to_string(im)+"/"+PT2_bs_data_triv[0][0].Tag[iens]+"/FF_triv/T2", "", "");
      
      //FF trivial T1
      if(im==0) { Corr.Tmin=18; Corr.Tmax=21 ; }
      else if(im==1) { Corr.Tmin=19; Corr.Tmax=23;}
      else if(im==2) { Corr.Tmin=17; Corr.Tmax=19;}
      else crash("im="+to_string(im)+" is invalid");

      FF_T1_trivial_list.distr_list.push_back( Corr.Fit_distr( T1));      
      FF_T2_trivial_list.distr_list.push_back( Corr.Fit_distr( T2));


     }

     Print_To_File( { }, {MBs_trivial.ave(), FF_T1_trivial_list.ave(), FF_T1_trivial_list.err(), FF_T2_trivial_list.ave(), FF_T2_trivial_list.err() }, "../data/Bs_phi_gamma/B64_trivial", "", "");
     Print_To_File( { }, {MBs_traj.ave(), FF_T1_traj_a_list.ave(), FF_T1_traj_a_list.err(), FF_T2_traj_a_list.ave(), FF_T2_traj_a_list.err() }, "../data/Bs_phi_gamma/B64_traj_a", "", "");
     Print_To_File( { }, {MBs_traj.ave(), FF_T1_traj_b_list.ave(), FF_T1_traj_b_list.err(), FF_T2_traj_b_list.ave(), FF_T2_traj_b_list.err() }, "../data/Bs_phi_gamma/B64_traj_b", "", "");
    
  }


  return;
}
