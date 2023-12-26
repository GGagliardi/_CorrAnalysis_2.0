#include "../include/HVP.h"


const double DTT = 0.5;
const double alpha = 1.0/137.035999;	     

using namespace std;

void Bounding_HVP(distr_t &amu_HVP, int &Tcut_opt,  const distr_t_list &V, const distr_t &a, string path,distr_t lowest_mass) {


  int Njacks=50;
  bool UseJack=true;
 
  int Simps_ord=3;
  double fm_to_inv_Gev= 1.0/0.197327;


  double T= V.size(); 

  auto LOG = [](double R_G, double t) { return log(fabs(R_G));};

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a , T/2);

  int T_ext_max= 300;

  distr_t_list Ker_extended= distr_t_list::f_of_distr(K, a , T_ext_max);

  int Tdatas_opt=-1;

  
  auto amu_HVP_func = [&V, &Ker, &UseJack, &Njacks,  &Simps_ord](double tcut) -> distr_t {
    
    distr_t ret_win(UseJack,UseJack?Njacks:800);
    for(int t=1;t<tcut;t++) ret_win = ret_win + 4.0*w(t,Simps_ord)*pow(alpha,2)*V.distr_list[t]*Ker.distr_list[t];
    return ret_win;
    
  };

  auto amu_HVP_int = [ &Ker_extended, &Simps_ord](double t) -> distr_t {

    return 4.0*w(t,Simps_ord)*pow(alpha,2)*Ker_extended.distr_list[t];

  };
  
 
 
  CorrAnalysis Corr(UseJack, Njacks,800);
  Corr.Nt = V.size();
  distr_t_list ratio_corr_V(UseJack);
  for(int t=0; t<Corr.Nt;t++) ratio_corr_V.distr_list.push_back( V.distr_list[t]/V.distr_list[(t+1)%Corr.Nt]);
  //decide whether to use m_eff(t) or log( V(t)/V(t+1)) (true = use log)  
  bool update_min_Tdata_from_log_ratio=true;    
  distr_t_list eff_mass_V = (update_min_Tdata_from_log_ratio)?distr_t_list::f_of_distr_list(LOG, ratio_corr_V):Corr.effective_mass_t(V, "");

  

  distr_t_list amu_HVP_min_Tdata(UseJack), amu_HVP_max_Tdata(UseJack), amu_HVP_T_2(UseJack);

  Vfloat TCUTS;

  Vfloat Is_T_data_opt;
  int slice_to_use_for_eff_mass=1;
    
    
  //loop over tcut
  for(int tcut=1; tcut<T/2;tcut++) {


    TCUTS.push_back( (double)tcut);
    distr_t amu_HVP_up_to_tcut = amu_HVP_func(tcut+1);
    amu_HVP_min_Tdata.distr_list.push_back(amu_HVP_up_to_tcut);
    amu_HVP_max_Tdata.distr_list.push_back(amu_HVP_up_to_tcut);
    amu_HVP_T_2.distr_list.push_back(amu_HVP_up_to_tcut);
    distr_t V_tcut = V.distr_list[tcut];
    Is_T_data_opt.push_back( 0.0);
      
    bool eff_mass_is_nan= isnan( eff_mass_V.ave(tcut));
    bool update_min_Tdata=true;
    if(eff_mass_is_nan || (eff_mass_V.err(tcut)/eff_mass_V.ave(tcut) > 0.05) || (eff_mass_V.ave(tcut) < 0 )) update_min_Tdata=false;

    if(update_min_Tdata) slice_to_use_for_eff_mass=tcut;
      
      

    for(int t=tcut+1; t < T_ext_max;t++) {

      //lambda function for lower and upper limit of single exp V(t)
      auto EXP_MIN = [&tcut, &t](double E) { return exp(-E*(t-tcut));};
	
      distr_t ker_val = V_tcut*amu_HVP_int(t);
      int size_min= amu_HVP_min_Tdata.size();
      int size_max= amu_HVP_max_Tdata.size();
      if(size_min != tcut || size_max != tcut) crash("size_min or size_max is different from tcut");
      distr_t lower_exp = distr_t::f_of_distr(EXP_MIN,eff_mass_V.distr_list[slice_to_use_for_eff_mass]);
      amu_HVP_min_Tdata.distr_list[size_min-1] = amu_HVP_min_Tdata.distr_list[size_min-1] + ker_val*lower_exp;
      amu_HVP_max_Tdata.distr_list[size_max-1] = amu_HVP_max_Tdata.distr_list[size_max-1] + ker_val*distr_t::f_of_distr(EXP_MIN, lowest_mass);
    }

      
  }


   


  //find interval where difference between amu_HVP_min_Tdata and amu_HVP_max_Tdata is smaller than 0.3 sigma
  //average
  bool Found_Tdata_opt=false;
  int tdata_opt=1;

   
  while(!Found_Tdata_opt && tdata_opt < T/2) {

    distr_t diff_max_min = amu_HVP_max_Tdata.distr_list[tdata_opt-1] - amu_HVP_min_Tdata.distr_list[tdata_opt-1];
    if( diff_max_min.ave()/min( amu_HVP_max_Tdata.err(tdata_opt-1) , amu_HVP_min_Tdata.err(tdata_opt-1)) < 0.3) Found_Tdata_opt=true;
    else tdata_opt++;
  }

      

  //if tdata_opt has not been found return tdata = -1 and amu_HVP =amu_HVP(T/2)
  if(!Found_Tdata_opt) {
	
    Tdatas_opt = -1;
    amu_HVP =  amu_HVP_func(T/2);
	
  }
  else { //tdata_opt has been found


    bool fit_to_constant=true;
    //average over min max over an interval of 0.5 fm
    if( amu_HVP_max_Tdata.size() != amu_HVP_min_Tdata.size()) crash("Size of amu_HVP_max_Tdata and amu_HVP_min_Tdata are not equal");
    int S_Tdata= amu_HVP_min_Tdata.size();
    int DT = min( (int)(DTT*fm_to_inv_Gev/a.ave()), (S_Tdata -tdata_opt));
    distr_t amu_HVP_ave_max= amu_HVP_max_Tdata.distr_list[tdata_opt-1];
    distr_t amu_HVP_ave_min= amu_HVP_min_Tdata.distr_list[tdata_opt-1];
    if(fit_to_constant) {
      amu_HVP_ave_max = amu_HVP_ave_max*1.0/pow(amu_HVP_max_Tdata.err(tdata_opt-1),2);
      amu_HVP_ave_min = amu_HVP_ave_min*1.0/pow(amu_HVP_min_Tdata.err(tdata_opt-1),2);
    }
    double amu_HVP_weight_max = (fit_to_constant)?1.0/(pow(amu_HVP_max_Tdata.err(tdata_opt-1),2)):1.0;
    double amu_HVP_weight_min = (fit_to_constant)?1.0/(pow(amu_HVP_min_Tdata.err(tdata_opt-1),2)):1.0;
    for(int t=tdata_opt; t< tdata_opt+DT; t++) {

      double weight_max_t = (fit_to_constant)?1.0/pow(amu_HVP_max_Tdata.err(t),2):1.0;
      double weight_min_t = (fit_to_constant)?1.0/pow(amu_HVP_min_Tdata.err(t),2):1.0;
	 
      amu_HVP_ave_max = amu_HVP_ave_max + amu_HVP_max_Tdata.distr_list[t]*weight_max_t;
      amu_HVP_ave_min = amu_HVP_ave_min + amu_HVP_min_Tdata.distr_list[t]*weight_min_t;
      amu_HVP_weight_max += weight_max_t;
      amu_HVP_weight_min += weight_min_t;

    }

    amu_HVP_ave_max= amu_HVP_ave_max/amu_HVP_weight_max;
    amu_HVP_ave_min= amu_HVP_ave_min/amu_HVP_weight_min;


    double den_weight= 1.0/(pow(amu_HVP_ave_max.err(),2))  + 1.0/(pow(amu_HVP_ave_min.err(),2));
    double fin_weight_max= (fit_to_constant)?(1.0/(pow(amu_HVP_ave_max.err(),2)*den_weight)):0.5;
    double fin_weight_min= (fit_to_constant)?(1.0/(pow(amu_HVP_ave_min.err(),2)*den_weight)):0.5;

    Tdatas_opt =tdata_opt;
    amu_HVP = fin_weight_min*amu_HVP_ave_min + fin_weight_max*amu_HVP_ave_max;
    Is_T_data_opt[tdata_opt-1] = 1.0;

  }

  Tcut_opt= Tdatas_opt;
      
  //Print to File
  Print_To_File({}, {TCUTS, amu_HVP_min_Tdata.ave(), amu_HVP_min_Tdata.err(), amu_HVP_max_Tdata.ave(), amu_HVP_max_Tdata.err(), amu_HVP_T_2.ave(), amu_HVP_T_2.err(), Is_T_data_opt}, path+".bound", "", "#tcut   lower    upper     T/2      Is_Tcut_opt.       Tcut_opt= "+to_string(Tdatas_opt));



  return;

}






void HVP() {





  int NB=2000;

  double GF=  1.1663787*1e-5;
  double hbar= 6.582119569*1e-25;
  double mt= 1.77686;
  double mk= 0.493677;
  double mpi=0.13957039;
  double mpi0=0.1349768;
  distr_t B_TK(false);
  distr_t B_TP(false);
  distr_t FK(false);
  distr_t FP(false);
  distr_t SEW(false);
  distr_t dR_TK(false);
  distr_t dR_TP(false);
  distr_t LT_tau(false);
  GaussianMersenne GV(645611);
  distr_t Vus(false);
  distr_t Vud(false);

  for(int iboot=0;iboot<NB;iboot++) {
    B_TK.distr.push_back(   (6.957 + GV()*0.096)*1e-3);
    B_TP.distr.push_back(   (0.10808 + GV()*0.00053));
    FK.distr.push_back( 0.1557 + GV()*0.0003);
    FP.distr.push_back( 0.1302 + GV()*0.0008);
    SEW.distr.push_back( 1.0232 + GV()*0.0003);
    LT_tau.distr.push_back( (2.903+ GV()*0.005)*1e-13);
    dR_TK.distr.push_back(  1 + (-0.15 + GV()*0.57)*1e-2);
    dR_TP.distr.push_back(  1 + (-0.24 + GV()*0.56)*1e-2);

  }
  
  Vus = GF*GF*(1.0/(16*M_PI*hbar))*FK*FK*LT_tau*pow(mt,3)*pow( 1 - pow(mk/mt,2),2)*SEW*dR_TK/B_TK;

  Vus= 1.0/Vus;
  Vus= SQRT_D(Vus);

  Vud = GF*GF*(1.0/(16*M_PI*hbar))*FP*FP*LT_tau*pow(mt,3)*pow( 1 - pow(mpi/mt,2),2)*SEW*dR_TP/B_TP;

  Vud= 1.0/Vud;
  Vud= SQRT_D(Vud);


  cout<<"Vus: "<<Vus.ave()<<" +- "<<Vus.err()<<endl;
  cout<<"Vud: "<<Vud.ave()<<" +- "<<Vud.err()<<endl;
   

  

  
  

  

  int Njacks=50;
  bool UseJack=true;
  double qu= 2.0/3.0;
  double qd= -1.0/3.0;
  double fm_to_inv_Gev= 1.0/0.197327;

  //create directories
  boost::filesystem::create_directory("../data/HVP");
  boost::filesystem::create_directory("../data/HVP/Bounding");
  boost::filesystem::create_directory("../data/HVP/Corr");

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
							       
		      
    string rA = A_bis.substr(A_bis.length()-5);
    string rB = B_bis.substr(B_bis.length()-5);
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(rA.substr(1,1));
      int n2 = stoi(rB.substr(1,1));
      if(rA == rB) {
	if(rA=="r0.h5" || rA=="r2.h5") return conf_num_A > conf_num_B;
	else if(rA=="r1.h5" || rA=="r3.h5") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    return A_bis<B_bis;
  };
  
  
  data_t Vk_data_tm, Vk_data_OS, P5P5_data_tm;

  Vk_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "VKVK", Sort_light_confs);
  P5P5_data_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  Vk_data_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "VKVK", Sort_light_confs);


  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;

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

  }
      

  //############################################################################################



























  //################################## E ENSEMBLE ANALYSIS ####################################


  bool Get_ASCII= true;

    if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../E112_RCs");
    

    vector<string> Ens_T1({"E.44.112"});
    vector<string> Ens_TT1({"cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_l_l",  "mix_l_s1", "mix_l_s2", "mix_s1_s1", "mix_s2_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../E112_RCs/"+channel);
	boost::filesystem::create_directory("../E112_RCs/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_A0P5", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_A0P5"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"channel: "<<channel<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../E112_RCs/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../E112_RCs/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	  else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	  PrintCorr.close();
	
	}

	fclose(stream);

	}
	
      }
    }
    }


    auto Sort_easy = [](string A, string B) {

      int conf_length_A= A.length();
      int conf_length_B= B.length();
      
      int pos_a_slash=-1;
      int pos_b_slash=-1;
      for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
      for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;
      
      string A_bis= A.substr(pos_a_slash+1);
      string B_bis= B.substr(pos_b_slash+1);

      return atoi( A_bis.c_str()) < atoi( B_bis.c_str());
      
    };
   

    D(1);

    data_t V_strange_L, V_strange_OS_L, V_strange_M, V_strange_OS_M;
    data_t corr_P5P5_strange, corr_P5P5_OS_strange, corr_P5P5_strange_heavy, corr_P5P5_OS_strange_heavy;
    data_t corr_A0P5_strange, corr_A0P5_OS_strange, corr_A0P5_strange_heavy, corr_A0P5_OS_strange_heavy;
    //strange
    //L
    V_strange_L.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
    V_strange_OS_L.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);
    D(2);
    //M
    V_strange_M.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKVK", "VKVK", Sort_easy); 
    V_strange_OS_M.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKVK", "VKVK", Sort_easy);
    D(3);
    //P5P5
    corr_P5P5_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_P5P5", "P5P5", Sort_easy);
    corr_P5P5_OS_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_P5P5", "P5P5", Sort_easy);
    corr_P5P5_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_P5P5", "P5P5", Sort_easy);
    corr_P5P5_OS_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_P5P5","P5P5",Sort_easy);
    D(4);
    //A0P5
    corr_A0P5_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_TM_A0P5", "A0P5", Sort_easy);
    corr_A0P5_OS_strange.Read("../E112_RCs/mix_s1_s1", "mes_contr_mix_s1_s1_OS_A0P5", "A0P5", Sort_easy);
    corr_A0P5_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_TM_A0P5", "A0P5", Sort_easy);
    corr_A0P5_OS_strange_heavy.Read("../E112_RCs/mix_s2_s2", "mes_contr_mix_s2_s2_OS_A0P5", "A0P5", Sort_easy);
    D(5);
    //VKVK light
    data_t V_light_tm, V_light_OS;
    V_light_tm.Read("../E112_RCs/mix_l_l", "mes_contr_mix_l_l_TM_VKVK", "VKVK", Sort_easy);
    V_light_OS.Read("../E112_RCs/mix_l_l", "mes_contr_mix_l_l_OS_VKVK", "VKVK", Sort_easy);
    int id_E=0;
    
    D(6);
  
  CorrAnalysis Corr(UseJack, Njacks,100);
  Corr.Nt = V_strange_L.nrows[id_E];

  boost::filesystem::create_directory("../data/RC_analysis_E");
  

  distr_t_list corr_s_V_L_tm = Corr.corr_t( V_strange_L.col(0)[id_E], "");
  distr_t_list corr_s_V_L_OS = Corr.corr_t( V_strange_OS_L.col(0)[id_E], "");

  distr_t_list corr_s_V_M_tm = Corr.corr_t( V_strange_M.col(0)[id_E], "");
  distr_t_list corr_s_V_M_OS = Corr.corr_t( V_strange_OS_M.col(0)[id_E], "");

  distr_t_list corr_s_P_L_tm = Corr.corr_t( corr_P5P5_strange.col(0)[id_E], "");
  distr_t_list corr_s_P_L_OS = Corr.corr_t( corr_P5P5_OS_strange.col(0)[id_E], "");

  distr_t_list corr_s_P_M_tm = Corr.corr_t( corr_P5P5_strange_heavy.col(0)[id_E], "");
  distr_t_list corr_s_P_M_OS = Corr.corr_t( corr_P5P5_OS_strange_heavy.col(0)[id_E], "");

  distr_t_list corr_s_AP_L_tm = Corr.corr_t( corr_A0P5_strange.col(0)[id_E], "");
  distr_t_list corr_s_AP_L_OS = Corr.corr_t( corr_A0P5_OS_strange.col(0)[id_E], "");

  distr_t_list corr_s_AP_M_tm = Corr.corr_t( corr_A0P5_strange_heavy.col(0)[id_E], "");
  distr_t_list corr_s_AP_M_OS = Corr.corr_t( corr_A0P5_OS_strange_heavy.col(0)[id_E], "");

  distr_t_list corr_V_tm = (5.0/9.0)*Corr.corr_t( V_light_tm.col(0)[id_E], "");
  distr_t_list corr_V_OS = (5.0/9.0)*Corr.corr_t( V_light_OS.col(0)[id_E], "");


  Corr.Tmin=62;
  Corr.Tmax=100;


  D(3);
  
  LatticeInfo L_info;
     
  L_info.LatInfo_new_ens(V_light_tm.Tag[id_E]);
     
  double aml= L_info.ml;
  double ams1= L_info.ms_L_new;
  double ams2= L_info.ms_M_new;

  distr_t Meta_tm_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_tm, "../data/RC_analysis_E/eta_L_tm"));
  distr_t Meta_tm_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_tm, "../data/RC_analysis_E/eta_M_tm"));
  distr_t Meta_OS_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_OS, "../data/RC_analysis_E/eta_L_OS"));
  distr_t Meta_OS_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_OS, "../data/RC_analysis_E/eta_M_OS"));
  
  distr_t_list RV_L = 2*ams1*corr_s_P_L_tm/(distr_t_list::derivative(corr_s_AP_L_tm, 0));  
  distr_t_list RV_M = 2*ams2*corr_s_P_M_tm/(distr_t_list::derivative(corr_s_AP_M_tm, 0));

  distr_t_list RA_L = 2*ams1*corr_s_P_L_OS/(distr_t_list::derivative(corr_s_AP_L_OS, 0)); 
  distr_t_list RA_M = 2*ams2*corr_s_P_M_OS/(distr_t_list::derivative(corr_s_AP_M_OS, 0));


  RA_L = RA_L*(Meta_OS_L/Meta_tm_L)*(SINH_D(Meta_OS_L)/SINH_D(Meta_tm_L))*Corr.matrix_element_t(corr_s_P_L_tm, "")/Corr.matrix_element_t(corr_s_P_L_OS, "");
  RA_M = RA_M*(Meta_OS_M/Meta_tm_M)*(SINH_D(Meta_OS_M)/SINH_D(Meta_tm_M))*Corr.matrix_element_t(corr_s_P_M_tm, "")/Corr.matrix_element_t(corr_s_P_M_OS, "");

  Print_To_File({}, {RV_L.ave(), RV_L.err(), RV_M.ave(), RV_M.err()} , "../data/RC_analysis_E/RV.dat", "", "");
  Print_To_File({}, {RA_L.ave(), RA_L.err(), RA_M.ave(), RA_M.err()} , "../data/RC_analysis_E/RA.dat", "", "");

  Corr.Tmin=47;
  Corr.Tmax=87;

  distr_t Zv_L = Corr.Fit_distr(RV_L);
  distr_t Zv_M = Corr.Fit_distr(RV_M);
  distr_t Za_L = Corr.Fit_distr(RA_L);
  distr_t Za_M = Corr.Fit_distr(RA_M);


   


  
  vector<distr_t> M2etas_fit({Meta_tm_L*Meta_tm_L/(a_E*a_E), Meta_tm_M*Meta_tm_M/(a_E*a_E)});

  distr_t m_etas_phys_distr;
  const double m_etas = 0.68989;
  const double m_etas_err= 0.00050;
  
 for(int ijack=0;ijack<Njacks;ijack++) m_etas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err/sqrt(Njacks-1.0));
  
  //Generate fake ms_distr
  distr_t ms_light_distr;
  distr_t ms_heavy_distr;
  for(int ijack=0;ijack<Njacks;ijack++) ms_light_distr.distr.push_back( ams1 );
  for(int ijack=0;ijack<Njacks;ijack++) ms_heavy_distr.distr.push_back( ams2 );
  
  //estrapolate ms phys
  vector<distr_t> ms_list( {ms_light_distr, ms_heavy_distr});
  distr_t ms_phys_extr = Obs_extrapolation_meson_mass(ms_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr,  "../data/RC_analysis_E", "ms_extrapolation_etas", UseJack, "SPLINE");

  vector<distr_t> Za_hadr_list, Zv_hadr_list;
  Za_hadr_list = {Za_L, Za_M};
  Zv_hadr_list = {Zv_L, Zv_M};

  distr_t Za = Obs_extrapolation_meson_mass( Za_hadr_list, ms_list, ms_phys_extr, "../data/RC_analysis_E", "Za_extr_quark_mass", UseJack, "SPLINE");
  distr_t Zv = Obs_extrapolation_meson_mass( Zv_hadr_list, ms_list, ms_phys_extr, "../data/RC_analysis_E", "Zv_extr_quark_mass", UseJack, "SPLINE");

  //evaluate correlator

  int T= Corr.Nt;

  distr_t agm2_SD_tm(UseJack, UseJack?Njacks:1000);
  distr_t agm2_SD_OS(UseJack, UseJack?Njacks:1000);

  double t0=  0.4*fm_to_inv_Gev;
  double Delta= 0.15*fm_to_inv_Gev;

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  auto th0 = [&t0, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a_E , T/2);


  //free correlator
  //Generate_free_corr_data();

  
  string Pt_free_oppor= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/OPPOR";
  Vfloat VV_free_oppor= Read_From_File(Pt_free_oppor, 1, 4);
  if(VV_free_oppor.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml w opposite r");
  //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
  string Pt_free_samer= "../Vkvk_cont/"+to_string(Corr.Nt/2)+"_m"+to_string_with_precision(L_info.ml,5)+"/SAMER";
  Vfloat VV_free_samer= Read_From_File(Pt_free_samer, 1, 4);
  if(VV_free_samer.size() != Corr.Nt) crash("Failed to read properly free VV correlator ml  w same r");
  //Insert electric charges
  for( auto & OP:VV_free_oppor) OP *= 5.0/9.0;
  for( auto & SA:VV_free_samer) SA *= 5.0/9.0;

   corr_V_tm = (1.0/(Za*Za))*(Za*Za*corr_V_tm + VV_free_oppor); //free_corr_log_art
   corr_V_OS = (1.0/(Zv*Zv))*(Zv*Zv*corr_V_OS + VV_free_samer); //free_corr_log_art

 

  for(int t=1; t <Corr.Nt/2; t++) {

    agm2_SD_tm = agm2_SD_tm + 4.0*pow(alpha,2)*w(t,3)*Za*Za*corr_V_tm[t]*Ker[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_E));
    agm2_SD_OS = agm2_SD_OS + 4.0*pow(alpha,2)*w(t,3)*Zv*Zv*corr_V_OS[t]*Ker[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_E));
  }


  cout<<"Za: "<<Za.ave()<<" "<<Za.err()<<endl;
  cout<<"Zv: "<<Zv.ave()<<" "<<Zv.err()<<endl;
  
  
  cout<<"tm: "<<agm2_SD_tm.ave()<<" "<<agm2_SD_tm.err()<<endl;
  cout<<"OS: "<<agm2_SD_OS.ave()<<" "<<agm2_SD_OS.err()<<endl;

  cout<<"a: "<<(a_E/fm_to_inv_Gev).ave()<<" "<<(a_E/fm_to_inv_Gev).err()<<endl;


  

















  //#############################################################################################




  int Nens= Vk_data_tm.size;


  for(int iens=0;iens<Nens;iens++) {

    cout<<"Analyzing ensemble: "<<Vk_data_tm.Tag[iens]<<endl; 

    CorrAnalysis Corr(UseJack, Njacks,800);
    Corr.Nt = Vk_data_tm.nrows[iens];
     
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);

    if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; }
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
    else crash("Ensemble not found");

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);

    int L= L_info.L;

      
       
    distr_t_list Vk_tm_distr = 1e10*Za*Za*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_tm.col(0)[iens] , "../data/HVP/Corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
    distr_t_list Vk_OS_distr = 1e10*Zv*Zv*(pow(qu,2)+pow(qd,2))*Corr.corr_t( Vk_data_OS.col(0)[iens] , "../data/HVP/Corr/Vk_OS_"+Vk_data_OS.Tag[iens]+".dat");
    distr_t_list P5_distr= Corr.corr_t( P5P5_data_tm.col(0)[iens] , "");

    if(Vk_data_tm.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=30; Corr.Tmax=70;}
    else if(Vk_data_tm.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=27; Corr.Tmax=50;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60;}
    else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80;}
    else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+Vk_data_tm.Tag[iens]);
    
    distr_t aMpi_distr= Corr.Fit_distr(Corr.effective_mass_t( P5_distr, ""));

    distr_t p2_rest= 2*aMpi_distr;

    distr_t p2_mot= 2*SQRT_D( aMpi_distr*aMpi_distr + pow( 2*M_PI/L,2)); 

    //B64 110 MeV
    distr_t C= (0.110 - 0.135)/(a_B*a_B);
    distr_t aMpi_N= (0.135 +  C*a_distr*a_distr)*a_distr;
    distr_t p2_mot_tm = 2*SQRT_D( (0.5*(aMpi_distr + aMpi_N))*(0.5*(aMpi_distr + aMpi_N)) + pow(2*M_PI/L,2));

    cout<<"pi neutral: "<< (aMpi_N/a_distr).ave()<< " +- "<< (aMpi_N/a_distr).err()<<endl;

   
    distr_t amu_HVP_tm(UseJack);
    distr_t amu_HVP_OS(UseJack);
    int Tcut_opt_tm, Tcut_opt_OS;

  
          
    Bounding_HVP(amu_HVP_tm, Tcut_opt_tm,  Vk_tm_distr, a_distr,"../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_tm" , p2_mot_tm);
    Bounding_HVP(amu_HVP_OS, Tcut_opt_OS, Vk_OS_distr, a_distr, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+"_OS" , p2_mot);

  

   
    vector<string> Tags({"tm", "OS"});
    distr_t_list amu_HVP;
    amu_HVP.distr_list.push_back(amu_HVP_tm);
    amu_HVP.distr_list.push_back(amu_HVP_OS);
    vector<double> Tmins({(double)Tcut_opt_tm, (double)Tcut_opt_OS});
    vector<double> Tmaxs({ Tcut_opt_tm + DTT*fm_to_inv_Gev/a_distr.ave(), Tcut_opt_OS + DTT*fm_to_inv_Gev/a_distr.ave()});

    Print_To_File(Tags, { amu_HVP.ave(), amu_HVP.err(), Tmins, Tmaxs}, "../data/HVP/Bounding/"+Vk_data_tm.Tag[iens]+".res", "", "");

    cout<<"#### "<<Vk_data_tm.Tag[iens]<<" ###"<<endl;
    cout<<"HVP tm: "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<" stat. "<< (amu_HVP_tm.err()*100/amu_HVP_tm.ave())<<"%"<<endl;
    cout<<"HVP OS: "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<" stat. "<< (amu_HVP_OS.err()*100/amu_HVP_OS.ave())<<"%"<<endl;
    cout<<"#######"<<endl;

    auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
    //distr_t_list Ker_el = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(0.0005/0.10565837) , Corr.Nt);
    //distr_t_list Ker_mu = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(1.77/0.1056837) , Corr.Nt);



    

   
       
  }


  return;

}
