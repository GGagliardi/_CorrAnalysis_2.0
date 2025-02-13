#include "../include/sea_quark_effects.h"
#include "Corr_analysis.h"
#include "input.h"
#include "numerics.h"


const double alpha = 1.0 / 137.035999;




using namespace std;



void valence_n_sea_quark_effects() {

  bool UseJack = true;

  sea_quark_effects();

  exit(-1);

  int Njacks = 100;
  
  double fm_to_inv_Gev= 1.0/0.197327;
  double qu = 2.0 / 3;
  double qd = -qu / 2;

  
   //generate lattice spacings and RCs
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
   
   

  auto Sort_light_confs = [](string A, string B) {

			   

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			    //A_bis=A;
			    //B_bis=B;

			     
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
  
 
  




  //get blinded VKVK correlator

  bool Get_ASCII_BLINDED=false;

  if(Get_ASCII_BLINDED) {
    //read binary
    vector<string> Corr_tags({"VkVk_OS.bin", "VkVk_tm.bin"});

    vector<string> CH({"OS", "TM"});
    
    for(int id=0; id<(signed)Corr_tags.size(); id++) {
      
      vector<string> C_list;
      ifstream CONFS_LIST;
      CONFS_LIST.open("../sea_quark_effects/HVP/cB211b.072.64/Conf_list.txt");
      while(!CONFS_LIST.eof()) {
	string a;
	CONFS_LIST>>a;
	if(!CONFS_LIST.eof()) C_list.push_back(a);
      }

      CONFS_LIST.close();
      
    
      
      
      
    
      string ch= Corr_tags[id];
      
      FILE *stream = fopen(("../../PEAKY_BLINDER/blinded_confs_RM3_test/B64/"+ch).c_str(), "rb");
      size_t Nconfs, TT, Nhits, Nsubs;
      bin_read(Nconfs, stream);
      bin_read(Nhits, stream);
      bin_read(TT, stream);
      bin_read(Nsubs, stream);
      
      
    
      
      cout<<"Nconfs: "<<Nconfs<<endl;
      cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
      cout<<"Nhits: "<<Nhits<<endl;
      cout<<"Nsubs: "<<Nsubs<<endl;
      
      
    
      
      for(size_t iconf=0;iconf<Nconfs;iconf++) {
	vector<double> C(TT/2+1,0.0);
	for(size_t t=0;t<TT/2+1;t++) {
	  for(size_t is=0;is<Nsubs;is++) {
	  
	    double c;
	    bin_read(c, stream);
	    C[t] += c/Nsubs;
	    
	  }
	}
	boost::filesystem::create_directory("../sea_quark_effects/HVP/cB211b.072.64/"+C_list[iconf]);
	ofstream PrintCorr("../sea_quark_effects/HVP/cB211b.072.64/"+C_list[iconf]+"/mes_contr_"+CH[id]+"_VKVK");
	PrintCorr.precision(16);
	PrintCorr<<"# VKVK"<<endl;
	for(size_t t=0;t<(TT/2+1);t++) PrintCorr<<C[t]<<endl;
	for(size_t t=TT/2+1; t<TT;t++) PrintCorr<<C[TT-t]<<endl;
	PrintCorr.close();

      }
    
      
      fclose(stream);
      
    }
  }
  

 
  //get ll,  ls and ss correlators
  
  bool Get_ASCII_Mix= false;
  
    if(Get_ASCII_Mix) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects");
      boost::filesystem::create_directory("../sea_quark_effects/Mix");

      vector<string> Ens_T1({"B.72.64"});
      vector<string> Ens_TT1({"cB211b.072.64"});
      
    
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({"mix_l_s1", "mix_s1_s1", "ll"});
	
	
	for(auto &channel : channels) {
	  boost::filesystem::create_directory("../sea_quark_effects/Mix/"+channel);
	  boost::filesystem::create_directory("../sea_quark_effects/Mix/"+channel+"/"+Ens_TT1[it]);
	}
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_VKVK", "TM_P5P5", "TM_V0P5", "OS_VKVK"});
	
	
	for(int id=0; id<(signed)Corr_tags.size(); id++) {
	  
	  
	  
	  for( auto &channel: channels) {
	    
	    
	    vector<string> C_list;
	    ifstream CONFS_LIST;
	    
	    if(channel != "ll") CONFS_LIST.open("../sea_quark_effects/confsListMix_ls_ss");
	    else CONFS_LIST.open("../sea_quark_effects/confsListll");
	    while(!CONFS_LIST.eof()) {
	      string a;
	      CONFS_LIST>>a;
	      if(!CONFS_LIST.eof()) C_list.push_back(a);
	    }
	    
	    
	    CONFS_LIST.close();
	    
	    FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	    size_t Nconfs, T, Nhits;
	    bin_read(Nconfs, stream);
	    bin_read(Nhits, stream);
	    bin_read(T, stream);
	    cout<<"Nconfs: "<<Nconfs<<endl;
	    if(Nconfs != C_list.size()) crash("number of configs different than expected: "+to_string(C_list.size()));
	    cout<<"T: "<<T<<" "<<T/2+1<<endl;
	    cout<<"Nhits: "<<Nhits<<endl;
	    for(size_t iconf=0;iconf<Nconfs;iconf++) {
	      vector<double> C(T/2+1);
	      for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	      boost::filesystem::create_directory("../sea_quark_effects/Mix/"+channel+"/"+Ens_TT1[it]+"/"+C_list[iconf]);
	      ofstream PrintCorr("../sea_quark_effects/Mix/"+channel+"/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	      PrintCorr.precision(16);
	      PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	      for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	      if(Corr_tags[id].substr(3,4) == "V0P5" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	      else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	      PrintCorr.close();
	      
	    }
	    
	    fclose(stream);
	    
	  }
	  
	}
      }
    }
    
    
   
    
  data_t data_bub_P;
  data_t data_bub_S;
  data_t data_bub_P_strange;
  data_t data_bub_S_strange;
  
  
  //scalar and pseudoscalar loop
  data_bub_P.Read("../sea_quark_effects/loops", "bub_P5_diff.txt", "P5P5", Sort_light_confs);
  data_bub_S.Read("../sea_quark_effects/loops", "bub_S.txt", "S0S0", Sort_light_confs);
  //strange
  data_bub_P_strange.Read("../sea_quark_effects/strange_loops", "bub_P5_diff.txt", "P5P5", Sort_light_confs);
  data_bub_S_strange.Read("../sea_quark_effects/strange_loops", "bub_S.txt", "S0S0", Sort_light_confs);
   
  
  //connected tm and OS correlators
  //data_t data_P5P5_OS, data_P5P5_tm;
  //data_P5P5_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs_h5);
  //data_P5P5_tm.Read("../R_ratio_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs_h5);

  //3pt with insertion of pseudoscalar density

  data_t P5P5_ll_dmcr_up, P5P5_ll_dmcr_down;
  data_t P5P5_ls_dmcr_up, P5P5_ls_dmcr_down;
  data_t V0P5_ll_dmcr_up, V0P5_ll_dmcr_down;
  data_t VKVK_ll_dmcr_up, VKVK_ll_dmcr_down;
  data_t VKVK_ss_dmcr_up, VKVK_ss_dmcr_down;

  P5P5_ll_dmcr_up.Read("../LIBE_data", "dmcr_plus_ll_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  P5P5_ll_dmcr_down.Read("../LIBE_data", "dmcr_minus_ll_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  P5P5_ls_dmcr_up.Read("../LIBE_data", "dmcr_plus_ls_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  P5P5_ls_dmcr_down.Read("../LIBE_data", "dmcr_minus_ls_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  V0P5_ll_dmcr_up.Read("../LIBE_data", "dmcr_plus_ll_muval_+7.2000e-04.dat",  "V0P5", Sort_light_confs);
  V0P5_ll_dmcr_down.Read("../LIBE_data", "dmcr_minus_ll_muval_+7.2000e-04.dat",  "V0P5", Sort_light_confs);
  VKVK_ll_dmcr_up.Read("../LIBE_data", "dmcr_plus_ll_muval_+7.2000e-04.dat",  "VKVK", Sort_light_confs);
  VKVK_ll_dmcr_down.Read("../LIBE_data", "dmcr_minus_ll_muval_+7.2000e-04.dat",  "VKVK", Sort_light_confs);
  VKVK_ss_dmcr_up.Read("../LIBE_data", "dmcr_plus_ll_muval_+1.8250e-02.dat",  "VKVK", Sort_light_confs);
  VKVK_ss_dmcr_down.Read("../LIBE_data", "dmcr_minus_ll_muval_+1.8250e-02.dat",  "VKVK", Sort_light_confs);

  //3pt with insertion of scalar density
  data_t P5P5_ll_dm_up, P5P5_ll_dm_down;
  P5P5_ll_dm_up.Read("../LIBE_data", "dm_plus_ll_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  P5P5_ll_dm_down.Read("../LIBE_data", "dm_minus_ll_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
  data_t P5P5_ls_dm_up;
  P5P5_ls_dm_up.Read("../LIBE_data", "dm_minus_ls_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
 

  //VKVK ll correlator in isoQCD
  data_t VKVK_ll_iso;
  VKVK_ll_iso.Read("../sea_quark_effects/HVP", "mes_contr_TM_VKVK",  "VKVK", Sort_light_confs);
  data_t VKVK_ll_iso_old;
  VKVK_ll_iso_old.Read("../sea_quark_effects/Mix/ll", "mes_contr_ll_TM_VKVK",  "VKVK", Sort_light_confs);

  //VKVK ss correlator in isoQCD
  data_t VKVK_ss_iso, VKVK_ss_OS_iso;
  VKVK_ss_iso.Read("../sea_quark_effects/Mix/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK",  "VKVK", Sort_light_confs);
  VKVK_ss_OS_iso.Read("../sea_quark_effects/Mix/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK",  "VKVK", Sort_light_confs);
  data_t V1V1_ss_iso_old, V2V2_ss_iso_old, V3V3_ss_iso_old;
  V1V1_ss_iso_old.Read("../sea_quark_effects/old_strange", "mes_contr_2pts_ll_1",  "V1V1", Sort_light_confs);
  V2V2_ss_iso_old.Read("../sea_quark_effects/old_strange", "mes_contr_2pts_ll_1",  "V2V2", Sort_light_confs);
  V3V3_ss_iso_old.Read("../sea_quark_effects/old_strange", "mes_contr_2pts_ll_1",  "V3V3", Sort_light_confs);
  
  //V0P5 ll correlator in isoQCD
  data_t V0P5_ll_iso;
  V0P5_ll_iso.Read("../sea_quark_effects/Mix/V0P5", "mes_contr_ll_TM_V0P5",  "V0P5", Sort_light_confs);
  
  
  //P5P5 ll correlator in isoQCD
  data_t P5P5_ll_iso;
  P5P5_ll_iso.Read("../sea_quark_effects/Mix/ll", "mes_contr_ll_TM_P5P5",  "P5P5", Sort_light_confs); 
  //P5P5 ls correlator in isoQCD
  data_t P5P5_ls_iso;
  P5P5_ls_iso.Read("../sea_quark_effects/Mix/mix_l_s1", "mes_contr_mix_l_s1_TM_P5P5",  "P5P5", Sort_light_confs);



  boost::filesystem::create_directory("../data/sea_quark_effects");
  boost::filesystem::create_directory("../data/sea_quark_effects/B64");

  int iens=0;
  distr_t a_distr= a_B;
  distr_t ZV= ZV_B;
  distr_t ZA= ZA_B;
  CorrAnalysis Corr(UseJack,Njacks,100);
  int T=128;
  double L = 64.0;
  double V= pow(L,3);
  Corr.Nt=T;
  

  //needed isoQCD correlators
  distr_t_list VKVK_ll_iso_distr= Corr.corr_t(VKVK_ll_iso.col(0)[iens],"../data/sea_quark_effects/B64/VKVK_ll_iso");
  distr_t_list VKVK_ll_iso_old_distr= Corr.corr_t(VKVK_ll_iso_old.col(0)[iens],"../data/sea_quark_effects/B64/VKVK_ll_iso");
  distr_t_list VKVK_ss_iso_distr= Corr.corr_t(VKVK_ss_iso.col(0)[iens], "../data/sea_quark_effects/B64/VKVK_ss_iso");
  distr_t_list VKVK_ss_iso_old_distr= (1.0/3.0)*Corr.corr_t(summ_master(V1V1_ss_iso_old.col(0)[iens], V2V2_ss_iso_old.col(0)[iens], V3V3_ss_iso_old.col(0)[iens]   ), "");
  distr_t_list P5P5_ll_iso_distr= Corr.corr_t(P5P5_ll_iso.col(0)[iens], "../data/sea_quark_effects/B64/P5P5_ll_iso");
  distr_t_list P5P5_ls_iso_distr= Corr.corr_t(P5P5_ls_iso.col(0)[iens], "../data/sea_quark_effects/B64/P5P5_ls_iso");
  Corr.Reflection_sign=1;
  distr_t_list V0P5_ll_iso_distr= Corr.corr_t(V0P5_ll_iso.col(0)[iens], "../data/sea_quark_effects/B64/V0P5_ll_iso");
  Corr.Reflection_sign=1;

  //correlators with insertion of dmcrit

  distr_t_list VKVK_ll_dmcr_distr= (2.0/V)*Corr.corr_t(VKVK_ll_dmcr_up.col(0)[iens], "");
  distr_t_list P5P5_ll_dmcr_distr= (2.0/V)*Corr.corr_t(P5P5_ll_dmcr_up.col(0)[iens], "");
  distr_t_list VKVK_ss_dmcr_distr= (2.0/V)*Corr.corr_t(VKVK_ss_dmcr_up.col(0)[iens], "");
  distr_t_list P5P5_ls_dmcr_distr= (1.0/V)*(Corr.corr_t(P5P5_ls_dmcr_up.col(0)[iens], "") + Corr.corr_t(P5P5_ls_dmcr_down.col(0)[iens],""));  //PAY ATTENTION TO THE SIGN
  Corr.Reflection_sign=-1;
  distr_t_list V0P5_ll_dmcr_distr= (2.0/V)*Corr.corr_t(V0P5_ll_dmcr_up.col(0)[iens], "");
  Corr.Reflection_sign=1;


  //correlators with insertion of dm
  distr_t_list P5P5_ll_dm_distr= (2.0/V)*Corr.corr_t(P5P5_ll_dm_up.col(0)[iens], "");
  distr_t_list P5P5_ls_dm_distr= (1.0/V)*Corr.corr_t(P5P5_ls_dm_up.col(0)[iens], "");


  //evaluate sea quark mass derivatives

  //define scalar and pseudoscalar bubble for each conf;

  int Nconfs_BUB=  data_bub_P.col(0)[0][0].size();

  if(Nconfs_BUB != (signed)VKVK_ll_iso.col(0)[iens][0].size()) crash("number of configs in loops and 2pts(deflated) is not the same");
  if(Nconfs_BUB != (signed)P5P5_ll_iso.col(0)[iens][0].size()) crash("number of configs in loops and 2pts is not the same");
  if(Nconfs_BUB != (signed)data_bub_P_strange.col(0)[iens][0].size()) crash("number of configs in light and strange loops is not the same");

  Vfloat S(Nconfs_BUB,0.0), P(Nconfs_BUB,0.0);
  Vfloat S_WI(Nconfs_BUB,0.0);
  Vfloat S_strange(Nconfs_BUB,0.0), P_strange(Nconfs_BUB,0.0);
  

  for(int iconf=0; iconf<Nconfs_BUB;iconf++)  
    for(int t=0;t<T;t++) {
      S[iconf] += 2.0*data_bub_S.col(1)[iens][t][iconf];
      P[iconf] += -2.0*data_bub_P.col(0)[iens][t][iconf];
      S_WI[iconf] += 2*0.00072*pow(L,3)*T*P5P5_ll_iso.col(0)[iens][t][iconf];
      S_strange[iconf] += data_bub_S_strange.col(1)[iens][t][iconf];
      P_strange[iconf] += -data_bub_P_strange.col(0)[iens][t][iconf];
    }

  
  cout<<"Printing strange condensate "<<endl;
  for(int iconf=0; iconf<Nconfs_BUB;iconf++) {
    cout<<"iconf: "<<iconf<<" , "<<S_strange[iconf]<<endl;
  }

  
  //correlate bubble with all observables of interest

  VVfloat VKVK_ll_dmcr_D(T), P5P5_ll_dmcr_D(T), VKVK_ss_dmcr_D(T), P5P5_ls_dmcr_D(T), V0P5_ll_dmcr_D(T);
  VVfloat VKVK_ll_dm_D(T), VKVK_ss_dm_D(T);
  VVfloat VKVK_ll_dm_strange_D(T), VKVK_ss_dm_strange_D(T);
  VVfloat VKVK_ss_OS_dm_D(T), VKVK_ss_OS_dm_strange_D(T);
  

  //to compute light-quark mass corrections to Mpi and MK
  VVfloat P5P5_ll_dm_D(T), P5P5_ls_dm_D(T), P5P5_ll_dm_strange_D(T),  P5P5_ls_dm_strange_D(T);

  VVfloat VKVK_ll_dm_WI_D(T), VKVK_ss_dm_WI_D(T);
													


  
  auto SEA_FUNC= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function SEA_FUNC expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[2];};


  for(int t=0; t<T;t++) {
    for(int iconf=0;iconf<Nconfs_BUB;iconf++) {
      VKVK_ll_dmcr_D[t].push_back( P[iconf]*VKVK_ll_iso.col(0)[iens][t][iconf] );
      VKVK_ll_dm_D[t].push_back( S[iconf]*VKVK_ll_iso.col(0)[iens][t][iconf] );
      VKVK_ss_dm_D[t].push_back( S[iconf]*VKVK_ss_iso.col(0)[iens][t][iconf] );
      VKVK_ss_OS_dm_D[t].push_back( S[iconf]*VKVK_ss_OS_iso.col(0)[iens][t][iconf] );
      VKVK_ss_dmcr_D[t].push_back( P[iconf]*VKVK_ss_iso.col(0)[iens][t][iconf] );
      P5P5_ll_dmcr_D[t].push_back( P[iconf]*P5P5_ll_iso.col(0)[iens][t][iconf] );
      P5P5_ls_dmcr_D[t].push_back( P[iconf]*P5P5_ls_iso.col(0)[iens][t][iconf] );
      V0P5_ll_dmcr_D[t].push_back( P[iconf]*V0P5_ll_iso.col(0)[iens][t][iconf] );

      
      P5P5_ll_dm_D[t].push_back( S[iconf]*P5P5_ll_iso.col(0)[iens][t][iconf] );
      P5P5_ls_dm_D[t].push_back( S[iconf]*P5P5_ls_iso.col(0)[iens][t][iconf] );


      //strange quark loops correlation
      VKVK_ll_dm_strange_D[t].push_back( S_strange[iconf]*VKVK_ll_iso.col(0)[iens][t][iconf] );
      VKVK_ss_dm_strange_D[t].push_back( S_strange[iconf]*VKVK_ss_iso.col(0)[iens][t][iconf] );
      VKVK_ss_OS_dm_strange_D[t].push_back( S_strange[iconf]*VKVK_ss_OS_iso.col(0)[iens][t][iconf] );
      P5P5_ls_dm_strange_D[t].push_back( S_strange[iconf]*P5P5_ls_iso.col(0)[iens][t][iconf] );
      P5P5_ll_dm_strange_D[t].push_back( S_strange[iconf]*P5P5_ll_iso.col(0)[iens][t][iconf] );

      //using condensate from Ward-Identity
      VKVK_ll_dm_WI_D[t].push_back( S_WI[iconf]*VKVK_ll_iso.col(0)[iens][t][iconf] );
      VKVK_ss_dm_WI_D[t].push_back( S_WI[iconf]*VKVK_ss_iso.col(0)[iens][t][iconf] );
    }
  }

  Jackknife J(10000,Njacks);
  distr_t_list VKVK_ll_dmcr_D_distr(UseJack), VKVK_ll_dm_D_distr(UseJack),  P5P5_ll_dmcr_D_distr(UseJack), VKVK_ss_dmcr_D_distr(UseJack), VKVK_ss_dm_D_distr(UseJack), P5P5_ls_dmcr_D_distr(UseJack), V0P5_ll_dmcr_D_distr(UseJack);
  distr_t_list P5P5_ll_dm_D_distr(UseJack), P5P5_ls_dm_D_distr(UseJack);
  //strange
  distr_t_list VKVK_ll_dm_strange_D_distr(UseJack), VKVK_ss_dm_strange_D_distr(UseJack);
  distr_t_list VKVK_ll_dm_WI_D_distr(UseJack), VKVK_ss_dm_WI_D_distr(UseJack);
  distr_t_list P5P5_ls_dm_strange_D_distr(UseJack), P5P5_ll_dm_strange_D_distr(UseJack);
  distr_t_list VKVK_ss_OS_dm_D_distr(UseJack), VKVK_ss_OS_dm_strange_D_distr(UseJack);
  
  for(int t=0;t<T;t++) {
    VKVK_ll_dmcr_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ll_dmcr_D[t], VKVK_ll_iso.col(0)[iens][t], P ));
    VKVK_ll_dm_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ll_dm_D[t], VKVK_ll_iso.col(0)[iens][t], S ));
    VKVK_ss_dmcr_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_dmcr_D[t], VKVK_ss_iso.col(0)[iens][t], P));
    VKVK_ss_dm_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_dm_D[t], VKVK_ss_iso.col(0)[iens][t], S));
    VKVK_ss_OS_dm_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_OS_dm_D[t], VKVK_ss_OS_iso.col(0)[iens][t], S));
    P5P5_ll_dmcr_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ll_dmcr_D[t], P5P5_ll_iso.col(0)[iens][t], P));
    P5P5_ls_dmcr_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ls_dmcr_D[t], P5P5_ls_iso.col(0)[iens][t], P));
    V0P5_ll_dmcr_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, V0P5_ll_dmcr_D[t], V0P5_ll_iso.col(0)[iens][t], P));

    P5P5_ll_dm_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ll_dm_D[t], P5P5_ll_iso.col(0)[iens][t], S));
    P5P5_ls_dm_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ls_dm_D[t], P5P5_ls_iso.col(0)[iens][t], S));

    //strange
    VKVK_ll_dm_strange_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ll_dm_strange_D[t], VKVK_ll_iso.col(0)[iens][t], S_strange ));
    VKVK_ss_dm_strange_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_dm_strange_D[t], VKVK_ss_iso.col(0)[iens][t], S_strange ));
    VKVK_ss_OS_dm_strange_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_OS_dm_strange_D[t], VKVK_ss_OS_iso.col(0)[iens][t], S_strange ));
    P5P5_ls_dm_strange_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ls_dm_strange_D[t], P5P5_ls_iso.col(0)[iens][t], S_strange));
    P5P5_ll_dm_strange_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, P5P5_ll_dm_strange_D[t], P5P5_ll_iso.col(0)[iens][t], S_strange));

    //using condensate from WI
    VKVK_ll_dm_WI_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ll_dm_WI_D[t], VKVK_ll_iso.col(0)[iens][t], S_WI ));
    VKVK_ss_dm_WI_D_distr.distr_list.push_back( J.DoJack( SEA_FUNC, 3, VKVK_ss_dm_WI_D[t], VKVK_ss_iso.col(0)[iens][t], S_WI ));

  }

  //symmetrize sea-quark correlators
  for(int t=0;t<=T/2;t++) {
    distr_t x(UseJack);
    x= 0.5*(VKVK_ll_dmcr_D_distr[t] + VKVK_ll_dmcr_D_distr[(T-t)%T]);  VKVK_ll_dmcr_D_distr.distr_list[t] = x; VKVK_ll_dmcr_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ll_dm_D_distr[t] + VKVK_ll_dm_D_distr[(T-t)%T]);  VKVK_ll_dm_D_distr.distr_list[t] = x; VKVK_ll_dm_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_dm_D_distr[t] + VKVK_ss_dm_D_distr[(T-t)%T]);  VKVK_ss_dm_D_distr.distr_list[t] = x; VKVK_ss_dm_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_OS_dm_D_distr[t] + VKVK_ss_OS_dm_D_distr[(T-t)%T]);  VKVK_ss_OS_dm_D_distr.distr_list[t] = x; VKVK_ss_OS_dm_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_dmcr_D_distr[t] + VKVK_ss_dmcr_D_distr[(T-t)%T]);  VKVK_ss_dmcr_D_distr.distr_list[t] = x; VKVK_ss_dmcr_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ll_dmcr_D_distr[t] + P5P5_ll_dmcr_D_distr[(T-t)%T]);  P5P5_ll_dmcr_D_distr.distr_list[t] = x; P5P5_ll_dmcr_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ls_dmcr_D_distr[t] + P5P5_ls_dmcr_D_distr[(T-t)%T]);  P5P5_ls_dmcr_D_distr.distr_list[t] = x; P5P5_ls_dmcr_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ll_dm_D_distr[t] + P5P5_ll_dm_D_distr[(T-t)%T]);  P5P5_ll_dm_D_distr.distr_list[t] = x; P5P5_ll_dm_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ls_dm_D_distr[t] + P5P5_ls_dm_D_distr[(T-t)%T]);  P5P5_ls_dm_D_distr.distr_list[t] = x; P5P5_ls_dm_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(V0P5_ll_dmcr_D_distr[t] + V0P5_ll_dmcr_D_distr[(T-t)%T]);  V0P5_ll_dmcr_D_distr.distr_list[t] = x; V0P5_ll_dmcr_D_distr.distr_list[(T-t)%T] = -1.0*x;

    //strange
    x= 0.5*(VKVK_ll_dm_strange_D_distr[t] + VKVK_ll_dm_strange_D_distr[(T-t)%T]);  VKVK_ll_dm_strange_D_distr.distr_list[t] = x; VKVK_ll_dm_strange_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_dm_strange_D_distr[t] + VKVK_ss_dm_strange_D_distr[(T-t)%T]);  VKVK_ss_dm_strange_D_distr.distr_list[t] = x; VKVK_ss_dm_strange_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_OS_dm_strange_D_distr[t] + VKVK_ss_OS_dm_strange_D_distr[(T-t)%T]);  VKVK_ss_OS_dm_strange_D_distr.distr_list[t] = x; VKVK_ss_OS_dm_strange_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ls_dm_strange_D_distr[t] + P5P5_ls_dm_strange_D_distr[(T-t)%T]);  P5P5_ls_dm_strange_D_distr.distr_list[t] = x; P5P5_ls_dm_strange_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(P5P5_ll_dm_strange_D_distr[t] + P5P5_ll_dm_strange_D_distr[(T-t)%T]);  P5P5_ll_dm_strange_D_distr.distr_list[t] = x; P5P5_ll_dm_strange_D_distr.distr_list[(T-t)%T] = x;

    //using condensate from WI
    x= 0.5*(VKVK_ll_dm_WI_D_distr[t] + VKVK_ll_dm_WI_D_distr[(T-t)%T]);  VKVK_ll_dm_WI_D_distr.distr_list[t] = x; VKVK_ll_dm_WI_D_distr.distr_list[(T-t)%T] = x;
    x= 0.5*(VKVK_ss_dm_WI_D_distr[t] + VKVK_ss_dm_WI_D_distr[(T-t)%T]);  VKVK_ss_dm_WI_D_distr.distr_list[t] = x; VKVK_ss_dm_WI_D_distr.distr_list[(T-t)%T] = x;
  }
      


  

  //determine dmcrit shift

  distr_t_list V0P5_ll_dmcr_TOT = -1.0*V0P5_ll_dmcr_distr + V0P5_ll_dmcr_D_distr;  //PAY ATTENTION TO THE RELATIVE SIGN
  distr_t_list P5P5_ll_dmcr_TOT = -1.0*P5P5_ll_dmcr_distr + P5P5_ll_dmcr_D_distr;  //PAY ATTENTION TO THE RELATIVE SIGN

  double mul = 0.00072;
  //mPCAC
  distr_t_list mPCAC_distr=  -0.5*distr_t_list::derivative(V0P5_ll_iso_distr,0)/(mul*P5P5_ll_iso_distr);

  //derivative  dmPCAC/dmcrit

  distr_t_list dmPCAC_distr= -0.5*(distr_t_list::derivative(V0P5_ll_dmcr_TOT,0)/P5P5_ll_iso_distr   -0.0*distr_t_list::derivative(V0P5_ll_iso_distr,0)*P5P5_ll_dmcr_TOT/(P5P5_ll_iso_distr*P5P5_ll_iso_distr));

 

  Print_To_File({}, {dmPCAC_distr.ave(), dmPCAC_distr.err()},  "../data/sea_quark_effects/B64/dmPCAC", "", "");
  Print_To_File({}, {mPCAC_distr.ave(), mPCAC_distr.err()},  "../data/sea_quark_effects/B64/mPCAC", "", "");

  Corr.Tmin=14;
  Corr.Tmax=45;

  //REQUIRE dmpcac = 0.05 *mul = 0.05*0.00072 -> dmcrit = (dmcrit/dmpcac)*0.05*mul = (dmcrit/dmpcac)*0.05*0.00072
  distr_t dmcrit= 0.05*0.00072/Corr.Fit_distr(dmPCAC_distr);

  double dm= (0.0006675 - 0.00072);

  cout<<"dmpcac/dmcrit: "<<Corr.Fit_distr(dmPCAC_distr).ave()<<" +- "<<Corr.Fit_distr(dmPCAC_distr).err()<<endl;
  cout<<"Typical size of dmcrit: "<<dmcrit.ave()<<" +- "<<dmcrit.err()<<endl;


  P5P5_ll_dmcr_TOT = dmcrit*P5P5_ll_dmcr_TOT;
  //correct all observables

  distr_t_list VKVK_ll_dmcr_TOT = dmcrit*(VKVK_ll_dmcr_distr + VKVK_ll_dmcr_D_distr);  //PAY ATTENTION TO THE RELATIVE SIGN
  distr_t_list VKVK_ss_dmcr_TOT = dmcrit*(VKVK_ss_dmcr_distr + VKVK_ss_dmcr_D_distr);  //PAY ATTENTION TO THE RELATIVE SIGN

  distr_t_list P5P5_ls_dmcr_TOT = dmcrit*(-1.0*P5P5_ls_dmcr_distr + P5P5_ls_dmcr_D_distr); //PAY ATTENTION TO THE RELATIVE SIGN
  distr_t_list P5P5_ls_dm_TOT = dm*(P5P5_ls_dm_distr + P5P5_ls_dm_D_distr); //PAY ATTENTION TO THE RELATIVE SIGN
  distr_t_list P5P5_ll_dm_TOT = dm*(-1.0*P5P5_ll_dm_distr + P5P5_ll_dm_D_distr);  //PAY ATTENTION TO THE RELATIVE SIGN

  Print_To_File({}, {(V0P5_ll_dmcr_distr).ave(), (V0P5_ll_dmcr_distr).err(), V0P5_ll_dmcr_D_distr.ave(), V0P5_ll_dmcr_D_distr.err() },  "../data/sea_quark_effects/B64/dV0P5_ll_dmcr", "", "");
  Print_To_File({}, {(P5P5_ll_dmcr_distr).ave(), (P5P5_ll_dmcr_distr).err(), P5P5_ll_dmcr_D_distr.ave(), P5P5_ll_dmcr_D_distr.err() },  "../data/sea_quark_effects/B64/dP5P5_ll_dmcr", "", "");
  Print_To_File({}, {(P5P5_ls_dmcr_distr).ave(), (P5P5_ls_dmcr_distr).err(), P5P5_ls_dmcr_D_distr.ave(), P5P5_ls_dmcr_D_distr.err() },  "../data/sea_quark_effects/B64/dP5P5_ls_dmcr", "", "");
  Print_To_File({}, {(VKVK_ll_dmcr_distr).ave(), (VKVK_ll_dmcr_distr).err(), VKVK_ll_dmcr_D_distr.ave(), VKVK_ll_dmcr_D_distr.err() },  "../data/sea_quark_effects/B64/dVKVK_ll_dmcr", "", "");
  Print_To_File({}, {(VKVK_ss_dmcr_distr).ave(), (VKVK_ss_dmcr_distr).err(), VKVK_ss_dmcr_D_distr.ave(), VKVK_ss_dmcr_D_distr.err() },  "../data/sea_quark_effects/B64/dVKVK_ss_dmcr", "", "");

  Corr.Tmin=35; Corr.Tmax=55;
  distr_t mP= Corr.Fit_distr(Corr.effective_mass_t( P5P5_ll_iso_distr, ""));
  Corr.Tmin=45; Corr.Tmax=59;
  distr_t mK= Corr.Fit_distr( Corr.effective_mass_t( P5P5_ls_iso_distr, ""));

  cout<<"mK/mP: "<<(mK/mP).ave()<<" +- "<<(mK/mP).err()<<endl;

  double dm_s = 0.10*0.0182782;

  distr_t_list dmPi = Corr.effective_slope_t( -1.0*P5P5_ll_dmcr_TOT/a_distr, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dmPi_eff_slope");
  distr_t_list dmK  = Corr.effective_slope_t( -1.0*P5P5_ls_dmcr_TOT/a_distr, P5P5_ls_iso_distr, "../data/sea_quark_effects/B64/dmK_eff_slope");
  distr_t_list dmK_s  = Corr.effective_slope_t( -1.0*dm_s*P5P5_ls_dm_strange_D_distr/a_distr, P5P5_ls_iso_distr, "../data/sea_quark_effects/B64/dmK_s_eff_slope");
  distr_t_list dmPi_s  = Corr.effective_slope_t( -1.0*dm_s*P5P5_ll_dm_strange_D_distr/a_distr, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dmPi_s_eff_slope");
  distr_t_list dmphi_s =  Corr.effective_slope_t( -1.0*dm_s*VKVK_ss_dm_strange_D_distr/a_distr, VKVK_ss_iso_distr, "../data/sea_quark_effects/B64/mphi_s_eff_slope");
  distr_t_list dA_pi = Corr.effective_ampl_rel_slope_t( P5P5_ll_dmcr_TOT, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dA_pi");
  distr_t_list dfPi_ov_fPi = Corr.decay_constant_rel_slope_t( P5P5_ll_dmcr_TOT, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dfPi_eff_slope");
  distr_t_list dfPi_ov_fPi_s = Corr.decay_constant_rel_slope_t( -1.0*dm_s*P5P5_ll_dm_strange_D_distr, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dfPi_s_eff_slope");
  distr_t_list dfK_ov_fK = Corr.decay_constant_rel_slope_t( P5P5_ls_dmcr_TOT, P5P5_ls_iso_distr, "../data/sea_quark_effects/B64/dfK_eff_slope");

  Corr.Tmin=18; Corr.Tmax=41;

  distr_t_list dmpi_ov_dmK=  Corr.effective_slope_t( -1.0*P5P5_ll_dm_TOT, P5P5_ll_iso_distr, "../data/sea_quark_effects/B64/dmPi_dm_eff_slope")/Corr.Fit_distr(Corr.effective_slope_t( -1.0*P5P5_ls_dm_TOT, P5P5_ls_iso_distr, "../data/sea_quark_effects/B64/dmK_dm_eff_slope")) ;

  distr_t_list ChPT_ratio= dmpi_ov_dmK*mP/mK;

  Print_To_File({}, {ChPT_ratio.ave(), ChPT_ratio.err() },  "../data/sea_quark_effects/B64/ChPT_ratio", "", "");
  
  //get HVP corrections
  
  distr_t_list damu_ll_dmcr_tcut_valence(UseJack), damu_ll_dmcr_tcut_sea(UseJack),  damu_ll_dm_sea_tcut(UseJack), damu_ss_dmcr_tcut_valence(UseJack), damu_ss_dmcr_tcut_sea(UseJack), damu_ss_dm_tcut(UseJack);
  distr_t_list damu_ll_dm_strange_sea_tcut(UseJack), damu_ss_dm_strange_sea_tcut(UseJack);
  distr_t_list damu_ss_OS_dm_tcut(UseJack), damu_ss_OS_dm_strange_sea_tcut(UseJack);

  //using condesate from WI
  distr_t_list damu_ll_dm_WI_sea_tcut(UseJack), damu_ss_dm_WI_tcut(UseJack);

  distr_t_list amu_ss_tcut(UseJack);
  distr_t_list amu_ss_old_tcut(UseJack);

 

    

  //define amu^HVP kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , T/2);
  

   
   
  for(int t=0;t<T/2;t++) {

    amu_ss_tcut.distr_list.push_back( (t==0)?(0.0*Get_id_jack_distr(Njacks)):( amu_ss_tcut[(t==0)?0:t-1] +  4.0*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_iso_distr[t]*Ker.distr_list[t] ));
    amu_ss_old_tcut.distr_list.push_back( (t==0)?(0.0*Get_id_jack_distr(Njacks)):( amu_ss_old_tcut[(t==0)?0:t-1] +  4.0*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_iso_old_distr[t]*Ker.distr_list[t] ));
    
    damu_ll_dmcr_tcut_valence.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ll_dmcr_tcut_valence[(t==0)?0:t-1] +  4.0*w(t,1)*dmcrit*(qu*qu + qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ll_dmcr_distr[t]*Ker.distr_list[t] ));
    damu_ss_dmcr_tcut_valence.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_dmcr_tcut_valence[(t==0)?0:t-1] +  4.0*w(t,1)*dmcrit*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_dmcr_distr[t]*Ker.distr_list[t] ));
    damu_ll_dmcr_tcut_sea.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ll_dmcr_tcut_sea[(t==0)?0:t-1] +  4.0*w(t,1)*dmcrit*(qu*qu + qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ll_dmcr_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_dmcr_tcut_sea.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_dmcr_tcut_sea[(t==0)?0:t-1] +  4.0*w(t,1)*dmcrit*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_dmcr_D_distr[t]*Ker.distr_list[t] ));
    damu_ll_dm_sea_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ll_dm_sea_tcut[(t==0)?0:t-1] +  4.0*dm*w(t,1)*(qu*qu + qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ll_dm_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_dm_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_dm_tcut[(t==0)?0:t-1] +  4.0*dm*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_dm_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_OS_dm_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_OS_dm_tcut[(t==0)?0:t-1] +  4.0*dm*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_ss_OS_dm_D_distr[t]*Ker.distr_list[t] ));
    //strange
    damu_ll_dm_strange_sea_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ll_dm_strange_sea_tcut[(t==0)?0:t-1] +  4.0*dm_s*w(t,1)*(qu*qu + qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ll_dm_strange_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_dm_strange_sea_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_dm_strange_sea_tcut[(t==0)?0:t-1] +  4.0*dm_s*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_dm_strange_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_OS_dm_strange_sea_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_OS_dm_strange_sea_tcut[(t==0)?0:t-1] +  4.0*dm_s*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_ss_OS_dm_strange_D_distr[t]*Ker.distr_list[t] ));

    //loop from WI
    damu_ll_dm_WI_sea_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ll_dm_WI_sea_tcut[(t==0)?0:t-1] +  4.0*dm*w(t,1)*(qu*qu + qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ll_dm_WI_D_distr[t]*Ker.distr_list[t] ));
    damu_ss_dm_WI_tcut.distr_list.push_back(  (t==0)?(0.0*Get_id_jack_distr(Njacks)):( damu_ss_dm_WI_tcut[(t==0)?0:t-1] +  4.0*dm*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_ss_dm_WI_D_distr[t]*Ker.distr_list[t] ));
  }

  distr_t_list damu_ll_dmcr_tcut = damu_ll_dmcr_tcut_valence + damu_ll_dmcr_tcut_sea;
  distr_t_list damu_ss_dmcr_tcut = damu_ss_dmcr_tcut_valence + damu_ss_dmcr_tcut_sea;
  
  Print_To_File({}, { damu_ll_dmcr_tcut.ave(), damu_ll_dmcr_tcut.err(),  damu_ll_dmcr_tcut_valence.ave(), damu_ll_dmcr_tcut_valence.err(), damu_ll_dmcr_tcut_sea.ave(), damu_ll_dmcr_tcut_sea.err(), damu_ll_dm_sea_tcut.ave(), damu_ll_dm_sea_tcut.err(), damu_ll_dm_strange_sea_tcut.ave(), damu_ll_dm_strange_sea_tcut.err(), damu_ll_dm_WI_sea_tcut.ave(), damu_ll_dm_WI_sea_tcut.err() } ,  "../data/sea_quark_effects/B64/damu_ll_tcut" , "", "# dmcr  dm");
  
  Print_To_File({}, { damu_ss_dmcr_tcut.ave(), damu_ss_dmcr_tcut.err(), damu_ss_dmcr_tcut_valence.ave(), damu_ss_dmcr_tcut_valence.err(), damu_ss_dmcr_tcut_sea.ave(), damu_ss_dmcr_tcut_sea.err(), damu_ss_dm_tcut.ave(), damu_ss_dm_tcut.err(), damu_ss_dm_strange_sea_tcut.ave(), damu_ss_dm_strange_sea_tcut.err(), damu_ss_dm_WI_tcut.ave(), damu_ss_dm_WI_tcut.err() } ,  "../data/sea_quark_effects/B64/damu_ss_tcut" , "", "# dmcr  dm");
  
  Print_To_File({}, {  damu_ss_OS_dm_tcut.ave(), damu_ss_OS_dm_tcut.err(), damu_ss_OS_dm_strange_sea_tcut.ave(), damu_ss_OS_dm_strange_sea_tcut.err() } ,  "../data/sea_quark_effects/B64/damu_ss_OS_tcut" , "", "# dmcr  dm");

  cout<<"amu(ss): "<<amu_ss_tcut.ave(50)<<" +- "<<amu_ss_tcut.err(50)<<endl;


  //HVP
  distr_t p2_mot= 2*SQRT_D( a_B*a_B*(0.140*0.140) + pow( 4*M_PI/64,2));
  distr_t amu_HVP_tm, amu_HVP_tm_old;
  int Tcut_opt_tm, Tcut_opt_tm_old;

  
  Bounding_HVP(amu_HVP_tm, Tcut_opt_tm,   1e10*( qu*qu + qd*qd   )*ZA_B*ZA_B*VKVK_ll_iso_distr, a_distr,"../data/HVP_blinded/B64/Bounding/tm" , p2_mot);
  Bounding_HVP(amu_HVP_tm_old, Tcut_opt_tm_old,   1e10*( qu*qu + qd*qd   )*ZA_B*ZA_B*VKVK_ll_iso_old_distr, a_distr,"../data/HVP_blinded/B64/Bounding/tm_old" , p2_mot);


  cout<<"amu(l,new): "<<amu_HVP_tm.ave()<<" "<<amu_HVP_tm.err()<<endl;
  cout<<"amu(l,old): "<<amu_HVP_tm_old.ave()<<" "<<amu_HVP_tm_old.err()<<endl;

  cout<<"amu(s,new): "<<amu_ss_tcut.ave(50)<<" +- "<<amu_ss_tcut.err(50)<<endl;
  cout<<"amu(s,old): "<<amu_ss_old_tcut.ave(50)<<" +- "<<amu_ss_old_tcut.err(50)<<endl;

  cout<<"amu(disco, new): "<<"-9.150955237099303 2.178609167198198"<<endl;
  cout<<"amu(disco, old): "<<"-11.51431935012435 4.691506508412073"<<endl;
  



  





  return;

}




void sea_quark_effects() {

  ofstream FAKE_L("../data/sea_quark_effects_NOVAL/table_Mpi.tex");
  FAKE_L.close();

  ofstream FAKE_K("../data/sea_quark_effects_NOVAL/table_MK.tex");
  FAKE_K.close();

  ofstream FAKE_Ds("../data/sea_quark_effects_NOVAL/table_MDs.tex");
  FAKE_Ds.close();

  
  ofstream der_FAKE_L("../data/sea_quark_effects_NOVAL/table_der_Mpi.tex");
  der_FAKE_L.close();

  ofstream der_FAKE_K("../data/sea_quark_effects_NOVAL/table_der_MK.tex");
  der_FAKE_K.close();

  ofstream der_FAKE_Ds("../data/sea_quark_effects_NOVAL/table_der_MDs.tex");
  der_FAKE_Ds.close();

  ofstream der_FAKE_FPI("../data/sea_quark_effects_NOVAL/table_der_Fpi.tex");
  der_FAKE_FPI.close();

  ofstream FAKE_FPI("../data/sea_quark_effects_NOVAL/table_Fpi.tex");
  FAKE_FPI.close();

  
  bool UseJack = true;
  int Njacks = 50;
  int Nboots=1000;
  
  double fm_to_inv_Gev= 1.0/0.197327;
  double qu = 2.0 / 3;
  double qd = -qu / 2;
  

  bool correlate_to_VKVK_strange=true;
  bool correlate_to_VKVK_charm=true;
  bool correlate_to_ls=true;
  bool correlate_to_sc=true;


  distr_t id_distr= Get_id_distr(UseJack?Njacks:Nboots, UseJack);

  double F= UseJack?(1.0/sqrt(Njacks-1.0)):1.0;
 


   //get ll,  ls and ss correlators

   bool Get_ASCII_l= false;
  
   if(Get_ASCII_l) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll");


      vector<string> Ens_T1({ "C.06.80", "B.72.64" , "D.54.96", "E.44.112"});
      vector<string> Ens_TT1({ "cC211a.06.80", "cB211b.072.64", "cD211a.054.96", "cE211a.044.112"});
    
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({"ll"});
	if(Ens_T1[it]=="E.44.112") channels = {"mix_l_l"};
	
		
	boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll/"+Ens_TT1[it]);
	
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_P5P5", "OS_P5P5", "OS_S0S0"});


	
   
          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"Ens: "<<Ens_T1[it]<<" (ll)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;

	string cc = channel;
	if(Ens_T1[it]=="E.44.112") cc= "ll";
	//read config file
	vector<string> C_list;
	ifstream CONFS_LIST;
	CONFS_LIST.open("../gm2_tau_rep_bin/"+Ens_T1[it]+"/confsListll.txt");
	while(!CONFS_LIST.eof()) {
	  string a;
	  CONFS_LIST>>a;
	  if(!CONFS_LIST.eof()) C_list.push_back(a);
	}

	if(Nconfs != (signed)C_list.size()) crash("Nconfs not matching in ll for ensemble: "+Ens_T1[it]+ "Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	CONFS_LIST.close();
	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL/ll/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+cc+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();
	  
	}
	
	fclose(stream);
	
	}
	
      }
    }
      
   }
  
	

 
   
  
  bool Get_ASCII_Mix_s= false;
  
  if(Get_ASCII_Mix_s) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s");


      vector<string> Ens_T1({"B.72.64" ,  "E.44.112"});
      vector<string> Ens_TT1({ "cB211b.072.64", "cE211a.044.112"});
    
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({ "mix_s1_s1", "mix_s2_s2"});
	
	
	boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]);
	
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_VKVK", "TM_P5P5", "OS_VKVK", "OS_P5P5", "TM_A0P5", "OS_A0P5"});


	
   
          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"Ens: "<<Ens_T1[it]<<" (Mix)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	//read config file
	vector<string> C_list;
	ifstream CONFS_LIST;
	CONFS_LIST.open("../gm2_tau_rep_bin/"+Ens_T1[it]+"/confsListMix");
	while(!CONFS_LIST.eof()) {
	  string a;
	  CONFS_LIST>>a;
	  if(!CONFS_LIST.eof()) C_list.push_back(a);
	}

	if(Nconfs != (signed)C_list.size()) crash("Nconfs not matching in Mix_s for ensemble: "+Ens_T1[it]+ " Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	    
	CONFS_LIST.close();
	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();
	  
	}
	
	fclose(stream);
	
	}
	
      }
    }
      
  }



  bool Get_ASCII_Mix_s_old= false;
  
  if(Get_ASCII_Mix_s_old) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s");


      vector<string> Ens_T1({"C.06.80", "D.54.96"});
      vector<string> Ens_TT1({ "cC211a.06.80", "cD211a.054.96"});
    
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({"ss1", "ss2"});
	
	
	boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]);
	
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_VKVK", "TM_P5P5", "OS_VKVK", "OS_P5P5", "TM_A0P5", "OS_A0P5"});


	
   
          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"Ens: "<<Ens_T1[it]<<" (Mix)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;

	string c, conf_file;
	if(channel=="ss1") { c="mix_s1_s1"; conf_file="confsListss1.txt"; }
	else if(channel=="ss2") { c="mix_s2_s2"; conf_file="confsListss2.txt";}
	
	//read config file
	vector<string> C_list;
	ifstream CONFS_LIST;
	CONFS_LIST.open("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+conf_file);
	while(!CONFS_LIST.eof()) {
	  string a;
	  CONFS_LIST>>a;
	  if(!CONFS_LIST.eof()) C_list.push_back(a);
	}

	if(Nconfs != (signed)C_list.size()) crash("Nconfs not matching in Mix_s for ensemble: "+Ens_T1[it]+ " Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	    
	CONFS_LIST.close();


	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+c+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();
	  
	}
	
	fclose(stream);
	
	}
	
      }
    }
      
  }

  
  //##########################################################################################



  bool Get_ASCII_new_format= false;
  
  if(Get_ASCII_new_format) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll_NF");


      vector<string> Ens_T1({"C.06.80", "D.54.96", "B.72.64", "E.44.112"});
      vector<string> Ens_TT1({ "cC211a.06.80", "cD211a.054.96", "cB211b.072.64", "cE211a.044.112"});

      //vector<string> Ens_T1({ "B.72.64"});
      //vector<string> Ens_TT1({  "cB211b.072.64"});
      
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({"gm2_run_l_l"});
	
	
	boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll_NF/"+Ens_TT1[it]);
	
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_P5P5", "TM_V0P5"});


	
   
          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin_new_format/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");

	BaseDataInfo Binfo;
	bin_read(Binfo,stream);
	size_t Nconfs=Binfo.nConfs;
	size_t T=Binfo.T;
	size_t Nhits=Binfo.nSources;

	string c = "ll";
	
       	cout<<"Ens: "<<Ens_T1[it]<<" (NF)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	cout<<"Version: "<<Binfo.version<<endl;
	cout<<"Git hash: "<<Binfo._gitHash<<endl;
	cout<<"sizeof(Binfo): "<<sizeof(Binfo)<<endl;

	//read confs list
	size_t confsListLength;
	bin_read(confsListLength,stream);
	cout<<"confListLength: "<<confsListLength<<endl;
	char confsList[confsListLength+1];
	fread(confsList,1,confsListLength,stream);
	confsList[confsListLength]='\0';
	cout<<confsList<<endl;
	size_t mesonsListLength;
	bin_read(mesonsListLength,stream);
	char mesonsListData[mesonsListLength];
	fread(mesonsListData,1,mesonsListLength,stream);
	//read config file
	vector<string> C_list(Nconfs,"");
	//assign configs
	int conf_csize= (confsListLength+1)/Nconfs;
	for(int iconf=0;iconf<Nconfs;iconf++) {
	  for(int a=0;a< conf_csize-1;a++) {
	    C_list[iconf] = C_list[iconf] + confsList[ iconf*conf_csize+a];
	  }
	}

	if(Nconfs != (signed)C_list.size()) crash("Nconfs not matching in NF for ensemble: "+Ens_T1[it]+ " Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	    
	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL/ll_NF/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL/ll_NF/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+c+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();
	  
	}
	
	fclose(stream);
	
	}
	
      }
    }
      
  }
  

  
 

  //##########################################################################################


  



  bool Get_ASCII_cc = false;
  if(Get_ASCII_cc) {

    vector<string> Ens_T1({ "C.06.80",  "D.54.96", "B.72.64", "E.44.112"});
    vector<string> Ens_TT1({ "cC211a.06.80", "cD211a.054.96",  "cB211b.072.64",  "cE211a.044.112" });


    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_c1_c1", "mix_c2_c2", "mix_s_c1"});
      if(Ens_T1[it] == "C.06.80" || Ens_T1[it] == "B.72.64")  channels={"mix_c1_c1", "mix_c2_c2", "mix_s1_c1"};

      boost::filesystem::create_directory("../sea_quark_effects_NOVAL_charm");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL_charm/cc");
           
      
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL_charm/cc/"+Ens_TT1[it]);
      
      //read binary
      vector<string> Corr_tags({"TM_P5P5", "TM_VKVK", "OS_VKVK"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../charm_E_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"Ens: "<<Ens_T1[it]<<" (cc)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	//read config file
	vector<string> C_list;
	string ch=channel;
	if(channel=="mix_s_c1") ch="mix_s1_c1";
	ifstream CONFS_LIST;
	CONFS_LIST.open("../charm_E_bin/"+Ens_T1[it]+"/confsList_lsc");
	while(!CONFS_LIST.eof()) {
	  string a;
	  CONFS_LIST>>a;
	  if(!CONFS_LIST.eof()) C_list.push_back(a);
	}

        if (Nconfs != (signed)C_list.size())
          crash("Nconfs not matching in cc for ensemble: "+Ens_T1[it]+ " Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	CONFS_LIST.close();
	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL_charm/cc/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL_charm/cc/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+ch+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();

	}

	fclose(stream);

	}
	
      }
    }
  }

  //################################################################################################################################################


  bool Get_ASCII_Mix_ls= false;
  
  if(Get_ASCII_Mix_ls) {
      //read binary files                                                                                                                                                                                                                                                          
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL_ls");
      boost::filesystem::create_directory("../sea_quark_effects_NOVAL_ls/Mix_s");

      vector<string> Ens_T1({ "C.06.80",  "D.54.96", "B.72.64", "E.44.112"});
      vector<string> Ens_TT1({ "cC211a.06.80", "cD211a.054.96",  "cB211b.072.64",  "cE211a.044.112" });

    
      
      for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	vector<string> channels({ "mix_l_s1"});
	
	
	boost::filesystem::create_directory("../sea_quark_effects_NOVAL_ls/Mix_s/"+Ens_TT1[it]);
	
	//read binary                                                                                                                                                                                                                                                              
	vector<string> Corr_tags({"TM_P5P5"});


	
   
          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"Ens: "<<Ens_T1[it]<<" (Mix)"<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	//read config file
	vector<string> C_list;
	ifstream CONFS_LIST;
	CONFS_LIST.open("../gm2_tau_rep_bin/"+Ens_T1[it]+"/confsListMix");
	while(!CONFS_LIST.eof()) {
	  string a;
	  CONFS_LIST>>a;
	  if(!CONFS_LIST.eof()) C_list.push_back(a);
	}

	if(Nconfs != (signed)C_list.size()) crash("Nconfs not matching in Mix_s for ensemble: "+Ens_T1[it]+ " Nconfs: "+to_string(Nconfs)+" expected: "+to_string( C_list.size()));
	    
	    
	CONFS_LIST.close();
	
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../sea_quark_effects_NOVAL_ls/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]);
	  ofstream PrintCorr("../sea_quark_effects_NOVAL_ls/Mix_s/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;
	  PrintCorr.close();
	  
	}
	
	fclose(stream);
	
	}
	
      }
    }
      
  }


  
  
  cout<<"Correlators read"<<endl;


  scale_setting_info SCALE_INFO = Get_scale_setting_info();

 
 

  auto Sort_light_confs = [](string A, string B) {

			   

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			    //A_bis=A;
			    //B_bis=B;

			     
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

   auto Sort_random = [](string A, string B) {

			   

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			    //A_bis=A;
			    //B_bis=B;

			     
			    string conf_num_A = A_bis.substr(0,4);
			    string conf_num_B = B_bis.substr(0,4);
							       
		      
			    string rA = A_bis.substr(A_bis.length()-2);
			    string rB = B_bis.substr(B_bis.length()-2);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(A_bis.substr(A_bis.length()-1));
			      int n2 = stoi(B_bis.substr(B_bis.length()-1));
			      if(rA == rB) {
			      if(rA=="r0" || rA=="r2") return conf_num_A < conf_num_B;
			      else if(rA=="r1" || rA=="r3") return conf_num_A > conf_num_B;
			      else crash("stream not recognized");
			      }
			      else return n1<n2;
			    }
			    return A_bis<B_bis;
  };

  

     
   
   
  //light
    data_t P5P5_TM_ll, P5P5_OS_ll;
    data_t V0P5_TM_ll;
  //strange-strange
    data_t P5P5_TM_ss1, VKVK_TM_ss1, VKVK_OS_ss1, P5P5_OS_ss1, A0P5_OS_ss1, A0P5_TM_ss1;
    data_t P5P5_TM_ss2, VKVK_TM_ss2, VKVK_OS_ss2;
    //light-strange
    data_t P5P5_TM_ls, P5P5_TM_sc;
    //charm-charm
    data_t P5P5_TM_cc1, VKVK_TM_cc1, VKVK_OS_cc1;
    data_t P5P5_TM_cc2, VKVK_TM_cc2, VKVK_OS_cc2;
    

    //read files

   
   
    P5P5_TM_ll.Read("../sea_quark_effects_NOVAL/ll_NF", "mes_contr_ll_TM_P5P5",  "P5P5", Sort_light_confs);


    

   
    //P5P5_OS_ll.Read("../sea_quark_effects_NOVAL/ll_NF", "mes_contr_ll_OS_P5P5",  "P5P5", Sort_light_confs);
    V0P5_TM_ll.Read("../sea_quark_effects_NOVAL/ll_NF", "mes_contr_ll_TM_V0P5",  "V0P5", Sort_light_confs);
   
    P5P5_TM_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_TM_P5P5",  "P5P5", Sort_light_confs);
    VKVK_TM_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_TM_VKVK",  "VKVK", Sort_light_confs);
    VKVK_OS_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_OS_VKVK",  "VKVK", Sort_light_confs);
    P5P5_OS_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_OS_P5P5",  "P5P5", Sort_light_confs);
    A0P5_TM_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_TM_A0P5",  "A0P5", Sort_light_confs);
    A0P5_OS_ss1.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s1_s1_OS_A0P5",  "A0P5", Sort_light_confs);
    
    P5P5_TM_ss2.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s2_s2_TM_P5P5",  "P5P5", Sort_light_confs);
    VKVK_TM_ss2.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s2_s2_TM_VKVK",  "VKVK", Sort_light_confs);
    VKVK_OS_ss2.Read("../sea_quark_effects_NOVAL/Mix_s", "mes_contr_mix_s2_s2_OS_VKVK",  "VKVK", Sort_light_confs);
    
    P5P5_TM_ls.Read("../sea_quark_effects_NOVAL_ls/Mix_s", "mes_contr_mix_l_s1_TM_P5P5",  "P5P5", Sort_light_confs);
    P5P5_TM_sc.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_s1_c1_TM_P5P5",  "P5P5", Sort_light_confs);
    
    P5P5_TM_cc1.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_c1_c1_TM_P5P5",  "P5P5", Sort_light_confs);
    VKVK_TM_cc1.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_c1_c1_TM_VKVK",  "VKVK", Sort_light_confs);
    VKVK_OS_cc1.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_c1_c1_OS_VKVK",  "VKVK", Sort_light_confs);
    
    P5P5_TM_cc2.Read("../sea_quark_effects_NOVAL_charm/cc" , "mes_contr_mix_c2_c2_TM_P5P5",  "P5P5", Sort_light_confs);
    VKVK_TM_cc2.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_c2_c2_TM_VKVK",  "VKVK", Sort_light_confs);
    VKVK_OS_cc2.Read("../sea_quark_effects_NOVAL_charm/cc", "mes_contr_mix_c2_c2_OS_VKVK",  "VKVK", Sort_light_confs);
    
    
    //read loops
    data_t loops_S_s, loops_S_c, loops_S_l, loops_P_l, loops_P_s;

    loops_S_l.Read("../sea_quark_effects_NOVAL/loop_ll", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_P_l.Read("../sea_quark_effects_NOVAL/loop_ll", "bub_P5_diff.txt",  "P5P5", Sort_light_confs);
    loops_P_s.Read("../sea_quark_effects_NOVAL/loop_ss", "bub_P5_diff.txt",  "P5P5", Sort_light_confs);
    loops_S_s.Read("../sea_quark_effects_NOVAL/loop_ss", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_S_c.Read("../sea_quark_effects_NOVAL/loop_cc", "bub_S.txt",  "S0S0", Sort_light_confs);


    //read loops_C
    data_t loops_C_S_s, loops_C_S_c, loops_C_S_l, loops_C_P_l;

    loops_C_S_l.Read("../sea_quark_effects_NOVAL_charm/loop_ll", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_C_P_l.Read("../sea_quark_effects_NOVAL_charm/loop_ll", "bub_P5_diff.txt",  "P5P5", Sort_light_confs);
    loops_C_S_s.Read("../sea_quark_effects_NOVAL_charm/loop_ss", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_C_S_c.Read("../sea_quark_effects_NOVAL_charm/loop_cc", "bub_S.txt",  "S0S0", Sort_light_confs);


    //read loops_ls
    data_t loops_ls_S_s, loops_ls_S_c, loops_ls_S_l, loops_ls_P_l;

    loops_ls_S_l.Read("../sea_quark_effects_NOVAL_ls/loop_ll", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_ls_P_l.Read("../sea_quark_effects_NOVAL_ls/loop_ll", "bub_P5_diff.txt",  "P5P5", Sort_light_confs);
    loops_ls_S_s.Read("../sea_quark_effects_NOVAL_ls/loop_ss", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_ls_S_c.Read("../sea_quark_effects_NOVAL_ls/loop_cc", "bub_S.txt",  "S0S0", Sort_light_confs);

    //read loops_sc
    data_t loops_sc_S_s, loops_sc_S_c, loops_sc_S_l, loops_sc_P_l;

    loops_sc_S_l.Read("../sea_quark_effects_NOVAL_charm/loop_ll", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_sc_P_l.Read("../sea_quark_effects_NOVAL_charm/loop_ll", "bub_P5_diff.txt",  "P5P5", Sort_light_confs);
    loops_sc_S_s.Read("../sea_quark_effects_NOVAL_charm/loop_ss", "bub_S.txt",  "S0S0", Sort_light_confs);
    loops_sc_S_c.Read("../sea_quark_effects_NOVAL_charm/loop_cc", "bub_S.txt",  "S0S0", Sort_light_confs);


    data_t V0P5_val, P5P5_val;

    V0P5_val.Read("../LIBE_data", "dmcr_plus_ll_muval_+7.2000e-04.dat",  "V0P5", Sort_light_confs);
    P5P5_val.Read("../LIBE_data", "dmcr_plus_ll_muval_+7.2000e-04.dat",  "P5P5", Sort_light_confs);
      
  

    int Nens= P5P5_TM_ll.size;

    distr_t_list damu_s_charm(UseJack);
    distr_t_list damu_W_s_charm(UseJack);
    distr_t_list damu_SD_s_charm(UseJack);
    distr_t_list Deltas(UseJack);
    vector<string> Ens_list;



    distr_t TOT_dmlc_for_E112(UseJack),  W_dmlc_for_E112(UseJack) ,  SD_dmlc_for_E112(UseJack);
    distr_t TOT_C_dmlc_for_E112(UseJack),  W_C_dmlc_for_E112(UseJack) ,  SD_C_dmlc_for_E112(UseJack);
    distr_t TOT_ZV_for_E112(UseJack), TOT_ZA_for_E112(UseJack);
    
    distr_t_list FPI_der_l(UseJack);
    distr_t_list FPI_der_s(UseJack);
    distr_t_list FPI_der_crit(UseJack);
    distr_t_list FPI_l(UseJack);
    distr_t_list FPI_s(UseJack);
    distr_t_list FPI_crit(UseJack);
    distr_t_list FPI_TRUE_list(UseJack);
    distr_t_list Daml_list(UseJack);
    distr_t_list Dams_list(UseJack);
    distr_t_list Damcr_list(UseJack);
    distr_t_list a_distr_list(UseJack);
    vector<string> EN_RED_list;
    
    for(int iens=0;iens<Nens;iens++) {

      string EN= P5P5_TM_ll.Tag[iens];

      int lens=-1; int sens=-1; int cens=-1;

      for(int b=0;b<(signed)SCALE_INFO.Ens_l.size();b++) { if(SCALE_INFO.Ens_l[b] == EN) lens=b;}
      for(int b=0;b<(signed)SCALE_INFO.Ens_c.size();b++) { if(SCALE_INFO.Ens_c[b] == EN) cens=b;}
      for(int b=0;b<(signed)SCALE_INFO.Ens.size();b++) { if(SCALE_INFO.Ens[b] == EN) sens=b;}

      if( (lens==-1) || (sens==-1) || (cens==-1)) crash("Cannot find ensemble: "+EN+" in SCALE_INFO");
   
      //int Njacks=  loops_P_l.col(0)[iens][0].size();

      if(P5P5_TM_ll.Tag[iens] != loops_S_c.Tag[iens]) crash("Ensembles do not match");
      if(VKVK_TM_ss1.Tag[iens] != P5P5_TM_ll.Tag[iens]) crash("Ensembles do not match");
      if(VKVK_TM_ss1.Tag[iens] != P5P5_TM_ls.Tag[iens]) crash("Ensembles do not match");
      if(VKVK_TM_ss1.Tag[iens] != P5P5_TM_sc.Tag[iens]) crash("Ensembles do not match");
      

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
      
      
      for(int ijack=0;ijack<(UseJack?Njacks:Nboots);ijack++) {
	
	a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(F)));
	a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(F)));
	a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(F)));
	a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(F)));
	a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(F)));
	a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(F)));
	ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(F));
	ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(F));
	ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(F));
	ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(F));
	ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(F));
	ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(F));
	ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(F));
	ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(F));
	ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err*(F));
	ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err*(F));
	
      }
      
      LatticeInfo L_info;
      string Ens = P5P5_TM_ll.Tag[iens];

      boost::filesystem::create_directory("../data/sea_quark_effects_NOVAL");
      boost::filesystem::create_directory("../data/sea_quark_effects_NOVAL/"+Ens);


      Ens_list.push_back(Ens);
      
      L_info.LatInfo_new_ens(Ens);
      
      //light and strange masses
      double aml= L_info.ml;
      double ams1= L_info.ms_L_new;
      double ams2= L_info.ms_M_new;
      //charm quark masses
      double amc1, amc2;
      double Zs_ov_Zp;
      if(Ens.substr(1,1)=="B") { amc1= 2.3000e-01 ; amc2= 2.4000e-01 ; Zs_ov_Zp = 1.2657137;}
      else if( Ens.substr(1,1)=="C")     { amc1= 0.18; amc2=0.19; Zs_ov_Zp = 1.214674 ;}
      else if(Ens.substr(1,1)=="D") { amc1= 0.15; amc2=0.16; Zs_ov_Zp = 1.1746661 ;} 
      else if(Ens.substr(1,1)=="E")  {amc1= 0.13; amc2=0.14; Zs_ov_Zp= 1.148525; }
      else crash("charm quark masses not found for Ensemble: "+Ens);
      //get lattice spacing and RCs
      distr_t a_distr(UseJack);
      distr_t ZV(UseJack), ZA(UseJack);
      if(Ens.substr(1,1)=="B") {a_distr=a_B; ZV = ZV_B; ZA = ZA_B;}
      else if(Ens.substr(1,1)=="C") {a_distr=a_C; ZV = ZV_C; ZA = ZA_C; }
      else if(Ens.substr(1,1)=="D") {a_distr=a_D; ZV = ZV_D; ZA = ZA_D; }
      else if(Ens.substr(1,1)=="E") {a_distr=a_E; ZV = ZV_E; ZA = ZA_E; }
      else crash("lattice spacing distribution for Ens: "+Ens+" not found");
      //lattice volume
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = P5P5_TM_ll.nrows[iens];
      int L= Corr.Nt/2;
      int T= Corr.Nt;
      double V4 = (double)pow(L,3)*T;
      double V= (double)pow(L,3);
      
      cout<<"####### Ensemble : "+Ens<<" ######"<<endl;
      

      distr_t_list P5P5_TM_ll_distr= Corr.corr_t( P5P5_TM_ll.col(0)[iens], "");
      //distr_t_list P5P5_OS_ll_distr= Corr.corr_t( P5P5_OS_ll.col(0)[iens], "");
      distr_t_list V0P5_TM_ll_distr= Corr.corr_t( V0P5_TM_ll.col(0)[iens], "");
      distr_t_list P5P5_TM_ss1_distr= Corr.corr_t( P5P5_TM_ss1.col(0)[iens], "");
      distr_t_list P5P5_OS_ss1_distr= Corr.corr_t( P5P5_OS_ss1.col(0)[iens], "");
      distr_t_list P5P5_TM_ss2_distr= Corr.corr_t( P5P5_TM_ss2.col(0)[iens], "");
      distr_t_list VKVK_TM_ss1_distr= Corr.corr_t( VKVK_TM_ss1.col(0)[iens], "");
      distr_t_list VKVK_TM_ss2_distr= Corr.corr_t( VKVK_TM_ss2.col(0)[iens], "");
      distr_t_list VKVK_OS_ss1_distr= Corr.corr_t( VKVK_OS_ss1.col(0)[iens], "");
      distr_t_list VKVK_OS_ss2_distr= Corr.corr_t( VKVK_OS_ss2.col(0)[iens], "");
      distr_t_list A0P5_TM_ss1_distr= Corr.corr_t( A0P5_TM_ss1.col(0)[iens], "");
      distr_t_list A0P5_OS_ss1_distr= Corr.corr_t( A0P5_OS_ss1.col(0)[iens], "");

      distr_t_list P5P5_TM_cc1_distr= Corr.corr_t( P5P5_TM_cc1.col(0)[iens], "");
      distr_t_list P5P5_TM_cc2_distr= Corr.corr_t( P5P5_TM_cc2.col(0)[iens], "");
      distr_t_list VKVK_TM_cc1_distr= Corr.corr_t( VKVK_TM_cc1.col(0)[iens], "");
      distr_t_list VKVK_TM_cc2_distr= Corr.corr_t( VKVK_TM_cc2.col(0)[iens], "");
      distr_t_list VKVK_OS_cc1_distr= Corr.corr_t( VKVK_OS_cc1.col(0)[iens], "");
      distr_t_list VKVK_OS_cc2_distr= Corr.corr_t( VKVK_OS_cc2.col(0)[iens], "");


      distr_t_list P5P5_TM_ls_distr=Corr.corr_t(P5P5_TM_ls.col(0)[iens], "");
      distr_t_list P5P5_TM_sc_distr=Corr.corr_t(P5P5_TM_sc.col(0)[iens], "");


      distr_t_list V0P5_val_distr(UseJack);
      distr_t_list P5P5_val_distr(UseJack);
      if(Ens=="cB211b.072.64") {
	Corr.Reflection_sign=-1;
	V0P5_val_distr = (2.0/V)*Corr.corr_t( V0P5_val.col(0)[0], "");
	Corr.Reflection_sign=1;
	P5P5_val_distr = (2.0/V)*Corr.corr_t( P5P5_val.col(0)[0], "");
      }
   
      
      
      //evaluate sea quark mass derivatives

      //define scalar and pseudoscalar bubble for each conf;

     int Nconfs=  loops_P_l.col(0)[iens][0].size();
     
     if(loops_P_l.Tag[iens] != loops_C_P_l.Tag[iens]) crash("Ens and Ens_C not match");

     int Nconfs_C = loops_C_P_l.col(0)[iens][0].size();

     int Nconfs_ls = loops_ls_P_l.col(0)[iens][0].size();
     int Nconfs_sc = loops_sc_P_l.col(0)[iens][0].size();

     
     //Vfloat Sl(Nconfs,0.0), Ss1(Nconfs,0.0), Ss2(Nconfs,0.0), Sc1(Nconfs,0.0), Sc2(Nconfs,0.0);
     
     Vfloat Sl(Nconfs,0.0), Ss(Nconfs,0.0), Sc(Nconfs,0.0), Pl(Nconfs,0.0), Ps(Nconfs, 0.0), Pl2(Nconfs,0.0);
     Vfloat Sl_WI(Nconfs,0.0);

     Vfloat Sl_C(Nconfs_C,0.0), Ss_C(Nconfs_C,0.0), Sc_C(Nconfs_C,0.0), Pl_C(Nconfs_C,0.0) ;

     Vfloat Sl_ls(Nconfs_ls,0.0), Ss_ls(Nconfs_ls,0.0), Sc_ls(Nconfs_ls,0.0), Pl_ls(Nconfs_ls,0.0) ;
     Vfloat Sl_sc(Nconfs_sc,0.0), Ss_sc(Nconfs_sc,0.0), Sc_sc(Nconfs_sc,0.0), Pl_sc(Nconfs_sc,0.0) ;

     double Dams= -0.1*ams1;
     double Daml = -0.1*aml ;
     double Damc= -0.1*amc1;
     //double Damcr = 0.05*aml/2.5;
     double amL1;

     string EN_RED;
     
     if( Ens =="cB211b.072.96")     {  amL1= 0.0006675; EN_RED ="B96";}
     else if(Ens =="cB211b.072.64") { amL1= 0.0006675;  EN_RED ="B64";}
     else if(Ens.substr(1,1)=="C")  { amL1= 0.000585;  EN_RED = "C80"; }
     else if(Ens.substr(1,1)=="D")  { amL1= 0.0004964;  EN_RED = "D96";}
     else if(Ens.substr(1,1)=="E")  { amL1= 0.000431;  EN_RED ="E112"; }
     else crash("Cannot recognize the ensemble: "+Ens+" in determining dmal ");

     EN_RED_list.push_back(EN_RED);
     
     Daml= amL1 - aml;
     
     if(Ens.substr(1,1)=="B") {
       Damc=  0.231567 - 0.22857550800 ;
       Dams=  0.0182782*0.0115*2;
     }
     else if(Ens.substr(1,1)=="C") {
       Damc= 0.198396 - 0.19479694565946 ;
       Dams= -0.0160609*0.0185*2;
     }
     else if(Ens.substr(1,1)=="D") {
       Damc= 0.164898 - 0.17241661418594 ;
       Dams= 0.012*0.013576;
     }
     else if(Ens.substr(1,1)=="E") {
       Damc= 0.141255 - 0.1426997693952;
       Dams= -0.0117933*0.0246;
     }
     else crash("Ens not found");

     Daml_list.distr_list.push_back( pow(10,3)*Daml/(a_distr));
     Dams_list.distr_list.push_back( pow(10,3)*Dams/(a_distr));
  

     a_distr_list.distr_list.push_back(a_distr);
         

     double Sl_ave=0; double Ss_ave=0.0; double  Sc_ave = 0.0; double Pl_ave=0.0; double Pl2_ave=0.0; double Ps_ave;

     double Sl_C_ave=0; double Ss_C_ave=0.0; double  Sc_C_ave = 0.0; double Pl_C_ave=0.0; double Pl2_C_ave=0.0;

     double Sl_ls_ave=0; double Ss_ls_ave=0.0; double  Sc_ls_ave = 0.0; double Pl_ls_ave=0.0; double Pl2_ls_ave=0.0;

     double Sl_sc_ave=0; double Ss_sc_ave=0.0; double  Sc_sc_ave = 0.0; double Pl_sc_ave=0.0; double Pl2_sc_ave=0.0;

     //Vfloat Sl_rew; Vfloat Ss_rew=0;Vfloat Sc_rew=0;Vfloat Pl_rew=0;

     double Sl_rew=0; double Ss_rew=0.0; double Sc_rew= 0.0; double Pl_rew=0.0;
     double ph_Sl=0; double ph_Ss=0; double ph_Sc=0; double ph_Pl=0.0;
     
     for(int iconf=0; iconf<Nconfs;iconf++)  {
       for(int t=0;t<T;t++) {
       	 Sl[iconf] += 2*loops_S_l.col(1)[iens][t][iconf];
	 Ss[iconf] += loops_S_s.col(1)[iens][t][iconf];
	 Sc[iconf] += loops_S_c.col(1)[iens][t][iconf];
	 Pl[iconf] += 2*loops_P_l.col(0)[iens][t][iconf];
	 Ps[iconf] += loops_P_s.col(0)[iens][t][iconf];
       }

       Sl_ave += Sl[iconf]/Nconfs;
       Ss_ave += Ss[iconf]/Nconfs;
       Sc_ave += Sc[iconf]/Nconfs;
       Pl_ave += Pl[iconf]/Nconfs;
       Ps_ave += Ps[iconf]/Nconfs;
       Pl2_ave += Pl[iconf]/Nconfs;
     }

     for(int iconf=0; iconf<Nconfs_C;iconf++)  {
       for(int t=0;t<T;t++) {
       
	 Sl_C[iconf] += 2*loops_C_S_l.col(1)[iens][t][iconf];
	 Ss_C[iconf] += loops_C_S_s.col(1)[iens][t][iconf];
	 Sc_C[iconf] += loops_C_S_c.col(1)[iens][t][iconf];
	 Pl_C[iconf] += 2*loops_C_P_l.col(0)[iens][t][iconf];
       }

       Sl_C_ave += Sl_C[iconf]/Nconfs_C;
       Ss_C_ave += Ss_C[iconf]/Nconfs_C;
       Sc_C_ave += Sc_C[iconf]/Nconfs_C;
       Pl_C_ave += Pl_C[iconf]/Nconfs_C;
      
     }

     for(int iconf=0; iconf<Nconfs_ls;iconf++)  {
       for(int t=0;t<T;t++) {
       
	 Sl_ls[iconf] += 2*loops_ls_S_l.col(1)[iens][t][iconf];
	 Ss_ls[iconf] += loops_ls_S_s.col(1)[iens][t][iconf];
	 Sc_ls[iconf] += loops_ls_S_c.col(1)[iens][t][iconf];
	 Pl_ls[iconf] += 2*loops_ls_P_l.col(0)[iens][t][iconf];
       }

       Sl_ls_ave += Sl_ls[iconf]/Nconfs_ls;
       Ss_ls_ave += Ss_ls[iconf]/Nconfs_ls;
       Sc_ls_ave += Sc_ls[iconf]/Nconfs_ls;
       Pl_ls_ave += Pl_ls[iconf]/Nconfs_ls;
      
     }

     for(int iconf=0; iconf<Nconfs_sc;iconf++)  {
       for(int t=0;t<T;t++) {
       
	 Sl_sc[iconf] += 2*loops_sc_S_l.col(1)[iens][t][iconf];
	 Ss_sc[iconf] += loops_sc_S_s.col(1)[iens][t][iconf];
	 Sc_sc[iconf] += loops_sc_S_c.col(1)[iens][t][iconf];
	 Pl_sc[iconf] += 2*loops_sc_P_l.col(0)[iens][t][iconf];
       }

       Sl_sc_ave += Sl_sc[iconf]/Nconfs_sc;
       Ss_sc_ave += Ss_sc[iconf]/Nconfs_sc;
       Sc_sc_ave += Sc_sc[iconf]/Nconfs_sc;
       Pl_sc_ave += Pl_sc[iconf]/Nconfs_sc;
      
     }


     

     /*

     for(int iconf=0;iconf<Nconfs;iconf++) {

       ph_Sl += exp(Daml*(Sl[iconf] - Sl_ave));
       ph_Ss += exp(Dams*(Ss[iconf] - Ss_ave));
       ph_Sc += exp(Damc*(Sc[iconf] - Sc_ave));
     }

     double checksum_l=0; double checksum_s=0; double checksum_c=0;

     for(int iconf=0; iconf<Nconfs;iconf++) {
       //Nconfs^eff/Nconfs

       double wl=exp(Daml*(Sl[iconf] - Sl_ave))/ph_Sl;
       double ws=exp(Dams*(Ss[iconf] - Ss_ave))/ph_Ss;
       double wc=exp(Damc*(Sc[iconf] - Sc_ave))/ph_Sc;

       checksum_l+=wl; checksum_s+=ws; checksum_c += wc;
      
       Sl_rew += pow(wl,2);
       Ss_rew += pow(ws,2);
       Sc_rew += pow(wc,2);
       //Sc_rew += 1/exp(Damcr*(Sc[iconf] - Pl_ave));

       cout<<"iconf "<<iconf<<" wl: "<< wl << " ws: "<<ws<<" wc: "<<wc<<endl;

       if(wc > 0.90) cout<<"WARNING : wc: "<<wc<<endl;
       if(wc > 0.95) cout<<"SUPER-WARNING : wc: "<<wc<<endl;
       if(wc > 0.99) cout<<"HYPER-WARNING : wc: "<<wc<<endl;

     }

     Sl_rew = 1.0/Sl_rew;
     Ss_rew = 1.0/Ss_rew;
     Sc_rew = 1.0/Sc_rew;


     cout<<"Nconfs: "<<Nconfs<<endl;
     cout<<"Eff_confs(l): "<<Sl_rew<<" checksum: "<<checksum_l<<endl;
     cout<<"Eff_confs(s): "<<Ss_rew<<" checksum: "<<checksum_s<<endl;
     cout<<"Eff_confs(c): "<<Sc_rew<<" checksum: "<<checksum_c<<endl;
     
     */

     

     double var_Sc=0; double var_Ss=0; double var_Sl=0.0;

     for(int iconf=0;iconf<Nconfs;iconf++) {

       var_Sc += pow( Sc[iconf] - Sc_ave,2)/(Nconfs-1.0);
       var_Ss += pow( Ss[iconf] - Ss_ave,2)/(Nconfs-1.0);
       var_Sl += pow( Sl[iconf] - Sl_ave,2)/(Nconfs-1.0);

            
       //Sl[iconf] = exp( ( Sl[iconf]-Sl_ave)*Daml);
       //Ss[iconf] = exp( (Ss[iconf]-Ss_ave)*Dams);
       //Sc[iconf] = exp( (Sc[iconf]-Sc_ave)*Damc);
       //Pl2[iconf] = exp( (Pl[iconf]-Pl_ave)*Damcr)*exp( 0.5*(Damcr/aml)*Damcr*(Pl2[iconf]-Pl2_ave));
       //Pl[iconf] = exp( (Pl[iconf]-Pl_ave)*Damcr);
       Pl2[iconf] = Pl[iconf] ;
     
     
     }


     
       
     cout<<"condensate computed! "<<endl;
     cout<<"sqrt{var} of <Sl>: "<<sqrt(var_Sl)<<", ave: "<<Sl_ave<<endl;
     cout<<"sqrt{var} of <Ss>: "<<sqrt(var_Ss)<<", ave: "<<Ss_ave<<endl;
     cout<<"sqrt{var} of <Sc>: "<<sqrt(var_Sc)<<", ave: "<<Sc_ave<<endl;

     VVfloat P5P5_TM_dmcrit(T), V0P5_TM_dmcrit(T), V0P5_TM_dmcrit_s(T);
     VVfloat P5P5_TM_dms(T), P5P5_TM_dmc(T), P5P5_TM_dml(T);
     VVfloat P5P5_TM_ls_dms(T), P5P5_TM_ls_dmcrit(T), P5P5_TM_ls_dml(T);
     VVfloat P5P5_TM_sc_dms(T), P5P5_TM_sc_dmcrit(T), P5P5_TM_sc_dml(T);

     //FOR ZV AND ZA
     VVfloat P5P5_TM_ss1_dmcrit(T), P5P5_TM_ss1_dms(T), P5P5_TM_ss1_dml(T);
     VVfloat P5P5_OS_ss1_dmcrit(T), P5P5_OS_ss1_dms(T), P5P5_OS_ss1_dml(T);
     VVfloat A0P5_TM_ss1_dmcrit(T), A0P5_TM_ss1_dms(T), A0P5_TM_ss1_dml(T);
     VVfloat A0P5_OS_ss1_dmcrit(T), A0P5_OS_ss1_dms(T), A0P5_OS_ss1_dml(T);
     
     VVfloat VKVK_TM_ss1_dms(T), VKVK_OS_ss1_dms(T), VKVK_TM_ss2_dms(T), VKVK_OS_ss2_dms(T);
     VVfloat VKVK_TM_ss1_dmc(T), VKVK_TM_ss1_dml(T), VKVK_OS_ss1_dmc(T), VKVK_OS_ss1_dml(T);

     VVfloat VKVK_TM_cc1_dmc(T), VKVK_OS_cc1_dmc(T), VKVK_TM_cc2_dmc(T), VKVK_OS_cc2_dmc(T);
     VVfloat VKVK_TM_cc1_dms(T), VKVK_TM_cc1_dml(T);


     VVfloat VKVK_TM_cc1_dmcrit(T), VKVK_OS_cc1_dmcrit(T), VKVK_TM_ss1_dmcrit(T), VKVK_OS_ss1_dmcrit(T);
     VVfloat VKVK_TM_cc1_d2mcrit(T), VKVK_OS_cc1_d2mcrit(T), VKVK_TM_ss1_d2mcrit(T), VKVK_OS_ss1_d2mcrit(T);
        
     
     auto SEA_FUNC_OR= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function SEA_FUNC expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[2];};

     auto SEA_FUNC= [&](const Vfloat& par) { if((signed)par.size() != 3) crash("Lambda function SEA_FUNC expects par[3], but par["+to_string((signed)par.size())+"] provided"); return par[0]/par[2] -par[1];};
     
     
  for(int t=0; t<T;t++) {
    for(int iconf=0;iconf<Nconfs;iconf++) {

      if(correlate_to_VKVK_strange) {
	VKVK_TM_ss1_dms[t].push_back( Ss[iconf]*VKVK_TM_ss1.col(0)[iens][t][iconf]);
	VKVK_OS_ss1_dms[t].push_back( Ss[iconf]*VKVK_OS_ss1.col(0)[iens][t][iconf]);
	VKVK_TM_ss2_dms[t].push_back( Ss[iconf]*VKVK_TM_ss2.col(0)[iens][t][iconf]);
	VKVK_OS_ss2_dms[t].push_back( Ss[iconf]*VKVK_OS_ss2.col(0)[iens][t][iconf]);
	VKVK_TM_ss1_dmc[t].push_back( Sc[iconf]*VKVK_TM_ss1.col(0)[iens][t][iconf]);
	VKVK_TM_ss1_dml[t].push_back( Sl[iconf]*VKVK_TM_ss1.col(0)[iens][t][iconf]);
	VKVK_OS_ss1_dmc[t].push_back( Sc[iconf]*VKVK_OS_ss1.col(0)[iens][t][iconf]);
	VKVK_OS_ss1_dml[t].push_back( Sl[iconf]*VKVK_OS_ss1.col(0)[iens][t][iconf]);
	VKVK_TM_ss1_dmcrit[t].push_back( Pl[iconf]*VKVK_TM_ss1.col(0)[iens][t][iconf]);
	VKVK_OS_ss1_dmcrit[t].push_back( Pl[iconf]*VKVK_OS_ss1.col(0)[iens][t][iconf]);
	P5P5_TM_dmcrit[t].push_back( Pl[iconf]*P5P5_TM_ll.col(0)[iens][t][iconf]);
	V0P5_TM_dmcrit[t].push_back( Pl[iconf]*V0P5_TM_ll.col(0)[iens][t][iconf]);
	V0P5_TM_dmcrit_s[t].push_back( Ps[iconf]*V0P5_TM_ll.col(0)[iens][t][iconf]);
	P5P5_TM_dms[t].push_back( Ss[iconf]*P5P5_TM_ll.col(0)[iens][t][iconf]);
	
	P5P5_TM_dml[t].push_back( Sl[iconf]*P5P5_TM_ll.col(0)[iens][t][iconf]);
	P5P5_TM_dmc[t].push_back( Sc[iconf]*P5P5_TM_ll.col(0)[iens][t][iconf]);


	//for ZV and ZA
	P5P5_TM_ss1_dmcrit[t].push_back( Pl[iconf]*P5P5_TM_ss1.col(0)[iens][t][iconf]);
	P5P5_TM_ss1_dms[t].push_back( Ss[iconf]*P5P5_TM_ss1.col(0)[iens][t][iconf]);
	P5P5_TM_ss1_dml[t].push_back( Sl[iconf]*P5P5_TM_ss1.col(0)[iens][t][iconf]);

	
	P5P5_OS_ss1_dmcrit[t].push_back( Pl[iconf]*P5P5_OS_ss1.col(0)[iens][t][iconf]);
	P5P5_OS_ss1_dms[t].push_back( Ss[iconf]*P5P5_OS_ss1.col(0)[iens][t][iconf]);
	P5P5_OS_ss1_dml[t].push_back( Sl[iconf]*P5P5_OS_ss1.col(0)[iens][t][iconf]);

	A0P5_OS_ss1_dmcrit[t].push_back( Pl[iconf]*A0P5_OS_ss1.col(0)[iens][t][iconf]);
	A0P5_OS_ss1_dms[t].push_back( Ss[iconf]*A0P5_OS_ss1.col(0)[iens][t][iconf]);
	A0P5_OS_ss1_dml[t].push_back( Sl[iconf]*A0P5_OS_ss1.col(0)[iens][t][iconf]);

	A0P5_TM_ss1_dmcrit[t].push_back( Pl[iconf]*A0P5_TM_ss1.col(0)[iens][t][iconf]);
	A0P5_TM_ss1_dms[t].push_back( Ss[iconf]*A0P5_TM_ss1.col(0)[iens][t][iconf]);
	A0P5_TM_ss1_dml[t].push_back( Sl[iconf]*A0P5_TM_ss1.col(0)[iens][t][iconf]);

	
      }
    }
  
    for(int iconf=0;iconf<Nconfs_C;iconf++) {
      if(correlate_to_VKVK_charm) {
	VKVK_TM_cc1_dmc[t].push_back( Sc_C[iconf]*VKVK_TM_cc1.col(0)[iens][t][iconf]);
	VKVK_OS_cc1_dmc[t].push_back( Sc_C[iconf]*VKVK_OS_cc1.col(0)[iens][t][iconf]);
	VKVK_TM_cc2_dmc[t].push_back( Sc_C[iconf]*VKVK_TM_cc2.col(0)[iens][t][iconf]);
	VKVK_OS_cc2_dmc[t].push_back( Sc_C[iconf]*VKVK_OS_cc2.col(0)[iens][t][iconf]);
	VKVK_TM_cc1_dms[t].push_back( Ss_C[iconf]*VKVK_TM_cc1.col(0)[iens][t][iconf]);
	VKVK_TM_cc1_dml[t].push_back( Sl_C[iconf]*VKVK_TM_cc1.col(0)[iens][t][iconf]);
	//dmcrit corrections
	VKVK_TM_cc1_dmcrit[t].push_back( Pl_C[iconf]*VKVK_TM_cc1.col(0)[iens][t][iconf]);
	VKVK_OS_cc1_dmcrit[t].push_back( Pl_C[iconf]*VKVK_OS_cc1.col(0)[iens][t][iconf]);
      }
 
      
    }

     for(int iconf=0;iconf<Nconfs_ls;iconf++) {
      if(correlate_to_ls) {
	P5P5_TM_ls_dmcrit[t].push_back( Pl_ls[iconf]*P5P5_TM_ls.col(0)[iens][t][iconf]);
	P5P5_TM_ls_dms[t].push_back( Ss_ls[iconf]*P5P5_TM_ls.col(0)[iens][t][iconf]);
	P5P5_TM_ls_dml[t].push_back( Sl_ls[iconf]*P5P5_TM_ls.col(0)[iens][t][iconf]);
      }
     }

     for(int iconf=0;iconf<Nconfs_sc;iconf++) {
      if(correlate_to_sc) {
	P5P5_TM_sc_dmcrit[t].push_back( Pl_sc[iconf]*P5P5_TM_sc.col(0)[iens][t][iconf]);
	P5P5_TM_sc_dms[t].push_back( Ss_sc[iconf]*P5P5_TM_sc.col(0)[iens][t][iconf]);
	P5P5_TM_sc_dml[t].push_back( Sl_sc[iconf]*P5P5_TM_sc.col(0)[iens][t][iconf]);
      }
     }


    
  }


   cout<<"Condensate correlated to Obs!"<<endl;

  
  
  Jackknife J(10000,Njacks);
  Jackknife J_C(10000,Njacks);
  Jackknife J_ls(10000,Njacks);
  Jackknife J_sc(10000,Njacks);

  Bootstrap B(Nboots,988453453,Nconfs);
  Bootstrap B_C(Nboots,988453453,Nconfs_C);
  Bootstrap B_ls(Nboots,988453453,Nconfs_ls);
  Bootstrap B_sc(Nboots,988453453,Nconfs_sc);

  
  distr_t_list VKVK_TM_ss1_dms_distr(UseJack), VKVK_OS_ss1_dms_distr(UseJack), VKVK_TM_ss2_dms_distr(UseJack), VKVK_OS_ss2_dms_distr(UseJack);
  distr_t_list VKVK_TM_ss1_dmc_distr(UseJack), VKVK_TM_ss1_dml_distr(UseJack), VKVK_OS_ss1_dmc_distr(UseJack), VKVK_OS_ss1_dml_distr(UseJack);

  distr_t_list VKVK_TM_cc1_dmc_distr(UseJack), VKVK_OS_cc1_dmc_distr(UseJack), VKVK_TM_cc2_dmc_distr(UseJack), VKVK_OS_cc2_dmc_distr(UseJack);
  distr_t_list VKVK_TM_cc1_dms_distr(UseJack), VKVK_TM_cc1_dml_distr(UseJack);


  distr_t_list VKVK_TM_ss1_dmcrit_distr(UseJack), VKVK_OS_ss1_dmcrit_distr(UseJack);
  distr_t_list VKVK_TM_cc1_dmcrit_distr(UseJack), VKVK_OS_cc1_dmcrit_distr(UseJack);

 
  distr_t_list P5P5_TM_dmcrit_distr(UseJack), V0P5_TM_dmcrit_distr(UseJack), V0P5_TM_dmcrit_s_distr(UseJack);
  distr_t_list P5P5_TM_dms_distr(UseJack), P5P5_TM_dmc_distr(UseJack), P5P5_TM_dml_distr(UseJack);

  distr_t_list P5P5_TM_ls_dms_distr(UseJack), P5P5_TM_ls_dmcrit_distr(UseJack), P5P5_TM_ls_dml_distr(UseJack);
  distr_t_list P5P5_TM_sc_dms_distr(UseJack), P5P5_TM_sc_dmcrit_distr(UseJack), P5P5_TM_sc_dml_distr(UseJack);

  //FOR ZV AND ZA
  distr_t_list P5P5_TM_ss1_dmcrit_distr(UseJack), P5P5_TM_ss1_dms_distr(UseJack), P5P5_TM_ss1_dml_distr(UseJack);
  distr_t_list P5P5_OS_ss1_dmcrit_distr(UseJack), P5P5_OS_ss1_dms_distr(UseJack), P5P5_OS_ss1_dml_distr(UseJack);
  distr_t_list A0P5_TM_ss1_dmcrit_distr(UseJack), A0P5_TM_ss1_dms_distr(UseJack), A0P5_TM_ss1_dml_distr(UseJack);
  distr_t_list A0P5_OS_ss1_dmcrit_distr(UseJack), A0P5_OS_ss1_dms_distr(UseJack), A0P5_OS_ss1_dml_distr(UseJack);
  

  
  for(int t=0;t<T;t++) {

    if(UseJack) {    //jackknife
      if(correlate_to_VKVK_strange) {
	VKVK_TM_ss1_dms_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_TM_ss1_dms[t], VKVK_TM_ss1.col(0)[iens][t], Ss ));
	VKVK_OS_ss1_dms_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_OS_ss1_dms[t], VKVK_OS_ss1.col(0)[iens][t], Ss ));
	VKVK_TM_ss2_dms_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_TM_ss2_dms[t], VKVK_TM_ss2.col(0)[iens][t], Ss ));
	VKVK_OS_ss2_dms_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_OS_ss2_dms[t], VKVK_OS_ss2.col(0)[iens][t], Ss ));
	VKVK_TM_ss1_dmc_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_TM_ss1_dmc[t], VKVK_TM_ss1.col(0)[iens][t], Sc ));
	VKVK_TM_ss1_dml_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_TM_ss1_dml[t], VKVK_TM_ss1.col(0)[iens][t], Sl ));
	VKVK_OS_ss1_dmc_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_OS_ss1_dmc[t], VKVK_OS_ss1.col(0)[iens][t], Sc ));
	VKVK_OS_ss1_dml_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_OS_ss1_dml[t], VKVK_OS_ss1.col(0)[iens][t], Sl ));
	VKVK_TM_ss1_dmcrit_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_TM_ss1_dmcrit[t], VKVK_TM_ss1.col(0)[iens][t], Pl ));
	VKVK_OS_ss1_dmcrit_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, VKVK_OS_ss1_dmcrit[t], VKVK_OS_ss1.col(0)[iens][t], Pl ));
	P5P5_TM_dmcrit_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_dmcrit[t], P5P5_TM_ll.col(0)[iens][t], Pl));
	V0P5_TM_dmcrit_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, V0P5_TM_dmcrit[t], V0P5_TM_ll.col(0)[iens][t], Pl));
	V0P5_TM_dmcrit_s_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, V0P5_TM_dmcrit_s[t], V0P5_TM_ll.col(0)[iens][t], Ps));
	P5P5_TM_dms_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_dms[t], P5P5_TM_ll.col(0)[iens][t], Ss));
	P5P5_TM_dmc_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_dmc[t], P5P5_TM_ll.col(0)[iens][t], Sc));
	P5P5_TM_dml_distr.distr_list.push_back( J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_dml[t], P5P5_TM_ll.col(0)[iens][t], Sl));

	//FOR ZV AND ZA
	P5P5_TM_ss1_dmcrit_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ss1_dmcrit[t], P5P5_TM_ss1.col(0)[iens][t], Pl )       );
	P5P5_OS_ss1_dmcrit_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_OS_ss1_dmcrit[t], P5P5_OS_ss1.col(0)[iens][t], Pl )       );
	A0P5_TM_ss1_dmcrit_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_TM_ss1_dmcrit[t], A0P5_TM_ss1.col(0)[iens][t], Pl )       );
	A0P5_OS_ss1_dmcrit_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_OS_ss1_dmcrit[t], A0P5_OS_ss1.col(0)[iens][t], Pl )       );

	P5P5_TM_ss1_dms_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ss1_dms[t], P5P5_TM_ss1.col(0)[iens][t], Ss )       );
	P5P5_OS_ss1_dms_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_OS_ss1_dms[t], P5P5_OS_ss1.col(0)[iens][t], Ss )       );
	A0P5_TM_ss1_dms_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_TM_ss1_dms[t], A0P5_TM_ss1.col(0)[iens][t], Ss )       );
	A0P5_OS_ss1_dms_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_OS_ss1_dms[t], A0P5_OS_ss1.col(0)[iens][t], Ss )       );

	P5P5_TM_ss1_dml_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ss1_dml[t], P5P5_TM_ss1.col(0)[iens][t], Sl )       );
	P5P5_OS_ss1_dml_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, P5P5_OS_ss1_dml[t], P5P5_OS_ss1.col(0)[iens][t], Sl )       );
	A0P5_TM_ss1_dml_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_TM_ss1_dml[t], A0P5_TM_ss1.col(0)[iens][t], Sl )       );
	A0P5_OS_ss1_dml_distr.distr_list.push_back(  J.DoJack( SEA_FUNC_OR, 3, A0P5_OS_ss1_dml[t], A0P5_OS_ss1.col(0)[iens][t], Sl )       );
	
      }
      
      if(correlate_to_VKVK_charm) {
	VKVK_TM_cc1_dmc_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_TM_cc1_dmc[t], VKVK_TM_cc1.col(0)[iens][t], Sc_C ));
	VKVK_OS_cc1_dmc_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_OS_cc1_dmc[t], VKVK_OS_cc1.col(0)[iens][t], Sc_C ));
	VKVK_TM_cc2_dmc_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_TM_cc2_dmc[t], VKVK_TM_cc2.col(0)[iens][t], Sc_C ));
	VKVK_OS_cc2_dmc_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_OS_cc2_dmc[t], VKVK_OS_cc2.col(0)[iens][t], Sc_C ));
	VKVK_TM_cc1_dms_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_TM_cc1_dms[t], VKVK_TM_cc1.col(0)[iens][t], Ss_C ));
	VKVK_TM_cc1_dml_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_TM_cc1_dml[t], VKVK_TM_cc1.col(0)[iens][t], Sl_C ));
	VKVK_TM_cc1_dmcrit_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_TM_cc1_dmcrit[t], VKVK_TM_cc1.col(0)[iens][t], Pl_C ));
	VKVK_OS_cc1_dmcrit_distr.distr_list.push_back( J_C.DoJack( SEA_FUNC_OR, 3, VKVK_OS_cc1_dmcrit[t], VKVK_OS_cc1.col(0)[iens][t], Pl_C ));
      }
      if(correlate_to_ls) {
	P5P5_TM_ls_dmcrit_distr.distr_list.push_back( J_ls.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ls_dmcrit[t], P5P5_TM_ls.col(0)[iens][t], Pl_ls));
	P5P5_TM_ls_dms_distr.distr_list.push_back( J_ls.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ls_dms[t], P5P5_TM_ls.col(0)[iens][t], Ss_ls));
	P5P5_TM_ls_dml_distr.distr_list.push_back( J_ls.DoJack( SEA_FUNC_OR, 3, P5P5_TM_ls_dml[t], P5P5_TM_ls.col(0)[iens][t], Sl_ls));
      }
      if(correlate_to_sc) {
	P5P5_TM_sc_dmcrit_distr.distr_list.push_back( J_sc.DoJack( SEA_FUNC_OR, 3, P5P5_TM_sc_dmcrit[t], P5P5_TM_sc.col(0)[iens][t], Pl_sc));
	P5P5_TM_sc_dms_distr.distr_list.push_back( J_sc.DoJack( SEA_FUNC_OR, 3, P5P5_TM_sc_dms[t], P5P5_TM_sc.col(0)[iens][t], Ss_sc));
	P5P5_TM_sc_dml_distr.distr_list.push_back( J_sc.DoJack( SEA_FUNC_OR, 3, P5P5_TM_sc_dml[t], P5P5_TM_sc.col(0)[iens][t], Sl_sc));
      }

      
    }
    else {    //bootstrap
        if(correlate_to_VKVK_strange) {
	VKVK_TM_ss1_dms_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_ss1_dms[t], VKVK_TM_ss1.col(0)[iens][t], Ss ));
	VKVK_OS_ss1_dms_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_ss1_dms[t], VKVK_OS_ss1.col(0)[iens][t], Ss ));
	VKVK_TM_ss2_dms_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_ss2_dms[t], VKVK_TM_ss2.col(0)[iens][t], Ss ));
	VKVK_OS_ss2_dms_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_ss2_dms[t], VKVK_OS_ss2.col(0)[iens][t], Ss ));
	VKVK_TM_ss1_dmc_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_ss1_dmc[t], VKVK_TM_ss1.col(0)[iens][t], Sc ));
	VKVK_TM_ss1_dml_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_ss1_dml[t], VKVK_TM_ss1.col(0)[iens][t], Sl ));
	VKVK_OS_ss1_dmc_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_ss1_dmc[t], VKVK_OS_ss1.col(0)[iens][t], Sc ));
	VKVK_OS_ss1_dml_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_ss1_dml[t], VKVK_OS_ss1.col(0)[iens][t], Sl ));
	VKVK_TM_ss1_dmcrit_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_ss1_dmcrit[t], VKVK_TM_ss1.col(0)[iens][t], Pl ));
	VKVK_OS_ss1_dmcrit_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_ss1_dmcrit[t], VKVK_OS_ss1.col(0)[iens][t], Pl ));
	P5P5_TM_dmcrit_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_dmcrit[t], P5P5_TM_ll.col(0)[iens][t], Pl));
	V0P5_TM_dmcrit_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, V0P5_TM_dmcrit[t], V0P5_TM_ll.col(0)[iens][t], Pl));
	P5P5_TM_dms_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_dms[t], P5P5_TM_ll.col(0)[iens][t], Ss));
	P5P5_TM_dmc_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_dmc[t], P5P5_TM_ll.col(0)[iens][t], Sc));
	P5P5_TM_dml_distr.distr_list.push_back( B.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_dml[t], P5P5_TM_ll.col(0)[iens][t], Sl));
      }
      
      if(correlate_to_VKVK_charm) {
	VKVK_TM_cc1_dmc_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_cc1_dmc[t], VKVK_TM_cc1.col(0)[iens][t], Sc_C ));
	VKVK_OS_cc1_dmc_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_cc1_dmc[t], VKVK_OS_cc1.col(0)[iens][t], Sc_C ));
	VKVK_TM_cc2_dmc_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_cc2_dmc[t], VKVK_TM_cc2.col(0)[iens][t], Sc_C ));
	VKVK_OS_cc2_dmc_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_cc2_dmc[t], VKVK_OS_cc2.col(0)[iens][t], Sc_C ));
	VKVK_TM_cc1_dms_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_cc1_dms[t], VKVK_TM_cc1.col(0)[iens][t], Ss_C ));
	VKVK_TM_cc1_dml_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_cc1_dml[t], VKVK_TM_cc1.col(0)[iens][t], Sl_C ));
	VKVK_TM_cc1_dmcrit_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_TM_cc1_dmcrit[t], VKVK_TM_cc1.col(0)[iens][t], Pl_C ));
	VKVK_OS_cc1_dmcrit_distr.distr_list.push_back( B_C.DoBoot( SEA_FUNC_OR, 3, VKVK_OS_cc1_dmcrit[t], VKVK_OS_cc1.col(0)[iens][t], Pl_C ));
      }

      if(correlate_to_ls) {
	P5P5_TM_ls_dmcrit_distr.distr_list.push_back( B_ls.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_ls_dmcrit[t], P5P5_TM_ls.col(0)[iens][t], Pl_ls));
	P5P5_TM_ls_dms_distr.distr_list.push_back( B_ls.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_ls_dms[t], P5P5_TM_ls.col(0)[iens][t], Ss_ls));
	P5P5_TM_ls_dml_distr.distr_list.push_back( B_ls.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_ls_dml[t], P5P5_TM_ls.col(0)[iens][t], Sl_ls));
      }
      if(correlate_to_sc) {
	P5P5_TM_sc_dmcrit_distr.distr_list.push_back( B_sc.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_sc_dmcrit[t], P5P5_TM_sc.col(0)[iens][t], Pl_sc));
	P5P5_TM_sc_dms_distr.distr_list.push_back( B_sc.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_sc_dms[t], P5P5_TM_sc.col(0)[iens][t], Ss_sc));
	P5P5_TM_sc_dml_distr.distr_list.push_back( B_sc.DoBoot( SEA_FUNC_OR, 3, P5P5_TM_sc_dml[t], P5P5_TM_sc.col(0)[iens][t], Sl_sc));
      }
      
      
    }
  }

  //symmetrize sea-quark correlators
  for(int t=0;t<=T/2;t++) {
    distr_t x(UseJack);


    if(correlate_to_VKVK_strange) {
      x= 0.5*(VKVK_TM_ss1_dms_distr[t] + VKVK_TM_ss1_dms_distr[(T-t)%T]);  VKVK_TM_ss1_dms_distr.distr_list[t] = x; VKVK_TM_ss1_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_OS_ss1_dms_distr[t] + VKVK_OS_ss1_dms_distr[(T-t)%T]);  VKVK_OS_ss1_dms_distr.distr_list[t] = x; VKVK_OS_ss1_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_ss2_dms_distr[t] + VKVK_TM_ss2_dms_distr[(T-t)%T]);  VKVK_TM_ss2_dms_distr.distr_list[t] = x; VKVK_TM_ss2_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_OS_ss2_dms_distr[t] + VKVK_OS_ss2_dms_distr[(T-t)%T]);  VKVK_OS_ss2_dms_distr.distr_list[t] = x; VKVK_OS_ss2_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_ss1_dmc_distr[t] + VKVK_TM_ss1_dmc_distr[(T-t)%T]);  VKVK_TM_ss1_dmc_distr.distr_list[t] = x; VKVK_TM_ss1_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_ss1_dml_distr[t] + VKVK_TM_ss1_dml_distr[(T-t)%T]);  VKVK_TM_ss1_dml_distr.distr_list[t] = x; VKVK_TM_ss1_dml_distr.distr_list[(T-t)%T] = x;

      x= 0.5*(VKVK_OS_ss1_dmc_distr[t] + VKVK_OS_ss1_dmc_distr[(T-t)%T]);  VKVK_OS_ss1_dmc_distr.distr_list[t] = x; VKVK_OS_ss1_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_OS_ss1_dml_distr[t] + VKVK_OS_ss1_dml_distr[(T-t)%T]);  VKVK_OS_ss1_dml_distr.distr_list[t] = x; VKVK_OS_ss1_dml_distr.distr_list[(T-t)%T] = x;

      //critical mass corrections
      x= 0.5*(VKVK_TM_ss1_dmcrit_distr[t] + VKVK_TM_ss1_dmcrit_distr[(T-t)%T]);  VKVK_TM_ss1_dmcrit_distr.distr_list[t] = x; VKVK_TM_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(VKVK_OS_ss1_dmcrit_distr[t] + VKVK_OS_ss1_dmcrit_distr[(T-t)%T]);  VKVK_OS_ss1_dmcrit_distr.distr_list[t] = x; VKVK_OS_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
   
      x= 0.5*(P5P5_TM_dmcrit_distr[t] + P5P5_TM_dmcrit_distr[(T-t)%T]); P5P5_TM_dmcrit_distr.distr_list[t] = x; P5P5_TM_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(V0P5_TM_dmcrit_distr[t] + V0P5_TM_dmcrit_distr[(T-t)%T]); V0P5_TM_dmcrit_distr.distr_list[t] = x; V0P5_TM_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(V0P5_TM_dmcrit_s_distr[t] + V0P5_TM_dmcrit_s_distr[(T-t)%T]); V0P5_TM_dmcrit_s_distr.distr_list[t] = x; V0P5_TM_dmcrit_s_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_dml_distr[t] + P5P5_TM_dml_distr[(T-t)%T]); P5P5_TM_dml_distr.distr_list[t] = x; P5P5_TM_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_dms_distr[t] + P5P5_TM_dms_distr[(T-t)%T]); P5P5_TM_dms_distr.distr_list[t] = x; P5P5_TM_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_dmc_distr[t] + P5P5_TM_dmc_distr[(T-t)%T]); P5P5_TM_dmc_distr.distr_list[t] = x; P5P5_TM_dmc_distr.distr_list[(T-t)%T] = x;


      //for ZV AND ZA
      x= 0.5*(P5P5_TM_ss1_dmcrit_distr[t] + P5P5_TM_ss1_dmcrit_distr[(T-t)%T]); P5P5_TM_ss1_dmcrit_distr.distr_list[t] = x; P5P5_TM_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_OS_ss1_dmcrit_distr[t] + P5P5_OS_ss1_dmcrit_distr[(T-t)%T]); P5P5_OS_ss1_dmcrit_distr.distr_list[t] = x; P5P5_OS_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_OS_ss1_dmcrit_distr[t] + A0P5_OS_ss1_dmcrit_distr[(T-t)%T]); A0P5_OS_ss1_dmcrit_distr.distr_list[t] = x; A0P5_OS_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_TM_ss1_dmcrit_distr[t] + A0P5_TM_ss1_dmcrit_distr[(T-t)%T]); A0P5_TM_ss1_dmcrit_distr.distr_list[t] = x; A0P5_TM_ss1_dmcrit_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(P5P5_TM_ss1_dms_distr[t] + P5P5_TM_ss1_dms_distr[(T-t)%T]); P5P5_TM_ss1_dms_distr.distr_list[t] = x; P5P5_TM_ss1_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_OS_ss1_dms_distr[t] + P5P5_OS_ss1_dms_distr[(T-t)%T]); P5P5_OS_ss1_dms_distr.distr_list[t] = x; P5P5_OS_ss1_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_OS_ss1_dms_distr[t] + A0P5_OS_ss1_dms_distr[(T-t)%T]); A0P5_OS_ss1_dms_distr.distr_list[t] = x; A0P5_OS_ss1_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_TM_ss1_dms_distr[t] + A0P5_TM_ss1_dms_distr[(T-t)%T]); A0P5_TM_ss1_dms_distr.distr_list[t] = x; A0P5_TM_ss1_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(P5P5_TM_ss1_dml_distr[t] + P5P5_TM_ss1_dml_distr[(T-t)%T]); P5P5_TM_ss1_dml_distr.distr_list[t] = x; P5P5_TM_ss1_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_OS_ss1_dml_distr[t] + P5P5_OS_ss1_dml_distr[(T-t)%T]); P5P5_OS_ss1_dml_distr.distr_list[t] = x; P5P5_OS_ss1_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_OS_ss1_dml_distr[t] + A0P5_OS_ss1_dml_distr[(T-t)%T]); A0P5_OS_ss1_dml_distr.distr_list[t] = x; A0P5_OS_ss1_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(A0P5_TM_ss1_dml_distr[t] + A0P5_TM_ss1_dml_distr[(T-t)%T]); A0P5_TM_ss1_dml_distr.distr_list[t] = x; A0P5_TM_ss1_dml_distr.distr_list[(T-t)%T] = x;
     
      
    }
    
    
    //charm
    if(correlate_to_VKVK_charm) {
      x= 0.5*(VKVK_TM_cc1_dmc_distr[t] + VKVK_TM_cc1_dmc_distr[(T-t)%T]);  VKVK_TM_cc1_dmc_distr.distr_list[t] = x; VKVK_TM_cc1_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_OS_cc1_dmc_distr[t] + VKVK_OS_cc1_dmc_distr[(T-t)%T]);  VKVK_OS_cc1_dmc_distr.distr_list[t] = x; VKVK_OS_cc1_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_cc2_dmc_distr[t] + VKVK_TM_cc2_dmc_distr[(T-t)%T]);  VKVK_TM_cc2_dmc_distr.distr_list[t] = x; VKVK_TM_cc2_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_OS_cc2_dmc_distr[t] + VKVK_OS_cc2_dmc_distr[(T-t)%T]);  VKVK_OS_cc2_dmc_distr.distr_list[t] = x; VKVK_OS_cc2_dmc_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_cc1_dms_distr[t] + VKVK_TM_cc1_dms_distr[(T-t)%T]);  VKVK_TM_cc1_dms_distr.distr_list[t] = x; VKVK_TM_cc1_dms_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_cc1_dml_distr[t] + VKVK_TM_cc1_dml_distr[(T-t)%T]);  VKVK_TM_cc1_dml_distr.distr_list[t] = x; VKVK_TM_cc1_dml_distr.distr_list[(T-t)%T] = x;
      
      x= 0.5*(VKVK_TM_cc1_dmcrit_distr[t] + VKVK_TM_cc1_dmcrit_distr[(T-t)%T]);  VKVK_TM_cc1_dmcrit_distr.distr_list[t] = x; VKVK_TM_cc1_dmcrit_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(VKVK_OS_cc1_dmcrit_distr[t] + VKVK_OS_cc1_dmcrit_distr[(T-t)%T]);  VKVK_OS_cc1_dmcrit_distr.distr_list[t] = x; VKVK_OS_cc1_dmcrit_distr.distr_list[(T-t)%T] = x;
      
   
    }

    if(correlate_to_ls) {
      x= 0.5*(P5P5_TM_ls_dml_distr[t] + P5P5_TM_ls_dml_distr[(T-t)%T]); P5P5_TM_ls_dml_distr.distr_list[t] = x; P5P5_TM_ls_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_ls_dms_distr[t] + P5P5_TM_ls_dms_distr[(T-t)%T]); P5P5_TM_ls_dms_distr.distr_list[t] = x; P5P5_TM_ls_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_ls_dmcrit_distr[t] + P5P5_TM_ls_dmcrit_distr[(T-t)%T]); P5P5_TM_ls_dmcrit_distr.distr_list[t] = x; P5P5_TM_ls_dmcrit_distr.distr_list[(T-t)%T] = x;
    }

    if(correlate_to_sc) {
      x= 0.5*(P5P5_TM_sc_dml_distr[t] + P5P5_TM_sc_dml_distr[(T-t)%T]); P5P5_TM_sc_dml_distr.distr_list[t] = x; P5P5_TM_sc_dml_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_sc_dms_distr[t] + P5P5_TM_sc_dms_distr[(T-t)%T]); P5P5_TM_sc_dms_distr.distr_list[t] = x; P5P5_TM_sc_dms_distr.distr_list[(T-t)%T] = x;
      x= 0.5*(P5P5_TM_sc_dmcrit_distr[t] + P5P5_TM_sc_dmcrit_distr[(T-t)%T]); P5P5_TM_sc_dmcrit_distr.distr_list[t] = x; P5P5_TM_sc_dmcrit_distr.distr_list[(T-t)%T] = x;
    }
    
    
    
  }



  cout<<"<OS> - <O><S> computed!"<<endl;
  
 

  //print corrections to dMpi and dfpi
  distr_t_list dmPi_l = Corr.effective_slope_t( P5P5_TM_dml_distr, P5P5_TM_ll_distr, "");
  distr_t_list dmPi_s = Corr.effective_slope_t( P5P5_TM_dms_distr, P5P5_TM_ll_distr, "");
  distr_t_list dmPi_c = Corr.effective_slope_t( P5P5_TM_dmc_distr, P5P5_TM_ll_distr, "");
  distr_t_list dmPi_crit = Corr.effective_slope_t( P5P5_TM_dmcrit_distr, P5P5_TM_ll_distr, "");

  distr_t_list dfPi_l = Corr.decay_constant_rel_slope_t( P5P5_TM_dml_distr, P5P5_TM_ll_distr, "");
  distr_t_list dfPi_s = Corr.decay_constant_rel_slope_t( P5P5_TM_dms_distr, P5P5_TM_ll_distr, "");
  distr_t_list dfPi_c = Corr.decay_constant_rel_slope_t( P5P5_TM_dmc_distr, P5P5_TM_ll_distr, "");
  distr_t_list dfPi_crit = Corr.decay_constant_rel_slope_t( P5P5_TM_dmcrit_distr, P5P5_TM_ll_distr, "");


  distr_t_list dmK_l= Corr.effective_slope_t( P5P5_TM_ls_dml_distr, P5P5_TM_ls_distr, "");
  distr_t_list dmK_s= Corr.effective_slope_t( P5P5_TM_ls_dms_distr, P5P5_TM_ls_distr, "");
  distr_t_list dmK_crit= Corr.effective_slope_t( P5P5_TM_ls_dmcrit_distr, P5P5_TM_ls_distr, "");

  distr_t_list dmDs_l= Corr.effective_slope_t( P5P5_TM_sc_dml_distr, P5P5_TM_sc_distr, "");
  distr_t_list dmDs_s= Corr.effective_slope_t( P5P5_TM_sc_dms_distr, P5P5_TM_sc_distr, "");
  distr_t_list dmDs_crit= Corr.effective_slope_t( P5P5_TM_sc_dmcrit_distr, P5P5_TM_sc_distr, "");

 

  distr_t_list damu_TM_ss1_dms(UseJack), damu_OS_ss1_dms(UseJack), damu_TM_ss2_dms(UseJack), damu_OS_ss2_dms(UseJack);
  distr_t_list damu_TM_ss1_dmc(UseJack), damu_TM_ss1_dml(UseJack), damu_OS_ss1_dmc(UseJack), damu_OS_ss1_dml(UseJack);

  distr_t_list damu_TM_cc1_dmc(UseJack), damu_OS_cc1_dmc(UseJack), damu_TM_cc2_dmc(UseJack), damu_OS_cc2_dmc(UseJack);
  distr_t_list damu_TM_cc1_dms(UseJack), damu_TM_cc1_dml(UseJack);

  //critical mass corrections
  distr_t_list damu_TM_ss1_dmcrit(UseJack), damu_OS_ss1_dmcrit(UseJack);
  distr_t_list damu_TM_cc1_dmcrit(UseJack), damu_OS_cc1_dmcrit(UseJack);





 

 
  distr_t_list amu_TM_ss1(UseJack);
  
  //get HVP corrections

  

  //get mpcac and dmpcac
  //mPCAC
  distr_t_list mPCAC_distr=  0.5*distr_t_list::derivative(V0P5_TM_ll_distr,0)/(aml*P5P5_TM_ll_distr);
 
  //derivative  dmPCAC/dmcrit
  distr_t_list dmPCAC_distr= -0.5*(distr_t_list::derivative(V0P5_TM_dmcrit_distr,0)/P5P5_TM_ll_distr   -0.0*distr_t_list::derivative(V0P5_TM_ll_distr,0)*P5P5_TM_dmcrit_distr/(P5P5_TM_ll_distr*P5P5_TM_ll_distr));

  distr_t_list dmPCAC_distr_s= -0.5*(distr_t_list::derivative(V0P5_TM_dmcrit_s_distr,0)/P5P5_TM_ll_distr);

 

  Print_To_File({}, {mPCAC_distr.ave(), mPCAC_distr.err(), (dmPCAC_distr).ave(), (dmPCAC_distr).err(), (dmPCAC_distr_s.ave()), (dmPCAC_distr_s.err())}, "../data/sea_quark_effects_NOVAL/"+Ens+"/mpcac", "", "");
  
  Corr.Tmin=20; //20
  Corr.Tmax=35; //40
  
  distr_t mPCAC = Corr.Fit_distr(mPCAC_distr);
  distr_t dmPCAC = Corr.Fit_distr(dmPCAC_distr); 
  cout<<"mpcac/mu: "<<mPCAC.ave()<<" +- "<<mPCAC.err()<<endl;
  cout<<"dmPCAC/dmcrit: "<<dmPCAC.ave()<<" +- "<<dmPCAC.err()<<" REN: "<<(dmPCAC*ZA*Zs_ov_Zp).ave()<<" +- "<<(dmPCAC*ZA*Zs_ov_Zp).err()<<endl;

  distr_t Damcr= aml*mPCAC/(1.1*dmPCAC);

  distr_t Damcr_new= Damcr.ave() + (Damcr-Damcr.ave())*sqrt( 1 + pow(0.1*Damcr.ave()/Damcr.err(),2));

  double kappa=0;
  if(Ens.substr(1,1)=="B") kappa=0.1394267;
  else if(Ens.substr(1,1)=="C") kappa=0.1387529;
  else if(Ens.substr(1,1)=="D") kappa=0.1379735;
  else if(Ens.substr(1,1)=="E") kappa=0.1374129;
  else crash("Ens: "+Ens+" not found");

  double m0_sim= (1.0/(2.0*kappa)) - 4.0;
  


  if(Ens.substr(1,1)=="E" ) {
    cout<<"Entering the loop: "<<"Ens: "<<Ens<<endl;
    mPCAC= -1.0*mPCAC;
    Damcr_new = aml*mPCAC/(1.1*1.1);
    Damcr_new= Damcr_new.ave() + (Damcr_new - Damcr_new.ave())*sqrt( 1.0 + pow(0.4*Damcr_new.ave()/Damcr_new.err(),2));
  }
  cout.precision(8);
  cout<<Ens<<": Damcr(nosyst): "<<Damcr.ave()<<" "<<Damcr.err()<<endl;
  cout<<Ens<<": Damcr: "<<Damcr_new.ave()<<" "<<Damcr_new.err()<<endl;
  cout<<Ens<<": m_cr^sim: "<<m0_sim<<endl;
  cout<<Ens<<": m_cr^iso(sum): "<<(m0_sim+Damcr_new.ave())<<" "<<Damcr_new.err()<<endl;
  cout<<Ens<<"; mcr^iso(diff): "<<(m0_sim-Damcr_new.ave())<<" "<<Damcr_new.err()<<endl;

  //ZV
  distr_t_list dZV_dmcrit = P5P5_TM_ss1_dmcrit_distr/P5P5_TM_ss1_distr   - distr_t_list::derivative( A0P5_TM_ss1_dmcrit_distr,0)/distr_t_list::derivative(A0P5_TM_ss1_distr,0);
  distr_t_list dZV_dms =  P5P5_TM_ss1_dms_distr/P5P5_TM_ss1_distr   - distr_t_list::derivative( A0P5_TM_ss1_dms_distr,0)/distr_t_list::derivative(A0P5_TM_ss1_distr,0);
  distr_t_list dZV_dml =  P5P5_TM_ss1_dml_distr/P5P5_TM_ss1_distr   - distr_t_list::derivative( A0P5_TM_ss1_dml_distr,0)/distr_t_list::derivative(A0P5_TM_ss1_distr,0);

  Print_To_File({}, {(a_distr*dZV_dmcrit).ave(), (a_distr*dZV_dmcrit).err(), (a_distr*dZV_dms).ave(), (a_distr*dZV_dms).err(), (a_distr*dZV_dml).ave(), (a_distr*dZV_dml).err() }, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dZV", "", "");
  Print_To_File({}, {(Damcr*dZV_dmcrit).ave(), (Damcr*dZV_dmcrit).err(), (Dams*dZV_dms).ave(), (Dams*dZV_dms).err(), (Daml*dZV_dml).ave(), (Daml*dZV_dml).err() }, "../data/sea_quark_effects_NOVAL/"+Ens+"/dZV", "", "");


  //ZA
  int Tmin_old= Corr.Tmin; int Tmax_old=Corr.Tmax;

   if(Ens.substr(1,1)=="B") { 
      if(Ens=="cB211b.072.64") {
	Corr.Tmin= 40; //40;
	Corr.Tmax= 57; //57;
      }
      else { 
	Corr.Tmin=40;
	Corr.Tmax=87;
      }
    }
    else if(Ens.substr(1,1)=="C") { 
      if(Ens=="cC211a.06.80")   {
	Corr.Tmin=44;
	Corr.Tmax=70;
      }
      else { 
	Corr.Tmin=40; Corr.Tmax=80;
      }
    }
    else if(Ens.substr(1,1)=="D") {Corr.Tmin=48; Corr.Tmax=80; }
    else if(Ens.substr(1,1)=="E") {Corr.Tmin=62; Corr.Tmax=100;}
    else crash("Ensemble not found");

 
   distr_t Meta_TM= Corr.Fit_distr( Corr.effective_mass_t(P5P5_TM_ss1_distr,""));
   distr_t Meta_OS= Corr.Fit_distr( Corr.effective_mass_t(P5P5_OS_ss1_distr, ""));
   
   distr_t_list dMeta_TM_dmcrit =-1.0*Corr.effective_slope_t(P5P5_TM_ss1_dmcrit_distr, P5P5_TM_ss1_distr, "");
   distr_t_list dMeta_TM_dms =-1.0*Corr.effective_slope_t(P5P5_TM_ss1_dms_distr, P5P5_TM_ss1_distr, "");
   distr_t_list dMeta_TM_dml =-1.0*Corr.effective_slope_t(P5P5_TM_ss1_dml_distr, P5P5_TM_ss1_distr, "");
  
   distr_t_list dMeta_OS_dmcrit =-1.0*Corr.effective_slope_t(P5P5_OS_ss1_dmcrit_distr, P5P5_OS_ss1_distr, "");
   distr_t_list dMeta_OS_dms =-1.0*Corr.effective_slope_t(P5P5_OS_ss1_dms_distr, P5P5_OS_ss1_distr, "");
   distr_t_list dMeta_OS_dml =-1.0*Corr.effective_slope_t(P5P5_OS_ss1_dml_distr, P5P5_OS_ss1_distr, "");

   Corr.Tmin=20;

   Corr.Tmax=40;
   
   distr_t_list dZA_dmcrit= P5P5_OS_ss1_dmcrit_distr/P5P5_OS_ss1_distr - distr_t_list::derivative( A0P5_OS_ss1_dmcrit_distr,0)/distr_t_list::derivative(A0P5_OS_ss1_distr,0) +
     dMeta_OS_dmcrit/Meta_OS -dMeta_TM_dmcrit/Meta_TM + dMeta_OS_dmcrit/TANH_D(Meta_OS) - dMeta_TM_dmcrit/TANH_D(Meta_TM) + Corr.effective_ampl_rel_slope_t(P5P5_TM_ss1_dmcrit_distr, P5P5_TM_ss1_distr, "") - Corr.effective_ampl_rel_slope_t(P5P5_OS_ss1_dmcrit_distr, P5P5_OS_ss1_distr,"");
   
   distr_t_list dZA_dms= P5P5_OS_ss1_dms_distr/P5P5_OS_ss1_distr - distr_t_list::derivative( A0P5_OS_ss1_dms_distr,0)/distr_t_list::derivative(A0P5_OS_ss1_distr,0) +
     dMeta_OS_dms/Meta_OS -dMeta_TM_dms/Meta_TM + dMeta_OS_dms/TANH_D(Meta_OS) - dMeta_TM_dms/TANH_D(Meta_TM) + Corr.effective_ampl_rel_slope_t(P5P5_TM_ss1_dms_distr, P5P5_TM_ss1_distr,"") - Corr.effective_ampl_rel_slope_t(P5P5_OS_ss1_dms_distr, P5P5_OS_ss1_distr, "");
   
   distr_t_list dZA_dml= P5P5_OS_ss1_dml_distr/P5P5_OS_ss1_distr - distr_t_list::derivative( A0P5_OS_ss1_dml_distr,0)/distr_t_list::derivative(A0P5_OS_ss1_distr,0) +
     dMeta_OS_dml/Meta_OS -dMeta_TM_dml/Meta_TM + dMeta_OS_dml/TANH_D(Meta_OS) - dMeta_TM_dml/TANH_D(Meta_TM) + Corr.effective_ampl_rel_slope_t(P5P5_TM_ss1_dml_distr, P5P5_TM_ss1_distr,"") - Corr.effective_ampl_rel_slope_t(P5P5_OS_ss1_dml_distr, P5P5_OS_ss1_distr,"");
   
   Print_To_File({}, {(a_distr*dZA_dmcrit).ave(), (a_distr*dZA_dmcrit).err(), (a_distr*dZA_dms).ave(), (a_distr*dZA_dms).err(), (a_distr*dZA_dml).ave(), (a_distr*dZA_dml).err() }, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dZA", "", "");
   Print_To_File({}, {(Damcr*dZA_dmcrit).ave(), (Damcr*dZA_dmcrit).err(), (Dams*dZA_dms).ave(), (Dams*dZA_dms).err(), (Daml*dZA_dml).ave(), (Daml*dZA_dml).err() }, "../data/sea_quark_effects_NOVAL/"+Ens+"/dZA", "", "");
   
   
   
   if(Ens.substr(1,1) == "C") {
     Corr.Tmin=9;
     Corr.Tmax=21;
   }

   
   if(Ens.substr(1,1) != "E") {
   
   distr_t dZA = Corr.Fit_distr(Daml*dZA_dml) + Corr.Fit_distr(Dams*dZA_dms) + Corr.Fit_distr(Damcr*dZA_dmcrit);
   dZA = dZA.ave() + (dZA-dZA.ave())*sqrt( 1 + pow( (Corr.Fit_distr(Dams*dZA_dms)*Damc/(12.0*Dams*dZA.err())).ave(),2));
   cout<<"dZA("<<Ens.substr(1,1)<<"): "<<1+dZA.ave()<<" "<<dZA.err()<<endl;
   distr_t dZV = Corr.Fit_distr(Daml*dZV_dml) + Corr.Fit_distr(Dams*dZV_dms) + Corr.Fit_distr(Damcr*dZV_dmcrit);
   dZV = dZV.ave() + (dZV-dZV.ave())*sqrt( 1 + pow( ( Corr.Fit_distr(Dams*dZV_dms)*Damc/(12.0*Dams*dZV.err())).ave(),2));
   cout<<"dZV("<<Ens.substr(1,1)<<"): "<<1+dZV.ave()<<" "<<dZV.err()<<endl;

   if(Ens.substr(1,1) == "C") {
     TOT_ZA_for_E112= Corr.Fit_distr((Damcr*dZA_dmcrit)) + Corr.Fit_distr(Daml*dZA_dml);
     TOT_ZV_for_E112= Corr.Fit_distr((Damcr*dZV_dmcrit)) + Corr.Fit_distr(Daml*dZV_dml);
   }
   
   }
   else {
     distr_t dZA = Corr.Fit_distr(Dams*dZA_dms);
     dZA = dZA.ave() + (dZA-dZA.ave())*sqrt( 1 + pow( (Corr.Fit_distr(Dams*dZA_dms)*Damc/(12.0*Dams*dZA.err())).ave(),2)  + pow(2e-5/dZA.err(),2));
     cout<<"dZA("<<Ens.substr(1,1)<<"): "<<1+dZA.ave()<<" "<<dZA.err()<<endl;
     distr_t dZV = Corr.Fit_distr(Dams*dZV_dms);
     dZV = dZV.ave() + (dZV-dZV.ave())*sqrt( 1 + pow( (Corr.Fit_distr(Dams*dZV_dms)*Damc/(12.0*Dams*dZV.err())).ave(),2) +  pow(2e-5/dZV.err(),2));
     cout<<"dZV("<<Ens.substr(1,1)<<"): "<<1+dZV.ave()<<" "<<dZV.err()<<endl;
     
   }
   
   



  Corr.Tmin=Tmin_old; Corr.Tmax=Tmax_old;

 


  if(Ens=="cB211b.072.64") {
    distr_t_list dmPCAC_val_distr = -0.5*(distr_t_list::derivative(-1.0*V0P5_val_distr,0)/P5P5_TM_ll_distr)  ;
    Print_To_File({}, {dmPCAC_val_distr.ave(), dmPCAC_val_distr.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/mpcac_val", "", "");
    distr_t_list Mpi_dmcr_der_val=  Corr.effective_slope_t( P5P5_val_distr, P5P5_TM_ll_distr, "");
    Corr.Tmin=20;
    Corr.Tmax=35;
    distr_t_list fpi_dmcr_der_val =  Corr.decay_constant_rel_slope_t( P5P5_val_distr, P5P5_TM_ll_distr, "");

    Print_To_File( {}, {(-1.0*Damcr*Mpi_dmcr_der_val/a_distr).ave(), (-1.0*Damcr*Mpi_dmcr_der_val/a_distr).err(), (Damcr*fpi_dmcr_der_val).ave(), (Damcr*fpi_dmcr_der_val).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_crit_val", "", "");
    
  }


  Damcr_list.distr_list.push_back( pow(10,3)*Damcr/(a_distr));


  distr_t DMK_l= Corr.Fit_distr(-1.0*Daml*dmK_l);
  distr_t DMK_crit= 0.5*Corr.Fit_distr( -1.0*Damcr*dmK_crit);
  distr_t DMK_s = Corr.Fit_distr(-1.0*Dams*dmK_s);

  distr_t DMDs_l = Corr.Fit_distr( -1.0*Daml*dmDs_l);
  distr_t DMDs_crit= Corr.Fit_distr(-1.0*Damcr*dmDs_crit);
  distr_t DMDs_s = Corr.Fit_distr( -1.0*Dams*dmDs_s);

  distr_t DMP_l = Corr.Fit_distr(-1.0*Daml*dmPi_l);
  distr_t DMP_crit = 0.5*Corr.Fit_distr( 1.0*Damcr*dmPi_crit);
  distr_t DMP_s = Corr.Fit_distr(-1.0*Dams*dmPi_s);


  distr_t der_DMK_l= DMK_l/Daml;
  distr_t der_DMK_crit= 2*DMK_crit/Damcr;
  distr_t der_DMK_s = DMK_s/Dams;

  distr_t der_DMDs_l= DMDs_l/Daml;
  distr_t der_DMDs_crit= DMDs_crit/Damcr;
  distr_t der_DMDs_s = DMDs_s/Dams;

  distr_t der_DMP_l = DMP_l/Daml;
  distr_t der_DMP_crit= 2.0*DMP_crit/Damcr;
  distr_t der_DMP_s = DMP_s/Dams;

  distr_t Mpi_TRUE= SCALE_INFO.Mpi[lens];
  distr_t fpi_TRUE= SCALE_INFO.fpi[lens];
  distr_t MK_TRUE= SCALE_INFO.MK1[sens];
  distr_t MDs_TRUE=SCALE_INFO.MDs1[cens];

  FPI_TRUE_list.distr_list.push_back( fpi_TRUE);

  distr_t_list DFPI_DMPCAC= fpi_TRUE*dfPi_crit/(1.1*dmPCAC);
  distr_t_list DMPI_DMPCAC=  dmPi_crit/(1.1*dmPCAC);

  Print_To_File({} , {DMPI_DMPCAC.ave(), DMPI_DMPCAC.err(), DFPI_DMPCAC.ave(), DFPI_DMPCAC.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_pcac", "", "");

  distr_t DFPI_l = fpi_TRUE*Daml*Corr.Fit_distr( dfPi_l);
  distr_t DFPI_s = fpi_TRUE*Dams*Corr.Fit_distr( dfPi_s);
  distr_t DFPI_crit = 0.5*fpi_TRUE*Damcr*Corr.Fit_distr( dfPi_crit);

  distr_t der_DFPI_l = DFPI_l/Daml;
  distr_t der_DFPI_s = DFPI_s/Dams;
  distr_t der_DFPI_crit = DFPI_crit/Damcr;


  FPI_der_l.distr_list.push_back( der_DFPI_l);
  FPI_der_s.distr_list.push_back( der_DFPI_s);
  FPI_der_crit.distr_list.push_back( der_DFPI_crit);

  FPI_l.distr_list.push_back( DFPI_l);
  FPI_s.distr_list.push_back( DFPI_s);
  FPI_crit.distr_list.push_back( DFPI_crit);


  //print results

  ofstream FL("../data/sea_quark_effects_NOVAL/table_Mpi.tex",ofstream::app);
  FL.precision(3);
  FL<<EN_RED<<" & "<<(pow(10,3)*DMP_l/a_distr).ave()<<"("<<(pow(10,3)*DMP_l/a_distr).err()<<") & "<<(DMP_l/Mpi_TRUE).ave()<<"("<<(DMP_l/Mpi_TRUE).err()<<") & "<<(DMP_l/Mpi_TRUE.err()).ave()<<"("<<(DMP_l/Mpi_TRUE.err()).err()<<") & ";
  FL<<(pow(10,3)*DMP_crit/a_distr).ave()<<"("<<(pow(10,3)*DMP_crit/a_distr).err()<<") & "<<(DMP_crit/Mpi_TRUE).ave()<<"("<<(DMP_crit/Mpi_TRUE).err()<<") & "<<(DMP_crit/Mpi_TRUE.err()).ave()<<"("<<(DMP_crit/Mpi_TRUE.err()).err()<<") & ";
  FL<<(pow(10,3)*DMP_s/a_distr).ave()<<"("<<(pow(10,3)*DMP_s/a_distr).err()<<") & "<<(DMP_s/Mpi_TRUE).ave()<<"("<<(DMP_s/Mpi_TRUE).err()<<") & "<<(DMP_s/Mpi_TRUE.err()).ave()<<"("<<(DMP_s/Mpi_TRUE.err()).err()<<") \\\\ \\hline"<<endl;
  FL.close();


  ofstream FK("../data/sea_quark_effects_NOVAL/table_MK.tex",ofstream::app);
  FK.precision(3);
  FK<<EN_RED<<" & "<<(pow(10,3)*DMK_l/a_distr).ave()<<"("<<(pow(10,3)*DMK_l/a_distr).err()<<") & "<<(DMK_l/MK_TRUE).ave()<<"("<<(DMK_l/MK_TRUE).err()<<") & "<<(DMK_l/MK_TRUE.err()).ave()<<"("<<(DMK_l/MK_TRUE.err()).err()<<") & ";
  FK<<(pow(10,3)*DMK_crit/a_distr).ave()<<"("<<(pow(10,3)*DMK_crit/a_distr).err()<<") & "<<(DMK_crit/MK_TRUE).ave()<<"("<<(DMK_crit/MK_TRUE).err()<<") & "<<(DMK_crit/MK_TRUE.err()).ave()<<"("<<(DMK_crit/MK_TRUE.err()).err()<<") & ";
  FK<<(pow(10,3)*DMK_s/a_distr).ave()<<"("<<(pow(10,3)*DMK_s/a_distr).err()<<") & "<<(DMK_s/MK_TRUE).ave()<<"("<<(DMK_s/MK_TRUE).err()<<") & "<<(DMK_s/MK_TRUE.err()).ave()<<"("<<(DMK_s/MK_TRUE.err()).err()<<") \\\\ \\hline"<<endl;
  FK.close();


  ofstream FDs("../data/sea_quark_effects_NOVAL/table_MDs.tex",ofstream::app);
  FDs.precision(3);
  FDs<<EN_RED<<" & "<<(pow(10,3)*DMDs_l/a_distr).ave()<<"("<<(pow(10,3)*DMDs_l/a_distr).err()<<") & "<<(DMDs_l/MDs_TRUE).ave()<<"("<<(DMDs_l/MDs_TRUE).err()<<") & "<<(DMDs_l/MDs_TRUE.err()).ave()<<"("<<(DMDs_l/MDs_TRUE.err()).err()<<") & ";
  FDs<<(pow(10,3)*DMDs_crit/a_distr).ave()<<"("<<(pow(10,3)*DMDs_crit/a_distr).err()<<") & "<<(DMDs_crit/MDs_TRUE).ave()<<"("<<(DMDs_crit/MDs_TRUE).err()<<") & "<<(DMDs_crit/MDs_TRUE.err()).ave()<<"("<<(DMDs_crit/MDs_TRUE.err()).err()<<") & ";
  FDs<<(pow(10,3)*DMDs_s/a_distr).ave()<<"("<<(pow(10,3)*DMDs_s/a_distr).err()<<") & "<<(DMDs_s/MDs_TRUE).ave()<<"("<<(DMDs_s/MDs_TRUE).err()<<") & "<<(DMDs_s/MDs_TRUE.err()).ave()<<"("<<(DMDs_s/MDs_TRUE.err()).err()<<") \\\\ \\hline"<<endl;
  FDs.close();


  ofstream FL_der("../data/sea_quark_effects_NOVAL/table_der_Mpi.tex",ofstream::app);
  FL_der.precision(3);
  FL_der<<EN_RED<<" & "<<der_DMP_l.ave()<<"("<<der_DMP_l.err()<<") & "<<der_DMP_crit.ave()<<"("<<der_DMP_crit.err()<<") & "<<der_DMP_s.ave()<<"("<<der_DMP_s.err()<<") \\\\ \\hline "<<endl;
  FL_der.close();

  ofstream FK_der("../data/sea_quark_effects_NOVAL/table_der_MK.tex",ofstream::app);
  FK_der.precision(3);
  FK_der<<EN_RED<<" & "<<der_DMK_l.ave()<<"("<<der_DMK_l.err()<<") & "<<der_DMK_crit.ave()<<"("<<der_DMK_crit.err()<<") & "<<der_DMK_s.ave()<<"("<<der_DMK_s.err()<<") \\\\ \\hline "<<endl;
  FK_der.close();

  ofstream FDs_der("../data/sea_quark_effects_NOVAL/table_der_MDs.tex",ofstream::app);
  FDs_der.precision(3);
  FDs_der<<EN_RED<<" & "<<der_DMDs_l.ave()<<"("<<der_DMDs_l.err()<<") & "<<der_DMDs_crit.ave()<<"("<<der_DMDs_crit.err()<<") & "<<der_DMDs_s.ave()<<"("<<der_DMDs_s.err()<<") \\\\ \\hline "<<endl;
  FDs_der.close();

  ofstream FPI_der("../data/sea_quark_effects_NOVAL/table_der_Fpi.tex",ofstream::app);
  FPI_der.precision(3);
  FPI_der<<EN_RED<<" & "<<der_DFPI_l.ave()<<"("<<der_DFPI_l.err()<<") & "<<der_DFPI_crit.ave()<<"("<<der_DFPI_crit.err()<<") & "<<der_DFPI_s.ave()<<"("<<der_DFPI_s.err()<<") \\\\ \\hline "<<endl;
  FPI_der.close();
	      


  
  

  //Damc=1.0;
  //Dams=1.0;
  //Daml=1.0;

 

  distr_t dfpi_mcr= Corr.Fit_distr(Damcr*dfPi_crit);

  cout<<"dfpi/fpi: "<<(0.5*dfpi_mcr).ave()<<" "<<(0.5*dfpi_mcr).err()<<" prediction from dfpi/fpi = sqrt(1 + (mpcac/mu)^2)-1: "<<(SQRT_D( 1 + ZA*ZA*mPCAC*mPCAC) -1).ave()<<" "<<(SQRT_D( 1 + ZA*ZA*mPCAC*mPCAC) -1).err()<<endl;

  Print_To_File( {}, {(-Daml*dmPi_l/(a_distr)).ave(), (-Daml*dmPi_l/a_distr).err(), (Daml*dfPi_l).ave(), (Daml*dfPi_l).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_l", "", "");
  Print_To_File( {}, {(-Dams*dmPi_s/(a_distr)).ave(), (-Dams*dmPi_s/a_distr).err(), (Dams*dfPi_s).ave(), (Dams*dfPi_s).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_s", "", "");
  Print_To_File( {}, {(-Damc*dmPi_c/(a_distr)).ave(), (-Damc*dmPi_c/a_distr).err(), (Damc*dfPi_c).ave(), (Damc*dfPi_c).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_c", "", "");
  Print_To_File( {}, {(-1.0*Damcr*dmPi_crit/(a_distr)).ave(), (-1.0*Damcr*dmPi_crit/a_distr).err(), (1.0*Damcr*dfPi_crit).ave(), (1.0*Damcr*dfPi_crit).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmfPi_crit", "", "");

  Print_To_File( {}, {(-Daml*dmK_l/(a_distr)).ave(), (-Daml*dmK_l/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmK_l", "", "");
  Print_To_File( {}, {(-Dams*dmK_s/(a_distr)).ave(), (-Dams*dmK_s/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmK_s", "", "");
  Print_To_File( {}, {(-1.0*Damcr*dmK_crit/(a_distr)).ave(), (-1.0*Damcr*dmK_crit/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmK_crit", "", "");

  Print_To_File( {}, {(-Daml*dmDs_l/(a_distr)).ave(), (-Daml*dmDs_l/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmDs_l", "", "");
  Print_To_File( {}, {(-Dams*dmDs_s/(a_distr)).ave(), (-Dams*dmDs_s/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmDs_s", "", "");
  Print_To_File( {}, {(-1.0*Damcr*dmDs_crit/(a_distr)).ave(), (-1.0*Damcr*dmDs_crit/a_distr).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/dmDs_crit", "", "");

  Print_To_File( {}, {dmPi_l.ave(), dmPi_l.err(), (a_distr*dfPi_l).ave(), (a_distr*dfPi_l).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmfPi_l", "", "");
  Print_To_File( {}, {dmPi_s.ave(), dmPi_s.err(), (a_distr*dfPi_s).ave(), (a_distr*dfPi_s).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmfPi_s", "", "");
  Print_To_File( {}, {dmPi_c.ave(), dmPi_c.err(), (a_distr*dfPi_c).ave(), (a_distr*dfPi_c).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmfPi_c", "", "");
  Print_To_File( {}, {dmPi_crit.ave(), dmPi_crit.err(), (a_distr*dfPi_crit).ave(), (a_distr*dfPi_crit).err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmfPi_crit", "", "");


  Print_To_File( {}, {dmK_l.ave(), dmK_l.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmK_l", "", "");
  Print_To_File( {}, {dmK_s.ave(), dmK_s.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmK_s", "", "");
  Print_To_File( {}, {dmK_crit.ave(), dmK_crit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmK_crit", "", "");

  Print_To_File( {}, {dmDs_l.ave(), dmDs_l.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmDs_l", "", "");
  Print_To_File( {}, {dmDs_s.ave(), dmDs_s.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmDs_s", "", "");
  Print_To_File( {}, {dmDs_crit.ave(), dmDs_crit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_dmDs_crit", "", "");



  

  //define amu^HVP kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , T/2);
  double t0=  0.4*fm_to_inv_Gev;
  double t1 = 1.0*fm_to_inv_Gev;
  double Delta= 0.15*fm_to_inv_Gev;
  auto th0 = [&t0, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
  auto th1 = [&t1, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};


  //ss Win and SD
  distr_t_list damu_W_TM_ss_dml(UseJack);
  distr_t_list damu_W_TM_ss_dms(UseJack);
  distr_t_list damu_W_TM_ss_dmcr(UseJack);
  distr_t_list damu_W_TM_ss_dmc(UseJack);
  distr_t_list damu_SD_TM_ss_dml(UseJack);
  distr_t_list damu_SD_TM_ss_dms(UseJack);
  distr_t_list damu_SD_TM_ss_dmcr(UseJack);
  distr_t_list damu_SD_TM_ss_dmc(UseJack);

  distr_t_list damu_W_TM_cc_dml(UseJack);
  distr_t_list damu_W_TM_cc_dms(UseJack);
  distr_t_list damu_W_TM_cc_dmcr(UseJack);
  distr_t_list damu_W_TM_cc_dmc(UseJack);

  distr_t_list damu_SD_TM_cc_dml(UseJack);
  distr_t_list damu_SD_TM_cc_dms(UseJack);
  distr_t_list damu_SD_TM_cc_dmcr(UseJack);
  distr_t_list damu_SD_TM_cc_dmc(UseJack);
  

  distr_t_list amu_W_TM_ss(UseJack);


  distr_t_list Integrand_damu_TM_ss_dml(UseJack);
  distr_t_list Integrand_damu_TM_ss_dms(UseJack);
  distr_t_list Integrand_damu_TM_ss_dmc(UseJack);
  distr_t_list Integrand_damu_TM_ss_dmcr(UseJack);

  
  for(int t=0;t<T/2;t++) {

    if(correlate_to_VKVK_strange) {
      //strange
      amu_TM_ss1.distr_list.push_back( (t==0)?(0.0*id_distr):( amu_TM_ss1[(t==0)?0:t-1] +  4.0*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_distr[t]*Ker.distr_list[t] ));
      
      damu_TM_ss1_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_ss1_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dms_distr[t]*Ker.distr_list[t] ));
      damu_OS_ss1_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_ss1_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_ss1_dms_distr[t]*Ker.distr_list[t] ));
      damu_TM_ss2_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_ss2_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss2_dms_distr[t]*Ker.distr_list[t] ));
      damu_OS_ss2_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_ss2_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_ss2_dms_distr[t]*Ker.distr_list[t] ));
      damu_TM_ss1_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_ss1_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmc_distr[t]*Ker.distr_list[t] ));
      damu_TM_ss1_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_ss1_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dml_distr[t]*Ker.distr_list[t] ));
      damu_OS_ss1_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_ss1_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_ss1_dmc_distr[t]*Ker.distr_list[t] ));
      damu_OS_ss1_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_ss1_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_ss1_dml_distr[t]*Ker.distr_list[t] ));
      

      Integrand_damu_TM_ss_dml.distr_list.push_back(  4.0*Daml*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dml_distr[t]*Ker.distr_list[t] );
      Integrand_damu_TM_ss_dms.distr_list.push_back(  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dms_distr[t]*Ker.distr_list[t] );
      Integrand_damu_TM_ss_dmc.distr_list.push_back(  4.0*Damc*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmc_distr[t]*Ker.distr_list[t] );
      Integrand_damu_TM_ss_dmcr.distr_list.push_back(  4.0*Damcr*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmcrit_distr[t]*Ker.distr_list[t] );
      
      //intermediate and short-distance window
      //W
      amu_W_TM_ss.distr_list.push_back(  (t==0)?(0.0*id_distr):( amu_W_TM_ss[(t==0)?0:t-1] +  4.0*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    ) ));
      damu_W_TM_ss_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_ss_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dml_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    )));
      damu_W_TM_ss_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_ss_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dms_distr[t]*Ker.distr_list[t]*(distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) )   ));
      damu_W_TM_ss_dmcr.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_ss_dmcr[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmcrit_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    )));
      damu_W_TM_ss_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_ss_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmc_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)  )  ));
      //SD
      damu_SD_TM_ss_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_ss_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dml_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) ); 
      damu_SD_TM_ss_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_ss_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dms_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) );  
      damu_SD_TM_ss_dmcr.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_ss_dmcr[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmcrit_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) ); 
      damu_SD_TM_ss_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_ss_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmc_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) );     
      
      
      //critical mass corrections
      damu_TM_ss1_dmcrit.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_ss1_dmcrit[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_ss1_dmcrit_distr[t]*Ker.distr_list[t] ));
      damu_OS_ss1_dmcrit.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_ss1_dmcrit[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qd*qd)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_ss1_dmcrit_distr[t]*Ker.distr_list[t] ));
      
      
           
    }
    
   
    //charm
    if(correlate_to_VKVK_charm) {
      damu_TM_cc1_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_cc1_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmc_distr[t]*Ker.distr_list[t] ));
      damu_OS_cc1_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_cc1_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_cc1_dmc_distr[t]*Ker.distr_list[t] ));
      damu_TM_cc2_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_cc2_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc2_dmc_distr[t]*Ker.distr_list[t] ));
      damu_OS_cc2_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_cc2_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_cc2_dmc_distr[t]*Ker.distr_list[t] ));
      damu_TM_cc1_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_cc1_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dms_distr[t]*Ker.distr_list[t] ));
      damu_TM_cc1_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_cc1_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dml_distr[t]*Ker.distr_list[t] ));
      damu_TM_cc1_dmcrit.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_TM_cc1_dmcrit[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmcrit_distr[t]*Ker.distr_list[t] ));
      damu_OS_cc1_dmcrit.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_OS_cc1_dmcrit[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZV*ZV*VKVK_OS_cc1_dmcrit_distr[t]*Ker.distr_list[t] ));
   
   
      damu_W_TM_cc_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_cc_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dml_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    )));
      damu_W_TM_cc_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_cc_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dms_distr[t]*Ker.distr_list[t]*(distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) )   ));
      damu_W_TM_cc_dmcr.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_cc_dmcr[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmcrit_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    )));
      damu_W_TM_cc_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_W_TM_cc_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmc_distr[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)  )  ));
      //SD
      damu_SD_TM_cc_dml.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_cc_dml[(t==0)?0:t-1] +  4.0*Daml*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dml_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) ); 
      damu_SD_TM_cc_dms.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_cc_dms[(t==0)?0:t-1] +  4.0*Dams*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dms_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) );  
      damu_SD_TM_cc_dmcr.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_cc_dmcr[(t==0)?0:t-1] +  4.0*Damcr*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmcrit_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) ); 
      damu_SD_TM_cc_dmc.distr_list.push_back( (t==0)?(0.0*id_distr):( damu_SD_TM_cc_dmc[(t==0)?0:t-1] +  4.0*Damc*w(t,1)*(qu*qu)*pow(alpha,2)*1e10*ZA*ZA*VKVK_TM_cc1_dmc_distr[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr))) );   




      
    }
  
  }




  distr_t DDmc= a_distr/Damc;
  distr_t DDms= a_distr/Dams;
  distr_t DDml= a_distr/Daml;
  distr_t DDmcrit = a_distr/(Damcr*dmPCAC*ZA);

  Deltas.distr_list.push_back( Damc/a_distr);

  if(correlate_to_VKVK_strange) {
    Print_To_File({}, {damu_TM_ss1_dms.ave(), damu_TM_ss1_dms.err(), damu_TM_ss2_dms.ave(), damu_TM_ss2_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_TM_dms", "", "");
    Print_To_File({}, {damu_OS_ss1_dms.ave(), damu_OS_ss1_dms.err(), damu_OS_ss2_dms.ave(), damu_OS_ss2_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_OS_dms", "", "");
    Print_To_File({}, {amu_TM_ss1.ave(), amu_TM_ss1.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_TM_dmc", "", "");
    Print_To_File({}, {amu_TM_ss1.ave(), amu_TM_ss1.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/amu_s_TM_dmc", "", "");
    Print_To_File({}, {damu_TM_ss1_dmc.ave(), damu_TM_ss1_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_TM_dmc", "", "");
    Print_To_File({}, {damu_TM_ss1_dml.ave(), damu_TM_ss1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_TM_dml", "", "");
    Print_To_File({}, {damu_OS_ss1_dmc.ave(), damu_OS_ss1_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_OS_dmc", "", "");
    Print_To_File({}, {damu_OS_ss1_dml.ave(), damu_OS_ss1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_OS_dml", "", "");
    //critical mass corrections
    //print correlator
    Print_To_File({}, {(VKVK_TM_ss1_distr+Daml*VKVK_TM_ss1_dml_distr).ave(), (VKVK_TM_ss1_distr+Daml*VKVK_TM_ss1_dml_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/C_s_TM_dml", "", "");
    Print_To_File({}, {(VKVK_TM_ss1_distr+Dams*VKVK_TM_ss1_dms_distr).ave(), (VKVK_TM_ss1_distr+Dams*VKVK_TM_ss1_dms_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/C_s_TM_dms", "", "");
    Print_To_File({}, {(VKVK_TM_ss1_distr+Damcr*VKVK_TM_ss1_dmcrit_distr).ave(), (VKVK_TM_ss1_distr+Damcr*VKVK_TM_ss1_dmcrit_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/C_s_TM_dmcr", "", "");
    Print_To_File({}, {(VKVK_TM_ss1_distr+Damc*VKVK_TM_ss1_dmc_distr).ave(), (VKVK_TM_ss1_distr+Damc*VKVK_TM_ss1_dmc_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/C_s_TM_dmc", "", "");

    Print_To_File({}, {(Daml*VKVK_TM_ss1_dml_distr).ave(), (Daml*VKVK_TM_ss1_dml_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dC_s_TM_dml", "", "");
    Print_To_File({}, {(Dams*VKVK_TM_ss1_dms_distr).ave(), (Dams*VKVK_TM_ss1_dms_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dC_s_TM_dms", "", "");
    Print_To_File({}, {(Damcr*VKVK_TM_ss1_dmcrit_distr).ave(), (Damcr*VKVK_TM_ss1_dmcrit_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dC_s_TM_dmcr", "", "");
    Print_To_File({}, {(Damc*VKVK_TM_ss1_dmc_distr).ave(), (Damc*VKVK_TM_ss1_dmc_distr).err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dC_s_TM_dmc", "", "");

    Print_To_File({}, {Integrand_damu_TM_ss_dml.ave(), Integrand_damu_TM_ss_dml.err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dInt_s_TM_dml", "", "");
    Print_To_File({}, {Integrand_damu_TM_ss_dms.ave(), Integrand_damu_TM_ss_dms.err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dInt_s_TM_dms", "", "");
    Print_To_File({}, {Integrand_damu_TM_ss_dmc.ave(), Integrand_damu_TM_ss_dmc.err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dInt_s_TM_dmc", "", "");
    Print_To_File({}, {Integrand_damu_TM_ss_dmcr.ave(), Integrand_damu_TM_ss_dmcr.err()},  "../data/sea_quark_effects_NOVAL/"+Ens+"/dInt_s_TM_dmcr", "", "");

    
    Print_To_File({}, {damu_TM_ss1_dmcrit.ave(), damu_TM_ss1_dmcrit.err(), damu_OS_ss1_dmcrit.ave(), damu_OS_ss1_dmcrit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_s_dmcrit", "", "");
  
  

    int tsep=(int)(2.5*fm_to_inv_Gev/a_distr.ave());
    int tsep_charm= (int)(1.5*fm_to_inv_Gev/a_distr.ave());

    if(Ens.substr(1,1)=="C") {
      TOT_dmlc_for_E112= 1.5*damu_TM_ss1_dmcrit[tsep]+damu_TM_ss1_dml[tsep];
      W_dmlc_for_E112= 1.5*damu_W_TM_ss_dmcr[tsep]+damu_W_TM_ss_dml[tsep];
      SD_dmlc_for_E112= 1.5*damu_SD_TM_ss_dmcr[tsep]+damu_SD_TM_ss_dml[tsep];
    }

    cout<<"intermediate-window: "<<amu_W_TM_ss.ave(Corr.Nt/2-10)<<" +- "<<amu_W_TM_ss.err(Corr.Nt/2-10)<<endl;

    cout<<"Damu(s)-ml+ms+mc"<<endl;

    double k_LD= -60*Damc/(12*a_distr.ave());
    double k_W = -14*Damc/(12*a_distr.ave());
    double k_SD = damu_SD_TM_ss_dms.ave(tsep)*(Damc/(12*Dams));

    cout<<"derivative amu^W :"<<(damu_W_TM_ss_dms[tsep]*DDms).ave()<<" "<<(damu_W_TM_ss_dms[tsep]*DDms).err()<<endl;
    cout<<"derivative amu^SD :"<<(damu_SD_TM_ss_dms[tsep]*DDms).ave()<<" "<<(damu_SD_TM_ss_dms[tsep]*DDms).err()<<endl;
    
    
    if(Ens.substr(1,1) != "E") {

      distr_t damu_HVP = damu_TM_ss1_dml[tsep]+damu_TM_ss1_dms[tsep]+damu_TM_ss1_dmcrit[tsep];
      distr_t damu_W= damu_W_TM_ss_dml[tsep]+damu_W_TM_ss_dms[tsep]+damu_W_TM_ss_dmcr[tsep];
      distr_t damu_SD= damu_SD_TM_ss_dml[tsep]+damu_SD_TM_ss_dms[tsep]+damu_SD_TM_ss_dmcr[tsep];
      distr_t damu_LD = damu_HVP-damu_W-damu_SD;

      cout<<"Ens: "<<Ens<<" Damu(s,HVP): "<<damu_HVP.ave()<<" +- "<<sqrt( pow(damu_HVP.err(),2) + pow(k_LD,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,W): "<<damu_W.ave()<<" +- "<<sqrt( pow(damu_W.err(),2) + pow(k_W,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,SD): "<<damu_SD.ave()<<" +- "<<sqrt( pow(damu_SD.err(),2) + pow(k_SD,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,LD): "<<damu_LD.ave()<<" +- "<<sqrt( pow(damu_LD.err(),2) + pow(k_LD-k_W-k_SD,2))<<endl;
      

     
    }
    else {

      distr_t damu_HVP = damu_TM_ss1_dms[tsep]+ TOT_dmlc_for_E112;
      distr_t damu_W = damu_W_TM_ss_dms[tsep]+ W_dmlc_for_E112;
      distr_t damu_SD = damu_SD_TM_ss_dms[tsep]+ SD_dmlc_for_E112;
      distr_t damu_LD = damu_HVP-damu_W-damu_SD;

      cout<<"Ens: "<<Ens<<" Damu(s,HVP): "<<damu_HVP.ave()<<" +- "<<sqrt( pow(damu_HVP.err(),2) + pow(k_LD,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,W): "<<damu_W.ave()<<" +- "<<sqrt( pow(damu_W.err(),2) + pow(k_W,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,SD): "<<damu_SD.ave()<<" +- "<<sqrt( pow(damu_SD.err(),2) + pow(k_SD,2))<<endl;
      cout<<"Ens: "<<Ens<<" Damu(s,LD): "<<damu_LD.ave()<<" +- "<<sqrt( pow(damu_LD.err(),2) + pow(k_LD-k_W-k_SD,2))<<endl;
        
    }

  
    

    damu_TM_ss1_dms = damu_TM_ss1_dms*DDms;
    damu_TM_ss2_dms = damu_TM_ss2_dms*DDms;
    damu_OS_ss1_dms = damu_OS_ss1_dms*DDms;
    damu_OS_ss2_dms = damu_OS_ss2_dms*DDms;
    damu_TM_ss1_dmc = damu_TM_ss1_dmc*DDmc;
    damu_OS_ss1_dmc = damu_OS_ss1_dmc*DDmc;
    damu_TM_ss1_dml = damu_TM_ss1_dml*DDml;
    damu_OS_ss1_dml = damu_OS_ss1_dml*DDml;
    damu_TM_ss1_dmcrit = damu_TM_ss1_dmcrit*DDmcrit;
    damu_OS_ss1_dmcrit = damu_OS_ss1_dmcrit*DDmcrit;
  
    damu_s_charm.distr_list.push_back( damu_TM_ss1_dmc[tsep_charm]);
    damu_W_s_charm.distr_list.push_back( damu_W_TM_ss_dmc[tsep_charm]);
    damu_SD_s_charm.distr_list.push_back( damu_SD_TM_ss_dmc[tsep_charm]);
  

    
    // derivatives
    Print_To_File({}, {damu_TM_ss1_dms.ave(), damu_TM_ss1_dms.err(), damu_TM_ss2_dms.ave(), damu_TM_ss2_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_TM_dms", "", "");
    Print_To_File({}, {damu_OS_ss1_dms.ave(), damu_OS_ss1_dms.err(), damu_OS_ss2_dms.ave(), damu_OS_ss2_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_OS_dms", "", "");
    Print_To_File({}, {damu_TM_ss1_dmc.ave(), damu_TM_ss1_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_TM_dmc", "", "");
    Print_To_File({}, {damu_TM_ss1_dml.ave(), damu_TM_ss1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_TM_dml", "", "");
    Print_To_File({}, {damu_OS_ss1_dmc.ave(), damu_OS_ss1_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_OS_dmc", "", "");
    Print_To_File({}, {damu_OS_ss1_dml.ave(), damu_OS_ss1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_OS_dml", "", "");
    //critical mass corrections
    Print_To_File({}, {damu_TM_ss1_dmcrit.ave(), damu_TM_ss1_dmcrit.err(), damu_OS_ss1_dmcrit.ave(), damu_OS_ss1_dmcrit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_s_dmcrit", "", "");
  


   

    

    
   
  }

  if(correlate_to_VKVK_charm) {
    Print_To_File({}, {damu_TM_cc1_dmc.ave(), damu_TM_cc1_dmc.err(), damu_TM_cc2_dmc.ave(), damu_TM_cc2_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_c_TM_dmc", "", "");
    Print_To_File({}, {damu_OS_cc1_dmc.ave(), damu_OS_cc1_dmc.err(), damu_OS_cc2_dmc.ave(), damu_OS_cc2_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_c_OS_dmc", "", "");
    Print_To_File({}, {damu_TM_cc1_dms.ave(), damu_TM_cc1_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_c_TM_dms", "", "");
    Print_To_File({}, {damu_TM_cc1_dml.ave(), damu_TM_cc1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_c_TM_dml", "", "");
    Print_To_File({}, {damu_TM_cc1_dmcrit.ave(), damu_TM_cc1_dmcrit.err(), damu_OS_cc1_dmcrit.ave(), damu_OS_cc1_dmcrit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/damu_c_dmcrit", "", "");




    int tsep=(int)(2.0*fm_to_inv_Gev/a_distr.ave());
    
    if(Ens.substr(1,1)=="C") {
      TOT_C_dmlc_for_E112= 1.5*damu_TM_cc1_dmcrit[tsep]+damu_TM_cc1_dml[tsep];
      W_C_dmlc_for_E112= 1.5*damu_W_TM_cc_dmcr[tsep]+damu_W_TM_cc_dml[tsep];
      SD_C_dmlc_for_E112= 1.5*damu_SD_TM_cc_dmcr[tsep]+damu_SD_TM_cc_dml[tsep];
    }

    cout<<"intermediate-window: "<<amu_W_TM_ss.ave(Corr.Nt/2-10)<<" +- "<<amu_W_TM_ss.err(Corr.Nt/2-10)<<endl;

    cout<<"Damu(s)-ml+ms+mc"<<endl;

    double k_LD= -0.0*60*Damc/(12*a_distr.ave());
    double k_W = -14*0.0*Damc/(12*a_distr.ave());
    double k_SD = 0.0*damu_SD_TM_cc_dms.ave(tsep)*(Damc/(12*Dams));

    cout<<"derivative amu^W(c) :"<<(damu_W_TM_cc_dms[tsep]*DDms).ave()<<" "<<(damu_W_TM_cc_dms[tsep]*DDms).err()<<endl;
    cout<<"derivative amu^SD(c) :"<<(damu_SD_TM_cc_dms[tsep]*DDms).ave()<<" "<<(damu_SD_TM_cc_dms[tsep]*DDms).err()<<endl;
    
    
    if(Ens.substr(1,1) != "E") {

      distr_t damu_HVP= damu_TM_cc1_dml[tsep]+damu_TM_cc1_dms[tsep]+damu_TM_cc1_dmcrit[tsep]+ ( damu_TM_cc1_dms[tsep] - damu_TM_cc1_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_W= damu_W_TM_cc_dml[tsep]+damu_W_TM_cc_dms[tsep]+damu_W_TM_cc_dmcr[tsep] + (damu_W_TM_cc_dms[tsep]- damu_W_TM_cc_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_SD = damu_SD_TM_cc_dml[tsep]+damu_SD_TM_cc_dms[tsep]+damu_SD_TM_cc_dmcr[tsep] + (damu_SD_TM_cc_dms[tsep] -damu_SD_TM_cc_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_LD= damu_HVP-damu_W-damu_SD;

      cout<<"Ens: "<<Ens<<" Damu(c,HVP): "<<damu_HVP.ave()<<" +- "<<damu_HVP.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,W): "<<damu_W.ave()<<" +- "<<damu_W.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,SD): "<<damu_SD.ave()<<" +- "<<damu_SD.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,LD): "<<damu_LD.ave()<<" +- "<<damu_LD.err()<<endl;
      
    }
    else {

      distr_t damu_HVP = damu_TM_cc1_dms[tsep]+ TOT_C_dmlc_for_E112 + (damu_TM_cc1_dms[tsep]-damu_TM_cc1_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_W = damu_W_TM_cc_dms[tsep]+ W_C_dmlc_for_E112 +  (damu_W_TM_cc_dms[tsep]-damu_W_TM_cc_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_SD = damu_SD_TM_cc_dms[tsep]+ SD_C_dmlc_for_E112 + (damu_SD_TM_cc_dms[tsep]-damu_SD_TM_cc_dms.ave(tsep))*(Damc/Dams)/12;
      distr_t damu_LD= damu_HVP-damu_W-damu_SD;

      cout<<"Ens: "<<Ens<<" Damu(c,HVP): "<<damu_HVP.ave()<<" +- "<<damu_HVP.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,W): "<<damu_W.ave()<<" +- "<<damu_W.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,SD): "<<damu_SD.ave()<<" +- "<<damu_SD.err()<<endl;
      cout<<"Ens: "<<Ens<<" Damu(c,LD): "<<damu_LD.ave()<<" +- "<<damu_LD.err()<<endl;
      
  
    }

  
    

    damu_TM_cc1_dmc = damu_TM_cc1_dmc*DDmc;
    damu_TM_cc2_dmc = damu_TM_cc2_dmc*DDmc;
    damu_OS_cc1_dmc = damu_OS_cc1_dmc*DDmc;
    damu_OS_cc2_dmc = damu_OS_cc2_dmc*DDmc;
    
    damu_TM_cc1_dms = damu_TM_cc1_dms*DDms;
    damu_TM_cc1_dml = damu_TM_cc1_dml*DDml;
    damu_TM_cc1_dmcrit = damu_TM_cc1_dmcrit*DDmcrit;
    damu_OS_cc1_dmcrit = damu_OS_cc1_dmcrit*DDmcrit;
  
      
    // derivatives
    Print_To_File({}, {damu_TM_cc1_dmc.ave(), damu_TM_cc1_dmc.err(), damu_TM_cc2_dmc.ave(), damu_TM_cc2_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_c_TM_dmc", "", "");
    Print_To_File({}, {damu_OS_cc1_dmc.ave(), damu_OS_cc1_dmc.err(), damu_OS_cc2_dmc.ave(), damu_OS_cc2_dmc.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_c_OS_dmc", "", "");
    Print_To_File({}, {damu_TM_cc1_dml.ave(), damu_TM_cc1_dml.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_c_TM_dml", "", "");
    Print_To_File({}, {damu_TM_cc1_dms.ave(), damu_TM_cc1_dms.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_c_TM_dms", "", "");
    //critical macc corrections
    Print_To_File({}, {damu_TM_cc1_dmcrit.ave(), damu_TM_cc1_dmcrit.err(), damu_OS_cc1_dmcrit.ave(), damu_OS_cc1_dmcrit.err()}, "../data/sea_quark_effects_NOVAL/"+Ens+"/derivative_damu_c_dmcrit", "", "");









    
 
    cout<<"HVP corrections computed!"<<endl;
  }



  

  
  
     
    }


    Vfloat wts;
    Vfloat wts_3;
    double sum_w3=0.0;
    double sum_w=0.0;
    for(int i=0;i<(signed)Ens_list.size();i++) {
      double w= 1.0/pow(FPI_der_s.err(i),2);
      sum_w +=w;
      wts.push_back(w);
      if(Ens_list[i] != "cC211a.06.80") {
	wts_3.push_back(w);
	sum_w3 +=w;
      }
      else { wts_3.push_back(0.0);}
      
    }
    //normalize
    for(int i=0;i<(signed)Ens_list.size();i++)  wts[i] /= sum_w;
    for(auto &ww: wts_3) ww /= sum_w3;
    
    distr_t FPI_s_der_FIT= 0.0*id_distr;
    distr_t FPI_s_der_FIT3=0.0*id_distr;
    //fit
    for(int i=0;i<(signed)Ens_list.size();i++) { FPI_s_der_FIT = FPI_s_der_FIT + wts[i]*FPI_der_s[i]; FPI_s_der_FIT3 = FPI_s_der_FIT3 + wts_3[i]*FPI_der_s[i];}
    //assign values of FPI_s
    distr_t_list FPI_s_FIT = FPI_s_der_FIT*Dams_list;
    FPI_s_FIT = FPI_s_der_FIT3*Dams_list;
    distr_t_list FPI_s_FIT_LAT = pow(10,-3)*FPI_s_FIT*a_distr_list;
   

    cout<<"dfpi/dms 3 beta: "<<FPI_s_der_FIT3.ave()<<" "<<FPI_s_der_FIT3.err()<<endl;
    cout<<"dfpi/dms 4 beta: "<<FPI_s_der_FIT.ave()<<" "<<FPI_s_der_FIT.err()<<endl;
    


    //constant fit to dfpi/ds and print
    for(int iens=0;iens<(signed)Ens_list.size();iens++) {

      ofstream FL("../data/sea_quark_effects_NOVAL/table_Fpi.tex",ofstream::app);
      FL.precision(3);
      FL<<EN_RED_list[iens]<<" & "<<(FPI_der_l*Daml_list).ave(iens)<<"("<<(FPI_der_l*Daml_list).err(iens)<<") & "<<(FPI_l/FPI_TRUE_list).ave(iens)<<"("<<(FPI_l/FPI_TRUE_list).err(iens)<<") & "<<(FPI_l/FPI_TRUE_list.err(iens)).ave(iens)<<"("<<(FPI_l/FPI_TRUE_list.err(iens)).err(iens)<<") & ";
      FL<<(FPI_der_crit*Damcr_list).ave(iens)<<"("<<(FPI_der_crit*Damcr_list).err(iens)<<") & "<<(FPI_crit/FPI_TRUE_list).ave(iens)<<"("<<(FPI_crit/FPI_TRUE_list).err(iens)<<") & "<<(FPI_crit/FPI_TRUE_list.err(iens)).ave(iens)<<"("<<(FPI_crit/FPI_TRUE_list.err(iens)).err(iens)<<") & ";
      FL<<(FPI_s_FIT).ave(iens)<<"("<<(FPI_s_FIT).err(iens)<<") & "<<(FPI_s_FIT_LAT/FPI_TRUE_list).ave(iens)<<"("<<(FPI_s_FIT_LAT/FPI_TRUE_list).err(iens)<<") & "<<(FPI_s_FIT_LAT/FPI_TRUE_list.err(iens)).ave(iens)<<"("<<(FPI_s_FIT_LAT/FPI_TRUE_list.err(iens)).err(iens)<<") \\\\ \\hline"<<endl;
      FL.close();
    }

   
    
    
  return;
}



