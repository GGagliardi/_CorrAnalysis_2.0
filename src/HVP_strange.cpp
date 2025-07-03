#include "../include/HVP.h"
#include "RC_WI_analysis.h"
#include "scale_setting_main.h"


const double DTT = 0.5;
const double alpha = 1.0/137.035999;
const double MK_FLAG = 0.494600000;

using namespace std;








void HVP_strange() {





 
  

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
  ZA_A_ave = a_info.Za_WI_FLAG;
  ZA_A_err = a_info.Za_WI_FLAG_err;
  ZV_A_ave = a_info.Zv_WI_FLAG;
  ZV_A_err = a_info.Zv_WI_FLAG_err;
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
      

  //############################################################################################

















  bool Get_ASCII_strange= true;

    if(Get_ASCII_strange) {
    //read binary files
    boost::filesystem::create_directory("../HVP_strange");
    

    vector<string> Ens_T1({"B.72.96", "B.72.64", "C.06.80", "C.06.112", "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cB211b.072.96", "cB211b.072.64", "cC211a.06.80", "cC211a.06.112", "cD211a.054.96", "cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_s1_s1",  "mix_s2_s2"});

      cout<<"Getting data for amu(s1-s2) on Ensemble: "<<Ens_T1[it]<<endl;

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../HVP_strange/"+channel);
	boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "OS_VKVK"});

          
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
	  boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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

    exit(-1);

    bool Get_ASCII_3rd_strange= true;

    if(Get_ASCII_3rd_strange) {
    //read binary files
    boost::filesystem::create_directory("../HVP_strange");
    

    vector<string> Ens_T1({"B.72.64", "C.06.80", "C.06.112", "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cB211b.072.64", "cC211a.06.80", "cC211a.06.112", "cD211a.054.96", "cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      cout<<"Getting data for amu(s3) on Ensemble: "<<Ens_T1[it]<<endl;

      vector<string> channels({"mix_s_s"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../HVP_strange/"+channel);
	boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "OS_VKVK"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../tau_decay_strange_bin_mu_corr/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
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
	  boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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
   



    boost::filesystem::create_directory("../data/HVP_strange");

    data_t Vk_data_tm, Vk_data_OS;
    data_t Vk_data_s2_tm, Vk_data_s2_OS;
    data_t Vk_data_s3_tm, Vk_data_s3_OS;
    
    Vk_data_tm.Read("../HVP_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
    Vk_data_OS.Read("../HVP_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);

    Vk_data_s2_tm.Read("../HVP_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKVK", "VKVK", Sort_easy);
    Vk_data_s2_OS.Read("../HVP_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKVK", "VKVK", Sort_easy);

    Vk_data_s3_tm.Read("../HVP_strange/mix_s_s", "mes_contr_mix_s_s_TM_VKVK", "VKVK", Sort_easy);
    Vk_data_s3_OS.Read("../HVP_strange/mix_s_s", "mes_contr_mix_s_s_OS_VKVK", "VKVK", Sort_easy);
    
    int Nens= Vk_data_tm.size;

    distr_t_list amu_TM_s(UseJack), amu_OS_s(UseJack), W_TM_s(UseJack), W_OS_s(UseJack), amu_TM_s_uncorr(UseJack), W_TM_s_uncorr(UseJack);
    distr_t_list a_distr_list(UseJack);


    scale_setting_info SCALE_INFO=Get_scale_setting_info();
    RCs_info RCs_INFO= Get_RCs("ls");
    RCs_info RCs_INFO_ss= Get_RCs("ss");

    GaussianMersenne GG(4387);

    cout<<"Nens: "<<Nens<<endl;


    for(int iens=0;iens<Nens;iens++) {
      
      cout<<"Analyzing ensemble: "<<Vk_data_tm.Tag[iens]<<endl;

      int m_ens=-1;
      for(int k=0;k<(signed)SCALE_INFO.Ens.size();k++) {
	if(SCALE_INFO.Ens[k] == Vk_data_tm.Tag[iens]) m_ens=k; 
      }
      if(m_ens==-1) crash("Cannot find info on ms and ms_OS in SCALE_INFO for Ens: "+Vk_data_tm.Tag[iens]);

      int RC_ens=-1;
      for(int k=0;k<(signed)RCs_INFO.Ens.size();k++) {
	if(RCs_INFO.Ens[k] == Vk_data_tm.Tag[iens]) RC_ens=k; 
      }
      if(RC_ens==-1) crash("Cannot find info on ms and ms_OS in SCALE_INFO for Ens: "+Vk_data_tm.Tag[iens]);

      distr_t ms_phys= SCALE_INFO.ms[m_ens];
      distr_t ms_OS_phys = SCALE_INFO.ms_OS[m_ens];
      distr_t ms_resampled(UseJack);
      if(Vk_data_tm.Tag[iens] != "cD211a.054.96") {
	for(int i=0;i<Njacks;i++) ms_resampled.distr.push_back( ms_phys.ave() + GG()*ms_phys.err()/sqrt(Njacks-1.0));
      }
      else {
	for(int i=0;i<Njacks;i++) ms_resampled.distr.push_back( 0.013542 + GG()*0.000028/sqrt(Njacks-1.0));
      }
	
      CorrAnalysis Corr(UseJack, Njacks,800);
      Corr.Nt = Vk_data_tm.nrows[iens];
     
      distr_t a_distr(UseJack);
      distr_t Zv(UseJack), Za(UseJack);

      distr_t MK1= SCALE_INFO.MK1[m_ens];
      distr_t MK2= SCALE_INFO.MK2[m_ens];

      double ams3;

      if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=SCALE_INFO.a_B; Zv = ZV_B; Za = ZA_B; ams3= 0.0182;}
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=SCALE_INFO.a_C; Zv = ZV_C; Za = ZA_C; ams3= 0.016067;}
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=SCALE_INFO.a_D; Zv = ZV_D; Za = ZA_D; ams3= 1.3557e-02;}
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="E") {a_distr=SCALE_INFO.a_E; Zv = ZV_E; Za = ZA_E; ams3= 1.1759e-02;}
      else crash("Ensemble not found");

      

      distr_t Zv_fl(UseJack), Za_fl(UseJack);

      Zv_fl = RCs_INFO_ss.Zv[RC_ens];
      Za_fl = RCs_INFO_ss.Za[RC_ens];
      
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);
      
      int L= L_info.L;

      double ams1= L_info.ms_L_new;
      double ams2= L_info.ms_M_new;

      D(1);
     
      
      distr_t_list Vk_tm_distr = 1e10*Za_fl*Za_fl*(pow(qd,2))*Corr.corr_t( Vk_data_tm.col(0)[iens] , "../data/HVP_strange/Corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
      distr_t_list Vk_OS_distr = 1e10*Zv_fl*Zv_fl*(pow(qd,2))*Corr.corr_t( Vk_data_OS.col(0)[iens] , "../data/HVP_strange/Corr/Vk_OS_"+Vk_data_OS.Tag[iens]+".dat");

      distr_t_list Vk_s2_tm_distr = 1e10*Za_fl*Za_fl*(pow(qd,2))*Corr.corr_t( Vk_data_s2_tm.col(0)[iens] , "../data/HVP_strange/Corr/Vk_s2_tm_"+Vk_data_tm.Tag[iens]+".dat");
      distr_t_list Vk_s2_OS_distr = 1e10*Zv_fl*Zv_fl*(pow(qd,2))*Corr.corr_t( Vk_data_s2_OS.col(0)[iens] , "../data/HVP_strange/Corr/Vk_s2_OS_"+Vk_data_OS.Tag[iens]+".dat");
           
      distr_t_list Vk_s3_tm_distr = 1e10*Za_fl*Za_fl*(pow(qd,2))*Corr.corr_t( Vk_data_s3_tm.col(0)[iens] , "../data/HVP_strange/Corr/Vk_s3_tm_"+Vk_data_tm.Tag[iens]+".dat");
      distr_t_list Vk_s3_OS_distr = 1e10*Zv_fl*Zv_fl*(pow(qd,2))*Corr.corr_t( Vk_data_s3_OS.col(0)[iens] , "../data/HVP_strange/Corr/Vk_s3_OS_"+Vk_data_OS.Tag[iens]+".dat");
      
      
      distr_t amu_HVP_tm(UseJack,Njacks), amu_HVP_OS(UseJack,Njacks);
      distr_t amu_s2_HVP_tm(UseJack,Njacks), amu_s2_HVP_OS(UseJack,Njacks);
      distr_t amu_s3_HVP_tm(UseJack,Njacks), amu_s3_HVP_OS(UseJack,Njacks);

      distr_t W_HVP_tm(UseJack,Njacks), W_HVP_OS(UseJack,Njacks);
      distr_t W_s2_HVP_tm(UseJack,Njacks), W_s2_HVP_OS(UseJack,Njacks);
      distr_t W_s3_HVP_tm(UseJack,Njacks), W_s3_HVP_OS(UseJack,Njacks);

      distr_t SD_HVP_tm(UseJack,Njacks), SD_HVP_OS(UseJack,Njacks);
      distr_t SD_s2_HVP_tm(UseJack,Njacks), SD_s2_HVP_OS(UseJack,Njacks);
      
      D(2);

      auto K = [&](double Mv, double t, int size) -> double { return kernel_K(t, Mv);};

    
      
      distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
      double t0=  0.4*fm_to_inv_Gev;
      double t1 = 1.0*fm_to_inv_Gev;
      double Delta= 0.15*fm_to_inv_Gev;
      auto th0 = [&t0, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
      auto th1 = [&t1, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

      for(int t=1;t<Corr.Nt/2;t++) {

	amu_HVP_tm = amu_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_tm_distr.distr_list[t]*Ker.distr_list[t];
	amu_HVP_OS = amu_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_OS_distr.distr_list[t]*Ker.distr_list[t];

	amu_s2_HVP_tm = amu_s2_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_tm_distr.distr_list[t]*Ker.distr_list[t];
	amu_s2_HVP_OS = amu_s2_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_OS_distr.distr_list[t]*Ker.distr_list[t];

	amu_s3_HVP_tm = amu_s3_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_s3_tm_distr.distr_list[t]*Ker.distr_list[t];
	amu_s3_HVP_OS = amu_s3_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_s3_OS_distr.distr_list[t]*Ker.distr_list[t];


	W_HVP_tm = W_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_tm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );
	W_HVP_OS = W_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_OS_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );

	SD_HVP_tm = SD_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_tm_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr)  );
	SD_HVP_OS = SD_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_OS_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr)  );

	SD_s2_HVP_tm = SD_s2_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_tm_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr) );
	SD_s2_HVP_OS = SD_s2_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_OS_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr) ); 


	W_s2_HVP_tm = W_s2_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_tm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );
	W_s2_HVP_OS = W_s2_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_s2_OS_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );

	W_s3_HVP_tm = W_s3_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_s3_tm_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );
	W_s3_HVP_OS = W_s3_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_s3_OS_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr)    );
	

      }

      D(3);
    cout<<"#### "<<Vk_data_tm.Tag[iens]<<" ###"<<endl;
    cout<<"HVP tm: "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<" stat. "<< (amu_HVP_tm.err()*100/amu_HVP_tm.ave())<<"%"<<endl;
    cout<<"HVP OS: "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<" stat. "<< (amu_HVP_OS.err()*100/amu_HVP_OS.ave())<<"%"<<endl;
    cout<<"#######"<<endl;


    //pp interpolation

    vector<distr_t> amu_HVP_tm_int({amu_HVP_tm, amu_s2_HVP_tm});
    vector<distr_t> amu_HVP_OS_int({amu_HVP_OS, amu_s2_HVP_OS});

    vector<distr_t> amu3_HVP_tm_int({amu_HVP_tm, amu_s2_HVP_tm, amu_s3_HVP_tm});
    vector<distr_t> amu3_HVP_OS_int({amu_HVP_OS, amu_s2_HVP_OS, amu_s3_HVP_OS});
    
    vector<distr_t> W_HVP_tm_int({W_HVP_tm, W_s2_HVP_tm});
    vector<distr_t> W_HVP_OS_int({W_HVP_OS, W_s2_HVP_OS});

    vector<distr_t> SD_HVP_tm_int({SD_HVP_tm, SD_s2_HVP_tm});
    vector<distr_t> SD_HVP_OS_int({SD_HVP_OS, SD_s2_HVP_OS});

    vector<distr_t> W3_HVP_tm_int({W_HVP_tm, W_s2_HVP_tm, W_s3_HVP_tm});
    vector<distr_t> W3_HVP_OS_int({W_HVP_OS, W_s2_HVP_OS, W_s3_HVP_OS});

    vector<distr_t> MKK({MK1,MK2});

    vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});

    vector<distr_t> MMS3({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2, Get_id_jack_distr(Njacks)*ams3});

    distr_t MK_FLAG_distr= MK_FLAG*Get_id_jack_distr(Njacks);


    cout<<"Correlation between masses is (TM): "<<(W_HVP_tm%W_s2_HVP_tm)/(W_HVP_tm.err()*W_s2_HVP_tm.err())<<endl;
    cout<<"Correlation between masses is (OS): "<<(W_HVP_OS%W_s2_HVP_OS)/(W_HVP_OS.err()*W_s2_HVP_OS.err())<<endl;
    cout<<"Correlation between s1-s3  is (OS): "<<(W_HVP_OS%W_s3_HVP_OS)/(W_HVP_OS.err()*W_s3_HVP_OS.err())<<endl;

    distr_t ms_phys_TM = ms_phys + (0.0181-0.0182782)*pow(a_distr.ave()/SCALE_INFO.a_B.ave(),2);

    distr_t SD_TM_s_distr=  Obs_extrapolation_meson_mass(SD_HVP_tm_int, MMS, ms_phys ,  "../data/HVP_strange"  , "SD_TM_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t SD_OS_s_distr=  Obs_extrapolation_meson_mass(SD_HVP_OS_int, MMS, ms_phys ,  "../data/HVP_strange"  , "SD_OS_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    
    

    distr_t W_TM_s_distr=  Obs_extrapolation_meson_mass(W_HVP_tm_int, MMS, ms_phys ,  "../data/HVP_strange"  , "W_TM_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t W_TM_s_distr3=  Obs_extrapolation_meson_mass(W3_HVP_tm_int, MMS3, ms_phys ,  "../data/HVP_strange"  , "W3_TM_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t W_OS_s_distr=  Obs_extrapolation_meson_mass(W_HVP_OS_int, MMS, ms_phys ,  "../data/HVP_strange"  , "W_OS_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t W_OS_s_distr_MK=  Obs_extrapolation_meson_mass(W_HVP_OS_int, MKK, MK_FLAG_distr ,  "../data/HVP_strange"  , "W_OS_MK_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t W_OS_s_distr3=  Obs_extrapolation_meson_mass(W3_HVP_OS_int, MMS3, ms_phys ,  "../data/HVP_strange"  , "W3_OS_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
   
  
    distr_t amu_TM_s_distr=  Obs_extrapolation_meson_mass(amu_HVP_tm_int, MMS, ms_phys ,  "../data/HVP_strange"  , "amu_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t amu_TM_s_distr3=  Obs_extrapolation_meson_mass(amu3_HVP_tm_int, MMS3, ms_phys ,  "../data/HVP_strange"  , "amu3_TM_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t amu_OS_s_distr=  Obs_extrapolation_meson_mass(amu_HVP_OS_int, MMS, ms_phys ,  "../data/HVP_strange"  , "amu_OS_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );
    distr_t amu_OS_s_distr3=  Obs_extrapolation_meson_mass(amu3_HVP_OS_int, MMS3, ms_phys ,  "../data/HVP_strange"  , "amu3_TOS_extrapolation_"+Vk_data_tm.Tag[iens]+".dat",  UseJack, "SPLINE" );

    distr_t W_OS_s_distr_bis= (  (W_HVP_OS - W_s2_HVP_OS)*ms_phys -W_HVP_OS*ams2 + W_s2_HVP_OS*ams1)/(ams1-ams2);

   
    

    cout<<"######## TEST INTERPOLATION ########"<<endl;

    cout<<"Linear interpolation in ms: "<<W_OS_s_distr.ave()<<" "<<W_OS_s_distr.err()<<endl;
    cout<<"Quadratic interpolation in ms: "<<W_OS_s_distr3.ave()<<" "<<W_OS_s_distr3.err()<<endl;
    cout<<"Linear interpolation in MK: "<<W_OS_s_distr_MK.ave()<<" "<<W_OS_s_distr_MK.err()<<endl;
    

    cout<<"####################################"<<endl;

    W_TM_s.distr_list.push_back( W_TM_s_distr);
    W_OS_s.distr_list.push_back( W_OS_s_distr);
    amu_TM_s.distr_list.push_back( amu_TM_s_distr);
    amu_OS_s.distr_list.push_back( amu_OS_s_distr);
    a_distr_list.distr_list.push_back(a_distr);
         
    }

    Print_To_File({}, {a_distr_list.ave(), W_TM_s.ave(), W_TM_s.err(), amu_TM_s.ave(), amu_TM_s.err()},  "../data/HVP_strange/amu_s_TM.dat", "", "");
    Print_To_File({}, {a_distr_list.ave(), W_OS_s.ave(), W_OS_s.err(), amu_OS_s.ave(), amu_OS_s.err()},  "../data/HVP_strange/amu_s_OS.dat", "", "");
    Print_To_File({}, {a_distr_list.ave(), (0.5*(W_TM_s+W_OS_s)).ave() , (0.5*(W_TM_s+W_OS_s)).err(), (0.5*(amu_TM_s+amu_OS_s)).ave() , (0.5*(amu_TM_s+amu_OS_s)).err()} ,  "../data/HVP_strange/amu_s_ave.dat", "", ""); 

  return;

}
