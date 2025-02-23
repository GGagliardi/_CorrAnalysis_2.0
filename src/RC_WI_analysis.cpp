#include "../include/RC_WI_analysis.h"

using namespace std;

const double Njacks = 50;
const bool UseJack=true;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const double m_phi= 1.019461;
const double m_phi_err = 0.000016;


void RC_WI_analysis() { Get_RCs("ss"); }

RCs_info Get_RCs(string AN_TYPE) {

  RCs_info RC_INFO;

  double m_etas = 0.68989;
  double m_etas_err = 0.00050;

  string mix_s1="mix_s1_s1";
  string mix_s2="mix_s2_s2";

  if (AN_TYPE =="ls") {
    mix_s1="mix_l_s1";
    mix_s2="mix_l_s2";
    m_etas= 0.494600000;
  }
  else if(AN_TYPE=="ss") {
    //do nothing
  }
  else crash("AN_TYPE: "+AN_TYPE+" not yet implemented");
  


   bool Get_ASCII= false;
   
   if(Get_ASCII) {
     //read binary files
     boost::filesystem::create_directory("../RC_WI");
     

     vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
     vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});
     
     for( int it=0; it<(signed)Ens_T1.size(); it++) {
       
       vector<string> channels({mix_s1,  mix_s2});
       
      for(auto &channel : channels) {
	boost::filesystem::create_directory("../RC_WI/"+channel);
	boost::filesystem::create_directory("../RC_WI/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_A0P5", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_A0P5"});
      
      
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	
	
	for( auto &channel: channels) {
	  
	  FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	  size_t Nconfs, T, Nhits;
	  bin_read(Nconfs, stream);
	  bin_read(Nhits, stream);
	  bin_read(T, stream);
	  
	  cout<<"Nconfs: "<<Nconfs<<endl;
	  cout<<"T: "<<T<<" "<<T/2+1<<endl;
	  cout<<"Nhits: "<<Nhits<<endl;
	  for(size_t iconf=0;iconf<Nconfs;iconf++) {
	    vector<double> C(T/2+1);
	    for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	    boost::filesystem::create_directory("../RC_WI/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	    ofstream PrintCorr("../RC_WI/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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


   cout<<"Producing ll correlators..."<<endl;
   bool Get_ASCII_l= false;
   
   if(Get_ASCII_l) {
     //read binary files
     boost::filesystem::create_directory("../RC_WI_light");
     

     vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
     vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});
     
     for( int it=0; it<(signed)Ens_T1.size(); it++) {
       
       vector<string> channels({"ll"});
       if(Ens_T1[it] == "E.44.112") channels= {"mix_l_l"};
       
       
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_A0P5", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_A0P5"});
      
      
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	
	
	for( auto &channel: channels) {

	  string c= "ll";
	  boost::filesystem::create_directory("../RC_WI_light/"+c);
	  boost::filesystem::create_directory("../RC_WI_light/"+c+"/"+Ens_TT1[it]);
	  
	  FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	  size_t Nconfs, T, Nhits;
	  bin_read(Nconfs, stream);
	  bin_read(Nhits, stream);
	  bin_read(T, stream);
	  
	  cout<<"Nconfs: "<<Nconfs<<endl;
	  cout<<"T: "<<T<<" "<<T/2+1<<endl;
	  cout<<"Nhits: "<<Nhits<<endl;
	  for(size_t iconf=0;iconf<Nconfs;iconf++) {
	    vector<double> C(T/2+1);
	    for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	    boost::filesystem::create_directory("../RC_WI_light/"+c+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	    ofstream PrintCorr("../RC_WI_light/"+c+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+c+"_"+Corr_tags[id]);
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
    

   data_t corr_P5P5_strange, corr_P5P5_OS_strange, corr_P5P5_strange_heavy, corr_P5P5_OS_strange_heavy;
   data_t corr_A0P5_strange, corr_A0P5_OS_strange, corr_A0P5_strange_heavy, corr_A0P5_OS_strange_heavy;

   //phi

   data_t corr_VKVK_strange, corr_VKVK_strange_heavy;


   //P5P5
   corr_P5P5_strange.Read("../RC_WI/"+mix_s1+"", "mes_contr_"+mix_s1+"_TM_P5P5", "P5P5", Sort_easy);

   corr_P5P5_OS_strange.Read("../RC_WI/"+mix_s1+"", "mes_contr_"+mix_s1+"_OS_P5P5", "P5P5", Sort_easy);
   corr_P5P5_strange_heavy.Read("../RC_WI/"+mix_s2+"", "mes_contr_"+mix_s2+"_TM_P5P5", "P5P5", Sort_easy);
   corr_P5P5_OS_strange_heavy.Read("../RC_WI/"+mix_s2+"", "mes_contr_"+mix_s2+"_OS_P5P5","P5P5",Sort_easy);
   corr_VKVK_strange.Read("../RC_WI/"+mix_s1+"", "mes_contr_"+mix_s1+"_TM_VKVK", "VKVK", Sort_easy);
   corr_VKVK_strange_heavy.Read("../RC_WI/"+mix_s2+"", "mes_contr_"+mix_s2+"_TM_VKVK", "VKVK", Sort_easy);
   //A0P5
   corr_A0P5_strange.Read("../RC_WI/"+mix_s1+"", "mes_contr_"+mix_s1+"_TM_A0P5", "A0P5", Sort_easy);
   corr_A0P5_OS_strange.Read("../RC_WI/"+mix_s1+"", "mes_contr_"+mix_s1+"_OS_A0P5", "A0P5", Sort_easy);
   corr_A0P5_strange_heavy.Read("../RC_WI/"+mix_s2+"", "mes_contr_"+mix_s2+"_TM_A0P5", "A0P5", Sort_easy);
   corr_A0P5_OS_strange_heavy.Read("../RC_WI/"+mix_s2+"", "mes_contr_"+mix_s2+"_OS_A0P5", "A0P5", Sort_easy);


   data_t corr_P5P5_light, corr_P5P5_OS_light;
   data_t corr_A0P5_light, corr_A0P5_OS_light;

   corr_P5P5_light.Read("../RC_WI_light/ll", "mes_contr_ll_TM_P5P5", "P5P5", Sort_easy);
   corr_P5P5_OS_light.Read("../RC_WI_light/ll", "mes_contr_ll_OS_P5P5", "P5P5", Sort_easy);

   corr_A0P5_light.Read("../RC_WI_light/ll", "mes_contr_ll_TM_A0P5", "A0P5", Sort_easy);
   corr_A0P5_OS_light.Read("../RC_WI_light/ll", "mes_contr_ll_OS_A0P5", "A0P5", Sort_easy);




   //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_E(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_E_ave, a_E_err;

  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;

  for(int ijack=0;ijack<Njacks;ijack++) {

  a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
  a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
  a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
  a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
  a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
  }


  boost::filesystem::create_directory("../data/RC_WI");
  boost::filesystem::create_directory("../data/RC_WI/masses");
  boost::filesystem::create_directory("../data/RC_WI/ZA");
  boost::filesystem::create_directory("../data/RC_WI/ZV");
  boost::filesystem::create_directory("../data/RC_WI/corr");
  

  int Nens= corr_P5P5_strange.size;

  distr_t_list Za_list(UseJack), Zv_list(UseJack), Zp_ov_Zs_list(UseJack);


  for(int iens=0;iens<Nens;iens++ ) {

    string Ens_Tag= corr_P5P5_strange.Tag[iens];
    cout<<"Analyzing ensemble: "<<Ens_Tag<<endl;

    assert( Ens_Tag == corr_P5P5_light.Tag[iens]);

    CorrAnalysis Corr(UseJack, Njacks,100);
    Corr.Nt = corr_P5P5_strange.nrows[iens];


    distr_t a_distr;
    if(Ens_Tag.substr(1,1)=="B") {  a_distr=a_B;
      if(Ens_Tag=="cB211b.072.64") {
	Corr.Tmin= 40; //40;
	Corr.Tmax= 57; //57;
      }
      else { // Ens_Tag = cB211b.072.96
	Corr.Tmin=40;
	Corr.Tmax=87;
      }

    }
    else if(Ens_Tag.substr(1,1)=="C") {  a_distr=a_C;
      if(Ens_Tag=="cC211a.06.80")   {
	Corr.Tmin=44;
	Corr.Tmax=70;
      }
      else { // Ens_Tag = cC211a.06.112
	Corr.Tmin=40; Corr.Tmax=80;
      }
    }
    else if(Ens_Tag.substr(1,1)=="D") {a_distr=a_D; Corr.Tmin=48; Corr.Tmax=80; }
    else if(Ens_Tag.substr(1,1)=="E") {a_distr=a_E; Corr.Tmin=62; Corr.Tmax=100;}
    else crash("Ensemble not found");

    distr_t_list corr_s_P_L_tm = Corr.corr_t( corr_P5P5_strange.col(0)[iens], "");
    distr_t_list corr_s_P_L_OS = Corr.corr_t( corr_P5P5_OS_strange.col(0)[iens], "");

    distr_t_list corr_s_P_M_tm = Corr.corr_t( corr_P5P5_strange_heavy.col(0)[iens], "");
    distr_t_list corr_s_P_M_OS = Corr.corr_t( corr_P5P5_OS_strange_heavy.col(0)[iens], "");
  
    distr_t_list corr_s_AP_L_tm = Corr.corr_t( corr_A0P5_strange.col(0)[iens], "../data/RC_WI/corr/A0P5_L_tm_"+Ens_Tag);
    distr_t_list corr_s_AP_L_OS = Corr.corr_t( corr_A0P5_OS_strange.col(0)[iens], "");

    distr_t_list corr_s_AP_M_tm = Corr.corr_t( corr_A0P5_strange_heavy.col(0)[iens], "../data/RC_WI/corr/A0P5_M_tm_"+Ens_Tag);
    distr_t_list corr_s_AP_M_OS = Corr.corr_t( corr_A0P5_OS_strange_heavy.col(0)[iens], "");


    distr_t_list corr_s_V_L_tm = Corr.corr_t( corr_VKVK_strange.col(0)[iens], "");
    distr_t_list corr_s_V_M_tm = Corr.corr_t( corr_VKVK_strange_heavy.col(0)[iens], "");

    //light
    distr_t_list corr_light_P_tm= Corr.corr_t( corr_P5P5_light.col(0)[iens], "");
    distr_t_list corr_light_P_OS= Corr.corr_t( corr_P5P5_OS_light.col(0)[iens], "");
    distr_t_list corr_light_AP_tm= Corr.corr_t( corr_A0P5_light.col(0)[iens], "");
    distr_t_list corr_light_AP_OS= Corr.corr_t( corr_A0P5_OS_light.col(0)[iens], "");
    
    

    LatticeInfo L_info;
     
    L_info.LatInfo_new_ens(Ens_Tag);
     
    double aml= L_info.ml;
    double ams1= L_info.ms_L_new;
    double ams2= L_info.ms_M_new;


    
    distr_t Meta_tm_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_tm, "../data/RC_WI/masses/eta_L_tm_"+Ens_Tag));
    distr_t Meta_tm_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_tm, "../data/RC_WI/masses/eta_M_tm_"+Ens_Tag));
    distr_t Meta_OS_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_L_OS, "../data/RC_WI/masses//eta_L_OS_"+Ens_Tag));
    distr_t Meta_OS_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_P_M_OS, "../data/RC_WI/masses/eta_M_OS_"+Ens_Tag));


    int Tmin_old= Corr.Tmin; int Tmax_old=Corr.Tmax;

    
    if(Ens_Tag.substr(1,1)=="B") { Corr.Tmin=27; Corr.Tmax=37;}
    else if(Ens_Tag.substr(1,1)=="C") { Corr.Tmin=43 ; Corr.Tmax = 60;}
    else if(Ens_Tag.substr(1,1)=="D") { Corr.Tmin= 57 ; Corr.Tmax = 78;}
    else if(Ens_Tag.substr(1,1)=="E") { Corr.Tmin=40 ; Corr.Tmax =80 ;}
    else crash("Ens not found");
    
    distr_t Mpi_OS =  Corr.Fit_distr( Corr.effective_mass_t(corr_light_P_OS, "../data/RC_WI/masses/Mpi_OS_"+Ens_Tag));
    distr_t Mpi_tm =  Corr.Fit_distr( Corr.effective_mass_t(corr_light_P_tm, "../data/RC_WI/masses/Mpi_tm_"+Ens_Tag));

    Corr.Tmin=Tmin_old; Corr.Tmax=Tmax_old;

    distr_t Mphi_tm_L = Corr.Fit_distr( Corr.effective_mass_t(corr_s_V_L_tm, "../data/RC_WI/masses/phi_L_tm_"+Ens_Tag));
    distr_t Mphi_tm_M = Corr.Fit_distr( Corr.effective_mass_t(corr_s_V_M_tm, "../data/RC_WI/masses/phi_M_tm_"+Ens_Tag));


    cout<<"Meta(tm,L): "<<Meta_tm_L.ave()<<" "<<Meta_tm_L.err()<<endl;
    cout<<"Meta(tm,M): "<<Meta_tm_M.ave()<<" "<<Meta_tm_M.err()<<endl;

    cout<<"Meta(OS,L): "<<Meta_OS_L.ave()<<" "<<Meta_OS_L.err()<<endl;
    cout<<"Meta(OS,M): "<<Meta_OS_M.ave()<<" "<<Meta_OS_M.err()<<endl;

    cout<<"Mphi(tm,L): "<<Mphi_tm_L.ave()<<" "<<Mphi_tm_L.err()<<endl;
    cout<<"Mphi(tm,M): "<<Mphi_tm_M.ave()<<" "<<Mphi_tm_M.err()<<endl;

    double F1= ams1; double F2=ams2;

    if(AN_TYPE=="ls") { F1 = 0.5*(aml+ams1) ; F2= 0.5*(aml+ams2);}
     
    distr_t_list RV_L = 2*F1*corr_s_P_L_tm/(distr_t_list::derivative(corr_s_AP_L_tm, 0));  
    distr_t_list RV_M = 2*F2*corr_s_P_M_tm/(distr_t_list::derivative(corr_s_AP_M_tm, 0));
    distr_t_list RA_L = 2*F1*corr_s_P_L_OS/(distr_t_list::derivative(corr_s_AP_L_OS, 0)); 
    distr_t_list RA_M = 2*F2*corr_s_P_M_OS/(distr_t_list::derivative(corr_s_AP_M_OS, 0));
    RA_L = RA_L*(Meta_OS_L/Meta_tm_L)*(SINH_D(Meta_OS_L)/SINH_D(Meta_tm_L))*Corr.matrix_element_t(corr_s_P_L_tm, "")/Corr.matrix_element_t(corr_s_P_L_OS, "");
    RA_M = RA_M*(Meta_OS_M/Meta_tm_M)*(SINH_D(Meta_OS_M)/SINH_D(Meta_tm_M))*Corr.matrix_element_t(corr_s_P_M_tm, "")/Corr.matrix_element_t(corr_s_P_M_OS, "");
    distr_t_list Zp_ov_Zs_L_distr = 1.0/(Corr.matrix_element_t(corr_s_P_L_tm, "")/Corr.matrix_element_t(corr_s_P_L_OS, "") );
    distr_t_list Zp_ov_Zs_M_distr = 1.0/(Corr.matrix_element_t(corr_s_P_M_tm, "")/Corr.matrix_element_t(corr_s_P_M_OS, "") );
    //light

    Tmin_old= Corr.Tmin; Tmax_old=Corr.Tmax;

    
    if(Ens_Tag.substr(1,1)=="B") { Corr.Tmin=27; Corr.Tmax=37;}
    else if(Ens_Tag.substr(1,1)=="C") { Corr.Tmin=43 ; Corr.Tmax = 60;}
    else if(Ens_Tag.substr(1,1)=="D") { Corr.Tmin= 57 ; Corr.Tmax = 78;}
    else if(Ens_Tag.substr(1,1)=="E") { Corr.Tmin=40 ; Corr.Tmax =80 ;}
    else crash("Ens not found");
    
    
    distr_t_list RV_ll = 2*aml*corr_light_P_tm/(distr_t_list::derivative(corr_light_AP_tm, 0));  
    distr_t_list RA_ll = 2*aml*corr_light_P_OS/(distr_t_list::derivative(corr_light_AP_OS, 0));
    RA_ll = RA_ll*(Mpi_OS/Mpi_tm)*(SINH_D(Mpi_OS)/SINH_D(Mpi_tm))*Corr.matrix_element_t(corr_light_P_tm, "")/Corr.matrix_element_t(corr_light_P_OS, "");
    distr_t_list Zp_ov_Zs_ll_distr = 1.0/(Corr.matrix_element_t(corr_light_P_tm, "")/Corr.matrix_element_t(corr_light_P_OS, "") );
    
    Corr.Tmin=Tmin_old; Corr.Tmax=Tmax_old;
    
    
    Print_To_File({}, {RV_L.ave(), RV_L.err(), RV_M.ave(), RV_M.err()} , "../data/RC_WI/ZV/RV_"+Ens_Tag+".dat", "", "");
    Print_To_File({}, {RA_L.ave(), RA_L.err(), RA_M.ave(), RA_M.err()} , "../data/RC_WI/ZA/RA_"+Ens_Tag+".dat", "", "");
    Print_To_File({}, {Zp_ov_Zs_L_distr.ave(), Zp_ov_Zs_L_distr.err(), Zp_ov_Zs_M_distr.ave(), Zp_ov_Zs_M_distr.err()} , "../data/RC_WI/Za/ZpZs_"+Ens_Tag+".dat", "", "");

    Print_To_File({}, {RV_ll.ave(), RV_ll.err()}, "../data/RC_WI/ZV/RV_ll"+Ens_Tag+".dat", "","");
    Print_To_File({}, {RA_ll.ave(), RA_ll.err()}, "../data/RC_WI/ZV/RV_ll"+Ens_Tag+".dat", "","");
    Print_To_File({}, {Zp_ov_Zs_ll_distr.ave(), Zp_ov_Zs_ll_distr.err()}, "../data/RC_WI/ZV/ZpZs_ll"+Ens_Tag+".dat", "","");

    //set time interval for Z factors
    if(Ens_Tag=="cB211b.072.64") {  Corr.Tmin= 40; Corr.Tmax= 55 ; }
    else if(Ens_Tag=="cB211b.072.96") {  Corr.Tmin= 41; Corr.Tmax= 78 ; }
    else if(Ens_Tag=="cC211a.06.80") {  Corr.Tmin= 40; Corr.Tmax= 70 ; }
    else if(Ens_Tag=="cC211a.06.112") {  Corr.Tmin= 45; Corr.Tmax= 90 ; }
    else if(Ens_Tag=="cD211a.054.96") {  Corr.Tmin= 48; Corr.Tmax= 80 ; }
    else if(Ens_Tag=="cE211a.044.112") {  Corr.Tmin= 55; Corr.Tmax= 100 ; }
    else crash(" Ensemble: "+Ens_Tag+" not found!");


    if(AN_TYPE=="ls") {

      if( Ens_Tag =="cB211b.072.96")     { Corr.Tmin=38; Corr.Tmax=59; }
      else if(Ens_Tag =="cB211b.072.64") { Corr.Tmin=45; Corr.Tmax=59; }
      else if(Ens_Tag.substr(1,1)=="C")  { Corr.Tmin=36; Corr.Tmax=63; }
      else if(Ens_Tag.substr(1,1)=="D")  { Corr.Tmin=37; Corr.Tmax=69; }
      else if(Ens_Tag.substr(1,1)=="E")  { Corr.Tmin=60; Corr.Tmax=85; }
      else crash("Cannot recognize the ensemble: "+Ens_Tag+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");
      
      
      
    }


    distr_t Zv_L = Corr.Fit_distr(RV_L);
    distr_t Zv_M = Corr.Fit_distr(RV_M);
    distr_t Za_L = Corr.Fit_distr(RA_L);
    distr_t Za_M = Corr.Fit_distr(RA_M);
    distr_t Zp_ov_Zs_L = Corr.Fit_distr(Zp_ov_Zs_L_distr);
    distr_t Zp_ov_Zs_M = Corr.Fit_distr(Zp_ov_Zs_M_distr);

    distr_t Zv_ll = Corr.Fit_distr(RV_ll);
    distr_t Za_ll = Corr.Fit_distr(RA_ll);
    distr_t Zp_ov_Zs_ll = Corr.Fit_distr(Zp_ov_Zs_ll_distr);

   
    

  
    vector<distr_t> M2etas_fit({Meta_tm_L*Meta_tm_L/(a_distr*a_distr), Meta_tm_M*Meta_tm_M/(a_distr*a_distr)});
    vector<distr_t> Metas_fit({Meta_tm_L/(a_distr), Meta_tm_M/(a_distr)});

    vector<distr_t> Mphi_fit({Mphi_tm_L*Mphi_tm_L/(a_distr*a_distr), Mphi_tm_M*Mphi_tm_M/(a_distr*a_distr)});

    distr_t m_etas_phys_distr, m_phi_phys_distr;
    
  
    for(int ijack=0;ijack<Njacks;ijack++) m_etas_phys_distr.distr.push_back( m_etas + GM()*1e-10*m_etas_err/sqrt(Njacks-1.0));
    for(int ijack=0;ijack<Njacks;ijack++) m_phi_phys_distr.distr.push_back( m_phi + GM()*m_phi_err/sqrt(Njacks-1.0));
  
    //Generate fake ms_distr
    distr_t ms_light_distr;
    distr_t ms_heavy_distr;
    for(int ijack=0;ijack<Njacks;ijack++) ms_light_distr.distr.push_back( ams1 );
    for(int ijack=0;ijack<Njacks;ijack++) ms_heavy_distr.distr.push_back( ams2 );
  
    //estrapolate ms phys
    vector<distr_t> ms_list( {ms_light_distr, ms_heavy_distr});
    distr_t ms_phys_extr = Obs_extrapolation_meson_mass(ms_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr,  "../data/RC_WI", "ms_extrapolation_etas_"+Ens_Tag, UseJack, "SPLINE");
    distr_t ms_phys_extr_phi = Obs_extrapolation_meson_mass(ms_list, Mphi_fit, m_phi_phys_distr*m_phi_phys_distr,  "../data/RC_WI", "ms_extrapolation_phi_"+Ens_Tag, UseJack, "SPLINE");

    cout<<"##### ms  ["<<Ens_Tag<<"] ####"<<endl;
    cout<<"ms1: "<<ams1<<" ms2: "<<ams2<<endl;
    cout<<"ms(etas): "<<ms_phys_extr.ave()<<" +- "<<ms_phys_extr.err()<<endl;
    cout<<"ms(phi): "<<ms_phys_extr_phi.ave()<<" +- "<<ms_phys_extr_phi.err()<<endl;

    vector<distr_t> Za_hadr_list, Zv_hadr_list, Zp_ov_Zs_hadr_list;
    Za_hadr_list = {Za_L, Za_M};
    Zv_hadr_list = {Zv_L, Zv_M};
    Zp_ov_Zs_hadr_list ={ Zp_ov_Zs_L, Zp_ov_Zs_M};

    //ms(phi): 0.0168083 +- 0.000960077
    //ms(phi): 0.0168319 +- 0.000936623


    distr_t Za = Obs_extrapolation_meson_mass( Za_hadr_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr, "../data/RC_WI", "Za_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");
    distr_t Zv = Obs_extrapolation_meson_mass( Zv_hadr_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr, "../data/RC_WI", "Zv_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");
    distr_t Zp_ov_Zs = Obs_extrapolation_meson_mass( Zp_ov_Zs_hadr_list, M2etas_fit, m_etas_phys_distr*m_etas_phys_distr, "../data/RC_WI", "ZpZs_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");

    distr_t Za_lin = Obs_extrapolation_meson_mass( Za_hadr_list, Metas_fit, m_etas_phys_distr, "../data/RC_WI", "Za_lin_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");
    distr_t Zv_lin = Obs_extrapolation_meson_mass( Zv_hadr_list, Metas_fit, m_etas_phys_distr, "../data/RC_WI", "Zv_lin_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");
    distr_t Zp_ov_Zs_lin = Obs_extrapolation_meson_mass( Zp_ov_Zs_hadr_list, Metas_fit, m_etas_phys_distr, "../data/RC_WI", "ZpZs_lin_extr_quark_mass_"+Ens_Tag, UseJack, "SPLINE");

    cout<<"#### RCs ["<<Ens_Tag<<"] ####"<<endl;
    cout.precision(8);
    cout<<"Za: "<<Za.ave()<<" "<<Za.err()<<" lin: "<<Za_lin.ave()<<" +- "<<Za_lin.err()<<endl;
    cout<<"Zv: "<<Zv.ave()<<" "<<Zv.err()<<" lin: "<<Zv_lin.ave()<<" +- "<<Zv_lin.err()<<endl;
    cout<<"Zp/Zs: "<<Zp_ov_Zs.ave()<<" "<<Zp_ov_Zs.err()<<" lin: "<<Zp_ov_Zs_lin.ave()<<" +- "<<Zp_ov_Zs_lin.err()<<endl;
    cout<<"Zs/Zp: "<<(1.0/Zp_ov_Zs).ave()<<" "<<(1.0/Zp_ov_Zs).err()<<" lin: "<<(1.0/Zp_ov_Zs_lin).ave()<<" +- "<<(1.0/Zp_ov_Zs_lin).err()<<endl;
    cout<<"From ll correlators: "<<endl;
    cout<<"Za: "<<Za_ll.ave()<<" +- "<<Za_ll.err()<<endl;
    cout<<"Zv: "<<Zv_ll.ave()<<" +- "<<Zv_ll.err()<<endl;
    cout<<"Zp/Zs: "<<Zp_ov_Zs_ll.ave()<<" +- "<<Zp_ov_Zs_ll.err()<<endl;
    cout<<"Zs/Zp: "<<(1.0/Zp_ov_Zs_ll).ave()<<" +- "<<(1.0/Zp_ov_Zs_ll).err()<<endl;


    cout<<"a2: "<<(a_distr*a_distr).ave()<<" "<<(Zp_ov_Zs/Zp_ov_Zs_ll -1.0).ave()<<" " <<(Zp_ov_Zs/Zp_ov_Zs_ll -1.0).err()<<endl;


    RC_INFO.Ens.push_back( Ens_Tag);
    Za_list.distr_list.push_back(Za);
    Zv_list.distr_list.push_back(Zv);
    Zp_ov_Zs_list.distr_list.push_back(Zp_ov_Zs);


  }


  RC_INFO.Za = Za_list;
  RC_INFO.Zv = Zv_list;
  RC_INFO.Zp_ov_Zs = Zp_ov_Zs_list;
  


  return RC_INFO;
}
