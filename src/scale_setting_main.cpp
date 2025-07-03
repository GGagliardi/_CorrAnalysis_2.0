#include "../include/scale_setting_main.h"
#include "Corr_analysis.h"
#include "Meson_mass_extrapolation.h"
#include "numerics.h"
#include "stat.h"



using namespace std;

//set constants

const bool UseJack=1;
const int Njacks=50; //50
const int Nboots=200;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const bool Use_three_finest_in_scale_setting_fp = false;

// FLAG VALUES FOR INPUT
const double MP_FLAG= 0.135;
const double MK_FLAG = 0.494600000;
const double MDs_FLAG = 1.967;
const double Metac_FLAG=2.980;

void Get_scale_setting() {

  scale_setting_info SC_INFO= Get_scale_setting_info();

  return;

}

scale_setting_info Get_scale_setting_info() {

  omp_set_num_threads(1);

  scale_setting_info SCALE_INFO;


   bool Get_ASCII= false;

    if(Get_ASCII) {
    //read binary files
    //boost::filesystem::create_directory("../E112_bin");
    

      vector<string> Ens_T1({"B.72.64", "B.72.96", "C.06.80","D.54.96","E.44.112"});
      vector<string> Ens_TT1({"cB211b.072.64", "cB211b.072.96", "cC211a.06.80",  "cD211a.054.96",  "cE211a.044.112"});


    boost::filesystem::create_directory("../corr_scale_setting/light");

    for( int it=0; it<(signed)Ens_T1.size(); it++) {


      
      
      boost::filesystem::create_directory("../corr_scale_setting/light/"+Ens_TT1[it]);
      //read binary
          
      cout<<"Analyzing Ens: "<<Ens_T1[it]<<endl<<flush;

      
      FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+((Ens_TT1[it]=="cE211a.044.112")?"mix_l_l_TM_P5P5":"ll_TM_P5P5")).c_str(), "rb");
      size_t Nconfs, T, Nhits;
      bin_read(Nconfs, stream);
      bin_read(Nhits, stream);
      bin_read(T, stream);

      cout<<"Nconfs: "<<Nconfs<<endl<<flush;
      cout<<"T: "<<T<<" "<<T/2+1<<endl<<flush;
      cout<<"Nhits: "<<Nhits<<endl<<flush;
      for(size_t iconf=0;iconf<Nconfs;iconf++) {
	
	string conf_id = to_string(iconf);
	if(conf_id.length() == 1) conf_id = "000"+conf_id;
	if(conf_id.length() == 2) conf_id = "00"+conf_id;
	if(conf_id.length() == 3) conf_id = "0"+conf_id;
	if(conf_id.length() > 4 ) crash("conf_id.size(): "+to_string(conf_id.size())+" conf: "+conf_id);

	conf_id = conf_id+"_r0";
	vector<double> C(T/2+1);
	for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	boost::filesystem::create_directory("../corr_scale_setting/light/"+Ens_TT1[it]+"/"+conf_id);
	ofstream PrintCorr("../corr_scale_setting/light/"+Ens_TT1[it]+"/"+conf_id+"/mes_contr_2pts_ll_1");
	PrintCorr.precision(16);
	PrintCorr<<"# "<<"P5P5"<<endl;
	for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<"   "<<"0.000000"<<endl;
	for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<"   "<<"0.000000"<<endl;
	PrintCorr.close();
	
      }
      
      fclose(stream);
      
	
	
	cout<<"Ens: "<<Ens_T1[it]<<" finished!"<<endl;
    }
    }

 


  


  //create directories g-2
  boost::filesystem::create_directory("../data/scale_setting");
  boost::filesystem::create_directory("../data/scale_setting/Mp");
  boost::filesystem::create_directory("../data/scale_setting/fp");
 

  
  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################

    


  
 
  //init Gaussian number generator
  GaussianMersenne GM(943832);


 
  
  data_t  pt2_pion, pt2_pion_B25, pt2_pion_A, pt2_pion_Z56;
    


  //Custom sorting for V_light to account for the two replica r0 and r1, ...., rn

 
  
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


 
   
  pt2_pion.Read("../corr_scale_setting/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_B25.Read("../corr_scale_setting/B25_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_A.Read("../corr_scale_setting/A_light_ens", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  pt2_pion_Z56.Read("../Aprime/mix_l2_l2", "mes_contr_mix_l2_l2_TM_P5P5", "P5P5", Sort_light_confs);

  
 

  //##################################################################################
  //###################    SCALE SETTING  USING fp, Mp    ########################


 
  distr_t a_from_fp_A(UseJack), a_from_fp_B(UseJack), a_from_fp_C(UseJack), a_from_fp_D(UseJack), a_from_fp_E(UseJack), a_from_fp_Z(UseJack);


  distr_t_list Mpi_scale_setting(UseJack), fpi_scale_setting(UseJack);
  distr_t_list Mpi_scale_setting_phys_point_ens(UseJack), fpi_scale_setting_phys_point_ens(UseJack);
  distr_t_list Mpi_scale_setting_Bens(UseJack), fpi_scale_setting_Bens(UseJack);
  distr_t_list Mpi_scale_setting_Aens(UseJack), fpi_scale_setting_Aens(UseJack);
  vector<string> Ensemble_scale_setting_tag_list;
  vector<string> Ensemble_phys_point_tag_list;
  vector<string> Ensemble_B_tag_list;
  vector<string> Ensemble_A_tag_list;
  Vfloat L_scale_setting_list;
  Vfloat L_phys_point, L_B_ens, L_A_ens;
  //physical point ensemble
  for(int iens=0;iens< pt2_pion.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion.nrows[iens];
      //Analyze correlators
      if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=40; Corr.Tmax=80;}
      else if(pt2_pion.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=35; Corr.Tmax=55;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="C") {Corr.Tmin=43; Corr.Tmax=65;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="D") {Corr.Tmin=45; Corr.Tmax=80;}
      else if(pt2_pion.Tag[iens].substr(1,1)=="E") {Corr.Tmin=60; Corr.Tmax=90;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion.Tag[iens]+".dat");


      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      Mpi_scale_setting_phys_point_ens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_phys_point_ens.distr_list.push_back(fpi_fit);
      Ensemble_phys_point_tag_list.push_back(pt2_pion.Tag[iens]);
      L_phys_point.push_back( L_info.L);
      if(pt2_pion.Tag[iens].substr(1,1) == "B") {
	Mpi_scale_setting_Bens.distr_list.push_back( Mpi_fit);
	fpi_scale_setting_Bens.distr_list.push_back( fpi_fit);
	Ensemble_B_tag_list.push_back(pt2_pion.Tag[iens]);
	L_B_ens.push_back(L_info.L);
      }
      
    }


    //B25 ensembles
    for(int iens=0;iens< pt2_pion_B25.size; iens++) {

      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_B25.Tag[iens]);
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_B25.nrows[iens];
      if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.24") {Corr.Tmin=15; Corr.Tmax=20;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.32") {Corr.Tmin=23; Corr.Tmax=30;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.25.48") {Corr.Tmin=25; Corr.Tmax=44;}
      else if(pt2_pion_B25.Tag[iens].substr(1,11)=="B211a.14.64") {Corr.Tmin = 39; Corr.Tmax=52;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_B25.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_B25.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_B25.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion_B25.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Ensemble_scale_setting_tag_list.push_back(pt2_pion_B25.Tag[iens]);
      Mpi_scale_setting_Bens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Bens.distr_list.push_back(fpi_fit);
      Ensemble_B_tag_list.push_back( pt2_pion_B25.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_B_ens.push_back( L_info.L);
     
    }

  

    //A ensembles
    for(int iens=0;iens< pt2_pion_A.size; iens++) {

      if(pt2_pion_A.Tag[iens].substr(1,1) == "A") {
      //Lattice info
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_A.Tag[iens]);
      double aml_A= L_info.ml;
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = pt2_pion_A.nrows[iens];
      if(pt2_pion_A.Tag[iens] == "cA211a.40.24") {Corr.Tmin=15; Corr.Tmax=21;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.53.24") {Corr.Tmin=15; Corr.Tmax=22;}
      else if(pt2_pion_A.Tag[iens] == "cA211ab.30.32") {Corr.Tmin=16; Corr.Tmax=29;}
      else if(pt2_pion_A.Tag[iens] == "cA211a.12.48") {Corr.Tmin=21; Corr.Tmax=34;}
      else crash("In scale setting analysis cannot find Tmin,Tmax for ensemble: "+pt2_pion_A.Tag[iens]);
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_A.col(0)[iens], "");
      distr_t_list Mpi_eff_distr = Corr.effective_mass_t(pt2_pion_A.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion_A.Tag[iens]+".dat");
      distr_t Mpi_fit = Corr.Fit_distr(Mpi_eff_distr);
      distr_t fpi_fit = Corr.Fit_distr(fpi_eff_distr);
      Mpi_scale_setting.distr_list.push_back(Mpi_fit);
      fpi_scale_setting.distr_list.push_back(fpi_fit);
      Mpi_scale_setting_Aens.distr_list.push_back(Mpi_fit);
      fpi_scale_setting_Aens.distr_list.push_back(fpi_fit);
      Ensemble_A_tag_list.push_back( pt2_pion_A.Tag[iens]);
      L_scale_setting_list.push_back( L_info.L);
      L_A_ens.push_back( L_info.L);
      cout<<pt2_pion_A.Tag[iens]<<" am_ud: "<<aml_A<<" aMpi: "<<Mpi_fit.ave()<<" +- "<<Mpi_fit.err()<<endl;
      }
      if(pt2_pion_A.Tag[iens].substr(1,1) == "A") {
	Ensemble_scale_setting_tag_list.push_back(pt2_pion_A.Tag[iens]);
      }

     
    }

    //Z56
    distr_t Mpi_Z56(UseJack), fpi_Z56(UseJack);

    for(int iens=0;iens< pt2_pion_Z56.size;iens++) {

      assert(iens==0);
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(pt2_pion_Z56.Tag[iens]);
      cout<<pt2_pion_Z56.Tag[iens]<<endl;
      double am= L_info.ml;
      cout<<"ml: "<<am<<endl;
      CorrAnalysis Corr(UseJack, Njacks, Nboots);
      Corr.Nt= pt2_pion_Z56.nrows[iens];
      cout<<"T:"<<Corr.Nt<<endl;
      Corr.Tmin=18; Corr.Tmax=30;
      distr_t_list pion_corr = Corr.corr_t(pt2_pion_Z56.col(0)[iens], "");
      distr_t_list Mpi_eff_distr =  Corr.effective_mass_t(pt2_pion_Z56.col(0)[iens], "../data/scale_setting/Mp/Mpi_"+pt2_pion_Z56.Tag[iens]+".dat");
      distr_t_list fpi_eff_distr = Corr.decay_constant_t(pow(2.0*am,2)*pion_corr, "../data/scale_setting/fp/fpi_"+pt2_pion_Z56.Tag[iens]+".dat");
      Mpi_Z56 = Corr.Fit_distr(Mpi_eff_distr);
      fpi_Z56 = Corr.Fit_distr(fpi_eff_distr);
      
    }
  
 
  
  
    Determine_scale_from_fp_FLAG(Mpi_scale_setting_phys_point_ens, fpi_scale_setting_phys_point_ens, Mpi_scale_setting_Bens, fpi_scale_setting_Bens, Mpi_scale_setting_Aens, fpi_scale_setting_Aens, L_phys_point, L_B_ens, L_A_ens, Ensemble_phys_point_tag_list, Ensemble_B_tag_list, Ensemble_A_tag_list, Mpi_Z56, fpi_Z56,  a_from_fp_A, a_from_fp_B, a_from_fp_C, a_from_fp_D, a_from_fp_E, a_from_fp_Z, UseJack, Use_three_finest_in_scale_setting_fp);


    //print pion masses
    cout<<"###############  PRINTING PION MASSES: "<<endl;
    for(int i=0;i<(signed)Ensemble_phys_point_tag_list.size();i++) {

      string T= Ensemble_phys_point_tag_list[i];
      distr_t M= Mpi_scale_setting_phys_point_ens.distr_list[i];
      distr_t F = fpi_scale_setting_phys_point_ens.distr_list[i];
      SCALE_INFO.Ens_l.push_back(T);
      distr_t a(UseJack);
      if(T.substr(1,1) == "B") { a= a_from_fp_B;}
      else if(T.substr(1,1) == "C" ) { a= a_from_fp_C;}
      else if(T.substr(1,1) == "D" ) { a= a_from_fp_D;}
      else if(T.substr(1,1) == "E" ) { a= a_from_fp_E;}
      cout<<"Mpi("<<T<<") : "<<(M/a).ave()<<" +- "<<(M/a).err()<<" [lu]: "<<M.ave()<<" +- "<<M.err()<<endl;
      
    }

    //push_back FPI and MPI
    SCALE_INFO.Mpi = Mpi_scale_setting_phys_point_ens;
    SCALE_INFO.fpi = fpi_scale_setting_phys_point_ens;


    // get lattice spacings [ Gev^-1 ]

    distr_t a_A = a_from_fp_A*1.00;
    distr_t a_B = a_from_fp_B*1.00;
    distr_t a_C = a_from_fp_C*1.00;
    distr_t a_D = a_from_fp_D*1.00;
    distr_t a_E = a_from_fp_E*1.00;
    

    //push-lattice spacing
    SCALE_INFO.a_A= a_A;
    SCALE_INFO.a_B= a_B;
    SCALE_INFO.a_C= a_C;
    SCALE_INFO.a_D= a_D;
    SCALE_INFO.a_E= a_E;


    //tune s quark mass on A ensemble

    Get_ms_A_ens(a_A);
    
    
    
    //determine strange and charm quark masses


    //convert data for ls correlators

    bool Get_ASCII_s= false;

    if(Get_ASCII_s) {
    //read binary files
    boost::filesystem::create_directory("../ls_correlator");

    vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});
    

    

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_l_l", "mix_l_s1", "mix_l_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../ls_correlator/"+channel);
	boost::filesystem::create_directory("../ls_correlator/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_P5P5", "OS_P5P5"});

          
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
	  boost::filesystem::create_directory("../ls_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../ls_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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



    bool Get_ASCII_s_dm= false;

    if(Get_ASCII_s_dm) {
    //read binary files
    boost::filesystem::create_directory("../ls_correlator");

    vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});
    

    

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_l1_s", "mix_l2_s"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../ls_correlator/"+channel);
	boost::filesystem::create_directory("../ls_correlator/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_P5P5", "OS_P5P5"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../tau_decay_strange_bin_mu_corr/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
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
	  boost::filesystem::create_directory("../ls_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../ls_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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




    //convert data for sc correlators

    bool Get_ASCII_c= true;

    if(Get_ASCII_c) {
    //read binary files
    boost::filesystem::create_directory("../sc_correlator");
    

    vector<string> Ens_T1({"C.06.80", "C.06.112",  "D.54.96", "E.44.112", "B.72.64"});
    vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cD211a.054.96", "cE211a.044.112", "cB211b.072.64"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_s_c1", "mix_s_c2", "mix_s_c3", "mix_c1_c1", "mix_c2_c2", "mix_c3_c3"});
      if(Ens_T1[it] == "C.06.80" || Ens_T1[it] =="B.72.64") channels = {"mix_s1_c1", "mix_s1_c2", "mix_s2_c1", "mix_s2_c2"};

      
      
      for(auto &channel : channels) {
	boost::filesystem::create_directory("../sc_correlator/"+channel);
	boost::filesystem::create_directory("../sc_correlator/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_P5P5", "TM_V0V0", "OS_V0V0", "OS_P5P5", "TM_VKVK", "OS_VKVK"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../charm_E_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
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
	  boost::filesystem::create_directory("../sc_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../sc_correlator/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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




    //determine strange quark masses on all the ensembles

    data_t ls_data_tm_P5P5, ls_H_data_tm_P5P5;
    data_t ls_data_OS_P5P5, ls_H_data_OS_P5P5;
   
    
    ls_data_tm_P5P5.Read("../ls_correlator/mix_l_s1", "mes_contr_mix_l_s1_TM_P5P5", "P5P5", Sort_easy);
    ls_H_data_tm_P5P5.Read("../ls_correlator/mix_l_s2", "mes_contr_mix_l_s2_TM_P5P5", "P5P5", Sort_easy);

    ls_data_OS_P5P5.Read("../ls_correlator/mix_l_s1", "mes_contr_mix_l_s1_OS_P5P5", "P5P5", Sort_easy);
    ls_H_data_OS_P5P5.Read("../ls_correlator/mix_l_s2", "mes_contr_mix_l_s2_OS_P5P5", "P5P5", Sort_easy);

    data_t l1s_data_tm_P5P5, l2s_data_tm_P5P5;
    
    l1s_data_tm_P5P5.Read("../ls_correlator/mix_l1_s", "mes_contr_mix_l1_s_TM_P5P5", "P5P5", Sort_easy);
    l2s_data_tm_P5P5.Read("../ls_correlator/mix_l2_s", "mes_contr_mix_l2_s_TM_P5P5", "P5P5", Sort_easy);

    data_t l1s_data_OS_P5P5, l2s_data_OS_P5P5;
    
    l1s_data_OS_P5P5.Read("../ls_correlator/mix_l1_s", "mes_contr_mix_l1_s_OS_P5P5", "P5P5", Sort_easy);
    l2s_data_OS_P5P5.Read("../ls_correlator/mix_l2_s", "mes_contr_mix_l2_s_OS_P5P5", "P5P5", Sort_easy);

    
    //pion masses in lattice units
    double amp_B_ave= 0.0565313; double amp_B_err= 1.438e-05;
    double amp_C_ave= 0.0472193; double amp_C_err= 3.45183e-05;
    double amp_D_ave= 0.0406214; double amp_D_err= 2.94047e-05;
    double amp_E_ave= 0.0338185; double amp_E_err= 3.18799e-05;
    
    distr_t aMp_B(UseJack), aMp_C(UseJack), aMp_D(UseJack), aMp_E(UseJack);
    
    for(int ijack=0;ijack<Njacks;ijack++ ) {
      aMp_B.distr.push_back( amp_B_ave + GM()*amp_B_err/sqrt(Njacks-1.0));
      aMp_C.distr.push_back( amp_C_ave + GM()*amp_C_err/sqrt(Njacks-1.0));
      aMp_D.distr.push_back( amp_D_ave + GM()*amp_D_err/sqrt(Njacks-1.0));
      aMp_E.distr.push_back( amp_E_ave + GM()*amp_E_err/sqrt(Njacks-1.0));
    }


    distr_t_list ams_list(UseJack), ams_corr_list(UseJack);
    distr_t_list ams_OS_list(UseJack);

    vector<string> Ens;


    int Nens = ls_data_tm_P5P5.size;

    for(int iens=0; iens<Nens;iens++) {


      Ens.push_back(ls_data_tm_P5P5.Tag[iens]);
      boost::filesystem::create_directory("../data/scale_setting/MK");
      boost::filesystem::create_directory("../data/scale_setting/MK/"+ls_data_tm_P5P5.Tag[iens]);
    
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = ls_data_tm_P5P5.nrows[iens];

      //get effective masses

      distr_t_list M_K = Corr.effective_mass_t( ls_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MK/"+ls_data_tm_P5P5.Tag[iens]+"/eff_mass_K");
      distr_t_list M_K_H = Corr.effective_mass_t( ls_H_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MK/"+ls_data_tm_P5P5.Tag[iens]+"/eff_mass_K_H");

      
      distr_t_list M_K_OS = Corr.effective_mass_t( ls_data_OS_P5P5.col(0)[iens], "../data/scale_setting/MK/"+ls_data_OS_P5P5.Tag[iens]+"/eff_mass_K_OS");
      distr_t_list M_K_H_OS = Corr.effective_mass_t( ls_H_data_OS_P5P5.col(0)[iens], "../data/scale_setting/MK/"+ls_data_OS_P5P5.Tag[iens]+"/eff_mass_K_H_OS");

      //distr_t_list dmK = Corr.effective_mass_t( summ_master( l1s_data_tm_P5P5.col(0)[iens], Multiply_Vvector_by_scalar(l2s_data_tm_P5P5.col(0)[iens], -1.0), ls_data_tm_P5P5.col(0)[iens]), "") - M_K;

      distr_t_list dmK = Corr.effective_mass_t(l1s_data_tm_P5P5.col(0)[iens],"") - Corr.effective_mass_t(l2s_data_tm_P5P5.col(0)[iens],"");
     
      distr_t_list dmK2 = POW_DL(Corr.effective_mass_t( l1s_data_tm_P5P5.col(0)[iens],""),2) - POW_DL(Corr.effective_mass_t(l2s_data_tm_P5P5.col(0)[iens],""),2);
      distr_t_list dmK2_OS = POW_DL(Corr.effective_mass_t( l1s_data_OS_P5P5.col(0)[iens],""),2) - POW_DL(Corr.effective_mass_t(l2s_data_OS_P5P5.col(0)[iens],""),2); 

      LatticeInfo L_info;
     
      L_info.LatInfo_new_ens(ls_data_tm_P5P5.Tag[iens]);
     
      double aml= L_info.ml;
      double ams1= L_info.ms_L_new;
      double ams2= L_info.ms_M_new;

      //get lattice spacing
     distr_t a_distr(UseJack);
     
     if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="B") {a_distr=a_B;  }
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="C") {a_distr=a_C; }
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="D") {a_distr=a_D; }
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="E") {a_distr=a_E; }
     else crash("lattice spacing distribution for Ens: "+ls_data_tm_P5P5.Tag[iens]+" not found");

  
    
     distr_t aMP;

     int Tmin_P5, Tmax_P5;

     if( ls_data_tm_P5P5.Tag[iens] =="cB211b.072.96")     { Tmin_P5=38; Tmax_P5=59; aMP=aMp_B;}
     else if(ls_data_tm_P5P5.Tag[iens] =="cB211b.072.64") { Tmin_P5=45; Tmax_P5=59; aMP=aMp_B;}
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="C")  { Tmin_P5=36; Tmax_P5=63; aMP=aMp_C;}
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="D")  { Tmin_P5=37; Tmax_P5=69; aMP=aMp_D;}
     else if(ls_data_tm_P5P5.Tag[iens].substr(1,1)=="E")  { Tmin_P5=60; Tmax_P5=85; aMP=aMp_E;}
     else crash("Cannot recognize the ensemble: "+ls_data_tm_P5P5.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");

     Corr.Tmin= Tmin_P5;
     Corr.Tmax= Tmax_P5;

     distr_t dmK2_fit = Corr.Fit_distr(dmK2)/(a_distr*a_distr);
     distr_t dmK2_fit_OS = Corr.Fit_distr(dmK2_OS)/(a_distr*a_distr);
    

     distr_t MK1= Corr.Fit_distr( M_K)/a_distr;
     distr_t MK2= Corr.Fit_distr( M_K_H )/a_distr;

     
     distr_t MK1_OS= Corr.Fit_distr( M_K_OS)/a_distr;
     distr_t MK2_OS= Corr.Fit_distr( M_K_H_OS )/a_distr;


     distr_t MK1_bis= Corr.Fit_distr( M_K)/a_distr;
     distr_t MK2_bis= Corr.Fit_distr( M_K_H)/a_distr;


     cout<<"lattice MK^2 mass for Ensemble: "<<ls_data_tm_P5P5.Tag[iens]<<endl;    
     cout<<"ams: "<<ams1<<" (aMK)^2: "<< (MK1*MK1*a_distr*a_distr).ave() <<"  "<<(MK1*MK1*a_distr*a_distr).err()<<endl;
     cout<<"ams: "<<ams2<<" (aMK)^2: "<< (MK2*MK2*a_distr*a_distr).ave()<<" " << (MK2*MK2*a_distr*a_distr).err()<<endl;
     cout<<"After correcting ml mistuning: "<<endl;
     cout<<"ams: "<<ams1<<" (aMK)^2: "<< (MK1_bis*MK1_bis*a_distr*a_distr).ave() <<"  "<<(MK1_bis*MK1_bis*a_distr*a_distr).err()<<endl;
     cout<<"ams: "<<ams2<<" (aMK)^2: "<< (MK2_bis*MK2_bis*a_distr*a_distr).ave()<<" " << (MK2_bis*MK2_bis*a_distr*a_distr).err()<<endl;
     cout<<"###############################"<<endl;
     
     //determine physical strange quark mass

     vector<distr_t> MMK2({MK1*MK1, MK2*MK2});
     vector<distr_t> MMK2_bis({MK1*MK1+ dmK2_fit, MK2*MK2 + dmK2_fit});
     vector<distr_t> MMK2_OS_bis({MK1_OS*MK1_OS+ dmK2_fit, MK2*MK2 + dmK2_fit});
     vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});

     distr_t MK_FLAG_corr= SQRT_D(  MK_FLAG*MK_FLAG + 0.5*( POW_D(aMP/a_distr,2)    - pow(MP_FLAG,2)));

     distr_t ams_phys = Obs_extrapolation_meson_mass(MMS, MMK2, MK_FLAG*MK_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MK"  , "ams_extrapolation_"+ls_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     distr_t ams_phys_corr = Obs_extrapolation_meson_mass(MMS, MMK2, MK_FLAG_corr*MK_FLAG_corr ,  "../data/scale_setting/MK"  , "ams_corr_extrapolation_"+ls_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     distr_t ams_phys_bis = Obs_extrapolation_meson_mass(MMS, MMK2_bis, MK_FLAG*MK_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MK"  , "ams_bis_extrapolation_"+ls_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     distr_t ams_phys_bis_OS = Obs_extrapolation_meson_mass(MMS, MMK2_OS_bis, MK_FLAG*MK_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MK"  , "ams_bis_OS_extrapolation_"+ls_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     
     ams_list.distr_list.push_back(ams_phys_bis);
     ams_OS_list.distr_list.push_back(ams_phys_bis_OS);
     ams_corr_list.distr_list.push_back(ams_phys_corr);

     SCALE_INFO.MK1.distr_list.push_back( SQRT_D(MK1*MK1+dmK2_fit));
     SCALE_INFO.MK2.distr_list.push_back( SQRT_D(MK2*MK2+dmK2_fit));
     
     cout<<"#### ENSEMBLE: "<<ls_data_tm_P5P5.Tag[iens]<<" ####"<<endl;
     cout<<"ams(lattice): "<<ams_phys_bis.ave()<<" +- "<<ams_phys_bis.err()<<endl;
     cout<<"ams(ChPT): "<<ams_phys_corr.ave()<<" +- "<<ams_phys_corr.err()<<endl;
     cout<<"ams(no-mistuning-corrections): "<<ams_phys.ave()<<" +- "<<ams_phys.err()<<endl;
     
    }


    //determine charm quark mass

    Get_mc_A_ens(a_A);

    //start from C80 and B64 ensemble [ 2 strange quark masses and 2 charm quark masses ]


    data_t s1c1_data_tm_P5P5, s1c2_data_tm_P5P5, s2c1_data_tm_P5P5, s2c2_data_tm_P5P5;

     
    s1c1_data_tm_P5P5.Read("../sc_correlator/mix_s1_c1", "mes_contr_mix_s1_c1_TM_P5P5", "P5P5", Sort_easy);
    s1c2_data_tm_P5P5.Read("../sc_correlator/mix_s1_c2", "mes_contr_mix_s1_c2_TM_P5P5", "P5P5", Sort_easy);
    s2c1_data_tm_P5P5.Read("../sc_correlator/mix_s2_c1", "mes_contr_mix_s2_c1_TM_P5P5", "P5P5", Sort_easy);
    s2c2_data_tm_P5P5.Read("../sc_correlator/mix_s2_c2", "mes_contr_mix_s2_c2_TM_P5P5", "P5P5", Sort_easy);
    
    int Nens_STRAT_1=s1c1_data_tm_P5P5.size;

    distr_t_list amc_phys_STRAT_1(UseJack);


    distr_t_list MDs1_list(UseJack);
    
    

    for(int iens=0; iens<Nens_STRAT_1; iens++) {

      boost::filesystem::create_directory("../data/scale_setting/MDs");
      boost::filesystem::create_directory("../data/scale_setting/MDs/"+s1c1_data_tm_P5P5.Tag[iens]);

      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = s1c1_data_tm_P5P5.nrows[iens];

      //get effective masses

      distr_t_list M_Ds_11 = Corr.effective_mass_t( s1c1_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+s1c1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_11");
      distr_t_list M_Ds_12 = Corr.effective_mass_t( s1c2_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+s1c1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_12");
      distr_t_list M_Ds_21 = Corr.effective_mass_t( s2c1_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+s1c1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_21");
      distr_t_list M_Ds_22 = Corr.effective_mass_t( s2c2_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+s1c1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_22");

      double ams1, ams2, amc1, amc2;
      string Ens=s1c1_data_tm_P5P5.Tag[iens];

      
      if(Ens == "cC211a.06.80") {
	ams1= 1.6000e-02;
	ams2= 1.7000e-02;
	amc1 = 1.8000e-01;
	amc2 = 1.9000e-01;
      }
      else if(Ens == "cB211b.072.64") {
	ams1= 1.7500e-02;
	ams2= 1.8500e-02;
	amc1 = 2.3000e-01;
	amc2 = 2.4000e-01;
      }
      else crash("Ensemble: "+Ens+" should not be analyzed with STRAT 1");

      int Tmin_P5, Tmax_P5;

     
      if(s1c1_data_tm_P5P5.Tag[iens].substr(1,1)=="C")  { Tmin_P5=32; Tmax_P5=43;}
      else if(s1c1_data_tm_P5P5.Tag[iens].substr(1,1)=="B")  { Tmin_P5=29; Tmax_P5=45;}
      else crash("Cannot recognize the ensemble: "+s1c1_data_tm_P5P5.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");
      
      Corr.Tmin= Tmin_P5;
      Corr.Tmax= Tmax_P5;



      distr_t a_distr = (Ens=="cC211a.06.80")?a_C:a_B;
      //get ams
      distr_t ams(UseJack);
      for(int bens=0; bens<(signed)ls_data_tm_P5P5.Tag.size(); bens++) {
	if(s1c1_data_tm_P5P5.Tag[iens] == ls_data_tm_P5P5.Tag[bens]) ams= ams_list[bens] ;
      }


      distr_t MDs_s1c1= Corr.Fit_distr( M_Ds_11)/a_distr;
      distr_t MDs_s1c2= Corr.Fit_distr( M_Ds_12 )/a_distr;
      distr_t MDs_s2c1= Corr.Fit_distr( M_Ds_21)/a_distr;
      distr_t MDs_s2c2= Corr.Fit_distr( M_Ds_22 )/a_distr;


      //determine physical strange quark mass

     vector<distr_t> MMDS1({MDs_s1c1, MDs_s2c1});
     vector<distr_t> MMDS2({MDs_s1c2, MDs_s2c2});
     vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});
     vector<distr_t> MMC({ Get_id_jack_distr(Njacks)*amc1, Get_id_jack_distr(Njacks)*amc2});

     SCALE_INFO.Ens_c.push_back(s1c1_data_tm_P5P5.Tag[iens]);
     MDs1_list.distr_list.push_back(Corr.Fit_distr(M_Ds_11));

     //interpolate MDs to the physical ms point
     distr_t MDs_1 = Obs_extrapolation_meson_mass(MMDS1, MMS, ams ,  "../data/scale_setting/MDs"  , "MDs_1_extrapolation_"+s1c1_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     distr_t MDs_2 = Obs_extrapolation_meson_mass(MMDS2, MMS, ams ,  "../data/scale_setting/MDs"  , "MDs_2_extrapolation_"+s1c1_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );

     //find amc

     vector<distr_t> MMDS({MDs_1, MDs_2});

     distr_t amc_phys=  Obs_extrapolation_meson_mass(MMC, MMDS, MDs_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MDs"  , "amc_extrapolation_"+s1c1_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );
     
      

     cout<<"#### ENSEMBLE: "<<s1c1_data_tm_P5P5.Tag[iens]<<" ####"<<endl;
     cout<<"amc: "<<amc_phys.ave()<<" +- "<<amc_phys.err()<<endl;
      
     amc_phys_STRAT_1.distr_list.push_back(amc_phys);

    }


    //##############################################################################################################################



    





























    //###############################################################################################################################

    


    //determine the charm quark mass on the other ensembles

    data_t sc1_data_tm_P5P5, sc2_data_tm_P5P5, sc3_data_tm_P5P5;

   

    data_t cc1_data_tm_P5P5, cc2_data_tm_P5P5, cc3_data_tm_P5P5;

    data_t cc1_data_OS_P5P5, cc2_data_OS_P5P5, cc3_data_OS_P5P5;

    data_t cc1_data_tm_V0V0, cc1_data_OS_V0V0;


    data_t cc1_data_tm_VKVK, cc1_data_OS_VKVK;

     
    sc1_data_tm_P5P5.Read("../sc_correlator/mix_s_c1", "mes_contr_mix_s_c1_TM_P5P5", "P5P5", Sort_easy);
    sc2_data_tm_P5P5.Read("../sc_correlator/mix_s_c2", "mes_contr_mix_s_c2_TM_P5P5", "P5P5", Sort_easy);
    sc3_data_tm_P5P5.Read("../sc_correlator/mix_s_c3", "mes_contr_mix_s_c3_TM_P5P5", "P5P5", Sort_easy);
    cc1_data_tm_P5P5.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_TM_P5P5", "P5P5", Sort_easy);
    cc2_data_tm_P5P5.Read("../sc_correlator/mix_c2_c2", "mes_contr_mix_c2_c2_TM_P5P5", "P5P5", Sort_easy);
    cc3_data_tm_P5P5.Read("../sc_correlator/mix_c3_c3", "mes_contr_mix_c3_c3_TM_P5P5", "P5P5", Sort_easy);
    cc1_data_tm_V0V0.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_TM_V0V0", "V0V0", Sort_easy);
    cc1_data_OS_V0V0.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_OS_V0V0", "V0V0", Sort_easy);

    cc1_data_OS_P5P5.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_OS_P5P5", "P5P5", Sort_easy);
    cc2_data_OS_P5P5.Read("../sc_correlator/mix_c2_c2", "mes_contr_mix_c2_c2_OS_P5P5", "P5P5", Sort_easy);
    cc3_data_OS_P5P5.Read("../sc_correlator/mix_c3_c3", "mes_contr_mix_c3_c3_OS_P5P5", "P5P5", Sort_easy);


    cc1_data_tm_VKVK.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_TM_VKVK", "VKVK", Sort_easy);
    cc1_data_OS_VKVK.Read("../sc_correlator/mix_c1_c1", "mes_contr_mix_c1_c1_OS_VKVK", "VKVK", Sort_easy);
    
    int Nens_c=sc1_data_tm_P5P5.size;

   

    distr_t_list amc_phys_list(UseJack);

    for(int iens=0; iens<Nens_c; iens++) {

      boost::filesystem::create_directory("../data/scale_setting/MDs");
      boost::filesystem::create_directory("../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]);

      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = sc1_data_tm_P5P5.nrows[iens];

      //get effective masses

      distr_t_list M_Ds_1 = Corr.effective_mass_t( sc1_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_1");
      distr_t_list M_Ds_2 = Corr.effective_mass_t( sc2_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_2");
      distr_t_list M_Ds_3 = Corr.effective_mass_t( sc3_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_Ds_3");


      distr_t_list M_etac_2 =  Corr.effective_mass_t( cc2_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_2");
      distr_t_list M_etac_1 =  Corr.effective_mass_t( cc1_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_1");
      distr_t_list M_etac_3 =  Corr.effective_mass_t( cc3_data_tm_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_3");

      distr_t_list M_etac_OS_2 =  Corr.effective_mass_t( cc2_data_OS_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_OS_2");
      distr_t_list M_etac_OS_1 =  Corr.effective_mass_t( cc1_data_OS_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_OS_1");
      distr_t_list M_etac_OS_3 =  Corr.effective_mass_t( cc3_data_OS_P5P5.col(0)[iens], "../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_etac_OS_3");

      distr_t_list M_V0_TM= Corr.effective_mass_t( cc1_data_tm_V0V0.col(0)[iens],"../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_V0_TM_1");
      distr_t_list M_V0_OS= Corr.effective_mass_t( cc1_data_OS_V0V0.col(0)[iens],"../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_V0_OS_1");


      distr_t_list M_Jpsi_TM_1 =  Corr.effective_mass_t( cc1_data_tm_VKVK.col(0)[iens],"../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_Jpsi_TM_1");
      distr_t_list M_Jpsi_OS_1 =  Corr.effective_mass_t( cc1_data_tm_VKVK.col(0)[iens],"../data/scale_setting/MDs/"+sc1_data_tm_P5P5.Tag[iens]+"/eff_mass_Jpsi_OS_1");

      double ams, amc1, amc2, amc3;
      distr_t a_distr(UseJack);

      if( sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="C")     { ams= 0.016067; amc1= 0.18; amc2=0.19; amc3=0.20; a_distr= a_C;}
      else if(sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="D") { ams= 1.3557e-02; amc1= 0.15; amc2=0.16; amc3=0.17; a_distr= a_D;} 
      else if(sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="E")  { ams= 1.1759e-02; amc1= 0.13; amc2=0.14; amc3=0.15; a_distr = a_E;} 
      else crash("Cannot recognize the ensemble: "+s1c1_data_tm_P5P5.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");


      int Tmin_P5, Tmax_P5;

      
      if(sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="C")  { Tmin_P5=40; Tmax_P5=48;}
      else if(sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="D")  { Tmin_P5=50; Tmax_P5=63;}
      else if(sc1_data_tm_P5P5.Tag[iens].substr(1,1)=="E")  { Tmin_P5=50; Tmax_P5=60;}
      else crash("Cannot recognize the ensemble: "+s1c1_data_tm_P5P5.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");
      
      Corr.Tmin= Tmin_P5;
      Corr.Tmax= Tmax_P5;



   
      //get ams
      distr_t ams_phys(UseJack);
      for(int bens=0; bens<(signed)ls_data_tm_P5P5.Tag.size(); bens++) {
	if(sc1_data_tm_P5P5.Tag[iens] == ls_data_tm_P5P5.Tag[bens]) ams_phys= ams_corr_list[bens] ;
      }
  
      cout<<"####ENSEMBLE: "<<sc1_data_tm_P5P5.Tag[iens]<<" ####"<<endl;
      cout<<"ams(phys): "<<ams_phys.ave()<<" +- "<<ams_phys.err()<<endl;
      cout<<"ams(simulated): "<<ams<<endl;


      distr_t MDs_c1= Corr.Fit_distr( M_Ds_1)/a_distr;
      distr_t MDs_c2= Corr.Fit_distr( M_Ds_2 )/a_distr;
      distr_t MDs_c3= Corr.Fit_distr( M_Ds_3)/a_distr;

      SCALE_INFO.Ens_c.push_back(sc1_data_tm_P5P5.Tag[iens]);
      MDs1_list.distr_list.push_back(Corr.Fit_distr(M_Ds_1));


      Corr.Tmin = 46;
      Corr.Tmax = 55;

      distr_t Metac1= Corr.Fit_distr(M_etac_1)/a_distr;
      distr_t Metac2= Corr.Fit_distr(M_etac_2)/a_distr;
      distr_t Metac3= Corr.Fit_distr(M_etac_3)/a_distr;

      vector<distr_t> MMETAC({Metac1, Metac2, Metac3});

      //determine physical strange quark mass

      vector<distr_t> MMDS({MDs_c1, MDs_c2, MDs_c3});
      
      vector<distr_t> MMC({ Get_id_jack_distr(Njacks)*amc1, Get_id_jack_distr(Njacks)*amc2, Get_id_jack_distr(Njacks)*amc3});

      


     
      
      distr_t amc_phys=  Obs_extrapolation_meson_mass(MMC, MMDS, MDs_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MDs"  , "amc_extrapolation_"+sc1_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );

      distr_t amc_phys_ETAC=  Obs_extrapolation_meson_mass(MMC, MMETAC, Metac_FLAG*Get_id_jack_distr(Njacks) ,  "../data/scale_setting/MDs"  , "amc_extrapolation_etac_"+sc1_data_tm_P5P5.Tag[iens]+".dat",  UseJack, "SPLINE" );

      distr_t etac_phys= Obs_extrapolation_meson_mass(MMETAC, MMC, amc_phys, "../data/scale_setting/MDs", "ametac_extrapolation_"+sc1_data_tm_P5P5.Tag[iens]+".dat", UseJack, "SPLINE");


      
     
      

      
      
      amc_phys_list.distr_list.push_back(amc_phys);

      cout<<"#### ENSEMBLE: "<<sc1_data_tm_P5P5.Tag[iens]<<" ####"<<endl;
      cout<<"amc: "<<amc_phys.ave()<<" +- "<<amc_phys.err()<<endl;
      cout<<"amc(from etac): "<<amc_phys_ETAC.ave()<<" +- "<<amc_phys_ETAC.err()<<endl;
      cout<<"am_etac: "<<etac_phys.ave()<<" +- "<<etac_phys.err()<<endl;

    }
    

    


    SCALE_INFO.ms = ams_list;
    SCALE_INFO.ms_OS = ams_OS_list;
    SCALE_INFO.Ens = Ens;
    SCALE_INFO.MDs1 = MDs1_list;
  
      
    //###################################################################################
      
      
    return SCALE_INFO ;

}


void Get_ms_A_ens(const distr_t &a_A) {

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



  //read data

  data_t l, s1, s2;


  l.Read("../A_ensembles_tuning/l", "mes_contr_l", "P5P5", Sort_light_confs);
  s1.Read("../A_ensembles_tuning/s1", "mes_contr_s1", "P5P5", Sort_light_confs);
  s2.Read("../A_ensembles_tuning/s2", "mes_contr_s2", "P5P5", Sort_light_confs); 


  boost::filesystem::create_directory("../data/scale_setting/MK");

  int Nens= l.size;

  cout<<"Number of ensembles: "<<Nens<<endl;

  for(int iens=0;iens<Nens;iens++) {


    

    string Ens= l.Tag[iens];

    double ams1, ams2;

   
    
    
    //set up statistical analysis
    CorrAnalysis Corr(UseJack,Njacks,100);
    Corr.Nt= l.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    cout<<"Analyzing ensemble: "<<Ens<<" T: "<<Corr.Nt<<endl;

    if(Ens=="cA211a.12.48") { ams1= 0.01760; ams2=0.02200; Corr.Tmin= 15; Corr.Tmax=30;}
    else if(Ens=="cA211a.30.32") { ams1=0.01760; ams2=0.02200; Corr.Tmin=15; Corr.Tmax=24;}
    else if(Ens=="cA211a.40.24") { ams1=0.01760; ams2=0.02200; Corr.Tmin=17; Corr.Tmax=22;}
    else crash("in Get_ms_A_ens() Ens: "+Ens+" not found");

    boost::filesystem::create_directory("../data/scale_setting/MK/"+Ens);

   
    distr_t_list Mpi_distr = Corr.effective_mass_t( l.col(0)[iens], "../data/scale_setting/MK/"+Ens+"/eff_mass_PI_unitary");
    distr_t_list MK_distr = Corr.effective_mass_t(s1.col(0)[iens], "../data/scale_setting/MK/"+Ens+"/eff_mass_K");
    distr_t_list MK_H_distr = Corr.effective_mass_t(s2.col(0)[iens], "../data/scale_setting/MK/"+Ens+"/eff_mass_K_H");


    distr_t Mpi = Corr.Fit_distr(Mpi_distr)/a_A;
    distr_t MK = Corr.Fit_distr(MK_distr)/a_A;
    distr_t MK_H = Corr.Fit_distr(MK_H_distr)/a_A;

    vector<distr_t> MMK2({ MK*MK , MK_H*MK_H });
    vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});


    distr_t MK_FLAG_corr= SQRT_D(  MK_FLAG*MK_FLAG + 0.5*( POW_D(Mpi,2)    - pow(MP_FLAG,2)));

    
    distr_t ams_phys_corr = Obs_extrapolation_meson_mass(MMS, MMK2, MK_FLAG_corr*MK_FLAG_corr ,  "../data/scale_setting/MK"  , "ams_corr_extrapolation_"+Ens+".dat",  UseJack, "SPLINE" );

    cout<<"Ens: "<<Ens<<" ams: "<<ams_phys_corr.ave()<<" +- "<<ams_phys_corr.err()<<endl;
    

  }



 
  return;

}



void Get_mc_A_ens(const distr_t &a_A) {

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



  //read data

  data_t c1, c2, c3;


  c1.Read("../A_ensembles_mc_tuning", "mes_contr_PT2_Ds1", "P5P5", Sort_light_confs);
  c2.Read("../A_ensembles_mc_tuning", "mes_contr_PT2_Ds2", "P5P5", Sort_light_confs);
  c3.Read("../A_ensembles_mc_tuning", "mes_contr_PT2_Ds3", "P5P5", Sort_light_confs); 


  boost::filesystem::create_directory("../data/scale_setting/MK");

  int Nens= 1;

  cout<<"Number of ensembles: "<<Nens<<endl;

  for(int iens=0;iens<Nens;iens++) {


    

    string Ens= c1.Tag[iens];

    double amc1= 0.255;
    double amc2= 0.265;
    double amc3= 0.275;

   
    
    
    //set up statistical analysis
    CorrAnalysis Corr(UseJack,Njacks,100);
    Corr.Nt= c1.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    cout<<"Analyzing ensemble: "<<Ens<<" T: "<<Corr.Nt<<endl;

    if(Ens=="cA211a.12.48") { Corr.Tmin= 15; Corr.Tmax=27;}
    else crash("in Get_mc_A_ens() Ens: "+Ens+" not found");

    boost::filesystem::create_directory("../data/scale_setting/MDs/"+Ens);

   
    distr_t_list MDs1_distr = Corr.effective_mass_t( c1.col(0)[iens], "../data/scale_setting/MDs/"+Ens+"/eff_mass_Ds_1.t");
    distr_t_list MDs2_distr = Corr.effective_mass_t( c2.col(0)[iens], "../data/scale_setting/MDs/"+Ens+"/eff_mass_Ds_2.t");
    distr_t_list MDs3_distr = Corr.effective_mass_t( c3.col(0)[iens], "../data/scale_setting/MDs/"+Ens+"/eff_mass_Ds_3.t");


   
    distr_t MDs1 = Corr.Fit_distr(MDs1_distr)/a_A;
    distr_t MDs2 = Corr.Fit_distr(MDs2_distr)/a_A;
    distr_t MDs3 = Corr.Fit_distr(MDs3_distr)/a_A;

    vector<distr_t> MMDs({ MDs1 , MDs2,MDs3 });
    vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*amc1, Get_id_jack_distr(Njacks)*amc2, Get_id_jack_distr(Njacks)*amc3});

    distr_t MDs_FLAG_corr = Get_id_jack_distr(Njacks)*MDs_FLAG;
    
        
    distr_t amc_phys = Obs_extrapolation_meson_mass(MMS, MMDs, MDs_FLAG_corr ,  "../data/scale_setting/MDs"  , "amc_corr_extrapolation_"+Ens+".dat",  UseJack, "SPLINE" );

    cout<<"Ens: "<<Ens<<" amc: "<<amc_phys.ave()<<" +- "<<amc_phys.err()<<endl;
    

  }


  exit(-1);
 
  return;

}
