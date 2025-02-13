#include "../include/weak_annihilation_inclusive.h"

using namespace std;


void compute_weak_annihilation() {

  int Njacks=50;
  bool UseJack=1;
  double fm_to_inv_Gev= 1.0/0.197327;
  int Nhits=1;
  int tsink=48;

  boost::filesystem::create_directory("../data/WA_inclusive_test");

  auto Sort_confs = [](string A, string B) {

			   

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
  

  vector<vector<data_t>> V_sou(Nhits), A_sou(Nhits);
  vector<vector<data_t>> V_sink(Nhits), A_sink(Nhits);
  vector<vector<data_t>> V_sou_opp(Nhits);
  vector<vector<data_t>> V_sink_opp(Nhits);
  
  for(int ihit=0;ihit<Nhits;ihit++) {
    V_sou[ihit].resize(4);
    A_sou[ihit].resize(4);
    V_sink[ihit].resize(4);
    A_sink[ihit].resize(4);
    V_sou_opp[ihit].resize(4);
    V_sink_opp[ihit].resize(4);
  }


  vector<data_t> P5_sou(Nhits), P5_sink(Nhits);
  
  for(int ihit=0;ihit<Nhits;ihit++) {

    cout<<"ihit: "<<ihit<<endl;

  
    for(int i=0;i<4; i++) {

      
      V_sou[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_H_S_PH_S0_H_th2Zp", "V"+to_string(i)+"P5", Sort_confs);
      A_sou[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_H_S_PH_S0_H_th2Zp", "A"+to_string(i)+"P5", Sort_confs);
      V_sou_opp[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_H_S_PH_S0_H_th2Zm", "V"+to_string(i)+"P5", Sort_confs);

        
      V_sink[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_PH_S48_H_th2Zp_H_S", "V"+to_string(i)+"P5", Sort_confs);
      A_sink[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_PH_S48_H_th2Zp_H_S", "A"+to_string(i)+"P5", Sort_confs);
      V_sink_opp[ihit][i].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_PH_S48_H_th2Zm_H_S", "V"+to_string(i)+"P5", Sort_confs);
        
      P5_sou[ihit].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_H_S_H_H_S0_H" , "P5P5", Sort_confs);
      P5_sink[ihit].Read("../weak_annihilation_test/out_disconnected_copy"+to_string(ihit), "mes_contr_H_S_H_H_S0_H", "P5P5", Sort_confs);

    }
  }

   


    int Nens=V_sou[0][0].Tag.size();



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






    boost::filesystem::create_directory("../data/WA_inclusive_test");

    

    for(int iens=0;iens<Nens;iens++) {

      cout<<"Analyzing ensemble: "<<P5_sou[0].Tag[iens]<<endl; 

      CorrAnalysis Corr(UseJack, Njacks,800);
      Corr.Nt = P5_sou[0].nrows[iens];
      int T=Corr.Nt;
      double V = (double)pow(T/2,3.0);
      double V4 = V*((double)T);

      Corr.Perform_Nt_t_average =0;
     
      distr_t a_distr(UseJack);
      distr_t Zv(UseJack), Za(UseJack);
      
      if(P5_sou[0].Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; }
      else if(P5_sou[0].Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
      else if(P5_sou[0].Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
      else crash("Ensemble not found");

      LatticeInfo L_info;
      
      L_info.LatInfo_new_ens(P5_sou[0].Tag[iens]);
      

      int Nconfs=V_sou[0][1].col(0)[iens][0].size();
      
      complex_distr_t_list V1V1(0, T*T, Nconfs);
      complex_distr_t_list V1V1_opp(0, T*T, Nconfs);
      complex_distr_t_list V1V1_a1(0, T*T, Nconfs);
      complex_distr_t_list V1V1_s1(0, T*T, Nconfs);
      complex_distr_t_list V1V1_a0(0, T*T, Nconfs);
      complex_distr_t_list V1V1_s0(0, T*T, Nconfs);

      Vfloat t1_s, t2_s;

      cout<<"V1V1(size after init.): "<<V1V1.size()<<endl;
    
      for(int ihit=0;ihit<Nhits;ihit++) {

	cout<<"Nconfs(ihit: "<<ihit<<"): "<<V_sou[ihit][1].col(0)[iens][0].size()<<endl;


	complex_distr_t_list A(0, V_sou[ihit][1].col(0)[iens], V_sou[ihit][1].col(1)[iens]);
	complex_distr_t_list A_opp(0, V_sou_opp[ihit][1].col(0)[iens], V_sou_opp[ihit][1].col(1)[iens]);
	complex_distr_t_list B(0, V_sink[ihit][1].col(0)[iens], V_sink[ihit][1].col(1)[iens]);
	complex_distr_t_list B_opp(0, V_sink_opp[ihit][1].col(0)[iens], V_sink_opp[ihit][1].col(1)[iens]);
	complex_distr_t_list B_dag=  B.dagger();
	complex_distr_t_list B_opp_dag = B_opp.dagger();

	int count=0;
	
	for(int t1=0;t1<Corr.Nt;t1++) {
	  for(int t2=0;t2<Corr.Nt;t2++) {
	    
	    if(ihit==0) { t1_s.push_back(t1); t2_s.push_back(t2); }
	    
	    complex_distr_t TT= pow(64,3.0)*A.distr_list[t1]*B_dag.distr_list[t2];
	    complex_distr_t TT_opp= pow(64,3.0)*A_opp.distr_list[t1]*B_opp_dag.distr_list[t2];
	    complex_distr_t TT_a1= TT + pow(64,3.0)*A.distr_list[t1]*B_opp_dag.distr_list[t2];
	    complex_distr_t TT_s1= TT - pow(64,3.0)*A.distr_list[t1]*B_opp_dag.distr_list[t2];
	    complex_distr_t TT_a0=  TT + pow(64,3.0)*A_opp.distr_list[t1]*B_dag.distr_list[t2];
	    complex_distr_t TT_s0=  TT - pow(64,3.0)*A_opp.distr_list[t1]*B_dag.distr_list[t2];

	    if(ihit==0) {
	      V1V1.distr_list[ count] = TT/Nhits;
	      V1V1_opp.distr_list[ count] = TT_opp/Nhits;

	      V1V1_a1.distr_list[ count] = TT_a1/Nhits;
	      V1V1_s1.distr_list[ count] = TT_s1/Nhits;
	      V1V1_a0.distr_list[ count] = TT_a0/Nhits;
	      V1V1_s0.distr_list[ count] = TT_s0/Nhits;
	    }
	    else {
	      V1V1.distr_list[ count] = V1V1.distr_list[count] + TT/Nhits;
	      V1V1_opp.distr_list[ count] = V1V1_opp.distr_list[count] + TT_opp/Nhits;
	      V1V1_a1.distr_list[ count] = V1V1_a1.distr_list[count] + TT_a1/Nhits;
	      V1V1_s1.distr_list[ count] = V1V1_s1.distr_list[count] + TT_s1/Nhits;
	      V1V1_a0.distr_list[ count] = V1V1_a0.distr_list[count] + TT_a0/Nhits;
	      V1V1_s0.distr_list[ count] = V1V1_s0.distr_list[count] + TT_s0/Nhits;
	    }

	    count++;
	  }

	}

      }


      //jackknife analysis of V1V1

      Corr.Nt= T*T;

      complex_distr_t_list V1V1_jack( Corr.corr_t( V1V1.Get_vvector(0) , ""), Corr.corr_t(V1V1.Get_vvector(1), ""));
      complex_distr_t_list V1V1_opp_jack( Corr.corr_t( V1V1_opp.Get_vvector(0) , ""), Corr.corr_t(V1V1_opp.Get_vvector(1), ""));

      complex_distr_t_list V1V1_a1_jack( Corr.corr_t( V1V1_a1.Get_vvector(0) , ""), Corr.corr_t(V1V1_a1.Get_vvector(1), ""));
      complex_distr_t_list V1V1_s1_jack( Corr.corr_t( V1V1_s1.Get_vvector(0) , ""), Corr.corr_t(V1V1_s1.Get_vvector(1), ""));
      complex_distr_t_list V1V1_a0_jack( Corr.corr_t( V1V1_a0.Get_vvector(0) , ""), Corr.corr_t(V1V1_a0.Get_vvector(1), ""));
      complex_distr_t_list V1V1_s0_jack( Corr.corr_t( V1V1_s0.Get_vvector(0) , ""), Corr.corr_t(V1V1_s0.Get_vvector(1), ""));

      complex_distr_t_list V1V1_a1a0_jack = V1V1_a1_jack + V1V1_s0_jack;
      complex_distr_t_list V1V1_s1s0_jack = V1V1_s1_jack + V1V1_a1_jack;

      complex_distr_t_list V1V1_ave_jack = 0.5*(V1V1_jack+ V1V1_opp_jack);

      //print

      Print_To_File( {} , { t1_s, t2_s, V1V1_jack.RE_ave(), V1V1_jack.RE_err(), V1V1_jack.IM_ave(), V1V1_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_tsink_"+to_string(tsink), "" , "");

      Print_To_File( {} , { t1_s, t2_s, V1V1_opp_jack.RE_ave(), V1V1_opp_jack.RE_err(), V1V1_opp_jack.IM_ave(), V1V1_opp_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_opp_tsink_"+to_string(tsink), "" , "");

      Print_To_File( {} , { t1_s, t2_s, V1V1_ave_jack.RE_ave(), V1V1_ave_jack.RE_err(), V1V1_ave_jack.IM_ave(), V1V1_ave_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_ave_tsink_"+to_string(tsink), "" , "");
      
      Print_To_File( {} , { t1_s, t2_s, V1V1_a1_jack.RE_ave(), V1V1_a1_jack.RE_err(), V1V1_a1_jack.IM_ave(), V1V1_a1_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_a1_tsink_"+to_string(tsink), "" , "");
      
      Print_To_File( {} , { t1_s, t2_s, V1V1_s1_jack.RE_ave(), V1V1_s1_jack.RE_err(), V1V1_s1_jack.IM_ave(), V1V1_s1_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_s1_tsink_"+to_string(tsink), "" , "");

      Print_To_File( {} , { t1_s, t2_s, V1V1_a0_jack.RE_ave(), V1V1_a0_jack.RE_err(), V1V1_a0_jack.IM_ave(), V1V1_a0_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_a0_tsink_"+to_string(tsink), "" , "");

      Print_To_File( {} , { t1_s, t2_s, V1V1_s0_jack.RE_ave(), V1V1_s0_jack.RE_err(), V1V1_s0_jack.IM_ave(), V1V1_s0_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_s0_tsink_"+to_string(tsink), "" , "");
      
      Print_To_File( {} , { t1_s, t2_s, V1V1_a1a0_jack.RE_ave(), V1V1_a1a0_jack.RE_err(), V1V1_a1a0_jack.IM_ave(), V1V1_a1a0_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_a1a0_tsink_"+to_string(tsink), "" , "");
      
      Print_To_File( {} , { t1_s, t2_s, V1V1_s1s0_jack.RE_ave(), V1V1_s1s0_jack.RE_err(), V1V1_s1s0_jack.IM_ave(), V1V1_s1s0_jack.IM_err() }, "../data/WA_inclusive_test/V1V1_s1s0_tsink_"+to_string(tsink), "" , "");
    }


    return;
    
    

} 
  
