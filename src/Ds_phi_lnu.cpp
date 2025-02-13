#include "../include/Ds_phi_lnu.h"
#include "numerics.h"
using namespace std;


const int Nboots = 800;
const bool UseJack=true;
const int Njacks=20;
const double Lambda_QCD= 0.3; //300 MeV
const double sign_kz = -1;
const double nqs = 5;
const int Nsou=3;



void Compute_Ds_phi_lnu() {


  
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

   
  vector<vector<data_t>> PT3_V1_data(nqs);
  vector<vector<data_t>> PT3_V2_data(nqs);
  vector<vector<data_t>> PT3_V3_data(nqs);
  vector<vector<data_t>> PT3_A0_data(nqs);
  vector<vector<data_t>> PT3_A1_data(nqs);
  vector<vector<data_t>> PT3_A2_data(nqs);
  vector<vector<data_t>> PT3_A3_data(nqs);
  vector<vector<data_t>> PT3_P5_data(nqs);

  vector<vector<data_t>> PT2_ss_data(nqs);
  vector<vector<data_t>> PT2_cs_data(nqs);
   
 
  for(int i=0;i<nqs;i++) {
    PT3_V1_data[i].resize(Nsou);
    PT3_V2_data[i].resize(Nsou);
    PT3_V3_data[i].resize(Nsou);
    PT3_A0_data[i].resize(Nsou);
    PT3_P5_data[i].resize(Nsou);
    PT3_A1_data[i].resize(Nsou);
    PT3_A2_data[i].resize(Nsou);
    PT3_A3_data[i].resize(Nsou);
    PT2_ss_data[i].resize(Nsou);
    PT2_cs_data[i].resize(Nsou);
  }
 

  for(int i=0;i<nqs;i++) {

    string read_path="";
    if(i==0) read_path="trivial";
    else if(i==1) read_path="traj_b";
    else if(i==2) read_path="traj_a";
    else if(i==3) read_path="rest";
    else if(i==4) read_path="rest";
    else crash("iqs="+to_string(i)+" not recognized");
    
    for(int j=0;j<Nsou;j++) {

      if(i!=3) {
      
	PT3_V1_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "V1P5", Sort_light_confs);
	PT3_V2_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "V2P5", Sort_light_confs);
	PT3_V3_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "V3P5", Sort_light_confs);
	PT3_A0_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "A0P5", Sort_light_confs);
	PT3_P5_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "P5P5", Sort_light_confs);
	PT3_A1_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "A1P5", Sort_light_confs);
	PT3_A2_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "A2P5", Sort_light_confs);
	PT3_A3_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_SOU_"+to_string(j), "A3P5", Sort_light_confs);
	PT2_ss_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT2_SS_0_SOU_"+to_string(j), "V1V1", Sort_light_confs);
	PT2_cs_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT2_BS_0_SOU_"+to_string(j), "P5P5", Sort_light_confs);

      }

      else {

	PT3_V1_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "V1P5", Sort_light_confs);
	PT3_V2_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "V2P5", Sort_light_confs);
	PT3_V3_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "V3P5", Sort_light_confs);
	PT3_A0_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "A0P5", Sort_light_confs);
	PT3_P5_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "P5P5", Sort_light_confs);
	PT3_A1_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "A1P5", Sort_light_confs);
	PT3_A2_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "A2P5", Sort_light_confs);
	PT3_A3_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT3_MB0_eps_SOU_"+to_string(j), "A3P5", Sort_light_confs);
	PT2_ss_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT2_SS_0_SOU_"+to_string(j), "V1V1", Sort_light_confs);
	PT2_cs_data[i][j].Read("../Bs_phi_gamma_data/"+read_path, "mes_contr_PT2_BS_0_SOU_"+to_string(j), "P5P5", Sort_light_confs);


      }
      
	
    }
    
  }

  
  int Nens= PT2_ss_data[0][0].Tag.size();

  for(int iens=0; iens<Nens;iens++) {

    distr_t_list FF_V_list(UseJack), FF_A0_list(UseJack), FF_A1_list(UseJack), FF_A2_list(UseJack);
    distr_t_list MDs_list(UseJack);

    Vfloat q2_list;

    //RCs
    distr_t Za, Zv, Z_T, a_distr;
    if(PT2_ss_data[0][0].Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A; Z_T=ZT_A; a_distr=a_A;}
    else if(PT2_ss_data[0][0].Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; Z_T=ZT_B; a_distr=a_B;}
    else if(PT2_ss_data[0][0].Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; Z_T=ZT_C; a_distr=a_C;}
    else if(PT2_ss_data[0][0].Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; Z_T=ZT_D; a_distr=a_D;}
    else crash("Ensemble: "+PT2_ss_data[0][0].Tag[iens]+" not recognised");


    Vfloat Thetas;
    vector<int> Tins;
    double mass_s, mass_c;
    if(PT2_ss_data[0][0].Tag[iens]=="cB211b.072.64") {
      Thetas = { 1.734914, 1.05335, 0.640817, -5e-5, 0.000000};
      mass_s=0.0184;
      mass_c= 0.237;
      Tins = {32, 32, 32, 25, 25};
    }
    else crash("Ensemble: "+PT2_ss_data[0][0].Tag[iens]+" not found");

  
    

    boost::filesystem::create_directory("../data/Ds_phi_lnu");
    
    cout<<"Analyzing ensemble: "<<PT2_ss_data[0][0].Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(PT2_ss_data[0][0].Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = PT2_ss_data[0][0].nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    distr_t M_phi_rest;

   

    for(int iq=4;iq>=0;iq--) {
      boost::filesystem::create_directory("../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]);
      boost::filesystem::create_directory("../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr");
      boost::filesystem::create_directory("../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/masses");
      boost::filesystem::create_directory("../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/FF");

           
      Vfloat TT, TT_ins;
      for(int t=0;t<Corr.Nt;t++) {TT_ins.push_back( Tins[iq]-t); };


      VVfloat PT3_Vx, PT3_Vyz, PT3_A0, PT3_Ax, PT3_Ayz, PT3_A0_IM, PT3_Ax_IM, PT3_Ayz_IM, PT3_P5;
      VVfloat PT2_SS, PT2_CS;
      
      for(int is=0;is<Nsou;is++) {

	if(is==0) {

	  if(iq != 3) {
	    PT3_Vx = PT3_V1_data[iq][is].col(1)[iens];
	    PT3_Vyz = summ_master( PT3_V2_data[iq][is].col(1)[iens]  ,  Multiply_Vvector_by_scalar( PT3_V3_data[iq][is].col(1)[iens], -1.0));
	    PT3_A0 = PT3_A0_data[iq][is].col(0)[iens];
	    PT3_A0_IM = PT3_A0_data[iq][is].col(1)[iens];
	    PT3_P5 = PT3_P5_data[iq][is].col(1)[iens];
	    PT3_Ax = PT3_A1_data[iq][is].col(0)[iens];
	    PT3_Ax_IM = PT3_A1_data[iq][is].col(1)[iens];
	    PT3_Ayz = summ_master( PT3_A2_data[iq][is].col(0)[iens], PT3_A3_data[iq][is].col(0)[iens]);
	    PT3_Ayz_IM = summ_master( PT3_A2_data[iq][is].col(1)[iens], PT3_A3_data[iq][is].col(1)[iens]);
	    PT2_SS= PT2_ss_data[iq][is].col(0)[iens];
	    PT2_CS= PT2_cs_data[iq][is].col(0)[iens];
	  }
	  else { //iq == 3, evaluate derivative

	    PT3_Vx = summ_master( PT3_V1_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_V1_data[4][is].col(1)[iens], -1.0));
	    PT3_Vyz = summ_master( PT3_V2_data[iq][is].col(1)[iens]  ,  Multiply_Vvector_by_scalar( PT3_V3_data[iq][is].col(1)[iens], -1.0) , Multiply_Vvector_by_scalar(PT3_V2_data[4][is].col(1)[iens], -1.0),  PT3_V3_data[4][is].col(1)[iens]);
	    PT3_A0 = summ_master( PT3_A0_data[iq][is].col(0)[iens], Multiply_Vvector_by_scalar(PT3_A0_data[4][is].col(0)[iens],-1.0));
	    PT3_A0_IM = summ_master( PT3_A0_data[iq][is].col(1)[iens] , Multiply_Vvector_by_scalar(PT3_A0_data[4][is].col(1)[iens], -1.0));
	    PT3_P5 = summ_master( PT3_P5_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_P5_data[4][is].col(1)[iens], -1.0));
	    PT3_Ax = PT3_A1_data[iq][is].col(0)[iens];
	    PT3_Ax_IM = PT3_A1_data[iq][is].col(1)[iens];
	    PT3_Ayz =  summ_master( PT3_A2_data[iq][is].col(0)[iens], PT3_A3_data[iq][is].col(0)[iens], Multiply_Vvector_by_scalar(PT3_A2_data[4][is].col(0)[iens], -1.0), Multiply_Vvector_by_scalar(PT3_A3_data[4][is].col(0)[iens], -1.0) );
	    PT3_Ayz_IM = summ_master( PT3_A2_data[iq][is].col(1)[iens], PT3_A3_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_A2_data[4][is].col(1)[iens], -1.0), Multiply_Vvector_by_scalar(PT3_A3_data[4][is].col(1)[iens], -1.0) );
	    PT2_SS= PT2_ss_data[iq][is].col(0)[iens];
	    PT2_CS= PT2_cs_data[iq][is].col(0)[iens];
	   
	  }
	    
	}
	else {

	  if(iq != 3) {
	  
	    PT3_Vx = summ_master( PT3_Vx, PT3_V1_data[iq][is].col(1)[iens]);
	    PT3_Vyz = summ_master( PT3_Vyz, summ_master( PT3_V2_data[iq][is].col(1)[iens]  ,  Multiply_Vvector_by_scalar( PT3_V3_data[iq][is].col(1)[iens], -1.0)));
	    PT3_A0 = summ_master( PT3_A0, PT3_A0_data[iq][is].col(0)[iens]);
	    PT3_A0_IM = summ_master( PT3_A0_IM, PT3_A0_data[iq][is].col(1)[iens]);
	    PT3_P5 = summ_master( PT3_P5, PT3_P5_data[iq][is].col(1)[iens]);
	    PT3_Ax = summ_master( PT3_Ax, PT3_A1_data[iq][is].col(0)[iens]);
	    PT3_Ax_IM = summ_master( PT3_Ax_IM, PT3_A1_data[iq][is].col(1)[iens]);
	    PT3_Ayz = summ_master( PT3_Ayz, summ_master( PT3_A2_data[iq][is].col(0)[iens], PT3_A3_data[iq][is].col(0)[iens]));
	    PT3_Ayz_IM = summ_master( PT3_Ayz_IM, summ_master( PT3_A2_data[iq][is].col(1)[iens], PT3_A3_data[iq][is].col(1)[iens]));
	    PT2_SS= summ_master( PT2_SS, PT2_ss_data[iq][is].col(0)[iens]);
	    PT2_CS= summ_master( PT2_CS, PT2_cs_data[iq][is].col(0)[iens]);

	  }
	  else { //iq==3, evaluate derivative
	    
	    PT3_Vx = summ_master( PT3_Vx, PT3_V1_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_V1_data[4][is].col(1)[iens], -1.0));
	    PT3_Vyz = summ_master( PT3_Vyz, summ_master( PT3_V2_data[iq][is].col(1)[iens]  ,  Multiply_Vvector_by_scalar( PT3_V3_data[iq][is].col(1)[iens], -1.0)),   summ_master( PT3_V3_data[4][is].col(1)[iens]  ,  Multiply_Vvector_by_scalar( PT3_V2_data[4][is].col(1)[iens], -1.0))  ) ;
	    PT3_A0 = summ_master( PT3_A0, PT3_A0_data[iq][is].col(0)[iens], Multiply_Vvector_by_scalar(PT3_A0_data[4][is].col(0)[iens], -1.0));
	    PT3_A0_IM = summ_master( PT3_A0_IM, PT3_A0_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_A0_data[4][is].col(1)[iens], -1.0));
	    PT3_P5 = summ_master( PT3_P5, PT3_P5_data[iq][is].col(1)[iens], Multiply_Vvector_by_scalar(PT3_P5_data[4][is].col(1)[iens], -1.0));
	    PT3_Ax = summ_master( PT3_Ax, PT3_A1_data[iq][is].col(0)[iens]);
	    PT3_Ax_IM = summ_master( PT3_Ax_IM, PT3_A1_data[iq][is].col(1)[iens]);
	    PT3_Ayz = summ_master( PT3_Ayz, summ_master( PT3_A2_data[iq][is].col(0)[iens], PT3_A3_data[iq][is].col(0)[iens]),  Multiply_Vvector_by_scalar( summ_master( PT3_A2_data[4][is].col(0)[iens], PT3_A3_data[4][is].col(0)[iens]) , -1.0)   );
	    PT3_Ayz_IM = summ_master( PT3_Ayz_IM, summ_master( PT3_A2_data[iq][is].col(1)[iens], PT3_A3_data[iq][is].col(1)[iens]), Multiply_Vvector_by_scalar( summ_master( PT3_A2_data[4][is].col(1)[iens], PT3_A3_data[4][is].col(1)[iens])   ,-1.0));
	    PT2_SS= summ_master( PT2_SS, PT2_ss_data[iq][is].col(0)[iens]);
	    PT2_CS= summ_master( PT2_CS, PT2_cs_data[iq][is].col(0)[iens]);
	    


	  }
	  
	}
	
      }


      PT3_Vx =  Multiply_Vvector_by_scalar(PT3_Vx, 1.0/Nsou);
      PT3_Vyz =  Multiply_Vvector_by_scalar(PT3_Vyz, 0.5/Nsou);
      PT3_A0 =  Multiply_Vvector_by_scalar(PT3_A0, 1.0/Nsou);
      PT3_A0_IM =  Multiply_Vvector_by_scalar(PT3_A0_IM, 1.0/Nsou);
      PT3_P5 =  Multiply_Vvector_by_scalar(PT3_P5, 1.0/Nsou);
      PT3_Ax =  Multiply_Vvector_by_scalar(PT3_Ax, 1.0/Nsou);
      PT3_Ax_IM =  Multiply_Vvector_by_scalar(PT3_Ax_IM, 1.0/Nsou);
      PT3_Ayz =  Multiply_Vvector_by_scalar(PT3_Ayz, 0.5/Nsou);
      PT3_Ayz_IM =  Multiply_Vvector_by_scalar(PT3_Ayz_IM, 0.5/Nsou);
      PT2_SS =  Multiply_Vvector_by_scalar(PT2_SS, 1.0/Nsou);
      PT2_CS =  Multiply_Vvector_by_scalar(PT2_CS, 1.0/Nsou);
      
      //analyze correlators
      Corr.Perform_Nt_t_average=0;
      distr_t_list PT3_Vx_distr = Corr.corr_t(PT3_Vx, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Vx_iq_"+to_string(iq));
      distr_t_list PT3_Vyz_distr = Corr.corr_t(PT3_Vyz, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Vyz_iq_"+to_string(iq));
      distr_t_list PT3_A0_distr = Corr.corr_t(PT3_A0, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_A0_iq_"+to_string(iq));
      distr_t_list PT3_A0_IM_distr = Corr.corr_t(PT3_A0_IM, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_A0_IM_iq_"+to_string(iq));
      distr_t_list PT3_P5_distr = Corr.corr_t(PT3_P5, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_P5_iq_"+to_string(iq));
      distr_t_list PT3_Ax_distr = Corr.corr_t(PT3_Ax, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Ax_iq_"+to_string(iq));
      distr_t_list PT3_Ax_IM_distr = Corr.corr_t(PT3_Ax_IM, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Ax_IM_iq_"+to_string(iq));
      distr_t_list PT3_Ayz_distr = Corr.corr_t(PT3_Ayz, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Ayz_iq_"+to_string(iq));
      distr_t_list PT3_Ayz_IM_distr = Corr.corr_t(PT3_Ayz_IM, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT3_Ayz_IM_iq_"+to_string(iq));\
      Corr.Perform_Nt_t_average=1;
      distr_t_list PT2_SS_distr = Corr.corr_t(PT2_SS, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT2_SS_iq_"+to_string(iq));
      distr_t_list PT2_CS_distr = Corr.corr_t(PT2_CS, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/corr/PT2_CS_iq_"+to_string(iq));
      distr_t_list Ds_eff_mass= Corr.effective_mass_t( PT2_CS_distr, "../data/Ds_phi_lnu/"+PT2_cs_data[0][0].Tag[iens]+"/masses/Ds_mass_iq_"+to_string(iq));
      distr_t_list phi_eff_mass= Corr.effective_mass_t( PT2_SS_distr, "../data/Ds_phi_lnu/"+PT2_ss_data[0][0].Tag[iens]+"/masses/phi_mass_iq_"+to_string(iq));



      //######
      

      int Tmin_Ds=0, Tmax_Ds=0;
      int Tmin_phi=0, Tmax_phi=0;
          
      if(PT2_ss_data[0][0].Tag[iens] =="cB211b.072.64") {
	if(iq==4 || iq==3) {Tmin_Ds=16; Tmax_Ds=30;}
	else {Tmin_Ds=20; Tmax_Ds=32; }
	Tmin_phi=15; Tmax_phi=20;
      }
      else crash("Ensemble: "+PT2_ss_data[0][0].Tag[iens]+" not found");

     
      Corr.Tmin=Tmin_Ds; Corr.Tmax=Tmax_Ds;
     
      distr_t M_Ds = Corr.Fit_distr(Ds_eff_mass);
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t M_phi = Corr.Fit_distr(phi_eff_mass);
      if(iq==3 || iq == 4) M_phi_rest=M_phi;
     
      MDs_list.distr_list.push_back( M_Ds/a_distr);
      
      //get momentum
      pt3_momenta pt3_mom(Vfloat({0.0,0.0,0.0}), Vfloat({0.0, 0.0, 0.0}), Vfloat({Thetas[iq], Thetas[iq], Thetas[iq]}), mass_c, mass_s, 0.0, L_info.L, L_info.T);
      distr_t q0= M_Ds - M_phi;

  

      double kz = pt3_mom.k()[2];
      //double q2 = ( (q0*q0 -3*kz*kz)/(M_Ds*M_Ds) ).ave();
      double q2 = ((M_Ds*M_Ds + M_phi_rest*M_phi_rest - 2*M_Ds*M_phi)/(M_Ds*M_Ds)).ave();
      if(iq==0) q2=0.0;
      distr_t aq2_distr=  (M_Ds*M_Ds + M_phi_rest*M_phi_rest - 2*M_Ds*M_phi); 

      cout<<"MDs: "<<(M_Ds/a_distr).ave()<<" +- "<<(M_Ds/a_distr).err()<<endl;
      cout<<"Mphi: "<<(M_phi_rest/a_distr).ave()<<" +- "<<(M_phi_rest/a_distr).err()<<endl;
      cout<<"akz: "<<kz<<endl;
      cout<<"aq0: "<<q0.ave()<<endl;
      cout<<"(a^2q^2): "<<aq2_distr.ave()<<endl;
      cout<<"iq "<<iq<<" q2/M2: "<<q2<<endl;

    
      
      Corr.Tmin=Tmin_Ds; Corr.Tmax=Tmax_Ds;
      distr_t F_D= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_CS_distr, ""))/2.0;
      Corr.Tmin=Tmin_phi; Corr.Tmax=Tmax_phi;
      distr_t F_S= Corr.Fit_distr( Corr.mel_ov_mass_t( PT2_SS_distr, ""))/2.0;
          
      distr_t_list Z_FACTOR= (1.0/(F_D*F_S))*EXPT_D(M_Ds, Corr.Nt)*EXPT_D(-1.0*M_phi, Corr.Nt)*EXP_D(M_phi*Tins[iq]);
 


      distr_t GAMMA= SQRT_D( 1.0 + kz*kz/(M_phi_rest*M_phi_rest));

      distr_t N3 = (kz/M_phi_rest)*(kz/M_phi_rest)/GAMMA;

         

    
      
      
      //GET FORM FACTORS

      distr_t_list V_distr(UseJack), A0_distr(UseJack), A1_distr(UseJack), A2_distr(UseJack), A0_distr_WI(UseJack);

      distr_t V(UseJack), A0(UseJack), A1(UseJack), A2(UseJack);


      if(fabs(kz)< 1e-14) { // phi at rest

	
	V_distr= Z_FACTOR*0.0; //set to 0
	A0_distr= Z_FACTOR*0.0; //set to 0
	A0_distr_WI = Z_FACTOR*0.0; //set to 0
	A2_distr = Z_FACTOR*0.0; //set to 0
	A1_distr= -1.0*Z_FACTOR*Za*PT3_Ax_distr/(M_Ds+M_phi_rest);
      }
      else {



	//define polarization vector
	distr_t_list eps(UseJack);
	eps.distr_list.push_back( N3*M_phi/fabs(kz) );
	eps.distr_list.push_back( N3*( (M_phi/kz)*(M_phi/kz) + 1));
	eps.distr_list.push_back( N3);
 
	double kg= kz;

	V_distr= Z_FACTOR*Zv*PT3_Vyz_distr*(M_Ds+M_phi_rest)*GAMMA/(2.0*M_Ds*kz);

	distr_t eps_dot_q = (kz>0?-1.0:1.0)*eps.distr_list[0]*M_Ds; // - M_phi) - (eps.distr_list[1] + 2*eps.distr_list[2])*kg;

	cout<<"eps * q: "<<eps_dot_q.ave()<<" , expected: "<<( (kz/M_phi_rest)*(kz/M_phi_rest)*M_Ds*M_phi/(GAMMA*fabs(kz))).ave()<<endl;

	vector<vector<distr_t>> CM(3);
	for(auto &c :CM ) {
	  for(int j=0;j<3;j++)   c.emplace_back(UseJack);
	}

	
	for(int ijack=0;ijack<Njacks;ijack++) {
	
	  Eigen::MatrixXd B(3,3);

	  double m= M_phi.distr[ijack];
	  double M= M_Ds.distr[ijack];
	  double m0 = M_phi_rest.distr[ijack];
	  double edq = eps_dot_q.distr[ijack];
	  double aq2 = aq2_distr.distr[ijack];
	  double aq0 = q0.distr[ijack];
	  double eps0 = eps.distr_list[0].distr[ijack];
	  double eps1 = eps.distr_list[1].distr[ijack];
	  double eps_yz = eps.distr_list[2].distr[ijack];

	  //kg is the photon momentum

	  
	  //########### DEFINE COEFFICIENT MATRIX #############

	  B(0,0) = (2*m0*edq/aq2)*aq0 ; B(0,1) = (M+m0)*(eps0 - edq*aq0/aq2) ; B(0,2) = -(edq/(M+m0))*( M+m - (M*M -m0*m0)*aq0/aq2);

	  B(1,0) = (2*m0*edq/aq2)*kg ; B(1,1) = (M+m0)*(eps1 - edq*kg/aq2)  ; B(1,2) = -(edq/(M+m0))*( -kg -  (M*M- m0*m0)*kg/aq2 ) ;

	  B(2,0) = (2*m0*edq/aq2)*kg ; B(2,1) = (M+m0)*(eps_yz - edq*kg/aq2) ; B(2,2) = -(edq/(M+m0))*( -kg - (M*M- m0*m0)*kg/aq2 );


	  //###################################################
	  

	  Eigen::MatrixXd Binv = B.inverse();

	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++)
	      CM[i][j].distr.push_back( Binv(i,j));
	}

	A0_distr_WI = Z_FACTOR*0.5*(mass_s+mass_c)*PT3_P5_distr/(eps_dot_q*M_phi_rest);
        A0_distr = Z_FACTOR*Za*( -1.0*CM[0][0]*PT3_A0_IM_distr  -CM[0][1]*PT3_Ax_distr -CM[0][2]*PT3_Ayz_distr);
	A1_distr = Z_FACTOR*Za*( -1.0*CM[1][0]*PT3_A0_IM_distr  -CM[1][1]*PT3_Ax_distr -CM[1][2]*PT3_Ayz_distr);
	A2_distr = Z_FACTOR*Za*( -1.0*CM[2][0]*PT3_A0_IM_distr  -CM[2][1]*PT3_Ax_distr -CM[2][2]*PT3_Ayz_distr);


      }


      //print form factors estimator to file

      Print_To_File( {}, {V_distr.ave(), V_distr.err(), A0_distr_WI.ave(), A0_distr_WI.err(),  A0_distr.ave(), A0_distr.err(), A1_distr.ave(), A1_distr.err(), A2_distr.ave(), A2_distr.err()}, "../data/Ds_phi_lnu/"+PT2_cs_data[0][0].Tag[iens]+"/FF/FF_distr_iq_"+to_string(iq), "", "# V0  A0  A1 A2");
      

      

      //fit form factors

      Corr.Tmin= 18; Corr.Tmax=25;

      if(iq==3 || iq==4) {Corr.Tmin=11; Corr.Tmax=16;}
      else if(iq==2)  {Corr.Tmin=20; Corr.Tmax=24;}    
      else if(iq==1)  {Corr.Tmin=19; Corr.Tmax=24;}
      else if(iq==0)  {Corr.Tmin=21; Corr.Tmax=24;}
      else crash("iq: "+to_string(iq)+" not found");
      
      V = Corr.Fit_distr(V_distr);

      if(iq==3 || iq == 4) {Corr.Tmin=10; Corr.Tmax=18;}
      else if(iq==2)  {Corr.Tmin=23; Corr.Tmax=28;}    
      else if(iq==1)  {Corr.Tmin=23; Corr.Tmax=28;}
      else if(iq==0)  {Corr.Tmin=21; Corr.Tmax=27;}
      else crash("iq: "+to_string(iq)+" not found");
      
      A0 = Corr.Fit_distr( A0_distr_WI);
      
      if(iq==3 || iq == 4) {Corr.Tmin=12; Corr.Tmax=17;}
      else if(iq==2)  {Corr.Tmin=18; Corr.Tmax=23;}    
      else if(iq==1)  {Corr.Tmin=19; Corr.Tmax=23;}
      else if(iq==0)  {Corr.Tmin=17; Corr.Tmax=21;}
      else crash("iq: "+to_string(iq)+" not found");
      
      A1 = Corr.Fit_distr( A1_distr);

      if(iq==3 || iq == 4) {Corr.Tmin=11; Corr.Tmax=16;}
      else if(iq==2)  {Corr.Tmin=24; Corr.Tmax=28;}    
      else if(iq==1)  {Corr.Tmin=24; Corr.Tmax=28;}
      else if(iq==0)  {Corr.Tmin=24; Corr.Tmax=28;}
      else crash("iq: "+to_string(iq)+" not found");
      
      A2 = Corr.Fit_distr( A2_distr);



      FF_V_list.distr_list.push_back(V);
      FF_A0_list.distr_list.push_back(A0);
      FF_A1_list.distr_list.push_back(A1);
      FF_A2_list.distr_list.push_back(A2);
      q2_list.push_back( q2);
      

      
    }

    Print_To_File( { }, {q2_list, FF_V_list.ave(), FF_V_list.err() ,  FF_A0_list.ave(), FF_A0_list.err(), FF_A1_list.ave(), FF_A1_list.err(),  FF_A2_list.ave(), FF_A2_list.err()}, "../data/Ds_phi_lnu/"+PT2_cs_data[0][0].Tag[iens]+"/FF/FF_list", "", "");
  
    
  }


  return;
}
