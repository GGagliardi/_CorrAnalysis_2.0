#include "../include/vph_Nissa.h"

using namespace std;

const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nboots= 100;
const bool UseJack=1;
const int nboots=150;
const int Njacks=30;
const double qu = 2.0/3.0; //electric charge of u-type quark
const double qd = -1.0/3.0; //electric charge of d-type quark
const string Meson="K";


void Compute_form_factors_Nissa() {

  //create directories
  boost::filesystem::create_directory("../data/vph_Nissa");
  boost::filesystem::create_directory("../data/vph_Nissa/C");
  boost::filesystem::create_directory("../data/vph_Nissa/mass");


  


  //axial
  vector<vector<data_t>> C_A_Fu_data(4), C_A_Fd_data(4), C_A_Bu_data(4), C_A_Bd_data(4);
  //vector
  vector<vector<data_t>> C_V_Fu_data(4), C_V_Fd_data(4), C_V_Bu_data(4), C_V_Bd_data(4);


  //2pts
  data_t data_2pts, data_2pts_test, data_2pts_pion;

  for(int mu=0;mu<4;mu++) {

    C_A_Fu_data[mu].resize(4);
    C_A_Fd_data[mu].resize(4);
    C_A_Bu_data[mu].resize(4);
    C_A_Bd_data[mu].resize(4);
    C_V_Fu_data[mu].resize(4);
    C_V_Fd_data[mu].resize(4);
    C_V_Bu_data[mu].resize(4);
    C_V_Bd_data[mu].resize(4);

  }

 
  //Read data

  //2pts function

  data_2pts.Read("../new_vph_gpu_data", "mes_contr_2pts_3", "P5P5");
  data_2pts_test.Read("../test_2pt_K_nissa", "mes_contr_2pts_3", "P5P5");
  data_2pts_pion.Read("../test_2pt_PI_nissa", "mes_contr_2pts_1", "P5P5");

  //vector

  
  //Fu
  C_V_Fu_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "V0P5");
  C_V_Fu_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "V0P5");
  C_V_Fu_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "V0P5");
  C_V_Fu_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "V0P5");
  C_V_Fu_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "V1P5");
  C_V_Fu_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "V1P5");
  C_V_Fu_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "V1P5");
  C_V_Fu_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "V1P5");
  C_V_Fu_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "V2P5");
  C_V_Fu_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "V2P5");
  C_V_Fu_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "V2P5");
  C_V_Fu_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "V2P5");
  C_V_Fu_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "V3P5");
  C_V_Fu_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "V3P5");
  C_V_Fu_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "V3P5");
  C_V_Fu_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "V3P5");

  //Fd
  C_V_Fd_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "V0P5");
  C_V_Fd_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "V0P5");
  C_V_Fd_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "V0P5");
  C_V_Fd_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "V0P5");
  C_V_Fd_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "V1P5");
  C_V_Fd_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "V1P5");
  C_V_Fd_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "V1P5");
  C_V_Fd_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "V1P5");
  C_V_Fd_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "V2P5");
  C_V_Fd_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "V2P5");
  C_V_Fd_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "V2P5");
  C_V_Fd_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "V2P5");
  C_V_Fd_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "V3P5");
  C_V_Fd_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "V3P5");
  C_V_Fd_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "V3P5");
  C_V_Fd_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "V3P5");

  //Bu
  C_V_Bu_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "V0P5");
  C_V_Bu_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "V0P5");
  C_V_Bu_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "V0P5");
  C_V_Bu_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "V0P5");
  C_V_Bu_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "V1P5");
  C_V_Bu_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "V1P5");
  C_V_Bu_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "V1P5");
  C_V_Bu_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "V1P5");
  C_V_Bu_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "V2P5");
  C_V_Bu_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "V2P5");
  C_V_Bu_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "V2P5");
  C_V_Bu_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "V2P5");
  C_V_Bu_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "V3P5");
  C_V_Bu_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "V3P5");
  C_V_Bu_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "V3P5");
  C_V_Bu_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "V3P5");

  //Bd
  C_V_Bd_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "V0P5");
  C_V_Bd_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "V0P5");
  C_V_Bd_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "V0P5");
  C_V_Bd_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "V0P5");
  C_V_Bd_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "V1P5");
  C_V_Bd_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "V1P5");
  C_V_Bd_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "V1P5");
  C_V_Bd_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "V1P5");
  C_V_Bd_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "V2P5");
  C_V_Bd_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "V2P5");
  C_V_Bd_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "V2P5");
  C_V_Bd_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "V2P5");
  C_V_Bd_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "V3P5");
  C_V_Bd_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "V3P5");
  C_V_Bd_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "V3P5");
  C_V_Bd_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "V3P5");


  //axial 

  
  //Fu
  C_A_Fu_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "A0P5");
  C_A_Fu_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "A0P5");
  C_A_Fu_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "A0P5");
  C_A_Fu_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "A0P5");
  C_A_Fu_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "A1P5");
  C_A_Fu_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "A1P5");
  C_A_Fu_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "A1P5");
  C_A_Fu_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "A1P5");
  C_A_Fu_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "A2P5");
  C_A_Fu_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "A2P5");
  C_A_Fu_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "A2P5");
  C_A_Fu_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "A2P5");
  C_A_Fu_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_FF_u", "A3P5");
  C_A_Fu_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_FF_u", "A3P5");
  C_A_Fu_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_FF_u", "A3P5");
  C_A_Fu_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_FF_u", "A3P5");

  //Fd
  C_A_Fd_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "A0P5");
  C_A_Fd_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "A0P5");
  C_A_Fd_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "A0P5");
  C_A_Fd_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "A0P5");
  C_A_Fd_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "A1P5");
  C_A_Fd_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "A1P5");
  C_A_Fd_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "A1P5");
  C_A_Fd_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "A1P5");
  C_A_Fd_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "A2P5");
  C_A_Fd_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "A2P5");
  C_A_Fd_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "A2P5");
  C_A_Fd_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "A2P5");
  C_A_Fd_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_FF_d", "A3P5");
  C_A_Fd_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_FF_d", "A3P5");
  C_A_Fd_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_FF_d", "A3P5");
  C_A_Fd_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_FF_d", "A3P5");

  //Bu
  C_A_Bu_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "A0P5");
  C_A_Bu_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "A0P5");
  C_A_Bu_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "A0P5");
  C_A_Bu_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "A0P5");
  C_A_Bu_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "A1P5");
  C_A_Bu_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "A1P5");
  C_A_Bu_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "A1P5");
  C_A_Bu_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "A1P5");
  C_A_Bu_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "A2P5");
  C_A_Bu_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "A2P5");
  C_A_Bu_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "A2P5");
  C_A_Bu_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "A2P5");
  C_A_Bu_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_BB_u", "A3P5");
  C_A_Bu_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_BB_u", "A3P5");
  C_A_Bu_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_BB_u", "A3P5");
  C_A_Bu_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_BB_u", "A3P5");

  //Bd
  C_A_Bd_data[0][0].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "A0P5");
  C_A_Bd_data[1][0].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "A0P5");
  C_A_Bd_data[2][0].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "A0P5");
  C_A_Bd_data[3][0].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "A0P5");
  C_A_Bd_data[0][1].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "A1P5");
  C_A_Bd_data[1][1].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "A1P5");
  C_A_Bd_data[2][1].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "A1P5");
  C_A_Bd_data[3][1].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "A1P5");
  C_A_Bd_data[0][2].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "A2P5");
  C_A_Bd_data[1][2].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "A2P5");
  C_A_Bd_data[2][2].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "A2P5");
  C_A_Bd_data[3][2].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "A2P5");
  C_A_Bd_data[0][3].Read("../new_vph_gpu_data", "C_mu_0_BB_d", "A3P5");
  C_A_Bd_data[1][3].Read("../new_vph_gpu_data", "C_mu_1_BB_d", "A3P5");
  C_A_Bd_data[2][3].Read("../new_vph_gpu_data", "C_mu_2_BB_d", "A3P5");
  C_A_Bd_data[3][3].Read("../new_vph_gpu_data", "C_mu_3_BB_d", "A3P5");


  int Nens = data_2pts.size;

  //define lambda function to combine FF and BB

 

  for(int iens=0;iens<Nens;iens++) {

    boost::filesystem::create_directory("../data/vph_Nissa/C/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/vph_Nissa/mass/"+data_2pts.Tag[iens]);

    cout<<"Analyzing ensemble: "<<data_2pts.Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = data_2pts.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;
    distr_t_list pt2_distr= Corr.corr_t(data_2pts.col(0)[iens], "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/corr_2pt.dat");
    distr_t_list pt2_test_distr= Corr.corr_t(data_2pts_test.col(0)[0], "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/test_corr_2pt.dat");
    distr_t_list pt2_pion_distr= Corr.corr_t(data_2pts_pion.col(0)[0], "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/pion_corr_2pt.dat");
    distr_t_list eff_mass = Corr.effective_mass_t(pt2_distr, "../data/vph_Nissa/mass/"+data_2pts.Tag[iens]+"/eff_mass.dat");

    vector<vector<distr_t_list>> Ax_tens(4);
    vector<vector<distr_t_list>> Vec_tens(4);
    

    //define vectors to combine FF and BB
    Vfloat th_FF(Corr.Nt);
    Vfloat th_BB(Corr.Nt);
    for(int t=0;t<Corr.Nt;t++) {
      if(t<Corr.Nt/2) {th_FF[t]=1.0; th_BB[t] = 0.0;}
      else {th_BB[t] = 1.0; th_FF[t] = 0.0;}
    }

    //loop over mu and nu
    for(int mu=0;mu<4;mu++) {
      for(int nu=0;nu<4;nu++) {

	cout<<"Analyzing mu: "<<mu<<" nu: "<<nu<<endl;
	int Im_Re;
	double parity;

	//vector
	
	Corr.Reflection_sign = -1;
	Im_Re=1;
	parity=1.0;

	Corr.Perform_Nt_t_average = 0;
	distr_t_list vec_F_u = qu*Corr.corr_t(C_V_Fu_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list vec_B_u = qu*Corr.corr_t(C_V_Bu_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list vec_F_d = qd*Corr.corr_t(C_V_Fd_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list vec_B_d = qd*Corr.corr_t(C_V_Bd_data[mu][nu].col(Im_Re)[iens],"");
	
	distr_t_list vec_u = th_FF*vec_F_u + th_BB*vec_B_u;
	distr_t_list vec_d = th_FF*vec_F_d + th_BB*vec_B_d;
	distr_t_list vec = vec_u + vec_d;
	distr_t_list vec_symm= vec;
	//symmetrize vec
	for(int t=0; t<Corr.Nt;t++) vec_symm.distr_list[t] = 0.5*(vec.distr_list[t] + Corr.Reflection_sign*vec.distr_list[( Corr.Nt -t)%Corr.Nt]);
	Corr.Perform_Nt_t_average=1;
	
	

	//axial
	if( (mu==0 || nu==0) && (mu != 0 || nu != 0)) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	else { Im_Re=0; Corr.Reflection_sign=1; parity=1.0;}

	Corr.Perform_Nt_t_average=0;
	distr_t_list ax_F_u = qu*Corr.corr_t(C_A_Fu_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list ax_B_u = qu*Corr.corr_t(C_A_Bu_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list ax_F_d = qd*Corr.corr_t(C_A_Fd_data[mu][nu].col(Im_Re)[iens],"");
	distr_t_list ax_B_d = qd*Corr.corr_t(C_A_Bd_data[mu][nu].col(Im_Re)[iens],"");

	
	distr_t_list ax_u = parity*(ax_F_u*th_FF  + ax_B_u*th_BB);
	distr_t_list ax_d = parity*(ax_F_d*th_FF  + ax_B_d*th_BB);
	distr_t_list ax = ax_u -ax_d;
	distr_t_list ax_symm=ax;
	//symmetrize ax
	for(int t=0; t<Corr.Nt;t++) ax_symm.distr_list[t] = 0.5*(ax.distr_list[t] + Corr.Reflection_sign*ax.distr_list[( Corr.Nt -t)%Corr.Nt]);
	Corr.Perform_Nt_t_average=1;

	//restore standard reflection sign
	Corr.Reflection_sign=1;

	//single contributions
	Print_To_File({}, {ax_u.ave() ,ax_u.err(), ax_d.ave(), ax_d.err()}, "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "#   ax_u    ax_d");
	Print_To_File({}, {vec_u.ave() ,vec_u.err(), vec_d.ave(), vec_d.err()}, "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "#   ax_u    ax_d");

	//total contribution without symmetrization

	Print_To_File({}, {ax.ave(), ax.err(), vec.ave(), vec.err()}, "../data/vph_Nissa/C/wo_Tsymm_"+data_2pts.Tag[iens]+"/mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat" , "", "#  t    A     V");

	//total contribution
	Print_To_File({}, {ax_symm.ave(), ax_symm.err(), vec_symm.ave(), vec_symm.err()}, "../data/vph_Nissa/C/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat" , "", "#  t    A     V");

	//push back
	Ax_tens[mu].push_back(ax_symm);
	Vec_tens[mu].push_back(vec_symm);
	
	

      }
    }




  }
  
  
 
  



  


  return;




}
