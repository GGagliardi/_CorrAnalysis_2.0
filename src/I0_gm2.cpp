#include "../include/I0_gm2.h"
#include "numerics.h"

const double alpha = 1.0 / 137.035999;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const double Njacks = 500;
const double DTT = 0.33;
const bool UseJack = true;
const vector<int> T_list({128,160});
const double qd = -1.0 / 3.0;
const double qu = 2.0 / 3.0;
const double qs= -1.0/3.0;



void Bounding_HVP_isoscalar(distr_t &amu_HVP, int &Tcut_opt,  const distr_t_list &V, const distr_t &a, string path,distr_t lowest_mass) {



  int Simps_ord=1;
  
  double T= V.size(); 

  auto LOG = [](double R_G, double t) { return log(fabs(R_G));};

  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};

  distr_t_list Ker = distr_t_list::f_of_distr(K, a , T/2);

  int T_ext_max= 300;

  distr_t_list Ker_extended= distr_t_list::f_of_distr(K, a , T_ext_max);

  int Tdatas_opt=-1;

  
  auto amu_HVP_func = [&V, &Ker, &Simps_ord](double tcut) -> distr_t {
    
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


   


  //find interval where difference between amu_HVP_min_Tdata and amu_HVP_max_Tdata is smaller than 0.2 sigma
  //average
  bool Found_Tdata_opt=false;
  int tdata_opt=1;

   
  while(!Found_Tdata_opt && tdata_opt < T/2) {

    distr_t diff_max_min = amu_HVP_max_Tdata.distr_list[tdata_opt-1] - amu_HVP_min_Tdata.distr_list[tdata_opt-1];
    if( fabs(diff_max_min.ave())/min( amu_HVP_max_Tdata.err(tdata_opt-1) , amu_HVP_min_Tdata.err(tdata_opt-1)) < 1) Found_Tdata_opt=true;
    else tdata_opt++;
  }

      

  //if tdata_opt has not been found return tdata = -1 and amu_HVP =amu_HVP(T/2)
  if(!Found_Tdata_opt) {
	
    Tdatas_opt = -1;
    amu_HVP =  amu_HVP_func(T/2);
	
  }
  else { //tdata_opt has been found


    bool fit_to_constant=true;
    //average over min max over an interval of 0.25 fm (DTT)
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

  distr_t_list amu_HVP_VEC(UseJack);
  vector<double> Is_in_fit_int;
  for(int i=0;i<(signed)TCUTS.size(); i++) {
    amu_HVP_VEC.distr_list.push_back( amu_HVP);
    Is_in_fit_int.push_back(( (i >= Tdatas_opt) && (i <= (Tdatas_opt+  DTT*fm_to_inv_Gev/a.ave())))?1.0:0.0);
  }

  Tcut_opt= Tdatas_opt;
      
  //Print to File
  Print_To_File({}, {TCUTS, amu_HVP_min_Tdata.ave(), amu_HVP_min_Tdata.err(), amu_HVP_max_Tdata.ave(), amu_HVP_max_Tdata.err(), amu_HVP_T_2.ave(), amu_HVP_T_2.err(), amu_HVP_VEC.ave(), amu_HVP_VEC.err(), Is_in_fit_int, Is_T_data_opt}, path+".bound", "", "#tcut   lower    upper     T/2  result   Is_in_fit_int?     Is_Tcut_opt.       Tcut_opt= "+to_string(Tdatas_opt));



  return;

}








void I0_gm2() {


  //define vectors containing mass parameters

  VVfloat m;
  m.push_back({0.00072, 2*0.00072, 3*0.00072, 5*0.00072, 0.0182782});
  m.push_back({0.00060, 2*0.00060, 3*0.00060, 5*0.00060, 0.016053});

  vector<string> Ens_List({"B64", "C80"});
  vector<string> Full_Tags({"cB211b.072.64", "cC211a.06.80"});
  vector<bool> Get_IB_HTP({true,true});
  vector<double> Mpi({0.1402,0.1367});
  
  VVfloat meps;
  meps.push_back({1e-6, 2*1e-6, 3*1e-6, 5*1e-6, 0.01825-0.0182782});
  meps.push_back({1e-6, 2e-6, 3e-6, 5e-6 , 0.016073 - 0.016053});  


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
  

  for(int ENS=0;ENS<(signed)Ens_List.size();ENS++) {

    //Read data

    ifstream Read_Confs("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/confs_list");
    
    vector<string> Confs_disco;
 
    while( !Read_Confs.eof()) {
      string a; 
      Read_Confs >> a;
      if(! Read_Confs.eof()) { Confs_disco.push_back(a) ; }
    }
    
    int Nconfs=Confs_disco.size();
    cout<<"Analyzing Ensemble: "<<Ens_List[ENS]<<endl;
    cout<<"Nconfs: "<<Nconfs<<endl;
    cout<<"T: "<<T_list[ENS]<<endl;
    //read isoQCD bubble

    VVVVfloat iso_ll2(3), iso_l2l3(3), iso_l3l4(3), iso_l4s(3);
    VVVVfloat IB_l(3), IB_l2(3), IB_l3(3), IB_l4(3), IB_s(3);
    
    for(int i=0;i<3;i++) {
      iso_ll2[i].resize(Confs_disco.size());
      iso_l2l3[i].resize(Confs_disco.size());
      iso_l3l4[i].resize(Confs_disco.size());
      iso_l4s[i].resize(Confs_disco.size());
      
      IB_l[i].resize(Confs_disco.size());
      IB_l2[i].resize(Confs_disco.size());
      IB_l3[i].resize(Confs_disco.size());
      IB_l4[i].resize(Confs_disco.size());
      IB_s[i].resize(Confs_disco.size());
    }
    
    vector<vector<int>> tsou_list_iso(Confs_disco.size());
    vector<vector<int>> tsou_list_IB(Confs_disco.size()); 
    
    
    for(int iconf=0;iconf<(signed)Confs_disco.size();iconf++) {
      
      ifstream Read_R0_ll2("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_ll2_R0");
      ifstream Read_R1_ll2("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_ll2_R1");
      while( !Read_R0_ll2.eof() && ! Read_R1_ll2.eof()) {
	string contr_R0, contr_R1, dummy;
	getline(Read_R0_ll2, contr_R0);
	getline(Read_R1_ll2, contr_R1);
	if(!Read_R0_ll2.eof() && !Read_R1_ll2.eof()) {
	  //extract tsou_list from contr_R0
	  unsigned pos_start= contr_R0.find_first_of("(");
	  string contr_R0_red= contr_R0.substr(pos_start+1);
	  tsou_list_iso[iconf].push_back( stoi( contr_R0_red.substr(0, contr_R0_red.find_first_of(","))));
	
		
	  for(int i=0;i<3;i++) {
	    getline(Read_R0_ll2,dummy);
	    getline(Read_R1_ll2, dummy);
	    Vfloat C;
	    for(int t=0;t<T_list[ENS];t++) {
	      double cr0,cr1;
	      
	      Read_R0_ll2 >> cr0; Read_R1_ll2 >> cr1;
	      if(!Read_R0_ll2.eof() && !Read_R1_ll2.eof()) C.push_back( -1.0*(m[ENS][1]+m[ENS][0])*(cr0+cr1));
	    }
	    iso_ll2[i][iconf].push_back(C);
	    getline(Read_R0_ll2, dummy); getline(Read_R1_ll2,dummy);
	  }
	}
      }
      Read_R0_ll2.close(); Read_R1_ll2.close();
      
      

      
      
      ifstream Read_R0_l2l3("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l2l3_R0");
      ifstream Read_R1_l2l3("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l2l3_R1");
      while( !Read_R0_l2l3.eof() && ! Read_R1_l2l3.eof()) {
	string contr_R0, contr_R1, dummy;
	
	getline(Read_R0_l2l3, contr_R0);
	getline(Read_R1_l2l3, contr_R1);
	if(!Read_R0_l2l3.eof() && !Read_R1_l2l3.eof()) {
	  for(int i=0;i<3;i++) {
	    getline(Read_R0_l2l3,dummy); getline(Read_R1_l2l3, dummy);
	    Vfloat C;
	    for(int t=0;t<T_list[ENS];t++) {
	      double cr0,cr1;
	      Read_R0_l2l3 >> cr0; Read_R1_l2l3 >> cr1;
	      if(!Read_R0_l2l3.eof() && !Read_R1_l2l3.eof()) C.push_back(  -1.0*(m[ENS][2]+m[ENS][1])*(cr0+cr1));
	    }
	    iso_l2l3[i][iconf].push_back(C);
	    getline(Read_R0_l2l3, dummy); getline(Read_R1_l2l3,dummy);
	  }
	}
      }
      Read_R0_l2l3.close(); Read_R1_l2l3.close();
      
      

    
      ifstream Read_R0_l3l4("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l3l4_R0");
      ifstream Read_R1_l3l4("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l3l4_R1");
      while( !Read_R0_l3l4.eof() && ! Read_R1_l3l4.eof()) {
      string contr_R0, contr_R1, dummy;
      
      getline(Read_R0_l3l4, contr_R0);
      getline(Read_R1_l3l4, contr_R1);
      if(!Read_R0_l3l4.eof() && !Read_R1_l3l4.eof()) {
	for(int i=0;i<3;i++) {
	  getline(Read_R0_l3l4,dummy); getline(Read_R1_l3l4, dummy);
	  Vfloat C;
	  for(int t=0;t<T_list[ENS];t++) {
	    double cr0,cr1;
	    Read_R0_l3l4 >> cr0; Read_R1_l3l4 >> cr1;
	    if(!Read_R0_l3l4.eof() && !Read_R1_l3l4.eof()) C.push_back(  -1.0*(m[ENS][3]+m[ENS][2])*(cr0+cr1));
	  }
	  iso_l3l4[i][iconf].push_back(C);
	  getline(Read_R0_l3l4, dummy); getline(Read_R1_l3l4,dummy);
	}
      }
      }
      Read_R0_l3l4.close(); Read_R1_l3l4.close();

   

    ifstream Read_R0_l4s("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l4s_R0");
    ifstream Read_R1_l4s("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_l4s_R1");
    while( !Read_R0_l4s.eof() && ! Read_R1_l4s.eof()) {
      string contr_R0, contr_R1, dummy;

      getline(Read_R0_l4s, contr_R0);
      getline(Read_R1_l4s, contr_R1);
      if(!Read_R0_l4s.eof() && !Read_R1_l4s.eof()) {
	for(int i=0;i<3;i++) {
	  getline(Read_R0_l4s,dummy); getline(Read_R1_l4s, dummy);
	  Vfloat C;
	  for(int t=0;t<T_list[ENS];t++) {
	    double cr0,cr1;
	    Read_R0_l4s >> cr0; Read_R1_l4s >> cr1;
	    if(!Read_R0_l4s.eof() && !Read_R1_l4s.eof()) C.push_back(  1.0*(m[ENS][4]+m[ENS][3])*(cr0+cr1));
	  }
	  iso_l4s[i][iconf].push_back(C);
	  getline(Read_R0_l4s, dummy); getline(Read_R1_l4s,dummy);
	}
      }
    }
    Read_R0_l4s.close(); Read_R1_l4s.close();

   

    ifstream Read_R0_IB_l("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll_R0");
    ifstream Read_R1_IB_l("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll_R1");
    while( !Read_R0_IB_l.eof() && ! Read_R1_IB_l.eof()) {
      string contr_R0, contr_R1, dummy;

      getline(Read_R0_IB_l, contr_R0);
      getline(Read_R1_IB_l, contr_R1);
      if(!Read_R0_IB_l.eof() && !Read_R1_IB_l.eof()) {
	for(int i=0;i<3;i++) {
	  getline(Read_R0_IB_l,dummy); getline(Read_R1_IB_l, dummy);
	  Vfloat C;
	  for(int t=0;t<T_list[ENS];t++) {
	    double cr0,cr1;
	    Read_R0_IB_l >> cr0; Read_R1_IB_l >> cr1;
	    if(!Read_R0_IB_l.eof() && !Read_R1_IB_l.eof()) C.push_back( ((2*m[ENS][0]+meps[ENS][0])/meps[ENS][0])*(cr0+cr1));
	  }
	  IB_l[i][iconf].push_back(C);
	  getline(Read_R0_IB_l, dummy); getline(Read_R1_IB_l,dummy);
	}
      }
    }
    Read_R0_IB_l.close(); Read_R1_IB_l.close();
    
   
    
    ifstream Read_R0_IB_s("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ss_R0");
    ifstream Read_R1_IB_s("../axial_WI_disco/isoQCD/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ss_R1");
    while( !Read_R0_IB_s.eof() && ! Read_R1_IB_s.eof()) {
      string contr_R0, contr_R1, dummy;

      getline(Read_R0_IB_s, contr_R0);
      getline(Read_R1_IB_s, contr_R1);
      if(!Read_R0_IB_s.eof() && !Read_R1_IB_s.eof()) {
	for(int i=0;i<3;i++) {
	  getline(Read_R0_IB_s,dummy); getline(Read_R1_IB_s, dummy);
	  Vfloat C;
	  for(int t=0;t<T_list[ENS];t++) {
	    double cr0,cr1;
	    Read_R0_IB_s >> cr0; Read_R1_IB_s >> cr1;
	    if(!Read_R0_IB_s.eof() && !Read_R1_IB_s.eof()) C.push_back(  ((2*m[ENS][4]+meps[ENS][4])/meps[ENS][4])*(cr0+cr1));
	  }
	  IB_s[i][iconf].push_back(C);
	  getline(Read_R0_IB_s, dummy); getline(Read_R1_IB_s,dummy);
	}
      }
    }
    Read_R0_IB_s.close(); Read_R1_IB_s.close();

    
    if(Get_IB_HTP[ENS]) {

      ifstream Read_R0_IB_l2("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll2_R0");
      ifstream Read_R1_IB_l2("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll2_R1");
      while( !Read_R0_IB_l2.eof() && ! Read_R1_IB_l2.eof()) {
	string contr_R0, contr_R1, dummy;
	getline(Read_R0_IB_l2, contr_R0);
	getline(Read_R1_IB_l2, contr_R1);

	if(!Read_R0_IB_l2.eof() && !Read_R1_IB_l2.eof()) {
	  //extract tsou_list from contr_R0
	  unsigned pos_start= contr_R0.find_first_of("(");
	  string contr_R0_red= contr_R0.substr(pos_start+1);
	  tsou_list_IB[iconf].push_back( stoi( contr_R0_red.substr(0, contr_R0_red.find_first_of(","))));
	  
	  for(int i=0;i<3;i++) {
	    getline(Read_R0_IB_l2,dummy); getline(Read_R1_IB_l2, dummy);
	    Vfloat C;
	    for(int t=0;t<T_list[ENS];t++) {
	      double cr0,cr1;
	      Read_R0_IB_l2 >> cr0; Read_R1_IB_l2 >> cr1;
	      if(!Read_R0_IB_l2.eof() && !Read_R1_IB_l2.eof()) C.push_back(  ((2*m[ENS][1]+meps[ENS][1])/meps[ENS][1])*(cr0+cr1));
	    }
	    IB_l2[i][iconf].push_back(C);
	    getline(Read_R0_IB_l2, dummy); getline(Read_R1_IB_l2,dummy);
	  }
	}
      }
      Read_R0_IB_l2.close(); Read_R1_IB_l2.close();
      
      

      
      ifstream Read_R0_IB_l3("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll3_R0");
      ifstream Read_R1_IB_l3("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll3_R1");
      while( !Read_R0_IB_l3.eof() && ! Read_R1_IB_l3.eof()) {
	string contr_R0, contr_R1, dummy;
	
	getline(Read_R0_IB_l3, contr_R0);
	getline(Read_R1_IB_l3, contr_R1);
	if(!Read_R0_IB_l3.eof() && !Read_R1_IB_l3.eof()) {
	  for(int i=0;i<3;i++) {
	    getline(Read_R0_IB_l3,dummy); getline(Read_R1_IB_l3, dummy);
	    Vfloat C;
	    for(int t=0;t<T_list[ENS];t++) {
	      double cr0,cr1;
	      Read_R0_IB_l3 >> cr0; Read_R1_IB_l3 >> cr1;
	      if(!Read_R0_IB_l3.eof() && !Read_R1_IB_l3.eof()) C.push_back(  ((2*m[ENS][2]+meps[ENS][2])/meps[ENS][2])*(cr0+cr1));
	    }
	    IB_l3[i][iconf].push_back(C);
	    getline(Read_R0_IB_l3, dummy); getline(Read_R1_IB_l3,dummy);
	  }
	}
      }
      Read_R0_IB_l3.close(); Read_R1_IB_l3.close();
      
            
      ifstream Read_R0_IB_l4("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll4_R0");
      ifstream Read_R1_IB_l4("../axial_WI_disco/IB/"+Ens_List[ENS]+"/"+Confs_disco[iconf]+"/mes_contr_OET_der_ll4_R1");
    
      while( !Read_R0_IB_l4.eof() && ! Read_R1_IB_l4.eof()) {
	string contr_R0, contr_R1, dummy;
	
	getline(Read_R0_IB_l4, contr_R0);
	getline(Read_R1_IB_l4, contr_R1);
	if(!Read_R0_IB_l4.eof() && !Read_R1_IB_l4.eof()) {
	  for(int i=0;i<3;i++) {
	    getline(Read_R0_IB_l4,dummy); getline(Read_R1_IB_l4, dummy);
	    Vfloat C;
	    for(int t=0;t<T_list[ENS];t++) {
	      double cr0,cr1;
	      Read_R0_IB_l4 >> cr0; Read_R1_IB_l4 >> cr1;
	      if(!Read_R0_IB_l4.eof() && !Read_R1_IB_l4.eof()) C.push_back(  ((2*m[ENS][3]+meps[ENS][3])/meps[ENS][3])*(cr0+cr1));
	    }
	    IB_l4[i][iconf].push_back(C);
	    getline(Read_R0_IB_l4, dummy); getline(Read_R1_IB_l4,dummy);
	  }
	}
      }
      Read_R0_IB_l4.close(); Read_R1_IB_l4.close();
    }

    }
    
    
  cout<<"correlator loaded"<<endl;


  //start the analysis

  VVVVfloat B_iso= summ_master( iso_ll2, iso_l2l3, iso_l3l4, iso_l4s) ;
  VVVVfloat B_l2_iso= summ_master( iso_l2l3, iso_l3l4, iso_l4s);
  VVVVfloat B_l3_iso = summ_master(iso_l3l4, iso_l4s);
  VVVVfloat B_l4_iso = iso_l4s;

  

  cout<<"convoluting..."<<endl;

  VVfloat C_disc(T_list[ENS]);
  VVfloat C_IB_l(T_list[ENS]);
  VVfloat C_IB_l2(T_list[ENS]);
  VVfloat C_IB_l3(T_list[ENS]);
  VVfloat C_IB_l4(T_list[ENS]);
  
  for(auto &C: C_disc) C.resize(Nconfs,0.0);
  for(auto &C: C_IB_l) C.resize(Nconfs,0.0);
  for(auto &C: C_IB_l2) C.resize(Nconfs,0.0);
  for(auto &C: C_IB_l3) C.resize(Nconfs,0.0);
  for(auto &C: C_IB_l4) C.resize(Nconfs,0.0);
  
    

  for(int k=0;k<3;k++) {
    C_disc = summ_master(C_disc, convolute(B_iso[k], B_iso[k], 1, tsou_list_iso));
    C_IB_l= summ_master( C_IB_l, convolute(B_iso[k], IB_l[k], 1, tsou_list_iso));
    C_IB_l2= summ_master( C_IB_l2, convolute_reduced(B_iso[k], IB_l2[k], 1,tsou_list_iso, tsou_list_IB));
    C_IB_l3= summ_master( C_IB_l3, convolute_reduced(B_iso[k], IB_l3[k], 1,tsou_list_iso, tsou_list_IB));
    C_IB_l4= summ_master( C_IB_l4, convolute_reduced(B_iso[k], IB_l4[k], 1,tsou_list_iso, tsou_list_IB));
  }

  cout<<"convolution done!"<<endl;

  CorrAnalysis Corr(UseJack,Njacks,100);
  Corr.Perform_Nt_t_average=0;
  Corr.Nt=T_list[ENS];

  double dmd=(2.56-0.196)*1e-4;
  double dmu=(-2.56-0.196)*1e-4;

  
  int T= Corr.Nt; 
  int L = Corr.Nt/2;
  double F= ( (double)( pow(L,3)*pow(T,2)))*(qs/2)*(qs/2);
  double F_libe= 2*( (double)( pow(L,3)*pow(T,2)))*(qs/2)*((qu*dmu + qd*dmd)/2); 

  distr_t_list C_disc_distr= Corr.corr_t(Multiply_Vvector_by_scalar(C_disc,F/3.0), "../data/axial_WI_disco/V_4OET_"+Ens_List[ENS]+"");
  distr_t_list C_IB_l_distr= Corr.corr_t(Multiply_Vvector_by_scalar(C_IB_l,F_libe/3.0), "../data/axial_WI_disco/V_IB_l_"+Ens_List[ENS]+"");
  distr_t_list C_IB_l2_distr= Corr.corr_t(Multiply_Vvector_by_scalar(C_IB_l2,F_libe/3.0), "../data/axial_WI_disco/V_IB_l2_"+Ens_List[ENS]+"");
  distr_t_list C_IB_l3_distr= Corr.corr_t(Multiply_Vvector_by_scalar(C_IB_l3,F_libe/3.0), "../data/axial_WI_disco/V_IB_l3_"+Ens_List[ENS]+"");
  distr_t_list C_IB_l4_distr= Corr.corr_t(Multiply_Vvector_by_scalar(C_IB_l4,F_libe/3.0), "../data/axial_WI_disco/V_IB_l4_"+Ens_List[ENS]+"");

  distr_t a_distr(UseJack);

  

  GaussianMersenne GM(43112);
  double a_ave, a_err,Zv;
  if(Ens_List[ENS] == "B64") {
    a_ave= 0.07948*fm_to_inv_Gev;
    a_err =0.00011*fm_to_inv_Gev;
    Zv= 0.70637654;
  }
  else if(Ens_List[ENS] == "C80") {
    a_ave= 0.06819*fm_to_inv_Gev;
    a_err =0.00014*fm_to_inv_Gev;
    Zv= 0.72540536;
  }
  else crash("Ens not found");
 
 
  for(int ijack=0;ijack<Njacks;ijack++) a_distr.distr.push_back( a_ave + GM()*a_err/sqrt(Njacks-1.0));

  //evaluate amu
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);


  distr_t_list amu_disc(UseJack), amu_IB_l(UseJack), amu_IB_l2(UseJack), amu_IB_l3(UseJack), amu_IB_l4(UseJack);


  distr_t p_sum_disc= 0.0*Get_id_jack_distr(Njacks);

  distr_t p_sum_IB_l =  0.0*Get_id_jack_distr(Njacks);
  distr_t p_sum_IB_l2 =  0.0*Get_id_jack_distr(Njacks);
  distr_t p_sum_IB_l3 =  0.0*Get_id_jack_distr(Njacks);
  distr_t p_sum_IB_l4 =  0.0*Get_id_jack_distr(Njacks);
  
 

  for(int t=0;t<T/2;t++) {

    p_sum_disc = p_sum_disc + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*C_disc_distr.distr_list[t]*Ker.distr_list[t];
    amu_disc.distr_list.push_back( p_sum_disc);

    
    p_sum_IB_l = p_sum_IB_l + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*C_IB_l_distr.distr_list[t]*Ker.distr_list[t];
    amu_IB_l.distr_list.push_back( p_sum_IB_l);

    
    p_sum_IB_l2 = p_sum_IB_l2 + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*C_IB_l2_distr.distr_list[t]*Ker.distr_list[t];
    amu_IB_l2.distr_list.push_back( p_sum_IB_l2);

    p_sum_IB_l3 = p_sum_IB_l3 + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*C_IB_l3_distr.distr_list[t]*Ker.distr_list[t];
    amu_IB_l3.distr_list.push_back( p_sum_IB_l3);

    p_sum_IB_l4 = p_sum_IB_l4 + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*C_IB_l4_distr.distr_list[t]*Ker.distr_list[t];
    amu_IB_l4.distr_list.push_back( p_sum_IB_l4);
    
  }


  Print_To_File({}, {amu_disc.ave(), amu_disc.err()}, "../data/axial_WI_disco/amu_disc_"+Ens_List[ENS]+"","","");
  Print_To_File({}, {amu_IB_l.ave(), amu_IB_l.err()}, "../data/axial_WI_disco/amu_IB_l_"+Ens_List[ENS]+"","","");

  
  Print_To_File({}, {amu_IB_l2.ave(), amu_IB_l2.err()}, "../data/axial_WI_disco/amu_IB_l2_"+Ens_List[ENS]+"","","");
  Print_To_File({}, {amu_IB_l3.ave(), amu_IB_l3.err()}, "../data/axial_WI_disco/amu_IB_l3_"+Ens_List[ENS]+"","","");
  Print_To_File({}, {amu_IB_l4.ave(), amu_IB_l4.err()}, "../data/axial_WI_disco/amu_IB_l4_"+Ens_List[ENS]+"","","");
  

  //do bounding

  //load amus and amu_l correlators
  VVfloat C_OS;
  vector<string> Corr_tags({"VkVk_OS.bin"});
  for(int id=0; id<(signed)Corr_tags.size(); id++) {
    string ch= Corr_tags[id];
    FILE *stream = fopen(("../../PEAKY_BLINDER/confs_to_blind/"+Ens_List[ENS]+"/"+ch).c_str(), "rb");
    size_t Nconfs_C, TT, Nhits, Nsubs;
    bin_read(Nconfs_C, stream);
    bin_read(Nhits, stream);
    bin_read(TT, stream);
    bin_read(Nsubs, stream);

    C_OS.resize(TT);
    for(auto & cc: C_OS) cc.resize(Nconfs_C,0);
   
    cout<<"Nconfs_C: "<<Nconfs_C<<endl;
    cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
    cout<<"Nhits: "<<Nhits<<endl;
    cout<<"Nsubs: "<<Nsubs<<endl;
      
    for(size_t iconf=0;iconf<Nconfs_C;iconf++) {
      for(size_t t=0;t<TT/2+1;t++) {
	for(size_t is=0;is<Nsubs;is++) {
	  double c;
	  bin_read(c, stream);
	  C_OS[t][iconf] += c/Nsubs;
	}
	//symmetrize
	if( t != 0) {
	  C_OS[TT-t][iconf] = C_OS[t][iconf];
	  
	}
      }
    }
    
    
    
    
    fclose(stream);
    
  }



  distr_t_list C_OS_l= Corr.corr_t( C_OS, "../data/axial_WI_disco/C_OS_conn_l_"+Ens_List[ENS]);

  data_t  Vk_s_data_OS;

  Vk_s_data_OS.Read("../HVP_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK");

  int iens=0;
  while( Vk_s_data_OS.Tag[iens] != Full_Tags[ENS]) {iens++;}

  distr_t_list C_OS_s = Corr.corr_t( Vk_s_data_OS.col(0)[iens], "../data/axial_WI_disco/C_OS_conn_s_"+Ens_List[ENS]);

  distr_t_list C_bound = 1e10*Zv*Zv*( ((qu*qu+qd*qd)/10.0)*C_OS_l + (qs*qs)*C_OS_s + C_disc_distr);

  
  double ampi= Mpi[ENS]*a_distr.ave();
  double k= 2*3.14159/L;

  distr_t Mpi3= (2*sqrt( pow(ampi,2) + k*k) + sqrt( pow(ampi,2) + 0.0*k*k))*Get_id_jack_distr(Njacks); 

  
  distr_t amu_disc_bounding(UseJack);

  int tcut_opt;
 
  Bounding_HVP_isoscalar(amu_disc_bounding, tcut_opt,  C_bound, a_distr, "../data/axial_WI_disco/isoscalar_"+Ens_List[ENS], Mpi3);


  cout<<"Ens: "<<Ens_List[ENS]<<" amu(I=0): "<<amu_disc_bounding.ave()<<" "<<amu_disc_bounding.err()<<endl;

  }
  
  
  
    
  


  return ;
}
