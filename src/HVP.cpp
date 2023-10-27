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
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;

  for(int ijack=0;ijack<Njacks;ijack++) {

  a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
  a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
  a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
  a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
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

    //cout<<"#### "<<Vk_data_tm.Tag[iens]<<" ###"<<endl;
    //cout<<"HVP tm: "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<" stat. "<< (amu_HVP_tm.err()*100/amu_HVP_tm.ave())<<"%"<<endl;
    //cout<<"HVP OS: "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<" stat. "<< (amu_HVP_OS.err()*100/amu_HVP_OS.ave())<<"%"<<endl;
    //cout<<"#######"<<endl;

    auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
    distr_t_list Ker_el = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(0.0005/0.10565837) , Corr.Nt);
    distr_t_list Ker_mu = 4*pow(alpha,2)*distr_t_list::f_of_distr(K, a_distr*(1.77/0.1056837) , Corr.Nt);

    for(int t=1;t<Corr.Nt/2;t++) cout<<t<<" "<<(Ker_el*Vk_tm_distr).ave(t)<<" "<<(Ker_el*Vk_tm_distr).err(t)<<" "<<(Ker_mu*Vk_tm_distr).ave(t)<<" "<<(Ker_mu*Vk_tm_distr).err(t)<<endl;



    exit(-1);
       
  }


  return;

}
