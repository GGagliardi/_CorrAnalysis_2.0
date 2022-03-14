#include "../include/gm2_fits.h"

using namespace std;

bool print_par_info=1;
const bool verb=0;
const double Mp_phys=0.135;
const double fp_phys= 0.1304;
const double csi_phys= pow(Mp_phys/(4.0*M_PI*fp_phys),2);
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;
const double w0_scale= 1.0/0.294;

class W_ipar {

public:
  W_ipar() : w_val(0.0), w_err(0.0){}

  
  double Mp, Mp_OS;
  double L;
  double w_val, w_err;
  double fp;
  double ibeta;
  double a;
  bool Is_tm;
   
};

class W_fpar {

public:
  W_fpar() {}
  W_fpar(const Vfloat &par)  {
    if((signed)par.size() != 12) crash("In class w_fpar constructor w_fpar(vector<double>) called with vector.size != 14 ");
    w0=par[0];
    Am=par[1];
    Plog=par[2];
    Al1=par[3];
    Al2_tm=par[4];
    Al2_OS=par[5];
    D_tm=par[6];
    D_OS=par[7];
    n=par[8];
    m=par[9];
    D4_tm= par[10];
    D4_OS= par[11];
       
    
  }

  double w0, Am, Plog,  Al1,Al2_tm, Al2_OS, D_tm, D_OS,n,m, D4_tm, D4_OS;
};







void Perform_Akaike_fits(const distr_t_list &meas_tm,const distr_t_list &meas_OS, const distr_t &a_A, const distr_t &a_B, const distr_t &a_C, const distr_t &a_D, Vfloat &L_list,const distr_t_list &a_distr_list, const distr_t_list &Mpi_fit, const distr_t_list &fp_fit, vector<string> &Ens_Tag, bool UseJack, int Njacks, int Nboots,  string W_type, string channel, vector<string> &Inc_a2_list, vector<string> &Inc_FSEs_list, vector<string> &Inc_a4_list, vector<string> &Inc_mass_extr_list, vector<string> &Inc_single_fit_list, VPfloat &Inc_log_corr_list, bool allow_a4_and_log, bool allow_only_finest, int vol_mult, int mass_mult, double cont_guess ) {


  //lambda functions to be used later
  auto FVE = [](double x) -> double {return (1.0/pow(x,1.5))*exp(-x);};
  auto LOG = [](double x) { return log(x);};

  //number of Ensembles
  int Nens = Ens_Tag.size();

  distr_t_list fp_CDH(UseJack), Mpi_CDH(UseJack), Csi_CDH(UseJack);

  cout<<"Fitting with Nens: "<<Nens<<endl;


  //compute CDH corrections to Mpi and fpi
  for(int iens=0; iens<Nens;iens++) {
    //GL and CDH correction for Mpi and fp
    distr_t csi_L = Mpi_fit.distr_list[iens]*Mpi_fit.distr_list[iens]/(pow(4.0*M_PI,2)*fp_fit.distr_list[iens]*fp_fit.distr_list[iens]);
    distr_t g1 = distr_t::f_of_distr(g1_l, Mpi_fit.distr_list[iens]*L_list[iens]);
    distr_t g2 = distr_t::f_of_distr(g2_l, Mpi_fit.distr_list[iens]*L_list[iens]);
    distr_t log_l = log(csi_phys) - distr_t::f_of_distr(LOG, csi_L);
    fp_CDH.distr_list.push_back(fp_fit.distr_list[iens]/(1.0 -2.0*csi_L*g1 +2.0*csi_L*csi_L*( (Cf1(l1ph,l2ph,l3ph,l4ph) + Sf1(s0,s1,s2,s3) + Cf1_log()*log_l)*g1 + (Cf2(l1ph,l2ph,l3ph,l4ph) + Sf2(s0,s1,s2,s3) + Cf2_log()*log_l)*g2)));
    Mpi_CDH.distr_list.push_back(Mpi_fit.distr_list[iens]/(1.0 + 0.5*csi_L*g1 - csi_L*csi_L*( (Cm1(l1ph,l2ph,l3ph,l4ph) + Sm1(s0,s1,s2,s3) + Cm1_log()*log_l)*g1 + (Cm2(l1ph,l2ph,l3ph,l4ph) + Sm2(s0,s1,s2,s3) + Cm2_log()*log_l)*g2)));
    Csi_CDH.distr_list.push_back( Mpi_CDH.distr_list[iens]*Mpi_CDH.distr_list[iens]/(16.0*M_PI*M_PI*fp_CDH.distr_list[iens]*fp_CDH.distr_list[iens]));
  }

 
  
  //total number of fits
  int Nfits = 0;

  Vfloat Aka_weight, Aka_weight_red, Ch2;
  Vint Npars, Ndof;

  //Contains the info on the continuum/thermodynamic and physical point extrapolated results for all fits
  Vfloat val, err;


  vector<string> which_lattice_spacings_list;
  if(allow_only_finest) which_lattice_spacings_list = {"all", "three_finest"};
  else which_lattice_spacings_list = {"all"};

  int num_ens_A_lat=0;
  for( auto &ens_t: Ens_Tag)
    if (ens_t.substr(1,1)=="A") num_ens_A_lat++;

  //store info about the allowed fits
  vector<string> FSEs_list, a2_list, a4_list, mass_extr_list, sin_fit_list, only_finest_list;
  VPfloat n_m_list;
  for (auto &which_lat: which_lattice_spacings_list)
    for (auto &t_a2: Inc_a2_list)
      for(auto &t_FSEs : Inc_FSEs_list)
	for(auto &t_a4: Inc_a4_list)
	  for(auto &t_mass: Inc_mass_extr_list)
	    for(auto &t_sin_fit: Inc_single_fit_list) {
	      
	      Vint m_single;
	      Vint n_single;
	      
	      for( auto& n_m_pair: Inc_log_corr_list) {
	      

		bool Fit_allowed= true; 

		/* check if a given combination of fit parameters is allowed
		   ##############################################################

		   WE DO NOT ALLOW:
	     
		   1) t_a4 != "off" if (n,m) != (0,0) if allow_a4_and_log = false

		   2) If single fit = "tm" (resp. "OS"), then a4 must be also also "tm" (resp. "OS") or  "off"

		   3) if single fit = "tm" (resp. "OS"), then FSEs cannot be "OS" (resp. "tm") 

		   4) If single fit = "tm" (resp. "OS") then m (resp. n) must be zero

		   5) Npars >= Nmeas

		   6) If single fit = "tm" (resp."OS") and FSEs = "off" then change Nmeas avoiding double counting of ensembles differing only by lattice volume

		   7) If single fit = "off" and FSEs = "off" then change Nmeas avoiding double counting of ensembles differing only by lattice volume

		   13) If single fit = "off" and FSEs = "tm" or "OS" then change Nmeas avoiding double counting of ensembles differing only by lattice volume

		   14) If single fit = "tm" (resp."OS") and t_mass = "off" then change Nmeas avoiding double counting of ensembles differing only by pion mass

		   15) If single fit = "off"  and t_mass = "off" then change Nmeas avoiding double counting of ensembles differing only by pion mass

		   8) if t_a2 = "off" && t_a4 != "off" 

		   9) if t_a2 = "off" && (n != 0 || m != 0)

		   10) if t_a2 = "tm" (resp. OS) && t_a4 = "OS" (resp. "tm")  || t_a4 = "on" 

		   11) if single_fit=="off" &&  t_a2 = "tm" (resp. OS) && m != 0 (resp. n != 0)

		   12) if single fit = "tm" (resp. "OS"), then a2 must be also "tm" (resp. "OS") or "off"

		   ############################################################## */



		//RULE 1
		if(!allow_a4_and_log) {
		  if( t_a4 != "off" && ( n_m_pair.first != 0 || n_m_pair.second != 0)) Fit_allowed=false;
		}

		//RULE 2
		if( t_sin_fit == "tm" && ( t_a4 != "tm" && t_a4 != "off")) Fit_allowed=false;
		if( t_sin_fit == "OS" && ( t_a4 != "OS" && t_a4 != "off")) Fit_allowed=false;

		//RULE 3
		if( t_sin_fit == "tm" && t_FSEs == "OS") Fit_allowed=false;
		if( t_sin_fit == "OS" && t_FSEs == "tm") Fit_allowed=false;

		//RULE 4 
		if( t_sin_fit == "tm" && (find(n_single.begin(), n_single.end(), n_m_pair.first) != n_single.end())) Fit_allowed = false; 
		if( t_sin_fit == "OS" && (find(m_single.begin(), m_single.end(), n_m_pair.second) != m_single.end())) Fit_allowed = false;

		//RULE 8
		if(t_a2 == "off" && t_a4 != "off") Fit_allowed=false;


		//RULE 9
		if(t_a2 == "off" && ( n_m_pair != make_pair(0.0,0.0))) Fit_allowed=false;

		//RULE 10
		if( t_a2 == "tm" && (t_a4 == "OS" || t_a4 == "on")) Fit_allowed=false;
		if( t_a2 == "OS" && (t_a4 == "tm" || t_a4 == "on")) Fit_allowed=false;

		//RULE 11
		if( t_sin_fit=="off" && (  t_a2=="tm" && n_m_pair.second != 0  )) Fit_allowed=false;
		if( t_sin_fit=="off" && (  t_a2=="OS" && n_m_pair.first != 0   )) Fit_allowed=false;

		//RULE 12
		if(t_sin_fit=="tm" && (t_a2 != "tm" && t_a2 != "off")) Fit_allowed=false;
		if(t_sin_fit=="OS" && (t_a2 != "OS" && t_a2 != "off")) Fit_allowed=false;


		//RULE 5
		int Nens_eff = (which_lat == "three_finest")?(Nens-num_ens_A_lat):Nens;
		int Nmeas = (t_sin_fit != "off")?Nens_eff:2*Nens_eff;
		int Tot_Nmeas = Nmeas;
		if(t_sin_fit != "off" && t_FSEs == "off") Nmeas -= vol_mult;  //RULE 6
		if(t_sin_fit == "off" && t_FSEs == "off") Nmeas -=2*vol_mult; //RULE 7
		if(t_sin_fit == "off" && (t_FSEs == "tm" || t_FSEs == "OS")) Nmeas -=vol_mult; //RULE 13
		if(t_sin_fit != "off" && t_mass == "off") Nmeas -= mass_mult; //RULE 14
		if(t_sin_fit == "off" && (t_mass == "off")) Nmeas -= 2*mass_mult; //RULE 15
	    
		int npars= 1; //for the continuum extrapolated value
		if( t_sin_fit != "off") {
		  npars += (t_a2 != "off"); //O(a^2)
		  npars += (t_mass == "on"); //(Mpi dep)
		  npars += (t_a4 != "off"); //a^4 term
		  npars += (t_FSEs != "off"); //FSEs
		}
		else {  // t_sin_fit == "off"
		  npars += 2*(t_a2 == "on"); //O(a^2) D_tm, D_OS
		  npars += (t_a2 == "tm" || t_a2 == "OS"); //O(a^2) on either tm or OS
		  npars += (t_mass == "on"); //(Mpi_dep)
		  npars += 2*(t_a4=="on"); //O(a^4) on both tm and OS
		  npars += (t_a4=="tm" || t_a4=="OS"); //O(a^4) on either tm or OS
		  npars += (t_FSEs != "off"); //FSEs; 
		}
		if(npars >= Nmeas) Fit_allowed= false;

	    


		if(Fit_allowed) {

	    

		  //push_back info about the allowed fit
		  //#####################################
		  if(t_sin_fit =="tm") n_single.push_back(n_m_pair.first);
		  if(t_sin_fit =="OS") m_single.push_back(n_m_pair.second);
		  FSEs_list.push_back(t_FSEs);
		  a2_list.push_back(t_a2);
		  a4_list.push_back(t_a4);
		  mass_extr_list.push_back(t_mass);
		  sin_fit_list.push_back(t_sin_fit);
		  n_m_list.push_back( n_m_pair);
		  only_finest_list.push_back( (which_lat=="three_finest")?"on":"off");
		  //#####################################
	      
		  //########################################################
  
		  //perform physical point + continuum + thermodynamic limit
  
		  bootstrap_fit<W_fpar,W_ipar> bf(UseJack?Njacks:Nboots);
		  bf.set_warmup_lev(1);
		  bf.Set_number_of_measurements(Tot_Nmeas);
		  bf.Set_verbosity(verb);
		  bf.Set_print_path("chi2_comb_W_gm2_"+W_type+"_"+channel+".out");
  
		  bf.Add_par("w0", cont_guess, cont_guess/100.0);
		  bf.Add_par("Am", -3.0, 0.1);
		  bf.Add_par("Plog", -1.0, 0.1);
		  bf.Add_par("Al1", 1.0, 10);
		  bf.Add_par("Al2_tm", -500.0, 10);
		  bf.Add_par("Al2_OS", -500.0, 10);
		  bf.Add_par("D_tm", -1.0, 0.1);
		  bf.Add_par("D_OS", 1.0, 0.1);
		  bf.Add_par("n", 3.0, 0.1);
		  bf.Add_par("m", 3.0, 0.1);
		  bf.Add_par("D4_tm", 1.0, 0.1);
		  bf.Add_par("D4_OS", 1.0, 0.1);

		  //Fix parameters depending on the fit type

 
		  bf.Fix_par("Plog", 0.0); //we never fit log corrections xi*log(xi)
		  if(t_mass == "off") bf.Fix_par("Am", 0.0);
		  if(t_FSEs != "on") bf.Fix_par("Al1", 0.0);
		  if(t_FSEs != "tm" || t_sin_fit == "OS") bf.Fix_par("Al2_tm", 0.0);
		  if(t_FSEs != "OS" || t_sin_fit == "tm") bf.Fix_par("Al2_OS", 0.0);
		  if( (t_a2 != "on" && t_a2 != "tm") || t_sin_fit == "OS") bf.Fix_par("D_tm", 0.0);
		  if(  (t_a2 != "on" && t_a2 != "OS") || t_sin_fit == "tm") bf.Fix_par("D_OS", 0.0);
		  if( (t_a4 != "on" && t_a4 != "tm") || t_sin_fit == "OS") bf.Fix_par("D4_tm", 0.0);
		  if( (t_a4 != "on" && t_a4 != "OS") || t_sin_fit == "tm") bf.Fix_par("D4_OS", 0.0);
		  bf.Fix_par("n", n_m_pair.first);
		  bf.Fix_par("m", n_m_pair.second);


		  int dof = Nmeas - npars;
		  Npars.push_back(npars);
		  Ndof.push_back(dof);
  

  
	      
		  bf.ansatz= [=](const W_fpar& X, const W_ipar& Y) {
			       double csi= Y.a==0?csi_phys:pow(Y.Mp/(4.0*M_PI*Y.fp),2); //if a=0 csi is set to its physical value
			       double MpL;
			       double Mp_dim = Y.a==0?Mp_phys:Y.Mp/Y.a; //if a=0 Mp is evaluated at the physical point
			       double n_val= Y.Is_tm?X.n:X.m;
			       MpL= Y.a==0?40.0:Y.Mp*Y.L; //if a=0 MpL is set to 40 (which ~ corresponds to thermodynamic limit)
			       double art;
			       if(n_val >=0) art = Y.a==0?0:pow(Y.a,2)/pow(log(w0_scale/Y.a),n_val);
			       else if(n_val == -1) art = Y.a==0?0:pow(Y.a,2)*log( fabs( log( pow(Y.a,2))));
			       else crash("log-parameter n: "+to_string(n_val)+" is not allowed");
			       double par_art = Y.Is_tm?X.D_tm:X.D_OS;
			       double par_art_a4 = Y.Is_tm?X.D4_tm:X.D4_OS;
			       double par_vol_art = Y.Is_tm?X.Al2_tm:X.Al2_OS;
			       return X.w0*(1.0 + X.Am*(Mp_dim-Mp_phys) + X.Plog*csi*log(csi/csi_phys)+ par_art*art + par_art_a4*pow(Y.a,4))*(1.0  + (X.Al1 +par_vol_art*pow(Y.a,2))*csi*(1.0/pow(MpL,1.5))*exp(-MpL)); // pion mass Y.Mp is in [GeV]
			     };
		  bf.measurement = [=](const W_fpar& X, const W_ipar& Y) {
				     return Y.w_val;
				   };
		  bf.error = [=](const W_fpar& X, const W_ipar& Y) {
			       return Y.w_err;
			     };

		  vector<vector<W_ipar>> ipar_all_ens(UseJack?Njacks:Nboots);
		  for(auto &ipar_jack: ipar_all_ens) ipar_jack.resize(Tot_Nmeas);

		  int id_obs=0;
	      
		  for(int iens=0; iens<Nens;iens++) {  
		    for(int ijack=0;ijack<(UseJack?Njacks:Nboots);ijack++) {
		      //tm_data
		      if(t_sin_fit != "OS" && (which_lat != "three_finest" || Ens_Tag[iens].substr(1,1) != "A") ) {
			if(Ens_Tag[iens].substr(1,1) == "A") {ipar_all_ens[ijack][id_obs].ibeta=0; ipar_all_ens[ijack][id_obs].a= a_A.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "B") {ipar_all_ens[ijack][id_obs].ibeta=1; ipar_all_ens[ijack][id_obs].a= a_B.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "C") {ipar_all_ens[ijack][id_obs].ibeta=2; ipar_all_ens[ijack][id_obs].a = a_C.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "D") {ipar_all_ens[ijack][id_obs].ibeta=3; ipar_all_ens[ijack][id_obs].a = a_D.distr[ijack];}
			else crash("In windows cont+therm+phys lim ensemble tag not listed");
			ipar_all_ens[ijack][id_obs].Mp = Mpi_CDH.distr_list[iens].distr[ijack]; 
			ipar_all_ens[ijack][id_obs].fp = fp_CDH.distr_list[iens].distr[ijack];
			ipar_all_ens[ijack][id_obs].L = L_list[iens];
			ipar_all_ens[ijack][id_obs].Is_tm = true;
			ipar_all_ens[ijack][id_obs].w_val = 1.0e10*meas_tm.distr_list[iens].distr[ijack];
			ipar_all_ens[ijack][id_obs].w_err = 1.0e10*meas_tm.err(iens);
		   
		      }
		    }
		    if(t_sin_fit != "OS" && (which_lat != "three_finest" || Ens_Tag[iens].substr(1,1) != "A")) id_obs++;
		  }

		  for(int iens=0; iens<Nens;iens++) {
		    for(int ijack=0;ijack<(UseJack?Njacks:Nboots);ijack++) {
		      //OS data
		      if(t_sin_fit != "tm" && (which_lat != "three_finest" || Ens_Tag[iens].substr(1,1) != "A") ) {
			if(Ens_Tag[iens].substr(1,1) == "A") {ipar_all_ens[ijack][id_obs].ibeta=0; ipar_all_ens[ijack][id_obs].a= a_A.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "B") {ipar_all_ens[ijack][id_obs].ibeta=1; ipar_all_ens[ijack][id_obs].a= a_B.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "C") {ipar_all_ens[ijack][id_obs].ibeta=2; ipar_all_ens[ijack][id_obs].a = a_C.distr[ijack];}
			else if(Ens_Tag[iens].substr(1,1) == "D") {ipar_all_ens[ijack][id_obs].ibeta=3; ipar_all_ens[ijack][id_obs].a = a_D.distr[ijack];}
			else crash("In windows cont+therm+phys lim ensemble tag not listed");
			ipar_all_ens[ijack][id_obs].Mp = Mpi_CDH.distr_list[iens].distr[ijack];
			ipar_all_ens[ijack][id_obs].fp = fp_CDH.distr_list[iens].distr[ijack];
			ipar_all_ens[ijack][id_obs].L = L_list[iens];
			ipar_all_ens[ijack][id_obs].Is_tm = false;
			ipar_all_ens[ijack][id_obs].w_val = 1.0e10*meas_OS.distr_list[iens].distr[ijack];
			ipar_all_ens[ijack][id_obs].w_err = 1.0e10*meas_OS.err(iens);
		  
		      }
		    }
		    if(t_sin_fit != "tm" && (which_lat != "three_finest" || Ens_Tag[iens].substr(1,1) != "A")) id_obs++;
		  }


		  //append
		  bf.Append_to_input_par(ipar_all_ens);
		  //fit
		  boot_fit_data<W_fpar> Bt_w_fit = bf.Perform_bootstrap_fit();

  
		  //retrieve parameter
		  distr_t w0, Am, Al1, Al2_tm, Al2_OS, D_tm, D_OS, Plog, n, m , D4_tm, D4_OS;

  
 
		  for(int ijack=0;ijack<(UseJack?Njacks:Nboots);ijack++) {
		    W_fpar my_w0_fit_pars = Bt_w_fit.par[ijack];
		    w0.distr.push_back(my_w0_fit_pars.w0);
		    Am.distr.push_back(my_w0_fit_pars.Am);
		    Al1.distr.push_back(my_w0_fit_pars.Al1);
		    Al2_tm.distr.push_back(my_w0_fit_pars.Al2_tm);
		    Al2_OS.distr.push_back(my_w0_fit_pars.Al2_OS);
		    D_tm.distr.push_back(my_w0_fit_pars.D_tm);
		    D_OS.distr.push_back(my_w0_fit_pars.D_OS);
		    Plog.distr.push_back(my_w0_fit_pars.Plog);
		    n.distr.push_back(my_w0_fit_pars.n);
		    m.distr.push_back(my_w0_fit_pars.m);
		    D4_tm.distr.push_back(my_w0_fit_pars.D4_tm);
		    D4_OS.distr.push_back(my_w0_fit_pars.D4_OS);
		  }
  
		  //print info
		  string only_fin_tag= (which_lat=="three_finest")?"on":"off";
		  if(print_par_info) {
		    cout<<"########  fit parameter for Obs: "<<W_type<<" channel: "<<channel<<" ##########"<<endl;
		    cout<<"single fit: "<<t_sin_fit<<endl;
		    cout<<"only_finest: "<<only_fin_tag<<endl;
		    cout<<"w0: ("<<w0.ave()<<" +- "<<w0.err()<<") x e-10"<<endl;
		    cout<<"Am: "<<Am.ave()<<" +- "<<Am.err()<<endl;
		    cout<<"Al1: "<<Al1.ave()<<" +- "<<Al1.err()<<endl;
		    cout<<"Al2 (tm): "<<Al2_tm.ave()<<" +- "<<Al2_tm.err()<<endl;
		    cout<<"Al2 (OS): "<<Al2_OS.ave()<<" +- "<<Al2_OS.err()<<endl;
		    cout<<"a^2 (tm): "<<D_tm.ave()<<" +- "<<D_tm.err()<<endl;
		    cout<<"a^2 (OS): "<<D_OS.ave()<<" +- "<<D_OS.err()<<endl;
		    cout<<"Plog: "<<Plog.ave()<<" +- "<<Plog.err()<<endl;
		    cout<<"a^4 (tm): "<<D4_tm.ave()<<" +- "<<D4_tm.err()<<endl;
		    cout<<"a^4 (OS): "<<D4_OS.ave()<<" +- "<<D4_OS.err()<<endl;
		    cout<<"n: "<<n.ave()<<endl;
		    cout<<"m: "<<m.ave()<<endl;
		    cout<<"#######################################################"<<endl;
		  }
  

		  //forward info on chi2, expectation values and Aka_weight

		  double small_sample_corr = (2.0*pow(npars,2) +2.0*npars)/(1.0*Nmeas-npars-1.0);
		  double ch2  = Bt_w_fit.get_ch2_ave();
		  double ak_weight= exp( -0.5*(ch2 + 2.0*npars - Nmeas));
		  double ak_weight_corr= exp(-0.5*(ch2 + 2.0*npars - Nmeas + small_sample_corr));
		  Ch2.push_back( ch2);
		  Aka_weight.push_back(ak_weight);
		  Aka_weight_red.push_back(ak_weight_corr);
 

   
		  //forward info on expectation value and error of cont extrapolated val
		  val.push_back(w0.ave());
		  err.push_back(w0.err());
   


		  //#################################     PLOT FITTING FUNCTION   ###############################################

		  double csi_points=1;
		  double mL_points=1;
		  double alat_points=40;
		  double tm_OS_points=2;
		  Vfloat csi, mL, alat, tm_OS;
		  Vfloat func_W_val, func_W_err;

		  for(int ir=0; ir<tm_OS_points;ir++) 
		    for(int i_csi=0;i_csi<csi_points;i_csi++)
		      for(int i_l=0;i_l<mL_points; i_l++)
			for(int i_a=0;i_a<alat_points;i_a++) {
			  double csi_i= csi_phys;
			  double li = 40.0;
			  double al;
			  al = i_a*2.0*a_B.ave()/alat_points;
			  csi.push_back(csi_i);
			  mL.push_back(li);
			  alat.push_back(al);
			  tm_OS.push_back(ir);
			  distr_t func_W(UseJack);
			  W_ipar Yi;
			  Yi.Mp = al*Mp_phys;
			  Yi.L = li/Yi.Mp;
			  Yi.a = al;
			  Yi.fp = al*fp_phys;
			  Yi.Is_tm = ir;
			  for(int ijack=0;ijack<(UseJack?Njacks:Nboots);ijack++) {
			    W_fpar Xi;
			    Xi.w0 = w0.distr[ijack];
			    Xi.Am = Am.distr[ijack];
			    Xi.Al1 = Al1.distr[ijack];
			    Xi.Al2_tm = Al2_tm.distr[ijack];
			    Xi.Al2_OS = Al2_OS.distr[ijack];
			    Xi.D_tm = D_tm.distr[ijack];
			    Xi.D_OS = D_OS.distr[ijack];
			    Xi.Plog = Plog.distr[ijack];
			    Xi.n = n.distr[ijack];
			    Xi.m = m.distr[ijack];
			    Xi.D4_tm = D4_tm.distr[ijack];
			    Xi.D4_OS = D4_OS.distr[ijack];
			    func_W.distr.push_back( bf.ansatz(Xi, Yi));
			  }
			  func_W_val.push_back( func_W.ave());
			  func_W_err.push_back( func_W.err());
			}



		  //###########################     END FITTING FUNCTION PLOT  #####################################################
 

    

  

		  //##########################   correct lattice data for FSEs and physical point extrapolation  ###################

		  distr_t_list data_extr_TM(UseJack), data_extr_OS(UseJack);

	     
		  for(int iens=0; iens<Nens;iens++) {

		    distr_t csi = Csi_CDH.distr_list[iens];
		    distr_t MpiL = Mpi_CDH.distr_list[iens]*L_list[iens];
		    distr_t a_iens = a_distr_list.distr_list[iens];
		    distr_t Mpi_dim = Mpi_CDH.distr_list[iens]/a_distr_list.distr_list[iens];
    
		    if(t_sin_fit != "OS") data_extr_TM.distr_list.push_back((meas_tm.distr_list[iens]/(1.0 + (Al1 + Al2_tm*a_iens*a_iens)*csi*distr_t::f_of_distr(FVE, MpiL))-1.0e-10*w0*Am*(Mpi_dim-Mp_phys)));
		    if(t_sin_fit != "tm") data_extr_OS.distr_list.push_back((meas_OS.distr_list[iens]/(1.0 + (Al1 + Al2_OS*a_iens*a_iens)*csi*distr_t::f_of_distr(FVE, MpiL))-1.0e-10*w0*Am*(Mpi_dim-Mp_phys)));
		  }

		  //#################################################################################################################

		  //Print fitting function
		  boost::filesystem::create_directory("../data/gm2/"+channel);
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_fit_func");
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_fit_func/"+W_type);
		  string only_finest_tag = (which_lat=="three_finest")?"on":"off";
		  string Fit_tag= "a2_"+t_a2+"_FSEs_"+t_FSEs+"_a4_"+t_a4+"_Mp_"+t_mass+"_sfit_"+t_sin_fit+"_only_finest_"+only_finest_tag+"_n_"+to_string(n_m_pair.first)+"_m_"+to_string(n_m_pair.second);
		  string Fit_info = "Ch2: "+to_string_with_precision(ch2,3)+" dof: "+to_string(dof)+" Nmeas: "+to_string(Nmeas);

		  Print_To_File({}, {csi,mL,alat, tm_OS, func_W_val, func_W_err}, "../data/gm2/"+channel+"/windows_fit_func/"+W_type+"/"+Fit_tag+".dat","", "#  xi   Mpi*L a(GeV^-1) tm/OS  val  err "+Fit_info);

		  //Print data corrected for FSEs and physical point extrapolation
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_extr_data");
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_extr_data/tm");
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_extr_data/OS");
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_extr_data/tm/"+W_type);
		  boost::filesystem::create_directory("../data/gm2/"+channel+"/windows_extr_data/OS/"+W_type);
  
  
		  if(t_sin_fit != "OS") {
		    Print_To_File(Ens_Tag, {L_list, a_distr_list.ave(), data_extr_TM.ave(), data_extr_TM.err()    }, "../data/gm2/"+channel+"/windows_extr_data/tm/"+W_type+"/"+Fit_tag+".dat", "", "#Ens L a val "+Fit_info);
		  }
		  if(t_sin_fit != "tm") {
		    Print_To_File(Ens_Tag, {L_list, a_distr_list.ave(), data_extr_OS.ave(), data_extr_OS.err()    }, "../data/gm2/"+channel+"/windows_extr_data/OS/"+W_type+"/"+Fit_tag+".dat", "", "#Ens L a val "+Fit_info);
		  }

		  Nfits++;
		}
	      }
	    }
  //########################## End correct lattice data for FSEs and physical point extrapolation  ###################





  
  //#######################################  COMPUTATION OF FINAL CENTRAL VALUES AND ERRORS #########################



  //############################    PRINT SUMMARY TABLES ##############################################

  boost::filesystem::create_directory("../data/gm2/"+channel+"/fit_summary_tables");
  

  //total Akaike weight
  double total_aka_weight =0.0;
  for(auto ak:Aka_weight) total_aka_weight += ak;
  
  //total Akaike weight reduced
  double total_aka_weight_red =0.0;
  for(auto ak_red:Aka_weight_red) total_aka_weight_red += ak_red;


  ofstream Print_table("../data/gm2/"+channel+"/fit_summary_tables/"+W_type+"_Nfits_"+to_string(Nfits)+".tab");
  Print_table<<"FSEs a2  a4   Mass   single_fit only_finest  n  m  Ch2  Npars  Ndof  Akaike   Akaike_reduced"<<endl;
  for(int ifit=0;ifit<Nfits;ifit++)
    Print_table<<FSEs_list[ifit]<<"\t"<<a2_list[ifit]<<"\t"<<a4_list[ifit]<<"\t"<<mass_extr_list[ifit]<<"\t"<<sin_fit_list[ifit]<<"\t"<<only_finest_list[ifit]<<"\t"<<n_m_list[ifit].first<<"\t"<<n_m_list[ifit].second<<"\t"<<Ch2[ifit]<<"\t"<<Npars[ifit]<<"\t"<<Ndof[ifit]<<"\t"<<Aka_weight[ifit]/total_aka_weight<<"\t"<<Aka_weight_red[ifit]/total_aka_weight_red<<endl;

  Print_table.close();
	    

  //###################################### END PRINTING SUMMARY TABLES ###########################


  //##################################### FINAL AVERAGE AND STD ERR ##############################


  //model 1  (standard Akaike)
  double fin_W_val(0.0), fin_W_stat(0.0), fin_W_sist(0.0);
 
  //model 2  (Akaike corrected for small sample)
  double fin_W_val_red(0.0), fin_W_stat_red(0.0), fin_W_sist_red(0.0);

  //model 3  (flat weight)
  double fin_W_val_flat(0.0), fin_W_stat_flat(0.0), fin_W_sist_flat(0.0);

  
  //compute <> and stds for model 1
  //#############################################################################################
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_val += (Aka_weight[ifit]/total_aka_weight)*val[ifit];
    fin_W_stat += (Aka_weight[ifit]/total_aka_weight)*pow(err[ifit],2);
  }
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_sist +=  (Aka_weight[ifit]/total_aka_weight)*pow(fin_W_val-val[ifit],2);
  }
  //#############################################################################################

    
  //compute <> and stds for model 2
  //#############################################################################################
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_val_red += (Aka_weight_red[ifit]/total_aka_weight_red)*val[ifit];
    fin_W_stat_red += (Aka_weight_red[ifit]/total_aka_weight_red)*pow(err[ifit],2);
  }
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_sist_red +=  (Aka_weight_red[ifit]/total_aka_weight_red)*pow(fin_W_val_red-val[ifit],2);
  }


  //compute <> and stds for model 3
  //#############################################################################################
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_val_flat += (1.0/Nfits)*val[ifit];
    fin_W_stat_flat += (1.0/Nfits)*pow(err[ifit],2);
  }
  for(int ifit=0; ifit<Nfits;ifit++) {
    fin_W_sist_flat +=  (1.0/Nfits)*pow(fin_W_val_flat-val[ifit],2);
  }
  //#############################################################################################


  //find smallest ch2, smallest reduced ch2, and largest Akaike weight (not reduced)
  int largest_Aka_id=0;
  int smallest_ch2_id=0;
  int smallest_ch2_red_id=0;
  double largest_Aka, smallest_ch2, smallest_ch2_red;

  for(int ifit=0;ifit<Nfits;ifit++) {
    if(ifit == 0 ) {
      largest_Aka = Aka_weight[ifit];
      smallest_ch2 = Ch2[ifit];
      smallest_ch2_red = Ch2[ifit]/Ndof[ifit];
    }
    else {
      if( Aka_weight[ifit] > largest_Aka) { largest_Aka = Aka_weight[ifit]; largest_Aka_id = ifit;}
      if( Ch2[ifit] < smallest_ch2) { smallest_ch2 = Ch2[ifit]; smallest_ch2_id = ifit;}
      if( Ch2[ifit]/Ndof[ifit] < smallest_ch2_red) { smallest_ch2_red = Ch2[ifit]/Ndof[ifit];  smallest_ch2_red_id = ifit;}

    }

  }


  //print the result

  cout<<"################# FINAL ERROR BUDGET ###################"<<endl;
  cout<<"Obs: "<<W_type<<", channel: "<<channel<<endl;
  cout<<"Performed "<<Nfits<<" total fits"<<endl;
  cout<<"(M1): ("<<fin_W_val<<" +- ["<<sqrt(fin_W_stat)<<"]_stat["<<sqrt(fin_W_sist)<<"]_sist)xe-10 = ("<<fin_W_val<<" +- "<<sqrt(fin_W_stat + fin_W_sist)<<")xe-10"<<endl;
  cout<<"(M2): ("<<fin_W_val_red<<" +- ["<<sqrt(fin_W_stat_red)<<"]_stat["<<sqrt(fin_W_sist_red)<<"]_sist)xe-10 = ("<<fin_W_val_red<<" +- "<<sqrt(fin_W_stat_red + fin_W_sist_red)<<")xe-10"<<endl;
  cout<<"(M3): ("<<fin_W_val_flat<<" +- ["<<sqrt(fin_W_stat_flat)<<"]_stat["<<sqrt(fin_W_sist_flat)<<"]_sist)xe-10 = ("<<fin_W_val_flat<<" +- "<<sqrt(fin_W_stat_flat + fin_W_sist_flat)<<")xe-10"<<endl;

  
  cout<<"# LARGEST AKAIKE WEIGHT OBTAINED: "<<largest_Aka/total_aka_weight<<endl;
  cout<<"# CORRESPONDING TO: "<<endl;
  cout<<"single fit: "<<sin_fit_list[largest_Aka_id]<<endl;
  cout<<"only finest: "<<only_finest_list[largest_Aka_id]<<endl;
  cout<<"FSEs mode: "<<FSEs_list[largest_Aka_id]<<endl;
  cout<<"a2 mode: "<<a2_list[largest_Aka_id]<<endl;
  cout<<"a4 mode: "<<a4_list[largest_Aka_id]<<endl;
  cout<<"mass extrapolation mode: "<<mass_extr_list[largest_Aka_id]<<endl;
  cout<<"(n,m): ("<<n_m_list[largest_Aka_id].first<<","<<n_m_list[largest_Aka_id].second<<")"<<endl;
  cout<<"Ndof: "<<Ndof[largest_Aka_id]<<endl;
  cout<<"Npars: "<<Npars[largest_Aka_id]<<endl;
  cout<<"Result: "<<val[largest_Aka_id]<<" +- "<<err[largest_Aka_id]<<endl;
  cout<<"# SMALLEST ch2 OBTAINED: "<<smallest_ch2<<endl;
  cout<<"# CORRESPONDING TO: "<<endl;
  cout<<"single fit: "<<sin_fit_list[smallest_ch2_id]<<endl;
  cout<<"only finest: "<<only_finest_list[smallest_ch2_id]<<endl;
  cout<<"FSEs mode: "<<FSEs_list[smallest_ch2_id]<<endl;
  cout<<"a2 mode: "<<a2_list[smallest_ch2_id]<<endl;
  cout<<"a4 mode: "<<a4_list[smallest_ch2_id]<<endl;
  cout<<"mass extrapolation mode: "<<mass_extr_list[smallest_ch2_id]<<endl;
  cout<<"(n,m): ("<<n_m_list[smallest_ch2_id].first<<","<<n_m_list[smallest_ch2_id].second<<")"<<endl;
  cout<<"Ndof: "<<Ndof[smallest_ch2_id]<<endl;
  cout<<"Npars: "<<Npars[smallest_ch2_id]<<endl;
  cout<<"Result: "<<val[smallest_ch2_id]<<" +- "<<err[smallest_ch2_id]<<endl;
  cout<<"# SMALLEST REDUCED ch2 OBTAINED: "<<smallest_ch2_red<<endl;
  cout<<"# CORRESPONDING TO: "<<endl;
  cout<<"single fit: "<<sin_fit_list[smallest_ch2_red_id]<<endl;
  cout<<"only finest: "<<only_finest_list[smallest_ch2_red_id]<<endl;
  cout<<"FSEs mode: "<<FSEs_list[smallest_ch2_red_id]<<endl;
  cout<<"a2 mode: "<<a2_list[smallest_ch2_red_id]<<endl;
  cout<<"a4 mode: "<<a4_list[smallest_ch2_red_id]<<endl;
  cout<<"mass extrapolation mode: "<<mass_extr_list[smallest_ch2_red_id]<<endl;
  cout<<"(n,m): ("<<n_m_list[smallest_ch2_red_id].first<<","<<n_m_list[smallest_ch2_red_id].second<<")"<<endl;
  cout<<"Ndof: "<<Ndof[smallest_ch2_red_id]<<endl;
  cout<<"Npars: "<<Npars[smallest_ch2_red_id]<<endl;
  cout<<"Result: "<<val[smallest_ch2_red_id]<<" +- "<<err[smallest_ch2_red_id]<<endl;
  cout<<"####################### Bye! ###########################"<<endl;
  
  return;  
}
