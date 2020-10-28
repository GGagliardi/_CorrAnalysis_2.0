#include "../include/PionMassAnalysis.h"


using namespace std;
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.139;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const double r0 = pow(0.672/0.197,2);
const int Nbranches = 1;
const int Nboots= 400;
const bool Use_JB_distribution= true;
const bool UseJack=1;
const int Njacks=15;
const int nboots=100;
const bool verbose=1;



class Y {

public:
  Y() {}
  double L, Mpi, Dm2_exch, Dm2_hand, Mpi_err, Dm2_exch_err, Dm2_hand_err, Cov_exch_hand;
  int ibeta;
};



class X {

public:
  X(const Vfloat &par) : ainv(3), Zv(3), Za(3) {
    if(par.size() != 17) {
      cout<<"In class X, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    this->chir=par[0];
    this->log=par[1];
    this->A_1=par[2];
    this->D=par[3];
    this->Dm=par[4];
    this->F_pi=par[5];
    this->log_a=par[6];
    this->f0=par[7];
    for(int ibeta=0;ibeta<3;ibeta++) {
      this->ainv[ibeta] = par[8+ibeta];
    }
    for(int ibeta=0; ibeta<3;ibeta++) {
      this->Zv[ibeta] = par[11+ibeta];
    }
    for(int ibeta=0; ibeta<3;ibeta++) {
      this->Za[ibeta] = par[14+ibeta];
    }
  }
  X() : ainv(3), Zv(3), Za(3) {}
  vector<double> ainv;
  vector<double> Zv;
  vector<double> Za;
  double f0;
  double D, Dm, chir, log, log_a, A_1, F_pi;
};





void Pion_mass_analysis(string CURRENT_TYPE, bool IncludeDisconnected) {
  
  data_t m_data, dm_exch_data, dm_hand_data, m_data_hand_run;
  if(CURRENT_TYPE=="CONSERVED")  {
    m_data.Read("../datasets", "mes_contr_00", "P5P5");
    dm_exch_data.Read("../datasets", "mes_contr_LL", "P5P5");
  }
  else if (CURRENT_TYPE=="LOCAL") { //current is local
    m_data.Read("../datasetslocal", "mes_contr_M0_R0_0_M0_R0_0", "P5P5");
    dm_exch_data.Read("../datasetslocal", "mes_contr_M0_R0_F_M0_R0_F", "P5P5");
    if(IncludeDisconnected) {
      m_data_hand_run.Read("../datasetslocal", "mes_contr_M0_R0_0_M0_R0_0_handcuffs", "P5P5");
      dm_hand_data.Read("../datasetslocal", "handcuffs", "P5P5");
    }
  }
  else crash("CURRENT_TYPE: "+CURRENT_TYPE+ " not yet implemented. Exiting...");
  
  
  //compute masses and effective masses

  
  
  int Nens = m_data.size;  //== dm_exch_data.size .....
  cout<<"N_ens: "<<Nens<<endl;



  distr_t_list Mpi_distr_list(UseJack), Dm2_exch_distr_list(UseJack), Dm2_hand_distr_list(UseJack);
  vector<Eigen::MatrixXd> CovMatrixPion(Nens);

  for(int i=0; i < Nens; i++) {
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(m_data.Tag[i].substr(0,1) == "A") Corr.Tmin = 13;
    else if(m_data.Tag[i].substr(0,1) =="B") Corr.Tmin = 16;
    else Corr.Tmin= 19;
    Corr.Tmax = m_data.nrows[i]/2 -5;
    Corr.Nt = m_data.nrows[i];
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/Mpi");
    string p= (CURRENT_TYPE=="CONSERVED")?"conserved":"local";
    boost::filesystem::create_directory("../data/Mpi/"+p);
   
    distr_t_list Mpi_eff_distr = Corr.effective_mass_t(m_data.col(0)[i], "../data/Mpi/"+p+"/mass."+m_data.Tag[i]);

    distr_t_list Dm_exch_eff_distr =  Corr.effective_slope_t(dm_exch_data.col(0)[i], m_data.col(0)[i], "../data/Mpi/"+p+"/dm_exch."+m_data.Tag[i]);
   
    distr_t_list Dm_hand_eff_distr, Mpi_eff_hand_run_distr;
    if(IncludeDisconnected) {
      Mpi_eff_hand_run_distr= Corr.effective_mass_t(m_data_hand_run.col(0)[i], "");
      Dm_hand_eff_distr= Corr.effective_slope_t(dm_hand_data.col(0)[i], m_data_hand_run.col(0)[i], "../data/Mpi/"+p+"/dm_hand."+m_data.Tag[i]);
    }

    
    Mpi_distr_list.distr_list.push_back(Corr.Fit_distr( Mpi_eff_distr));
    Dm2_exch_distr_list.distr_list.push_back ( 2.0*Corr.Fit_distr(Dm_exch_eff_distr)*Corr.Fit_distr(Mpi_eff_distr));
    if(IncludeDisconnected) Dm2_hand_distr_list.distr_list.push_back( 2.0*Corr.Fit_distr(Dm_hand_eff_distr)*Corr.Fit_distr(Mpi_eff_hand_run_distr));
    else Dm2_hand_distr_list.distr_list.emplace_back(UseJack, Dm2_exch_distr_list.distr_list[i].size());

    if(!Use_JB_distribution) { //resample from Gaussian distribution
      if(IncludeDisconnected) Compute_covariance_matrix(UseJack, CovMatrixPion[i], 3, Mpi_distr_list.distr_list[i].distr, Dm2_exch_distr_list.distr_list[i].distr, Dm2_hand_distr_list.distr_list[i].distr);
      else Compute_covariance_matrix(UseJack, CovMatrixPion[i],2, Mpi_distr_list.distr_list[i].distr, Dm2_exch_distr_list.distr_list[i].distr);
      
    }
  }


  
  //everything is set up. Start fitting!
  LatticeInfo L_info(CURRENT_TYPE);


  bootstrap_fit<X,Y> bf(nboots);


  bf.Set_number_of_measurements(Nens);
  bf.Set_verbosity(verbose);
  

  //Add parameters
  bf.Add_par("chir", 4e-5, 1e-6);
  bf.Add_par("log", 3.7 , 0.7);
  bf.Add_par("A_1", -5.7, 0.02);
  bf.Add_par("D", 2e-3, 1e-5);
  bf.Add_par("Dm", 0.5, 1e-3);
  bf.Add_par("F_pi", 2.5, 0.1);
  bf.Add_par("log_a", 0.1, 1e-3);
  bf.Add_prior_par("f0", 0.121, 1e-3);
  bf.Add_prior_pars({"ainv0", "ainv1", "ainv2", "Zv0", "Zv1", "Zv2", "Za0", "Za1", "Za2"});
    
  //Fix some parameters to make test
  //bf.Fix_par("Dm", 0.0);
  bf.Fix_par("log_a", 0.0);

  //Add List of parameters to be released after first minimization
  bf.Fix_n_release({"ainv0", "ainv1", "ainv2", "f0"});

  if(CURRENT_TYPE=="CONSERVED") {
    bf.Fix_par("Zv0", 1.0);
    bf.Fix_par("Zv1",1.0);
    bf.Fix_par("Zv2",1.0 );
    bf.Fix_par("Za0",0.0);
    bf.Fix_par("Za1",0.0);
    bf.Fix_par("Za2",0.0);
    
  }
  else if(CURRENT_TYPE=="LOCAL") {
    bf.Fix_n_release("Zv0");
    bf.Fix_n_release("Zv1");
    bf.Fix_n_release("Zv2");
    if(IncludeDisconnected) {
      bf.Fix_n_release("Za0");
      bf.Fix_n_release("Za1");
      bf.Fix_n_release("Za2");
    }
    else {  bf.Fix_par("Za0",0.0);  bf.Fix_par("Za1",0.0);  bf.Fix_par("Za2",0.0);}
  }
  
  else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");


  
 
  Vfloat L, m_l, T;
  Read_pars_from_ensemble_tag(m_data.Tag, m_l, L, T);
  


  //define fitting function

  //############################################################################################################

  bf.ansatz =  [=](const X &p, const Y &ip) -> double {


    double ainv = p.ainv[ip.ibeta];
    double L = ip.L;

     
    double Mp = ainv*ip.Mpi;
    double Mp2 = pow(Mp,2);

     
    double SD_FVE = (L>0.0)?(e2/3.0)*ip.Mpi*(1.0/pow(L/ainv,3))*r0*(1 + p.F_pi/(pow(ainv,2))):0.0;

     
    double FVE_universal = (L>0.0)?(kappa*alpha/L)*( ip.Mpi + 2.0/L):0.0;

     
    double fitted_value = 16*M_PI*alpha*p.chir*(1.0/pow(p.f0,2))
    -(1.0/(4*M_PI))*alpha*(p.log +p.log_a/pow(ainv,2))*Mp2*log( Mp2/(pow(4*M_PI*p.f0,2)))
    + alpha*Mp2*(1.0/(4.0*M_PI))*p.A_1
    + p.D/pow(ainv,2)
    + SD_FVE - FVE_universal*pow(ainv,2)
    + p.Dm*Mp2/pow(ainv,2);


    
     
    fitted_value /= pow(ainv,2);
    

    return fitted_value;
  };
  

  bf.measurement = [=](const X& p,const Y& ip) -> double {

    double Zv = p.Zv[ip.ibeta];
    double Za = p.Za[ip.ibeta];
    
        
    double m1 = (e2/2)*pow(Zv,2)*ip.Dm2_exch;
    double m2 = (e2/2)*pow(Za,2)*ip.Dm2_hand;
    
    return m1-m2;
  };
  

  bf.error =  [=](const X& p,const  Y &ip) -> double {

    double Zv = p.Zv[ip.ibeta];
    double Za = p.Za[ip.ibeta];

       
    double m1 = (e2/2)*pow(Zv,2)*ip.Dm2_exch_err;
    double m2 = (e2/2)*pow(Za,2)*ip.Dm2_hand_err;
    double m3 = pow( e2/2,2)*pow(Za*Zv,2);
    return sqrt( pow(m1,2) +
		 pow(m2,2)
		 -2.0*m3*ip.Cov_exch_hand);
  };

 
  
  
  //############################################################################################################

  //init random Number Generators

  GaussianMersenne GM(54353);
  RandomMersenne RM(23423, UseJack?(Njacks-1):(Nboots-1));

  Eigen::MatrixXd CovMatrixInput(9,9); //covariance matrix of input parameters
  Eigen::VectorXd Ave_input_parameters(9);

  vector<boot_fit_data<X>> Bt_fit;


  //we store here all resampled data. 
  vector<vector<vector<Y>>> all_ens(Nbranches);
  for(auto &all_ens_br : all_ens) {
    all_ens_br.resize(nboots);
    for(auto & all_ens_ib: all_ens_br) {
      all_ens_ib.resize(Nens);
      for(int iens=0; iens<Nens;iens++) {
	all_ens_ib[iens].L=L[iens];
	if(m_data.Tag[iens].substr(0,1) == "A") all_ens_ib[iens].ibeta=0;
	else if(m_data.Tag[iens].substr(0,1) == "B") all_ens_ib[iens].ibeta=1;
	else all_ens_ib[iens].ibeta=2;
	all_ens_ib[iens].Mpi_err = Mpi_distr_list.err(iens);
	all_ens_ib[iens].Dm2_hand_err = Dm2_hand_distr_list.err(iens);
	all_ens_ib[iens].Dm2_exch_err = Dm2_exch_distr_list.err(iens);
	all_ens_ib[iens].Cov_exch_hand = Dm2_exch_distr_list[iens]%Dm2_hand_distr_list[iens]; 
      }
    }
  }


  for(int ibranch=0; ibranch < Nbranches; ibranch++) {


    
    ReadBranch(ibranch, CovMatrixInput, Ave_input_parameters );


    //clear all measurements
    bf.Clear_priors();
    bf.Clear_input_pars();
    
    

    //generate bootstrap data
    for(int iboot=0; iboot<nboots; iboot++) {
      bf.ib= &iboot;
      bool IsLoc = CURRENT_TYPE != "CONSERVED";
      Vfloat lat_input = Covariate(CovMatrixInput, Ave_input_parameters, GM);
      double Zv0= gauss(L_info.Retrieve_Zv("A",ibranch) , GM);
      double Zv1= gauss(L_info.Retrieve_Zv("B", ibranch), GM);
      double Zv2= gauss(L_info.Retrieve_Zv("D", ibranch), GM);
      double Za0= gauss(L_info.Retrieve_Za("A", ibranch),GM);
      double Za1= gauss(L_info.Retrieve_Za("B", ibranch),GM);
      double Za2= gauss(L_info.Retrieve_Za("D", ibranch) ,GM);
      bf.Append_to_prior("ainv0", lat_input[3], sqrt(CovMatrixInput(3,3)));
      bf.Append_to_prior("ainv1", lat_input[4], sqrt(CovMatrixInput(4,4)));
      bf.Append_to_prior("ainv2", lat_input[5], sqrt(CovMatrixInput(5,5)));
      bf.Append_to_prior("f0", lat_input[2], sqrt(CovMatrixInput(2,2)));
      bf.Append_to_prior("Zv0", IsLoc?Zv0:1.0, L_info.Retrieve_Zv("A", ibranch).second);
      bf.Append_to_prior("Zv1", IsLoc?Zv1:1.0, L_info.Retrieve_Zv("B", ibranch).second);
      bf.Append_to_prior("Zv2", IsLoc?Zv2:1.0, L_info.Retrieve_Zv("D", ibranch).second);
      bf.Append_to_prior("Za0", IncludeDisconnected?Za0:0.0, L_info.Retrieve_Za("A", ibranch).second);
      bf.Append_to_prior("Za1", IncludeDisconnected?Za1:0.0, L_info.Retrieve_Za("B", ibranch).second);
      bf.Append_to_prior("Za2", IncludeDisconnected?Za2:0.0, L_info.Retrieve_Za("D", ibranch).second);


      int k= RM();
      for(int imeas=0; imeas <Nens; imeas++) {
	if(!Use_JB_distribution) {
	  Eigen::VectorXd vec(CovMatrixPion[imeas].rows());
	  if(IncludeDisconnected)  vec<<Mpi_distr_list.ave(imeas),Dm2_exch_distr_list.ave(imeas),Dm2_hand_distr_list.ave(imeas);
	  else vec<<Mpi_distr_list.ave(imeas),Dm2_exch_distr_list.ave(imeas);
	  Vfloat res_meas(3,0.0);
	 
	  res_meas= Covariate(CovMatrixPion[imeas],vec, GM);
	  all_ens[ibranch][iboot][imeas].Mpi= res_meas[0];
	  all_ens[ibranch][iboot][imeas].Dm2_hand = IncludeDisconnected?res_meas[2]:0.0;
	  all_ens[ibranch][iboot][imeas].Dm2_exch = res_meas[1];
	 	 
	}
      
	else {
	  all_ens[ibranch][iboot][imeas].Mpi = Mpi_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].Dm2_exch = Dm2_exch_distr_list.distr_list[imeas].distr[k];
	  all_ens[ibranch][iboot][imeas].Dm2_hand = Dm2_hand_distr_list.distr_list[imeas].distr[k];
	}
      }

    }
    
    bf.Append_to_input_par(all_ens[ibranch]);
    Bt_fit.push_back(bf.Perform_bootstrap_fit());
    
  }
  
 
  //print the data
  
  //define lambda functions

  auto FVE = [](double L, double Mp) {
    return (kappa*alpha/L)*( Mp + 2.0/L);
  };

  auto SD_FVE = [](double L, double Mp, double ainv,  double F_p) {
    double SDE =(e2/3.0)*Mp*(1.0/pow(L,3))*pow( ainv,3)*r0*(1 + F_p*(1.0/pow(ainv,2)));
    return SDE/pow(ainv,2);
  };

  auto cont_ansatz = [](X p, double M) {

    double Mp2 = pow(M,2);
    return  16*M_PI*alpha*p.chir*(1.0/pow(p.f0,2))
    -(1.0/(4*M_PI))*alpha*p.log*Mp2*log( Mp2/(pow(4*M_PI*p.f0,2)))
    + alpha*Mp2*(1.0/(4.0*M_PI))*p.A_1;
  };


  
  
  Vfloat MP(600);
  Vfloat MP_measured(bf.Get_number_of_measurements());
  for(unsigned int m=0; m<MP.size();m++) MP[m] = 0.0010*m;
  VVVfloat SD_subtracted_data, Mp_measured, SD_subtracted_data_dim;
  VVVVfloat Fitted_func_adim, Fitted_func_dim;
  VVfloat Physical_point;
  cascade_resize(SD_subtracted_data, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(SD_subtracted_data_dim, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Mp_measured, Vint{bf.Get_number_of_measurements(), Nbranches, nboots});
  cascade_resize(Fitted_func_adim, Vint{3, (int)MP.size(), Nbranches, nboots});
  cascade_resize(Fitted_func_dim, Vint{4, (int)MP.size(), Nbranches, nboots});
  cascade_resize(Physical_point, Vint{Nbranches, nboots});

  

  
  for(int ibranch=0; ibranch<Nbranches;ibranch++) {
    for(int iboot=0; iboot<nboots;iboot++) {
      X P = Bt_fit[ibranch].par[iboot];
      for(int imeas=0; imeas< bf.Get_number_of_measurements(); imeas ++) {
	Y E = all_ens[ibranch][iboot][imeas];
	SD_subtracted_data[imeas][ibranch][iboot] = bf.measurement(P, E) + FVE(E.L, E.Mpi) - SD_FVE(E.L, E.Mpi, P.ainv[E.ibeta], P.F_pi);
	SD_subtracted_data_dim[imeas][ibranch][iboot] = SD_subtracted_data[imeas][ibranch][iboot]*pow(P.ainv[E.ibeta],2);
	Mp_measured[imeas][ibranch][iboot] = E.Mpi*P.ainv[E.ibeta];
      }
      for(int ibeta=0; ibeta < 4;ibeta++) {
	for(unsigned int m=0; m<MP.size(); m++){
	  if(ibeta<3) {
	    Y point;
	    point.L = -1.0;
	    point.ibeta=ibeta;
	    point.Mpi= MP[m]/P.ainv[ibeta];
	    Fitted_func_adim[ibeta][m][ibranch][iboot] = bf.ansatz(P,point);
	    Fitted_func_dim[ibeta][m][ibranch][iboot] = bf.ansatz(P, point)*pow(P.ainv[ibeta],2);
	  }
	  else Fitted_func_dim[ibeta][m][ibranch][iboot] = cont_ansatz(P, MP[m]) ;
	}
      }
      //Add physical point
      Physical_point[ibranch][iboot] = cont_ansatz(P, MPiPhys);
    }
  }


   

  Vfloat SD_subtracted_val, SD_subtracted_err, MMP, SD_subtracted_dim_val, SD_subtracted_dim_err;
  for(int imeas=0;imeas<bf.Get_number_of_measurements();imeas++) {
    SD_subtracted_val.push_back( Boot_ave(SD_subtracted_data[imeas]));
    SD_subtracted_err.push_back( Boot_err(SD_subtracted_data[imeas]));
    SD_subtracted_dim_val.push_back( Boot_ave( SD_subtracted_data_dim[imeas]));
    SD_subtracted_dim_err.push_back( Boot_err( SD_subtracted_data_dim[imeas]));
    MMP.push_back( Boot_ave(Mp_measured[imeas]));
  }
  
  VVfloat Fitted_func_adim_val, Fitted_func_adim_err, Fitted_func_dim_val, Fitted_func_dim_err;
  cascade_resize(Fitted_func_adim_val, {3,(int)MP.size()});
  cascade_resize(Fitted_func_adim_err, {3,(int)MP.size()});
  cascade_resize(Fitted_func_dim_val, {4, (int)MP.size()});
  cascade_resize(Fitted_func_dim_err, {4, (int)MP.size()});

  for(int ibeta=0;ibeta<4;ibeta++) {
    for(unsigned int m=0; m<MP.size();m++) {
      if(ibeta<3) {
	Fitted_func_adim_val[ibeta][m] = Boot_ave(Fitted_func_adim[ibeta][m]);
	Fitted_func_adim_err[ibeta][m] = Boot_err(Fitted_func_adim[ibeta][m]);
	Fitted_func_dim_val[ibeta][m] = Boot_ave(Fitted_func_dim[ibeta][m]);
	Fitted_func_dim_err[ibeta][m] = Boot_err(Fitted_func_dim[ibeta][m]);
      }
      else {
	Fitted_func_dim_val[ibeta][m] = Boot_ave(Fitted_func_dim[ibeta][m]);
	Fitted_func_dim_err[ibeta][m] = Boot_err(Fitted_func_dim[ibeta][m]);
      }
    }
  }

  //plot the result

  string print_path = (CURRENT_TYPE=="CONSERVED")?"conserved":"local";
  string exch_or_tot = IncludeDisconnected?"tot":"exch";
  // Set the size of output image to 1200x780 pixels
  plt::figure_size(1200*1.2, 780*1.2);
  plt::xlabel("$M_{\pi}$   $[GeV]$");
  plt::ylabel("$\Delta M_{\pi}^{2}$   $[lattice units]$");
  plt::xlim(0.0, MP[MP.size()-1]);
  
  plt::errorbar(MP, Fitted_func_adim_val[0], Fitted_func_adim_err[0], { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "Fit to A"}});
  plt::errorbar(MP, Fitted_func_adim_val[1], Fitted_func_adim_err[1], { {"c", "blue"}, {"marker", "."} , {"ls" , "-"}, {"label", "Fit to B"}});
  plt::errorbar(MP, Fitted_func_adim_val[2], Fitted_func_adim_err[2], { {"c", "green"}, {"marker", "."} , {"ls" , "-"}, {"label", "Fit to D"}});
  plt::errorbar(MMP, SD_subtracted_val, SD_subtracted_err, { {"c", "black"}, {"marker", "o"}, {"ls", ""}, {"label", "Data all FVE subtracted"}});
  // Enable legend.
  plt::legend();
  // Save the image (file format is determined by the extension)
  boost::filesystem::create_directory("../plots");
  boost::filesystem::create_directory("../plots/Mpi");
  string figure_path= "../plots/Mpi/fit_dimless_"+print_path+"_"+exch_or_tot;
  plt::save(figure_path.c_str());
    
  //Print data to file
 
  boost::filesystem::create_directory("../data");
  boost::filesystem::create_directory("../data/Mpi/"+print_path);
  string type = IncludeDisconnected?"total":"exchange";
  for(int ibeta=0;ibeta<3;ibeta++) {
    Print_To_File({}, {MP, Fitted_func_adim_val[ibeta], Fitted_func_adim_err[ibeta]}, "../data/Mpi/"+print_path+"/Fit_function"+type+"_ibeta_"+to_string(ibeta)+".data", "", "#Mpi     #Mpi^2_+ - Mpi^2_0     #err" );
  }
  Print_To_File(m_data.Tag, { MMP, SD_subtracted_val, SD_subtracted_err}, "../data/Mpi/"+print_path+"/sd_subtracted_data_"+type+".dat", "", "#Ens        #Mpi        #Mpi^2_+ - Mpi^2_0      #err");



  plt::clf();
  plt::xlim(0.0, MP[MP.size()-1]);
  plt::xlabel("$M_{\pi}$   $[GeV]$");
  plt::ylabel("$ \Delta M_{pi}^{2}$    $[GeV^{-2}]$");
  plt::errorbar(MP, Fitted_func_dim_val[0], Fitted_func_dim_err[0], { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "A"}});
  plt::errorbar(MP, Fitted_func_dim_val[1], Fitted_func_dim_err[1], { {"c", "blue"}, {"marker", "."} , {"ls" , "-"}, {"label", "B"}});
  plt::errorbar(MP, Fitted_func_dim_val[2], Fitted_func_dim_err[2], { {"c", "green"}, {"marker", "."} , {"ls" , "-"}, {"label", "D"}});
  plt::errorbar(MP, Fitted_func_dim_val[3], Fitted_func_dim_err[3], { {"c", "yellow"}, {"marker", "."} , {"ls" , "-"}, {"label", "cont."}});
  plt::errorbar(MMP, SD_subtracted_dim_val, SD_subtracted_dim_err, { {"c", "black"}, {"marker", "o"}, {"ls", ""}, {"label", "Data, all FVE subtracted"}});
  plt::errorbar(Vfloat{MPiPhys}, Vfloat{1137.0*1e-6}, Vfloat{63.0*1e-6}, {{"c", "red"}, {"marker", "o"}, {"ls", ""}, {"label", "phys from 1707"}});
  plt::plot(MP, Vfloat(MP.size(), 1261.2*1e-6), { {"c","grey"}, {"label", "phys."}});

  plt::legend();

  figure_path= "../plots/Mpi/fit_"+print_path+"_"+exch_or_tot;
  plt::save(figure_path.c_str());


  cout<<"Physical Pion mass difference squared [Mev2]: "<< Boot_ave(Physical_point)*1e+6<<"    "<<Boot_err(Physical_point)*1e+6<<endl;
  cout<<"Average chi2:"<<endl;
  for(int ibr=0;ibr<Nbranches;ibr++) cout<<"Branch: "<<ibr<<" chi2 = "<<accumulate(Bt_fit[ibr].chi2.begin(), Bt_fit[ibr].chi2.end(), 0.0)/Bt_fit[ibr].chi2.size()<<endl;

  
    

    
  return;
}
 
					
  






















