#include "../include/virtual_ff_fit.h"

using namespace std;

const int nsteps= 100;
const double step_size= 1.0/(double)nsteps;
double rk_guess_V = pow(493.677/775.4,2)  ; //Mk^2 /Mrho^2
double rq_guess_V = pow(493.677/892 , 2); //Mk^2/Mk*^2(892)
double rk_guess_A = pow(493.677/775.4,2) ; //Mk^2/Mrho^2
double rq_guess_A = pow(493.677/1253 ,2); //Mk^2/Mk*^2(1270)


class Yff{

public:
  Yff() : Za_ov_Zv(1.0), Za_ov_Zv_err(0.0) {}
  double Mpi, fp, ff, xk, xq, ff_err, Za_ov_Zv, Za_ov_Zv_err;
  
};

class Xff_VMD{

public:
  Xff_VMD() {}
  Xff_VMD(const Vfloat &par) {
    if(par.size() != 5) {
      cout<<"In class Xff_VMD, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    a0=par[0];
    ampl=par[1];
    rk= par[2];
    rq=par[3];
    Za_ov_Zv=par[4];
  }

  double a0, ampl, rk,rq, Za_ov_Zv;
  


};

class Xff_ChPT{
  
public:
  Xff_ChPT() {}
  Xff_ChPT(const Vfloat &par) {
    if(par.size() != 6) {
      cout<<"In class Xff_ChPT, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    a0=par[0];
    ak= par[1];
    aq= par[2];
    akq= par[3];
    a2kq=par[4];
    Za_ov_Zv=par[5];
  }
  double a0, ak, aq, akq, a2kq, Za_ov_Zv;
};

 

void Fit_virtual_FF_VMD( vector<function<double(double, double)>> &fit_func, const vector<distr_t> &FF, distr_t &f_p,distr_t &m_p, distr_t& Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag, string Meson, bool UseJack, bool ConstFit, int t_fit) {

  int Nmeas= FF.size();
  int njacks = f_p.size();
  if(f_p.size() != m_p.size()) crash(" Number of jackknives are not the same for fp and mp");
  bool verbose=1;


  cout<<"Nmeas: "<<Nmeas<<endl;
  cout<<"Njacks: "<<njacks<<endl;
  //perform bootstrap fit

  double L, m_l, T;
  Read_pars_from_ensemble_tag(Ens_tag, m_l, L, T);
  LatticeInfo L_info("LOCAL"); //CURRENT_TYPE == LOCAL


  bootstrap_fit<Xff_VMD,Yff> bf(njacks);


  bf.Set_number_of_measurements(FF.size());
  bf.Set_verbosity(verbose);
  

  //Add parameters
  double a0, a0_err;
  if(ff_type=="H1" || ff_type == "H2") {a0=0.2; a0_err=0.01;} 
  else if(ff_type=="FA")  {a0=0.03; a0_err=0.003;} 
  else if(ff_type=="FV") {a0=0.08*Za_ov_Zv.ave(); a0_err=0.002*Za_ov_Zv.ave();}
  else crash("string : "+ff_type+" does not name a form factor");
  bf.Add_par("a0",a0, a0_err);
  if(W != "H2") bf.Add_par("ampl",a0/100  , a0_err/100);
  else bf.Add_par("ampl",-a0/100  , -a0_err/100);
  if(W=="A") {
  bf.Add_prior_par("rk", rk_guess_A, rk_guess_A/10);
  bf.Add_prior_par("rq", rq_guess_A, rq_guess_A/10);
  }
  else if(W=="V") {
    
    bf.Add_prior_par("rk", rk_guess_V, rk_guess_V/10);
    bf.Add_prior_par("rq", rq_guess_V, rq_guess_V/10);
  }
  else crash("string W: "+W+" is not axial or vector");


    //fix n release
  bf.Add_prior_par("Za_ov_Zv", Za_ov_Zv.ave(), Za_ov_Zv.err());
  if(W != "V") bf.Fix_par("Za_ov_Zv", 1.0);
  bf.Fix_n_release("rk");
  bf.Fix_n_release("rq");
  if(W=="V") bf.Fix_n_release("Za_ov_Zv");

  
  bf.ansatz =  [=](const Xff_VMD &p, const Yff &ip) -> double {
    return  p.a0 + p.ampl/((1.0-pow(ip.xk,2)*p.rk)*(1.0- pow(ip.xq,2)*p.rq));
  };
  bf.measurement = [=](const Xff_VMD& p,const Yff& ip) -> double {
    return p.Za_ov_Zv*ip.ff;
  };
  bf.error =  [=](const Xff_VMD& p,const  Yff &ip) -> double {
    return p.Za_ov_Zv*ip.ff_err;
  };


  //where to store fitted parameters
  boot_fit_data<Xff_VMD> Bt_fit;

  //input parameters
  vector<vector<Yff>> data(njacks);
    for(auto &data_iboot: data) data_iboot.resize(Nmeas);
  
    for(int ijack=0;ijack<njacks;ijack++) {
      for(int imeas=0;imeas<Nmeas;imeas++) {
        data[ijack][imeas].Mpi = m_p.distr[ijack];
	data[ijack][imeas].fp = f_p.distr[ijack];
	data[ijack][imeas].ff = FF[imeas].distr[ijack];
	data[ijack][imeas].ff_err = FF[imeas].err();
	data[ijack][imeas].xk = fabs(mom[imeas].virt()/m_p.distr[ijack]) ;
	data[ijack][imeas].xq = sqrt( 1.0 + pow(data[ijack][imeas].xk,2) - mom[imeas].x_gamma(m_p).distr[ijack]) ;
      } 
    }
   
 


  //add prior 
  for(int ijack=0;ijack<njacks;ijack++) {
     bf.ib= &ijack;
     bf.Append_to_prior("rk", rk_guess_A, rk_guess_A );
     if(W=="A")  { bf.Append_to_prior("rq", rq_guess_A, rq_guess_A ); bf.Append_to_prior("Za_ov_Zv", 1.0, 1.0);}
     else if(W=="V") { bf.Append_to_prior("rq", rq_guess_V, rq_guess_V ) ; bf.Append_to_prior("Za_ov_Zv", Za_ov_Zv.distr[ijack], Za_ov_Zv.err());}
     else crash("string W: "+W+" is not axial nor vector");
  }


  //fit
  bf.Append_to_input_par(data);
  Bt_fit =bf.Perform_bootstrap_fit();


  //define lambda functions
  for(int ijack=0;ijack<njacks;ijack++) {
    auto F = [=](double xk, double xq) -> double { return Bt_fit.par[ijack].a0 + Bt_fit.par[ijack].ampl/((1.0-pow(xk,2)*Bt_fit.par[ijack].rk)*(1.0 -pow(xq,2)*Bt_fit.par[ijack].rq));};
    fit_func.push_back(F);
  }

  //print data
  //set a grid in the xk, xq plane

  Vfloat xk_list, xq_list, FF_list, FF_err_list;
  for(int dxk=0; dxk < nsteps; dxk++) {
    for(int dxq=0; dxq < nsteps; dxq++) {
      double xxk = step_size*dxk;
      double xxq = step_size*dxq;
      distr_t func(UseJack);
      for(int ijack=0;ijack<njacks;ijack++) {
	func.distr.push_back( fit_func[ijack](xxk,xxq));
      }
      xk_list.push_back(xxk);
      xq_list.push_back(xxq);
      FF_list.push_back(func.ave());
      FF_err_list.push_back(func.err());
    }
  }


 
  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag);
  if(ConstFit) {
  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag+"/"+ff_type+".dat", "", "#xk   xq  val    err" );
  }
  else  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag+"/"+ff_type+"time_"+to_string(t_fit)+".dat", "", "#xk   xq  val    err" );
  
  
 

  return;

}







void Fit_virtual_FF_ChPT(vector<function<double(double, double)>> &fit_func, const vector<distr_t> &FF, distr_t &f_p, distr_t &m_p, distr_t &Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag, string Meson, bool UseJack, bool ConstFit, int t_fit) {

  int Nmeas= FF.size();
  int njacks = f_p.size();
  if((f_p.size() != m_p.size()) || (f_p.size() != Za_ov_Zv.size())) crash(" Number of jackknives are not the same for fp,mp and RC ZA/ZV");
  bool verbose=1;

 cout<<"Nmeas: "<<Nmeas<<endl;
 cout<<"Njacks: "<<njacks<<endl;
  //perform bootstrap fit

 


  bootstrap_fit<Xff_ChPT,Yff> bf(njacks);


  bf.Set_number_of_measurements(FF.size());
  bf.Set_verbosity(verbose);
  

  //Add parameters
  if(ff_type=="H1")  { bf.Add_par("a0", 0.2, 0.004); bf.Add_par("ak", 0.12 , 0.008); bf.Add_par("aq", 0.09 ,0.001);}
  else if(ff_type == "H2") { bf.Add_par("a0", 0.2, 0.005); bf.Add_par("ak",0.23 ,0.003); bf.Add_par("aq",-0.1 ,0.002);}
  else if(ff_type=="FA") {bf.Add_par("a0",0.032, 0.002); bf.Add_par("ak",0.036 ,0.001); bf.Add_par("aq",-0.0002 ,0.00001);}
  else if(ff_type=="FV") {bf.Add_par("a0",0.074*Za_ov_Zv.ave(), 0.001*Za_ov_Zv.ave()); bf.Add_par("ak",0.045*Za_ov_Zv.ave() ,0.0001*Za_ov_Zv.ave()); bf.Add_par("aq",0.02*Za_ov_Zv.ave() ,0.001*Za_ov_Zv.ave());}
  else crash("string : "+ff_type+" does not name a form factor");
  bf.Add_par("akq", 0.1, 0.01);
  bf.Add_par("a2kq",0.1,0.01);
  bf.Fix_par("akq",0.0);
  //bf.Fix_par("akq", 0.0);
  //bf.Fix_par("aq",0.0);
  //bf.Fix_par("akq",0.0);
  //bf.Fix_par("ak",0.0);

  //fix n release Za/Zv
  bf.Add_prior_par("Za_ov_Zv", Za_ov_Zv.ave(), Za_ov_Zv.err());
  if(W != "V") bf.Fix_par("Za_ov_Zv", 1.0);
  if(W =="V") bf.Fix_n_release("Za_ov_Zv");

  

  // if(ff_type=="H1" || ff_type=="H2") {bf.Fix_par("ak",0.0);bf.Fix_par("aq",0.0); bf.Fix_par("akq",0.0);}
  bf.ansatz =  [=](const Xff_ChPT &p, const Yff &ip) -> double {
    return  p.a0 + p.ak*pow(ip.xk,2) + p.aq*pow(ip.xq,2)+ p.akq*ip.xq*ip.xk +  p.a2kq*pow(ip.xq*ip.xk,2)  ;
  };
  bf.measurement = [=](const Xff_ChPT& p,const Yff& ip) -> double {
    return p.Za_ov_Zv*ip.ff;
  };
  bf.error =  [=](const Xff_ChPT& p,const  Yff &ip) -> double {
    return p.Za_ov_Zv*ip.ff_err;
  };


  //where to store fitted parameters
  boot_fit_data<Xff_ChPT> Bt_fit;

  //input parameters
  vector<vector<Yff>> data(njacks);
  for(auto &data_iboot: data) data_iboot.resize(Nmeas);

  
    for(int ijack=0;ijack<njacks;ijack++) {
      for(int imeas=0;imeas<Nmeas;imeas++) {
        data[ijack][imeas].Mpi = m_p.distr[ijack];
	data[ijack][imeas].fp = f_p.distr[ijack];
	data[ijack][imeas].ff = FF[imeas].distr[ijack];
	data[ijack][imeas].ff_err = FF[imeas].err();
	data[ijack][imeas].xk = fabs(mom[imeas].virt()/m_p.distr[ijack]) ;
	data[ijack][imeas].xq = sqrt( 1.0 + pow(data[ijack][imeas].xk,2) - mom[imeas].x_gamma(m_p).distr[ijack]) ;
      } 
    }


     //add prior 
  for(int ijack=0;ijack<njacks;ijack++) {
     bf.ib= &ijack;
     if(W=="A")  {  bf.Append_to_prior("Za_ov_Zv", 1.0, 1.0);}
     else if(W=="V") { bf.Append_to_prior("Za_ov_Zv", Za_ov_Zv.distr[ijack], Za_ov_Zv.err());}
     else crash("string W: "+W+" is not axial nor vector");
  }

    
  //fit
  bf.Append_to_input_par(data);
  Bt_fit =bf.Perform_bootstrap_fit();


  //define lambda functions
  for(int ijack=0;ijack<njacks;ijack++) {
    auto F = [=](double xk, double xq) -> double { return Bt_fit.par[ijack].a0 + Bt_fit.par[ijack].ak*pow(xk,2) + Bt_fit.par[ijack].aq*pow(xq,2) +Bt_fit.par[ijack].akq*xk*xq + Bt_fit.par[ijack].a2kq*pow(xk*xq,2);};
    fit_func.push_back(F);
  }

  //print data
  //set a grid in the xk, xq plane

  Vfloat xk_list, xq_list, FF_list, FF_err_list;
  for(int dxk=0; dxk < nsteps; dxk++) {
    for(int dxq=0; dxq < nsteps; dxq++) {
      double xxk = step_size*dxk;
      double xxq = step_size*dxq;
      distr_t func(UseJack);
      for(int ijack=0;ijack<njacks;ijack++) {
	func.distr.push_back( fit_func[ijack](xxk,xxq));
      }
      xk_list.push_back(xxk);
      xq_list.push_back(xxq);
      FF_list.push_back(func.ave());
      FF_err_list.push_back(func.err());
    }
  }


  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag);
   
  if(ConstFit) {
  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag+"/"+ff_type+".dat", "", "#xk   xq  val    err" );
  }
  else  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag+"/"+ff_type+"time_"+to_string(t_fit)+".dat", "", "#xk   xq  val    err" );
  
 

  return;

}