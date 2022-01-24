#include "../include/Meson_mass_extrapolation.h"

const bool verb=true;
using namespace std;

class Ov {

public:
  Ov()  {}
  
  double Y,Y_err, X;
};

class fpar {

public:
  fpar() {}
  fpar(const Vfloat &par) {
    if((signed)par.size() != 2) crash("In Obs_extrapolation_meson_mass, class fpar has size different from two.");
    a0=par[0];
    D=par[1];
  }

  double a0,D;

};



distr_t Obs_extrapolation_meson_mass( vector<distr_t> &Meas, vector<distr_t> &Op, double Op_phys, string Print_dir, string Tag) {

  //Get number of jackknives
  double Njacks= Meas[0].size();
  double Nmeas = Meas.size();
  cout<<"Njacks: "<<Njacks<<endl;
  cout<<"Nmeas: "<<Nmeas<<endl;
  cout<<"Print directory: "<<Print_dir<<endl;
  cout<<"Tag: "<<Tag<<endl;
  cout<<"##########"<<endl;
  cout<<"Meson obs values (Exp: "<<Op_phys<<" )"<<endl;
  for(int imeas=0;imeas<Nmeas;imeas++) cout<<Op[imeas].ave()<<" +- "<<Op[imeas].err()<<" ; "<<Meas[imeas].ave()<<" +- "<<Meas[imeas].err()<<endl;
  cout<<"##########"<<endl;
  bootstrap_fit<fpar,Ov> bf(Njacks);
  bf.set_warmup_lev(2); //sets warmup

  bf.Set_number_of_measurements(Meas.size());
  bf.Set_verbosity(verb);
  double a0_guess =0;
  for(auto &ms: Meas) a0_guess += ms.ave()/Nmeas;
  bf.Add_par("a0", a0_guess, 1e-2*a0_guess);
  bf.Add_par("D", -1.0e-1, 1e-3);

  bf.ansatz =  [&](const fpar &p, const Ov &ip) -> double {
		 return p.a0*(1.0 + p.D*(ip.X - Op_phys));
	       };

  bf.measurement = [&](const fpar &p, const Ov &ip) ->double {
		     return ip.Y;
		   };
  bf.error = [&](const fpar &p, const Ov &ip) -> double {
	       return ip.Y_err;
	     };


  //meas and errs

  
  //fill the data
  vector<vector<Ov>> data(Njacks);
  for(auto &data_ij: data) data_ij.resize(Nmeas);
  //allocate space for output result
  boot_fit_data<fpar> Bt_fit;

  for(auto &data_iboot: data) data_iboot.resize(Nmeas);
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int i=0;i<Nmeas;i++) {
      data[ijack][i].Y = Meas[i].distr[ijack];
      data[ijack][i].Y_err = Meas[i].err();
      data[ijack][i].X = Op[i].distr[ijack];
    } 
  }


    
  //append
  bf.Append_to_input_par(data);

  //fit
  Bt_fit= bf.Perform_bootstrap_fit();
  Bt_fit.ch2_ave();


  //retrieve fit parameter a0

  distr_t a0_fit, D_fit;

  for(int ijack=0;ijack<Njacks;ijack++) {a0_fit.distr.push_back( Bt_fit.par[ijack].a0); D_fit.distr.push_back(Bt_fit.par[ijack].D);}


  if(Print_dir != "") {
    boost::filesystem::create_directory(Print_dir+"/extr");
    vector<double> X, Xerr, Y,Yerr;
    for(int imeas=0;imeas<Nmeas;imeas++) {
      X.push_back( Op[imeas].ave());
      Xerr.push_back( Op[imeas].err());
      Y.push_back( Meas[imeas].ave());
      Yerr.push_back( Meas[imeas].err());
    }
   
    X.push_back(Op_phys);
    Xerr.push_back(0.0);
    Y.push_back(a0_fit.ave());
    Yerr.push_back(a0_fit.err());

    Print_To_File({}, {X,Xerr, Y, Yerr}, Print_dir+"/extr/"+Tag+".meas", "", "# X Xerr Y Y err");

    //print fitting function
    int Npoints = 300;
    double xmin= 0.5*Op_phys;
    double xmax= 1.5*Op_phys;
    Vfloat Xpoints;
    distr_t_list Fit_func;
    for(int ipoint=0;ipoint<Npoints;ipoint++) {
      double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
      Fit_func.distr_list.push_back(a0_fit*(1.0 + D_fit*(xp-Op_phys)));
      Xpoints.push_back(xp);
    }

    Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},Print_dir+"/extr/"+Tag+".fit_func", "", "");
    


  }

  return a0_fit;
}
