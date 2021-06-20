#ifndef __bootstrap_fit__
#define __bootstrap_fit__


#include "numerics.h"

using namespace std;

template <class T>
class boot_fit_data {

public:
  
  boot_fit_data() {}
  
  Vfloat chi2;
  Vfloat EDM;
  Vfloat N_it;
  vector<bool> IsValid;
  vector<T> par;
};





template <class T1, class T2> 

class bootstrap_fit : public ROOT::Minuit2::FCNBase  {
  
public:
 bootstrap_fit():  theErrorDef(1.0), NumberOfMeasurements(0) {
    this->PATH="chi2.out"; this->verbose=0;
  }
 bootstrap_fit(int nb) :  Input_pars(nb), theErrorDef(1.0), NumberOfMeasurements(0), nboots(nb) { this->PATH="chi2.out"; this->verbose=0;}
 
  virtual ~bootstrap_fit() {}
  function<double(const T1& p, const T2 &ip)> ansatz;
  function<double(const T1& p, const T2 &ip)> measurement;
  function<double(const T1& p, const T2 &ip)> error;
  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const Vfloat& par) const;

  boot_fit_data<T1> Perform_bootstrap_fit();
  
  
  void setErrorDef(double def) {theErrorDef = def;}
  void Add_par(string Name, double val, double err);
  void Add_prior_par(string Name, double val, double err);
  void Add_prior_pars(const vector<string>& Names);
  void Append_to_prior(string Name, double val, double err);
  void Fix_par(string Name, double val);
  void Fix_par(string Name) { Fix_par(Name, nan("1"));};
  void Append_to_input_par(vector<vector<T2>> values) {this->Input_pars = values;}
  void Append_to_input_par(vector<T2> val) { this->Input_pars[*ib] = val;}
  void Append_to_input_par(T2 val);
  void Clear_input_pars() {Input_pars.clear(); Input_pars.resize(nboots);}
  void Clear_prior(string A);
  void Clear_priors() { for (auto &pr: Priors) pr.second.clear();}
  void Delete_priors() {Priors.clear();}
  void Set_number_of_measurements(int Nens) { this->NumberOfMeasurements = Nens;}
  int Get_number_of_measurements() {return this->NumberOfMeasurements;}
  void Set_verbosity(int a) { this->verbose=a;}
  void Set_nboots(int nboots) { this->nboots=nboots;}
  void Fix_n_release(string A, double val) {this->To_release.insert(make_pair(A, val));}
  void Fix_n_release(string A) {Fix_n_release(A, nan("1"));}
  void Fix_n_release(const vector<string> &A) { for( auto &p:A) Fix_n_release(p);}
  int P(string A) const ;


  
 
  int* ib;

 private:
  ROOT::Minuit2::MnUserParameters PList;
  
  map<string, double> To_release;
  map<string, int> PNames;
  map<string, double> Fixed_pars;
  vector<vector<T2>> Input_pars;
  map<string, VPfloat> Priors;
  double theErrorDef;
  int NumberOfMeasurements;
  int nboots;
  int verbose;
  bool PRINT;
  string TAG;
  string PATH;
 
  
};


template <class T1, class T2> 
double bootstrap_fit<T1, T2>::operator()(const Vfloat& par) const {
    int pos=0;
    double chi2=0.0;
    ofstream PrintChi;
    if(PRINT) {
      PrintChi.open(PATH, ofstream::app); PrintChi.precision(10);
      PrintChi<<"#############################################"<<endl;
      PrintChi<<this->TAG<<endl;
    }
    while(pos<NumberOfMeasurements)  {
      T1 pp(par);
      // for(unsigned int pi=0;pi<par.size();pi++) if(isnan(par[pi])) crash("par nr: "+to_string(pi)+" is nan, bootstrap: "+to_string(*ib));
      double ansatz = this->ansatz(pp, Input_pars[*ib][pos]);
      double measurement = this->measurement(pp, Input_pars[*ib][pos]);
      double error = this->error(pp, Input_pars[*ib][pos]);
    
     
      
      double res = pow( (ansatz-measurement)/error,2);
      if(PRINT) PrintChi<<pos<<setw(20)<<ansatz<<setw(20)<<measurement<<setw(20)<<error<<setw(20)<<res<<endl;
      chi2+=res;
      pos++;
    }

    if(PRINT) PrintChi<<"chi2 w.o. priors: "<<chi2<<endl;
    //Add Gaussian Prior
    for(auto const& [key, val] : this->Priors) {
      if(PRINT) PrintChi<<key<<setw(20)<<par[P(key)]<<setw(20)<<val[*ib].first<<setw(20)<<val[*ib].second<<endl;
      chi2 += pow((par[P(key)]- val[*ib].first)/val[*ib].second,2);
    }
    if(PRINT) PrintChi<<"chi_2: "<<chi2<<endl;
    if(PRINT) PrintChi.close();
    return chi2;
}



template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Append_to_input_par(T2 val) {

  if(Input_pars[*ib].size() >= this->NumberOfMeasurements) crash("In bootstrap_fit::Append_to_input_par, data size larger than NumberOfMeasurements");
  else Input_pars[*ib].push_back(val);
  return;
}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Clear_prior(string A) {
 

  for(auto & [key,val] : Priors) {
    if(key == A) val.clear();
  }
  return;
}

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Append_to_prior(string Name, double val, double err) {

  map<string, VPfloat>::iterator it;
  it = Priors.find(Name);
  if(it == Priors.end()) crash("Could not find the prior parameter "+Name+" in Prior map");
  else {
    it->second.push_back(make_pair(val,err));
  }
  return;
  

}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_par(string Name, double val, double err) {
    PList.Add(Name, val, err);
    int Position= PNames.size();
    PNames.insert( pair<string,int>(Name, Position));
    return; 
  }

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_prior_pars(const vector<string> &Names) {

  
    for(unsigned int i=0; i<Names.size(); i++) {
      this->PList.Add(Names[i], 0.0, 1.0);
      int Position=PNames.size();
      PNames.insert( pair<string, int>(Names[i], Position));
      VPfloat start;
      Priors.insert(pair<string, VPfloat>(Names[i], start));
    }
    return;
}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_prior_par(string Name, double val, double err) {

   
  PList.Add(Name, val, err);
  int Position=PNames.size();
  PNames.insert( pair<string, int>(Name, Position));
  VPfloat start;
  Priors.insert(pair<string, VPfloat>(Name, start));
  
  return;
}

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Fix_par(string Name, double val) {
  map<string,double>::iterator it;
  it = Fixed_pars.find(Name);
  if(it != Fixed_pars.end()) {
    it->second = val;
    return;
  }
  
  Fixed_pars.insert(make_pair(Name, val)); 
  return;  
}


template <class T1, class T2> 
boot_fit_data<T1> bootstrap_fit<T1,T2>::Perform_bootstrap_fit() {

  boot_fit_data<T1> boot_result;

  int boot_index=0;
  ib = &boot_index;

  ROOT::Minuit2::MnMigrad migrad(*this, PList, 2);
  
  
  
  while(boot_index<nboots) {

     
    //fix parameters that must be fixed initially
    for(auto &name: To_release) {
      if(!isnan(name.second)) migrad.SetValue(name.first.c_str(), name.second);
      map<string, VPfloat>::iterator it;
      it= Priors.find(name.first);
      if( it != Priors.end()) migrad.SetValue(name.first.c_str(), it->second[*ib].first);
      migrad.Fix(name.first.c_str());
    }
    for(auto &name: Fixed_pars) {
      if(!isnan(name.second)) migrad.SetValue(name.first.c_str(), name.second);
      migrad.Fix(name.first.c_str());
    }
    

    
    //fit
    ROOT::Minuit2::FunctionMinimum chi2 = migrad();

    if(verbose) { //print chi^2 details
      this->PRINT= true;
      this->TAG= "Fixed priors";
      Vfloat PP(PNames.size());
      for(auto const& [key,val]: PNames) PP[val] = chi2.UserState().Value(key);
      this->operator()(PP);
      this->PRINT = false;
    }



    //release them
    for(auto & [key, val] : To_release) migrad.Release(key.c_str());

    
    //refit
    chi2 = migrad();
    boot_result.chi2.push_back(chi2.Fval());
    boot_result.EDM.push_back(chi2.Edm());
    boot_result.N_it.push_back(chi2.NFcn());
    boot_result.IsValid.push_back(chi2.IsValid());
    Vfloat PP(PNames.size());
    for(auto const& [key,val]: PNames) PP[val] = chi2.UserState().Value(key);
    T1 P(PP);
    boot_result.par.push_back(P);

    if(verbose) { //print chi^2 details
      this->PRINT= true;
      this->TAG= "Priors released";
      this->operator()(PP);
      this->PRINT = false;
    }
    

    
    

    cout<<"Chi2 for bootstrap: "<<boot_index<<" = "<<chi2.Fval()<<endl;
    cout<<"Parameters value:"<<endl;
    for(unsigned int par_id=0; par_id<PNames.size();par_id++) {
      for(auto const& [key,val] : PNames) {
	if(val==(signed)par_id) cout<<key<<"\t"<<chi2.UserState().Value(key)<<endl;
      }
    }
    if(verbose) {
      cout<<"Exit status: "<<chi2.IsValid()<<endl;
      cout<<"EDM :"<<chi2.Edm()<<endl;
      cout<<"N_iterations :"<<chi2.NFcn()<<endl;
    }
    cout<<"#################"<<endl;


      


    boot_index++;
  }

  
  return boot_result;
}

template <class T1, class T2> 
int bootstrap_fit<T1,T2>::P(string A) const {

  int occurrences=0;
  int pos=0;

  for(auto const& [key, val] : this->PNames)
    if(key== A) { pos = val; occurrences++;}

  if(occurrences!= 1) crash("In bootstrap_fit<T1,T2>::P(string A), either par_name: "+A+" does not exist or it's present multiple times");

  return pos;
  
}
 





double Boot_ave(Vfloat& A);
Pfloat Boot_ave_err(Vfloat& A);
double Boot_err(Vfloat& A) ;
Pfloat Boot_ave_err(VVfloat& A) ;
double Boot_ave(VVfloat& A) ;
double Boot_err(VVfloat& A);




































#endif

