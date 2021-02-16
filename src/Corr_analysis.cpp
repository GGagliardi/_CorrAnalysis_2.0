#include "../include/Corr_analysis.h"


using namespace std;

Vfloat CorrAnalysis::ASymm(const VVfloat &data, int t ) {

  if((signed)data.size() <= t) crash("In call to CorrAnalysis::ASymm(VVfloat, int t), index t is larger than VVfloat.size(). Exiting..."); 

  
   
  if(!Perform_Nt_t_average) return data[t];
  
  return Multiply_vector_by_scalar(Sum_vectors(data[t], Multiply_vector_by_scalar(data[(Nt-t)%Nt], (double)this->Reflection_sign)), 0.5);
 
}




distr_t_list  CorrAnalysis::corr_t(const VVfloat &corr_A, string Obs) {

  distr_t_list correlator(UseJack);
  VPfloat result;

  
  if((signed)corr_A.size() != Nt) {
    crash("corr_t in CorrAnalysis called with vector of size != Nt="+to_string(Nt));
    cout<<"Printing data..."<<endl;
    for(unsigned int t=0; t< corr_A.size();t++) printV(corr_A[t], "t="+to_string(t), 0);
  }
 
  if(UseJack) {
    Jackknife J(10, Njacks);
  
    
    for(int t=0; t<Nt;t++) 
      { correlator.distr_list.push_back(J.DoJack(1,ASymm(corr_A,t))); result.push_back( correlator[t].ave_err()); }
    
  }

  else { //use bootstrap
    Bootstrap B(Nboots, seed, corr_A[0].size());
    for(unsigned int t=0; t<corr_A.size(); t++) {correlator.distr_list.push_back(B.DoBoot(1, ASymm(corr_A,t))); result.push_back(correlator.ave_err(t));}     
  }
    
  


  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }

  return correlator;
}


distr_t_list CorrAnalysis::effective_mass_t(const VVfloat &corr_A, string Obs) {

  distr_t_list effective_mass_t(UseJack);
  VPfloat result;
   if((signed)corr_A.size() != Nt) {
    crash("effective_mass_t in CorrAnalysis called with vector of size != Nt="+to_string(Nt));
    cout<<"Printing data..."<<endl;
    for(unsigned int t=0; t< corr_A.size();t++) printV(corr_A[t], "t="+to_string(t), 0);
   }
 

  if(UseJack) { //use jackknife
    Jackknife J(10, Njacks);
    distr_t_list JackDistr_t(UseJack);
    for(int t=0; t<Nt;t++) JackDistr_t.distr_list.push_back(J.DoJack(1,ASymm(corr_A,t))); 
    
    for(int t=0; t<Nt;t++) {
      effective_mass_t.distr_list.emplace_back(UseJack);
      for(int ijack=0; ijack<Njacks;ijack++) effective_mass_t.distr_list[t].distr.push_back( Root_Brent( JackDistr_t.distr_list[t].distr[ijack]/JackDistr_t.distr_list[(t+1)%Nt].distr[ijack], t, Nt));
    }
  }
  
    
  else { //use bootstrap
    Bootstrap B(Nboots, seed, corr_A[0].size());
    distr_t_list BootDistr_t(UseJack);
    for( int t=0; t<Nt; t++) BootDistr_t.distr_list.push_back(B.DoBoot(1,ASymm(corr_A,t)));
    for(int t=0; t<Nt; t++) {
      effective_mass_t.distr_list.emplace_back(UseJack);
      for(int iboot=0; iboot<Nboots; iboot++) effective_mass_t.distr_list[t].distr.push_back( Root_Brent( BootDistr_t.distr_list[t].distr[iboot]/BootDistr_t.distr_list[(t+1)%Nt].distr[iboot], t, Nt));
    }
  }
    

  result = effective_mass_t.ave_err();


  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
  return effective_mass_t;
}


//overloading function effective_mass_t

distr_t_list CorrAnalysis::effective_mass_t(const distr_t_list &corr_A_distr, string Obs) {

 
  distr_t_list effective_mass_t(corr_A_distr.UseJack, corr_A_distr.size());

  vector<Pfloat> result;
 

   for(int t=0; t<corr_A_distr.size();t++) {
    if(corr_A_distr.distr_list[t].size() != corr_A_distr.distr_list[(t+1)%corr_A_distr.size()].size()) crash("Call to distr_t_list effective_mass_t(distr_t_list&) is invalid, distributions in distr_t_list do not have same size"); 
    for(int is=0; is<corr_A_distr.distr_list[t].size();is++) effective_mass_t.distr_list[t].distr.push_back( Root_Brent( corr_A_distr.distr_list[t].distr[is]/corr_A_distr.distr_list[(t+1)%corr_A_distr.size()].distr[is], t, corr_A_distr.size()));
  }

  result = effective_mass_t.ave_err();

 
  
  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(int t=0; t<corr_A_distr.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }

   
  return effective_mass_t;
}


distr_t_list CorrAnalysis::effective_slope_t(const VVfloat &corr_A, const VVfloat &corr_B, string Obs) {

  vector<Pfloat> result;

  distr_t_list effective_slope_t(UseJack);
  if((signed)corr_A.size() != Nt || (signed)corr_B.size() != Nt) crash("effective_slope_t in CorrAnalysis called with vectors of size != Nt");

 
  if(UseJack) {
    Jackknife J(10, Njacks);
    distr_t_list DM_Jack_t(UseJack), M_Jack_t(UseJack);
    for(int t=0; t<Nt;t++) {
      DM_Jack_t.distr_list.push_back( J.DoJack(1,ASymm(corr_A,t)));
      M_Jack_t.distr_list.push_back(J.DoJack(1,ASymm(corr_B,t)));
    }

    for( int t=0; t < Nt; t++) {
    effective_slope_t.distr_list.emplace_back(UseJack);
   
    for(int ijack=0; ijack<Njacks;ijack++) {

      double M_j = Root_Brent(M_Jack_t.distr_list[(t-1+Nt)%Nt].distr[ijack]/M_Jack_t.distr_list[t].distr[ijack], t-1, Nt);
      double der_back = DM_Jack_t[t].distr[ijack]/M_Jack_t[t].distr[ijack] - DM_Jack_t[(t-1+Nt)%Nt].distr[ijack]/M_Jack_t[(t-1+Nt)%Nt].distr[ijack];
      
      effective_slope_t.distr_list[t].distr.push_back(der_back/(FTAN(M_j, t,Nt)));

    }
    }
  }
  
 
  else { //use bootstrap
    Bootstrap B(Nboots, seed, corr_A[0].size());
    distr_t_list DM_Boot_t(UseJack), M_Boot_t(UseJack);
    for(int t=0; t<Nt; t++) {
      DM_Boot_t.distr_list.push_back(B.DoBoot(1,ASymm(corr_A,t)));
      M_Boot_t.distr_list.push_back(B.DoBoot(1,ASymm(corr_B,t)));
    }
    for(int t=0; t<Nt;t++) {
      effective_slope_t.distr_list.emplace_back(UseJack);
      for(int iboot=0;iboot<Nboots;iboot++) {
	double M_j = Root_Brent(M_Boot_t.distr_list[(t-1+Nt)%Nt].distr[iboot]/M_Boot_t.distr_list[t].distr[iboot],t-1,Nt);
	
	double der_back = DM_Boot_t[t].distr[iboot]/M_Boot_t[t].distr[iboot] - DM_Boot_t[(t-1+Nt)%Nt].distr[iboot]/M_Boot_t[(t-1+Nt)%Nt].distr[iboot];
	effective_slope_t.distr_list[t].distr.push_back( der_back/(FTAN(M_j,t,Nt)));
      }
    }     
  }


  result = effective_slope_t.ave_err();


  


  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
  return effective_slope_t;
}


//overloading function effective_slope_t


distr_t_list CorrAnalysis::effective_slope_t(const distr_t_list &corr_A_distr, const distr_t_list &corr_B_distr, string Obs) {

  distr_t_list effective_slope_t(corr_A_distr.UseJack, corr_A_distr.size());

  VPfloat result;
 
  if(corr_A_distr.size() != corr_B_distr.size()) crash("Call to  distr_t_list effective_slope_t(distr_t_list&, distr_t_list&) is invalid, the two distributions list do not have same size");

  for(int t=0; t<corr_A_distr.size();t++) {


    if(corr_A_distr.distr_list[t].size() != corr_B_distr.distr_list[t].size()) crash("Call to  distr_t_list effective_slope_t(distr_t_list&, distr_t_list&) is invalid, distributions in distr_t_list do not have same size");
    
    for(int is=0; is<corr_A_distr.distr_list[t].size();is++) {

      double M_j = Root_Brent(corr_B_distr.distr_list[(t-1+corr_B_distr.size())%corr_B_distr.size()].distr[is]/corr_B_distr.distr_list[t].distr[is], t-1, corr_B_distr.size());
      double der_back = corr_A_distr.distr_list[t].distr[is]/corr_B_distr.distr_list[t].distr[is] - corr_A_distr.distr_list[(t-1+Nt)%Nt].distr[is]/corr_B_distr.distr_list[(t-1+Nt)%Nt].distr[is];
      effective_slope_t.distr_list[t].distr.push_back( der_back/(FTAN(M_j, t, Nt)));
      
    }

    
  }

  result = effective_slope_t.ave_err();

    if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(int t=0; t<corr_A_distr.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
   
  return effective_slope_t;
}

distr_t_list CorrAnalysis::effective_slope_t_2nd_ord(const distr_t_list& corr_A_distr,const distr_t_list& corr_B_distr,const distr_t_list& corr_C_distr, string Obs) {

  
  distr_t effective_slope_t_2nd_ord(corr_A_distr.UseJack);

  double Mj=0;

  distr_t m_eff = Fit_distr(effective_mass_t(corr_B_distr,""));

  distr_t z(corr_A_distr.UseJack);

 


 
  if(corr_A_distr.size() != corr_B_distr.size()) crash("Call to  distr_t_list effective_slope_t(distr_t_list&, distr_t_list&) is invalid, the  distributions list do not have same size");
  if(corr_A_distr.size() != corr_C_distr.size()) crash("Call to  distr_t_list effective_slope_t(distr_t_list&, distr_t_list&) is invalid, the  distributions list do not have same size");


  distr_t_list ratio_A = corr_A_distr/corr_B_distr;
  distr_t_list ratio_C = corr_C_distr/corr_B_distr;

  Vfloat Y_conn_err, Y_disc_err;

  for(int t=10;t<=Nt/2;t++) { Y_conn_err.push_back(ratio_A.err()[t]); Y_disc_err.push_back(ratio_C.err()[t]);}

  for(int is=0;is<corr_A_distr.distr_list[0].size();is++) {
    Vfloat X, Y_conn, Y_disc ;
    for(int t=10;t<=Nt/2;t++) {
      if(corr_A_distr.distr_list[t].size() != corr_A_distr.distr_list[0].size()) crash("In CorrAnalysis::effective_slope_t_2nd_ord jackknife samples in corr_A are not valid");
      if(corr_C_distr.distr_list[t].size() != corr_A_distr.distr_list[0].size()) crash("In CorrAnalysis::effective_slope_t_2nd_ord jackknife samples in corr_C are not valid");
      X.push_back(t);
      Y_conn.push_back( ratio_A.distr_list[t].distr[is]);
      Y_disc.push_back( ratio_C.distr_list[t].distr[is]);
    }
    
    
    Mj = m_eff.distr[is];
    auto anz = [&](const Vfloat &par, double t) -> double { return par[0] + par[1]*pow((Nt/2 -t),2) + par[2]*(Nt/2 -t)*tanh( (Nt/2-t)*Mj);};
    T_fit conn(X,Y_conn, Y_conn_err);    
    conn.ansatz = anz;
    T_fit disc(X, Y_disc, Y_disc_err);
    disc.ansatz = anz;
    conn.add_pars({10,-50,100});
    disc.add_pars({10,-50,100});
    fit_t_res conn_res= conn.fit();
    fit_t_res disc_res= disc.fit();
    z.distr.push_back( conn_res.pars[1]/disc_res.pars[1]);
    //effective_slope_t_2nd_ord.distr.push_back( -1*conn_res.pars[2]  + conn_res.pars[1]*disc_res.pars[2]/disc_res.pars[1]);
    //cout<<"ch2 : "<< conn_res.chi2/(Nt/2-(10) - 3)<<" "<<disc_res.chi2/(Nt/2-(10)-3)<<endl;
  }

  return effective_slope_t( corr_A_distr- z*corr_C_distr, corr_B_distr, Obs);


}




distr_t_list CorrAnalysis::effective_slope2_t(const VVfloat &corr_A, const VVfloat &corr_B, string Obs) {

  vector<Pfloat> result;
  
  distr_t_list effective_slope2_t(UseJack);
  if((signed)corr_A.size() != Nt || (signed)corr_B.size() != Nt) crash("effective_slope_t in CorrAnalysis called with vectors of size != Nt");

  
  if(UseJack) {
    Jackknife J(10, Njacks);
    distr_t_list DM_Jack_t(UseJack), M_Jack_t(UseJack);
    for(int t=0; t<Nt;t++) {
      DM_Jack_t.distr_list.push_back( J.DoJack(1,ASymm(corr_A,t)));
      M_Jack_t.distr_list.push_back(J.DoJack(1,ASymm(corr_B,t)));
    }
    
    for(int t=0; t < Nt; t++) {
      effective_slope2_t.distr_list.emplace_back(UseJack);
      for(int ijack=0; ijack<Njacks;ijack++) {

	double M_j = Root_Brent(M_Jack_t.distr_list[(t-1+Nt)%Nt].distr[ijack]/M_Jack_t.distr_list[t].distr[ijack], t-1, Nt);
	double der_back = DM_Jack_t[t].distr[ijack]/M_Jack_t[t].distr[ijack] - DM_Jack_t[(t-1+Nt)%Nt].distr[ijack]/M_Jack_t[(t-1+Nt)%Nt].distr[ijack];
	effective_slope2_t.distr_list[t].distr.push_back( 2.0*M_j*der_back/FTAN(M_j, t, Nt));
      }
    }
  }
  
    
  else { //use bootstrap
    Bootstrap B(Nboots, seed, corr_A[0].size());
    distr_t_list DM_Boot_t(UseJack), M_Boot_t(UseJack);
    for(int t=0; t<Nt; t++) {
      DM_Boot_t.distr_list.push_back(B.DoBoot(1,ASymm(corr_A,t)));
      M_Boot_t.distr_list.push_back(B.DoBoot(1,ASymm(corr_B,t)));
    }
    for(int t=0; t<Nt;t++) {
      effective_slope2_t.distr_list.emplace_back(UseJack);
      for(int iboot=0;iboot<Nboots;iboot++) {

	double M_j = Root_Brent(M_Boot_t.distr_list[(t-1+Nt)%Nt].distr[iboot]/M_Boot_t.distr_list[t].distr[iboot], t-1, Nt);
	double der_back = DM_Boot_t[t].distr[iboot]/M_Boot_t[t].distr[iboot] - DM_Boot_t[(t-1+Nt)%Nt].distr[iboot]/M_Boot_t[(t-1+Nt)%Nt].distr[iboot];
	effective_slope2_t.distr_list[t].distr.push_back( 2.0*M_j*der_back/FTAN(M_j, t, Nt));
	
      }
    }     
  }
    
  result = effective_slope2_t.ave_err();


  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
  return effective_slope2_t;
}





distr_t_list CorrAnalysis::residue_t(const VVfloat &corr_A, string Obs) {

  distr_t_list residue(UseJack);
  VPfloat result;
  if((signed)corr_A.size() != Nt) crash("residue_t in CorrAnalysis called with vector of size != Nt");

  distr_t_list corr_distr = corr_t(corr_A, "");
  
  distr_t effective_mass_fit_distr= Fit_distr(effective_mass_t(corr_distr, ""));
  

  distr_t_list analytic_factor(UseJack,Nt);
  for(int t=0;t<Nt;t++) {
    for(int is=0; is < effective_mass_fit_distr.size(); is++) {
      double el = effective_mass_fit_distr.distr[is];
      analytic_factor.distr_list[t].distr.push_back( (exp(-el*t) +Reflection_sign*exp(-el*(Nt-t)))/(2*el));
    }
  }

  
  residue = corr_distr/analytic_factor;

  result = residue.ave_err();

  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }

  return residue;
}
 


//overloading function residue_t

distr_t_list CorrAnalysis::residue_t(const distr_t_list &corr_A_distr, string Obs) {

  distr_t_list residue(corr_A_distr.UseJack, corr_A_distr.size());
  distr_t effective_mass_fit_distr = Fit_distr(effective_mass_t(corr_A_distr, ""));


  VPfloat result;

  distr_t_list analytic_factor(UseJack,Nt);
  for(int t=0;t<Nt;t++) {
    for(int is=0; is < effective_mass_fit_distr.size(); is++) {
      double el = effective_mass_fit_distr.distr[is];
      analytic_factor.distr_list[t].distr.push_back( (exp(-el*t) +Reflection_sign*exp(-el*(Nt-t)))/(2*el));
    }
  }

  residue = corr_A_distr/analytic_factor;

  

  result = residue.ave_err();

  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for( int t=0; t<corr_A_distr.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }


  
  return residue;
 


}



distr_t_list CorrAnalysis::decay_constant_t(const VVfloat &corr_A, string Obs) {

  distr_t_list decay_constant_t(UseJack,Nt);
  VPfloat result;

  distr_t_list residue = residue_t(corr_A, "");
  distr_t effective_mass_fit_distr = Fit_distr( effective_mass_t(corr_A, ""));

     
  
  for(int t=0;t<Nt;t++) {
    for(int is=0; is < effective_mass_fit_distr.size(); is++) {
      double el = effective_mass_fit_distr.distr[is];
      double Z = residue.distr_list[t].distr[is];
      decay_constant_t.distr_list[t].distr.push_back( sqrt(Z)/(el*sinh(el)));
    }
  }


  
  result = decay_constant_t.ave_err();
  
  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(unsigned int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
  
  return decay_constant_t;
}

//overloading function decay_constant_t

distr_t_list CorrAnalysis::decay_constant_t(const distr_t_list &corr_A, string Obs) {

  distr_t_list decay_constant_t(UseJack,Nt);
  VPfloat result;


  distr_t_list residue = residue_t(corr_A, "");
  distr_t effective_mass_fit_distr = Fit_distr( effective_mass_t(corr_A, ""));
  
  
  for(int t=0;t<Nt;t++) {
    for(int is=0; is < effective_mass_fit_distr.size(); is++) {
      double el = effective_mass_fit_distr.distr[is];
      double Z = residue.distr_list[t].distr[is];
      decay_constant_t.distr_list[t].distr.push_back( sqrt(Z)/(el*sinh(el)));
    }
  }
  

  result = decay_constant_t.ave_err();
  
  if(Obs != "") {
    ofstream Print(Obs+".t", ofstream::out);
    Print.precision(10);
    for(int t=0; t<corr_A.size(); t++) Print<<t<<setw(20)<<result[t].first<<setw(20)<<result[t].second<<endl;
    Print.close();
  }
  
  return decay_constant_t;
}








Pfloat CorrAnalysis::Fit_(const distr_t_list& M_distr) {

  
  Vfloat err;
  for(int t=Tmin; t<=Tmax; t++) err.push_back(M_distr.err(t));

  distr_t Fit_distr(M_distr.UseJack);

  if(M_distr.size() ==0) crash("Fit_distr called with an empty distr_t_list element");
  int distr_size= M_distr.distr_list[0].size();
  
  for(int i_distr=0; i_distr<distr_size;i_distr++) {
    Vfloat M_distr_i;
    for(int t=Tmin;t<=Tmax;t++) {
      if(M_distr.distr_list[t].size() != distr_size) crash("In Fit_distr, the elements of distr_t_list do not have same sizes");
      M_distr_i.push_back( M_distr.distr_list[t].distr[i_distr]);
    }
    Fit_distr.distr.push_back(DoConstantFit(M_distr_i,err));
  }

  return Fit_distr.ave_err();
}
 




distr_t CorrAnalysis::Fit_distr(const distr_t_list& M_distr) {

  Vfloat err;
  for(int t=Tmin; t<=Tmax; t++) err.push_back(M_distr.err(t));

  distr_t Fit_distr(M_distr.UseJack);

  if(M_distr.size() ==0) crash("Fit_distr called with an empty distr_t_list element");
  int distr_size= M_distr.distr_list[0].size();
  for(int i_distr=0; i_distr<distr_size;i_distr++) {
    Vfloat M_distr_i;
    for(int t=Tmin;t<=Tmax;t++) {
      if(M_distr.distr_list[t].size() != distr_size) crash("In Fit_distr, the elements of distr_t_list do not have same sizes");
      M_distr_i.push_back( M_distr.distr_list[t].distr[i_distr]);
    }
    Fit_distr.distr.push_back(DoConstantFit(M_distr_i,err));

  }

  
  return Fit_distr;
}
