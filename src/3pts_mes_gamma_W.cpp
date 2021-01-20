#include "../include/3pts_mes_gamma_W.h"



using namespace std;
namespace plt = matplotlibcpp;

const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nbranches = 8;
const int Nboots= 100;
const bool UseJack=1;
const int nboots=150;
const int Njacks=5;
const int sm_lev=0;
bool Include_k0_noise=1;

class ExpFit {

 
public:
  ExpFit(const Vfloat &par)  { if((signed)par.size() != 7) crash("par size != 7");
    a1=par[0];
    a2=par[1];
    a3=par[2];
    a4=par[3];
    a5=par[4];
    a6=par[5];
    a7=par[6];
  }

  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
  double a6;
  double a7;
};

class Val{
public:
  Val() {}
  double t;
  double meas;
  double err;
};

void Get_Tmin_Tmax(string corr_type, string Ens_tag, CorrAnalysis &corr, double xg, string W) {

  if(corr_type=="2pt") {
    corr.Tmax = corr.Nt/2 -4;
    if(Ens_tag.substr(0,1) =="A" || Ens_tag.substr(0,1) == "C" || Ens_tag.substr(0,1) == "B" || Ens_tag.substr(0,1)== "G" || Ens_tag.substr(0,1)=="D" || Ens_tag.substr(0,1) =="E" ) corr.Tmin = 15;
    else if(Ens_tag.substr(0,1) =="B") corr.Tmin = 16;
    else corr.Tmin = 19;
  }

  if(corr_type=="3pt") {
    if(Ens_tag.substr(0,6)=="A40.32") {
      if(W=="A") { 
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 17; corr.Tmin=13;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if (xg >0.09 && xg < 0.11) { corr.Tmax = 29; corr.Tmin= 23;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 27; corr.Tmin= 22;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 24; corr.Tmin= 16;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 25; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 23; corr.Tmin= 12;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 17; corr.Tmin= 13;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      
    }
    else  if(Ens_tag.substr(0,6)=="B40.32") {
      if(W=="A") { 
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin=12;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if (xg >0.09 && xg < 0.11) { corr.Tmax = 20; corr.Tmin= 15;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 20; corr.Tmin= 15;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 20; corr.Tmin= 12;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 16; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 16; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      
    }
    else if(Ens_tag.substr(0,6)=="A40.24" || Ens_tag.substr(0,6)=="B40.24" || Ens_tag.substr(0,6)=="C40.24") {
      if(W=="A") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin= 9;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 9;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 14; corr.Tmin= 8;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 8;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 14; corr.Tmin= 8;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11) { corr.Tmax = 18; corr.Tmin= 14;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 18; corr.Tmin= 14;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 12;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 13; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 10;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 13; corr.Tmin= 10;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else crash("cannot find the ensemble");
    }
      else if(Ens_tag.substr(0,6)=="A40.40") {
      if(W=="A") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin= 9;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 9;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 14; corr.Tmin= 8;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 8;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 14; corr.Tmin= 8;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 20; corr.Tmin= 15;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 12;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 15; corr.Tmin= 11;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 10;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 14; corr.Tmin= 9;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else crash("cannot find the ensemble");
    }

    
    else { corr.Tmin= 9; corr.Tmax=16;}
  }

  return;

}


void Add_to_mom_list(pt3_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt3_momenta> A;

  for(auto & c: header.comb) A.emplace_back(Vfloat{c.th0[0], c.th0[1], c.th0[2]},Vfloat{c.ths[0], c.ths[1], c.ths[2]},Vfloat{c.tht[0], c.tht[1], c.tht[2] }, c.mu1, c.mu2, c.off,L, header.tmax); 
  M.mom.push_back(A);
  return;
}


void Add_to_mom_list(pt2_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt2_momenta> A;

  for(auto & c: header.comb) A.emplace_back( Vfloat{c.th0[0], c.th0[1], c.th0[2]}, Vfloat{c.ths[0], c.ths[1], c.ths[2]}, c.mu1, c.mu2, L, c.i0, c.is);

  M.mom.push_back(A);

  return;
}


distr_t_list V_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, distr_t& Meson_mass, pt3_momenta& Mom,int twall) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};

  



  if(distr_mom_k.size() != 4 || distr_mom_0.size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas=4");
  int size= distr_mom_k[0][0].size();
  if(size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has size zero, exiting....");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has zero distr_size, exiting....");


  distr_t_list V_ave(UseJack, size, distr_size);

  Vfloat exp_Nt_t;
  for(int t=0; t<Mom.Nt();t++) exp_Nt_t.push_back( exp( fabs((twall-t))*Mom.Egamma()));
  //for(int t=0;t<Mom.Nt(); t++) exp_Nt_t.push_back( exp( ((Mom.Nt()/2.0)-t)*Mom.Egamma()) - exp((t- Mom.Nt()/2.0)*Mom.Egamma()));

  bool is_mom_zero=false;
  Vfloat theta_t_old(Mom.Theta(2));
  if(Mom.k_mod() < eps(5)) { is_mom_zero=true; Vfloat theta_t = {0.0, 0.0, 1.0}; Mom.Set_theta(theta_t, 2);}

  
  
  for(int alpha=1; alpha <= 2; alpha ++) {
   
    if(distr_mom_k[alpha].size() != 4 || distr_mom_0[alpha].size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
    for(int mu=0; mu<4;mu++) {

      if(mu != alpha) {
    
	distr_t_list den1(3, Meson_mass);
	distr_t_list den2(3, Meson_mass);
	
	
      den1 = den1*Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3),Mom.k()),-1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3), Mom.p()), Mom.Egamma());
    
      den2 = den2*Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3),Mom.k()), -1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3), Mom.p()), Mom.Egamma());
    
      V_ave = V_ave+ e[1][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den1[alpha-1];

      V_ave = V_ave+ e[2][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den2[alpha-1];


      }

     
    }
  }


 
  
  
  cout<<"V##########"<<endl;
  distr_t_list DIFF2 = 2.0*( distr_mom_k[2][1]*exp_Nt_t -((double)Include_k0_noise)*distr_mom_0[2][1] - distr_mom_k[1][2]*exp_Nt_t + ((double)Include_k0_noise)*distr_mom_0[1][2]);
  distr_t_list DIFF= DIFF2/(Mom.k()[2]*Meson_mass);
  //for(int t=0; t < V_ave.size(); t++) cout<<V_ave.distr_list[t].ave()<<" "<<V_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
  cout<<"############"<<endl;

  if(is_mom_zero) Mom.Set_theta(theta_t_old, 2);
  
  return V_ave;
  //return DIFF2;


}

distr_t_list A_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};



  if(distr_mom_k.size() != 4) crash("A_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas= 4");
  for(auto & distr_mom_k_alpha : distr_mom_k) if(distr_mom_k_alpha.size() != 4) crash("A_ave unpolarized called with a with a vv<distr_t_list>> having an element of size != 4");
  int size = distr_mom_k[0][0].size();
  if(size==0) crash(" A_ave_unpolarized called with Nt=0");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size ==0) crash(" A_ave unpolarized called with a vv<distr_t_list>> having zero statistics");
  distr_t_list A_ave(UseJack, size, distr_size);

 
  for(int alpha=1; alpha <= 2; alpha ++) {
 
    if(distr_mom_k[alpha].size() != 4) crash(" A_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
 
    for(int mu=0; mu< 4; mu ++) {

 
      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[1][mu]/e[1][alpha]);

      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[2][mu]/e[2][alpha]);
      

   
    }
  }


  /*
  cout<<"A##########"<<endl;
  distr_t_list DIFF = 2.0*( distr_mom_k[1][1] +distr_mom_k[2][2]);
  for(int t=0; t < A_ave.size(); t++) cout<<A_ave.distr_list[t].ave()<<" "<<A_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
  cout<<"############"<<endl;
  */

  //return distr_mom_k[1][1];
  return A_ave;
  
}


distr_t_list H_V(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, pt3_momenta& Mom) {


  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  
  if(distr_mom_k.size() != 4 || distr_mom_0.size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas=4");
  int size= distr_mom_k[0][0].size();
  if(size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has size zero, exiting....");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has zero distr_size, exiting....");
  
  
  distr_t_list H_ave(UseJack, size, distr_size);

  Vfloat exp_Nt_t;
  for(int t=0; t<Mom.Nt();t++) exp_Nt_t.push_back( exp( fabs((Mom.Nt()/2.0-t))*Mom.Egamma()));
  //for(int t=0;t<Mom.Nt(); t++) exp_Nt_t.push_back( exp( ((Mom.Nt()/2.0)-t)*Mom.Egamma()) - exp((t- Mom.Nt()/2.0)*Mom.Egamma()));

  
  for(int alpha=1; alpha < 2; alpha ++) {
   
    if(distr_mom_k[alpha].size() != 4 || distr_mom_0[alpha].size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
    for(int mu=0; mu<4;mu++) {

      if(mu != alpha) {
    
    
      H_ave = H_ave+ e[1][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu]);

      H_ave = H_ave+ e[2][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu]);

      }

     
    }
  }

  

  return H_ave;

}













void Compute_form_factors(string Meson) {




  string dir= "../3pt_data/"+Meson;
  int Nens=0;
  vector<string> Ens_tag;
  

  //read number of ensembles and store the paths

  boost::filesystem::directory_iterator end_itr;

  for (boost::filesystem::directory_iterator itr(dir); itr != end_itr; ++itr) {
  

      if(!boost::filesystem::is_regular_file(itr->path())) {

	Nens++;
	Ens_tag.push_back(itr->path().string().substr(dir.length()+1));
        }
  }



  
  pt3_momenta_list mom3_l;
  pt2_momenta_list mom2_l;


  vector<vector<distr_t_list>> bar_R_A_distr(Nens);
  vector<vector<distr_t_list>> bar_R_V_distr(Nens);
  vector<vector<distr_t_list>> bar_R_A_no_bar_distr(Nens);
  vector<vector<distr_t_list>> bar_R_V_no_bar_distr(Nens);
  vector<vector<distr_t_list>> HV_k(Nens);
  vector<vector<distr_t>> bar_R_A_fit_distr(Nens);
  vector<vector<distr_t>> bar_R_V_fit_distr(Nens);
  vector<vector<distr_t_list>> fp_distr_list(Nens);
  vector<vector<distr_t_list>> m_distr_list(Nens);
  vector<vector<distr_t>> fp_fit_distr(Nens);
  vector<vector<distr_t>> m_fit_distr(Nens);
  vector<vector<distr_t_list>> overlap_distr_list(Nens);
  vector<vector<distr_t_list>> overlap_smeared_distr_list(Nens);
  vector<vector<distr_t>> overlap_fit_distr(Nens);
  vector<vector<distr_t>> overlap_smeared_fit_distr(Nens);
  
  

  for (int i_ens=0; i_ens<Nens;i_ens++) {


    FILE *stream_3pt;
    FILE *stream_2pt;

    struct header_virph header_3pt;
    struct header_virph header_2pt;

   

    string Path3pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat";
    string Path2pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat2";



    stream_3pt= fopen(Path3pt.c_str(), "rb");
    stream_2pt= fopen(Path2pt.c_str(), "rb");



    read_header_bin(stream_3pt, header_3pt);
    read_header_bin(stream_2pt, header_2pt);


    //print size of header file
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 2pt (Bytes): "<<header_2pt.header_size<<endl;
    cout<<"Number of configs (from 2pt): "<<Get_number_of_configs_2pt(stream_2pt, header_2pt)<<endl;
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 3pt (Bytes): "<<header_3pt.header_size<<endl;
    cout<<"Number of configs (from 3pt): "<<Get_number_of_configs_3pt(stream_3pt, header_3pt)<<endl;
  


    //Forward informations to 3pt_momenta_list
    double ml, l, Nt;
    Read_pars_from_ensemble_tag(Ens_tag[i_ens], ml, l , Nt);
    Add_to_mom_list(mom3_l, header_3pt, l);
    Add_to_mom_list(mom2_l, header_2pt, l);
    mom2_l.Nt.push_back(Nt);
    mom3_l.Nt.push_back(Nt);

    //init CorrAnalysis class
    CorrAnalysis corr(UseJack, Get_number_of_configs_3pt(stream_3pt, header_3pt), nboots);
    //CorrAnalysis corr(UseJack, Njacks, nboots);
    corr.Nt=Nt;

    

    //2pts
    Get_Tmin_Tmax("2pt", Ens_tag[i_ens], corr, 0, "V");
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/form_factors");
    boost::filesystem::create_directory("../data/form_factors/"+Meson);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]);

    int ncomb2pt= header_2pt.ncomb;

    
    cout<<"number of combs in 2pt: "<<ncomb2pt<<endl;

    for(int icomb2pt=0; icomb2pt<ncomb2pt;icomb2pt++) {

      auto c = header_2pt.comb[icomb2pt];
      string tag=to_string(icomb2pt);
      overlap_distr_list[i_ens].push_back ( corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,0), ""));
      overlap_smeared_distr_list[i_ens].push_back( corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,sm_lev), ""));
      overlap_fit_distr[i_ens].push_back( corr.Fit_distr(overlap_distr_list[i_ens][icomb2pt]));
      overlap_smeared_fit_distr[i_ens].push_back( corr.Fit_distr( overlap_smeared_distr_list[i_ens][icomb2pt]));
      fp_distr_list[i_ens].push_back( (c.mu1+c.mu2)*corr.decay_constant_t( Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, 0), "../data/form_factors/"+Meson+"/fp_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));
      m_distr_list[i_ens].push_back( corr.effective_mass_t(Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, sm_lev), "../data/form_factors/"+Meson+"/meson_mass_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));
      fp_fit_distr[i_ens].push_back(corr.Fit_distr(fp_distr_list[i_ens][icomb2pt]));
      m_fit_distr[i_ens].push_back(corr.Fit_distr(m_distr_list[i_ens][icomb2pt]));
    }

    //end anaysis 2pts

    //beginning analysis 3pts

    int ncomb3pt= header_3pt.ncomb;



   
    //loop over kinematics

      

    for(int icomb3pt=0; icomb3pt<ncomb3pt;icomb3pt++) {

      vector<vector<distr_t_list>> distr_V(4), distr_V_k0(4), distr_A(4), distr_A_k0(4);
      vector<distr_t_list> distr_corr;
      for(int i=0;i<4;i++) { distr_V[i].resize(4); distr_V_k0[i].resize(4); distr_A[i].resize(4); distr_A_k0[i].resize(4);}
      
      auto c = header_3pt.comb[icomb3pt];
      
      int icomb_k0= Get_comb_k0(header_3pt, icomb3pt);

      int symmetric_comb=Get_symmetric_comb(header_3pt, icomb3pt);
      int symmetric_comb_k0 = Get_symmetric_comb(header_3pt, icomb_k0);

      cout<<"###BEG###"<<endl;
      cout<<"icomb: "<<icomb3pt<<" th_t: "<<c.tht[2]<<endl;
      cout<<"icomb_k0: "<<icomb_k0<<" th_t: "<<header_3pt.comb[icomb_k0].tht[2]<<endl;
      cout<<"icomb_symm: "<<symmetric_comb<<" th_t: "<<header_3pt.comb[symmetric_comb].tht[2]<<endl;
      cout<<"icomb_symm_k0: "<<symmetric_comb_k0<<" th_t: "<<header_3pt.comb[symmetric_comb_k0].tht[2]<<endl;
      cout<<"##############"<<endl;






     

      //print correlators in the polarization basis

      distr_t_list V_pol_1_alpha_1, V_pol_1_alpha_2, V_pol_2_alpha_1, V_pol_2_alpha_2;
      distr_t_list A_pol_1_alpha_1, A_pol_1_alpha_2, A_pol_2_alpha_1, A_pol_2_alpha_2;

      corr.Perform_Nt_t_average = 1 ;
      corr.Reflection_sign= -1;

      Vfloat ph_exp;
      double Eg = mom3_l.mom[i_ens][icomb3pt].Egamma();
      double EgT= sinh(Eg)*(1.0 - exp(Nt*Eg));
      for(int t=0; t<Nt;t++) ph_exp.push_back( -1.0*exp( fabs((Nt/2.0-t))*Eg));

      V_pol_1_alpha_1 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, 1, 1, "V", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, 1, 2, "V", 0), ""));
      V_pol_2_alpha_1 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, 1, 1, "V", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, 1, 2, "V", 0), ""));
      V_pol_1_alpha_2 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, 2, 1, "V", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1, icomb3pt, 2, 2, "V", 0), ""));
      V_pol_2_alpha_2 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1, icomb3pt, 2, 1, "V", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1, icomb3pt, 2, 2, "V", 0), ""));

      corr.Reflection_sign=1;
      A_pol_1_alpha_1 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, 1, 1, "A", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0, icomb3pt, 1, 2, "A", 0), ""));
      A_pol_2_alpha_1 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0, icomb3pt, 1, 1, "A", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0, icomb3pt, 1, 2, "A", 0), ""));
      A_pol_1_alpha_2 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, 2, 1, "A", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, 2, 2, "A", 0), ""));
      A_pol_2_alpha_2 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, 2, 1, "A", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,0, icomb3pt, 2, 2, "A", 0), ""));

      corr.Reflection_sign=-1;

      distr_t_list symm_V_pol_1_alpha_1 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, 1, 1, "V", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, 1, 2, "V", 0), ""));
      distr_t_list symm_V_pol_2_alpha_1 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,symmetric_comb, 1, 1, "V", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,  symmetric_comb, 1, 2, "V", 0), ""));
      distr_t_list symm_V_pol_1_alpha_2 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,symmetric_comb, 2, 1, "V", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1,  symmetric_comb, 2, 2, "V", 0), ""));
      distr_t_list symm_V_pol_2_alpha_2 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1,symmetric_comb, 2, 1, "V", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  1, symmetric_comb, 2, 2, "V", 0), ""));

      corr.Reflection_sign=1;
      distr_t_list symm_A_pol_1_alpha_1 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0,symmetric_comb, 1, 1, "A", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0,  symmetric_comb, 1, 2, "A", 0), ""));
      distr_t_list symm_A_pol_2_alpha_1 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0,symmetric_comb, 1, 1, "A", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,  0,  symmetric_comb, 1, 2, "A", 0), ""));
      distr_t_list  symm_A_pol_1_alpha_2 = (1.0/sqrt(2))*(-1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0,symmetric_comb, 2, 1, "A", 0), "") -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, symmetric_comb, 2, 2, "A", 0), ""));
      distr_t_list symm_A_pol_2_alpha_2 = (1.0/sqrt(2))*(corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0,symmetric_comb, 2, 1, "A", 0), "")  -1.0*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,0,  symmetric_comb, 2, 2, "A", 0), ""));

      
      distr_t_list sum_V = ph_exp*( -1.0*( 1.0*V_pol_1_alpha_1/4.0 + +1.0*V_pol_2_alpha_1/4.0 -1.0*V_pol_1_alpha_2/4.0 + V_pol_2_alpha_2/4.0));
      distr_t_list sum_A = ph_exp*( -1.0*A_pol_1_alpha_1/4.0 +1.0*A_pol_2_alpha_1/4.0 -1.0*A_pol_1_alpha_2/4.0 - 1.0*A_pol_2_alpha_2/4.0);
      distr_t_list symm_sum_V = ph_exp*( -1.0*( 1.0*symm_V_pol_1_alpha_1/4.0 + +1.0*symm_V_pol_2_alpha_1/4.0 -1.0*symm_V_pol_1_alpha_2/4.0 + symm_V_pol_2_alpha_2/4.0));
      distr_t_list symm_sum_A = ph_exp*( -1.0*symm_A_pol_1_alpha_1/4.0 +1.0*symm_A_pol_2_alpha_1/4.0 -1.0*symm_A_pol_1_alpha_2/4.0 - 1.0*symm_A_pol_2_alpha_2/4.0);
      distr_t_list sum_V_double_ins, sum_A_double_ins;
      sum_V_double_ins = (2.0/3.0)*sum_V -(1.0/3.0)*symm_sum_V;
      sum_A_double_ins = (2.0/3.0)*sum_A +(1.0/3.0)*symm_sum_A; 
      //if( Eg != 0) { sum_V_alpha_1 = sum_V_alpha_1/(2.0*EgT); sum_A_alpha_1/(2.0*EgT);}

      Print_To_File({}, {sum_V.ave(), sum_V.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/summ_V_icomb_"+to_string(icomb3pt)+".dat","", "");
      Print_To_File({}, {sum_A.ave(), sum_A.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/summ_A_icomb_"+to_string(icomb3pt)+".dat","", "");
      Print_To_File({}, {sum_A_double_ins.ave(), sum_A_double_ins.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/ins_A_icomb_"+to_string(icomb3pt)+".dat","", "");
      Print_To_File({}, {sum_V_double_ins.ave(), sum_V_double_ins.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/ins_V_icomb_"+to_string(icomb3pt)+".dat","", "");

      

      corr.Perform_Nt_t_average = 1;

      Print_To_File({}, {V_pol_1_alpha_1.ave(), V_pol_1_alpha_1.err(), V_pol_1_alpha_2.ave(), V_pol_1_alpha_2.err(), V_pol_2_alpha_1.ave(), V_pol_2_alpha_1.err(), V_pol_2_alpha_2.ave(), V_pol_2_alpha_2.err() }, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/polariz_basis_V_"+to_string(icomb3pt)+".dat","", "# pol_1_alpha1               pol_1_alpha2           pol_2_alpha1        pol_2_alpha_2");
       Print_To_File({}, {A_pol_1_alpha_1.ave(), A_pol_1_alpha_1.err(), A_pol_1_alpha_2.ave(), A_pol_1_alpha_2.err(), A_pol_2_alpha_1.ave(), A_pol_2_alpha_1.err(), A_pol_2_alpha_2.ave(), A_pol_2_alpha_2.err() }, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/polariz_basis_A_"+to_string(icomb3pt)+".dat","", "# pol_1_alpha1               pol_1_alpha2           pol_2_alpha1        pol_2_alpha_2");


        //loop over alpha and mu
        
      for(int alpha=0; alpha<4;alpha++) {

	
	for(int mu=0; mu<4;mu++) {

	  //double e_f1 =2.0/3.0;
	  //double e_f2 =-1.0/3.0;

	  double e_f1 = 2.0/3.0;
	  double e_f2 = -1.0/3.0;



	  corr.Reflection_sign = -1;
   
        
	  
	  distr_V[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, alpha, mu, "V", sm_lev), "");
	  
	  distr_V_k0[alpha][mu] =e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,icomb_k0, alpha, mu, "V", sm_lev), "");

	  
	  corr.Reflection_sign = 1;
	  
	  distr_A[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, alpha, mu, "A", sm_lev), "");
	  
	  distr_A_k0[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,0, icomb_k0, alpha, mu, "A", sm_lev), "");



	  
	  Print_To_File({}, {distr_A[alpha][mu].ave(), distr_A[alpha][mu].err(), distr_V[alpha][mu].ave(), distr_V[alpha][mu].err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/no_symm"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");




	  
	  corr.Reflection_sign= -1;
	  distr_V[alpha][mu] = distr_V[alpha][mu] + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, alpha, mu, "V", sm_lev), "");
	  distr_V_k0[alpha][mu] = distr_V_k0[alpha][mu]  + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb_k0, alpha, mu, "V", sm_lev), "");
	  corr.Reflection_sign = 1;
	  distr_A[alpha][mu] = distr_A[alpha][mu]  - e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, symmetric_comb, alpha, mu, "A", sm_lev), "");
	  distr_A_k0[alpha][mu] = distr_A_k0[alpha][mu]   - e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, symmetric_comb_k0, alpha, mu, "A", sm_lev), "");

	  auto F_min_new = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};
	  int pt2_k0p0_new = Get_2pt_k0p0(header_2pt, c.mu1, c.mu2);
	  int pt2_p_new = Get_2pt_p(header_2pt, c.i0, c.is);
	  distr_t_list E_meson_new = distr_t_list::f_of_distr(F_min_new, m_fit_distr[i_ens][pt2_p_new], Nt);
	  Vfloat eeexp;
	  double Egamma_new = mom3_l.mom[i_ens][icomb3pt].Egamma();
	  double xg_new = mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0_new]).ave();
	  for(int t=0; t<Nt;t++) eeexp.push_back( exp( fabs((Nt/2-t))*Egamma_new));

	  distr_t_list distr_V_exp = distr_V[alpha][mu]*E_meson_new*eeexp;

	  distr_t_list distr_A_exp = distr_A[alpha][mu]*E_meson_new*eeexp;

	  if(alpha==mu) distr_corr.push_back(distr_A_exp);

	  if(xg_new < -1 && (c.mu1 > c.mu2) && (alpha==mu)) {


	    cout<<"printing correlator"<<endl;
	    for(int tt=0;tt<Nt;tt++) cout<<distr_A_exp.ave()[tt]<<" "<<distr_A_exp.err()[tt]<<endl;
      
	    cout<<"fitting contaminations for xg: "<<xg_new<<endl;
	    int nboot_exp=1;
	    bootstrap_fit<ExpFit, Val> fit(nboot_exp);
	    fit.Set_number_of_measurements(Nt/2 - 12);
	    fit.Set_verbosity(1);
	    /*  if(Nt==64) {
	      fit.Add_par("a1", -0.001,0.001/10); }
	    else fit.Add_par("a1", -0.003, 0.005);
	    fit.Add_par("a2", 0.03,0.01);
	    fit.Add_par("a3", 0.3,(0.3/10));
	    fit.Add_par("a4", 0.03,0.01);
	    fit.Add_par("a5", 0.4,(0.4/10));  */

	 

	    if(alpha==0) {
	      if(xg_new > 0.098 && xg_new < 0.11) {fit.Add_par("a1", -0.6,0.1);}
	      else if(xg_new > 0.198 && xg_new < 0.21) {fit.Add_par("a1", -0.4,0.1); }
	      else if(xg_new > 0.398 && xg_new < 0.41) {fit.Add_par("a1",-0.2, 0.1);}
	      else if(xg_new > 0.598 && xg_new < 0.61) {fit.Add_par("a1", -0.1,0.05);}
	      else if(xg_new > 0.798 && xg_new < 0.81) {fit.Add_par("a1", -0.05,0.05);}
	      else if(xg_new > 0.99 && xg_new < 1.1) {fit.Add_par("a1", -0.05, 0.05);}
	    }

	    else  {

	       fit.Add_par("a1", -0.05,0.005);
	    }
	    fit.Add_par("a2", -0.04,0.01);
	    fit.Add_par("a3", 0.45,0.03);
	    
	    if(Nt==64) { 
	      fit.Add_par("a4", 0.02, 0.01);
	    }
	    if(Nt==48) {
	      fit.Add_par("a4", 0.01, 0.001);
	    }
	    fit.Add_par("a5", 0.4,0.1);
	    fit.Add_par("a6", 0.1, 0.01);
	    fit.Add_par("a7", 0.1,0.01);




	    

	    

    
	    //fit.Fix_n_release("a4",0);

	    auto exp_sub = [&](double gap, double ti) ->double {
	      double res=0;
	      for(int n=0; n<= 20; n++) res += pow(-1*ti*gap,n)/(fact(n)*(n+1));

	      return res*ti;
	    };

	    fit.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
	    fit.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
	    fit.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) +par.a6*exp_sub(par.a7,val.t)+  par.a4*exp(-1*par.a5*(Nt/2-val.t));};

	    int a=0;
	    fit.ib=&a;
	    for(int i=2;i<Nt/2 -10;i++) {
	      double t= (double)i;
	      double meas= distr_A_exp.ave()[i];
	      double err = distr_A_exp.err()[i];
	      cout<<"meas: "<<meas<<" "<<err<<endl;
	      Val X;
	      X.t= t;
	      X.meas=meas;
	      X.err=err;
	      fit.Append_to_input_par(X);
	    }

	    fit.Fix_par("a4",0.0);
	    fit.Fix_par("a5",0.0);
	    fit.Fix_par("a6",0.0);
	    fit.Fix_par("a7",0.0);
	    boot_fit_data<ExpFit> Fit_output1 = fit.Perform_bootstrap_fit();


	    bootstrap_fit<ExpFit, Val> fit2(nboot_exp);
	    fit2.Set_number_of_measurements(Nt/2-2 - (Nt/2 - 10));
	    fit2.Set_verbosity(1);
	    /*  if(Nt==64) {
	      fit2.Add_par("a1", -0.001,0.001/10); }
	    else fit2.Add_par("a1", -0.003, 0.005);
	    fit2.Add_par("a2", 0.03,0.01);
	    fit2.Add_par("a3", 0.1,(0.3/10));
	    fit2.Add_par("a4", 0.03,0.01);
	    fit2.Add_par("a5", 0.2,(0.4/10));  */

	    if(alpha==0) {
	      if(xg_new > 0.098 && xg_new < 0.11) {fit2.Add_par("a1", -0.6,0.1);}
	      else if(xg_new > 0.198 && xg_new < 0.21) {fit2.Add_par("a1", -0.4,0.1); }
	      else if(xg_new > 0.398 && xg_new < 0.41) {fit2.Add_par("a1",-0.2, 0.1);}
	      else if(xg_new > 0.598 && xg_new < 0.61) {fit2.Add_par("a1", -0.1,0.05);}
	      else if(xg_new > 0.798 && xg_new < 0.81) {fit2.Add_par("a1", -0.05,0.05);}
	      else if(xg_new > 0.99 && xg_new < 1.1) {fit2.Add_par("a1", -0.05, 0.05);}
	    }
	    
	    else  {

	       fit2.Add_par("a1", -0.05,0.005);
	    }
	      fit2.Add_par("a2", -0.2,0.1);
	      fit2.Add_par("a3", 0.4,0.1);
	      if(Nt==64) {
	      fit2.Add_par("a4", 0.02, 0.01);
	      }
	      if(Nt==48) {
	      fit2.Add_par("a4", 0.1, 0.01);
	      }
	      fit2.Add_par("a5", 0.4,0.1);
	    fit2.Add_par("a6",0.1,0.01);
	    fit2.Add_par("a7",0.001,0.0001);

	    fit2.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
	    fit2.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
	    fit2.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) + par.a6*exp_sub(par.a7,val.t)+ par.a4*exp(-1*par.a5*(Nt/2-val.t));};


	    int b=0;
	    fit2.ib=&b;

	    for(int i=Nt/2-10;i<Nt/2 - 2;i++) {
	      double t= (double)i;
	      double meas= distr_A_exp.ave()[i];
	      double err = distr_A_exp.err()[i];
	      Val X;
	      X.t= t;
	      X.meas=meas;
	      X.err=err;
	      fit2.Append_to_input_par(X);
	    }

	    fit2.Fix_par("a2",0.0);
	    fit2.Fix_par("a3",0.0);
	    fit2.Fix_par("a6",0.0);
	    fit2.Fix_par("a7",0.0);
	    boot_fit_data<ExpFit> Fit_output2 = fit2.Perform_bootstrap_fit();

	


	    //total fit

	    nboot_exp=100;
	    bootstrap_fit<ExpFit, Val> fit_tot(nboot_exp);
	    fit_tot.Set_number_of_measurements(Nt/2 - 4);
	    fit_tot.Set_verbosity(1);

	    fit_tot.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
	    fit_tot.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
	    fit_tot.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) +  par.a6*exp_sub(par.a7,val.t) + par.a4*exp(-1*par.a5*(Nt/2-val.t));};

	    fit_tot.Add_par("a1", Fit_output1.par[0].a1, Fit_output1.par[0].a1/100);

           
	    fit_tot.Add_par("a2", Fit_output1.par[0].a2,Fit_output1.par[0].a2/100);
	    fit_tot.Add_par("a3", Fit_output1.par[0].a3,Fit_output1.par[0].a3/100);
	    fit_tot.Add_par("a4", Fit_output2.par[0].a4,Fit_output2.par[0].a4/100);
	    fit_tot.Add_par("a5", Fit_output2.par[0].a5,Fit_output2.par[0].a5/100);
	    //fit_tot.Add_par("a6", Fit_output1.par[0].a6, Fit_output1.par[0].a6/100);
	    //fit_tot.Add_par("a7", Fit_output1.par[0].a7, Fit_output1.par[0].a7/100);
	    fit_tot.Add_par("a6", 0.001, 0.001);
	    fit_tot.Add_par("a7",1.5*Egamma_new,1.5*Egamma_new/10);

	    if(xg_new >= 0.0) { fit_tot.Fix_par("a6",0.0); fit_tot.Fix_par("a7",0.0);}

             
	    GaussianMersenne GM(645645);
	    for(int iboot=0;iboot<nboot_exp;iboot++) {
	      fit_tot.ib =&iboot;
	      for(int i=2;i<Nt/2 -2;i++) {
		double t= (double)i;
		double meas= distr_A_exp.ave()[i];
		double err = distr_A_exp.err()[i];
		Val X;
		X.t= t;
		X.meas=meas+ GM()*err/10.0 ;
		X.err=err;
		fit_tot.Append_to_input_par(X);
	      }
	    }

          

      
	    boot_fit_data<ExpFit> Fit_output= fit_tot.Perform_bootstrap_fit();
	    //compute average
	    double ave_ch2=0;
	    for(auto &c : Fit_output.chi2) ave_ch2+=c/nboot_exp;
	    Vfloat a1_ave;
	    Vfloat a2_ave;
	    Vfloat a4_ave;
	    Vfloat a3_ave;
	    Vfloat a5_ave;
	    Vfloat a6_ave;
	    Vfloat a7_ave;
	    for(int iboot=0;iboot<nboot_exp;iboot++) {a1_ave.push_back(Fit_output.par[iboot].a1); a2_ave.push_back(Fit_output.par[iboot].a2); a3_ave.push_back(Fit_output.par[iboot].a3); a4_ave.push_back(Fit_output.par[iboot].a4); a5_ave.push_back(Fit_output.par[iboot].a5); a6_ave.push_back(Fit_output.par[iboot].a6); a7_ave.push_back(Fit_output.par[iboot].a7); }

	    cout<<"alpha= "<<alpha<<"  "<<"mu= "<<mu<<" "<<"xg: "<<xg_new<<endl;
	    cout<<"printing final values of params"<<endl;
	    cout<<"C_\alpha^{mu}: "<<xg_new<<"   ENS: "<<Ens_tag[i_ens]<<endl;
	    cout<<"a1: "<<Boot_ave(a1_ave)<<" "<<10*Boot_err(a1_ave)<<endl;
	    cout<<"a2: "<<Boot_ave(a2_ave)<<" "<<10*Boot_err(a2_ave)<<endl;
	    cout<<"a3: "<<Boot_ave(a3_ave)<<" "<<10*Boot_err(a3_ave)<<endl;
	    cout<<"a4: "<<Boot_ave(a4_ave)<<" "<<10*Boot_err(a4_ave)<<endl;
	    cout<<"a5: "<<Boot_ave(a5_ave)<<" "<<10*Boot_err(a5_ave)<<endl;
	    cout<<"a6: "<<Boot_ave(a6_ave)<<" "<<10*Boot_err(a6_ave)<<endl;
	    cout<<"a7: "<<Boot_ave(a7_ave)<<" "<<10*Boot_err(a7_ave)<<endl;
	    cout<<"ch2: "<<ave_ch2<<"/"<<Nt/2 -4 - 5  <<endl;

     
	  }

	 

	  Print_To_File({}, {distr_A_exp.ave(), distr_A_exp.err(), distr_V_exp.ave(), distr_V_exp.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/symm_comb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");


	 
	}
      }

      
      int pt2_k0p0 = Get_2pt_k0p0(header_2pt, c.mu1, c.mu2);
      double xg3 = mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0]).ave();
      cout<<"printing ratio corr:"<<endl;
      cout<<"Ens: "<<Ens_tag[i_ens]<<"icomb: "<<icomb3pt<<"  xg: "<<xg3<<endl;
      cout<<"alpha=1"<<endl;
      if(xg3 > 0.005 && Ens_tag[i_ens].substr(0,1)=="E") {
      distr_t_list ratio_1 = distr_corr[1]/distr_corr[2];
      distr_t_list ratio_2 = distr_corr[3]/distr_corr[2];
      for(int i=0;i<Nt;i++) cout<<i<<" "<<ratio_1.ave()[i]<<"   "<<ratio_1.err()[i]<<endl;
      cout<<"alpha=3"<<endl;
      for(int i=0;i<Nt;i++) cout<<i<<" "<<ratio_2.ave()[i]<<"   "<<ratio_2.err()[i]<<endl;
      }


      //get \bar{R}_{A/V}

      double Egamma = mom3_l.mom[i_ens][icomb3pt].Egamma();
      cout<<"Egamma["<<icomb3pt<<"] :"<<Egamma<<endl;
     
      int pt2_p = Get_2pt_p(header_2pt, c.i0, c.is);
      cout<<"Mass: "<<m_fit_distr[i_ens][pt2_k0p0].ave()<<endl;
      cout<<"x_gamma: "<<mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_p]).ave()<<endl;
     
      
      int twall= ( (Ens_tag[i_ens].substr(0,1) == "T") || (Ens_tag[i_ens].substr(0,1) == "L"))?20:Nt/2;
      distr_t_list V_ave = V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[i_ens][pt2_p], mom3_l.mom[i_ens][icomb3pt], twall);
      distr_t_list A_ave = A_ave_unpolarized(distr_A);
      distr_t_list A_ave_zero_mom = A_ave_unpolarized(distr_A_k0);
      distr_t_list HV = H_V(distr_V, distr_V_k0, mom3_l.mom[i_ens][icomb3pt]);

      
      Vfloat exp_Nt_t;
      auto F_min = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};
      distr_t_list E_meson = distr_t_list::f_of_distr(F_min, m_fit_distr[i_ens][pt2_p], Nt);
      for(int t=0; t<Nt;t++) exp_Nt_t.push_back( exp( fabs((twall-t))*Egamma));
      distr_t F_A_distr_factor= (2.0*fp_fit_distr[i_ens][pt2_k0p0]/(m_fit_distr[i_ens][pt2_k0p0]*mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0])));


      auto sq= [](double x)->double {return sqrt(x);};

      distr_t matrix_el_PS_smeared = overlap_smeared_fit_distr[i_ens][pt2_p];
      distr_t matrix_el_PS_smeared_sqrt = distr_t::f_of_distr(sq, overlap_smeared_fit_distr[i_ens][pt2_p]);
      distr_t matrix_el_PS = distr_t::f_of_distr(sq, overlap_fit_distr[i_ens][pt2_p]);
      distr_t matrix_el_PS_smeared_new = matrix_el_PS_smeared;
      matrix_el_PS_smeared = matrix_el_PS_smeared/matrix_el_PS;
     

      cout<<"matrix_el_PS: "<<matrix_el_PS.ave()<<"  "<<matrix_el_PS.err()<<endl;
      cout<<"matrix_el_PS_smeared: "<<matrix_el_PS_smeared.ave()<<" "<<matrix_el_PS_smeared.err()<<endl;

      //for(int i=0; i<matrix_el_PS.size();i++) cout<<matrix_el_PS.distr[i]<<" "<<matrix_el_PS_smeared.distr[i]<<endl;
      // bar_R_A_distr[i_ens].push_back( (exp_Nt_t*(A_ave/A_ave_zero_mom)) -1.0);
      bar_R_A_distr[i_ens].push_back( (exp_Nt_t*(A_ave/A_ave_zero_mom))-1.0);
      // bar_R_V_distr[i_ens].push_back( (V_ave/A_ave_zero_mom)); //FRA NORMALIZATION FOR VECTOR FORM FACTOR
      //bar_R_V_distr[i_ens].push_back( (V_ave/A_ave_zero_mom)*fp_fit_distr[i_ens][pt2_k0p0]/mom3_l.mom[i_ens][icomb3pt].k()[2]);
      // bar_R_V_distr[i_ens].push_back( (V_ave/A_ave_zero_mom)*fp_fit_distr[i_ens][pt2_k0p0]*m_fit_distr[i_ens][pt2_k0p0]);
      bar_R_V_distr[i_ens].push_back( (V_ave/A_ave_zero_mom)*fp_fit_distr[i_ens][pt2_k0p0]*m_fit_distr[i_ens][pt2_k0p0]);
      //bar_R_A_no_bar_distr[i_ens].push_back(((A_ave/(matrix_el_PS*mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0])))*E_meson*exp_Nt_t) - 2.0*fp_fit_distr[i_ens][pt2_p]/(mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0])*m_fit_distr[i_ens][pt2_k0p0]));
      bar_R_A_no_bar_distr[i_ens].push_back( (A_ave/matrix_el_PS_smeared_new)*E_meson*exp_Nt_t);
      bar_R_V_no_bar_distr[i_ens].push_back(V_ave*E_meson);
      //bar_R_V_no_bar_distr[i_ens].push_back((V_ave*m_fit_distr[i_ens][pt2_p]/matrix_el_PS_smeared)*(E_meson*m_fit_distr[i_ens][pt2_k0p0]/2.0));
      HV_k[i_ens].push_back((HV*m_fit_distr[i_ens][pt2_p]/matrix_el_PS_smeared)*(E_meson*m_fit_distr[i_ens][pt2_k0p0]/2.0));
      distr_t_list F_A_distr = bar_R_A_distr[i_ens][icomb3pt]*F_A_distr_factor;
      
      VVfloat To_print({bar_R_A_distr[i_ens][icomb3pt].ave(), bar_R_A_distr[i_ens][icomb3pt].err(), bar_R_V_distr[i_ens][icomb3pt].ave(), bar_R_V_distr[i_ens][icomb3pt].err(), bar_R_A_no_bar_distr[i_ens][icomb3pt].ave(), bar_R_A_no_bar_distr[i_ens][icomb3pt].err(), bar_R_V_no_bar_distr[i_ens][icomb3pt].ave(), bar_R_V_no_bar_distr[i_ens][icomb3pt].err(), HV_k[i_ens][icomb3pt].ave(), HV_k[i_ens][icomb3pt].err()});

      Print_To_File({}, To_print, "../data/form_factors/"+Meson+"/bar_R_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       RA    RA_err     RV      RV_err    RA_no_bar      RA_no_bar_err        RV_no_bar      RV_no_bar_err");
      Print_To_File({}, {F_A_distr.ave(), F_A_distr.err()} , "../data/form_factors/"+Meson+"/F_A_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       RA    RA_err     RV      RV_err");


      double xg = mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0]).ave();
      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "V");
      distr_t fit_result_V = corr.Fit_distr(bar_R_V_distr[i_ens][icomb3pt]);
      if(c.mu1 > c.mu2) {  
      boost::filesystem::create_directory("../plots/form_factors");
      plt::clf();
      plt::xlabel("$t/a$");
      plt::ylabel("$F_{V}$");
      plt::xlim(max(0,corr.Tmin - 20), corr.Tmax+12);
      plt::grid(4);
      Vfloat TT;
      Vfloat TT_tot;
      for(int t=0;t<corr.Nt;t++) TT_tot.push_back(t);
      //generate data_point for fit;
      Vfloat FV, FV_err;
      for(int t=corr.Tmin;t<=corr.Tmax;t++) {
      FV.push_back( fit_result_V.ave());
      FV_err.push_back(fit_result_V.err());
      TT.push_back(t);
      
      }
      plt::errorbar(TT_tot, bar_R_V_distr[i_ens][icomb3pt].ave(), bar_R_V_distr[i_ens][icomb3pt].err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "xg="+to_string_with_precision(xg,4)}});
      plt::errorbar(TT, FV, FV_err, { {"c", "red"}, {"marker", "."}, {"ls", "-"}, {"label", "fitted value"}});
      // Enable legend.
      plt::legend();
      string figure_path_V= "../plots/form_factors/"+Meson+"/FV_"+Ens_tag[i_ens]+"_xg_"+to_string_with_precision(xg,4)+".png";
      plt::save(figure_path_V.c_str());
      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "A");
      distr_t fit_result_A = corr.Fit_distr(F_A_distr);
      plt::clf();
      TT.clear();
      TT_tot.clear();
      plt::xlabel("$t/a$");
      plt::ylabel("$F_{A}$");
      plt::grid(4);
      plt::xlim(max(corr.Tmin - 20,0), corr.Tmax+12);
      plt::ylim(-0.05, 0.09);
      for(int t= 0;t<corr.Nt;t++) TT_tot.push_back(t);
      //generate data_point for fit;
      Vfloat FA, FA_err;
      for(int t=corr.Tmin;t<=corr.Tmax;t++) {
	FA.push_back(fit_result_A.ave());
	FA_err.push_back(fit_result_A.err());
	TT.push_back(t);
      }
      plt::errorbar(TT_tot, F_A_distr.ave(), F_A_distr.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "xg="+to_string_with_precision(xg,4)}});
      plt::errorbar(TT, FA, FA_err, { {"c", "red"}, {"marker", "."}, {"ls", "-"}, {"label", "fitted value"}});
    
    

      // Enable legend.
      plt::legend();
      string figure_path_A= "../plots/form_factors/"+Meson+"/FA_"+Ens_tag[i_ens]+"_xg_"+to_string_with_precision(xg,4)+".png";
      plt::save(figure_path_A.c_str());
      cout<<"IRUN:"<<Ens_tag[i_ens]<<" ICOMB:"<<icomb3pt<<" "<<mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_p]).ave()<<"\t"<<fit_result_V.ave()<<"\t"<<fit_result_V.err()<<"\t"<<fit_result_A.ave()<<"\t"<<fit_result_A.err()<<endl;
      

    //plot

      if(xg > 0.01) {
      
      cout<<"fitting contaminations for xg: "<<xg<<endl;
      int nboot_exp=1;
      bootstrap_fit<ExpFit, Val> fit(nboot_exp);
      fit.Set_number_of_measurements(Nt/2 - 16);
      fit.Set_verbosity(1);
      /*     if(Nt==64) {
	fit.Add_par("a1", -0.7,0.7/10); }
      else fit.Add_par("a1", -0.7, 0.7/10);
      fit.Add_par("a2", 0.3,0.3);
      fit.Add_par("a3", 0.3,(0.3/10));
      fit.Add_par("a4", 0.03,0.03);
      fit.Add_par("a5", 0.4,(0.4/10));  */


      //bar_R_A
      /*

      

       fit.Add_par("a1", 0.02,0.02);
       if(Ens_tag[i_ens].substr(0,1)=="D") fit.Add_par("a2", 0.12,0.02);
       else fit.Add_par("a2", 0.12,0.02);
       fit.Add_par("a3", 0.2,0.1);
       
       if(Nt==64) { 
      	fit.Add_par("a4", -0.2, 0.1);
       }
      if(Nt==48) {
      	fit.Add_par("a4", -0.2, 0.1);
       }
       fit.Add_par("a5", 0.4,0.2); 
      fit.Add_par("a6", 0.1, 0.01);
      fit.Add_par("a7", 0.1,0.01); 
      */

      //C_A

      /*
          fit.Add_par("a1", -1.17,0.03);
       if(Ens_tag[i_ens].substr(0,1)=="D") fit.Add_par("a2", -0.2,0.02);
       else fit.Add_par("a2", -0.2,0.02);
       fit.Add_par("a3", 0.3,0.05);
       
       if(Nt==64) { 
      	fit.Add_par("a4", 0.4, 0.1);
       }
      if(Nt==48) {
      	fit.Add_par("a4", 0.4, 0.1);
       }
       fit.Add_par("a5", 0.4,0.1); 
      fit.Add_par("a6", 0.1, 0.01);
      fit.Add_par("a7", 0.1,0.01);

      */
      


      
      //bar_R_V
      
      
      fit.Add_par("a1", fit_result_A.ave() ,fit_result_A.err());

        
      if(Ens_tag[i_ens].substr(0,1)=="D") fit.Add_par("a2", 0.1,0.03);
      else fit.Add_par("a2",0.1,0.03);
      fit.Add_par("a3", 0.3,0.1);
      fit.Add_par("a4", -0.1,0.03);
      fit.Add_par("a5", 0.3,0.2); 
      fit.Add_par("a6",0.1,0.01);
      fit.Add_par("a7",0.1,0.0001);

    
      //fit.Fix_n_release("a4",0);

      auto exp_sub = [&](double gap, double ti) ->double {
	double res=0;
	for(int n=0; n<= 20; n++) res += pow(-1*ti*gap,n)/(fact(n)*(n+1));

	return res*ti;
      };

      fit.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
      fit.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
      fit.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) +par.a6*exp_sub(par.a7,val.t)+  par.a4*exp(-1*par.a5*(Nt/2-val.t));};

      int a=0;
      fit.ib=&a;
      	for(int i=2;i<Nt/2 -14;i++) {
	double t= (double)i;
	double meas=  F_A_distr.ave()[i];
	double err =  F_A_distr.err()[i];
	Val X;
	X.t= t;
	X.meas=meas;
	X.err=err;
	fit.Append_to_input_par(X);
	}

	fit.Fix_par("a4",0.0);
	fit.Fix_par("a5",0.0);
	fit.Fix_par("a6",0.0);
	fit.Fix_par("a7",0.0);
	boot_fit_data<ExpFit> Fit_output1 = fit.Perform_bootstrap_fit();


      bootstrap_fit<ExpFit, Val> fit2(nboot_exp);
      fit2.Set_number_of_measurements(Nt/2-2 - (Nt/2 - 12));
      fit2.Set_verbosity(1);


      /*
      //bar_R_A
      
      fit2.Add_par("a1", 0.02,0.02);
      if(Ens_tag[i_ens].substr(0,1)=="D") fit2.Add_par("a2", 0.1,0.01);
      else fit2.Add_par("a2", 0.1,0.01);
      fit2.Add_par("a3", 0.3,0.1);
      fit2.Add_par("a4", -0.2,0.1);
      fit2.Add_par("a5", 0.4,(0.4/10));  */

      //C_A

      /*
      
      fit2.Add_par("a1", -1.17,0.03);
       if(Ens_tag[i_ens].substr(0,1)=="D") fit2.Add_par("a2", -0.2,0.1);
       else fit2.Add_par("a2", -0.2,0.1);
       fit2.Add_par("a3", 0.3,0.1);
       
       if(Nt==64) { 
      	fit2.Add_par("a4", 0.4, 0.1);
       }
      if(Nt==48) {
      	fit2.Add_par("a4", 0.4, 0.1);
       }
      fit2.Add_par("a5", 0.4,0.1); 
      fit2.Add_par("a6", 0.1, 0.01);
      fit2.Add_par("a7", 0.1,0.01);

      */

      // bar_R_V
       
     
      fit2.Add_par("a1",  fit_result_A.ave(), fit_result_A.err());
        

      
      if(Ens_tag[i_ens].substr(0,1)=="D") fit2.Add_par("a2", 0.1,0.03);
      else fit2.Add_par("a2",0.1,0.03);
      fit2.Add_par("a3", 0.3,0.1);
      fit2.Add_par("a4", -0.1,0.03);
      fit2.Add_par("a5", 0.3,0.2); 
      fit2.Add_par("a6",0.1,0.01);
      fit2.Add_par("a7",0.1,0.0001);

      /* fit2.Add_par("a1", -0.2,0.05);
      fit2.Add_par("a2", -0.2,0.1);
      fit2.Add_par("a3", 0.4,0.1);
      if(Nt==64) {
      	fit2.Add_par("a4", 0.02, 0.01);
      }
      if(Nt==48) {
      	fit2.Add_par("a4", 0.1, 0.01);
      }
      fit2.Add_par("a5", 0.4,0.1); */
      fit2.Add_par("a6",0.1,0.01);
      fit2.Add_par("a7",0.001,0.0001);

      fit2.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
      fit2.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
      fit2.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) + par.a6*exp_sub(par.a7,val.t)+ par.a4*exp(-1*par.a5*(Nt/2-val.t));};


      int b=0;
      fit2.ib=&b;

      	for(int i=Nt/2-12;i<Nt/2 - 2;i++) {
	double t= (double)i;
	double meas=  F_A_distr.ave()[i];
	double err =  F_A_distr.err()[i];
	Val X;
	X.t= t;
	X.meas=meas;
	X.err=err;
	fit2.Append_to_input_par(X);
	}

	fit2.Fix_par("a2",0.0);
	fit2.Fix_par("a3",0.0);
	fit2.Fix_par("a6",0.0);
	fit2.Fix_par("a7",0.0);
	boot_fit_data<ExpFit> Fit_output2 = fit2.Perform_bootstrap_fit();

	


	//total fit

      nboot_exp=100;
      bootstrap_fit<ExpFit, Val> fit_tot(nboot_exp);
      fit_tot.Set_number_of_measurements(Nt/2 - 4);
      fit_tot.Set_verbosity(1);

      fit_tot.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
      fit_tot.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
      fit_tot.ansatz = [=](const ExpFit& par, const Val& val) -> double { return par.a1 + par.a2*exp(-1*par.a3*val.t) +  par.a6*exp_sub(par.a7,val.t) + par.a4*exp(-1*par.a5*(Nt/2-val.t));};

      fit_tot.Add_par("a1", Fit_output1.par[0].a1, Fit_output1.par[0].a1/100);

           
      fit_tot.Add_par("a2", Fit_output1.par[0].a2,Fit_output1.par[0].a2/100);
      fit_tot.Add_par("a3", Fit_output1.par[0].a3,Fit_output1.par[0].a3/100);
      fit_tot.Add_par("a4", Fit_output2.par[0].a4,Fit_output2.par[0].a4/100);
      fit_tot.Add_par("a5", Fit_output2.par[0].a5,Fit_output2.par[0].a5/100);
      //fit_tot.Add_par("a6", Fit_output1.par[0].a6, Fit_output1.par[0].a6/100);
      //fit_tot.Add_par("a7", Fit_output1.par[0].a7, Fit_output1.par[0].a7/100);
      fit_tot.Add_par("a6", 0.001, 0.001);
      fit_tot.Add_par("a7",1.5*Egamma,1.5*Egamma/10);

      if(xg >= 0.0) { fit_tot.Fix_par("a6",0.0); fit_tot.Fix_par("a7",0.0);}

             
      GaussianMersenne GM(645645);
      for(int iboot=0;iboot<nboot_exp;iboot++) {
	fit_tot.ib =&iboot;
	for(int i=2;i<Nt/2 -2;i++) {
	double t= (double)i;
	double meas=  F_A_distr.ave()[i];
	double err =  F_A_distr.err()[i];
	Val X;
	X.t= t;
	X.meas=meas+ GM()*err/10.0 ;
	X.err=err;
	fit_tot.Append_to_input_par(X);
      }
      }

          

      
      boot_fit_data<ExpFit> Fit_output= fit_tot.Perform_bootstrap_fit();
      //compute average
      double ave_ch2=0;
      for(auto &c : Fit_output.chi2) ave_ch2+=c/nboot_exp;
      Vfloat a1_ave;
      Vfloat a2_ave;
      Vfloat a4_ave;
      Vfloat a3_ave;
      Vfloat a5_ave;
      Vfloat a6_ave;
      Vfloat a7_ave;
      for(int iboot=5;iboot<nboot_exp;iboot++) {a1_ave.push_back(Fit_output.par[iboot].a1); a2_ave.push_back(Fit_output.par[iboot].a2); a3_ave.push_back(Fit_output.par[iboot].a3); a4_ave.push_back(Fit_output.par[iboot].a4); a5_ave.push_back(Fit_output.par[iboot].a5); a6_ave.push_back(Fit_output.par[iboot].a6); a7_ave.push_back(Fit_output.par[iboot].a7); }
      cout<<"printing final values of params"<<endl;
      cout<<"XG: "<<xg<<"   ENS: "<<Ens_tag[i_ens]<<endl;
      cout<<"a1: "<<Boot_ave(a1_ave)<<" "<<10*Boot_err(a1_ave)<<endl;
      cout<<"a2: "<<Boot_ave(a2_ave)<<" "<<10*Boot_err(a2_ave)<<endl;
      cout<<"a3: "<<Boot_ave(a3_ave)<<" "<<10*Boot_err(a3_ave)<<endl;
      cout<<"a4: "<<Boot_ave(a4_ave)<<" "<<10*Boot_err(a4_ave)<<endl;
      cout<<"a5: "<<Boot_ave(a5_ave)<<" "<<10*Boot_err(a5_ave)<<endl;
      cout<<"a6: "<<Boot_ave(a6_ave)<<" "<<10*Boot_err(a6_ave)<<endl;
      cout<<"a7: "<<Boot_ave(a7_ave)<<" "<<10*Boot_err(a7_ave)<<endl;
      cout<<"ch2: "<<ave_ch2<<"/"<<Nt/2 -4 - 5  <<endl;

     
      }
      
      }


    }
    
    
    
    fclose(stream_3pt);
    fclose(stream_2pt);
    
  }

  return;
}

