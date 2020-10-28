#include "../include/3pts_mes_gamma_W.h"



using namespace std;

constexpr double kappa=2.837297;
const double M2PiPhys=pow(0.139,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const double r0 = pow(0.672/0.197,2);
const int Nbranches = 8;
const int Nboots= 100;
const bool UseJack=1;
const int nboots=100;
const int Njacks=40;
const int sm_lev=0;

void Get_Tmin_Tmax(string corr_type, string Ens_tag, CorrAnalysis &corr) {

  if(corr_type=="2pt") {
    corr.Tmax = corr.Nt/2 - 4;
    if(Ens_tag.substr(0,1) =="A") corr.Tmin = 15;
    else if(Ens_tag.substr(0,1) =="B") corr.Tmin = 16;
    else corr.Tmin = 19;
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


distr_t_list V_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, distr_t& Meson_mass, pt3_momenta& Mom) {

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
  for(int t=0; t<Mom.Nt();t++) exp_Nt_t.push_back( exp( fabs((Mom.Nt()/2.0-t))*Mom.Egamma()));

  if(Mom.k_mod() < eps(15)) return V_ave; 
  
  
  for(int alpha=1; alpha <= 2; alpha ++) {
   
    if(distr_mom_k[alpha].size() != 4 || distr_mom_0[alpha].size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
    for(int mu=0; mu<4;mu++) {
    
      distr_t_list den1(3, Meson_mass);
      distr_t_list den2(3, Meson_mass);

	
      den1 = den1*Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3),Mom.k()),-1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3), Mom.p()), Mom.Egamma());
    
      den2 = den2*Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3),Mom.k()), -1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3), Mom.p()), Mom.Egamma());
    
      V_ave = V_ave+ e[1][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - distr_mom_0[alpha][mu])/den1[alpha-1];

      V_ave = V_ave+ e[2][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - distr_mom_0[alpha][mu])/den2[alpha-1];

     
    }
  }
  

  return V_ave;

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

  
  return A_ave;
  
}








void Compute_form_factors(string Meson) {




  string dir= "../3pt_data";
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
  vector<vector<distr_t>> bar_R_A_fit_distr(Nens);
  vector<vector<distr_t>> bar_R_V_fit_distr(Nens);
  vector<vector<distr_t_list>> fp_distr_list(Nens);
  vector<vector<distr_t_list>> m_distr_list(Nens);
  vector<vector<distr_t>> fp_fit_distr(Nens);
  vector<vector<distr_t>> m_fit_distr(Nens);
  

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


    //Forward informations to 3pt_momenta_list
    double ml, l, Nt;
    Read_pars_from_ensemble_tag(Ens_tag[i_ens], ml, l , Nt);
    Add_to_mom_list(mom3_l, header_3pt, l);
    Add_to_mom_list(mom2_l, header_2pt, l);
    mom2_l.Nt.push_back(Nt);
    mom3_l.Nt.push_back(Nt);

    //init CorrAnalysis class
    CorrAnalysis corr(UseJack, Njacks, nboots);
    corr.Nt=Nt;

    

    //2pts
    Get_Tmin_Tmax("2pt", Ens_tag[i_ens], corr);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/form_factors");
    boost::filesystem::create_directory("../data/form_factors/"+Meson);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/H");

    int ncomb2pt= header_2pt.ncomb;
;

    for(int icomb2pt=0; icomb2pt<ncomb2pt;icomb2pt++) {

      auto c = header_2pt.comb[icomb2pt];
      string tag=to_string(c.i0)+"_"+to_string(c.is)+"_"+to_string(c.it);
      fp_distr_list[i_ens].push_back( (c.mu1+c.mu2)*corr.decay_constant_t( Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, sm_lev), "../data/form_factors/"+Meson+"/fp_"+Ens_tag[i_ens]+"_"+tag));
      m_distr_list[i_ens].push_back( corr.effective_mass_t(Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, sm_lev), "../data/form_factors/"+Meson+"/meson_mass_"+Ens_tag[i_ens]+"_"+tag));
      fp_fit_distr[i_ens].push_back(corr.Fit_distr(fp_distr_list[i_ens][icomb2pt]));
      m_fit_distr[i_ens].push_back(corr.Fit_distr(m_distr_list[i_ens][icomb2pt]));
    }

    //end anaysis 2pts

    //beginning analysis 3pts

    int ncomb3pt= header_3pt.ncomb;




    //loop over kinematics
    

    for(int icomb3pt=0; icomb3pt<ncomb3pt;icomb3pt++) {

      vector<vector<distr_t_list>> distr_V(4), distr_V_k0(4), distr_A(4), distr_A_k0(4);
      for(int i=0;i<4;i++) { distr_V[i].resize(4); distr_V_k0[i].resize(4); distr_A[i].resize(4); distr_A_k0[i].resize(4);}

      auto c = header_3pt.comb[icomb3pt];
      
      int icomb_k0= Get_comb_k0(header_3pt, icomb3pt);




      //loop over alpha and mu

        
      for(int alpha=0; alpha<4;alpha++) {
	
	for(int mu=0; mu<4;mu++) {

	  double e_f1 = 2.0/3;
	  double e_f2 = -1.0/3;

	  int symmetric_comb=Get_symmetric_comb(header_3pt, icomb3pt);
	  int symmetric_comb_k0 = Get_symmetric_comb(header_3pt, icomb_k0);

	  corr.Reflection_sign = -1;
	  
	  distr_V[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, alpha, mu, "V", sm_lev), "") + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, alpha, mu, "V", sm_lev), "");
	  
	  distr_V_k0[alpha][mu] =e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,icomb_k0, alpha, mu, "V", sm_lev), "") + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb_k0, alpha, mu, "V", sm_lev), "");

	  
	  corr.Reflection_sign = 1;
	  
	  distr_A[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, icomb3pt, alpha, mu, "A", sm_lev), "") - e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, symmetric_comb, alpha, mu, "A", sm_lev), "");
	  
	  distr_A_k0[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,0, icomb_k0, alpha, mu, "A", sm_lev), "") - e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 0, symmetric_comb_k0, alpha, mu, "A", sm_lev), "");

	  
	  Print_To_File({}, {distr_A[alpha][mu].ave(), distr_A[alpha][mu].err()}, "../data/form_factors/"+Meson+"/H/comb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu),"", "#");

	}
      }


      //get \bar{R}_{A/V}

      double Egamma = mom3_l.mom[i_ens][icomb3pt].Egamma();
      int pt2_k0p0 = Get_2pt_k0p0(header_2pt, c.mu1, c.mu2);
      int pt2_p = Get_2pt_p(header_2pt, c.i0, c.is);

      distr_t_list V_ave = V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[i_ens][pt2_p], mom3_l.mom[i_ens][icomb3pt]);
      distr_t_list A_ave = A_ave_unpolarized(distr_A);
      distr_t_list A_ave_zero_mom = A_ave_unpolarized(distr_A_k0);

      
      Vfloat exp_Nt_t;
      for(int t=0; t<Nt;t++) exp_Nt_t.push_back( exp( fabs((Nt/2.0-t))*Egamma));
      distr_t F_A_distr_factor= (2.0*fp_fit_distr[i_ens][pt2_k0p0]/(m_fit_distr[i_ens][pt2_k0p0]*mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[i_ens][pt2_k0p0])));

      
      bar_R_A_distr[i_ens].push_back( (A_ave/A_ave_zero_mom)*exp_Nt_t);
      bar_R_V_distr[i_ens].push_back( (V_ave/A_ave_zero_mom)*fp_fit_distr[i_ens][pt2_k0p0]*m_fit_distr[i_ens][pt2_k0p0]);
      distr_t_list F_A_distr = F_A_distr_factor*bar_R_A_distr[i_ens][icomb3pt];
      
      VVfloat To_print({bar_R_A_distr[i_ens][icomb3pt].ave(), bar_R_A_distr[i_ens][icomb3pt].err(), bar_R_V_distr[i_ens][icomb3pt].ave(), bar_R_V_distr[i_ens][icomb3pt].err()});

      Print_To_File({}, To_print, "../data/form_factors/"+Meson+"/bar_R_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+".dat", "", "#t       RA    RA_err     RV      RV_err");
      Print_To_File({}, {F_A_distr.ave(), F_A_distr.err()} , "../data/form_factors/"+Meson+"/F_A_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+".dat", "", "#t       RA    RA_err     RV      RV_err");
    }

    
    fclose(stream_3pt);
    fclose(stream_2pt);
  }


  return;
}
