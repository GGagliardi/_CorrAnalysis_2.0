#include "../include/Bs_mumu_gamma.h"
#include "Corr_analysis.h"
#include "numerics.h"
using namespace std;

const int Nboots= 800;
const int NJ=20;
const double Qu = -1.0/3.0; //electric charge of d-type quark
const double Qd = -1.0/3.0; //electric charge of d-type quark
const double qqu = 2.0/3.0; //electric charge of u-type quark
const double qqd = -1.0/3.0; //electric charge of d-type quark
string MESON = "Bs";
const double Lambda_QCD= 0.3; //300 MeV
const bool Is_reph=true;
Vfloat Bs_xg_t_list;
Vfloat Bs_xg_to_spline;
Vfloat Bs_xg_to_spline_VMD;
Vfloat Bs_a_to_print;
const double fT_Jpsi= 0.3927; //https://arxiv.org/pdf/2008.02024.pdf  [MSbar 2GeV]
const double fT_Jpsi_err= 0.0027; //https://arxiv.org/pdf/2008.02024.pdf [MSbar 2GeV]
const int shift=0;
const string ph_type = Is_reph ? "rph" : "vph";
const int Nmasses = 5;
const Vfloat ratio_mh({1, 1.0/1.4835, 1.0/2.02827, 1.0/2.53531, 1.0/3.04236});

 
double evolutor_ZT_MS_bar( double mu1, double mu2)  { //evolves the Z_T RC from Z_T(mu1) to Z_T(mu2) in the MSbar scheme at three loops, namely the function returns Z_T(mu2)/Z_T(mu1)



  
  double Nf= 4;
  double CF= 4.0/3;
  double CA= 3;
  double TF= 0.5;
  double zeta_3=  riemann_zeta(3);

  
  
  auto Int= [&Nf, &CF, &CA, &TF, &zeta_3 ](double mu) {

    auto alphas = [&Nf](double x) { return Get_4l_alpha_s(x, Nf);};
    
    double a =  alphas(mu)/(4.0*M_PI);

    double anomalous_dim = (4.0/3.0)*a  - 2.0*( 26*Nf - 543)*a*a/27.0 - pow(a,3)*(36*pow(Nf,2) + 1440*zeta_3*Nf + 5240*Nf + 2784*zeta_3 - 52555)/81.0;
    
   
    return anomalous_dim/mu;
    
  };


  double val, err;

 
  double prec=1e-6;
  gsl_function_pp<decltype(Int)> integrand(Int);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function *G = static_cast<gsl_function*>(&integrand);
  gsl_integration_qags(G, mu1, mu2,  0.0, prec, 1000, w, &val, &err);
  gsl_integration_workspace_free(w);
  if(fabs(err/val) > 5*prec) crash("In Compute_differential_decay_rate, cannot reach target precision: "+to_string_with_precision(prec,5));


  return exp(val);

  
		      

}

pair<double,double> WC(int id, double q2, double G1, double G2) {

  if( id == 1) return make_pair(0.241,0);
  if( id == 2) return make_pair(-1.1,0);
  if( id == 7 ) return make_pair(0.312,0);
  if( id == 10) return make_pair(4.64,0);
  if( id == 9) {

    double RE_C9 = -4.21;
    double IM_C9 = 0.0;
    double Cbar= 0.359/3.0;

    

    vector<string> Res({"Psi_2S", "Psi_3770", "Psi_4040", "Psi_4160", "Psi_4415"});
    double K= M_PI/180.0;
    Vfloat Bs({8e-3, 9.6e-6, 10.7e-6,6.9e-6,2.0e-5 });
    Vfloat Deltas({ 0, 0, 133*K, 301*K, 246*K});
    Vfloat Gammas({ 2.94e-4, 2.72e-2, 8e-2, 7e-2, 6.2e-2}); //GeV
    Vfloat Masses({3.686, 3.774, 4.039, 4.191, 4.421}); //GeV
    for(int i=0; i < (signed)Res.size() ; i++) {
      RE_C9 += ((9.0*M_PI)/(pow(alpha_em,2)))*Cbar*(1+G1)*(Masses[i]*Bs[i]*Gammas[i]/( pow( q2- pow(Masses[i],2),2) + pow( Masses[i]*Gammas[i],2)))*( cos(Deltas[i]+G2)*(q2 - pow(Masses[i],2)) + sin(Deltas[i]+G2)*Masses[i]*Gammas[i]);
      
      IM_C9 +=  ((9.0*M_PI)/(pow(alpha_em,2)))*Cbar*(1+G1)*(Masses[i]*Bs[i]*Gammas[i]/( pow( q2- pow(Masses[i],2),2) + pow( Masses[i]*Gammas[i],2)))*( sin(Deltas[i]+G2)*(q2 - pow(Masses[i],2)) -cos(Deltas[i]+G2)*Masses[i]*Gammas[i]);
    }
    
    return make_pair(RE_C9, IM_C9);
  }

  crash("In WC, you should not be here!");
  return make_pair(0,0);
}

//####################################################



void Get_Bs_xg_t_list(int num_xg) {

  for(int ixg=1;ixg<num_xg;ixg++) Bs_xg_t_list.push_back( ixg*0.1);
  return;
}



void Get_Bs_lattice_spacings_to_print() {

  int Nlat=300;
  double sx= 0.08*1.5/(Nlat-1.0); //fm
  for(int a=0;a<Nlat;a++) Bs_a_to_print.push_back( sx*a);


}

void Get_Bs_xg_to_spline() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { Bs_xg_to_spline.push_back(0.0+sx*ixg);}
    
  return;
}

void Get_Bs_xg_to_spline_VMD() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { Bs_xg_to_spline_VMD.push_back(0.0+sx*ixg);}
    
  return;
}


void Get_Bs_Mh_Tmin_Tmax(string obs, string mh,  int &Tmin, int &Tmax, int ixg, string Ens) {

  if(ixg==0) { Tmin=10; Tmax=15; return;}


  bool use_old=true;

  if(mh == "B0s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")) {

      
      
      if(obs=="FAu") {
	if(ixg==1) {  Tmin=16; Tmax=36;
	  if(Ens=="cA211a.12.48") { Tmin += 12; Tmax+=8;}
	  if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	  if(Ens=="cD211a.054.96") { Tmin -=6; Tmax=Tmin+14;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin -= 5; Tmax-=7;}
	  if(Ens=="cD211a.054.96") { Tmax-=5;}
	}
	else if(ixg==3) {  Tmin=21; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin -= 7; Tmax-=7;}
	}
	else if(ixg==4) {  Tmin=22; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=8;}
	  if(Ens=="cD211a.054.96") { Tmin -=1; Tmax-=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {	
	if(ixg==1) {  Tmin=18; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin += 7; Tmax+=7;}
	}
	else if(ixg==2) {  Tmin=19; Tmax=29;	}
	else if(ixg==3) {  Tmin=19; Tmax=29;  }
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=39; Tmax=49;
	  if(Ens=="cA211a.12.48") { Tmin -= 8; Tmax-=8;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin -=10; Tmax-=10;}
	}
	else if(ixg==3) {  Tmin=38; Tmax=48;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=10;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin -=10; Tmax-=10;}
	}
	else if(ixg==4) {  Tmin=38; Tmax=48;
	   if(Ens=="cA211a.12.48") { Tmin -= 12; Tmax-=12;}
	   if(Ens=="cB211b.072.64") { Tmin -=11; Tmax-=11;}
	   if(Ens=="cD211a.054.96") { Tmin -=12; Tmax-=12;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=35; Tmax=47;  }
	else if(ixg==2) {  Tmin=27; Tmax=38;
	  if(Ens=="cA211a.12.48") { Tmin += 7; Tmax+=7;}
	  if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") { Tmin +=9; Tmax+=9;}
	}
	else if(ixg==3) {  Tmin=29; Tmax=39;
	   if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	   if(Ens=="cD211a.054.96") { Tmin +=10; Tmax+=10;}
	}
	else if(ixg==4) {  Tmin=28; Tmax=40;
	  if(Ens=="cB211b.072.64") { Tmin +=9; Tmax+=9;}
	   if(Ens=="cD211a.054.96") { Tmin +=10; Tmax+=10;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=19; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin += 3;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 4; }
	}
	else if(ixg==3) {  Tmin=20; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 4; }
	}
	else if(ixg==4) {  Tmin=16; Tmax=30;
	  if(Ens=="cB211b.072.64") { Tmin -= 5;  Tmax -=5;}
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax+=5;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=23; Tmax=35;  }
	else if(ixg==2) {  Tmin=14; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin += 9; Tmax+=9;}
	}
	else if(ixg==3) {  Tmin=18; Tmax=31;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;}
	}
	else if(ixg==4) {  Tmin=18; Tmax=29;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=41; Tmax=50;	}
	else if(ixg==3) {  Tmin=38; Tmax=46;
	  if(Ens=="cB211b.072.64") { Tmin -= 12; Tmax-=12;}
	}
	else if(ixg==4) {  Tmin=38; Tmax=44;
	  if(Ens=="cB211b.072.64") { Tmin -= 15; Tmax-=11;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {  //quite bad
	if(ixg==1) {  Tmin=40; Tmax=51;  }
	else if(ixg==2) {  Tmin=40; Tmax=51;
	  if(Ens=="cA211a.12.48") { Tmin -= 8; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=38; Tmax=48;
	  if(Ens=="cA211a.12.48") { Tmin -= 9; Tmax-=9;}
	}
	else if(ixg==4) {  Tmin=37; Tmax=45;
	  if(Ens=="cA211a.12.48") { Tmin -= 10; Tmax-=5;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=16; Tmax=36;  }
	else if(ixg==2) {  Tmin=21; Tmax=32;  }
	else if(ixg==3) {  Tmin=21; Tmax=30;  }
	else if(ixg==4) {  Tmin=22; Tmax=32;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {	
	if(ixg==1) {  Tmin=18; Tmax=31;  }
	else if(ixg==2) {  Tmin=14; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=24;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=38; Tmax=55; use_old=false;  }
	else if(ixg==2) {  Tmin=35; Tmax=54; use_old=false;  }
	else if(ixg==3) {  Tmin=38; Tmax=48;  }
	else if(ixg==4) {  Tmin=27; Tmax=37;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=35; Tmax=47;  }
	else if(ixg==2) {  Tmin=36; Tmax=47;  }
	else if(ixg==3) {  Tmin=29; Tmax=39;  }
	else if(ixg==4) {  Tmin=28; Tmax=40;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=21; Tmax=40; use_old=false;  }
	else if(ixg==2) {  Tmin=25; Tmax=40; use_old=false;  }
	else if(ixg==3) {  Tmin=28; Tmax=40; use_old=false;  }
	else if(ixg==4) {  Tmin=27; Tmax=41;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=24; Tmax=36; use_old=false; }
	else if(ixg==2) {  Tmin=28; Tmax=39; use_old=false;  }
	else if(ixg==3) {  Tmin=32; Tmax=40; use_old=false;  }
	else if(ixg==4) {  Tmin=32; Tmax=51; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=39; Tmax=50;  }
	else if(ixg==2) {  Tmin=41; Tmax=50;  }
	else if(ixg==3) {  Tmin=35; Tmax=54; use_old=false;  }
	else if(ixg==4) {  Tmin=32; Tmax=54; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") { 
	if(ixg==1) {  Tmin=40; Tmax=57; use_old=false;  }
	else if(ixg==2) {  Tmin=35; Tmax=55; use_old=false;  }
	else if(ixg==3) {  Tmin=35; Tmax=55; use_old=false; }
	else if(ixg==4) {  Tmin=32; Tmax=55; use_old=false; }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }


  else if(mh == "B1s") {

    if( (Ens=="cA211a.12.48") || (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96") ) {

      if(obs=="FAu") {
       	if(ixg==1) {  Tmin=24; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin+=3;}
	  if(Ens=="cB211b.072.64") {Tmin+=4;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=29;
	  if(Ens=="cA211a.12.48") {Tmin-=4; Tmax -=4;}
	  if(Ens=="cD211a.054.96") {Tmin-=8; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=18; Tmax=26;	}
	else if(ixg==4) {  Tmin=15; Tmax=26;	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=25; Tmax=30;
	   if(Ens=="cA211a.12.48") {Tmin-=3;}
	   if(Ens=="cB211b.072.64") {Tmin-=3; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=17; Tmax=25;
	  if(Ens=="cD211a.054.96") {Tmin+=8; Tmax+=8;}
	}
	else if(ixg==3) {  Tmin=16; Tmax=29;
	  if(Ens=="cD211a.054.96") {Tmin+=8; Tmax+=8;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=29;
	  if(Ens=="cB211b.072.64") {Tmin-=3;}
	  if(Ens=="cD211a.054.96") {Tmin+=9; Tmax+=9;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=34; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	}
	else if(ixg==2) {  Tmin=34; Tmax=49;
	  if(Ens=="cA211a.12.48") {Tmin-=7;}
	}
	else if(ixg==3) {  Tmin=33; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=37;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=33; Tmax=44;  }
	else if(ixg==2) {  Tmin=27; Tmax=36;
	  if(Ens=="cA211a.12.48") {Tmin+=10; Tmax+=10;}
	}
	else if(ixg==3) {  Tmin=26; Tmax=37;  }
	else if(ixg==4) {  Tmin=24; Tmax=29;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax+=4;}
	  if(Ens=="cB211b.072.64") {Tmax+=6;}
	  if(Ens=="cD211a.054.96") {Tmin+=11; Tmax+=11;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=20; Tmax=28;  }
	else if(ixg==2) {  Tmin=18; Tmax=26;  }
	else if(ixg==3) {  Tmin=18; Tmax=26;  }
	else if(ixg==4) {  Tmin=18; Tmax=30;
	  if(Ens=="cD211a.054.96") {Tmin-=5; Tmax-=7;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=24; Tmax=35;  }
	else if(ixg==2) {  Tmin=24; Tmax=31;
	    if(Ens=="cA211a.12.48") {Tmin-=2; Tmax+=2;}
	}
	else if(ixg==3) {  Tmin=18; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin+=3; Tmax+=3;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin+=4; Tmax+=4;}
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=1;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=37; Tmax=50;
	  if(Ens=="cA211a.12.48") {Tmin-=6;}
	  if(Ens=="cD211a.054.96") {Tmin-=5; Tmax-=5;}
	}
	else if(ixg==2) {  Tmin=33; Tmax=50;
	  if(Ens=="cA211a.12.48") {Tmin-=5; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin-=3;}
	}
	else if(ixg==3) {  Tmin=25; Tmax=37;
	  if(Ens=="cB211b.072.64") {Tmax-=2; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=33;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=34; Tmax=45;
	  if(Ens=="cD211a.054.96") {Tmin-=3; Tmax+=5;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=7;}
	}
	else if(ixg==3) {  Tmin=26; Tmax=33;
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax+=8;}
	}
	else if(ixg==4) {  Tmin=23; Tmax=31;
	  if(Ens=="cA211a.12.48") {Tmin-=7;}
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=9; Tmax+=9;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {
    
      if(obs=="FAu") {
       	if(ixg==1) {  Tmin=24; Tmax=34;  }
	else if(ixg==2) {  Tmin=17; Tmax=29;  }
	else if(ixg==3) {  Tmin=16; Tmax=26;  }
	else if(ixg==4) {  Tmin=15; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=19; Tmax=37; use_old=false;  }
	else if(ixg==2) {  Tmin=15; Tmax=22;  }
	else if(ixg==3) {  Tmin=13; Tmax=29;  }
	else if(ixg==4) {  Tmin=12; Tmax=29;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=38; Tmax=51; use_old=false;  }
	else if(ixg==2) {  Tmin=34; Tmax=49;  }
	else if(ixg==3) {  Tmin=33; Tmax=45;  }
	else if(ixg==4) {  Tmin=24; Tmax=33; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=42; Tmax=53;  }
	else if(ixg==2) {  Tmin=27; Tmax=36;  }
	else if(ixg==3) {  Tmin=22; Tmax=31;  }
	else if(ixg==4) {  Tmin=22; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=20; Tmax=28;  }
	else if(ixg==2) {  Tmin=18; Tmax=26;  }
	else if(ixg==3) {  Tmin=18; Tmax=26;  }
	else if(ixg==4) {  Tmin=18; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=28; Tmax=47; use_old=false;  }
	else if(ixg==2) {  Tmin=26; Tmax=37; use_old=false;  }
	else if(ixg==3) {  Tmin=20; Tmax=32; use_old=false;  }
	else if(ixg==4) {  Tmin=19; Tmax=29; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=37; Tmax=50;  }
	else if(ixg==2) {  Tmin=28; Tmax=42;  }
	else if(ixg==3) {  Tmin=30; Tmax=47; use_old=false;  }
	else if(ixg==4) {  Tmin=24; Tmax=39; use_old=false;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=40; Tmax=45; }
	else if(ixg==2) {  Tmin=30; Tmax=44; use_old=false;}
	else if(ixg==3) {  Tmin=22; Tmax=33;  }
	else if(ixg==4) {  Tmin=17; Tmax=31;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else crash("Ensemble: "+Ens+" not yet implemented");
  }


  else if(mh == "B2s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=26; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=4;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin-=3;}
	}
	else if(ixg==3) {  Tmin=25; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin-=14; Tmax-=14;}
	}
	else if(ixg==4) {  Tmin=26; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=2;}
	  if(Ens=="cD211a.054.96") {Tmin-=11; Tmax-=11;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") {Tmin-=8;}
	  if(Ens=="cD211a.054.96") {Tmin-=12;}
	}
	else if(ixg==2) {  Tmin=24; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin-=5; Tmax-=7;}
	  if(Ens=="cB211b.072.64") {Tmin-=6; Tmax-=6;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=35;
	   if(Ens=="cA211a.12.48") {Tmin-=4;}
	   if(Ens=="cB211b.072.64") {Tmin-=6; Tmax-=6;}
	}
	else if(ixg==4) {  Tmin=24; Tmax=35;
	   if(Ens=="cA211a.12.48") {Tmin-=4; Tmax-=4;}
	   if(Ens=="cB211b.072.64") {Tmin-=4; Tmax-=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=31; Tmax=39;
	  if(Ens=="cA211a.12.48") {Tmin-=3; Tmax-=4;}
	  if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=33; Tmax=43;
	  if(Ens=="cA211a.12.48") {Tmin-=6; Tmax-=8;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=32;
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=32;
	  if(Ens=="cA211a.12.48") {Tmin-=3;}
	  if(Ens=="cB211b.072.64") {Tmin-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=31; Tmax=45;
	  if(Ens=="cA211a.12.48") {Tmin-=1;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=36;  }
	else if(ixg==3) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=4;}
	  if(Ens=="cB211b.072.64") {Tmin+=5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=8; Tmax-=8;}
	  if(Ens=="cB211b.072.64") {Tmin-=7; Tmax-=8;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=18; Tmax=25;  }
	else if(ixg==2) {  Tmin=18; Tmax=24;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=28;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=2;}
	}
	else if(ixg==2) {  Tmin=20; Tmax=30;  }
	else if(ixg==3) {  Tmin=20; Tmax=30;  }
	else if(ixg==4) {  Tmin=19; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=1; Tmax-=3;}
	  if(Ens=="cB211b.072.64") { Tmin -= 4; Tmax-=4;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=31; Tmax=41;
	  if(Ens=="cA211a.12.48") {Tmin-=1; Tmax-=1;}
	}
	else if(ixg==2) {  Tmin=31; Tmax=39;
	  if(Ens=="cA211a.12.48") {Tmin-=7; Tmax-=7;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=35;
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=4;   }
	}
	else if(ixg==4) {  Tmin=24; Tmax=38;
	  if(Ens=="cA211a.12.48") {Tmin-=4;}
	  if(Ens=="cB211b.072.64") { Tmin -= 3; Tmax-=3;   }
	  if(Ens=="cD211a.054.96") { Tmin += 3; Tmax+=3;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=29; Tmax=35;
	  if(Ens=="cA211a.12.48") {Tmin-=1; Tmax-=1;}
	  if(Ens=="cB211b.072.64") { Tmin += 1; Tmax+=3;   }
	  if(Ens=="cD211a.054.96") { Tmin += 3; Tmax+=3;   }
	}
	else if(ixg==2) {  Tmin=27; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=2; Tmax-=2;}
	}
	else if(ixg==3) {  Tmin=24; Tmax=34;
	  if(Ens=="cA211a.12.48") {Tmin-=5;}
	}
	else if(ixg==4) {  Tmin=20; Tmax=30;
	  if(Ens=="cA211a.12.48") {Tmin-=5;}
	  if(Ens=="cB211b.072.64") { Tmin -=6; Tmax-=6;   }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
	if(ixg==1) {  Tmin=26; Tmax=34;  }
	else if(ixg==2) {  Tmin=23; Tmax=30;  }
	else if(ixg==3) {  Tmin=15; Tmax=24;  }
	else if(ixg==4) {  Tmin=17; Tmax=25;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") {
	if(ixg==1) {  Tmin=18; Tmax=32;  }
	else if(ixg==2) {  Tmin=12; Tmax=33;  }
	else if(ixg==3) {  Tmin=14; Tmax=33;  }
	else if(ixg==4) {  Tmin=14; Tmax=33;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=31; Tmax=43;  }
	else if(ixg==2) {  Tmin=33; Tmax=43;  }
	else if(ixg==3) {  Tmin=28; Tmax=37;  }
	else if(ixg==4) {  Tmin=19; Tmax=35;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=33; Tmax=45;  }
	else if(ixg==2) {  Tmin=24; Tmax=36;  }
	else if(ixg==3) {  Tmin=24; Tmax=30;  }
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=16; Tmax=22;  }
	else if(ixg==2) {  Tmin=18; Tmax=23;  }
	else if(ixg==3) {  Tmin=20; Tmax=28;  }
	else if(ixg==4) {  Tmin=19; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=27; Tmax=34;  }
	else if(ixg==2) {  Tmin=20; Tmax=30;  }
	else if(ixg==3) {  Tmin=17; Tmax=27;  }
	else if(ixg==4) {  Tmin=12; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") {
	if(ixg==1) {  Tmin=31; Tmax=41;  }
	else if(ixg==2) {  Tmin=29; Tmax=39;  }
	else if(ixg==3) {  Tmin=24; Tmax=38;  }
	else if(ixg==4) {  Tmin=24; Tmax=38;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=27; Tmax=34;  }
	else if(ixg==2) {  Tmin=25; Tmax=34;  }
	else if(ixg==3) {  Tmin=24; Tmax=34;  }
	else if(ixg==4) {  Tmin=20; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }



  else if(mh == "B3s") {

    if( (Ens=="cA211a.12.48") || (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") {
       if(ixg==1) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	 if(Ens=="cD211a.054.96") { Tmin += 6; Tmax +=6 ;}
       }
       else if(ixg==2) {  Tmin=20; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
       }
       else if(ixg==3) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
       }
       else if(ixg==4) {  Tmin=21; Tmax=27;
	 if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	 if(Ens=="cD211a.054.96") { Tmin -= 3; Tmax -=3 ;}
	 
       }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //non bellissimi
	if(ixg==1) {  Tmin=17; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin+=2; Tmax+=1;}
	}
	else if(ixg==2) {  Tmin=24; Tmax=31;
	  if(Ens=="cA211a.12.48") { Tmin-=12; Tmax-=12;}
	  if(Ens=="cB211b.072.64") { Tmin-=2; Tmax-=2;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=33;
	  if(Ens=="cA211a.12.48") { Tmin-=12; Tmax-=16;}
	  if(Ens=="cB211b.072.64") { Tmin-=2; Tmax-=2;}
	}
	else if(ixg==4) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=1;}
	  if(Ens=="cD211a.054.96") { Tmin -= 4; Tmax -=4 ;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=28; Tmax=40;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	}
	else if(ixg==2) {  Tmin=27; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	  if(Ens=="cD211a.054.96") { Tmin -= 1; Tmax +=4 ;}
	}
	else if(ixg==3) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=3; Tmax+=3;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=26; Tmax=32;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=4;}
	  if(Ens=="cD211a.054.96") { Tmin+=3; Tmax+=3;}
	}
	else if(ixg==2) {  Tmin=23; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4;}
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax +=4 ;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") { Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") { Tmin += 2; Tmax +=2 ;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=17; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	}
	else if(ixg==2) {  Tmin=13; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin+=3; Tmax-=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=5; Tmax+=5;}
	}
	else if(ixg==3) {  Tmin=13; Tmax=24;
	   if(Ens=="cA211a.12.48") { Tmin+=3; }
	   if(Ens=="cB211b.072.64") { Tmin+=3;}
	   if(Ens=="cD211a.054.96") { Tmin += 3;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=20; Tmax=27;
	   if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=2;}
	}
	else if(ixg==2) {  Tmin=17; Tmax=27;
	  if(Ens=="cA211a.12.48") { Tmin+=4; Tmax+=2;}
	}
	else if(ixg==3) {  Tmin=16; Tmax=27;
	  if(Ens=="cB211b.072.64") { Tmin-=5; Tmax-=5;}
	}
	else if(ixg==4) {  Tmin=14; Tmax=22;
	  if(Ens=="cB211b.072.64") { Tmin-=4; Tmax-=7;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //two plataux for small xg
	if(ixg==1) {  Tmin=31; Tmax=38;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") { Tmin-=6; Tmax-=6;}
	}
	else if(ixg==2) {  Tmin=25; Tmax=30;
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax +=5 ;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=28;
	  if(Ens=="cB211b.072.64") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cD211a.054.96") { Tmin += 8; Tmax +=8 ;}
	}
	else if(ixg==4) {  Tmin=14; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin+=2;}
	  if(Ens=="cB211b.072.64") { Tmin+=6; Tmax+=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=40;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=1;}
	  if(Ens=="cB211b.072.64") { Tmin+=2; Tmax+=2;}
	  if(Ens=="cD211a.054.96") {Tmin +=3;}
	}
	else if(ixg==2) {  Tmin=25; Tmax=34;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=9;}
	  if(Ens=="cD211a.054.96") {Tmin +=3; Tmax+=3;}
	}
	else if(ixg==3) {  Tmin=20; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4;}
	  if(Ens=="cB211b.072.64") { Tmin-=2; Tmax-=2;}
	  if(Ens=="cD211a.054.96") {Tmin +=6; Tmax+=6;}
	}
	else if(ixg==4) {  Tmin=17; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=4;}
	  if(Ens=="cB211b.072.64") { Tmin-=3;}
	  if(Ens=="cD211a.054.96") {Tmin +=2; Tmax+=2;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }

    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") {
       if(ixg==1) {  Tmin=21; Tmax=27;  }
       else if(ixg==2) {  Tmin=20; Tmax=25;  }
       else if(ixg==3) {  Tmin=21; Tmax=27;  }
       else if(ixg==4) {  Tmin=21; Tmax=27;  }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //non bellissimi
	if(ixg==1) {  Tmin=17; Tmax=24;  }
	else if(ixg==2) {  Tmin=24; Tmax=31;  }
	else if(ixg==3) {  Tmin=22; Tmax=32;  }
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=30; Tmax=40;  }
	else if(ixg==2) {  Tmin=26; Tmax=33;  }
	else if(ixg==3) {  Tmin=22; Tmax=30;  }
	else if(ixg==4) {  Tmin=16; Tmax=30;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=26; Tmax=32;  }
	else if(ixg==2) {  Tmin=24; Tmax=31;  }
	else if(ixg==3) {  Tmin=20; Tmax=25;  }
	else if(ixg==4) {  Tmin=16; Tmax=25;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=17; Tmax=26;  }
	else if(ixg==2) {  Tmin=19; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=20; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=20; Tmax=27;  }
	else if(ixg==2) {  Tmin=21; Tmax=31;  }
	else if(ixg==3) {  Tmin=14; Tmax=27;  }
	else if(ixg==4) {  Tmin=14; Tmax=22;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //two plataux for small xg
	if(ixg==1) {  Tmin=31; Tmax=38;  }
	else if(ixg==2) {  Tmin=26; Tmax=30;  }
	else if(ixg==3) {  Tmin=25; Tmax=33;  }
	else if(ixg==4) {  Tmin=23; Tmax=31;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=24; Tmax=40;  }
	else if(ixg==2) {  Tmin=21; Tmax=30;  }
	else if(ixg==3) {  Tmin=20; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    
    else crash("Ensemble: "+Ens+" not yet implemented");
  }

  
  else if(mh == "B4s") {

    if( (Ens=="cA211a.12.48") ||  (Ens=="cB211b.072.64") || (Ens=="cB211b.072.96") || (Ens=="cD211a.054.96")  ) {

      if(obs=="FAu") { //not very nice
       if(ixg==1) {  Tmin=18; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=1;}
	 if(Ens=="cB211b.072.64") {Tmin+=1;}
       }
       else if(ixg==2) {  Tmin=20; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=3;}
	 if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=1;}
       }
       else if(ixg==3) {  Tmin=21; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=4;}
	 if(Ens=="cB211b.072.64") {Tmin-=1; Tmax+=1;}
       }
       else if(ixg==4) {  Tmin=21; Tmax=25;
	 if(Ens=="cA211a.12.48") { Tmin-=5;}
	 if(Ens=="cB211b.072.64") {Tmin-=2; Tmax+=2;}
       }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //not very nice
	if(ixg==1) {  Tmin=21; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=2; Tmax+=6;}
	  if(Ens=="cD211a.054.96") { Tmin += 2; Tmax+=2;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") {Tmin+=4; Tmax+=4;}
	  if(Ens=="cD211a.054.96") { Tmin += 4; Tmax+=9;}
	}
	else if(ixg==3) {  Tmin=21; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	  if(Ens=="cB211b.072.64") {Tmin+=4; Tmax+=4;}
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax+=5;}
	}
	else if(ixg==4) {  Tmin=21; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=2;}
	  if(Ens=="cB211b.072.64") {Tmin-=2; Tmax-=2;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=28; Tmax=33;
	  if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=5;}
	}
	else if(ixg==2) {  Tmin=22; Tmax=28;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=2;}
	  if(Ens=="cB211b.072.64") {Tmin+=3; Tmax+=3;}
	  if(Ens=="cD211a.054.96") { Tmin += 5; Tmax+=5;}
	}
	else if(ixg==3) {  Tmin=19; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else if(ixg==4) {  Tmin=16; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=2;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=22; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=2;}
	  if(Ens=="cD211a.054.96") { Tmin += 7; Tmax+=4;}
	}
	else if(ixg==2) {  Tmin=21; Tmax=27;
	   if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	   if(Ens=="cD211a.054.96") { Tmin += 3;}
	}
	else if(ixg==3) {  Tmin=17; Tmax=24;
	  if(Ens=="cB211b.072.64") {Tmin+=1; Tmax-=1;}
	  if(Ens=="cD211a.054.96") {Tmin+=4; Tmax += 4; }
	}
	else if(ixg==4) {  Tmin=14; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax-=6;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=15; Tmax=24; if(Ens=="cB211b.072.96") Tmin -=4;  }
	else if(ixg==2) {  Tmin=17; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=18; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=18; Tmax=27;
	  if(Ens=="cA211a.12.48") { Tmin-=1; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=1; Tmax += 4; }
	}
	else if(ixg==2) {  Tmin=18; Tmax=28;
	  if(Ens=="cA211a.12.48") { Tmin-=2; Tmax=Tmin+4;}
	  if(Ens=="cD211a.054.96") {Tmin-=2; Tmax += 3; }
	}
	else if(ixg==3) {  Tmin=15; Tmax=23;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=7;}
	}
	else if(ixg==4) {  Tmin=12; Tmax=21;
	  if(Ens=="cD211a.054.96") {Tmin+=10; Tmax += 10; }
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //bad plateaux at  xg=0.1
	if(ixg==1) {  Tmin=26; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=5; Tmax-=5;}
	  if(Ens=="cB211b.072.64") {Tmin-=1; Tmax += 2; }
	  if(Ens=="cD211a.054.96") {Tmin+=2; Tmax += 2; }
	}
	else if(ixg==2) {  Tmin=21; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=2;}
	  if(Ens=="cB211b.072.64") {Tmin+=4; Tmax += 4; }
	  if(Ens=="cD211a.054.96") {Tmin+=8; Tmax += 8; }
	}
	else if(ixg==3) {  Tmin=22; Tmax=25;
	  if(Ens=="cA211a.12.48") { Tmin-=4; Tmax-=2;}
	  if(Ens=="cD211a.054.96") {Tmin+=4; Tmax += 4; }
	}
	else if(ixg==4) {  Tmin=20; Tmax=26;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=3;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=30;
	  if(Ens=="cA211a.12.48") { Tmin-=7; Tmax-=7;}
	  if(Ens=="cD211a.054.96") {Tmin+=3; Tmax += 6; }
	}
	else if(ixg==2) {  Tmin=24; Tmax=30;
	    if(Ens=="cA211a.12.48") { Tmin-=6; Tmax-=6;}
	    if(Ens=="cD211a.054.96") {Tmin+=3; Tmax += 6; }
	}
	else if(ixg==3) {  Tmin=16; Tmax=24;
	  if(Ens=="cA211a.12.48") { Tmin-=3; Tmax-=5;}
	  if(Ens=="cD211a.054.96") {Tmin+=5; Tmax += 3; }
	}
	else if(ixg==4) {  Tmin=11; Tmax=17;
	  if(Ens=="cA211a.12.48") { Tmin+=1; Tmax+=1;}
	}
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else if(Ens=="cC211a.06.80") {

      if(obs=="FAu") { //not very nice
       if(ixg==1) {  Tmin=18; Tmax=25;  }
       else if(ixg==2) {  Tmin=18; Tmax=25;  }
       else if(ixg==3) {  Tmin=18; Tmax=25;  }
       else if(ixg==4) {  Tmin=21; Tmax=25;  }
       else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FAd") { //not very nice
	if(ixg==1) {  Tmin=21; Tmax=24;  }
	else if(ixg==2) {  Tmin=27; Tmax=33;  }
	else if(ixg==3) {  Tmin=27; Tmax=34;  }
	else if(ixg==4) {  Tmin=21; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVu") {
	if(ixg==1) {  Tmin=25; Tmax=33;  }
	else if(ixg==2) {  Tmin=28; Tmax=33;  }
	else if(ixg==3) {  Tmin=21; Tmax=24;  }
	else if(ixg==4) {  Tmin=16; Tmax=24;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FVd") {
	if(ixg==1) {  Tmin=24; Tmax=30;  }
	else if(ixg==2) {  Tmin=24; Tmax=27;  }
	else if(ixg==3) {  Tmin=18; Tmax=27;  }
	else if(ixg==4) {  Tmin=15; Tmax=23;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Tu") {
	if(ixg==1) {  Tmin=15; Tmax=24;  }
	else if(ixg==2) {  Tmin=17; Tmax=24;  }
	else if(ixg==3) {  Tmin=17; Tmax=24;  }
	else if(ixg==4) {  Tmin=18; Tmax=28;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FA_Td") {
	if(ixg==1) {  Tmin=18; Tmax=27;  }
	else if(ixg==2) {  Tmin=21; Tmax=29;  }
	else if(ixg==3) {  Tmin=11; Tmax=20;  }
	else if(ixg==4) {  Tmin=12; Tmax=21;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Tu") { //bad plateaux at  xg=0.1
	if(ixg==1) {  Tmin=28; Tmax=32;  }
	else if(ixg==2) {  Tmin=29; Tmax=33;  }
	else if(ixg==3) {  Tmin=26; Tmax=29;  }
	else if(ixg==4) {  Tmin=20; Tmax=26;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
      else if(obs=="FV_Td") {
	if(ixg==1) {  Tmin=26; Tmax=30;  }
	else if(ixg==2) {  Tmin=24; Tmax=30;  }
	else if(ixg==3) {  Tmin=16; Tmax=24;  }
	else if(ixg==4) {  Tmin=11; Tmax=17;  }
	else crash("ixg: "+to_string(ixg)+" not yet implemented");
      }
    }
    else crash("Ensemble: "+Ens+" not yet implemented");
  }

  else crash("mh: "+mh+" not yet implemented");
  

  
  if(use_old && Ens=="cC211a.06.80") { Tmin= (int)(Tmin*0.07957/0.0682); Tmax= (int)( Tmax*0.07957/0.0682);}

  if(use_old && Ens=="cD211a.054.96") { Tmin= (int)(Tmin*0.07957/0.0569); Tmax= (int)( Tmax*0.07957/0.0569);}

  if(use_old && Ens=="cA211a.12.48") { Tmin=(int)(Tmin*0.07957/0.0907593); Tmax= (int)(Tmax*0.07957/0.0907593);}
 

  return;

}











void Compute_Bs_mumu_gamma() {

  




  
  Get_Bs_xg_to_spline();
  Get_Bs_xg_to_spline_VMD();
  Get_Bs_lattice_spacings_to_print();

  GaussianMersenne GM(61111223);

  //read masses
  data_t M1_SMSM, M2_SMSM, M3_SMSM, M4_SMSM;
  data_t M1_SM, M2_SM, M3_SM, M4_SM;

  M1_SM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M1_SM_2", "P5P5");
  M2_SM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M2_SM_2", "P5P5");
  M3_SM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M3_SM_2", "P5P5");
  M4_SM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M4_SM_2", "P5P5");

  M1_SMSM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M1_SMSM_2", "P5P5");
  M2_SMSM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M2_SMSM_2", "P5P5");
  M3_SMSM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M3_SMSM_2", "P5P5");
  M4_SMSM.Read("../Bs_mumu_gamma_data/masses", "mes_contr_2pts_M4_SMSM_2", "P5P5");

  int Ns= M1_SM.size;

  for(int i=0;i<Ns;i++) {

    boost::filesystem::create_directory( "../data/mass_list/"+M1_SM.Tag[i]);
    CorrAnalysis Corr;
    Corr.Reflection_sign=1;
    Corr.UseJack=1;
    Corr.Njacks=20;
    Corr.Nt = M1_SM.nrows[i];
    distr_t_list eff_M1 = Corr.effective_mass_t( M1_SM.col(0)[i], "../data/mass_list/"+M1_SM.Tag[i]+"/eff_mass_M1_SM");
    distr_t_list eff_M2 = Corr.effective_mass_t( M2_SM.col(0)[i], "../data/mass_list/"+M1_SM.Tag[i]+"/eff_mass_M2_SM");
    distr_t_list eff_M3 = Corr.effective_mass_t( M3_SM.col(0)[i], "../data/mass_list/"+M1_SM.Tag[i]+"/eff_mass_M3_SM");
    distr_t_list eff_M4 = Corr.effective_mass_t( M4_SM.col(0)[i], "../data/mass_list/"+M1_SM.Tag[i]+"/eff_mass_M4_SM");

    distr_t_list eff_M1_SMSM = Corr.effective_mass_t( M1_SMSM.col(0)[i], "../data/mass_list/"+M1_SMSM.Tag[i]+"/eff_mass_M1_SMSM");
    distr_t_list eff_M2_SMSM = Corr.effective_mass_t( M2_SMSM.col(0)[i], "../data/mass_list/"+M1_SMSM.Tag[i]+"/eff_mass_M2_SMSM");
    distr_t_list eff_M3_SMSM = Corr.effective_mass_t( M3_SMSM.col(0)[i], "../data/mass_list/"+M1_SMSM.Tag[i]+"/eff_mass_M3_SMSM");
    distr_t_list eff_M4_SMSM = Corr.effective_mass_t( M4_SMSM.col(0)[i], "../data/mass_list/"+M1_SMSM.Tag[i]+"/eff_mass_M4_SMSM");
   
  }


  vector<rt_FF_Bs> FF_ret_list;

  
 

  vector<bool> Use_three_finest_list({0,0,0,0,0});
  vector<bool> Include_a4_list({0,0,0,0,0});
  vector<bool> UseJack_list({1,1,1,1,1});
  vector<int>  num_xg_list({5,5,5,5,5});
  vector<double> Perform_continuum_extrapolation_list({1,1,1,1,1});
  vector<string> Corr_path_list({ "../Bs_mumu_gamma_data/mh0",  "../Bs_mumu_gamma_data/mh1", "../Bs_mumu_gamma_data/mh2", "../Bs_mumu_gamma_data/mh3", "../Bs_mumu_gamma_data/mh4"});
  vector<string> out_tag_list({"mh0", "mh1", "mh2", "mh3", "mh4"});
  vector<string> Meson_list({"B0s", "B1s", "B2s", "B3s", "B4s"});

  int N= UseJack_list.size();
  for(int i=0;i<N;i++) {
    Bs_xg_t_list.clear();
    Get_Bs_xg_t_list(num_xg_list[i]);
    MESON=Meson_list[i];
    //create output directories
    boost::filesystem::create_directory("../data/ph_emission");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/mass");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/C");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/H");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/C_T");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/H_T");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC"); 
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/per_kin");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_u");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_u/per_kin");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_d");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+MESON+"/FF_d/per_kin");

    
    string Fit_tag= ( (Include_a4_list[i]==true)?"wa4_":"");
    Fit_tag = ( (Use_three_finest_list[i]==true)?"wtf_":Fit_tag);
    FF_ret_list.push_back(Get_Bs_mumu_gamma_form_factors(num_xg_list[i], Perform_continuum_extrapolation_list[i], Use_three_finest_list[i], Include_a4_list[i], UseJack_list[i], Fit_tag, Corr_path_list[i], out_tag_list[i]));
  }


  //Fit F_Bs

  

  


  cout<<"Generating data to be used for meson mass extrapolation"<<endl;
  //print continuum extrapolated results for each contribution as a function of the mass
  vector<string> contribs({"FA", "FA_u", "FA_d", "FV",  "FV_u", "FV_d",  "FA_T", "FA_T_u", "FA_T_d", "FV_T", "FV_T_u", "FV_T_d"});
  for(int c=0; c<(signed)contribs.size(); c++) {

    distr_t_list FF_Bs_list;
    
    distr_t_list masses(1), fps(1);
    vector<distr_t_list> FF_M;
    for(int ixg=1; ixg<num_xg_list[0]; ixg++) FF_M.emplace_back(1);
    for(int m=0; m<(signed)FF_ret_list.size(); m++) {
      masses.distr_list.push_back( FF_ret_list[m].mass );
      fps.distr_list.push_back( FF_ret_list[m].fp);
      for(int ixg=1; ixg<num_xg_list[0];ixg++) { FF_M[ixg-1].distr_list.push_back( FF_ret_list[m].Get_FF(c).distr_list[ixg-1]);}
    }

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+"Bs_extr");
    for(int ixg=1;ixg<(signed)num_xg_list[0];ixg++) {
     
      Print_To_File({}, { ratio_mh, (1.0/masses).ave(), (1.0/masses).err(), fps.ave(), fps.err(), FF_M[ixg-1].ave(), FF_M[ixg-1].err(), (FF_M[ixg-1]/fps).ave(), (FF_M[ixg-1]/fps).err() }, "../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_ixg_"+to_string(ixg)+".dat", "", "");


      //perform extrapolation in the heavy quark mass
      if(contribs[c] != "FA" && contribs[c] != "FV" && contribs[c] != "FV_T" && contribs[c] != "FA_T" ) {
	
	 
	class ipar_MEX {
	public:
	  ipar_MEX() : FF(0.0), FF_err(0.0) {}
	  double FF, FF_err, M, xg;
	};
  
	class fpar_MEX {
	public:
	  fpar_MEX() {}
	  fpar_MEX(const Vfloat &par) {
	    if((signed)par.size() != 3) crash("In class fpar_MEX  class constructor Vfloat par has size != 3");
	    A=par[0];
	    B=par[1];
	    C=par[2];
	  }
	  double A,B,C;
	};
	
	
	//init bootstrap fit
	bootstrap_fit<fpar_MEX,ipar_MEX> bf_MEX(NJ);
	bf_MEX.set_warmup_lev(0); //sets warmup
	bf_MEX.Set_number_of_measurements(Nmasses);
	bf_MEX.Set_verbosity(1);
	
	double guess_B= sqrt(2*0.24/Bs_xg_t_list[ixg-1]);
	if(contribs[c].substr(1,1)=="A")  guess_B=  sqrt(2*0.5/Bs_xg_t_list[ixg-1]);
	
	bf_MEX.Add_par("A", (FF_M[ixg-1]/fps).ave(0), (FF_M[ixg-1]/fps).ave(0)/10);
	bf_MEX.Add_par("B", guess_B , guess_B/10);
	bf_MEX.Add_par("C", 0.1, 0.001);
	//fit on mean values to get ch2
	bootstrap_fit<fpar_MEX,ipar_MEX> bf_MEX_ch2(1);
	bf_MEX_ch2.set_warmup_lev(2); //sets warmup
	bf_MEX_ch2.Set_number_of_measurements(Nmasses);
	bf_MEX_ch2.Set_verbosity(1);
	bf_MEX_ch2.Add_par("A", (FF_M[ixg-1]/fps).ave(0), (FF_M[ixg-1]/fps).ave(0)/10);
	bf_MEX_ch2.Add_par("B", guess_B, guess_B/10); 
	bf_MEX_ch2.Add_par("C", 0.1, 0.001); 
	bool C_fixed=false;
	if(contribs[c]=="FA_d") C_fixed=true;
	if(C_fixed) {
	  bf_MEX.Fix_par("C", 0.0);
	  bf_MEX_ch2.Fix_par("C", 0.0);
	}
	
	//ansatz
	bf_MEX.ansatz=  [&contribs, &c ](const fpar_MEX &p, const ipar_MEX &ip) {

	  bool is_axial= (contribs[c].substr(1,1) == "A");

	  double x= ip.M*ip.M;
	  
	  if(is_axial) x= ip.M;

	  double func= (p.A/( 1 + p.B*p.B/x ))*( 1 + p.C/ip.M) ; 

	  if(contribs[c].substr( contribs[c].length() -1, 1) == "u" ) func/= ip.M;

	  return func;
	  
	  
	};
	bf_MEX.measurement=  [ ](const fpar_MEX &p, const ipar_MEX &ip) {
	  
	  return ip.FF;
	};
	bf_MEX.error=  [ ](const fpar_MEX &p, const ipar_MEX &ip) {
	  
	  return ip.FF_err;
	};
	
	bf_MEX_ch2.ansatz= bf_MEX.ansatz;
	bf_MEX_ch2.measurement = bf_MEX.measurement;
	bf_MEX_ch2.error = bf_MEX.error;
      
	//start fitting
	//fill the data
	vector<vector<ipar_MEX>> data_MEX(NJ);
	vector<vector<ipar_MEX>> data_MEX_ch2(1);
	//allocate space for output result
	boot_fit_data<fpar_MEX> Bt_fit_MEX;
	boot_fit_data<fpar_MEX> Bt_fit_MEX_ch2;
	for(auto &data_iboot: data_MEX) data_iboot.resize(Nmasses);
	for(auto &data_iboot: data_MEX_ch2) data_iboot.resize(Nmasses);
	for(int ijack=0;ijack<NJ;ijack++) {
	  for(int im=0;im<Nmasses;im++) {
	    data_MEX[ijack][im].FF = (FF_M[ixg-1].distr_list[im]/fps.distr_list[im]).distr[ijack];
	    data_MEX[ijack][im].FF_err= (FF_M[ixg-1]/fps).err(im);
	    data_MEX[ijack][im].xg= Bs_xg_t_list[ixg-1];
	    data_MEX[ijack][im].M= masses.distr_list[im].distr[ijack];
	    if(ijack==0) {
	      data_MEX_ch2[ijack][im].FF = (FF_M[ixg-1]/fps).ave(im);
	      data_MEX_ch2[ijack][im].FF_err= (FF_M[ixg-1]/fps).err(im);
	      data_MEX_ch2[ijack][im].xg= Bs_xg_t_list[ixg-1];
	      data_MEX_ch2[ijack][im].M = masses.ave(im);
	      
	    }
	  }
	}
	
	//append
	bf_MEX.Append_to_input_par(data_MEX);
	bf_MEX_ch2.Append_to_input_par(data_MEX_ch2);
	//fit
	cout<<"Fitting "<<contribs[c]<<" to extrapolate at the Bs-meson mass"<<endl;
	Bt_fit_MEX= bf_MEX.Perform_bootstrap_fit();
	Bt_fit_MEX_ch2= bf_MEX_ch2.Perform_bootstrap_fit();
	int Npars= C_fixed?2:3;
	double ch2_red_MEX= Bt_fit_MEX_ch2.get_ch2_ave()/( Nmasses -Npars);
	
	//retrieve params
	distr_t A(1), B(1), C(1);
	for(int ijack=0;ijack<NJ;ijack++) {
	  A.distr.push_back( Bt_fit_MEX.par[ijack].A);
	  B.distr.push_back( Bt_fit_MEX.par[ijack].B);
	  C.distr.push_back( Bt_fit_MEX.par[ijack].C);
	}
	
	//print fit function
	distr_t_list F_MEX_fit(1);
	Vfloat Inv_Masses_to_print;

	//### Determine FF at Bs meson mass
	ipar_MEX pp_Bs;
	double MBs=5.367;
	pp_Bs.M= MBs;
	distr_t FF_Bs(1);
	for(int ijack=0;ijack<NJ;ijack++) FF_Bs.distr.push_back( bf_MEX.ansatz( Bt_fit_MEX.par[ijack], pp_Bs));
	FF_Bs_list.distr_list.push_back( FF_Bs);
	//#################################

	//######### PRINT FIT FUNCTION ##########
	for(int i=0;i<1000;i++) { Inv_Masses_to_print.push_back( 0.7*2*i/999  ) ;}
	//plot fit function
	for(auto &Minv: Inv_Masses_to_print) {
	  ipar_MEX pp_MEX;
	  pp_MEX.M=1.0/Minv;
	  distr_t F_MEX_M(1);
	  for(int ijack=0;ijack<NJ;ijack++) {
	    F_MEX_M.distr.push_back( bf_MEX.ansatz( Bt_fit_MEX.par[ijack], pp_MEX));
	  }
	  
	  F_MEX_fit.distr_list.push_back( F_MEX_M);
	}
	//print
	string out_path="../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_ixg_"+to_string(ixg)+".fit_func";
	Print_To_File({},{ Inv_Masses_to_print, F_MEX_fit.ave(), F_MEX_fit.err()} , out_path , "", "ch2/dof: "+to_string_with_precision(ch2_red_MEX,5));
	
      }
    }
    //print FF_Bs
    if(contribs[c] != "FA" && contribs[c] != "FV" && contribs[c] != "FV_T" && contribs[c] != "FA_T" ) {
      Print_To_File({}, {Bs_xg_t_list, FF_Bs_list.ave(), FF_Bs_list.err()}, "../data/ph_emission/"+ph_type+"/Bs_extr/"+contribs[c]+"_Bs.dat", "", "");
    }
  }
  
  exit(-1);
  MESON="Bs";
  
  
 
  
  Bs_xg_t_list.clear();
  int nxg=11;
  Get_Bs_xg_t_list(nxg);


  vector<distr_t_list> FF_final;

  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> F_interpolated;

  distr_t_list Coupling_1_list(1), Pole_1_list(1), Coupling_2_list(1), Pole_2_list(1), Bg_list(1);

  

  for( int i=0; i<(signed)contribs.size(); i++) {

 

    distr_t_list FF_jack(1), FF_boot(0);
    distr_t_list FF_jack_AIC(1), FF_boot_AIC(0);

    for(int ixg=1; ixg<nxg;ixg++) {

      rt_FF_Bs Fit_A_jack = FF_ret_list[1];
      rt_FF_Bs Fit_B_jack = FF_ret_list[2];
      //rt_FF_Bs Fit_A_boot = FF_ret_list[1];
      //rt_FF_Bs Fit_B_boot = FF_ret_list[4];

      double wA = 1.0;
      double wB = 1.0;
      double sum=wA+wB;
      wA /= sum;
      wB /= sum;

      double wA_AIC= exp(-0.5*(Fit_A_jack.Get_ch2(i)[ixg-1]*Fit_A_jack.Ndof + 2*Fit_A_jack.Npars - 2*Fit_A_jack.Nmeas));
      double wB_AIC =exp(-0.5*(Fit_B_jack.Get_ch2(i)[ixg-1]*Fit_B_jack.Ndof + 2*Fit_B_jack.Npars - 2*Fit_B_jack.Nmeas));
      double sum_AIC= wA_AIC+wB_AIC;
      wA_AIC /= sum_AIC;
      wB_AIC /= sum_AIC;

    
      

      distr_t D_FF_jack= wA*Fit_A_jack.Get_FF(i).distr_list[ixg-1] + wB*Fit_B_jack.Get_FF(i).distr_list[ixg-1];
      //distr_t D_FF_boot= wA*Fit_A_boot.Get_FF(i).distr_list[ixg-1] + wB*Fit_B_boot.Get_FF(i).distr_list[ixg-1];

      distr_t D_FF_jack_AIC= wA_AIC*Fit_A_jack.Get_FF(i).distr_list[ixg-1] + wB_AIC*Fit_B_jack.Get_FF(i).distr_list[ixg-1];
      //distr_t D_FF_boot_AIC= wA_AIC*Fit_A_boot.Get_FF(i).distr_list[ixg-1] + wB_AIC*Fit_B_boot.Get_FF(i).distr_list[ixg-1];


      

      double syst= wA*pow( Fit_A_jack.Get_FF(i).ave(ixg-1) - D_FF_jack.ave(),2) + wB*pow( Fit_B_jack.Get_FF(i).ave(ixg-1) - D_FF_jack.ave(),2);
      syst = sqrt(syst);
      distr_t syst_jack(1), syst_boot(0);

      double syst_AIC= wA_AIC*pow( Fit_A_jack.Get_FF(i).ave(ixg-1) - D_FF_jack_AIC.ave(),2) + wB_AIC*pow( Fit_B_jack.Get_FF(i).ave(ixg-1) - D_FF_jack_AIC.ave(),2);
      syst_AIC = sqrt(syst_AIC);
      distr_t syst_AIC_jack(1), syst_AIC_boot(0);

      for(int iboot=0;iboot<Nboots;iboot++) { double rn= GM();  syst_boot.distr.push_back( rn*syst); syst_AIC_boot.distr.push_back( rn*syst_AIC);}
      for(int ijack=0;ijack<NJ;ijack++)     { double rn= GM();  syst_jack.distr.push_back( rn*syst/sqrt(NJ-1.0)); syst_AIC_jack.distr.push_back( rn*syst_AIC/sqrt(NJ-1.0));}

      D_FF_jack = D_FF_jack + syst_jack;
      //D_FF_boot = D_FF_boot + syst_boot;

      D_FF_jack_AIC = D_FF_jack_AIC + syst_AIC_jack;
      //D_FF_boot_AIC = D_FF_boot_AIC + syst_AIC_boot;
      
      FF_jack.distr_list.push_back(  D_FF_jack );
      //FF_boot.distr_list.push_back(  D_FF_boot );
      FF_jack_AIC.distr_list.push_back( D_FF_jack_AIC );
      //FF_boot_AIC.distr_list.push_back( D_FF_boot_AIC );



      cout<<"#########################"<<endl;
      cout<<"Performing BMA for contrib: "<<contribs[i]<<" ixg: "<<0.1*ixg<<endl;
      cout<<"Four lattice spacings (Meas A): "<<endl;
      cout<<"w(U): "<<wA<<" w(AIC): "<<wA_AIC<<endl;
      cout<<"Nmeas: "<<Fit_A_jack.Nmeas<<endl;
      cout<<"Ndof: "<<Fit_A_jack.Ndof<<endl;
      cout<<"Npars: "<<Fit_A_jack.Npars<<endl;
      cout<<"ch2/dof: "<<Fit_A_jack.Get_ch2(i)[ixg-1]<<endl;
      cout<<"Three finest (Meas B): "<<endl;
      cout<<"w(U): "<<wB<<" w(AIC): "<<wB_AIC<<endl;
      cout<<"Nmeas: "<<Fit_B_jack.Nmeas<<endl;
      cout<<"Ndof: "<<Fit_B_jack.Ndof<<endl;
      cout<<"Npars: "<<Fit_B_jack.Npars<<endl;
      cout<<"ch2/dof: "<<Fit_B_jack.Get_ch2(i)[ixg-1]<<endl;
      cout<<"Meas A: "<<Fit_A_jack.Get_FF(i).ave(ixg-1)<<" +- "<<Fit_A_jack.Get_FF(i).err(ixg-1)<<endl;
      cout<<"Meas B: "<<Fit_B_jack.Get_FF(i).ave(ixg-1)<<" +- "<<Fit_B_jack.Get_FF(i).err(ixg-1)<<endl;
      
      cout<<"Combined(U): "<<D_FF_jack.ave()<<" +- "<<D_FF_jack.err()<<endl;
      cout<<"Expected error stat(U): "<< sqrt( wA*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"Expected error syst(U): "<< syst<<endl;
      cout<<"stat+syst(U): "<< sqrt( syst*syst + wA*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      
      cout<<"Combined(AIC): "<<D_FF_jack_AIC.ave()<<" +- "<<D_FF_jack_AIC.err()<<endl;
      cout<<"Expected error stat(AIC): "<< sqrt( wA_AIC*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB_AIC*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"Expected error syst(AIC): "<< syst_AIC<<endl;
      cout<<"stat+syst(AIC): "<< sqrt( syst_AIC*syst_AIC + wA_AIC*pow(Fit_A_jack.Get_FF(i).err(ixg-1),2) + wB_AIC*pow(Fit_B_jack.Get_FF(i).err(ixg-1),2))<<endl;
      cout<<"#########################"<<endl;
      
      
    }


    //print form factors

    Print_To_File({}, {Bs_xg_t_list, FF_jack.ave(), FF_jack.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.dat", "", "");
    //Print_To_File({}, {Bs_xg_t_list, FF_boot.ave(), FF_boot.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.dat", "", "");
    
    Print_To_File({}, {Bs_xg_t_list, FF_jack_AIC.ave(), FF_jack_AIC.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.dat", "", "");
    //Print_To_File({}, {Bs_xg_t_list, FF_boot_AIC.ave(), FF_boot_AIC.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.dat", "", "");


    //print covariance matrix

    //Eigen::MatrixXd Cov_boot(nxg-1, nxg-1);
    Eigen::MatrixXd Cov_jack(nxg-1, nxg-1);
    //Eigen::MatrixXd Cov_boot_AIC(nxg-1, nxg-1);
    Eigen::MatrixXd Cov_jack_AIC(nxg-1, nxg-1);

    //Eigen::MatrixXd Corr_boot(nxg-1, nxg-1);
    Eigen::MatrixXd Corr_jack(nxg-1, nxg-1);
    //Eigen::MatrixXd Corr_boot_AIC(nxg-1, nxg-1);
    Eigen::MatrixXd Corr_jack_AIC(nxg-1, nxg-1);
 
    for(int x_xg=1; x_xg<nxg;x_xg++) {
      for(int y_xg=1; y_xg<nxg;y_xg++) {

      
	//Cov_boot(x_xg-1, y_xg-1) = FF_boot.distr_list[x_xg-1]%FF_boot.distr_list[y_xg-1];
	Cov_jack(x_xg-1, y_xg-1) = FF_jack.distr_list[x_xg-1]%FF_jack.distr_list[y_xg-1];
	//Cov_boot_AIC(x_xg-1, y_xg-1) = FF_boot_AIC.distr_list[x_xg-1]%FF_boot_AIC.distr_list[y_xg-1];
	Cov_jack_AIC(x_xg-1, y_xg-1) = FF_jack_AIC.distr_list[x_xg-1]%FF_jack_AIC.distr_list[y_xg-1];

      
	//Corr_boot(x_xg-1, y_xg-1) = Cov_boot(x_xg-1,y_xg-1)/(FF_boot.err(x_xg-1)*FF_boot.err(y_xg-1));
	Corr_jack(x_xg-1, y_xg-1) = Cov_jack(x_xg-1,y_xg-1)/(FF_jack.err(x_xg-1)*FF_jack.err(y_xg-1));
	//Corr_boot_AIC(x_xg-1, y_xg-1) = Cov_boot_AIC(x_xg-1,y_xg-1)/(FF_boot_AIC.err(x_xg-1)*FF_boot_AIC.err(y_xg-1));
	Corr_jack_AIC(x_xg-1, y_xg-1) = Cov_jack_AIC(x_xg-1,y_xg-1)/(FF_jack_AIC.err(x_xg-1)*FF_jack_AIC.err(y_xg-1));
      }
    }

    //Print To File
  
    //ofstream Print_Cov_boot("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.cov");
    ofstream Print_Cov_jack("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.cov");
    //ofstream Print_Cov_boot_AIC("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.cov");
    ofstream Print_Cov_jack_AIC("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.cov");

    //ofstream Print_Corr_boot("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_boot.corr");
    ofstream Print_Corr_jack("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_jack.corr");
    //ofstream Print_Corr_boot_AIC("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_boot.corr");
    ofstream Print_Corr_jack_AIC("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_jack.corr");


    //Print_Cov_boot<<Cov_boot<<endl;
    //Print_Cov_boot.close();
    Print_Cov_jack<<Cov_jack<<endl;
    Print_Cov_jack.close();
    //Print_Cov_boot_AIC<<Cov_boot_AIC<<endl;
    //Print_Cov_boot_AIC.close();
    Print_Cov_jack_AIC<<Cov_jack_AIC<<endl;
    Print_Cov_jack_AIC.close();

  
    //Print_Corr_boot<<Corr_boot<<endl;
    //Print_Corr_boot.close();
    Print_Corr_jack<<Corr_jack<<endl;
    Print_Corr_jack.close();
    //Print_Corr_boot_AIC<<Corr_boot_AIC<<endl;
    //Print_Corr_boot_AIC.close();
    Print_Corr_jack_AIC<<Corr_jack_AIC<<endl;
    Print_Corr_jack_AIC.close();


    
    //interpolate the form factors

    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FF_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FF_AIC_cont_interpol_jacks;

    //interpolate continuum extrapolated form factors
    for(int ijack=0;ijack<NJ;ijack++) {

      Vfloat FF_ij, FF_AIC_ij;
  
      for(int ixg=1;ixg<nxg;ixg++) {
	FF_ij.push_back( FF_jack.distr_list[ixg-1].distr[ijack]);
	FF_AIC_ij.push_back( FF_jack_AIC.distr_list[ixg-1].distr[ijack]);
      }

      FF_cont_interpol_jacks.emplace_back( FF_ij.begin(), FF_ij.end(), 0.1, 0.1);
      FF_AIC_cont_interpol_jacks.emplace_back( FF_AIC_ij.begin(), FF_AIC_ij.end(), 0.1, 0.1);

    }

    //push_back
    F_interpolated.push_back( FF_AIC_cont_interpol_jacks);


    auto FF_cont_interpol_distr= [&FF_cont_interpol_jacks](double xg) -> distr_t {
      distr_t return_distr(1);
      for(int ijack=0; ijack<NJ;ijack++) { return_distr.distr.push_back( FF_cont_interpol_jacks[ijack](xg));}
      return return_distr;
    };

    auto FF_AIC_cont_interpol_distr= [&FF_AIC_cont_interpol_jacks](double xg) -> distr_t {
      distr_t return_distr(1);
      for(int ijack=0; ijack<NJ;ijack++) { return_distr.distr.push_back( FF_AIC_cont_interpol_jacks[ijack](xg));}
      return return_distr;
    };
  
 
    distr_t_list FF_interpol_cont_to_print_distr(1);
    distr_t_list FF_AIC_interpol_cont_to_print_distr(1);

  
  
    for(auto &X: Bs_xg_to_spline) {
      FF_interpol_cont_to_print_distr.distr_list.push_back( FF_cont_interpol_distr(X));
      FF_AIC_interpol_cont_to_print_distr.distr_list.push_back( FF_AIC_cont_interpol_distr(X));
    }

    Print_To_File({}, {Bs_xg_to_spline, FF_interpol_cont_to_print_distr.ave(), FF_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/uniform/"+contribs[i]+"_interpol.dat", "", "#xg "+contribs[i]+" "+contribs[i]+"_err");

    Print_To_File({}, {Bs_xg_to_spline, FF_AIC_interpol_cont_to_print_distr.ave(), FF_AIC_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_interpol.dat", "", "#xg "+contribs[i]+" "+contribs[i]+"_err");



 

    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################
    //################                          #################
    //################                          #################
    //################      POLE-FITS           #################
    //################                          #################
    //################                          #################
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################

    class ipar_VMD {
    public:
      ipar_VMD() : FF(0.0), FF_err(0.0) {}
      double FF, FF_err, xg;
    };
  
    class fpar_VMD {
    public:
      fpar_VMD() {}
      fpar_VMD(const Vfloat &par) {
	if((signed)par.size() != 5) crash("In class fpar_VMD  class constructor Vfloat par has size != 5");
	A=par[0];
	M=par[1];
	A2=par[2];
	M2=par[3];
	B=par[4];
      }
      double A, M, A2, M2,B;
    };

  
    //init bootstrap fit
    bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD(NJ);
    bf_VMD.set_warmup_lev(0); //sets warmup
    bf_VMD.Set_number_of_measurements(Bs_xg_t_list.size());
    bf_VMD.Set_verbosity(1);
  
    //################################################
    //define initial condition for params
    double A_guess, M_guess, M2_guess;
    if( contribs[i].substr(1,1) == "V") { M_guess=1.073; M2_guess= 1.38 ; }
    else if( contribs[i].substr(1,1) == "A") { M_guess=1.25; M2_guess= 1.29 ;}
    else crash("FF: "+contribs[i]+" is not V or A");
    A_guess= FF_jack_AIC.ave(0)*M_guess*(M_guess-1.0);
    //###############################################

  
    bf_VMD.Add_par("A", A_guess, A_guess/10);
    bf_VMD.Add_par("M", M_guess, M_guess/10);
    bf_VMD.Add_par("A2", A_guess, A_guess/10);
    bf_VMD.Add_par("M2", M2_guess, M2_guess/10);
    bf_VMD.Add_par("B", A_guess, A_guess/10);
    //fit on mean values to get ch2
    bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD_ch2(1);
    bf_VMD_ch2.set_warmup_lev(0); //sets warmup
    bf_VMD_ch2.Set_number_of_measurements(Bs_xg_t_list.size());
    bf_VMD_ch2.Set_verbosity(1);
    bf_VMD_ch2.Add_par("A", A_guess, A_guess/10);
    bf_VMD_ch2.Add_par("M", M_guess, M_guess/10);
    bf_VMD_ch2.Add_par("A2", A_guess, A_guess/10);
    bf_VMD_ch2.Add_par("M2", M2_guess, M2_guess/10);
    bf_VMD_ch2.Add_par("B", A_guess, A_guess/10);
    //bf_VMD.Fix_par("M",  M_guess);
    //bf_VMD_ch2.Fix_par("M", M_guess);
    bf_VMD.Fix_par("A2", 0.0);
    bf_VMD_ch2.Fix_par("A2", 0.0);
    bf_VMD.Fix_par("M2", M2_guess);
    bf_VMD_ch2.Fix_par("M2", M2_guess);
    //bf_VMD.Fix_par("B",0);
    //bf_VMD_ch2.Fix_par("B",0);

    //ansatz
    bf_VMD.ansatz=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

      double E_res= sqrt( p.M*p.M + ip.xg*ip.xg/4.0);
      double E_res_2= sqrt( p.M2*p.M2 + ip.xg*ip.xg/4.0);
      return ( p.A/( E_res*( E_res - (1.0 - ip.xg/2.0)) )  + p.A2/(E_res_2*( E_res_2 - (1.0 - ip.xg/2.0)) ) + p.B)   ;
    };
    bf_VMD.measurement=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

      return ip.FF;
    };
    bf_VMD.error=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

      return ip.FF_err;
    };

    bf_VMD_ch2.ansatz= bf_VMD.ansatz;
    bf_VMD_ch2.measurement = bf_VMD.measurement;
    bf_VMD_ch2.error = bf_VMD.error;


    //start fitting
    //fill the data
    vector<vector<ipar_VMD>> data_VMD(NJ);
    vector<vector<ipar_VMD>> data_VMD_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_VMD> Bt_fit_VMD;
    boot_fit_data<fpar_VMD> Bt_fit_VMD_ch2;
    for(auto &data_iboot: data_VMD) data_iboot.resize(Bs_xg_t_list.size());
    for(auto &data_iboot: data_VMD_ch2) data_iboot.resize(Bs_xg_t_list.size());
    for(int ijack=0;ijack<NJ;ijack++) {
      for(int ix=0;ix<(signed)Bs_xg_t_list.size();ix++) {
	data_VMD[ijack][ix].FF = FF_jack_AIC.distr_list[ix].distr[ijack];
	data_VMD[ijack][ix].FF_err= FF_jack_AIC.err(ix);
	data_VMD[ijack][ix].xg= Bs_xg_t_list[ix];
	if(ijack==0) {
	  data_VMD_ch2[ijack][ix].FF = FF_jack_AIC.ave(ix);
	  data_VMD_ch2[ijack][ix].FF_err= FF_jack_AIC.err(ix);
	  data_VMD_ch2[ijack][ix].xg= Bs_xg_t_list[ix];

	}
      }
    }

    //append
    bf_VMD.Append_to_input_par(data_VMD);
    bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
    //fit
    cout<<"Fitting "<<contribs[i]<<" using VMD ansatz"<<endl;
    Bt_fit_VMD= bf_VMD.Perform_bootstrap_fit();
    Bt_fit_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();
    double ch2_red_VMD= Bt_fit_VMD_ch2.get_ch2_ave()/( Bs_xg_t_list.size() -3.0);

    //retrieve params
    distr_t Ampl_F(1), pole_F(1), Ampl_F2(1), pole_F2(1), Bg_F(1);
    for(int ijack=0;ijack<NJ;ijack++) {
      Ampl_F.distr.push_back( Bt_fit_VMD.par[ijack].A); pole_F.distr.push_back( Bt_fit_VMD.par[ijack].M);
      Ampl_F2.distr.push_back( Bt_fit_VMD.par[ijack].A2); pole_F2.distr.push_back( Bt_fit_VMD.par[ijack].M2);
      Bg_F.distr.push_back( Bt_fit_VMD.par[ijack].B);

    }


    distr_t_list F_VMD_fit(1);
    //plot fit function
    for(auto &X: Bs_xg_to_spline_VMD) {
      ipar_VMD pp_VMD;
      pp_VMD.xg=X;
      distr_t F_VMD_xg(1);
      for(int ijack=0;ijack<NJ;ijack++) {
	F_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_VMD.par[ijack], pp_VMD));
      }

      F_VMD_fit.distr_list.push_back( F_VMD_xg);
    }

    //print 
    string header_F= "Ampl1: "+to_string_with_precision(Ampl_F.ave(),5)+" +- "+to_string_with_precision(Ampl_F.err(),5)+" M^res1/Mp: "+to_string_with_precision(pole_F.ave(), 5)+" +- "+to_string_with_precision(pole_F.err(), 5)+" Ampl2: "+to_string_with_precision(Ampl_F2.ave(),5)+" +- "+to_string_with_precision(Ampl_F2.err(),5)+" M^res2/Mp: "+to_string_with_precision(pole_F2.ave(), 5)+" +- "+to_string_with_precision(pole_F2.err(), 5)+" B: "+to_string_with_precision(Bg_F.ave(),5)+" +- "+to_string_with_precision(Bg_F.err(),5)+"   ch2/dof: "+to_string_with_precision(ch2_red_VMD  ,5)+" dof: "+to_string(Bs_xg_t_list.size()-2);
  
    Print_To_File({},{ Bs_xg_to_spline_VMD, F_VMD_fit.ave(), F_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/"+contribs[i]+"_VMD.fit" , "", header_F);

    //push_back lambda function

    Coupling_1_list.distr_list.push_back(Ampl_F);
    Pole_1_list.distr_list.push_back(pole_F);
    Coupling_2_list.distr_list.push_back(Ampl_F2);
    Pole_2_list.distr_list.push_back(pole_F2);
    Bg_list.distr_list.push_back(Bg_F);
  
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################
    //###########################################################

    //push_back_analysis
    FF_final.push_back( FF_jack_AIC);

  
  }

  //print ratio between tensor and vector form factors
  distr_t_list R_A_u= FF_final[1]/FF_final[7];
  distr_t_list R_A_d= -1.0*FF_final[2]/FF_final[8];
  distr_t_list R_V_u= FF_final[4]/FF_final[10];
  distr_t_list R_V_d= FF_final[5]/FF_final[11];

  Print_To_File( {}, {Bs_xg_t_list, R_A_u.ave(), R_A_u.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/ratio_Au_jack.dat"  ,   "", "");
  Print_To_File( {}, {Bs_xg_t_list, R_A_d.ave(), R_A_d.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/ratio_Ad_jack.dat"  ,   "", "");
  Print_To_File( {}, {Bs_xg_t_list, R_V_u.ave(), R_V_u.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/ratio_Vu_jack.dat"  ,   "", "");
  Print_To_File( {}, {Bs_xg_t_list, R_V_d.ave(), R_V_d.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/ratio_Vd_jack.dat"  ,   "", "");

  //mass ratio mb/mc
  double mc_to_mb= 3.29;


  //define rescaled form factors at B mass

  distr_t_list FA_Bs =    pow(mc_to_mb,-1.5)*FF_final[1]  - pow(mc_to_mb,-0.5)*FF_final[2];
  distr_t_list FV_Bs =    pow(mc_to_mb,-1.5)*FF_final[4]  + pow(mc_to_mb,-0.5)*FF_final[5];
  distr_t_list FA_T_Bs =  pow(mc_to_mb,-1.5)*FF_final[7]  + pow(mc_to_mb,-0.5)*FF_final[8];
  distr_t_list FV_T_Bs =  pow(mc_to_mb,-1.5)*FF_final[10] + pow(mc_to_mb,-0.5)*FF_final[11];

  distr_t_list FA_BDs =    FF_final[1]  - FF_final[2];
  distr_t_list FV_BDs =    FF_final[4]  + FF_final[5];
  distr_t_list FA_T_BDs =  FF_final[7]  + FF_final[8];
  distr_t_list FV_T_BDs =  FF_final[10] + FF_final[11];

  //print
  Print_To_File({}, {Bs_xg_t_list, FA_Bs.ave(), FA_Bs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_Bs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FV_Bs.ave(), FV_Bs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_Bs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FA_T_Bs.ave(), FA_T_Bs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_T_Bs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FV_T_Bs.ave(), FV_T_Bs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_T_Bs_jack.dat", "", "");

  Print_To_File({}, {Bs_xg_t_list, FA_BDs.ave(), FA_BDs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_BDs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FV_BDs.ave(), FV_BDs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_BDs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FA_T_BDs.ave(), FA_T_BDs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_T_BDs_jack.dat", "", "");
  Print_To_File({}, {Bs_xg_t_list, FV_T_BDs.ave(), FV_T_BDs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_T_BDs_jack.dat", "", "");

   

  auto pole_ansatz= [&](int i, int jack,  double x) {

    double M= Pole_1_list.distr_list[i].distr[jack];
    double A= Coupling_1_list.distr_list[i].distr[jack];
    double E_res= sqrt( M*M +  x*x/4.0);
    double res_1 = A/(E_res*(E_res -(1.0 - x/2)));

    double M2= Pole_2_list.distr_list[i].distr[jack];
    double A2= Coupling_2_list.distr_list[i].distr[jack];
    double E_res_2= sqrt( M2*M2 +  x*x/4.0);
    double res_2 = A2/(E_res_2*(E_res_2 -(1.0 - x/2)));

    double bg = Bg_list.distr_list[i].distr[jack];
     
    return res_1 + res_2 + bg;

  };

  distr_t_list FA_Bs_interpolated(true, Bs_xg_to_spline.size()), FV_Bs_interpolated(true, Bs_xg_to_spline.size()), FA_T_Bs_interpolated(true, Bs_xg_to_spline.size()), FV_T_Bs_interpolated(true, Bs_xg_to_spline.size());
  distr_t_list FA_BDs_interpolated(true, Bs_xg_to_spline.size()), FV_BDs_interpolated(true, Bs_xg_to_spline.size()), FA_T_BDs_interpolated(true, Bs_xg_to_spline.size()), FV_T_BDs_interpolated(true, Bs_xg_to_spline.size());


    
  for(int ix=0; ix < (signed)Bs_xg_to_spline.size(); ix++) {
    double x= Bs_xg_to_spline[ix];
    for(int ijack=0; ijack<NJ;ijack++) {
      //Bs
      FA_Bs_interpolated.distr_list[ix].distr.push_back(  pow(mc_to_mb, -1.5)*pole_ansatz(1,ijack,x) -pow(mc_to_mb,-0.5)*pole_ansatz(2,ijack,x)  ) ;
      FV_Bs_interpolated.distr_list[ix].distr.push_back( pow(mc_to_mb,-1.5)*pole_ansatz(4,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(5,ijack,x) ) ;
      FA_T_Bs_interpolated.distr_list[ix].distr.push_back(  pow(mc_to_mb, -1.5)*pole_ansatz(7,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(8,ijack,x) ) ;
      FV_T_Bs_interpolated.distr_list[ix].distr.push_back(  pow(mc_to_mb,-1.5)*pole_ansatz(10,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(11,ijack,x) ) ;
      //BDs
      FA_BDs_interpolated.distr_list[ix].distr.push_back(  pole_ansatz(1,ijack,x) -pole_ansatz(2,ijack,x)  ) ;
      FV_BDs_interpolated.distr_list[ix].distr.push_back( pole_ansatz(4,ijack,x) + pole_ansatz(5,ijack,x) ) ;
      FA_T_BDs_interpolated.distr_list[ix].distr.push_back(  pole_ansatz(7,ijack,x) + pole_ansatz(8,ijack,x) ) ;
      FV_T_BDs_interpolated.distr_list[ix].distr.push_back(  pole_ansatz(10,ijack,x) + pole_ansatz(11,ijack,x) ) ;
    }
  }

  //Bs
  Print_To_File( {}, {Bs_xg_to_spline, FA_Bs_interpolated.ave(), FA_Bs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_Bs_interpol.dat", "", "#xg FA_Bs  FA_Bs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FV_Bs_interpolated.ave(), FV_Bs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_Bs_interpol.dat", "", "#xg FV_Bs  FV_Bs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FA_T_Bs_interpolated.ave(), FA_T_Bs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_T_Bs_interpol.dat", "", "#xg FA_T_Bs  FA_T_Bs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FV_T_Bs_interpolated.ave(), FV_T_Bs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_T_Bs_interpol.dat", "", "#xg FV_T_Bs  FV_T_Bs_err");
  //BDs
  Print_To_File( {}, {Bs_xg_to_spline, FA_BDs_interpolated.ave(), FA_BDs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_BDs_interpol.dat", "", "#xg FA_BDs  FA_BDs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FV_BDs_interpolated.ave(), FV_BDs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_BDs_interpol.dat", "", "#xg FV_BDs  FV_BDs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FA_T_BDs_interpolated.ave(), FA_T_BDs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FA_T_BDs_interpol.dat", "", "#xg FA_T_BDs  FA_T_BDs_err");
  Print_To_File( {}, {Bs_xg_to_spline, FV_T_BDs_interpolated.ave(), FV_T_BDs_interpolated.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/BMA/AIC/FV_T_BDs_interpol.dat", "", "#xg FV_T_BDs  FV_T_BDs_err");
    


  auto fnull= [](double x) { return 0;};
 
   

  int Nxgs=201;
  distr_t_list rate_Bs(true, Nxgs), rate_Bs_diff(true,Nxgs), AFB_Bs(true, Nxgs);
  distr_t_list rate_BDs(true, Nxgs), rate_BDs_diff(true,Nxgs), AFB_BDs(true, Nxgs);
  vector<double> xg_min_list(Nxgs), xg_max_list(Nxgs), qmin_list(Nxgs);
  double Mb=5.36692;
  for(int it=0;it<Nxgs;it++) { xg_min_list[it] = 0.0; xg_max_list[it] =  1 - pow( 3.8 + it*(Mb-3.8)/(Nxgs-1.0)  ,2)/pow(Mb,2); qmin_list[it] = Mb*sqrt(1-xg_max_list[it]);  }
   
  double evolutor_ZT = evolutor_ZT_MS_bar(2.0,5.0);

  cout<<"ZT(5GeV, MSbar)/ZT(2GeV, MSbar): "<<evolutor_ZT<<endl;
  cout<<"Riemann zeta(3): "<<riemann_zeta(3)<<endl;

  distr_t g_Bs_phi(1);
  for(int i=0;i<NJ;i++) g_Bs_phi.distr.push_back( 0.27 + 0.27*GM()/sqrt(NJ-1.0));

  
   
  for(int Is_BDs=0; Is_BDs <= 1; Is_BDs++ ) {

    mc_to_mb= Is_BDs?1.0:3.29;

    for(int it=0; it < Nxgs;it++) {

      for(int ijack=0;ijack<NJ;ijack++) {


	 
	auto FA = [&](double x) {
	  return pow(mc_to_mb, -1.5)*pole_ansatz(1,ijack,x) -pow(mc_to_mb,-0.5)*pole_ansatz(2,ijack,x);
	};
	
	auto FV = [&](double x) {
	  return pow(mc_to_mb,-1.5)*pole_ansatz(4,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(5,ijack,x);
	};


	
	 
	auto FA_T = [&](double x) {

	  double Mbs=5.36692;
	  double q= sqrt((1-x))*Mbs;
	  double mphi=1.019461;
	  double f_phi=0.241;
	  double Gamma_phi= 0.00425;
	  double mq= mphi/q;
	  double gq= Gamma_phi/q;
	  
	  double real=  (pow(mc_to_mb, -1.5)*pole_ansatz(7,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(8,ijack,x))*evolutor_ZT;
	  double virt=  (pow(mc_to_mb, -1.5)*pole_ansatz(7,ijack,1.0) + pow(mc_to_mb,-0.5)*pole_ansatz(8,ijack,1.0))*evolutor_ZT;
	  
	  double phi_pole_re = -(2.0/3.0)*g_Bs_phi.distr[ijack]*(f_phi/mphi)*(1-mq*mq)/( pow(1.0 - mq*mq,2) + pow(mq*gq,2) );
	  
	  return real+virt+phi_pole_re;
	};

	auto FV_T = [&](double x) {

	  double Mbs=5.36692;
	  double q= sqrt((1-x))*Mbs;
	  double mphi=1.019461;
	  double f_phi=0.241;
	  double Gamma_phi= 0.00425;
	  double mq= mphi/q;
	  double gq= Gamma_phi/q;
	  
	  double real = (pow(mc_to_mb,-1.5)*pole_ansatz(10,ijack,x) + pow(mc_to_mb,-0.5)*pole_ansatz(11,ijack,x))*evolutor_ZT;
	  double virt=  (pow(mc_to_mb, -1.5)*pole_ansatz(10,ijack,1.0) + pow(mc_to_mb,-0.5)*pole_ansatz(11,ijack,1.0))*evolutor_ZT;
	  
	  double phi_pole_re = -(2.0/3.0)*g_Bs_phi.distr[ijack]*(f_phi/mphi)*(1-mq*mq)/( pow(1.0 - mq*mq,2) + pow(mq*gq,2) );
	  
	  return real+virt+phi_pole_re;
	  
	};
	
	auto FA_T_IM = [&](double x) {
	  
	  double Mbs=5.36692;
	  double q= sqrt((1-x))*Mbs;
	  double mphi=1.019461;
	  double f_phi=0.241;
	  double Gamma_phi= 0.00425;
	  double mq= mphi/q;
	  double gq= Gamma_phi/q;
	  
	  double pr= (2.0/3.0)*g_Bs_phi.distr[ijack]*(f_phi/mphi)*mq*gq/( pow(1.0 - mq*mq,2) + pow(mq*gq,2) ) ; 
	  return pr;
	};
	
	auto FV_T_IM = [&](double x) {
	  
	  double Mbs=5.36692;
	  double q= sqrt((1-x))*Mbs;
	  double mphi=1.019461;
	  double f_phi=0.241;
	  double Gamma_phi= 0.00425;
	  double mq= mphi/q;
	  double gq= Gamma_phi/q;
	  
	  double pr= (2.0/3.0)*g_Bs_phi.distr[ijack]*(f_phi/mphi)*mq*gq/( pow(1.0 - mq*mq,2) + pow(mq*gq,2) ) ; 
	  return pr;
	  
	};

	
	if(!Is_BDs) {
	  rate_Bs.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate( ijack, NJ, FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM , xg_max_list[it], xg_min_list[it]  ));
	  rate_Bs_diff.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(ijack, NJ, xg_max_list[it], FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM ));
	  AFB_Bs.distr_list[it].distr.push_back( Compute_AFB(ijack, NJ, xg_max_list[it], FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM ));
	}
	else {
	  rate_BDs.distr_list[it].distr.push_back(  Compute_Bs_mumugamma_decay_rate( ijack, NJ, FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM , xg_max_list[it], xg_min_list[it]  ));
	  rate_BDs_diff.distr_list[it].distr.push_back( Compute_Bs_mumugamma_differential_decay_rate(ijack, NJ, xg_max_list[it], FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM ));
	  AFB_BDs.distr_list[it].distr.push_back( Compute_AFB(ijack, NJ, xg_max_list[it], FV, FA, FV_T, FA_T, FV_T_IM, FA_T_IM ));
	}
      }
    }
  }
  

  //Bs
  Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_Bs.ave(), rate_Bs.err(), rate_Bs_diff.ave(), rate_Bs_diff.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/rate_Bs.dat", "", "#qmin[GeV] xg_min xg_max  rate   rate_diff    ");
  Print_To_File( {}, { qmin_list, xg_max_list, AFB_Bs.ave(), AFB_Bs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/AFB_Bs.dat", "", "#q[GeV] xg   AFG    ");
  //BDs
  Print_To_File( {}, { qmin_list, xg_min_list, xg_max_list, rate_BDs.ave(), rate_BDs.err(), rate_BDs_diff.ave(), rate_BDs_diff.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/rate_BDs.dat", "", "#qmin[GeV] xg_min xg_max  rate   rate_diff    ");
  Print_To_File( {}, { qmin_list, xg_max_list, AFB_BDs.ave(), AFB_BDs.err()},  "../data/ph_emission/"+ph_type+"/"+MESON+"/AFB_BDs.dat", "", "#q[GeV] xg   AFG    ");
   
  cout<<"Done! Bye"<<endl;
    
  return;
}


rt_FF_Bs Get_Bs_mumu_gamma_form_factors(int num_xg, int Perform_continuum_extrapolation, bool Use_three_finest, bool Include_a4, bool UseJack, string Fit_tag, string path_list, string out_tag) {


  cout<<"CALLING Get_Bs_mumu_gamma_form_factors"<<endl;

  int size_mu_nu= Is_reph?2:4;

  int size_mu_nu_T=2;

  //return type
  rt_FF_Bs rt_fit(UseJack);
  rt_fit.Use_three_finest=Use_three_finest;
  rt_fit.Include_a4 = Include_a4;
  rt_fit.num_xg=num_xg;

  
  
  string _Fit_tag = "";

  if(Fit_tag.size() > 0) { _Fit_tag= "_"+Fit_tag.substr(0, Fit_tag.size()-1); } 

  const int Njacks=((UseJack==true)?NJ:Nboots);

   
  string ph_type_mes=ph_type+"/"+MESON;
  
  




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

  
    if(A_bis.length() <= 4) return A_bis < B_bis;
     
    string rA = A_bis.substr(A_bis.length()-2);
    string rB = B_bis.substr(B_bis.length()-2);

        
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(rA.substr(1,1));
      int n2 = stoi(rB.substr(1,1));
      if(rA == rB) {
	if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
	else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    return A_bis<B_bis;
  };
  
  
 
  int off_i = (Is_reph?1:0);
  int off_T = 1;
  


  //2pts
  data_t data_2pts, data_2pts_SM, data_2pts_SMSM;
  data_t data_2pts_V1, data_2pts_V2, data_2pts_V3;
  data_t data_2pts_u_V1, data_2pts_u_V2, data_2pts_u_V3;
  data_t data_2pts_d_V1, data_2pts_d_V2, data_2pts_d_V3;
  data_t data_2pts_u_T1, data_2pts_u_T2, data_2pts_u_T3;
  data_t data_2pts_u_T1_SMSM, data_2pts_u_T2_SMSM, data_2pts_u_T3_SMSM;
  data_t data_2pts_u_B1, data_2pts_u_B2, data_2pts_u_B3;

  
  //Read data

  //2pts function

  data_2pts.Read(path_list, "mes_contr_2pts_SM_3", "P5P5", Sort_light_confs);
  data_2pts_SM.Read(path_list, "mes_contr_2pts_SM_3", "P5P5", Sort_light_confs);
  data_2pts_SMSM.Read(path_list, "mes_contr_2pts_SMSM_3", "P5P5", Sort_light_confs);
  data_2pts_V1.Read(path_list, "mes_contr_2pts_SM_3", "V1V1", Sort_light_confs);
  data_2pts_V2.Read(path_list, "mes_contr_2pts_SM_3", "V2V2", Sort_light_confs);
  data_2pts_V3.Read(path_list, "mes_contr_2pts_SM_3", "V3V3", Sort_light_confs);
  data_2pts_d_V1.Read(path_list, "mes_contr_2pts_SM_1", "V1V1", Sort_light_confs);
  data_2pts_d_V2.Read(path_list, "mes_contr_2pts_SM_1", "V2V2", Sort_light_confs);
  data_2pts_d_V3.Read(path_list, "mes_contr_2pts_SM_1", "V3V3", Sort_light_confs);
  data_2pts_u_V1.Read(path_list, "mes_contr_2pts_SM_2", "V1V1" , Sort_light_confs);
  data_2pts_u_V2.Read(path_list, "mes_contr_2pts_SM_2", "V2V2", Sort_light_confs);
  data_2pts_u_V3.Read(path_list, "mes_contr_2pts_SM_2", "V3V3", Sort_light_confs);
  data_2pts_u_T1.Read(path_list, "mes_contr_2pts_SM_2", "T1T1" , Sort_light_confs);
  data_2pts_u_T2.Read(path_list, "mes_contr_2pts_SM_2", "T2T2", Sort_light_confs);
  data_2pts_u_T3.Read(path_list, "mes_contr_2pts_SM_2", "T3T3", Sort_light_confs);
  data_2pts_u_T1_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T1T1" , Sort_light_confs);
  data_2pts_u_T2_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T2T2", Sort_light_confs);
  data_2pts_u_T3_SMSM.Read(path_list, "mes_contr_2pts_SMSM_2", "T3T3", Sort_light_confs);
  data_2pts_u_B1.Read(path_list, "mes_contr_2pts_SM_2", "B1B1" , Sort_light_confs);
  data_2pts_u_B2.Read(path_list, "mes_contr_2pts_SM_2", "B2B2", Sort_light_confs);
  data_2pts_u_B3.Read(path_list, "mes_contr_2pts_SM_2", "B3B3", Sort_light_confs);



  //axial
  vector<vector<vector<data_t>>> C_A_Fu_data(size_mu_nu), C_A_Fd_data(size_mu_nu), C_A_Bu_data(size_mu_nu), C_A_Bd_data(size_mu_nu);
  //vector
  vector<vector<vector<data_t>>> C_V_Fu_data(size_mu_nu), C_V_Fd_data(size_mu_nu), C_V_Bu_data(size_mu_nu), C_V_Bd_data(size_mu_nu);
  //TK
  vector<vector<vector<data_t>>> C_TK_Fu_data(size_mu_nu), C_TK_Fd_data(size_mu_nu), C_TK_Bu_data(size_mu_nu), C_TK_Bd_data(size_mu_nu);
  //BK
  vector<vector<vector<data_t>>> C_BK_Fu_data(size_mu_nu), C_BK_Fd_data(size_mu_nu), C_BK_Bu_data(size_mu_nu), C_BK_Bd_data(size_mu_nu);

  for(int mu=0;mu<size_mu_nu;mu++) {

    C_A_Fu_data[mu].resize(size_mu_nu);
    C_A_Fd_data[mu].resize(size_mu_nu);
    C_A_Bu_data[mu].resize(size_mu_nu);
    C_A_Bd_data[mu].resize(size_mu_nu);
    C_V_Fu_data[mu].resize(size_mu_nu);
    C_V_Fd_data[mu].resize(size_mu_nu);
    C_V_Bu_data[mu].resize(size_mu_nu);
    C_V_Bd_data[mu].resize(size_mu_nu);

    C_BK_Fu_data[mu].resize(size_mu_nu_T);
    C_BK_Fd_data[mu].resize(size_mu_nu_T);
    C_BK_Bu_data[mu].resize(size_mu_nu_T);
    C_BK_Bd_data[mu].resize(size_mu_nu_T);
    C_TK_Fu_data[mu].resize(size_mu_nu_T);
    C_TK_Fd_data[mu].resize(size_mu_nu_T);
    C_TK_Bu_data[mu].resize(size_mu_nu_T);
    C_TK_Bd_data[mu].resize(size_mu_nu_T);

    
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_A_Fu_data[mu][nu].resize(num_xg);
      C_A_Fd_data[mu][nu].resize(num_xg);
      C_A_Bu_data[mu][nu].resize(num_xg);
      C_A_Bd_data[mu][nu].resize(num_xg);
      C_V_Fu_data[mu][nu].resize(num_xg);
      C_V_Fd_data[mu][nu].resize(num_xg);
      C_V_Bu_data[mu][nu].resize(num_xg);
      C_V_Bd_data[mu][nu].resize(num_xg);

     


    }

    for(int nu=0;nu<size_mu_nu_T;nu++) {

      C_BK_Fu_data[mu][nu].resize(num_xg);
      C_BK_Fd_data[mu][nu].resize(num_xg);
      C_BK_Bu_data[mu][nu].resize(num_xg);
      C_BK_Bd_data[mu][nu].resize(num_xg);
      C_TK_Fu_data[mu][nu].resize(num_xg);
      C_TK_Fd_data[mu][nu].resize(num_xg);
      C_TK_Bu_data[mu][nu].resize(num_xg);
      C_TK_Bd_data[mu][nu].resize(num_xg);

    }

  }

  //read data
  for(int ixg=0;ixg<num_xg;ixg++) {
    for(int mu=0;mu<size_mu_nu;mu++) {
      for(int nu=0;nu<size_mu_nu;nu++) {

	//vector
	//Fu
	C_V_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_V_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_V_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_V_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5", Sort_light_confs );

	//axial
	//Fu
	C_A_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Fd
	C_A_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bu
	C_A_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
	//Bd
	C_A_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5", Sort_light_confs );
      }
      for(int nu=0;nu<size_mu_nu_T;nu++) {

	//tensor electric part
	//Fu
	C_TK_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Fd
	C_TK_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bu
	C_TK_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bd
	C_TK_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "T"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//tensor magnetic part
	//Fu
	C_BK_Fu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Fd
	C_BK_Fd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bu
	C_BK_Bu_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
	//Bd
	C_BK_Bd_data[mu][nu][ixg].Read(path_list, "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "B"+to_string(nu+off_T)+"P5", Sort_light_confs );
      }
    }
  }

  //#######################################################################################################################################



 

  int Nens = data_2pts.size;
  GaussianMersenne GM(543543);
  GaussianMersenne GM_T(5432);

  //vectors to store FV and FA for all ensembles and values of gamma
  vector<distr_t_list> FA_per_ens(Nens), FV_per_ens(Nens), FA_u_per_ens(Nens), FV_u_per_ens(Nens), FA_d_per_ens(Nens), FV_d_per_ens(Nens);
  vector<distr_t_list> FA_per_kin(num_xg-1), FV_per_kin(num_xg-1), FA_u_per_kin(num_xg-1), FV_u_per_kin(num_xg-1), FA_d_per_kin(num_xg-1), FV_d_per_kin(num_xg-1); //num_xg -1 is to exclude xg=0 which is undefined
  vector<distr_t_list> FA_T_per_ens(Nens), FV_T_per_ens(Nens), FA_T_u_per_ens(Nens), FV_T_u_per_ens(Nens), FA_T_d_per_ens(Nens), FV_T_d_per_ens(Nens);
  vector<distr_t_list> FA_T_per_kin(num_xg-1), FV_T_per_kin(num_xg-1), FA_T_u_per_kin(num_xg-1), FV_T_u_per_kin(num_xg-1), FA_T_d_per_kin(num_xg-1), FV_T_d_per_kin(num_xg-1);
  distr_t_list a_distr_list(UseJack);
  distr_t_list a_distr_list_red(UseJack);
  Vfloat L_list, L_list_red;
  //M_Ds, F_Ds
  distr_t_list MP_list(UseJack);
  distr_t_list MP_list_red(UseJack);
  distr_t_list FP_list(UseJack);
  distr_t_list ZT_list(UseJack);
  distr_t_list MP_ov_FP_list(UseJack);
  //resize
  distr_t_list FP_3pt(UseJack);

  
  
  
  //define lambda function to combine FF and BB

  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  //resample Jpsi tensor decay constant
  distr_t fT_Jpsi_distr(UseJack);

  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

    fT_Jpsi_distr.distr.push_back( fT_Jpsi + fT_Jpsi_err*GM_T()/((UseJack==true)?sqrt(Njacks-1.0):1.0));

    
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM_T()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM_T()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM_T()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM_T()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

  }

 
  for(int iens=0;iens<Nens;iens++) {

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/C/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/H/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results");

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/C_T/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/H_T/"+data_2pts.Tag[iens]);
      
    
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/mass/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/decay_const/"+data_2pts.Tag[iens]);
    
    cout<<"Analyzing ensemble: "<<data_2pts.Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts.Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = data_2pts.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    
    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d, virts;

    thetas= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 1 , 5);
    virts=  Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 2 , 5);
    masses_u= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File(path_list+"/"+data_2pts.Tag[iens]+"/pars_list.dat", 4 , 5);

    cout<<"pars_list.dat: Read!"<<endl;

    //if((signed)thetas.size() != num_xg) crash("Number of rows in pars_list.dat does not match num_xg"); 

    int num_xg_per_ens= (signed)thetas.size();
       
    //RCs
    distr_t Za, Zv, Z_T, a_distr;
    if(data_2pts.Tag[iens].substr(1,1)=="A") {      Za= ZA_A; Zv=ZV_A; Z_T=ZT_A; a_distr=a_A;}
    else if(data_2pts.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; Z_T=ZT_B; a_distr=a_B;}
    else if(data_2pts.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; Z_T=ZT_C; a_distr=a_C;}
    else if(data_2pts.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; Z_T=ZT_D; a_distr=a_D;}
    else crash("Ensemble: "+data_2pts.Tag[iens]+" not recognised");


    a_distr_list.distr_list.push_back(a_distr);
    L_list.push_back(L_info.L);
    
    if(data_2pts.Tag[iens] != "cB211b.072.96") {
      a_distr_list_red.distr_list.push_back(a_distr);
      L_list_red.push_back(L_info.L);
    }
   
    

    //read masses
    double mu= masses_u[0];
    double md= masses_d[0];

    
    cout<<"ZA: "<<Za.ave()<<" +- "<<Za.err()<<endl;
    cout<<"ZV: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
    cout<<"mu: "<<mu<<endl;
    cout<<"md: "<<md<<endl;
    cout<<"Number of kinematic present: "<<num_xg_per_ens<<endl;
    cout<<"Number of kinematic analyzed: "<<num_xg<<endl;
    
    
   
    
    //set time interval for eff_mass_fit SM
    if(data_2pts.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=14; Corr.Tmax= 24;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 24;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=20;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=19;}
      else { Corr.Tmin=14; Corr.Tmax=27;}
      
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=24; Corr.Tmax=38;}
      else if(MESON == "B2s") {Corr.Tmin=22; Corr.Tmax=34;}
      else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=23;}
      else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=22;}
      else {Corr.Tmin=23; Corr.Tmax=33;}
    }
    else if(data_2pts.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON == "B2s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=31;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=30;}
      else {Corr.Tmin=29; Corr.Tmax=42;}
    }
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=32; Corr.Tmax=53;}
      else if(MESON == "B2s") { Corr.Tmin=28; Corr.Tmax= 37;}
      else if(MESON=="B3s") {Corr.Tmin=28; Corr.Tmax=37;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=37;}
      else {Corr.Tmin=38; Corr.Tmax=52;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SM.dat");
    distr_t_list pt2_distr_SMSM= Corr.corr_t(data_2pts_SMSM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SMSM.dat");
    distr_t_list eff_mass_SMSM = Corr.effective_mass_t(pt2_distr_SMSM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SMSM.dat");
    distr_t mel_SMSM= Corr.Fit_distr(Corr.mel_ov_mass_t(pt2_distr_SMSM, ""))/2.0;


    distr_t_list pt2_T_distr= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_T1.col(0)[iens], data_2pts_u_T2.col(0)[iens], data_2pts_u_T3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/corr_T_2pt.dat");
    distr_t_list pt2_T_distr_SMSM= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_T1_SMSM.col(0)[iens], data_2pts_u_T2_SMSM.col(0)[iens], data_2pts_u_T3_SMSM.col(0)[iens]), "");
    distr_t_list pt2_B_distr= -1*(1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_B1.col(0)[iens], data_2pts_u_B2.col(0)[iens], data_2pts_u_B3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/corr_B_2pt.dat");
    distr_t_list eff_mass_T= Corr.effective_mass_t(pt2_T_distr,"../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_T.dat");
    distr_t_list eff_mass_B= Corr.effective_mass_t(pt2_B_distr,"../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_B.dat"); 
    distr_t_list ft_distr = Corr.residue_t( pt2_T_distr, "")/(Corr.Fit_distr(eff_mass_T)*Corr.matrix_element_t(pt2_T_distr_SMSM, ""));

  


   

    
    
    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    distr_t M_P_SMSM = Corr.Fit_distr(eff_mass_SMSM);
    distr_t M_P= 0.5*(M_P_SM+M_P_SMSM);
    auto SINH= [](double x) { return sinh(x);};
    distr_t_list FP_SM_distr_list = (mu + md )*Corr.residue_t( pt2_distr_SM, "")/(M_P*distr_t::f_of_distr(SINH, M_P)*Corr.matrix_element_t(pt2_distr_SMSM, ""));
    Print_To_File({}, {FP_SM_distr_list.ave(), FP_SM_distr_list.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const_SM.dat.t", "", "");
    distr_t FP_SM= Corr.Fit_distr( FP_SM_distr_list  );

    
    distr_t F_P= FP_SM;
    distr_t F_T=Corr.Fit_distr(ft_distr);
       
    MP_list.distr_list.push_back(M_P/a_distr);
    FP_list.distr_list.push_back(F_P/a_distr);
    ZT_list.distr_list.push_back(Z_T);
    MP_ov_FP_list.distr_list.push_back( M_P/F_P);
    if(data_2pts.Tag[iens] != "cB211b.072.96") {
      MP_list_red.distr_list.push_back(M_P/a_distr);
    }

    distr_t ZZT = fT_Jpsi_distr*a_distr/F_T;
   
    distr_t RF= Get_id_jack_distr(NJ); // (a_distr/FP_SM);

    cout<<"aM_P: "<<M_P.ave()<<" +- "<<M_P.err()<<" -> "<< (M_P/(a_distr)).ave()<<" +- "<<(M_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P: "<<F_P.ave()<<" +- "<<F_P.err()<<" -> "<< (F_P/(a_distr)).ave()<<" +- "<<(F_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P(SM): " << FP_SM.ave() << " +- " << FP_SM.err() << " -> "<< (FP_SM/(a_distr)).ave()<<" +- "<<(FP_SM/a_distr).err()<<" GeV"<< endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"aM_P(SM): "<<M_P_SM.ave()<<" +- "<<M_P_SM.err()<<" -> "<< (M_P_SM/(a_distr)).ave()<<" +- "<<(M_P_SM/a_distr).err()<<" GeV"<<endl;
    cout<<"aM_P(SMSM): "<<M_P_SMSM.ave()<<" +- "<<M_P_SMSM.err()<<" -> "<<(M_P_SMSM/a_distr).ave()<<" +- "<<(M_P_SMSM/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"Ensemble: "<<data_2pts.Tag[iens]<<" ZT(ETMC): "<<Z_T.ave()<<" +- "<<Z_T.err()<<" from F_T: "<<ZZT.ave()<<" +- "<<ZZT.err()<<endl;


    //vector ud meson
    distr_t_list pt2_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_V1.col(0)[iens], data_2pts_V2.col(0)[iens], data_2pts_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_ud.dat");
    distr_t_list eff_mass_V = Corr.effective_mass_t(pt2_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_ud.dat");
    //vector u meson
    distr_t_list uu_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_u_V1.col(0)[iens], data_2pts_u_V2.col(0)[iens], data_2pts_u_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_uu.dat");
    distr_t_list eff_mass_V_uu= Corr.effective_mass_t(uu_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_uu.dat");
    //vector d meson
    distr_t_list dd_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_d_V1.col(0)[iens], data_2pts_d_V2.col(0)[iens], data_2pts_d_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V_dd.dat");
    distr_t_list eff_mass_V_dd= Corr.effective_mass_t(dd_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V_dd.dat");

    //#######################
    
    //#######################


    //define meson mass exponential to be removed
    auto EXP_MES_FUNC = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};

    distr_t_list EXP_MES= distr_t_list::f_of_distr( EXP_MES_FUNC, M_P, Corr.Nt);

    //for FV and FA
    vector<vector<vector<distr_t_list>>> Ax_glb, Vec_glb, Ax_u_glb, Ax_d_glb, Vec_u_glb, Vec_d_glb, Ax_Ds_glb, Ax_Ds_u_glb, Ax_Ds_d_glb;
    distr_t_list FV(UseJack), FA(UseJack), FV_u(UseJack), FV_d(UseJack), FA_u(UseJack), FA_d(UseJack);
    //for FV_T and FA_T
    vector<vector<vector<distr_t_list>>> Ax_T_glb, Vec_T_glb, Ax_T_u_glb, Ax_T_d_glb, Vec_T_u_glb, Vec_T_d_glb;
    distr_t_list FV_T(UseJack), FA_T(UseJack), FV_T_u(UseJack), FV_T_d(UseJack), FA_T_u(UseJack), FA_T_d(UseJack);


    distr_t_list xg_list(UseJack);

    //for FV and FA
    Vfloat ax_Tmin, ax_Tmax;
    Vfloat vec_Tmin, vec_Tmax;
    Vfloat ax_u_Tmin, ax_u_Tmax;
    Vfloat vec_u_Tmin, vec_u_Tmax;
    Vfloat ax_d_Tmin, ax_d_Tmax;
    Vfloat vec_d_Tmin, vec_d_Tmax;


    //for FV_T and FA_T
    Vfloat ax_T_Tmin, ax_T_Tmax;
    Vfloat vec_T_Tmin, vec_T_Tmax;
    Vfloat ax_T_u_Tmin, ax_T_u_Tmax;
    Vfloat vec_T_u_Tmin, vec_T_u_Tmax;
    Vfloat ax_T_d_Tmin, ax_T_d_Tmax;
    Vfloat vec_T_d_Tmin, vec_T_d_Tmax;
    

    for(int ixg=0;ixg<num_xg;ixg++) {

      //for FV and FA
      vector<vector<distr_t_list>> Ax_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_d(size_mu_nu);
      //for FV_T and FA_T
      vector<vector<distr_t_list>> Ax_T_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_T_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_T_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_T_tens_d(size_mu_nu);
      //for FA of the Ds meson
      vector<vector<distr_t_list>> Ax_Ds_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_Ds_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_Ds_tens_d(size_mu_nu);


      //get xg, Eg, kz from thetas

      double theta=thetas[ixg];

      pt3_momenta pt3_mom(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], virts[ixg], L_info.L, L_info.T);

      double Eg= pt3_mom.Egamma();
      distr_t xg= pt3_mom.x_gamma(M_P);
      if(xg.ave() > 1e-10 ) xg_list.distr_list.push_back(xg);
      double kz = pt3_mom.k()[2];
    
   

      //define photon exponential to be removed
      Vfloat EXP_PH(Corr.Nt,0.0);
      Vfloat EXP_PH_RED(Corr.Nt, 0.0);
      for(int t=0; t < Corr.Nt;t++) EXP_PH[t] = exp( Eg*abs( Corr.Nt/2 - t));
      for(int t=0; t < Corr.Nt;t++) EXP_PH_RED[t] = exp( -Eg*t);
     


      cout<<"##### Considering kinematic with......"<<endl;
      cout<<"Eg: "<<Eg<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;



      

      //define vectors to combine FF and BB
      Vfloat th_FF(Corr.Nt);
      Vfloat th_BB(Corr.Nt);
      for(int t=0;t<Corr.Nt;t++) {
	if(t<Corr.Nt/2) {th_FF[t]=1.0; th_BB[t] = 0.0;}
	else {th_BB[t] = 1.0; th_FF[t] = 0.0;}
      }

      //loop over mu and nu
      for(int mu=0;mu<size_mu_nu;mu++) {

	
	for(int nu=0;nu<size_mu_nu;nu++) {

	  cout<<"Analyzing FA and FV for  mu: "<<mu+off_i<<" nu: "<<nu+off_i<<endl;
	  int Im_Re;
	  double parity;

	  //vector
	
	  Corr.Reflection_sign = -1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_F_u = Qu*Corr.corr_t(C_V_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_u = Qu*Corr.corr_t(C_V_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_F_d = Qd*Corr.corr_t(C_V_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_d = Qd*Corr.corr_t(C_V_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
	  distr_t_list vec_u = th_FF*vec_F_u + th_BB*vec_B_u;
	  distr_t_list vec_d = th_FF*vec_F_d + th_BB*vec_B_d;
	  distr_t_list vec = vec_u + vec_d;
	  distr_t_list vec_symm= vec;
	  distr_t_list vec_u_symm= vec_u;
	  distr_t_list vec_d_symm= vec_d;
	  //symmetrize vec
	  for(int t=0; t<Corr.Nt;t++) {
	    vec_symm.distr_list[t] = 0.5*(vec.distr_list[t] + Corr.Reflection_sign*vec.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_u_symm.distr_list[t] = 0.5*(vec_u.distr_list[t] + Corr.Reflection_sign*vec_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_d_symm.distr_list[t] = 0.5*(vec_d.distr_list[t] + Corr.Reflection_sign*vec_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;
	
	

	  //axial
	  if( (mu==0 || nu==0) && (mu != 0 || nu != 0)) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	  else { Im_Re=0; Corr.Reflection_sign=1; parity=1.0;}

	  Corr.Perform_Nt_t_average=0;
	  distr_t_list ax_F_u = Qu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_u = Qu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_F_d = Qd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_d = Qd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	  //For Ds meson
	  distr_t_list ax_Ds_F_u = qqu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_B_u = qqu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_F_d = qqd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_Ds_B_d = qqd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_u = parity*(ax_F_u*th_FF  + ax_B_u*th_BB);
	  distr_t_list ax_d = parity*(ax_F_d*th_FF  + ax_B_d*th_BB);
	  distr_t_list ax = ax_u -ax_d;
	  distr_t_list ax_symm=ax;
	  distr_t_list ax_u_symm= ax_u;
	  distr_t_list ax_d_symm= ax_d;


	  //for Ds meson
	  distr_t_list ax_Ds_u = parity*(ax_Ds_F_u*th_FF  + ax_Ds_B_u*th_BB);
	  distr_t_list ax_Ds_d = parity*(ax_Ds_F_d*th_FF  + ax_Ds_B_d*th_BB);
	  distr_t_list ax_Ds = ax_Ds_u -ax_Ds_d;
	  distr_t_list ax_Ds_symm=ax_Ds;
	  distr_t_list ax_Ds_u_symm= ax_Ds_u;
	  distr_t_list ax_Ds_d_symm= ax_Ds_d;
	  
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_symm.distr_list[t] = 0.5*(ax.distr_list[t] + Corr.Reflection_sign*ax.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_u_symm.distr_list[t] = 0.5*(ax_u.distr_list[t] + Corr.Reflection_sign*ax_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_d_symm.distr_list[t] = 0.5*(ax_d.distr_list[t] + Corr.Reflection_sign*ax_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    //for Ds meson
	    ax_Ds_symm.distr_list[t] = 0.5*(ax_Ds.distr_list[t] + Corr.Reflection_sign*ax_Ds.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_Ds_u_symm.distr_list[t] = 0.5*(ax_Ds_u.distr_list[t] + Corr.Reflection_sign*ax_Ds_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_Ds_d_symm.distr_list[t] = 0.5*(ax_Ds_d.distr_list[t] + Corr.Reflection_sign*ax_Ds_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;

	  //restore standard reflection sign
	  Corr.Reflection_sign=1;


	

	  //get H tensor
	  distr_t_list HA_u= ax_u*EXP_MES*EXP_PH;
	  distr_t_list HA_d= ax_d*EXP_MES*EXP_PH;
	  distr_t_list HV_u= vec_u*EXP_MES*EXP_PH;
	  distr_t_list HV_d= vec_d*EXP_MES*EXP_PH;
	  distr_t_list HA_u_symm= ax_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_d_symm= ax_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_u_symm= vec_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_d_symm= vec_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_symm= ax_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_symm= vec_symm*EXP_MES*EXP_PH;



	  //########################   C TENSOR ###########################


	  //print forward and backward contribution to C
	  Print_To_File({}, {ax_F_u.ave(), ax_F_u.err(), ax_F_d.ave(), ax_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_F_u    ax_F_d");
	  Print_To_File({}, {ax_B_u.ave(), ax_B_u.err(), ax_B_d.ave(), ax_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_B_u    ax_B_d");
	  Print_To_File({}, {vec_F_u.ave(), vec_F_u.err(), vec_F_d.ave(), vec_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_F_u    vec_F_d");
	  Print_To_File({}, {vec_B_u.ave(), vec_B_u.err(), vec_B_d.ave(), vec_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_B_u    vec_B_d");
	

	  //single contributions to C
	  Print_To_File({}, {ax_u.ave() ,ax_u.err(), ax_d.ave(), ax_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u.ave() ,vec_u.err(), vec_d.ave(), vec_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //single contributions to C with symmetrization
	  Print_To_File({}, {ax_u_symm.ave() ,ax_u_symm.err(), ax_d_symm.ave(), ax_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u_symm.ave() ,vec_u_symm.err(), vec_d_symm.ave(), vec_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //error scaling of C^munu(t)
	  for(int t=1;t< Corr.Nt/2;t++) {

	    if(ixg==0) {
	      ofstream Print_scaling_A("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat");
	      Print_scaling_A<<xg.ave()<<" "<<(ax_u_symm*EXP_PH).err(t)<<" "<<(ax_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_A.close();

	      ofstream Print_scaling_V("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat");
	      Print_scaling_V<<xg.ave()<<" "<<(vec_u_symm*EXP_PH).err(t)<<" "<<(vec_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_V.close();
	    }
	    else {
	      ofstream Print_scaling_A("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat", ofstream::app);
	      Print_scaling_A<<xg.ave()<<" "<<(ax_u_symm*EXP_PH).err(t)<<" "<<(ax_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_A.close();

	      ofstream Print_scaling_V("../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/C_scaling_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_t_"+to_string(t)+".dat", ofstream::app);
	      Print_scaling_V<<xg.ave()<<" "<<(vec_u_symm*EXP_PH).err(t)<<" "<<(vec_d_symm*EXP_PH).err(t)<<endl;
	      Print_scaling_V.close();
	    }

	  }

	  //total contribution without symmetrization to C

	  Print_To_File({}, {ax.ave(), ax.err(), vec.ave(), vec.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");

	  //total contribution to C
	  Print_To_File({}, {ax_symm.ave(), ax_symm.err(), vec_symm.ave(), vec_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");


	  //###############################################################


	  //########################   H TENSOR ###########################


	  //single contributions to H without symmetrization
	   Print_To_File({}, {HA_u.ave() ,HA_u.err(), HA_d.ave(), HA_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u.ave() ,HV_u.err(), HV_d.ave(), HV_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");
	  

	  
	  //single contributions to H with symmetrization
	  Print_To_File({}, {HA_u_symm.ave() ,HA_u_symm.err(), HA_d_symm.ave(), HA_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u_symm.ave() ,HV_u_symm.err(), HV_d_symm.ave(), HV_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");

	  //total contribution to H
	  Print_To_File({}, {HA_symm.ave(), HA_symm.err(), HV_symm.ave(), HV_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    HA     HV");

	  //################################################################

	  //push back
	  Ax_tens[mu].push_back(ax_symm);
	  Ax_tens_u[mu].push_back(ax_u_symm);
	  Ax_tens_d[mu].push_back(ax_d_symm);
	  Vec_tens[mu].push_back(vec_symm);
	  Vec_tens_u[mu].push_back(vec_u_symm);
	  Vec_tens_d[mu].push_back(vec_d_symm);
	  //Ds meson
	  Ax_Ds_tens[mu].push_back(ax_Ds_symm);
	  Ax_Ds_tens_u[mu].push_back(ax_Ds_u_symm);
	  Ax_Ds_tens_d[mu].push_back(ax_Ds_d_symm);
	

	}


	for(int nu=0; nu<size_mu_nu_T; nu++) {


	  cout<<"Analyzing FA_T and FV_T for  mu: "<<mu+off_i<<" nu: "<<nu+off_T<<endl;
	  int Im_Re;
	  double parity;

	  //tensor- electric part
	  Corr.Reflection_sign= -1;
	  if( (mu+off_i) != (nu+off_T)) Corr.Reflection_sign = 1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_T_F_u = Qu*Corr.corr_t(C_TK_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_B_u = Qu*Corr.corr_t(C_TK_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_F_d = Qd*Corr.corr_t(C_TK_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_T_B_d = Qd*Corr.corr_t(C_TK_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
	  distr_t_list vec_T_u = th_FF*vec_T_F_u + th_BB*vec_T_B_u;
	  distr_t_list vec_T_d = th_FF*vec_T_F_d + th_BB*vec_T_B_d;
	  distr_t_list vec_T = vec_T_u + vec_T_d;
	  distr_t_list vec_T_symm= vec_T;
	  distr_t_list vec_T_u_symm= vec_T_u;
	  distr_t_list vec_T_d_symm= vec_T_d;

	  //symmetrize vec
	  for(int t=0; t<Corr.Nt;t++) {
	    vec_T_symm.distr_list[t] = 0.5*(vec_T.distr_list[t] + Corr.Reflection_sign*vec_T.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_T_u_symm.distr_list[t] = 0.5*(vec_T_u.distr_list[t] + Corr.Reflection_sign*vec_T_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_T_d_symm.distr_list[t] = 0.5*(vec_T_d.distr_list[t] + Corr.Reflection_sign*vec_T_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;


	  //tensor - magnetic part
	  if(  (  ((mu+off_i) ==0) || ((nu+off_T)==0)) && ( ((mu+off_i) != 0) || ( (nu+off_T) != 0))) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	  else {
	    Im_Re=0;
	    Corr.Reflection_sign=1;
	    if(  (mu+off_i) == (nu + off_T)) Corr.Reflection_sign=-1; 
	    parity=1.0;
	  }


	  Corr.Perform_Nt_t_average=0;
	  distr_t_list ax_T_F_u = Qu*Corr.corr_t(C_BK_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_B_u = Qu*Corr.corr_t(C_BK_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_F_d = Qd*Corr.corr_t(C_BK_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_T_B_d = Qd*Corr.corr_t(C_BK_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_T_u = parity*(ax_T_F_u*th_FF  + ax_T_B_u*th_BB);
	  distr_t_list ax_T_d = parity*(ax_T_F_d*th_FF  + ax_T_B_d*th_BB);
	  distr_t_list ax_T = ax_T_u +ax_T_d;
	  distr_t_list ax_T_symm=ax_T;

	  distr_t_list ax_T_u_symm= ax_T_u;
	  distr_t_list ax_T_d_symm= ax_T_d;
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_T_symm.distr_list[t] = 0.5*(ax_T.distr_list[t] + Corr.Reflection_sign*ax_T.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_T_u_symm.distr_list[t] = 0.5*(ax_T_u.distr_list[t] + Corr.Reflection_sign*ax_T_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_T_d_symm.distr_list[t] = 0.5*(ax_T_d.distr_list[t] + Corr.Reflection_sign*ax_T_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;

	  //restore standard reflection sign
	  Corr.Reflection_sign=1;
	  
	  //get H tensor
	  distr_t_list H_T_A_u= ax_T_u*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_d= ax_T_d*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_u= vec_T_u*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_d= vec_T_d*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_u_symm= ax_T_u_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_d_symm= ax_T_d_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_u_symm= vec_T_u_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_d_symm= vec_T_d_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_A_symm= ax_T_symm*EXP_MES*EXP_PH;
	  distr_t_list H_T_V_symm= vec_T_symm*EXP_MES*EXP_PH;



	  //print forward and backward contribution to C
	  Print_To_File({}, {ax_T_F_u.ave(), ax_T_F_u.err(), ax_T_F_d.ave(), ax_T_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_F_u    ax_F_d");
	  Print_To_File({}, {ax_T_B_u.ave(), ax_T_B_u.err(), ax_T_B_d.ave(), ax_T_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_B_u    ax_B_d");
	  Print_To_File({}, {vec_T_F_u.ave(), vec_T_F_u.err(), vec_T_F_d.ave(), vec_T_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_F_u    vec_F_d");
	  Print_To_File({}, {vec_T_B_u.ave(), vec_T_B_u.err(), vec_T_B_d.ave(), vec_T_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_B_u    vec_B_d");
	

	  //single contributions to C
	  Print_To_File({}, {ax_T_u.ave() ,ax_T_u.err(), ax_T_d.ave(), ax_T_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_T_u.ave() ,vec_T_u.err(), vec_T_d.ave(), vec_T_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //single contributions to C with symmetrization
	  Print_To_File({}, {ax_T_u_symm.ave() ,ax_T_u_symm.err(), ax_T_d_symm.ave(), ax_T_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_T_u_symm.ave() ,vec_T_u_symm.err(), vec_T_d_symm.ave(), vec_T_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");


	  //total contribution without symmetrization to C

	  Print_To_File({}, {ax_T.ave(), ax_T.err(), vec_T.ave(), vec_T.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C_T/"+data_2pts.Tag[iens]+"/wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");

	  //total contribution to C
	  Print_To_File({}, {ax_T_symm.ave(), ax_T_symm.err(), vec_T_symm.ave(), vec_T_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");


	  //###############################################################


	  //########################   H TENSOR ###########################


	  //single contributions to H without symmetrization
	   Print_To_File({}, {H_T_A_u.ave() ,H_T_A_u.err(), H_T_A_d.ave(), H_T_A_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {H_T_V_u.ave() ,H_T_V_u.err(), H_T_V_d.ave(), H_T_V_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");
	  

	  
	  //single contributions to H with symmetrization
	  Print_To_File({}, {H_T_A_u_symm.ave() ,H_T_A_u_symm.err(), H_T_A_d_symm.ave(), H_T_A_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {H_T_V_u_symm.ave() ,H_T_V_u_symm.err(), H_T_V_d_symm.ave(), H_T_V_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");

	  //total contribution to H
	  Print_To_File({}, {H_T_A_symm.ave(), H_T_A_symm.err(), H_T_V_symm.ave(), H_T_V_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H_T/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_T)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    HA     HV");


	  
	  //push back
	  Ax_T_tens[mu].push_back(ax_T_symm);
	  Ax_T_tens_u[mu].push_back(ax_T_u_symm);
	  Ax_T_tens_d[mu].push_back(ax_T_d_symm);
	  Vec_T_tens[mu].push_back(vec_T_symm);
	  Vec_T_tens_u[mu].push_back(vec_T_u_symm);
	  Vec_T_tens_d[mu].push_back(vec_T_d_symm);
	  

	}
      }


      //######################################################################################################################
      //####################     STANDARD     FV     AND    FA     COMPUTATION     ###########################

      //######

      //push_back Ax_tens and Vec_tens
      Ax_glb.push_back(Ax_tens);
      Ax_u_glb.push_back(Ax_tens_u);
      Ax_d_glb.push_back(Ax_tens_d);
      
      Vec_glb.push_back(Vec_tens);
      Vec_u_glb.push_back(Vec_tens_u);
      Vec_d_glb.push_back(Vec_tens_d);


      //For Ds
      Ax_Ds_glb.push_back(Ax_Ds_tens);
      Ax_Ds_u_glb.push_back(Ax_Ds_tens_u);
      Ax_Ds_d_glb.push_back(Ax_Ds_tens_d);

      distr_t_list FA0_distr= 0.5*(Ax_Ds_glb[0][1-off_i][1-off_i] + Ax_Ds_glb[0][2-off_i][2-off_i]);
      
      //define exponential mass factor for each form factor
      //RIGHT ONE IS sign_kz = -1; //relative sign in tensor FF contributions
      double sign_kz= -1;
      vector<distr_t_list> EXP_MASS_CH;
      auto F_EXP = [](double m, double t) { return exp(m*t) ;};
      vector<string> ff_tags({"Au", "Vu", "Ad", "Vd", "ATu", "VTu", "ATd", "VTd"});
      //######
         
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[1-off_i][1-off_i])*EXP_PH  -Qu*FA0_distr, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH, "")));
      
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[1-off_i][1-off_i])*EXP_PH  -Qd*FA0_distr, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH, "")));

      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_T_tens_u[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*EXP_PH, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t(  ((Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u[2-off_i][2-off_T]*(xg/2)))*EXP_PH, "") ));

      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t( (Ax_T_tens_d[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*EXP_PH, "")));
      EXP_MASS_CH.push_back( distr_t_list::f_of_distr_list(F_EXP, Corr.effective_mass_t(  ((Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d[2-off_i][2-off_T]*(xg/2)))*EXP_PH, "") ));

          
    

      

      //Compute FV and FA
    
      distr_t_list FA_distr = -1.0*RF*(0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH - 0.0*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH)*(1.0/(Eg*mel_SMSM))*EXP_MES; 
      distr_t_list FV_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_new_distr = (0.5*RF*Za*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH - (Vec_glb[0][1-off_i][2-off_i] + Vec_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (up-component)
      distr_t_list FA0_u_distr= 0.5*RF*(Ax_u_glb[0][1-off_i][1-off_i] + Ax_u_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_u_distr = -1.0*RF*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_u_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MES;
      distr_t_list FA_u_var= -1.0*RF*Zv*(0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - Qu*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MASS_CH[0];
      distr_t_list FV_u_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_u_new_distr = (0.5*RF*Za*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_u_var = (0.5*RF*Za*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MASS_CH[1];
      distr_t_list FV_u_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH - (Vec_u_glb[0][1-off_i][2-off_i] + Vec_u_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (d-component)
      distr_t_list FA0_d_distr= 0.5*RF*(Ax_d_glb[0][1-off_i][1-off_i] + Ax_d_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_d_distr = -1.0*RF*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/Eg)*(FP_SM/(-1.0*FA0_distr));
      distr_t_list FA_d_new_distr= -1.0*RF*Zv*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MES;
      distr_t_list FA_d_var= -1.0*RF*Zv*(0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - Qd*FA0_distr)*(1.0/(Eg*mel_SMSM))*EXP_MASS_CH[2];
      distr_t_list FV_d_distr = 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_d_new_distr=   (0.5*RF*Za*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MES;
      distr_t_list FV_d_var=   (0.5*RF*Za*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/(kz*mel_SMSM))*EXP_MASS_CH[3];
      distr_t_list FV_d_sub_distr= 0.5*RF*(Za/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH - (Vec_d_glb[0][1-off_i][2-off_i] + Vec_d_glb[0][2-off_i][1-off_i])*1.0)/kz;

      
      
      //Print FV and FA
      Print_To_File({}, {FA_distr.ave(), FA_distr.err(), FA0_distr.ave(), FA0_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_distr.ave(), FV_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_u_distr.ave(), FA_u_distr.err(), FA_u_new_distr.ave(), FA_u_new_distr.err(),  FA_u_var.ave(), FA_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_distr.ave(), FV_u_distr.err(), FV_u_new_distr.ave(), FV_u_new_distr.err(),  FV_u_var.ave(), FV_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_sub_distr.ave(), FV_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_d_distr.ave(), FA_d_distr.err(), FA_d_new_distr.ave(), FA_d_new_distr.err(),  FA_d_var.ave(), FA_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_d_distr).ave(), (FV_d_distr).err(), FV_d_new_distr.ave(), FV_d_new_distr.err(),  FV_d_var.ave(), FV_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_d_sub_distr).ave(), (FV_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors


      
      //set interval depending on xg
      int Tmin_V, Tmax_V, Tmin_A, Tmax_A;
      if(xg.ave() > 1e-10) {
      //set time intervals for each heavy mass
      //FA_u
      Get_Bs_Mh_Tmin_Tmax("FAu",MESON, Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_A; Corr.Tmax= Tmax_A;
      distr_t FA_u_fit= Corr.Fit_distr(FA_u_new_distr);
      FA_u.distr_list.push_back(FA_u_fit);
      FA_u_per_kin[ixg-1].distr_list.push_back(FA_u_fit);
      ax_u_Tmin.push_back( Corr.Tmin); ax_u_Tmax.push_back( Corr.Tmax);
      //FA_d
      Get_Bs_Mh_Tmin_Tmax("FAd",MESON, Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_A; Corr.Tmax= Tmax_A;
      distr_t FA_d_fit= Corr.Fit_distr(FA_d_new_distr);
      FA_d.distr_list.push_back(FA_d_fit);
      FA_d_per_kin[ixg-1].distr_list.push_back(FA_d_fit);
      ax_d_Tmin.push_back( Corr.Tmin); ax_d_Tmax.push_back( Corr.Tmax);
      ax_Tmin.push_back( Corr.Tmin); ax_Tmax.push_back( Corr.Tmax);
      FA.distr_list.push_back( FA_u_fit+ FA_d_fit);
      FA_per_kin[ixg-1].distr_list.push_back( FA_u_fit + FA_d_fit);

      //FV_u
      Get_Bs_Mh_Tmin_Tmax("FVu",MESON, Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_V; Corr.Tmax= Tmax_V;
      distr_t FV_u_fit= Corr.Fit_distr(FV_u_new_distr);
      FV_u.distr_list.push_back(FV_u_fit);
      FV_u_per_kin[ixg-1].distr_list.push_back(FV_u_fit);
      vec_u_Tmin.push_back( Corr.Tmin); vec_u_Tmax.push_back( Corr.Tmax);
      //FV_d
      Get_Bs_Mh_Tmin_Tmax("FVd",MESON, Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_V; Corr.Tmax= Tmax_V;
      distr_t FV_d_fit= Corr.Fit_distr(FV_d_new_distr);
      FV_d.distr_list.push_back(FV_d_fit);
      FV_d_per_kin[ixg-1].distr_list.push_back(FV_d_fit);
      vec_d_Tmin.push_back( Corr.Tmin); vec_d_Tmax.push_back( Corr.Tmax);
      vec_Tmin.push_back( Corr.Tmin); vec_Tmax.push_back( Corr.Tmax);
      FV.distr_list.push_back( FV_u_fit+ FV_d_fit);
      FV_per_kin[ixg-1].distr_list.push_back( FV_u_fit + FV_d_fit);
      }





      //######################################################################################################################



      //######################################################################################################################
      //####################     FV_T     AND    FA_T     COMPUTATION     ###########################

      //push_back Ax_tens and Vec_tens
      Ax_T_glb.push_back(Ax_T_tens);
      Ax_T_u_glb.push_back(Ax_T_tens_u);
      Ax_T_d_glb.push_back(Ax_T_tens_d);
      
      Vec_T_glb.push_back(Vec_T_tens);
      Vec_T_u_glb.push_back(Vec_T_tens_u);
      Vec_T_d_glb.push_back(Vec_T_tens_d);
    
      //Compute FV and FA
      distr_t_list FA_T_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens[1-off_i][2-off_T]*(xg/2) + Ax_T_tens[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens[2-off_i][2-off_T]*(xg/2)))*(1.0/Eg)*EXP_PH;
      //improved estimator
      distr_t_list FA_T_new_distr =(0.5*RF*(Z_T)*( Ax_T_tens[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens[1-off_i][2-off_T]*(xg/2) + Ax_T_tens[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      
      //compute FV and FA (up-component)
      distr_t_list FA_T_u_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens_u[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_u_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u[2-off_i][2-off_T]*(xg/2)))*(1.0/Eg)*EXP_PH;
      //improved estimator
      distr_t_list FA_T_u_new_distr =(0.5*RF*(Z_T)*( Ax_T_tens_u[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_u_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      //another variation
      distr_t_list FA_T_u_var =(0.5*RF*(Z_T)*( Ax_T_tens_u[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_u[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_u[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_u[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[4];
      distr_t_list FV_T_u_var = (0.5*RF*(Z_T)*( (Vec_T_tens_u[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_u[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_u[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_u[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[5];
      

      //compute FV and FA (d-component)
      distr_t_list FA_T_d_distr =0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( Ax_T_tens_d[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/Eg)*EXP_PH;
      distr_t_list FV_T_d_distr = 0.5*RF*(Z_T/Zv)*(FP_SM/(-1.0*FA0_distr))*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d[2-off_i][2-off_T]*(xg/2)))*(1.0/Eg)*EXP_PH;
      //improved estimator
      distr_t_list FA_T_d_new_distr =(0.5*RF*(Z_T)*( Ax_T_tens_d[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      distr_t_list FV_T_d_new_distr = (0.5*RF*(Z_T)*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MES;
      //another variation
      distr_t_list FA_T_d_var =(0.5*RF*(Z_T)*( Ax_T_tens_d[1-off_i][1-off_T]*(1-xg/2) - sign_kz*Vec_T_tens_d[1-off_i][2-off_T]*(xg/2) + Ax_T_tens_d[2-off_i][2-off_T]*(1-xg/2) +sign_kz*Vec_T_tens_d[2-off_i][1-off_T]*(xg/2))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[6];
      distr_t_list FV_T_d_var = (0.5*RF*(Z_T)*( (Vec_T_tens_d[1-off_i][2-off_T]*(1-xg/2) -sign_kz*Ax_T_tens_d[1-off_i][1-off_T]*(xg/2)) -(Vec_T_tens_d[2-off_i][1-off_i]*(1-xg/2) + sign_kz*Ax_T_tens_d[2-off_i][2-off_T]*(xg/2)))*(1.0/(mel_SMSM*Eg))*EXP_PH)*EXP_MASS_CH[7];
   
      //Print FV and FA
      Print_To_File({}, {FA_T_distr.ave(), FA_T_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_T_distr.ave(), FV_T_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      //Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_T_u_distr.ave(), FA_T_u_distr.err(), FA_T_u_new_distr.ave(), FA_T_u_new_distr.err(), FA_T_u_var.ave(), FA_T_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_T_u_distr.ave(), FV_T_u_distr.err(), FV_T_u_new_distr.ave(), FV_T_u_new_distr.err(), FV_T_u_var.ave(), FV_T_u_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      //Print_To_File({}, {FV_T_u_sub_distr.ave(), FV_T_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_T_d_distr.ave(), FA_T_d_distr.err(), FA_T_d_new_distr.ave(), FA_T_d_new_distr.err(), FA_T_d_var.ave(), FA_T_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(FV_T_d_distr).ave(), (FV_T_d_distr).err(), FV_T_d_new_distr.ave(), FV_T_d_new_distr.err(), FV_T_d_var.ave(), FV_T_d_var.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_T_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      // Print_To_File({}, {(FV_T_d_sub_distr).ave(), (FV_T_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_T_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors

      //set interval depending on xg
      int Tmin_T_V, Tmax_T_V, Tmin_T_A, Tmax_T_A;
      
      if(xg.ave() > 1e-10) {
      //set time intervals for each heavy mass
      //FA_T_u
      Get_Bs_Mh_Tmin_Tmax("FA_Tu",MESON, Tmin_T_A, Tmax_T_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_A; Corr.Tmax= Tmax_T_A;
      distr_t FA_T_u_fit= Corr.Fit_distr(FA_T_u_new_distr);
      FA_T_u.distr_list.push_back(FA_T_u_fit);
      FA_T_u_per_kin[ixg-1].distr_list.push_back(FA_T_u_fit);
      ax_T_u_Tmin.push_back( Corr.Tmin); ax_T_u_Tmax.push_back( Corr.Tmax);
      //FA_T_d
      Get_Bs_Mh_Tmin_Tmax("FA_Td", MESON, Tmin_T_A, Tmax_T_A, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_A; Corr.Tmax= Tmax_T_A;
      distr_t FA_T_d_fit= Corr.Fit_distr(FA_T_d_new_distr);
      FA_T_d.distr_list.push_back(FA_T_d_fit);
      FA_T_d_per_kin[ixg-1].distr_list.push_back(FA_T_d_fit);
      ax_T_d_Tmin.push_back( Corr.Tmin); ax_T_d_Tmax.push_back( Corr.Tmax);
      ax_T_Tmin.push_back( Corr.Tmin); ax_T_Tmax.push_back( Corr.Tmax);
      FA_T.distr_list.push_back( FA_T_u_fit+ FA_T_d_fit);
      FA_T_per_kin[ixg-1].distr_list.push_back( FA_T_u_fit + FA_T_d_fit);

      //FV_T_u
      Get_Bs_Mh_Tmin_Tmax("FV_Tu", MESON, Tmin_T_V, Tmax_T_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_V; Corr.Tmax= Tmax_T_V;
      distr_t FV_T_u_fit= Corr.Fit_distr(FV_T_u_new_distr);
      FV_T_u.distr_list.push_back(FV_T_u_fit);
      FV_T_u_per_kin[ixg-1].distr_list.push_back(FV_T_u_fit);
      vec_T_u_Tmin.push_back( Corr.Tmin); vec_T_u_Tmax.push_back( Corr.Tmax);
      //FV_T_d
      Get_Bs_Mh_Tmin_Tmax("FV_Td", MESON, Tmin_T_V, Tmax_T_V, ixg, data_2pts.Tag[iens]);
      Corr.Tmin = Tmin_T_V; Corr.Tmax= Tmax_T_V;
      distr_t FV_T_d_fit= Corr.Fit_distr(FV_T_d_new_distr);
      FV_T_d.distr_list.push_back(FV_T_d_fit);
      FV_T_d_per_kin[ixg-1].distr_list.push_back(FV_T_d_fit);
      vec_T_d_Tmin.push_back( Corr.Tmin); vec_T_d_Tmax.push_back( Corr.Tmax);
      vec_T_Tmin.push_back( Corr.Tmin); vec_T_Tmax.push_back( Corr.Tmax);
      FV_T.distr_list.push_back( FV_T_u_fit+ FV_T_d_fit);
      FV_T_per_kin[ixg-1].distr_list.push_back( FV_T_u_fit + FV_T_d_fit);
      }

      




      //######################################################################################################################
      
      
    
      
    }

    //for FV and FA
    FA_per_ens[iens] = FA;
    FV_per_ens[iens] = FV;
    FA_u_per_ens[iens] = FA_u;
    FA_d_per_ens[iens] = FA_d;
    FV_u_per_ens[iens] = FV_u;
    FV_d_per_ens[iens] = FV_d;
    //for FV_T and FA_T
    FA_T_per_ens[iens] = FA_T;
    FV_T_per_ens[iens] = FV_T;
    FA_T_u_per_ens[iens] = FA_T_u;
    FA_T_d_per_ens[iens] = FA_T_d;
    FV_T_u_per_ens[iens] = FV_T_u;
    FV_T_d_per_ens[iens] = FV_T_d;
    


    //Print fitted form factors
    //FV and FA
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV.ave(), FV.err(), vec_Tmin, vec_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA.ave(), FA.err(), ax_Tmin, ax_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_u.ave(), FV_u.err(), vec_u_Tmin, vec_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_u.ave(), FA_u.err(), ax_u_Tmin, ax_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (FV_d).ave(), (FV_d).err(), vec_d_Tmin, vec_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_d.ave(), FA_d.err(), ax_d_Tmin, ax_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");


    //FV_T and FA_T
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_T.ave(), FV_T.err(), vec_T_Tmin, vec_T_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T.ave(), FA_T.err(), ax_T_Tmin, ax_T_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_T_u.ave(), FV_T_u.err(), vec_T_u_Tmin, vec_T_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T_u.ave(), FA_T_u.err(), ax_T_u_Tmin, ax_T_u_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (FV_T_d).ave(), (FV_T_d).err(), vec_T_d_Tmin, vec_T_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV_T.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_T_d.ave(), FA_T_d.err(), ax_T_d_Tmin, ax_T_d_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA_T.dat", "", "#xg Dxg  FA DFA");


  }


  


  if(num_xg==1) crash("Exiting...No interpolation nor continuum extrapolation to perform, only xg available is zero");


  //Print all ensembles for given xg
  for(int ixg=1;ixg<num_xg;ixg++) {

    //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_per_kin[ixg-1].ave(), FA_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_per_kin[ixg-1].ave(), FV_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_per_kin[ixg-1].ave(), FA_T_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_per_kin[ixg-1].ave(), FV_T_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");


    //########### u contribution ############
    //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_u_per_kin[ixg-1].ave(), FA_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_u_per_kin[ixg-1].ave(), FV_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_u_per_kin[ixg-1].ave(), FA_T_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_u_per_kin[ixg-1].ave(), FV_T_u_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_u/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");


    //#######################################



    //########### d contribution ############
     //for FA and FV
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_d_per_kin[ixg-1].ave(), FA_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_d_per_kin[ixg-1].ave(), FV_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //for FA_T and FV_T
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_T_d_per_kin[ixg-1].ave(), FA_T_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FA_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_T_d_per_kin[ixg-1].ave(), FV_T_d_per_kin[ixg-1].err(), L_list}, "../data/ph_emission/"+ph_type_mes+"/FF_d/per_kin/FV_T_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err L/a");

    //#######################################



    
    
  }

  //################################################################################################

  //interpolate form factors for each ensemble

  if(num_xg > 5) {

  cout<<"Interpolating vector and tensor form factors!"<<endl;

  //for FV and FA
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FA_interpol_jacks(Nens);
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FV_interpol_jacks(Nens);

  //for FV_T and FA_T
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FA_T_interpol_jacks(Nens);
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FV_T_interpol_jacks(Nens);
  

  for(int iens=0;iens<Nens;iens++) {

    int num_xg_iens= FA_per_ens[iens].size();

    for(int ijack=0;ijack<Njacks;ijack++) {

      //For FV and FA
      Vfloat FA_jacks, FV_jacks;
      for(int ixg=1;ixg<num_xg_iens;ixg++) { FA_jacks.push_back( FA_per_ens[iens].distr_list[ixg-1].distr[ijack]); FV_jacks.push_back( FV_per_ens[iens].distr_list[ixg-1].distr[ijack]);}
      FA_interpol_jacks[iens].emplace_back( FA_jacks.begin(), FA_jacks.end(), 0.1, 0.1);
      FV_interpol_jacks[iens].emplace_back( FV_jacks.begin(), FV_jacks.end(), 0.1, 0.1);
      //For FV_T and FA_T
      Vfloat FA_T_jacks, FV_T_jacks;
      for(int ixg=1;ixg<num_xg_iens;ixg++) { FA_T_jacks.push_back( FA_T_per_ens[iens].distr_list[ixg-1].distr[ijack]); FV_T_jacks.push_back( FV_T_per_ens[iens].distr_list[ixg-1].distr[ijack]);}
      FA_T_interpol_jacks[iens].emplace_back( FA_T_jacks.begin(), FA_T_jacks.end(), 0.1, 0.1);
      FV_T_interpol_jacks[iens].emplace_back( FV_T_jacks.begin(), FV_T_jacks.end(), 0.1, 0.1);
    }
  }

  auto FA_interpol_distr = [&FA_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {
			     distr_t return_distr(UseJack);
			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_interpol_jacks[iens][ijack](xg));}
			     return return_distr;
  };
  auto FV_interpol_distr = [&FV_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {
			     distr_t return_distr(UseJack);
			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_interpol_jacks[iens][ijack](xg));}
			     return return_distr;
  } ;

  auto FA_T_interpol_distr = [&FA_T_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {
			     distr_t return_distr(UseJack);
			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_T_interpol_jacks[iens][ijack](xg));}
			     return return_distr;
  };
  auto FV_T_interpol_distr = [&FV_T_interpol_jacks, &UseJack, &Njacks](double xg, int iens) -> distr_t {
			     distr_t return_distr(UseJack);
			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_T_interpol_jacks[iens][ijack](xg));}
			     return return_distr;
  } ;



  for(int iens=0;iens<Nens;iens++) {

    //for FV and FA
    distr_t_list FA_interpol_to_print_distr(UseJack);
    distr_t_list FV_interpol_to_print_distr(UseJack);
    for(auto &X: Bs_xg_to_spline) {
      FA_interpol_to_print_distr.distr_list.push_back( FA_interpol_distr(X, iens));
      FV_interpol_to_print_distr.distr_list.push_back( FV_interpol_distr(X, iens));
    }

    Print_To_File({}, {Bs_xg_to_spline, FA_interpol_to_print_distr.ave(), FA_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/"+data_2pts.Tag[iens]+"/FA_interpol.dat", "", "#xg FA FA_err");
    Print_To_File({}, {Bs_xg_to_spline, FV_interpol_to_print_distr.ave(), FV_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/"+data_2pts.Tag[iens]+"/FV_interpol.dat", "", "#xg FV FV_err");


    //for FV_T and FA_T
    distr_t_list FA_T_interpol_to_print_distr(UseJack);
    distr_t_list FV_T_interpol_to_print_distr(UseJack);
    for(auto &X: Bs_xg_to_spline) {
      FA_T_interpol_to_print_distr.distr_list.push_back( FA_T_interpol_distr(X, iens));
      FV_T_interpol_to_print_distr.distr_list.push_back( FV_T_interpol_distr(X, iens));
    }

    Print_To_File({}, {Bs_xg_to_spline, FA_T_interpol_to_print_distr.ave(), FA_T_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/"+data_2pts.Tag[iens]+"/FA_T_interpol.dat", "", "#xg FA FA_err");
    Print_To_File({}, {Bs_xg_to_spline, FV_T_interpol_to_print_distr.ave(), FV_T_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/"+data_2pts.Tag[iens]+"/FV_T_interpol.dat", "", "#xg FV FV_err");


  }

  }


  //Print MP, FP, MP_ov_FP
  Print_To_File({}, {a_distr_list_red.ave(), MP_list_red.ave(), MP_list_red.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/mass/masses.list", "", "#a MP MP_err");
  Print_To_File({}, {a_distr_list.ave(), FP_list.ave(), FP_list.err(), MP_ov_FP_list.ave(), MP_ov_FP_list.err(), ZT_list.ave(), ZT_list.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/fP.list", "", "#a FP FP_err MP/FP MP/FP_err ZT ");


  //#################################################################################################

  cout<<"Estimate finite-volume effects"<<endl;


  //If ensemble B96 is present estimate FSEs
  int iens_B64=-1;
  int iens_B96=-1;
  for(int iens=0;iens<Nens;iens++) {
    if(data_2pts.Tag[iens] == "cB211b.072.64") iens_B64=iens;
    else if(data_2pts.Tag[iens] == "cB211b.072.96") iens_B96=iens;
  }

  if( (iens_B64 > -1)  && (iens_B96 > -1)) { //compute FSE

    //############################################################################################################################
    //###########################################################################################################################
    //###########################################################################################################################
    //FOR FV and FA
    ofstream Print_FVE("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_per_ens[iens_B96].err(ixg-1),2) +pow(FA_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_per_ens[iens_B96].ave(ixg-1)-FA_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_per_ens[iens_B96].err(ixg-1),2) +pow(FV_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_per_ens[iens_B96].ave(ixg-1)-FV_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;

      if(ixg==7) {
	 Print_FVE<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_per_ens[iens_B64].ave(ixg-1))/FA_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_per_ens[iens_B64].ave(ixg-1))/FV_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
    Print_FVE.close();



    //u contribution
    ofstream Print_FVE_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_u_per_ens[iens_B96].err(ixg-1),2) +pow(FA_u_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_u_per_ens[iens_B96].ave(ixg-1)-FA_u_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_u_per_ens[iens_B96].err(ixg-1),2) +pow(FV_u_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_u_per_ens[iens_B96].ave(ixg-1)-FV_u_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_u_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_u_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_u<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;

      if(ixg==7) {
	 Print_FVE_u<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_u<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_u<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_u_per_ens[iens_B64].ave(ixg-1))/FA_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_u_per_ens[iens_B64].ave(ixg-1))/FV_u_per_ens[iens_B64].err(ixg-1))<<endl;
      }
      
      
    }
 
    Print_FVE_u.close();



    //d contribution
    ofstream Print_FVE_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_d_per_ens[iens_B96].err(ixg-1),2) + pow(FA_d_per_ens[iens_B64].err(ixg-1),2));
      double FA_corr= fabs(FA_d_per_ens[iens_B96].ave(ixg-1)-FA_d_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_d_per_ens[iens_B96].err(ixg-1),2) + pow(FV_d_per_ens[iens_B64].err(ixg-1),2));
      double FV_corr= fabs(FV_d_per_ens[iens_B96].ave(ixg-1)-FV_d_per_ens[iens_B64].ave(ixg-1));
      FA_corr= fabs((FA_corr/FA_d_per_ens[iens_B64].ave(ixg-1))*erf( FA_corr/(sqrt(2.0)*sA_max)));
      FV_corr= fabs((FV_corr/FV_d_per_ens[iens_B64].ave(ixg-1))*erf( FV_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_d<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;

      if(ixg==7) {
	 Print_FVE_d<<"0.80"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_d<<"0.90"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_d<<"1.00"<<" "<<FA_corr<<" "<<FV_corr<<" "<<(FA_corr*fabs(FA_d_per_ens[iens_B64].ave(ixg-1))/FA_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_corr*fabs(FV_d_per_ens[iens_B64].ave(ixg-1))/FV_d_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
    Print_FVE_d.close();
    //############################################################################################################################
    //###########################################################################################################################
    //###########################################################################################################################








    //############################################################################################################################
    //###########################################################################################################################
    //###########################################################################################################################
    //FOR FV_T and FA_T
    ofstream Print_FVE_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_T_per_ens[iens_B96].err(ixg-1),2) +pow(FA_T_per_ens[iens_B64].err(ixg-1),2));
      double FA_T_corr= fabs(FA_T_per_ens[iens_B96].ave(ixg-1)-FA_T_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_T_per_ens[iens_B96].err(ixg-1),2) +pow(FV_T_per_ens[iens_B64].err(ixg-1),2));
      double FV_T_corr= fabs(FV_T_per_ens[iens_B96].ave(ixg-1)-FV_T_per_ens[iens_B64].ave(ixg-1));
      FA_T_corr= fabs((FA_T_corr/FA_T_per_ens[iens_B64].ave(ixg-1))*erf( FA_T_corr/(sqrt(2.0)*sA_max)));
      FV_T_corr= fabs((FV_T_corr/FV_T_per_ens[iens_B64].ave(ixg-1))*erf( FV_T_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_T<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_per_ens[iens_B64].ave(ixg-1))/FA_T_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_per_ens[iens_B64].ave(ixg-1))/FV_T_per_ens[iens_B64].err(ixg-1))<<endl;

      if(ixg==7) {
	 Print_FVE_T<<"0.80"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_per_ens[iens_B64].ave(ixg-1))/FA_T_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_per_ens[iens_B64].ave(ixg-1))/FV_T_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_T<<"0.90"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_per_ens[iens_B64].ave(ixg-1))/FA_T_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_per_ens[iens_B64].ave(ixg-1))/FV_T_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_T<<"1.00"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_per_ens[iens_B64].ave(ixg-1))/FA_T_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_per_ens[iens_B64].ave(ixg-1))/FV_T_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
    Print_FVE_T.close();



    //u contribution
    ofstream Print_FVE_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_T_u_per_ens[iens_B96].err(ixg-1),2) +pow(FA_T_u_per_ens[iens_B64].err(ixg-1),2));
      double FA_T_corr= fabs(FA_T_u_per_ens[iens_B96].ave(ixg-1)-FA_T_u_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_T_u_per_ens[iens_B96].err(ixg-1),2) +pow(FV_T_u_per_ens[iens_B64].err(ixg-1),2));
      double FV_T_corr= fabs(FV_T_u_per_ens[iens_B96].ave(ixg-1)-FV_T_u_per_ens[iens_B64].ave(ixg-1));
      FA_T_corr= fabs((FA_T_corr/FA_T_u_per_ens[iens_B64].ave(ixg-1))*erf( FA_T_corr/(sqrt(2.0)*sA_max)));
      FV_T_corr= fabs((FV_T_corr/FV_T_u_per_ens[iens_B64].ave(ixg-1))*erf( FV_T_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_T_u<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_u_per_ens[iens_B64].ave(ixg-1))/FA_T_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_u_per_ens[iens_B64].ave(ixg-1))/FV_T_u_per_ens[iens_B64].err(ixg-1))<<endl;


      if(ixg==7) {
	 Print_FVE_T_u<<"0.80"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_u_per_ens[iens_B64].ave(ixg-1))/FA_T_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_u_per_ens[iens_B64].ave(ixg-1))/FV_T_u_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_T_u<<"0.90"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_u_per_ens[iens_B64].ave(ixg-1))/FA_T_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_u_per_ens[iens_B64].ave(ixg-1))/FV_T_u_per_ens[iens_B64].err(ixg-1))<<endl;
	 Print_FVE_T_u<<"1.00"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_u_per_ens[iens_B64].ave(ixg-1))/FA_T_u_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_u_per_ens[iens_B64].ave(ixg-1))/FV_T_u_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
   
    Print_FVE_T_u.close();



    //d contribution
    ofstream Print_FVE_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {

      double sA_max= sqrt( pow(FA_T_d_per_ens[iens_B96].err(ixg-1),2) + pow(FA_T_d_per_ens[iens_B64].err(ixg-1),2));
      double FA_T_corr= fabs(FA_T_d_per_ens[iens_B96].ave(ixg-1)-FA_T_d_per_ens[iens_B64].ave(ixg-1));
      double sV_max= sqrt( pow(FV_T_d_per_ens[iens_B96].err(ixg-1),2) + pow(FV_T_d_per_ens[iens_B64].err(ixg-1),2));
      double FV_T_corr= fabs(FV_T_d_per_ens[iens_B96].ave(ixg-1)-FV_T_d_per_ens[iens_B64].ave(ixg-1));
      FA_T_corr= fabs((FA_T_corr/FA_T_d_per_ens[iens_B64].ave(ixg-1))*erf( FA_T_corr/(sqrt(2.0)*sA_max)));
      FV_T_corr= fabs((FV_T_corr/FV_T_d_per_ens[iens_B64].ave(ixg-1))*erf( FV_T_corr/(sqrt(2.0)*sV_max)));
      Print_FVE_T_d<<to_string_with_precision(ixg*0.1,2)<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_d_per_ens[iens_B64].ave(ixg-1))/FA_T_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_d_per_ens[iens_B64].ave(ixg-1))/FV_T_d_per_ens[iens_B64].err(ixg-1))<<endl;



      if(ixg==7) {
	Print_FVE_T_d<<"0.80"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_d_per_ens[iens_B64].ave(ixg-1))/FA_T_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_d_per_ens[iens_B64].ave(ixg-1))/FV_T_d_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_T_d<<"0.90"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_d_per_ens[iens_B64].ave(ixg-1))/FA_T_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_d_per_ens[iens_B64].ave(ixg-1))/FV_T_d_per_ens[iens_B64].err(ixg-1))<<endl;
	Print_FVE_T_d<<"1.00"<<" "<<FA_T_corr<<" "<<FV_T_corr<<" "<<(FA_T_corr*fabs(FA_T_d_per_ens[iens_B64].ave(ixg-1))/FA_T_d_per_ens[iens_B64].err(ixg-1))<<" "<<(FV_T_corr*fabs(FV_T_d_per_ens[iens_B64].ave(ixg-1))/FV_T_d_per_ens[iens_B64].err(ixg-1))<<endl;
      }

      
    }
    
    Print_FVE_T_d.close();
    //############################################################################################################################
    //###########################################################################################################################
    //###########################################################################################################################


    
  }
   
  //continuum extrapolation


 
 

  if(Perform_continuum_extrapolation) {


   
  
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_cont_interpol_jacks;

    //u contribution
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_u_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_u_cont_interpol_jacks;

    //d contribution
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_d_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_d_cont_interpol_jacks;



    //load FVE
    //VECTOR FORM FACTORS

    ifstream Read_FVE("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A;
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE.eof()) crash("Read_FVE reached eof before expected");
      Read_FVE>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_per_ens[iens].distr_list[ixg-1] = FA_per_ens[iens].distr_list[ixg-1]*FVE_A;
      FV_per_ens[iens].distr_list[ixg-1] = FV_per_ens[iens].distr_list[ixg-1]*FVE_V;
    }
    }
    Read_FVE.close();



    //u contribution
    ifstream Read_FVE_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A;
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE_u.eof()) crash("Read_FVE_u reached eof before expected");
      Read_FVE_u>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_u_per_ens[iens].distr_list[ixg-1] = FA_u_per_ens[iens].distr_list[ixg-1]*FVE_A;
      FV_u_per_ens[iens].distr_list[ixg-1] = FV_u_per_ens[iens].distr_list[ixg-1]*FVE_V;
    }
    }
    Read_FVE_u.close();


    //d contribution
    ifstream Read_FVE_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V(UseJack), FVE_A(UseJack);
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE_d.eof()) crash("Read_FVE_d reached eof before expected");
      Read_FVE_d>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_d_per_ens[iens].distr_list[ixg-1] = FA_d_per_ens[iens].distr_list[ixg-1]*FVE_A;
      FV_d_per_ens[iens].distr_list[ixg-1] = FV_d_per_ens[iens].distr_list[ixg-1]*FVE_V;
    }
    }
    Read_FVE_d.close();



    //TENSOR FORM FACTORS
    ifstream Read_FVE_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V_T(UseJack), FVE_A_T;
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE_T.eof()) crash("Read_FVE_T reached eof before expected");
      Read_FVE_T>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A_T.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V_T.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_T_per_ens[iens].distr_list[ixg-1] = FA_per_ens[iens].distr_list[ixg-1]*FVE_A_T;
      FV_T_per_ens[iens].distr_list[ixg-1] = FV_per_ens[iens].distr_list[ixg-1]*FVE_V_T;
    }
    }
    Read_FVE_T.close();



    //u contribution
    ifstream Read_FVE_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T_u.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V_T(UseJack), FVE_A_T;
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE_T_u.eof()) crash("Read_FVE_T_u reached eof before expected");
      Read_FVE_T_u>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A_T.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V_T.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_T_u_per_ens[iens].distr_list[ixg-1] = FA_T_u_per_ens[iens].distr_list[ixg-1]*FVE_A_T;
      FV_T_u_per_ens[iens].distr_list[ixg-1] = FV_T_u_per_ens[iens].distr_list[ixg-1]*FVE_V_T;
    }
    }
    Read_FVE_T_u.close();


    //d contribution
    ifstream Read_FVE_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/FVE_T_d.dat");
    for(int ixg=1;ixg<num_xg;ixg++) {
      distr_t FVE_V_T(UseJack), FVE_A_T(UseJack);
      double xx, sV, sA, srel_V, srel_A;
      if(Read_FVE_T_d.eof()) crash("Read_FVE_T_d reached eof before expected");
      Read_FVE_T_d>>xx>>sA>>sV>>srel_A>>srel_V;
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_A_T.distr.push_back( 1.0 + sA*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
      for(int ijack=0;ijack<Njacks;ijack++) { FVE_V_T.distr.push_back( 1.0 + sV*GM()/((UseJack==true)?sqrt(Njacks -1.0):1.0));   }
    for(int iens=0;iens<Nens;iens++) {
      FA_T_d_per_ens[iens].distr_list[ixg-1] = FA_T_d_per_ens[iens].distr_list[ixg-1]*FVE_A_T;
      FV_T_d_per_ens[iens].distr_list[ixg-1] = FV_T_d_per_ens[iens].distr_list[ixg-1]*FVE_V_T;
    }
    }
    Read_FVE_T_d.close();
    
    

  class ipar_FF_Nissa {
  public:
    ipar_FF_Nissa() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, a;
    int is;
  };
  
  class fpar_FF_Nissa {
  public:
    fpar_FF_Nissa() {}
    fpar_FF_Nissa(const Vfloat &par) {
      if((signed)par.size() != 3) crash("In class fpar_FF_Nissa  class constructor Vfloat par has size != 3");
      F0=par[0];
      D1=par[1];
      D2=par[2];
    }
    double F0, D1,D2;
  };

  
  //init bootstrap fit
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF(Njacks);
  bf_FF.set_warmup_lev(0); //sets warmup
  bf_FF.Set_number_of_measurements(Nens);
  bf_FF.Set_verbosity(1);
  bf_FF.Add_par("F0", 0.07, 0.001);
  bf_FF.Add_par("D1", 1.0, 0.1);
  bf_FF.Add_par("D2", 1.0, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF_ch2(1);
  bf_FF_ch2.set_warmup_lev(0); //sets warmup
  bf_FF_ch2.Set_number_of_measurements(Nens);
  bf_FF_ch2.Set_verbosity(1);
  bf_FF_ch2.Add_par("F0", 0.07, 0.001);
  bf_FF_ch2.Add_par("D1", 1.0, 0.1);
  bf_FF_ch2.Add_par("D2", 1.0, 0.1);

  //count number of dof
  int dof= Nens-2;
  if(Include_a4==false) { bf_FF.Fix_par("D2", 0.0); bf_FF_ch2.Fix_par("D2",0.0);}
  else dof-- ;

  if( Use_three_finest) dof--;


  //push_back info
  rt_fit.Nmeas=Nens-(Use_three_finest==true);
  rt_fit.Ndof= dof;
  rt_fit.Npars= rt_fit.Nmeas-dof;
  

  //ansatz
  bf_FF.ansatz=  [&Use_three_finest ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

    if( (ip.is==0) && Use_three_finest) return 0.0;
		  
    return p.F0 + p.D1*pow(ip.a*Lambda_QCD,2) + p.D2*pow(ip.a*Lambda_QCD,4);

  };

  
  bf_FF.measurement=  [&Use_three_finest ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {
    
    if( (ip.is==0) && Use_three_finest) return 0.0;
    
    return ip.FF;
       
  };
  bf_FF.error=  [ ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

		 return ip.FF_err;
		 };

  bf_FF_ch2.ansatz= bf_FF.ansatz;
  bf_FF_ch2.measurement = bf_FF.measurement;
  bf_FF_ch2.error = bf_FF.error;



  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                           FITTING FV AND FA                     ##########################################
  //###########################################                     INCLUDING SINGLE CONTRIBUTIONS              ##########################################
  //###########################################                        FROM STRANGE AND "CHARM"                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################




  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA for all xg                             ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  Vfloat ch2_FA;
  distr_t_list F0_A_list(UseJack);
  distr_t_list D1_A_list(UseJack);
  distr_t_list D2_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_A_list.distr_list.push_back(F0);
  D1_A_list.distr_list.push_back(D1);
  D2_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV for all xg                             ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  
  Vfloat ch2_FV;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_V_list(UseJack);
  distr_t_list D1_V_list(UseJack);
  distr_t_list D2_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_V_list.distr_list.push_back(F0);
  D1_V_list.distr_list.push_back(D1);
  D2_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_A_list.ave(), F0_A_list.err(), (D1_A_list/F0_A_list).ave(), (D1_A_list/F0_A_list).err(), (D2_A_list/F0_A_list).ave(), (D2_A_list/F0_A_list).err(), ch2_FA}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_V_list.ave(), F0_V_list.err(), (D1_V_list/F0_V_list).ave(), (D1_V_list/F0_V_list).err(), (D2_V_list/F0_V_list).ave(), (D2_V_list/F0_V_list).err(), ch2_FV}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


  //Print covariance matrix
  string Analysis_tag= ( (UseJack==true)?"jack":"boot")+_Fit_tag;
  Eigen::MatrixXd Cov_FA(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA(x_xg-1, y_xg-1) = F0_A_list.distr_list[x_xg-1]%F0_A_list.distr_list[y_xg-1];
      Cov_FV(x_xg-1, y_xg-1) = F0_V_list.distr_list[x_xg-1]%F0_V_list.distr_list[y_xg-1];
      Corr_FA(x_xg-1, y_xg-1) = Cov_FA(x_xg-1,y_xg-1)/(F0_A_list.err(x_xg-1)*F0_A_list.err(y_xg-1));
      Corr_FV(x_xg-1, y_xg-1) = Cov_FV(x_xg-1,y_xg-1)/(F0_V_list.err(x_xg-1)*F0_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Analysis_tag+".corr");

  Print_Cov_FA<<Cov_FA<<endl;
  Print_Cov_FV<<Cov_FV<<endl;
  Print_Corr_FA<<Corr_FA<<endl;
  Print_Corr_FV<<Corr_FV<<endl;


  Print_Cov_FA.close();
  Print_Cov_FV.close();
  Print_Corr_FA.close();
  Print_Corr_FV.close();


  
      
	

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ Up-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA(u) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################


  
  
  Vfloat ch2_FA_u;
  distr_t_list F0_u_A_list(UseJack);
  distr_t_list D1_u_A_list(UseJack);
  distr_t_list D2_u_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_u_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_u_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_A_list.distr_list.push_back(F0);
  D1_u_A_list.distr_list.push_back(D1);
  D2_u_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV(u) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_u;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_u_V_list(UseJack);
  distr_t_list D1_u_V_list(UseJack);
  distr_t_list D2_u_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_u_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_u_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_V_list.distr_list.push_back(F0);
  D1_u_V_list.distr_list.push_back(D1);
  D2_u_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_u_A_list.ave(), F0_u_A_list.err(), (D1_u_A_list/F0_u_A_list).ave(), (D1_u_A_list/F0_u_A_list).err(), (D2_u_A_list/F0_u_A_list).ave(), (D2_u_A_list/F0_u_A_list).err(), ch2_FA_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_u_V_list.ave(), F0_u_V_list.err(), (D1_u_V_list/F0_u_V_list).ave(), (D1_u_V_list/F0_u_V_list).err(), (D2_u_V_list/F0_u_V_list).ave(), (D2_u_V_list/F0_u_V_list).err(), ch2_FV_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Eigen::MatrixXd Cov_FA_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_u(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_u(x_xg-1, y_xg-1) = F0_u_A_list.distr_list[x_xg-1]%F0_u_A_list.distr_list[y_xg-1];
      Cov_FV_u(x_xg-1, y_xg-1) = F0_u_V_list.distr_list[x_xg-1]%F0_u_V_list.distr_list[y_xg-1];
      Corr_FA_u(x_xg-1, y_xg-1) = Cov_FA_u(x_xg-1,y_xg-1)/(F0_u_A_list.err(x_xg-1)*F0_u_A_list.err(y_xg-1));
      Corr_FV_u(x_xg-1, y_xg-1) = Cov_FV_u(x_xg-1,y_xg-1)/(F0_u_V_list.err(x_xg-1)*F0_u_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Analysis_tag+".corr");

  Print_Cov_FA_u<<Cov_FA_u<<endl;
  Print_Cov_FV_u<<Cov_FV_u<<endl;
  Print_Corr_FA_u<<Corr_FA_u<<endl;
  Print_Corr_FV_u<<Corr_FV_u<<endl;


  Print_Cov_FA_u.close();
  Print_Cov_FV_u.close();
  Print_Corr_FA_u.close();
  Print_Corr_FV_u.close();





  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ down-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA(d) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
 
  Vfloat ch2_FA_d;
  distr_t_list F0_d_A_list(UseJack);
  distr_t_list D1_d_A_list(UseJack);
  distr_t_list D2_d_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_d_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_d_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_A_list.distr_list.push_back(F0);
  D1_d_A_list.distr_list.push_back(D1);
  D2_d_A_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV(d) for all xg                          ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_d;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_d_V_list(UseJack);
  distr_t_list D1_d_V_list(UseJack);
  distr_t_list D2_d_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_d_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_d_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV(d), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_V_list.distr_list.push_back(F0);
  D1_d_V_list.distr_list.push_back(D1);
  D2_d_V_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_d_A_list.ave(), F0_d_A_list.err(), (D1_d_A_list/F0_d_A_list).ave(), (D1_d_A_list/F0_d_A_list).err(), (D2_d_A_list/F0_d_A_list).ave(), (D2_d_A_list/F0_d_A_list).err(), ch2_FA_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_d_V_list.ave(), F0_d_V_list.err(), (D1_d_V_list/F0_d_V_list).ave(), (D1_d_V_list/F0_d_V_list).err(), (D2_d_V_list/F0_d_V_list).ave(), (D2_d_V_list/F0_d_V_list).err(), ch2_FV_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  

  Eigen::MatrixXd Cov_FA_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_d(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_d(x_xg-1, y_xg-1) = F0_d_A_list.distr_list[x_xg-1]%F0_d_A_list.distr_list[y_xg-1];
      Cov_FV_d(x_xg-1, y_xg-1) = F0_d_V_list.distr_list[x_xg-1]%F0_d_V_list.distr_list[y_xg-1];
      Corr_FA_d(x_xg-1, y_xg-1) = Cov_FA_d(x_xg-1,y_xg-1)/(F0_d_A_list.err(x_xg-1)*F0_d_A_list.err(y_xg-1));
      Corr_FV_d(x_xg-1, y_xg-1) = Cov_FV_d(x_xg-1,y_xg-1)/(F0_d_V_list.err(x_xg-1)*F0_d_V_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Analysis_tag+".corr");

  Print_Cov_FA_d<<Cov_FA_d<<endl;
  Print_Cov_FV_d<<Cov_FV_d<<endl;
  Print_Corr_FA_d<<Corr_FA_d<<endl;
  Print_Corr_FV_d<<Corr_FV_d<<endl;


  Print_Cov_FA_d.close();
  Print_Cov_FV_d.close();
  Print_Corr_FA_d.close();
  Print_Corr_FV_d.close();




  //###################################################################################################################################################################


  

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                           FITTING FV_T AND FA_T                  ##########################################
  //###########################################                     INCLUDING SINGLE CONTRIBUTIONS              ##########################################
  //###########################################                        FROM STRANGE AND "CHARM"                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################




  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  Vfloat ch2_FA_T;
  distr_t_list F0_A_T_list(UseJack);
  distr_t_list D1_A_T_list(UseJack);
  distr_t_list D2_A_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_T_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA_T, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_A_T_list.distr_list.push_back(F0);
  D1_A_T_list.distr_list.push_back(D1);
  D2_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T for all xg                           ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

  
  Vfloat ch2_FV_T;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_V_T_list(UseJack);
  distr_t_list D1_V_T_list(UseJack);
  distr_t_list D2_V_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_T_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV_T, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_V_T_list.distr_list.push_back(F0);
  D1_V_T_list.distr_list.push_back(D1);
  D2_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_A_T_list.ave(), F0_A_T_list.err(), (D1_A_T_list/F0_A_T_list).ave(), (D1_A_T_list/F0_A_T_list).err(), (D2_A_T_list/F0_A_T_list).ave(), (D2_A_T_list/F0_A_T_list).err(), ch2_FA_T}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_V_T_list.ave(), F0_V_T_list.err(), (D1_V_T_list/F0_V_T_list).ave(), (D1_V_T_list/F0_V_T_list).err(), (D2_V_T_list/F0_V_T_list).ave(), (D2_V_T_list/F0_V_T_list).err(), ch2_FV_T}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2  ch2/dof");


  //Print covariance matrix
  Analysis_tag= ( (UseJack==true)?"jack":"boot")+_Fit_tag;
  Eigen::MatrixXd Cov_FA_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T(x_xg-1, y_xg-1) = F0_A_T_list.distr_list[x_xg-1]%F0_A_T_list.distr_list[y_xg-1];
      Cov_FV_T(x_xg-1, y_xg-1) = F0_V_T_list.distr_list[x_xg-1]%F0_V_T_list.distr_list[y_xg-1];
      Corr_FA_T(x_xg-1, y_xg-1) = Cov_FA_T(x_xg-1,y_xg-1)/(F0_A_T_list.err(x_xg-1)*F0_A_T_list.err(y_xg-1));
      Corr_FV_T(x_xg-1, y_xg-1) = Cov_FV_T(x_xg-1,y_xg-1)/(F0_V_T_list.err(x_xg-1)*F0_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_"+Analysis_tag+".corr");

  Print_Cov_FA_T<<Cov_FA_T<<endl;
  Print_Cov_FV_T<<Cov_FV_T<<endl;
  Print_Corr_FA_T<<Corr_FA_T<<endl;
  Print_Corr_FV_T<<Corr_FV_T<<endl;


  Print_Cov_FA_T.close();
  Print_Cov_FV_T.close();
  Print_Corr_FA_T.close();
  Print_Corr_FV_T.close();


  
      
	

  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ Up-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################


  
  
  Vfloat ch2_FA_T_u;
  distr_t_list F0_u_A_T_list(UseJack);
  distr_t_list D1_u_A_T_list(UseJack);
  distr_t_list D2_u_A_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_T_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_u_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_u_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_A_T_list.distr_list.push_back(F0);
  D1_u_A_T_list.distr_list.push_back(D1);
  D2_u_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T(u) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_T_u;
  bf_FF.Set_par_val("F0", 0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", 0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_u_V_T_list(UseJack);
  distr_t_list D1_u_V_T_list(UseJack);
  distr_t_list D2_u_V_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_T_u_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_u_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_u_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_u_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_u_V_T_list.distr_list.push_back(F0);
  D1_u_V_T_list.distr_list.push_back(D1);
  D2_u_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T_u.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_u_A_T_list.ave(), F0_u_A_T_list.err(), (D1_u_A_T_list/F0_u_A_T_list).ave(), (D1_u_A_T_list/F0_u_A_T_list).err(), (D2_u_A_T_list/F0_u_A_T_list).ave(), (D2_u_A_T_list/F0_u_A_T_list).err(), ch2_FA_T_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_u_V_T_list.ave(), F0_u_V_T_list.err(), (D1_u_V_T_list/F0_u_V_T_list).ave(), (D1_u_V_T_list/F0_u_V_T_list).err(), (D2_u_V_T_list/F0_u_V_T_list).ave(), (D2_u_V_T_list/F0_u_V_T_list).err(), ch2_FV_T_u}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  Eigen::MatrixXd Cov_FA_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T_u(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T_u(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T_u(x_xg-1, y_xg-1) = F0_u_A_T_list.distr_list[x_xg-1]%F0_u_A_T_list.distr_list[y_xg-1];
      Cov_FV_T_u(x_xg-1, y_xg-1) = F0_u_V_T_list.distr_list[x_xg-1]%F0_u_V_T_list.distr_list[y_xg-1];
      Corr_FA_T_u(x_xg-1, y_xg-1) = Cov_FA_T_u(x_xg-1,y_xg-1)/(F0_u_A_T_list.err(x_xg-1)*F0_u_A_T_list.err(y_xg-1));
      Corr_FV_T_u(x_xg-1, y_xg-1) = Cov_FV_T_u(x_xg-1,y_xg-1)/(F0_u_V_T_list.err(x_xg-1)*F0_u_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_u_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T_u("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_u_"+Analysis_tag+".corr");

  Print_Cov_FA_T_u<<Cov_FA_T_u<<endl;
  Print_Cov_FV_T_u<<Cov_FV_T_u<<endl;
  Print_Corr_FA_T_u<<Corr_FA_T_u<<endl;
  Print_Corr_FV_T_u<<Corr_FV_T_u<<endl;


  Print_Cov_FA_T_u.close();
  Print_Cov_FV_T_u.close();
  Print_Corr_FA_T_u.close();
  Print_Corr_FV_T_u.close();





  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################

                                                  //############ down-quark CONTRIBUTION ##################



  
    
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FA_T(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
 
  Vfloat ch2_FA_T_d;
  distr_t_list F0_d_A_T_list(UseJack);
  distr_t_list D1_d_A_T_list(UseJack);
  distr_t_list D2_d_A_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_T_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_T_d_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_T_d_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_T_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA_T(u), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_A_T_list.distr_list.push_back(F0);
  D1_d_A_T_list.distr_list.push_back(D1);
  D2_d_A_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FA_T_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

 
  //print fit func
  distr_t_list FA_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FA_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FA_T_xg_to_print.ave(), FA_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //###########################################                                                                 ##########################################
  //###########################################                                                                 ##########################################
  //###########################################                   FIT FV_T(d) for all xg                        ##########################################
  //###########################################                                                                 ##########################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  //######################################################################################################################################################
  Vfloat ch2_FV_T_d;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t_list F0_d_V_T_list(UseJack);
  distr_t_list D1_d_V_T_list(UseJack);
  distr_t_list D2_d_V_T_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_T_d_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_T_d_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data[ijack][iens].a = a_A.distr[ijack]; data[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data[ijack][iens].a = a_B.distr[ijack]; data[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data[ijack][iens].a = a_C.distr[ijack]; data[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data[ijack][iens].a = a_D.distr[ijack]; data[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_T_d_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_T_d_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") { data_ch2[ijack][iens].a = a_A.distr[ijack]; data_ch2[ijack][iens].is=0; }
		else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1;}
		else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_ch2[ijack][iens].a = a_B.distr[ijack]; data_ch2[ijack][iens].is=1; }
		else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_ch2[ijack][iens].a = a_C.distr[ijack]; data_ch2[ijack][iens].is=2; }
		else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_ch2[ijack][iens].a = a_D.distr[ijack]; data_ch2[ijack][iens].is=3; }
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV_T(d), xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack), D2(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1); D2.distr.push_back( Bt_fit.par[ijack].D2);}
  //push_back retrieved parameters
  F0_d_V_T_list.distr_list.push_back(F0);
  D1_d_V_T_list.distr_list.push_back(D1);
  D2_d_V_T_list.distr_list.push_back(D2);
  //push_back ch2
  ch2_FV_T_d.push_back( Bt_fit_ch2.get_ch2_ave()/dof);

  //print fit func
  distr_t_list FV_T_xg_to_print(UseJack);
  for(auto &a: Bs_a_to_print) FV_T_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2) + D2*pow(a*fmTGeV*Lambda_QCD,4));
  Print_To_File({}, {Bs_a_to_print, FV_T_xg_to_print.ave(), FV_T_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Fit_tag+"xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }

  //Print continuum extrapolated form factors
  Print_To_File({}, {Bs_xg_t_list, F0_d_A_T_list.ave(), F0_d_A_T_list.err(), (D1_d_A_T_list/F0_d_A_T_list).ave(), (D1_d_A_T_list/F0_d_A_T_list).err(), (D2_d_A_T_list/F0_d_A_T_list).ave(), (D2_d_A_T_list/F0_d_A_T_list).err(), ch2_FA_T_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  Print_To_File({}, {Bs_xg_t_list, F0_d_V_T_list.ave(), F0_d_V_T_list.err(), (D1_d_V_T_list/F0_d_V_T_list).ave(), (D1_d_V_T_list/F0_d_V_T_list).err(), (D2_d_V_T_list/F0_d_V_T_list).ave(), (D2_d_V_T_list/F0_d_V_T_list).err(), ch2_FV_T_d}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Fit_tag+"cont.dat", "", "#xg  F0   D1   D2   ch2/dof");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  

  Eigen::MatrixXd Cov_FA_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Cov_FV_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FA_T_d(num_xg-1, num_xg-1);
  Eigen::MatrixXd Corr_FV_T_d(num_xg-1, num_xg-1);

  for(int x_xg=1; x_xg<num_xg;x_xg++) {
    for(int y_xg=1; y_xg<num_xg;y_xg++) {
      Cov_FA_T_d(x_xg-1, y_xg-1) = F0_d_A_T_list.distr_list[x_xg-1]%F0_d_A_T_list.distr_list[y_xg-1];
      Cov_FV_T_d(x_xg-1, y_xg-1) = F0_d_V_T_list.distr_list[x_xg-1]%F0_d_V_T_list.distr_list[y_xg-1];
      Corr_FA_T_d(x_xg-1, y_xg-1) = Cov_FA_T_d(x_xg-1,y_xg-1)/(F0_d_A_T_list.err(x_xg-1)*F0_d_A_T_list.err(y_xg-1));
      Corr_FV_T_d(x_xg-1, y_xg-1) = Cov_FV_T_d(x_xg-1,y_xg-1)/(F0_d_V_T_list.err(x_xg-1)*F0_d_V_T_list.err(y_xg-1));
    }
  }

  //Print To File
  ofstream Print_Cov_FA_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Analysis_tag+".cov");
  ofstream Print_Cov_FV_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Analysis_tag+".cov");
  ofstream Print_Corr_FA_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_T_d_"+Analysis_tag+".corr");
  ofstream Print_Corr_FV_T_d("../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_T_d_"+Analysis_tag+".corr");

  Print_Cov_FA_T_d<<Cov_FA_T_d<<endl;
  Print_Cov_FV_T_d<<Cov_FV_T_d<<endl;
  Print_Corr_FA_T_d<<Corr_FA_T_d<<endl;
  Print_Corr_FV_T_d<<Corr_FV_T_d<<endl;


  Print_Cov_FA_T_d.close();
  Print_Cov_FV_T_d.close();
  Print_Corr_FA_T_d.close();
  Print_Corr_FV_T_d.close();




  //###################################################################################################################################################################

  








  

 
  //push_back FIT INFO
  rt_fit.Ch2_FA= ch2_FA;
  rt_fit.Ch2_FV= ch2_FV;
  rt_fit.FA = F0_A_list;
  rt_fit.FV = F0_V_list;
  rt_fit.Ch2_FA_T= ch2_FA_T;
  rt_fit.Ch2_FV_T= ch2_FV_T;
  rt_fit.FA_T = F0_A_T_list;
  rt_fit.FV_T = F0_V_T_list;

  rt_fit.Ch2_FA_u= ch2_FA_u;
  rt_fit.Ch2_FV_u= ch2_FV_u;
  rt_fit.FA_u = F0_u_A_list;
  rt_fit.FV_u = F0_u_V_list;
  rt_fit.Ch2_FA_T_u= ch2_FA_T_u;
  rt_fit.Ch2_FV_T_u= ch2_FV_T_u;
  rt_fit.FA_T_u = F0_u_A_T_list;
  rt_fit.FV_T_u = F0_u_V_T_list;


  rt_fit.Ch2_FA_d= ch2_FA_d;
  rt_fit.Ch2_FV_d= ch2_FV_d;
  rt_fit.FA_d = F0_d_A_list;
  rt_fit.FV_d = F0_d_V_list;
  rt_fit.Ch2_FA_T_d= ch2_FA_T_d;
  rt_fit.Ch2_FV_T_d= ch2_FV_T_d;
  rt_fit.FA_T_d = F0_d_A_T_list;
  rt_fit.FV_T_d = F0_d_V_T_list;
  

  //##########################################




  /*
  //VMD-like fit for FV and FA
  //###############################################################################################################

  class ipar_VMD {
  public:
    ipar_VMD() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, xg;
  };

  class fpar_VMD {
  public:
    fpar_VMD() {}
    fpar_VMD(const Vfloat &par) {
      if((signed)par.size() != 2) crash("In class fpar_VMD  class constructor
Vfloat par has size != 2"); A=par[0]; M=par[1];
    }
    double A, M;
  };


  //init bootstrap fit
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD(Njacks);
  bf_VMD.set_warmup_lev(0); //sets warmup
  bf_VMD.Set_number_of_measurements(Bs_xg_t_list.size());
  bf_VMD.Set_verbosity(1);
  bf_VMD.Add_par("A", 0.07, 0.001);
  bf_VMD.Add_par("M", 1.3, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD_ch2(1);
  bf_VMD_ch2.set_warmup_lev(0); //sets warmup
  bf_VMD_ch2.Set_number_of_measurements(Bs_xg_t_list.size());
  bf_VMD_ch2.Set_verbosity(1);
  bf_VMD_ch2.Add_par("A", 0.07, 0.001);
  bf_VMD_ch2.Add_par("M", 1.3, 0.1);
  //bf_VMD.Fix_par("M", 1.25);
  //bf_VMD_ch2.Fix_par("M", 1.25);


  //ansatz
  bf_VMD.ansatz=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

                    double E_res= sqrt( p.M*p.M + ip.xg*ip.xg/4.0);
                    return p.A/( E_res*( E_res - (1.0 - ip.xg/2.0)) );
                 };
  bf_VMD.measurement=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

                 return ip.FF;
                 };
  bf_VMD.error=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

                 return ip.FF_err;
                 };

  bf_VMD_ch2.ansatz= bf_VMD.ansatz;
  bf_VMD_ch2.measurement = bf_VMD.measurement;
  bf_VMD_ch2.error = bf_VMD.error;


  //start fitting FA
  //fill the data
  vector<vector<ipar_VMD>> data_VMD(Njacks);
  vector<vector<ipar_VMD>> data_VMD_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD_ch2;
  for(auto &data_iboot: data_VMD) data_iboot.resize(Bs_xg_t_list.size());
  for(auto &data_iboot: data_VMD_ch2) data_iboot.resize(Bs_xg_t_list.size());
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)Bs_xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_A_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_A_list.err(ix);
      data_VMD[ijack][ix].xg= Bs_xg_t_list[ix];
      if(ijack==0) {
        data_VMD_ch2[ijack][ix].FF = F0_A_list.ave(ix);
        data_VMD_ch2[ijack][ix].FF_err= F0_A_list.err(ix);
        data_VMD_ch2[ijack][ix].xg= Bs_xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FA using VMD ansatz"<<endl;
  Bt_fit_FA_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FA_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();
  double ch2_red_FA_VMD= Bt_fit_FA_VMD_ch2.get_ch2_ave()/( Bs_xg_t_list.size()
-2.0);

  //retrieve params
  distr_t Ampl_FA(UseJack), pole_FA(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FA.distr.push_back(
Bt_fit_FA_VMD.par[ijack].A); pole_FA.distr.push_back(
Bt_fit_FA_VMD.par[ijack].M);}


  //start fitting FV

  bf_VMD.Set_par_val("A", -0.07, 0.001);
  bf_VMD_ch2.Set_par_val("A", -0.07, 0.001);
  //bf_VMD.Fix_par("M", 1.073);
  //bf_VMD_ch2.Fix_par("M", 1.073);

  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD_ch2;

  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)Bs_xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_V_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_V_list.err(ix);
      data_VMD[ijack][ix].xg= Bs_xg_t_list[ix];
      if(ijack==0) {
        data_VMD_ch2[ijack][ix].FF = F0_V_list.ave(ix);
        data_VMD_ch2[ijack][ix].FF_err= F0_V_list.err(ix);
        data_VMD_ch2[ijack][ix].xg= Bs_xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FV using VMD ansatz"<<endl;
  Bt_fit_FV_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FV_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();

  double ch2_red_FV_VMD= Bt_fit_FV_VMD_ch2.get_ch2_ave()/( Bs_xg_t_list.size()
-2.0);


  //retrieve params
  distr_t Ampl_FV(UseJack), pole_FV(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FV.distr.push_back(
Bt_fit_FV_VMD.par[ijack].A); pole_FV.distr.push_back(
Bt_fit_FV_VMD.par[ijack].M);}


  distr_t_list FA_VMD_fit(UseJack), FV_VMD_fit(UseJack);
  //plot fit function
   for(auto &X: Bs_xg_to_spline_VMD) {
     ipar_VMD pp_VMD;
     pp_VMD.xg=X;
     distr_t FA_VMD_xg(UseJack), FV_VMD_xg(UseJack);
     for(int ijack=0;ijack<Njacks;ijack++) {
       FA_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FA_VMD.par[ijack],
pp_VMD)); FV_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FV_VMD.par[ijack],
pp_VMD));
     }

     FA_VMD_fit.distr_list.push_back( FA_VMD_xg);
     FV_VMD_fit.distr_list.push_back( FV_VMD_xg);

   }

   //print

   string header_FA= "Ampl: "+to_string_with_precision(Ampl_FA.ave(),5)+" +-
"+to_string_with_precision(Ampl_FA.err(),5)+" M^res/Mp:
"+to_string_with_precision(pole_FA.ave(), 5)+" +-
"+to_string_with_precision(pole_FA.err(), 5)+" ch2/dof:
"+to_string_with_precision(ch2_red_FA_VMD  ,5)+" dof:
"+to_string(Bs_xg_t_list.size()-2); string header_FV= "Ampl:
"+to_string_with_precision(Ampl_FV.ave(),5)+" +-
"+to_string_with_precision(Ampl_FV.err(),5)+" M^res/Mp:
"+to_string_with_precision(pole_FV.ave(), 5)+" +-
"+to_string_with_precision(pole_FV.err(), 5)+" ch2/dof:
"+to_string_with_precision(ch2_red_FV_VMD  ,5)+" dof:
"+to_string(Bs_xg_t_list.size()-2); Print_To_File({},{ Bs_xg_to_spline_VMD,
FA_VMD_fit.ave(), FA_VMD_fit.err()} ,
"../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"VMD.fit" ,
"", header_FA); Print_To_File({},{ Bs_xg_to_spline_VMD, FV_VMD_fit.ave(),
FV_VMD_fit.err()} ,
"../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"VMD.fit" ,
"", header_FV);







  //#################################################################################################################

  */

  
  
  
  
  
  //interpolate continuum extrapolated form factors
  if(num_xg > 5) {
  for(int ijack=0;ijack<Njacks;ijack++) {

    Vfloat FV_cont_ij, FA_cont_ij;
    Vfloat FV_u_cont_ij, FA_u_cont_ij;
    Vfloat FV_d_cont_ij, FA_d_cont_ij;
    for(int ixg=1;ixg<num_xg;ixg++) {
      FA_cont_ij.push_back( F0_A_list.distr_list[ixg-1].distr[ijack]);
      FV_cont_ij.push_back( F0_V_list.distr_list[ixg-1].distr[ijack]);

      //u contribution
      FA_u_cont_ij.push_back( F0_u_A_list.distr_list[ixg-1].distr[ijack]);
      FV_u_cont_ij.push_back( F0_u_V_list.distr_list[ixg-1].distr[ijack]);

      //d contribution
      FA_d_cont_ij.push_back( F0_d_A_list.distr_list[ixg-1].distr[ijack]);
      FV_d_cont_ij.push_back( F0_d_V_list.distr_list[ixg-1].distr[ijack]);
    }

    FA_cont_interpol_jacks.emplace_back( FA_cont_ij.begin(), FA_cont_ij.end(), 0.1, 0.1);
    FV_cont_interpol_jacks.emplace_back( FV_cont_ij.begin(), FV_cont_ij.end(), 0.1, 0.1);

    //u contribution
    FA_u_cont_interpol_jacks.emplace_back( FA_u_cont_ij.begin(), FA_u_cont_ij.end(), 0.1, 0.1);
    FV_u_cont_interpol_jacks.emplace_back( FV_u_cont_ij.begin(), FV_u_cont_ij.end(), 0.1, 0.1);

    //d contribution
    FA_d_cont_interpol_jacks.emplace_back( FA_d_cont_ij.begin(), FA_d_cont_ij.end(), 0.1, 0.1);
    FV_d_cont_interpol_jacks.emplace_back( FV_d_cont_ij.begin(), FV_d_cont_ij.end(), 0.1, 0.1);

  }


  auto FA_cont_interpol_distr= [&FA_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_cont_interpol_distr= [&FV_cont_interpol_jacks,  &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };


  //u contribution
  auto FA_u_cont_interpol_distr= [&FA_u_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_u_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_u_cont_interpol_distr= [&FV_u_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_u_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  

  //d contribution
  auto FA_d_cont_interpol_distr= [&FA_d_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_d_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };
  auto FV_d_cont_interpol_distr= [&FV_d_cont_interpol_jacks, &UseJack, &Njacks](double xg) -> distr_t {
    distr_t return_distr(UseJack);
    for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_d_cont_interpol_jacks[ijack](xg));}
    return return_distr;
  };

 
  distr_t_list FA_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_interpol_cont_to_print_distr(UseJack);

  //u contribution
  distr_t_list FA_u_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_u_interpol_cont_to_print_distr(UseJack);

  //d contribution
  distr_t_list FA_d_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_d_interpol_cont_to_print_distr(UseJack);
  
  for(auto &X: Bs_xg_to_spline) {
    FA_interpol_cont_to_print_distr.distr_list.push_back( FA_cont_interpol_distr(X));
    FV_interpol_cont_to_print_distr.distr_list.push_back( FV_cont_interpol_distr(X));

    //u contribution
    FA_u_interpol_cont_to_print_distr.distr_list.push_back( FA_u_cont_interpol_distr(X));
    FV_u_interpol_cont_to_print_distr.distr_list.push_back( FV_u_cont_interpol_distr(X));

    //d contribution
    FA_d_interpol_cont_to_print_distr.distr_list.push_back( FA_d_cont_interpol_distr(X));
    FV_d_interpol_cont_to_print_distr.distr_list.push_back( FV_d_cont_interpol_distr(X));
  }

  Print_To_File({}, {Bs_xg_to_spline, FA_interpol_cont_to_print_distr.ave(), FA_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {Bs_xg_to_spline, FV_interpol_cont_to_print_distr.ave(), FV_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");

  //u contribution
  Print_To_File({}, {Bs_xg_to_spline, FA_u_interpol_cont_to_print_distr.ave(), FA_u_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_u_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {Bs_xg_to_spline, FV_u_interpol_cont_to_print_distr.ave(), FV_u_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_u_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");

  //d contribution
  Print_To_File({}, {Bs_xg_to_spline, FA_d_interpol_cont_to_print_distr.ave(), FA_d_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FA_d_"+Fit_tag+"interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {Bs_xg_to_spline, FV_d_interpol_cont_to_print_distr.ave(), FV_d_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/FF/continuum/FV_d_"+Fit_tag+"interpol.dat", "", "#xg FV FV_err");

  }
  //#################################################################################






  //##########################################################################################
  //##########################################################################################

  //perform continuum extrapolation of MP, fP, and MP/fP
  //MP
  bf_FF.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  double ch2_MP;
  distr_t F0_MP(UseJack);
  distr_t D1_MP(UseJack);
  distr_t D2_MP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ch2;
  for(auto &data_iboot: data_MP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_MP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP[ijack][iens].FF= MP_list.distr_list[iens].distr[ijack];
      data_MP[ijack][iens].FF_err= MP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP[ijack][iens].a = a_A.distr[ijack]; data_MP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP[ijack][iens].a = a_B.distr[ijack]; data_MP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP[ijack][iens].a = a_B.distr[ijack]; data_MP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP[ijack][iens].a = a_C.distr[ijack]; data_MP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP[ijack][iens].a = a_D.distr[ijack]; data_MP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ch2[ijack][iens].FF= MP_list.ave(iens);
	data_MP_ch2[ijack][iens].FF_err= MP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ch2[ijack][iens].a = a_A.distr[ijack]; data_MP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ch2[ijack][iens].a = a_C.distr[ijack]; data_MP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ch2[ijack][iens].a = a_D.distr[ijack]; data_MP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP);
  bf_FF_ch2.Append_to_input_par(data_MP_ch2);
  //fit
  cout<<"Fitting MP"<<endl;
  Bt_fit_MP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ch2 = bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP.distr.push_back( Bt_fit_MP.par[ijack].F0 ); D1_MP.distr.push_back( Bt_fit_MP.par[ijack].D1 ); D2_MP.distr.push_back( Bt_fit_MP.par[ijack].D2 );}
  //get ch2
  ch2_MP= Bt_fit_MP_ch2.get_ch2_ave();


  //FP
  double ch2_FP;
  bf_FF.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t F0_FP(UseJack);
  distr_t D1_FP(UseJack);
  distr_t D2_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP_ch2;
  for(auto &data_iboot: data_FP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_FP[ijack][iens].FF= FP_list.distr_list[iens].distr[ijack];
      data_FP[ijack][iens].FF_err= FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP[ijack][iens].a = a_A.distr[ijack]; data_FP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP[ijack][iens].a = a_B.distr[ijack]; data_FP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP[ijack][iens].a = a_C.distr[ijack]; data_FP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP[ijack][iens].a = a_D.distr[ijack]; data_FP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_FP_ch2[ijack][iens].FF= FP_list.ave(iens);
	data_FP_ch2[ijack][iens].FF_err= FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_FP);
  bf_FF_ch2.Append_to_input_par(data_FP_ch2);
  //fit
  cout<<"Fitting FP"<<endl;
  Bt_fit_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_FP.distr.push_back( Bt_fit_FP.par[ijack].F0); D1_FP.distr.push_back( Bt_fit_FP.par[ijack].D1); D2_FP.distr.push_back( Bt_fit_FP.par[ijack].D2);}
  //get ch2
  ch2_FP = Bt_fit_FP_ch2.get_ch2_ave();



  //MP/FP
  double ch2_MP_ov_FP;
  bf_FF.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF.Set_par_val("D2", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  if(Include_a4) bf_FF_ch2.Set_par_val("D2", 1.0, 0.1);
  distr_t F0_MP_ov_FP(UseJack);
  distr_t D1_MP_ov_FP(UseJack);
  distr_t D2_MP_ov_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP_ch2;
  for(auto &data_iboot: data_MP_ov_FP) data_iboot.resize(Nens);
   for(auto &data_iboot: data_MP_ov_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP_ov_FP[ijack][iens].FF= MP_ov_FP_list.distr_list[iens].distr[ijack];
      data_MP_ov_FP[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP[ijack][iens].is=0; }
      else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1;}
      else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP[ijack][iens].is=1; }
      else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP[ijack][iens].is=2; }
      else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP[ijack][iens].is=3; }
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ov_FP_ch2[ijack][iens].FF= MP_ov_FP_list.ave(iens);
	data_MP_ov_FP_ch2[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") { data_MP_ov_FP_ch2[ijack][iens].a = a_A.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=0; }
	else if(data_2pts.Tag[iens] == "cB211b.072.64") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1;}
	else if(data_2pts.Tag[iens] == "cB211b.072.96") { data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=1; }
	else if(data_2pts.Tag[iens] == "cC211a.06.80")  {data_MP_ov_FP_ch2[ijack][iens].a = a_C.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=2; }
	else if(data_2pts.Tag[iens] == "cD211a.054.96")  {data_MP_ov_FP_ch2[ijack][iens].a = a_D.distr[ijack]; data_MP_ov_FP_ch2[ijack][iens].is=3; }
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP_ov_FP);
  bf_FF_ch2.Append_to_input_par(data_MP_ov_FP_ch2);
  //fit
  cout<<"Fitting MP/FP"<<endl;
  Bt_fit_MP_ov_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ov_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].F0); D1_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].D1); D2_MP_ov_FP.distr.push_back(Bt_fit_MP_ov_FP.par[ijack].D2);}
  //get ch2
  ch2_MP_ov_FP = Bt_fit_MP_ov_FP_ch2.get_ch2_ave();


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list MP_to_print(UseJack);
  distr_t_list FP_to_print(UseJack);
  distr_t_list MP_ov_FP_to_print(UseJack);
  for(auto &a: Bs_a_to_print) {
    MP_to_print.distr_list.push_back( F0_MP + D1_MP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP*pow(a*fmTGeV*Lambda_QCD,4));
    FP_to_print.distr_list.push_back( F0_FP + D1_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_FP*pow(a*fmTGeV*Lambda_QCD,4));
    MP_ov_FP_to_print.distr_list.push_back( F0_MP_ov_FP + D1_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,2) + D2_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,4));
  }
  Print_To_File({}, {Bs_a_to_print, MP_to_print.ave(), MP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/mass/masses"+_Fit_tag+".fit_func", "", "#a[fm]  MP MP_err");
  Print_To_File({}, {Bs_a_to_print, FP_to_print.ave(), FP_to_print.err(), MP_ov_FP_to_print.ave(), MP_ov_FP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+MESON+"/decay_const/fP"+_Fit_tag+".fit_func", "", "#a[fm]  FP FP_err   MP/FP  MP/FP_err");


  //push back the result for the mass and decay constant
  rt_fit.mass = F0_MP;
  rt_fit.fp = F0_FP;


  //print ch2 summary
  cout<<"#######   SUMMARY OF Ch2 ##########"<<endl;
  cout<<"### FA u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_u[ixg-1]<<endl; }
  cout<<"### FA d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_d[ixg-1]<<endl; }
  cout<<"### FV u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_u[ixg-1]<<endl; }
  cout<<"### FV d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_d[ixg-1]<<endl; }
  cout<<"### FAT u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_T_u[ixg-1]<<endl; }
  cout<<"### FAT d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FA_T_d[ixg-1]<<endl; }
  cout<<"### FVT u ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_T_u[ixg-1]<<endl; }
  cout<<"### FVT d ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<Bs_xg_t_list[ixg-1]<<"): "<<ch2_FV_T_d[ixg-1]<<endl; }
  cout<<"##########"<<endl;
  cout<<"ch2(MP): "<< ch2_MP<<endl;
  cout<<"ch2(FP): "<< ch2_FP<<endl;
  cout<<"ch2(MP/FP): "<<ch2_MP_ov_FP<<endl;
  cout<<"###################################"<<endl;

  



  
  
  //##########################################################################################
  //##########################################################################################
  
 

  }

  
 
 
  return rt_fit;




}
