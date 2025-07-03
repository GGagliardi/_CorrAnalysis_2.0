#include "../include/tau_decay_LIBE_ISO.h"


#define assertm(exp, msg) assert(((void)msg, exp))

const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const bool UseJack = 1;
const int Nboots= 2000;
const double fm_to_inv_Gev= 1.0/0.197327;
const int prec = 128;
const bool tau_LIBE_ISO_verbosity_lev=1;
const double GF= 1.1663787*1e-5; //[GeV^-2]
//CKM matrix elements
const double m_tau = 1.77686;
const double m_kappa = 0.494600000; //old is 0.4942
const double m_kappa_err = 1e-14;
const double E0_l = 0.9*m_kappa;
const double E0_sp = 0.9*m_kappa; // 0.9*(m_kappa+MPiPhys); //0.9 * (m_kappa + MPiPhys); //E0_l
const double E0_A_sp = 0.9*m_kappa; //0.9*(m_kappa+ 2*MPiPhys); // 0.9*(m_kappa+2*MPiPhys);
const double C_V = 2*M_PI/(pow(m_tau,3));
const double GAMMA_FACT= 12*M_PI; //12*M_PI*pow(Vud*GF,2);
const string MODE="TANT";
const bool Use_t_up_to_T_half_strange=false;
const int sm_func_mode= 0;
const string SM_TYPE_0= "KL_"+to_string(sm_func_mode);
const string SM_TYPE_1= "KT_"+to_string(sm_func_mode);
const bool Use_Customized_plateaux_strange=true;
using namespace std;


double Customized_plateaux_tau_spectre_LIBE_ISO( double alpha, double Emax, string channel, string reg, double s, string Ens ) {
 
  double Ra0=-1;
  int alpha_m= (int)(alpha+1);

  if(alpha_m < 3) crash("Customized_plateaux_tau_spectre called with alpha = "+to_string_with_precision(alpha,2));

  if(Ens=="cZ211a.077.64") Ens = "cB211b.072.64";

  if(alpha_m > 3) {

    if(s< 0.125) {
  
      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	  else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	  else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e4; if( s > 0.11) Ra0=9e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5; if( s >0.11) Ra0=1.6e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=6e5; if (s > 0.11) Ra0 = 3e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e5;  if(s > 0.11) Ra0=6e4;  }
	  else if(Ens == "cD211a.054.96") { Ra0=7.5e4; if(s > 0.11) Ra0=3e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=9e5; if( s > 0.11) Ra0=1e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=2e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=1.3e4;  if(s > 0.11) Ra0= 6e3; }
	  else if(Ens == "cB211b.072.96") { Ra0=4e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=2e5;  if(s > 0.11) Ra0=2e4; }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e5; if(s > 0.11) Ra0=3e4;  }
	  else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5; if( s > 0.11) Ra0=5e4;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
    }
    else { // s > 0.12


      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	  else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	  else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e4;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=3e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e4;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=2e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=2e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=2e4;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e4;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e5;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
      
      
    }
    
  }
  else {  // alpha=3

    if(s < 0.125) {

     if( reg=="tm") {
       if(channel=="Aii") {
	 if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	 else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	 else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	 else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	 else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	 else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="A0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	 else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=8e3; if(s>0.11) Ra0=4.8e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e4; if(s>0.11) Ra0=1e4;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e4; if(s>0.11) Ra0=1e4;  }
	 else if(Ens == "cD211a.054.96") { Ra0=4e4;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e5; if(s > 0.11) Ra0=1e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=8e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1.1e5;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2.5e5;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=6e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
       
     }
     else if(reg=="OS") {
       
       if(channel=="Aii") {
	 if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   }  //increase E_th
	 else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	 else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	 else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	 else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="A0") {
	 if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
    }
       else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	 else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
    }
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=4.5e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=2e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e4; if( s > 0.11) Ra0= 7e3;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e4; if( s > 0.11) Ra0=7e3;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3.5e4;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=4e4;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1.5e6;   }
	 else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
       
     }
     
     else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");

    }

    else { // sigma > 0.12

      
      if( reg=="tm") {
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") {  Ra0 =3e6;  }
	  else if(Ens == "cB211b.072.96") { Ra0= 4e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0= 6e6;   } 
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=8e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=4e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0= 2.5e7;    } 
	 else if(Ens == "cC211a.06.80") {  Ra0=3e7;  } 
	 else if(Ens == "cC211a.06.112") { Ra0= 2e6;   }
	 else if(Ens == "cD211a.054.96") {  Ra0= 2e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=3e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=1e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1e7;   }
	  else if(Ens == "cC211a.06.80") {  Ra0=1e6;  }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e6;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e7;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=3.5e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	 if(Ens == "cB211b.072.64") { Ra0=1e6;   }
	 else if(Ens == "cB211b.072.96") { Ra0=1e6;   }
	 else if(Ens == "cC211a.06.80") {  Ra0=3e6;  }
	 else if(Ens == "cC211a.06.112") { Ra0= 4e7;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2.5e7;   }
	 else if(Ens == "cE211a.044.112") {  Ra0=1e7;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
       else if(channel=="T") {
	 if(Ens == "cB211b.072.64") { Ra0=1e3;   }
	 else if(Ens == "cB211b.072.96") { Ra0=4e4;   }
	 else if(Ens == "cC211a.06.80") { Ra0=2e3;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 2e3;   }
	 else if(Ens == "cD211a.054.96") { Ra0=9e3;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e4;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else if(channel=="L") {
	 if(Ens == "cB211b.072.64") { Ra0=8e5;   }
	 else if(Ens == "cB211b.072.96") { Ra0=8e5;   }
	 else if(Ens == "cC211a.06.80") { Ra0=1e5;   }
	 else if(Ens == "cC211a.06.112") { Ra0= 1e5;   }
	 else if(Ens == "cD211a.054.96") { Ra0=2e6;  }
	 else if(Ens == "cE211a.044.112") {  Ra0=1.5e5;  }
	 else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
       }
       else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
      }
      else if(reg=="OS") {
	
	if(channel=="Aii") {
	  if(Ens == "cB211b.072.64") { Ra0=1.7e6;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.5e7;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1.2e7;   }  //increase E_th
	  else if(Ens == "cC211a.06.112") { Ra0= 6e6;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=5.5e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=8e6;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="Vii") {
	  if(Ens == "cB211b.072.64") { Ra0=2.5e6; } // Ra0= 6e5;   }
	  else if(Ens == "cB211b.072.96") {  Ra0= 9e7;    }
	  else if(Ens == "cC211a.06.80") {   Ra0= 2e7;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 2e7;   }
	  else if(Ens == "cD211a.054.96") {  Ra0=2e7;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="A0") {
	  if(Ens == "cB211b.072.64") { Ra0=2e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=1.6e6;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e6;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 1e7;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;   }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e8;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="V0") {
	  if(Ens == "cB211b.072.64") { Ra0=5e7;   }
	  else if(Ens == "cB211b.072.96") { Ra0=2e8;   }
	  else if(Ens == "cC211a.06.80") { Ra0=1e9;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 7e8;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e9;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=5e3;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="T") {
	  if(Ens == "cB211b.072.64") { Ra0=1e3;   }
	  else if(Ens == "cB211b.072.96") { Ra0=4e4;   }
	  else if(Ens == "cC211a.06.80") { Ra0=3e3;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 3e3;   }
	  else if(Ens == "cD211a.054.96") { Ra0=5e3;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=1e4;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else if(channel=="L") {
	  if(Ens == "cB211b.072.64") { Ra0=1e5;   }
	  else if(Ens == "cB211b.072.96") { Ra0=3e5;   }
	  else if(Ens == "cC211a.06.80") { Ra0=5e5;   }
	  else if(Ens == "cC211a.06.112") { Ra0= 5e5;   }
	  else if(Ens == "cD211a.054.96") { Ra0=3e6;  }
	  else if(Ens == "cE211a.044.112") {  Ra0=2e7;  }
	  else crash("In Customized_plateaux_spectre, ensemble: "+Ens+" not recognized");
	}
	else crash("In Customized_plateaux_tau_spectre, channel: "+channel+" not yet implemented");
	
     }
      
      else crash("In Customized_plateaux_tau_spectre, reg: "+reg+" not yet implemented");
    }
  }
  
  assertm(Ra0 > 0, "Assertion Ra0 > 0 failed");
  return Ra0;
}



void Generate_data() {


  bool Get_ASCII= true;

  if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../tau_decay_strange");
    

    vector<string> Ens_T1({"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"});
    
    for( int it=0; it<(signed)Ens_T1.size(); it++) {
      
      vector<string> channels({"mix_l_l",  "mix_l_s1", "mix_l_s2", "mix_s1_s1", "mix_s2_s2"});
      
      for(auto &channel : channels) {
	boost::filesystem::create_directory("../tau_decay_strange/"+channel);
	boost::filesystem::create_directory("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_S0S0", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_S0S0"});
      
      
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	
	for( auto &channel: channels) {
	  
	  FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	  size_t Nconfs, T, Nhits;
	  bin_read(Nconfs, stream);
	  bin_read(Nhits, stream);
	  bin_read(T, stream);
	  
	  cout<<"Nconfs: "<<Nconfs<<endl;
	  cout<<"T: "<<T<<" "<<T/2+1<<endl;
	  cout<<"Nhits: "<<Nhits<<endl;
	  for(size_t iconf=0;iconf<Nconfs;iconf++) {
	    vector<double> C(T/2+1);
	    for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	    boost::filesystem::create_directory("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	    ofstream PrintCorr("../tau_decay_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	    PrintCorr.precision(16);
	    PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	    for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	    if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	    else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	    PrintCorr.close();
	    
	  }
	  
	  fclose(stream);
	  
	}
	
      }
    }
  }

  bool Get_ASCII_dm=true;
  
  if(Get_ASCII_dm) {

    cout<<"Analyzing light-quark mass corrections"<<endl;
    
    //mass correction
    boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction");
    boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/mix_l1_s");
    boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/mix_l2_s");
    
    
    vector<string> Ens_T1 = {"C.06.80", "C.06.112", "B.72.64", "B.72.96" , "D.54.96", "E.44.112"};
    vector<string> Ens_TT1 = {"cC211a.06.80", "cC211a.06.112", "cB211b.072.64", "cB211b.072.96" , "cD211a.054.96", "cE211a.044.112"};
    
    
    for( int it=0; it<(signed)Ens_T1.size(); it++) {
      
      vector<string> channels({"mix_l1_s", "mix_l2_s"});
      
      for(auto &channel : channels) {
	boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel);
	boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "TM_AKAK", "TM_A0A0", "TM_V0V0", "TM_P5P5", "TM_S0S0", "OS_VKVK", "OS_AKAK", "OS_A0A0", "OS_V0V0", "OS_P5P5", "OS_S0S0"});
      
      
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	
	

	for( auto &channel: channels) {
	  
	  FILE *stream = fopen( ("../tau_decay_strange_bin_mu_corr/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	  size_t Nconfs, T, Nhits;
	  bin_read(Nconfs, stream);
	  bin_read(Nhits, stream);
	  bin_read(T, stream);
	  
	  cout<<"Nconfs: "<<Nconfs<<endl;
	  cout<<"T: "<<T<<" "<<T/2+1<<endl;
	  cout<<"Nhits: "<<Nhits<<endl;
	  for(size_t iconf=0;iconf<Nconfs;iconf++) {
	    vector<double> C(T/2+1);
	    for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	    boost::filesystem::create_directory("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	    ofstream PrintCorr("../tau_decay_strange/light_mass_correction/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	    PrintCorr.precision(16);
	    PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	    for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	    if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	    else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	    PrintCorr.close();
	    
	  }
	  
	  fclose(stream);
	  
	}
	
      }
    }
    
   

    
  }


  return;
  
}


LIBE_tau_ret Compute_tau_decay_width_LIBE_ISO(bool Is_Emax_Finite, double Emax, double alpha, double sigma, int Njacks, string Ens) {


  double beta=alpha;

  LIBE_tau_ret RET;
 
  PrecFloat::setDefaultPrecision(prec);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int _hostname_len;
  char _hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(_hostname, &_hostname_len);

   
  string Tag_reco_type="Beta_"+to_string_with_precision(beta,2);
  Tag_reco_type+="_Emax_"+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1));
  string alpha_Emax_Tag= "{"+to_string_with_precision(beta,2)+","+((Is_Emax_Finite==0)?"inf":to_string_with_precision(Emax,1))+"}";


  cout<<"STARTING COMPUTATION OF: {alpha,Emax} : "<<alpha_Emax_Tag<<endl;
  cout<<"Creating output directories...";
  
  cout.precision(5);


  //create directories

  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO");
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type);
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange");
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/Br");
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr");
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/mass");
  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance");
  

  cout<<"done!"<<endl;


  //Read data

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

   auto Sort_easy = [](string A, string B) {

      int conf_length_A= A.length();
      int conf_length_B= B.length();
      
      int pos_a_slash=-1;
      int pos_b_slash=-1;
      for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
      for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;
      
      string A_bis= A.substr(pos_a_slash+1);
      string B_bis= B.substr(pos_b_slash+1);

      return atoi( A_bis.c_str()) < atoi( B_bis.c_str());
      
  };
   

  
  //light strange mass

  data_t ll_data_tm_P5P5, ss_data_tm_P5P5, ls_data_tm_P5P5;
  data_t ll_data_tm_AKAK, ss_data_tm_AKAK, ls_data_tm_AKAK;
  data_t ll_data_tm_A0A0, ss_data_tm_A0A0, ls_data_tm_A0A0;
  data_t ll_data_tm_V0V0, ss_data_tm_V0V0, ls_data_tm_V0V0;
  data_t ll_data_tm_VKVK, ss_data_tm_VKVK, ls_data_tm_VKVK;

  data_t ll_data_OS_P5P5, ss_data_OS_P5P5, ls_data_OS_P5P5;
  data_t ll_data_OS_AKAK, ss_data_OS_AKAK, ls_data_OS_AKAK;
  data_t ll_data_OS_A0A0, ss_data_OS_A0A0, ls_data_OS_A0A0;
  data_t ll_data_OS_V0V0, ss_data_OS_V0V0, ls_data_OS_V0V0;
  data_t ll_data_OS_VKVK, ss_data_OS_VKVK, ls_data_OS_VKVK;
  
  ll_data_tm_P5P5.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_P5P5", "P5P5", Sort_easy);
  ss_data_tm_P5P5.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_P5P5", "P5P5", Sort_easy);
  ls_data_tm_P5P5.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_P5P5", "P5P5", Sort_easy);
  ll_data_tm_AKAK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_AKAK", "AKAK", Sort_easy);
  ss_data_tm_AKAK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_AKAK", "AKAK", Sort_easy);
  ls_data_tm_AKAK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_AKAK", "AKAK", Sort_easy);

  ll_data_tm_A0A0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_A0A0", "A0A0", Sort_easy);
  ss_data_tm_A0A0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_A0A0", "A0A0", Sort_easy);
  ls_data_tm_A0A0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_A0A0", "A0A0", Sort_easy);

  ll_data_tm_V0V0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_V0V0", "V0V0", Sort_easy);
  ss_data_tm_V0V0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_V0V0", "V0V0", Sort_easy);
  ls_data_tm_V0V0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_V0V0", "V0V0", Sort_easy);
  
  ll_data_tm_VKVK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_TM_VKVK", "VKVK", Sort_easy);
  ss_data_tm_VKVK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
  ls_data_tm_VKVK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_TM_VKVK", "VKVK", Sort_easy);

  ll_data_OS_P5P5.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_P5P5", "P5P5", Sort_easy);
  ss_data_OS_P5P5.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_P5P5", "P5P5", Sort_easy);
  ls_data_OS_P5P5.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_P5P5", "P5P5", Sort_easy);
  ll_data_OS_AKAK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_AKAK", "AKAK", Sort_easy);
  ss_data_OS_AKAK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_AKAK", "AKAK", Sort_easy);
  ls_data_OS_AKAK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_AKAK", "AKAK", Sort_easy);

  ll_data_OS_A0A0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_A0A0", "A0A0", Sort_easy);
  ss_data_OS_A0A0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_A0A0", "A0A0", Sort_easy);
  ls_data_OS_A0A0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_A0A0", "A0A0", Sort_easy);

  ll_data_OS_V0V0.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_V0V0", "V0V0", Sort_easy);
  ss_data_OS_V0V0.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_V0V0", "V0V0", Sort_easy);
  ls_data_OS_V0V0.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_V0V0", "V0V0", Sort_easy);
  
  ll_data_OS_VKVK.Read("../tau_decay_strange/mix_l_l", "mes_contr_mix_l_l_OS_VKVK", "VKVK", Sort_easy);
  ss_data_OS_VKVK.Read("../tau_decay_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);
  ls_data_OS_VKVK.Read("../tau_decay_strange/mix_l_s1", "mes_contr_mix_l_s1_OS_VKVK", "VKVK", Sort_easy);

  

  //heavier strange mass

  data_t  ss_H_data_tm_P5P5, ls_H_data_tm_P5P5;
  data_t  ss_H_data_tm_AKAK, ls_H_data_tm_AKAK;
  data_t  ss_H_data_tm_A0A0, ls_H_data_tm_A0A0;
  data_t  ss_H_data_tm_V0V0, ls_H_data_tm_V0V0;
  data_t  ss_H_data_tm_VKVK, ls_H_data_tm_VKVK;
  

  data_t ss_H_data_OS_P5P5, ls_H_data_OS_P5P5;
  data_t ss_H_data_OS_AKAK, ls_H_data_OS_AKAK;
  data_t ss_H_data_OS_A0A0, ls_H_data_OS_A0A0;
  data_t ss_H_data_OS_V0V0, ls_H_data_OS_V0V0;
  data_t ss_H_data_OS_VKVK, ls_H_data_OS_VKVK;


  ss_H_data_tm_P5P5.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_P5P5", "P5P5", Sort_easy);
  ls_H_data_tm_P5P5.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_P5P5", "P5P5", Sort_easy);
  ss_H_data_tm_AKAK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_AKAK", "AKAK", Sort_easy);
  ls_H_data_tm_AKAK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_AKAK", "AKAK", Sort_easy);

  ss_H_data_tm_A0A0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_A0A0", "A0A0", Sort_easy);
  ls_H_data_tm_A0A0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_A0A0", "A0A0", Sort_easy);
  
  ss_H_data_tm_V0V0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_V0V0", "V0V0", Sort_easy);
  ls_H_data_tm_V0V0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_V0V0", "V0V0", Sort_easy);
  
  ss_H_data_tm_VKVK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_TM_VKVK", "VKVK", Sort_easy);
  ls_H_data_tm_VKVK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_TM_VKVK", "VKVK", Sort_easy);


  ss_H_data_OS_P5P5.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_P5P5", "P5P5", Sort_easy);
  ls_H_data_OS_P5P5.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_P5P5", "P5P5", Sort_easy);
  ss_H_data_OS_AKAK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_AKAK", "AKAK", Sort_easy);
  ls_H_data_OS_AKAK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_AKAK", "AKAK", Sort_easy);

  ss_H_data_OS_A0A0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_A0A0", "A0A0", Sort_easy);
  ls_H_data_OS_A0A0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_A0A0", "A0A0", Sort_easy);

  
  ss_H_data_OS_V0V0.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_V0V0", "V0V0", Sort_easy);
  ls_H_data_OS_V0V0.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_V0V0", "V0V0", Sort_easy);
  
  ss_H_data_OS_VKVK.Read("../tau_decay_strange/mix_s2_s2", "mes_contr_mix_s2_s2_OS_VKVK", "VKVK", Sort_easy);
  ls_H_data_OS_VKVK.Read("../tau_decay_strange/mix_l_s2", "mes_contr_mix_l_s2_OS_VKVK", "VKVK", Sort_easy);




  //light-quark mass correction

  data_t ls_ph_data_tm_VKVK, ls_ph_data_tm_V0V0, ls_ph_data_tm_AKAK, ls_ph_data_tm_A0A0;
  data_t ls_ph_data_OS_VKVK, ls_ph_data_OS_V0V0, ls_ph_data_OS_AKAK, ls_ph_data_OS_A0A0;
  data_t ls_ph_data_tm_P5P5;


  //light 
  //tm
  ls_ph_data_tm_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_VKVK", "VKVK", Sort_easy);
  ls_ph_data_tm_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_V0V0", "V0V0", Sort_easy);
  ls_ph_data_tm_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_AKAK", "AKAK", Sort_easy);
  ls_ph_data_tm_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_A0A0", "A0A0", Sort_easy);
  ls_ph_data_tm_P5P5.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_TM_P5P5", "P5P5", Sort_easy);
  //OS
  ls_ph_data_OS_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_VKVK", "VKVK", Sort_easy);
  ls_ph_data_OS_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_V0V0", "V0V0", Sort_easy);
  ls_ph_data_OS_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_AKAK", "AKAK", Sort_easy);
  ls_ph_data_OS_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l2_s", "mes_contr_mix_l2_s_OS_A0A0", "A0A0", Sort_easy);

  


  data_t ls_uni_data_tm_VKVK, ls_uni_data_tm_V0V0, ls_uni_data_tm_AKAK, ls_uni_data_tm_A0A0;
  data_t ls_uni_data_OS_VKVK, ls_uni_data_OS_V0V0, ls_uni_data_OS_AKAK, ls_uni_data_OS_A0A0;
  data_t ls_uni_data_tm_P5P5;



  //light

  //tm
  ls_uni_data_tm_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_VKVK", "VKVK", Sort_easy);
  ls_uni_data_tm_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_V0V0", "V0V0", Sort_easy);
  ls_uni_data_tm_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_AKAK", "AKAK", Sort_easy);
  ls_uni_data_tm_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_A0A0", "A0A0", Sort_easy);
  ls_uni_data_tm_P5P5.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_TM_P5P5", "P5P5", Sort_easy);
  //OS
  ls_uni_data_OS_VKVK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_VKVK", "VKVK", Sort_easy);
  ls_uni_data_OS_V0V0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_V0V0", "V0V0", Sort_easy);
  ls_uni_data_OS_AKAK.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_AKAK", "AKAK", Sort_easy);
  ls_uni_data_OS_A0A0.Read("../tau_decay_strange/light_mass_correction/mix_l1_s", "mes_contr_mix_l1_s_OS_A0A0", "A0A0", Sort_easy);

    


  int Nens = ll_data_tm_P5P5.size;

  int Nens_mcorr= ls_ph_data_OS_A0A0.size;

  cout<<"Nens_mcorr: "<<Nens_mcorr<<endl;

  boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO_strange");
  
  



  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
   
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cB211b.072.96");
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cC211a.06.112");
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;

  
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
      a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err));
    }
  }




  //############################################################################################




  //pion masses in lattice units
  double amp_B_ave= 0.0565313; double amp_B_err= 1.438e-05;
  double amp_C_ave= 0.0472193; double amp_C_err= 3.45183e-05;
  double amp_D_ave= 0.0406214; double amp_D_err= 2.94047e-05;
  double amp_E_ave= 0.0338185; double amp_E_err= 3.18799e-05;

  distr_t aMp_B(UseJack), aMp_C(UseJack), aMp_D(UseJack), aMp_E(UseJack);

  for(int ijack=0;ijack<Njacks;ijack++ ) {

    aMp_B.distr.push_back( amp_B_ave + GM()*amp_B_err/sqrt(Njacks-1.0));
    aMp_C.distr.push_back( amp_C_ave + GM()*amp_C_err/sqrt(Njacks-1.0));
    aMp_D.distr.push_back( amp_D_ave + GM()*amp_D_err/sqrt(Njacks-1.0));
    aMp_E.distr.push_back( amp_E_ave + GM()*amp_E_err/sqrt(Njacks-1.0));

  }
  
 
  
  distr_t Mk_iso(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Mk_iso.distr.push_back( m_kappa + GM()*m_kappa_err/sqrt(Njacks-1.0));}
 
  //loop over the ensembles
  for(int iens=0; iens<Nens;iens++) {
    
    if(ll_data_tm_P5P5.Tag[iens]== Ens) {
    
      boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses");
      boost::filesystem::create_directory("../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]);
    
      CorrAnalysis Corr(UseJack, Njacks,Nboots);
      Corr.Nt = ls_data_tm_VKVK.nrows[iens];
      
   
      //effective masses
      //tm
      //pseudoscalar
      distr_t_list M_pi_tm = Corr.effective_mass_t( ll_data_tm_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_Mpi_tm");
      distr_t_list M_K_tm = Corr.effective_mass_t( ls_data_tm_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K_tm");
      distr_t_list M_K_H_tm = Corr.effective_mass_t( ls_H_data_tm_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_tm_P5P5.Tag[iens]+"/eff_mass_K_H_tm");
      
      
      //OS
      //pseudoscalar
      distr_t_list M_pi_OS = Corr.effective_mass_t( ll_data_OS_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_Mpi_OS");
      distr_t_list M_K_OS = Corr.effective_mass_t( ls_data_OS_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_OS");
      distr_t_list M_K_H_OS = Corr.effective_mass_t( ls_H_data_OS_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_H_OS");
     
      LatticeInfo L_info;
     
      L_info.LatInfo_new_ens(ls_data_tm_VKVK.Tag[iens]);
      
      double aml= L_info.ml;
      double ams1= L_info.ms_L_new;
      double ams2= L_info.ms_M_new;
      
    
      CorrAnalysis Corr_block_1(1, ls_data_tm_VKVK.Nconfs[iens], Nboots);
      Corr_block_1.Nt= ls_data_tm_VKVK.nrows[iens];
      int T = Corr.Nt;
      
     
   
     
      //get lattice spacing
      distr_t a_distr(UseJack);
      if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="B") {a_distr=a_B;}
      else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="C") {a_distr=a_C;}
      else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="D") {a_distr=a_D; ; }
      else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="E") {a_distr=a_E; }
      else crash("lattice spacing distribution for Ens: "+ls_data_tm_VKVK.Tag[iens]+" not found");


    
     //determine Mpi and Mpi^OS

     int Tmin_pion_tm=-1; int Tmin_pion_OS=-1;
     int Tmax_pion_tm=-1; int Tmax_pion_OS=-1;

     if(ll_data_OS_P5P5.Tag[iens] == "cB211b.072.64" || ll_data_OS_P5P5.Tag[iens] == "cB211b.072.96" ) {
       Tmin_pion_tm=40;   Tmax_pion_tm=55;
       Tmin_pion_OS=28;   Tmax_pion_OS=49;
     }
     else if(ll_data_OS_P5P5.Tag[iens] == "cC211a.06.80") {
       Tmin_pion_tm=43;   Tmax_pion_tm=65;
       Tmin_pion_OS=36;   Tmax_pion_OS=65;
     }
     else if(ll_data_OS_P5P5.Tag[iens] == "cC211a.06.112") {
       Tmin_pion_tm=43;   Tmax_pion_tm=65;
       Tmin_pion_OS=36;   Tmax_pion_OS=65;
     }
     else if(ll_data_OS_P5P5.Tag[iens] == "cD211a.054.96") {
       Tmin_pion_tm=45;   Tmax_pion_tm=80;
       Tmin_pion_OS=36;   Tmax_pion_OS=61;
     }
     else if(ll_data_OS_P5P5.Tag[iens] == "cE211a.044.112") {
       Tmin_pion_tm=60;   Tmax_pion_tm=90;
       Tmin_pion_OS=56;   Tmax_pion_OS=90;
     }


     Corr.Tmin= Tmin_pion_tm; Corr.Tmax= Tmax_pion_tm;
     
     distr_t MPI_tm = Corr.Fit_distr( M_pi_tm)/a_distr;

     Corr.Tmin= Tmin_pion_OS; Corr.Tmax= Tmax_pion_OS;

     distr_t MPI_OS = Corr.Fit_distr( M_pi_OS)/a_distr;


     cout<<"MPI_tm("<<ls_data_tm_VKVK.Tag[iens]<<") : "<<MPI_tm.ave()<<" +- "<<MPI_tm.err()<<endl;
     cout<<"MPI_OS("<<ls_data_tm_VKVK.Tag[iens]<<") : "<<MPI_OS.ave()<<" +- "<<MPI_OS.err()<<endl;
    
 
    //############# LIGHTER MASS ################//
    
    //tm
    distr_t_list Vk_tm_distr, Ak_tm_distr, A0_tm_distr, V0_tm_distr;
    //tm block1
    distr_t_list Vk_tm_block_1_distr, Ak_tm_block_1_distr, A0_tm_block_1_distr, V0_tm_block_1_distr;
    //OS
    distr_t_list Vk_OS_distr, Ak_OS_distr, A0_OS_distr, V0_OS_distr;
    //OS block1
    distr_t_list Vk_OS_block_1_distr,  Ak_OS_block_1_distr, A0_OS_block_1_distr, V0_OS_block_1_distr;
    distr_t_list P5_tm_distr;

      
    //light-tm sector
    Vk_tm_distr = Corr.corr_t(ls_data_tm_VKVK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Vk_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    Ak_tm_distr = Corr.corr_t(ls_data_tm_AKAK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Ak_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    A0_tm_distr = Corr.corr_t(ls_data_tm_A0A0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    V0_tm_distr = Corr.corr_t(ls_data_tm_V0V0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    P5_tm_distr = Corr.corr_t(ls_data_tm_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/P5_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat");

    //light-OS sector
    Vk_OS_distr = Corr.corr_t(ls_data_OS_VKVK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Vk_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    Ak_OS_distr = Corr.corr_t(ls_data_OS_AKAK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Ak_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    A0_OS_distr = Corr.corr_t(ls_data_OS_A0A0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");
    V0_OS_distr = Corr.corr_t(ls_data_OS_V0V0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat");


    distr_t_list MK_OS_A0_distr = Corr.effective_mass_t(A0_OS_distr, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/masses/"+ll_data_OS_P5P5.Tag[iens]+"/eff_mass_K_OS_A0");
    
    
    //light-tm sector
    Vk_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_VKVK.col(0)[iens], "");
    Ak_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_AKAK.col(0)[iens], "");
    A0_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_A0A0.col(0)[iens], "");
    V0_tm_block_1_distr = Corr_block_1.corr_t(ls_data_tm_V0V0.col(0)[iens], "");

    //light-OS sector
    Vk_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_VKVK.col(0)[iens], "");
    Ak_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_AKAK.col(0)[iens], "");
    A0_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_A0A0.col(0)[iens], "");
    V0_OS_block_1_distr = Corr_block_1.corr_t(ls_data_OS_V0V0.col(0)[iens], "");

    //############# LIGHTER MASS ################//

    //############# HEAVIER MASS ################//
    
    //tm
    distr_t_list Vk_H_tm_distr, Ak_H_tm_distr, A0_H_tm_distr, V0_H_tm_distr;
    //tm block1
    distr_t_list Vk_H_tm_block_1_distr, Ak_H_tm_block_1_distr, A0_H_tm_block_1_distr, V0_H_tm_block_1_distr;
    //OS
    distr_t_list Vk_H_OS_distr, Ak_H_OS_distr, A0_H_OS_distr, V0_H_OS_distr;
    //OS block1
    distr_t_list Vk_H_OS_block_1_distr,  Ak_H_OS_block_1_distr, A0_H_OS_block_1_distr, V0_H_OS_block_1_distr;
    distr_t_list P5_H_tm_distr;
   
    //light-tm sector
    Vk_H_tm_distr = Corr.corr_t(ls_H_data_tm_VKVK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Vk_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    Ak_H_tm_distr = Corr.corr_t(ls_H_data_tm_AKAK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Ak_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    A0_H_tm_distr = Corr.corr_t(ls_H_data_tm_A0A0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    V0_H_tm_distr = Corr.corr_t(ls_H_data_tm_V0V0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    P5_H_tm_distr = Corr.corr_t(ls_H_data_tm_P5P5.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/P5_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");

    //light-OS sector
    Vk_H_OS_distr = Corr.corr_t(ls_H_data_OS_VKVK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Vk_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    Ak_H_OS_distr = Corr.corr_t(ls_H_data_OS_AKAK.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/Ak_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    A0_H_OS_distr = Corr.corr_t(ls_H_data_OS_A0A0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
    V0_H_OS_distr = Corr.corr_t(ls_H_data_OS_V0V0.col(0)[iens], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");

  

    //analyze data with Njacks=Nconfs
    //light-tm sector
    Vk_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_VKVK.col(0)[iens], "");
    Ak_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_AKAK.col(0)[iens], "");
    A0_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_A0A0.col(0)[iens], "");
    V0_H_tm_block_1_distr = Corr_block_1.corr_t(ls_H_data_tm_V0V0.col(0)[iens], "");

   

    //light-OS sector
    Vk_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_VKVK.col(0)[iens], "");
    Ak_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_AKAK.col(0)[iens], "");
    A0_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_A0A0.col(0)[iens], "");
    V0_H_OS_block_1_distr = Corr_block_1.corr_t(ls_H_data_OS_V0V0.col(0)[iens], "");

  

    //############# HEAVIER MASS ################//


     //light-quark mass corrections

     distr_t_list V0V0_tm_dm_distr(UseJack), VKVK_tm_dm_distr(UseJack),  A0A0_tm_dm_distr(UseJack), AKAK_tm_dm_distr(UseJack);

     distr_t_list V0V0_OS_dm_distr(UseJack), VKVK_OS_dm_distr(UseJack),  A0A0_OS_dm_distr(UseJack), AKAK_OS_dm_distr(UseJack);

     distr_t_list P5P5_tm_L1(UseJack), P5P5_tm_L2(UseJack);


     //light- quark mass corrections bootstrap

     distr_t_list V0V0_tm_block_1_dm_distr(false), VKVK_tm_block_1_dm_distr(false),  A0A0_tm_block_1_dm_distr(false), AKAK_tm_block_1_dm_distr(false);

     distr_t_list V0V0_OS_block_1_dm_distr(false), VKVK_OS_block_1_dm_distr(false),  A0A0_OS_block_1_dm_distr(false), AKAK_OS_block_1_dm_distr(false);

   
     bool read_dm=false;
     int iens_dm=-1;
     for(int j=0;j<Nens_mcorr;j++) {
       if( ls_ph_data_OS_A0A0.Tag[j] == ll_data_OS_P5P5.Tag[iens]) { read_dm=true; iens_dm=j;}
     }

     P5P5_tm_L1 = Corr.corr_t( ls_ph_data_tm_P5P5.col(0)[iens_dm], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/P5P5_L1_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");
     P5P5_tm_L2 = Corr.corr_t( ls_uni_data_tm_P5P5.col(0)[iens_dm], "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/P5P5_L2_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");


     if(!read_dm) crash("Cannot find dm-corrections for ensemble: "+ll_data_OS_P5P5.Tag[iens]);

     if(read_dm) {

       cout<<"Appling dm corrections..."<<endl;
       
       V0V0_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_V0V0.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_dm_tm_"+ls_H_data_tm_V0V0.Tag[iens]+".dat");
       VKVK_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_VKVK.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/VK_dm_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat");
       A0A0_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_A0A0.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_dm_tm_"+ls_H_data_tm_A0A0.Tag[iens]+".dat");
       AKAK_tm_dm_distr = Corr.corr_t( summ_master( ls_ph_data_tm_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_AKAK.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/AK_dm_tm_"+ls_H_data_tm_AKAK.Tag[iens]+".dat");


       V0V0_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_V0V0.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/V0_dm_OS_"+ls_H_data_OS_V0V0.Tag[iens]+".dat");
       VKVK_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_VKVK.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/VK_dm_OS_"+ls_H_data_OS_VKVK.Tag[iens]+".dat");
       A0A0_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_A0A0.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/A0_dm_OS_"+ls_H_data_OS_A0A0.Tag[iens]+".dat");
       AKAK_OS_dm_distr = Corr.corr_t( summ_master( ls_ph_data_OS_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_AKAK.col(0)[iens_dm], -1.0))  , "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/corr/AK_dm_OS_"+ls_H_data_OS_AKAK.Tag[iens]+".dat");


       V0V0_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_V0V0.col(0)[iens_dm], -1.0)), "");
       VKVK_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_VKVK.col(0)[iens_dm], -1.0)), "");
       A0A0_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_A0A0.col(0)[iens_dm], -1.0)), "");
       AKAK_tm_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_tm_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_tm_AKAK.col(0)[iens_dm], -1.0)), "");


       V0V0_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_V0V0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_V0V0.col(0)[iens_dm], -1.0)), "");
       VKVK_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_VKVK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_VKVK.col(0)[iens_dm], -1.0)), "");
       A0A0_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_A0A0.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_A0A0.col(0)[iens_dm], -1.0)), "");
       AKAK_OS_block_1_dm_distr = Corr_block_1.corr_t( summ_master( ls_ph_data_OS_AKAK.col(0)[iens_dm],  Multiply_Vvector_by_scalar(ls_uni_data_OS_AKAK.col(0)[iens_dm], -1.0)), "");


    
       
       V0_tm_distr = V0_tm_distr + V0V0_tm_dm_distr;
       Vk_tm_distr = Vk_tm_distr + VKVK_tm_dm_distr;
       A0_tm_distr = A0_tm_distr + A0A0_tm_dm_distr;
       Ak_tm_distr = Ak_tm_distr + AKAK_tm_dm_distr;

       V0_OS_distr = V0_OS_distr + V0V0_OS_dm_distr;
       Vk_OS_distr = Vk_OS_distr + VKVK_OS_dm_distr;
       A0_OS_distr = A0_OS_distr + A0A0_OS_dm_distr;
       Ak_OS_distr = Ak_OS_distr + AKAK_OS_dm_distr;


       V0_tm_block_1_distr = V0_tm_block_1_distr + V0V0_tm_block_1_dm_distr;
       Vk_tm_block_1_distr = Vk_tm_block_1_distr + VKVK_tm_block_1_dm_distr;
       A0_tm_block_1_distr = A0_tm_block_1_distr + A0A0_tm_block_1_dm_distr;
       Ak_tm_block_1_distr = Ak_tm_block_1_distr + AKAK_tm_block_1_dm_distr;

       V0_OS_block_1_distr = V0_OS_block_1_distr + V0V0_OS_block_1_dm_distr;
       Vk_OS_block_1_distr = Vk_OS_block_1_distr + VKVK_OS_block_1_dm_distr;
       A0_OS_block_1_distr = A0_OS_block_1_distr + A0A0_OS_block_1_dm_distr;
       Ak_OS_block_1_distr = Ak_OS_block_1_distr + AKAK_OS_block_1_dm_distr;


       V0_H_tm_distr = V0_H_tm_distr + V0V0_tm_dm_distr;
       Vk_H_tm_distr = Vk_H_tm_distr + VKVK_tm_dm_distr;
       A0_H_tm_distr = A0_H_tm_distr + A0A0_tm_dm_distr;
       Ak_H_tm_distr = Ak_H_tm_distr + AKAK_tm_dm_distr;

       V0_H_OS_distr = V0_H_OS_distr + V0V0_OS_dm_distr;
       Vk_H_OS_distr = Vk_H_OS_distr + VKVK_OS_dm_distr;
       A0_H_OS_distr = A0_H_OS_distr + A0A0_OS_dm_distr;
       Ak_H_OS_distr = Ak_H_OS_distr + AKAK_OS_dm_distr;


       V0_H_tm_block_1_distr = V0_H_tm_block_1_distr + V0V0_tm_block_1_dm_distr;
       Vk_H_tm_block_1_distr = Vk_H_tm_block_1_distr + VKVK_tm_block_1_dm_distr;
       A0_H_tm_block_1_distr = A0_H_tm_block_1_distr + A0A0_tm_block_1_dm_distr;
       Ak_H_tm_block_1_distr = Ak_H_tm_block_1_distr + AKAK_tm_block_1_dm_distr;

       V0_H_OS_block_1_distr = V0_H_OS_block_1_distr + V0V0_OS_block_1_dm_distr;
       Vk_H_OS_block_1_distr = Vk_H_OS_block_1_distr + VKVK_OS_block_1_dm_distr;
       A0_H_OS_block_1_distr = A0_H_OS_block_1_distr + A0A0_OS_block_1_dm_distr;
       Ak_H_OS_block_1_distr = Ak_H_OS_block_1_distr + AKAK_OS_block_1_dm_distr;
       

     }
    
     
    
    


    int Tmin_P5=0;
    int Tmax_P5=0;
    distr_t aMP(UseJack);
    if( ls_data_tm_VKVK.Tag[iens] =="cB211b.072.96")     { Tmin_P5=38; Tmax_P5=59;   aMP=aMp_B;}
    else if(ls_data_tm_VKVK.Tag[iens] =="cB211b.072.64") { Tmin_P5=45; Tmax_P5=59;   aMP=aMp_B;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="C")  { Tmin_P5=36; Tmax_P5=63;   aMP=aMp_C;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="D")  { Tmin_P5=37; Tmax_P5=69;   aMP=aMp_D;}
    else if(ls_data_tm_VKVK.Tag[iens].substr(1,1)=="E")  { Tmin_P5=60; Tmax_P5=85;   aMP=aMp_E;}
    else crash("Cannot recognize the ensemble: "+ls_data_tm_VKVK.Tag[iens]+" in assigning Tmin_P5,Tmax_P5 for ensemble: ");

    Corr.Tmin = Tmin_P5; Corr.Tmax= Tmax_P5;

    distr_t MK1= Corr.Fit_distr( M_K_tm )/a_distr;
    distr_t MKOS_1 = Corr.Fit_distr(M_K_OS)/a_distr;

    distr_t MK2= Corr.Fit_distr( M_K_H_tm )/a_distr;
    distr_t MKOS_2 = Corr.Fit_distr(M_K_H_OS)/a_distr;


    vector<distr_t> MMK2({MK1*MK1, MK2*MK2});
    vector<distr_t> MMS({ Get_id_jack_distr(Njacks)*ams1, Get_id_jack_distr(Njacks)*ams2});

    distr_t Mk_iso_corr= SQRT_D(  Mk_iso*Mk_iso + 0.5*( POW_D(aMP/a_distr,2)  - pow(0.135,2)));

    distr_t ams_phys = Obs_extrapolation_meson_mass(MMS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "ams_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+".dat",  UseJack, "SPLINE" );

    vector<distr_t> MMKOS_2({MKOS_1*MKOS_1, MKOS_2*MKOS_2});
    distr_t MK_OS_phys =  SQRT_D(Obs_extrapolation_meson_mass(MMKOS_2, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "MK2_OS_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+".dat",  UseJack, "SPLINE" ));

    cout<<"MK_OS_PHYS("<<ls_data_tm_VKVK.Tag[iens]<<") : "<<MK_OS_phys.ave()<<" +- "<<MK_OS_phys.err()<<endl;
    
    double resc_GeV = C_V*GAMMA_FACT/(pow(a_distr.ave(),3));
    distr_t resc_GeV_distr= resc_GeV*Get_id_distr(Njacks,UseJack);
         
    //print covariance matrix

    //########### LIGHTER MASS #############

    Vfloat cov_A0_tm, cov_V0_tm,  cov_Ak_tm, cov_Vk_tm, cov_A0_OS, cov_V0_OS,  cov_Ak_OS, cov_Vk_OS, TT, RR;
    Vfloat corr_m_A0_tm, corr_m_V0_tm,  corr_m_Ak_tm, corr_m_Vk_tm, corr_m_A0_OS, corr_m_V0_OS, corr_m_Ak_OS, corr_m_Vk_OS;
    
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt;rr++) {
	TT.push_back(tt);
	RR.push_back(rr);

	double err_resc_A0_tm= A0_tm_distr.err(tt)*A0_tm_distr.err(rr)/(A0_tm_block_1_distr.err(tt)*A0_tm_block_1_distr.err(rr));
	double err_resc_V0_tm= V0_tm_distr.err(tt)*V0_tm_distr.err(rr)/(V0_tm_block_1_distr.err(tt)*V0_tm_block_1_distr.err(rr));
	double err_resc_Ak_tm= Ak_tm_distr.err(tt)*Ak_tm_distr.err(rr)/(Ak_tm_block_1_distr.err(tt)*Ak_tm_block_1_distr.err(rr));
	double err_resc_Vk_tm= Vk_tm_distr.err(tt)*Vk_tm_distr.err(rr)/(Vk_tm_block_1_distr.err(tt)*Vk_tm_block_1_distr.err(rr));


	double err_resc_A0_OS= A0_OS_distr.err(tt)*A0_OS_distr.err(rr)/(A0_OS_block_1_distr.err(tt)*A0_OS_block_1_distr.err(rr));
	double err_resc_V0_OS= V0_OS_distr.err(tt)*V0_OS_distr.err(rr)/(V0_OS_block_1_distr.err(tt)*V0_OS_block_1_distr.err(rr));
	double err_resc_Ak_OS= Ak_OS_distr.err(tt)*Ak_OS_distr.err(rr)/(Ak_OS_block_1_distr.err(tt)*Ak_OS_block_1_distr.err(rr));
	double err_resc_Vk_OS= Vk_OS_distr.err(tt)*Vk_OS_distr.err(rr)/(Vk_OS_block_1_distr.err(tt)*Vk_OS_block_1_distr.err(rr));

	
	cov_A0_tm.push_back( (A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr])*err_resc_A0_tm);
	cov_V0_tm.push_back( (V0_tm_block_1_distr.distr_list[tt]%V0_tm_block_1_distr.distr_list[rr])*err_resc_V0_tm);
	cov_Ak_tm.push_back( (Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr])*err_resc_Ak_tm);
	cov_Vk_tm.push_back( (Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr])*err_resc_Vk_tm);
	cov_A0_OS.push_back( (A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr])*err_resc_A0_OS);
	cov_V0_OS.push_back( (V0_OS_block_1_distr.distr_list[tt]%V0_OS_block_1_distr.distr_list[rr])*err_resc_V0_OS);
	cov_Ak_OS.push_back( (Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr])*err_resc_Ak_OS);
	cov_Vk_OS.push_back( (Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr])*err_resc_Vk_OS);


	corr_m_A0_tm.push_back( (A0_tm_block_1_distr.distr_list[tt]%A0_tm_block_1_distr.distr_list[rr])/(A0_tm_block_1_distr.err(tt)*A0_tm_block_1_distr.err(rr)));
	corr_m_V0_tm.push_back( (V0_tm_block_1_distr.distr_list[tt]%V0_tm_block_1_distr.distr_list[rr])/(V0_tm_block_1_distr.err(tt)*V0_tm_block_1_distr.err(rr)));
	corr_m_Ak_tm.push_back( (Ak_tm_block_1_distr.distr_list[tt]%Ak_tm_block_1_distr.distr_list[rr])/( Ak_tm_block_1_distr.err(tt)*Ak_tm_block_1_distr.err(rr)));
	corr_m_Vk_tm.push_back( (Vk_tm_block_1_distr.distr_list[tt]%Vk_tm_block_1_distr.distr_list[rr])/( Vk_tm_block_1_distr.err(tt)*Vk_tm_block_1_distr.err(rr)));
	corr_m_A0_OS.push_back( (A0_OS_block_1_distr.distr_list[tt]%A0_OS_block_1_distr.distr_list[rr])/( A0_OS_block_1_distr.err(tt)*A0_OS_block_1_distr.err(rr)));
	corr_m_V0_OS.push_back( (V0_OS_block_1_distr.distr_list[tt]%V0_OS_block_1_distr.distr_list[rr])/( V0_OS_block_1_distr.err(tt)*V0_OS_block_1_distr.err(rr)));
	corr_m_Ak_OS.push_back( (Ak_OS_block_1_distr.distr_list[tt]%Ak_OS_block_1_distr.distr_list[rr])/( Ak_OS_block_1_distr.err(tt)*Ak_OS_block_1_distr.err(rr)));
	corr_m_Vk_OS.push_back( (Vk_OS_block_1_distr.distr_list[tt]%Vk_OS_block_1_distr.distr_list[rr])/( Vk_OS_block_1_distr.err(tt)*Vk_OS_block_1_distr.err(rr)));
	
      }

    Print_To_File({}, {TT,RR,cov_A0_tm, corr_m_A0_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/A0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_tm, corr_m_V0_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/V0_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_tm, corr_m_Ak_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Aii_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_tm, corr_m_Vk_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Vii_tm_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_OS, corr_m_A0_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/A0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_OS, corr_m_V0_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/V0_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_OS, corr_m_Ak_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Aii_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_OS, corr_m_Vk_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Vii_OS_"+ls_data_tm_VKVK.Tag[iens]+".dat", "", "");

    

    
    //########### HEAVIER MASS #############

    Vfloat cov_A0_H_tm, cov_V0_H_tm,  cov_Ak_H_tm, cov_Vk_H_tm, cov_A0_H_OS, cov_V0_H_OS,  cov_Ak_H_OS, cov_Vk_H_OS;
    Vfloat corr_m_A0_H_tm, corr_m_V0_H_tm,  corr_m_Ak_H_tm, corr_m_Vk_H_tm, corr_m_A0_H_OS, corr_m_V0_H_OS, corr_m_Ak_H_OS, corr_m_Vk_H_OS;

       
    for(int tt=0;tt<Corr.Nt;tt++)
      for(int rr=0;rr<Corr.Nt; rr++) {


	double err_resc_A0_H_tm= A0_H_tm_distr.err(tt)*A0_H_tm_distr.err(rr)/(A0_H_tm_block_1_distr.err(tt)*A0_H_tm_block_1_distr.err(rr));
	double err_resc_V0_H_tm= V0_H_tm_distr.err(tt)*V0_H_tm_distr.err(rr)/(V0_H_tm_block_1_distr.err(tt)*V0_H_tm_block_1_distr.err(rr));
	double err_resc_Ak_H_tm= Ak_H_tm_distr.err(tt)*Ak_H_tm_distr.err(rr)/(Ak_H_tm_block_1_distr.err(tt)*Ak_H_tm_block_1_distr.err(rr));
	double err_resc_Vk_H_tm= Vk_H_tm_distr.err(tt)*Vk_H_tm_distr.err(rr)/(Vk_H_tm_block_1_distr.err(tt)*Vk_H_tm_block_1_distr.err(rr));


	double err_resc_A0_H_OS= A0_H_OS_distr.err(tt)*A0_H_OS_distr.err(rr)/(A0_H_OS_block_1_distr.err(tt)*A0_H_OS_block_1_distr.err(rr));
	double err_resc_V0_H_OS= V0_H_OS_distr.err(tt)*V0_H_OS_distr.err(rr)/(V0_H_OS_block_1_distr.err(tt)*V0_H_OS_block_1_distr.err(rr));
	double err_resc_Ak_H_OS= Ak_H_OS_distr.err(tt)*Ak_H_OS_distr.err(rr)/(Ak_H_OS_block_1_distr.err(tt)*Ak_H_OS_block_1_distr.err(rr));
	double err_resc_Vk_H_OS= Vk_H_OS_distr.err(tt)*Vk_H_OS_distr.err(rr)/(Vk_H_OS_block_1_distr.err(tt)*Vk_H_OS_block_1_distr.err(rr));
	

	cov_A0_H_tm.push_back( (A0_H_tm_block_1_distr.distr_list[tt]%A0_H_tm_block_1_distr.distr_list[rr])*err_resc_A0_H_tm);
	cov_V0_H_tm.push_back( (V0_H_tm_block_1_distr.distr_list[tt]%V0_H_tm_block_1_distr.distr_list[rr])*err_resc_V0_H_tm);
	cov_Ak_H_tm.push_back( (Ak_H_tm_block_1_distr.distr_list[tt]%Ak_H_tm_block_1_distr.distr_list[rr])*err_resc_Ak_H_tm);
	cov_Vk_H_tm.push_back( (Vk_H_tm_block_1_distr.distr_list[tt]%Vk_H_tm_block_1_distr.distr_list[rr])*err_resc_Vk_H_tm);
	cov_A0_H_OS.push_back( (A0_H_OS_block_1_distr.distr_list[tt]%A0_H_OS_block_1_distr.distr_list[rr])*err_resc_A0_H_OS);
	cov_V0_H_OS.push_back( (V0_H_OS_block_1_distr.distr_list[tt]%V0_H_OS_block_1_distr.distr_list[rr])*err_resc_V0_H_OS);
	cov_Ak_H_OS.push_back( (Ak_H_OS_block_1_distr.distr_list[tt]%Ak_H_OS_block_1_distr.distr_list[rr])*err_resc_Ak_H_OS);
	cov_Vk_H_OS.push_back( (Vk_H_OS_block_1_distr.distr_list[tt]%Vk_H_OS_block_1_distr.distr_list[rr])*err_resc_Vk_H_OS);


	corr_m_A0_H_tm.push_back( (A0_H_tm_block_1_distr.distr_list[tt]%A0_H_tm_block_1_distr.distr_list[rr])/(A0_H_tm_block_1_distr.err(tt)*A0_H_tm_block_1_distr.err(rr)));
	corr_m_V0_H_tm.push_back( (V0_H_tm_block_1_distr.distr_list[tt]%V0_H_tm_block_1_distr.distr_list[rr])/(V0_H_tm_block_1_distr.err(tt)*V0_H_tm_block_1_distr.err(rr)));
	corr_m_Ak_H_tm.push_back( (Ak_H_tm_block_1_distr.distr_list[tt]%Ak_H_tm_block_1_distr.distr_list[rr])/( Ak_H_tm_block_1_distr.err(tt)*Ak_H_tm_block_1_distr.err(rr)));
	corr_m_Vk_H_tm.push_back( (Vk_H_tm_block_1_distr.distr_list[tt]%Vk_H_tm_block_1_distr.distr_list[rr])/( Vk_H_tm_block_1_distr.err(tt)*Vk_H_tm_block_1_distr.err(rr)));
	corr_m_A0_H_OS.push_back( (A0_H_OS_block_1_distr.distr_list[tt]%A0_H_OS_block_1_distr.distr_list[rr])/( A0_H_OS_block_1_distr.err(tt)*A0_H_OS_block_1_distr.err(rr)));
	corr_m_V0_H_OS.push_back( (V0_H_OS_block_1_distr.distr_list[tt]%V0_H_OS_block_1_distr.distr_list[rr])/( V0_H_OS_block_1_distr.err(tt)*V0_H_OS_block_1_distr.err(rr)));
	corr_m_Ak_H_OS.push_back( (Ak_H_OS_block_1_distr.distr_list[tt]%Ak_H_OS_block_1_distr.distr_list[rr])/( Ak_H_OS_block_1_distr.err(tt)*Ak_H_OS_block_1_distr.err(rr)));
	corr_m_Vk_H_OS.push_back( (Vk_H_OS_block_1_distr.distr_list[tt]%Vk_H_OS_block_1_distr.distr_list[rr])/( Vk_H_OS_block_1_distr.err(tt)*Vk_H_OS_block_1_distr.err(rr)));

	
      }


    //set to zero correlation and covariance if < 0.1
   

    Print_To_File({}, {TT,RR,cov_A0_H_tm, corr_m_A0_H_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/A0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_H_tm, corr_m_V0_H_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/V0_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_H_tm, corr_m_Ak_H_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Aii_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_H_tm, corr_m_Vk_H_tm}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Vii_H_tm_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_A0_H_OS, corr_m_A0_H_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/A0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_V0_H_OS, corr_m_V0_H_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/V0_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Ak_H_OS, corr_m_Ak_H_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Aii_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
    Print_To_File({}, {TT,RR,cov_Vk_H_OS, corr_m_Vk_H_OS}, "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange/covariance/Vii_H_OS_"+ls_H_data_tm_VKVK.Tag[iens]+".dat", "", "");
       

    distr_t_list A0_tm, V0_tm,  Aii_tm, A0_OS, V0_OS, Aii_OS, Vii_tm, Vii_OS;
    distr_t_list A0_H_tm, V0_H_tm,  Aii_H_tm, A0_H_OS, V0_H_OS, Aii_H_OS, Vii_H_tm, Vii_H_OS;

    
    //######### DEFINE 0th and ii component of C^munu ###########
    //lighter
    //tm
    A0_tm = A0_tm_distr;
    V0_tm = V0_tm_distr;
    Aii_tm =Ak_tm_distr;
    Vii_tm = Vk_tm_distr;

    //OS
    A0_OS = A0_OS_distr;
    V0_OS = V0_OS_distr;
    Aii_OS = Ak_OS_distr;
    Vii_OS = Vk_OS_distr;

    //heavier
    //tm
    A0_H_tm = A0_H_tm_distr;
    V0_H_tm = V0_H_tm_distr;
    Aii_H_tm =Ak_H_tm_distr;
    Vii_H_tm = Vk_H_tm_distr;
    //OS
    A0_H_OS = A0_H_OS_distr;
    V0_H_OS = V0_H_OS_distr;
    Aii_H_OS = Ak_H_OS_distr;
    Vii_H_OS = Vk_H_OS_distr;
    //###########################################################


    //SET [tmin-tmax]
    //################   LIGHTER   ##############################

    bool Found_error_less_x_percent=false;
    double x=10;
    //tm
    int tmax_tm_A0=1;
    while(!Found_error_less_x_percent && tmax_tm_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_tm.distr_list[tmax_tm_A0]).err()/fabs( (A0_tm.distr_list[tmax_tm_A0]).ave()) <  0.01*x) tmax_tm_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_tm_V0=1;
    while(!Found_error_less_x_percent && tmax_tm_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_tm.distr_list[tmax_tm_V0]).err()/fabs( (V0_tm.distr_list[tmax_tm_V0]).ave()) <  0.01*x) tmax_tm_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    

    int tmax_tm_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_tm_1_Aii < Corr.Nt/2 -1 ) {
   
      if( (Aii_tm.distr_list[tmax_tm_1_Aii]).err()/fabs( (Aii_tm.distr_list[tmax_tm_1_Aii]).ave()) <  0.01*x) tmax_tm_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_tm_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_tm_1_Vii < Corr.Nt/2 -1 ) {
   
      if( (Vii_tm.distr_list[tmax_tm_1_Vii]).err()/fabs( (Vii_tm.distr_list[tmax_tm_1_Vii]).ave()) <  0.01*x) tmax_tm_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    //#########################################################################
    //OS

    int tmax_OS_A0=1;
      while(!Found_error_less_x_percent && tmax_OS_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_OS.distr_list[tmax_OS_A0]).err()/fabs( (A0_OS.distr_list[tmax_OS_A0]).ave()) <  0.01*x) tmax_OS_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_OS_V0=1;
    while(!Found_error_less_x_percent && tmax_OS_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_OS.distr_list[tmax_OS_V0]).err()/fabs( (V0_OS.distr_list[tmax_OS_V0]).ave()) <  0.01*x) tmax_OS_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_OS_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_OS_1_Aii < Corr.Nt/2 -1) {
   
      if( (Aii_OS.distr_list[tmax_OS_1_Aii]).err()/fabs( (Aii_OS.distr_list[tmax_OS_1_Aii]).ave()) <  0.01*x) tmax_OS_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_OS_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_OS_1_Vii < Corr.Nt/2 -1) {
   
      if( (Vii_OS.distr_list[tmax_OS_1_Vii]).err()/fabs( (Vii_OS.distr_list[tmax_OS_1_Vii]).ave()) <  0.01*x) tmax_OS_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    

    //################   HEAVIER   ##############################
    //tm
    int tmax_H_tm_A0=1;
    while(!Found_error_less_x_percent && tmax_H_tm_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_H_tm.distr_list[tmax_H_tm_A0]).err()/fabs( (A0_H_tm.distr_list[tmax_H_tm_A0]).ave()) <  0.01*x) tmax_H_tm_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
    int tmax_H_tm_V0=1;
    while(!Found_error_less_x_percent && tmax_H_tm_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_H_tm.distr_list[tmax_H_tm_V0]).err()/fabs( (V0_H_tm.distr_list[tmax_H_tm_V0]).ave()) <  0.01*x) tmax_H_tm_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    

    int tmax_H_tm_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_H_tm_1_Aii < Corr.Nt/2 -1 ) {
   
      if( (Aii_H_tm.distr_list[tmax_H_tm_1_Aii]).err()/fabs( (Aii_H_tm.distr_list[tmax_H_tm_1_Aii]).ave()) <  0.01*x) tmax_H_tm_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_tm_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_H_tm_1_Vii < Corr.Nt/2 -1 ) {
   
      if( (Vii_H_tm.distr_list[tmax_H_tm_1_Vii]).err()/fabs( (Vii_H_tm.distr_list[tmax_H_tm_1_Vii]).ave()) <  0.01*x) tmax_H_tm_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    //OS

    int tmax_H_OS_A0=1;
    while(!Found_error_less_x_percent && tmax_H_OS_A0 < Corr.Nt/2 -1 ) {
   
      if( (A0_H_OS.distr_list[tmax_H_OS_A0]).err()/fabs( (A0_H_OS.distr_list[tmax_H_OS_A0]).ave()) <  0.01*x) tmax_H_OS_A0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;


    int tmax_H_OS_V0=1;
    while(!Found_error_less_x_percent && tmax_H_OS_V0 < Corr.Nt/2 -1 ) {
   
      if( (V0_H_OS.distr_list[tmax_H_OS_V0]).err()/fabs( (V0_H_OS.distr_list[tmax_H_OS_V0]).ave()) <  0.01*x) tmax_H_OS_V0++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_OS_1_Aii=1;
    while(!Found_error_less_x_percent && tmax_H_OS_1_Aii < Corr.Nt/2 -1) {
   
      if( (Aii_H_OS.distr_list[tmax_H_OS_1_Aii]).err()/fabs( (Aii_H_OS.distr_list[tmax_H_OS_1_Aii]).ave()) <  0.01*x) tmax_H_OS_1_Aii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;

    int tmax_H_OS_1_Vii=1;
    while(!Found_error_less_x_percent && tmax_H_OS_1_Vii < Corr.Nt/2 -1) {
   
      if( (Vii_H_OS.distr_list[tmax_H_OS_1_Vii]).err()/fabs( (Vii_H_OS.distr_list[tmax_H_OS_1_Vii]).ave()) <  0.01*x) tmax_H_OS_1_Vii++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    
   
      
    const auto K0 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

      
		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      PrecFloat XX= E/(m_tau*a_distr.ave());
		      if(X_ave < E0) return 0.0;
		      
		      //if(X_ave > 1e4) return 0.0;

		  

                      if(ijack == -1)  X=X_ave;
		      else 	X= E/(m_tau*a_distr.distr[ijack]);
		      
		      
		      PrecFloat sm_theta;

		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
						 
		      return (1/XX)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    auto K0_dub = [&a_distr](double &E, double &s) -> double {

      
		      
      double X = E/(m_tau*a_distr.ave());
      double XX= E/(m_tau*a_distr.ave());
		     
          
      double sm_theta;
      
      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
      
      return (1/XX)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
    };

    const auto K1 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*a_distr.ave());
		      PrecFloat XX= E/(m_tau*a_distr.ave());
		      if( X_ave < E0) return 0.0;

		      //if(X_ave > 1e4) return 0.0;

		  

		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);

		      
		      PrecFloat sm_theta;
		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
		      
		      return (1 + 2*pow(X,2))*(1/(XX))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };

    auto K1_dub = [&a_distr](double E, double s) -> double {

      
		      
      double X = E/(m_tau*a_distr.ave());
      double XX= E/(m_tau*a_distr.ave());
      
      
      double sm_theta;

           
      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");

      return (1 + 2*pow(X,2))*(1/(XX))*pow(( 1 -pow(X,2)),2)*sm_theta;
     		   
    };



    const auto K0_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

      
		      
		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
		      if(X_ave < E0) return 0.0;
		      
		      //if(X_ave > 1e4) return 0.0;

		      PrecFloat FF= pow(a_distr.ave()/(a_distr.ave()+a_distr.err()),3);


                      if(ijack == -1)  X=X_ave;
		      else 	X= E/(m_tau*a_distr.distr[ijack]);
		      
		      
		      PrecFloat sm_theta;

		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
						 
		      return FF*(1/X)*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		 };

    const auto K1_shifted = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

		      PrecFloat X;
		      PrecFloat X_ave = E/(m_tau*(a_distr.ave() + a_distr.err()));
		      if( X_ave < E0) return 0.0;

		      //if(X_ave > 1e4) return 0.0;

		      PrecFloat FF= pow(a_distr.ave()/(a_distr.ave()+a_distr.err()),3);


		      if(ijack==-1) {
			X=X_ave;
		      }
		      else X= E/(m_tau*a_distr.distr[ijack]);

		      
		      PrecFloat sm_theta;
		      if(sm_func_mode==0) sm_theta= 1/(1+ exp(-(1-X)/s));
		      else if(sm_func_mode==1) sm_theta= 1/(1+ exp(-sinh((1-X)/s)));
		      else if(sm_func_mode==2) sm_theta= (1+erf((1-X)/s))/2;
		      else crash("sm_func_mode: "+to_string(sm_func_mode)+" not yet implemented");
		      
		      return FF*(1 + 2*pow(X,2))*(1/(X))*pow(( 1 -pow(X,2)),2)*sm_theta;
		   
		    };
    


    if(Use_t_up_to_T_half_strange) {
      tmax_tm_1_Aii = tmax_tm_1_Vii= tmax_tm_A0 = tmax_OS_A0 = tmax_tm_V0 = tmax_OS_V0 = tmax_OS_1_Aii = tmax_OS_1_Vii= Corr.Nt/2 -1;
      tmax_H_tm_1_Aii = tmax_H_tm_1_Vii= tmax_H_tm_A0 = tmax_H_OS_A0 = tmax_H_tm_V0 = tmax_H_OS_V0 = tmax_H_OS_1_Aii = tmax_H_OS_1_Vii= Corr.Nt/2 -1;

    }

  
    double s= sigma;

    
    //#####    LIGHTER  ##########
    distr_t Br_sigma_A0_tm;
    distr_t Br_sigma_V0_tm;
    distr_t Br_sigma_Aii_tm;
    distr_t Br_sigma_Vii_tm;
    distr_t Br_sigma_A0_OS;
    distr_t Br_sigma_V0_OS;
    distr_t Br_sigma_Aii_OS;
    distr_t Br_sigma_Vii_OS;

    
    distr_t Br_s_sigma_A0_tm;
    distr_t Br_s_sigma_V0_tm;
    distr_t Br_s_sigma_Aii_tm;
    distr_t Br_s_sigma_Vii_tm;
    distr_t Br_s_sigma_A0_OS;
    distr_t Br_s_sigma_V0_OS;
    distr_t Br_s_sigma_Aii_OS;
    distr_t Br_s_sigma_Vii_OS;
    
    //int tmax= T/2 -4;
    double lA0_tm, lAii_tm, lVii_tm, lV0_tm;
    double lA0_OS, lAii_OS, lVii_OS, lV0_OS;
    
    double slA0_tm, slAii_tm, slVii_tm, slV0_tm;
    double slA0_OS, slAii_OS, slVii_OS, slV0_OS;
          
    double syst_A0_tm, syst_Aii_tm, syst_Vii_tm, syst_V0_tm;
    double syst_A0_OS, syst_Aii_OS, syst_Vii_OS, syst_V0_OS;
    
    double syst_s_A0_tm, syst_s_Aii_tm, syst_s_Vii_tm, syst_s_V0_tm;
    double syst_s_A0_OS, syst_s_Aii_OS, syst_s_Vii_OS, syst_s_V0_OS;
    
    
    //#####    HEAVIER  ##########
    distr_t Br_sigma_A0_H_tm;
    distr_t Br_sigma_V0_H_tm;
    distr_t Br_sigma_Aii_H_tm;
    distr_t Br_sigma_Vii_H_tm;
    distr_t Br_sigma_A0_H_OS;
    distr_t Br_sigma_V0_H_OS;
    distr_t Br_sigma_Aii_H_OS;
    distr_t Br_sigma_Vii_H_OS;
    
    distr_t Br_s_sigma_A0_H_tm;
    distr_t Br_s_sigma_V0_H_tm;
    distr_t Br_s_sigma_Aii_H_tm;
    distr_t Br_s_sigma_Vii_H_tm;
    distr_t Br_s_sigma_A0_H_OS;
    distr_t Br_s_sigma_V0_H_OS;
    distr_t Br_s_sigma_Aii_H_OS;
    distr_t Br_s_sigma_Vii_H_OS;
    
    //int tmax= T/2 -4;
    double lA0_H_tm, lAii_H_tm, lVii_H_tm, lV0_H_tm;
    double lA0_H_OS, lAii_H_OS, lVii_H_OS, lV0_H_OS;
    
    double slA0_H_tm, slAii_H_tm, slVii_H_tm, slV0_H_tm;
    double slA0_H_OS, slAii_H_OS, slVii_H_OS, slV0_H_OS;
          
    double syst_A0_H_tm, syst_Aii_H_tm, syst_Vii_H_tm, syst_V0_H_tm;
    double syst_A0_H_OS, syst_Aii_H_OS, syst_Vii_H_OS, syst_V0_H_OS;
    
    double syst_s_A0_H_tm, syst_s_Aii_H_tm, syst_s_Vii_H_tm, syst_s_V0_H_tm;
    double syst_s_A0_H_OS, syst_s_Aii_H_OS, syst_s_Vii_H_OS, syst_s_V0_H_OS;
          
    double mult=1e4;
    if( (beta > 2) && (s < 0.15) ) mult=1e5;

    
    //###########################################################//
    //###########################################################//
    //###########################################################//
    //###########################################################//
    //##############                              ###############//
    //##############                              ###############//
    //##############                              ###############//
    //##############         LIGHTER              ###############//
    //##############                              ###############//
    //##############                              ###############//
    //##############                              ###############//
      
      auto start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax,  "Aii", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_tm, syst_Aii_tm, mult, lAii_tm, MODE, "tm", "Aii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_tm_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_tm, syst_s_Aii_tm, mult, slAii_tm, MODE, "tm", "Aii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      
      syst_s_Aii_tm = fabs( Br_s_sigma_Aii_tm.ave() - Br_sigma_Aii_tm.ave());
      Br_sigma_Aii_tm= Br_sigma_Aii_tm.ave() + (Br_sigma_Aii_tm-Br_sigma_Aii_tm.ave())*sqrt( pow(Br_sigma_Aii_tm.err(),2) + pow(syst_s_Aii_tm,2))/Br_sigma_Aii_tm.err();
      auto end = chrono::system_clock::now();
      cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush;
      chrono::duration<double> elapsed_seconds = end-start;
      double time_Aii_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

      
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "Aii", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_OS, syst_Aii_OS, mult, lAii_OS, MODE, "OS", "Aii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_OS_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_OS, syst_s_Aii_OS, mult, slAii_OS, MODE, "OS", "Aii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);

      syst_s_Aii_OS = fabs( Br_s_sigma_Aii_OS.ave() - Br_sigma_Aii_OS.ave());
      Br_sigma_Aii_OS= Br_sigma_Aii_OS.ave() + (Br_sigma_Aii_OS-Br_sigma_Aii_OS.ave())*sqrt( pow(Br_sigma_Aii_OS.err(),2) + pow(syst_s_Aii_OS,2))/Br_sigma_Aii_OS.err();

           
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Aii_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "Vii", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_tm, syst_Vii_tm, mult, lVii_tm, MODE, "tm", "Vii_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_tm_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_tm, syst_s_Vii_tm, mult, slVii_tm, MODE, "tm", "Vii_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
     
      syst_s_Vii_tm = fabs( Br_s_sigma_Vii_tm.ave() - Br_sigma_Vii_tm.ave());
      Br_sigma_Vii_tm= Br_sigma_Vii_tm.ave() + (Br_sigma_Vii_tm-Br_sigma_Vii_tm.ave())*sqrt( pow(Br_sigma_Vii_tm.err(),2) + pow(syst_s_Vii_tm,2))/Br_sigma_Vii_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      
      
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax,  "Vii", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_OS, syst_Vii_OS, mult, lVii_OS, MODE, "OS", "Vii_strange_"+ls_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_OS, fake_func,0, fake_func_d  ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_OS_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_OS, syst_s_Vii_OS, mult, slVii_OS, MODE, "OS", "Vii_s_strange_"+ls_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);

      syst_s_Vii_OS = fabs( Br_s_sigma_Vii_OS.ave() - Br_sigma_Vii_OS.ave());
      Br_sigma_Vii_OS= Br_sigma_Vii_OS.ave() + (Br_sigma_Vii_OS-Br_sigma_Vii_OS.ave())*sqrt( pow(Br_sigma_Vii_OS.err(),2) + pow(syst_s_Vii_OS,2))/Br_sigma_Vii_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "A0", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_A0, prec, SM_TYPE_0,K0, A0_tm, syst_A0_tm, mult, lA0_tm, MODE, "tm", "A0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_A0, prec, SM_TYPE_0,K0_shifted, A0_tm, syst_s_A0_tm, mult, slA0_tm, MODE, "tm", "A0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_tm = fabs( Br_s_sigma_A0_tm.ave() - Br_sigma_A0_tm.ave());
      Br_sigma_A0_tm= Br_sigma_A0_tm.ave() + (Br_sigma_A0_tm-Br_sigma_A0_tm.ave())*sqrt( pow(Br_sigma_A0_tm.err(),2) + pow(syst_s_A0_tm,2))/Br_sigma_A0_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "A0", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_A0, prec, SM_TYPE_0,K0, A0_OS, syst_A0_OS, mult, lA0_OS, MODE, "OS", "A0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_A0, prec, SM_TYPE_0,K0_shifted, A0_OS, syst_s_A0_OS, mult, slA0_OS, MODE, "OS", "A0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_OS = fabs( Br_s_sigma_A0_OS.ave() - Br_sigma_A0_OS.ave());
      Br_sigma_A0_OS= Br_sigma_A0_OS.ave() + (Br_sigma_A0_OS-Br_sigma_A0_OS.ave())*sqrt( pow(Br_sigma_A0_OS.err(),2) + pow(syst_s_A0_OS,2))/Br_sigma_A0_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "V0", "tm" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_V0, prec, SM_TYPE_0,K0, V0_tm, syst_V0_tm, mult, lV0_tm, MODE, "tm", "V0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_tm_V0, prec, SM_TYPE_0,K0_shifted, V0_tm, syst_s_V0_tm, mult, slV0_tm, MODE, "tm", "V0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );

      syst_s_V0_tm = fabs( Br_s_sigma_V0_tm.ave() - Br_sigma_V0_tm.ave());
      Br_sigma_V0_tm= Br_sigma_V0_tm.ave() + (Br_sigma_V0_tm-Br_sigma_V0_tm.ave())*sqrt( pow(Br_sigma_V0_tm.err(),2) + pow(syst_s_V0_tm,2))/Br_sigma_V0_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_tm, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "V0", "OS" , s, ls_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_V0, prec, SM_TYPE_0,K0, V0_OS, syst_V0_OS, mult, lV0_OS, MODE, "OS", "V0_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_OS_V0, prec, SM_TYPE_0,K0_shifted, V0_OS, syst_s_V0_OS, mult, slV0_OS, MODE, "OS", "V0_s_strange_"+ls_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_OS = fabs( Br_s_sigma_V0_OS.ave() - Br_sigma_V0_OS.ave());
      Br_sigma_V0_OS= Br_sigma_V0_OS.ave() + (Br_sigma_V0_OS-Br_sigma_V0_OS.ave())*sqrt( pow(Br_sigma_V0_OS.err(),2) + pow(syst_s_V0_OS,2))/Br_sigma_V0_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_OS, sigma: "<<s<<", Ens: "<<ls_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      
      start = chrono::system_clock::now();

      //###########################################################//
      //###########################################################//
      //###########################################################//
      //###########################################################//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############         HEAVIER              ###############//
      //##############                              ###############//
      //##############                              ###############//
      //##############                              ###############//
  
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax,  "Aii", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_tm_1_Aii, prec, SM_TYPE_1,K1, Aii_H_tm, syst_Aii_H_tm, mult, lAii_H_tm, MODE, "tm", "Aii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_tm_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_H_tm, syst_s_Aii_H_tm, mult, slAii_H_tm, MODE, "tm", "Aii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_H_tm = fabs( Br_s_sigma_Aii_H_tm.ave() - Br_sigma_Aii_H_tm.ave());
      Br_sigma_Aii_H_tm= Br_sigma_Aii_H_tm.ave() + (Br_sigma_Aii_H_tm-Br_sigma_Aii_H_tm.ave())*sqrt( pow(Br_sigma_Aii_H_tm.err(),2) + pow(syst_s_Aii_H_tm,2))/Br_sigma_Aii_H_tm.err();
      end = chrono::system_clock::now();
      cout<<"node: "<<_hostname<<", rank: "<<rank<<", thread_id: "<<omp_get_thread_num()<<" core-id: "<<sched_getcpu()<<endl<<flush;
      elapsed_seconds = end-start;
      double time_Aii_H_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

        
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "Aii", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Aii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_OS_1_Aii, prec, SM_TYPE_1,K1, Aii_H_OS, syst_Aii_H_OS, mult, lAii_H_OS, MODE, "OS", "Aii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Aii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_A_sp*a_distr.ave(),  T, tmax_H_OS_1_Aii, prec, SM_TYPE_1,K1_shifted, Aii_H_OS, syst_s_Aii_H_OS, mult, slAii_H_OS, MODE, "OS", "Aii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0,resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Ak_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Aii_H_OS = fabs( Br_s_sigma_Aii_H_OS.ave() - Br_sigma_Aii_H_OS.ave());
      Br_sigma_Aii_H_OS= Br_sigma_Aii_H_OS.ave() + (Br_sigma_Aii_H_OS-Br_sigma_Aii_H_OS.ave())*sqrt( pow(Br_sigma_Aii_H_OS.err(),2) + pow(syst_s_Aii_H_OS,2))/Br_sigma_Aii_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Aii_H_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Aii_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Aii_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
        
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "Vii", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_1_Vii, prec, SM_TYPE_1,K1, Vii_H_tm, syst_Vii_H_tm, mult, lVii_H_tm, MODE, "tm", "Vii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_H_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_H_tm = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_tm_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_H_tm, syst_s_Vii_H_tm, mult, slVii_H_tm, MODE, "tm", "Vii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_H_tm, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, beta);
      syst_s_Vii_H_tm = fabs( Br_s_sigma_Vii_H_tm.ave() - Br_sigma_Vii_H_tm.ave());
      Br_sigma_Vii_H_tm= Br_sigma_Vii_H_tm.ave() + (Br_sigma_Vii_H_tm-Br_sigma_Vii_H_tm.ave())*sqrt( pow(Br_sigma_Vii_H_tm.err(),2) + pow(syst_s_Vii_H_tm,2))/Br_sigma_Vii_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_H_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;

     
      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax,  "Vii", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_Vii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_1_Vii, prec, SM_TYPE_1,K1, Vii_H_OS, syst_Vii_H_OS, mult, lVii_H_OS, MODE, "OS", "Vii_strange_H_"+ls_H_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_H_OS, fake_func,0, fake_func_d  ,  Is_Emax_Finite, Emax, beta);
      Br_s_sigma_Vii_H_OS = Get_Laplace_transfo(  0.0,  s, E0_sp*a_distr.ave(),  T, tmax_H_OS_1_Vii, prec, SM_TYPE_1,K1_shifted, Vii_H_OS, syst_s_Vii_H_OS, mult, slVii_H_OS, MODE, "OS", "Vii_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens],1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_Vk_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta);
      syst_s_Vii_H_OS = fabs( Br_s_sigma_Vii_H_OS.ave() - Br_sigma_Vii_H_OS.ave());
      Br_sigma_Vii_H_OS= Br_sigma_Vii_H_OS.ave() + (Br_sigma_Vii_H_OS-Br_sigma_Vii_H_OS.ave())*sqrt( pow(Br_sigma_Vii_H_OS.err(),2) + pow(syst_s_Vii_H_OS,2))/Br_sigma_Vii_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_Vii_H_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[Vii_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_Vii_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "A0", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_A0, prec, SM_TYPE_0,K0, A0_H_tm, syst_A0_H_tm, mult, lA0_H_tm, MODE, "tm", "A0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_A0, prec, SM_TYPE_0,K0_shifted, A0_H_tm, syst_s_A0_H_tm, mult, slA0_H_tm, MODE, "tm", "A0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_H_tm = fabs( Br_s_sigma_A0_H_tm.ave() - Br_sigma_A0_H_tm.ave());
      Br_sigma_A0_H_tm= Br_sigma_A0_H_tm.ave() + (Br_sigma_A0_H_tm-Br_sigma_A0_H_tm.ave())*sqrt( pow(Br_sigma_A0_H_tm.err(),2) + pow(syst_s_A0_H_tm,2))/Br_sigma_A0_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_H_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "A0", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_A0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_A0, prec, SM_TYPE_0,K0, A0_H_OS, syst_A0_H_OS, mult, lA0_H_OS, MODE, "OS", "A0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_A0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_A0, prec, SM_TYPE_0,K0_shifted, A0_H_OS, syst_s_A0_H_OS, mult, slA0_H_OS, MODE, "OS", "A0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_A0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_A0_H_OS = fabs( Br_s_sigma_A0_H_OS.ave() - Br_sigma_A0_H_OS.ave());
      Br_sigma_A0_H_OS= Br_sigma_A0_H_OS.ave() + (Br_sigma_A0_H_OS-Br_sigma_A0_H_OS.ave())*sqrt( pow(Br_sigma_A0_H_OS.err(),2) + pow(syst_s_A0_H_OS,2))/Br_sigma_A0_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_A0_H_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[A0_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_A0_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;



      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "V0", "tm" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_V0, prec, SM_TYPE_0,K0, V0_H_tm, syst_V0_H_tm, mult, lV0_H_tm, MODE, "tm", "V0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_H_tm = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_tm_V0, prec, SM_TYPE_0,K0_shifted, V0_H_tm, syst_s_V0_H_tm, mult, slV0_H_tm, MODE, "tm", "V0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_H_tm, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_H_tm = fabs( Br_s_sigma_V0_H_tm.ave() - Br_sigma_V0_H_tm.ave());
      Br_sigma_V0_H_tm= Br_sigma_V0_H_tm.ave() + (Br_sigma_V0_H_tm-Br_sigma_V0_H_tm.ave())*sqrt( pow(Br_sigma_V0_H_tm.err(),2) + pow(syst_s_V0_H_tm,2))/Br_sigma_V0_H_tm.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_H_tm= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_H_tm, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_H_tm<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;
      

      start = chrono::system_clock::now();
      if( (beta > 2) && Use_Customized_plateaux_strange) mult=  Customized_plateaux_tau_spectre_LIBE_ISO( beta, Emax, "V0", "OS" , s, ls_H_data_tm_VKVK.Tag[iens] );
      Br_sigma_V0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_V0, prec, SM_TYPE_0,K0, V0_H_OS, syst_V0_H_OS, mult, lV0_H_OS, MODE, "OS", "V0_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      Br_s_sigma_V0_H_OS = Get_Laplace_transfo(  0.0,  s, E0_l*a_distr.ave(),  T, tmax_H_OS_V0, prec, SM_TYPE_0,K0_shifted, V0_H_OS, syst_s_V0_H_OS, mult, slV0_H_OS, MODE, "OS", "V0_s_strange_H_"+ls_H_data_tm_VKVK.Tag[iens], 1e-3,0, resc_GeV, 0.0, "tau_decay_LIBE_ISO", cov_V0_H_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, beta );
      syst_s_V0_H_OS = fabs( Br_s_sigma_V0_H_OS.ave() - Br_sigma_V0_H_OS.ave());
      Br_sigma_V0_H_OS= Br_sigma_V0_H_OS.ave() + (Br_sigma_V0_H_OS-Br_sigma_V0_H_OS.ave())*sqrt( pow(Br_sigma_V0_H_OS.err(),2) + pow(syst_s_V0_H_OS,2))/Br_sigma_V0_H_OS.err();
      end = chrono::system_clock::now();
      elapsed_seconds = end-start;
      double time_V0_H_OS= elapsed_seconds.count();
      if(tau_LIBE_ISO_verbosity_lev) {
	cout<<endl<<flush;
	cout<<"Elapsed time[V0_H_OS, sigma: "<<s<<", Ens: "<<ls_H_data_tm_VKVK.Tag[iens]<<", rank: "<<rank<<", #thread="<<omp_get_thread_num()<<"] : "<<time_V0_H_OS<<" s"<<endl<<flush;
      }
      else cout<<"."<<flush;


      

      start = chrono::system_clock::now();
      

      //########### EXTRAPOLATE TO THE PHYSICAL STRANGE POINT ################
      
      //tm
      vector<distr_t> BR_LIST_Vii_tm({ Br_sigma_Vii_tm, Br_sigma_Vii_H_tm});
      vector<distr_t> BR_LIST_Aii_tm({ Br_sigma_Aii_tm, Br_sigma_Aii_H_tm});
      vector<distr_t> BR_LIST_A0_tm({ Br_sigma_A0_tm, Br_sigma_A0_H_tm});
      vector<distr_t> BR_LIST_V0_tm({ Br_sigma_V0_tm, Br_sigma_V0_H_tm});

           
     
      //OS
      vector<distr_t> BR_LIST_Vii_OS({ Br_sigma_Vii_OS, Br_sigma_Vii_H_OS});
      vector<distr_t> BR_LIST_Aii_OS({ Br_sigma_Aii_OS, Br_sigma_Aii_H_OS});
      vector<distr_t> BR_LIST_A0_OS({ Br_sigma_A0_OS, Br_sigma_A0_H_OS});
      vector<distr_t> BR_LIST_V0_OS({ Br_sigma_V0_OS, Br_sigma_V0_H_OS});

      
      //tm
      distr_t Br_sigma_Vii_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_Vii_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "Vii_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_Aii_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_Aii_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "Aii_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_A0_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_A0_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "A0_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_V0_tm_extr=  Obs_extrapolation_meson_mass(BR_LIST_V0_tm, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "V0_tm_extrapolation_"+ls_data_tm_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );

      //OS
      distr_t Br_sigma_Vii_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_Vii_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "Vii_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_Aii_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_Aii_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "Aii_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_A0_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_A0_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "A0_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
      distr_t Br_sigma_V0_OS_extr=  Obs_extrapolation_meson_mass(BR_LIST_V0_OS, MMK2, Mk_iso*Mk_iso ,  "../data/tau_decay_LIBE_ISO/"+Tag_reco_type+"/strange"  , "V0_OS_extrapolation_"+ls_data_OS_VKVK.Tag[iens]+"_sigma_"+to_string_with_precision(s,3)+".dat",  UseJack, "SPLINE" );
           
      double syst_tm_A0 = Br_sigma_A0_tm_extr.ave()*max(    syst_A0_tm/Br_sigma_A0_tm.ave(),  syst_A0_H_tm/Br_sigma_A0_H_tm.ave() ); // syst_A0_tm;
      double syst_tm_V0 =  Br_sigma_V0_tm_extr.ave()*max(   syst_V0_tm/Br_sigma_V0_tm.ave(), syst_V0_H_tm/Br_sigma_V0_H_tm.ave() ); // syst_V0_tm;
      double syst_tm_Ak =  Br_sigma_Aii_tm_extr.ave()*max(   syst_Aii_tm/Br_sigma_Aii_tm.ave(), syst_Aii_H_tm/Br_sigma_Aii_H_tm.ave() ); // syst_Aii_tm;
      double syst_tm_Vk =  Br_sigma_Vii_tm_extr.ave()*max(  syst_Vii_tm/Br_sigma_Vii_tm.ave(),  syst_Vii_H_tm/Br_sigma_Vii_H_tm.ave() ); // syst_Vii_tm;

      
      double syst_OS_A0 = Br_sigma_A0_OS_extr.ave()*max(  syst_A0_OS/Br_sigma_A0_OS.ave(),  syst_A0_H_OS/Br_sigma_A0_H_OS.ave() ); // syst_A0_OS;
      double syst_OS_V0 =  Br_sigma_V0_OS_extr.ave()*max( syst_V0_OS/Br_sigma_V0_OS.ave(),  syst_V0_H_OS/Br_sigma_V0_H_OS.ave() ); // syst_V0_OS;
      double syst_OS_Ak =  Br_sigma_Aii_OS_extr.ave()*max( syst_Aii_OS/Br_sigma_Aii_OS.ave(), syst_Aii_H_OS/Br_sigma_Aii_H_OS.ave() ); // syst_Aii_OS;
      double syst_OS_Vk =   Br_sigma_Vii_OS_extr.ave()*max(  syst_Vii_OS/Br_sigma_Vii_OS.ave(),  syst_Vii_H_OS/Br_sigma_Vii_H_OS.ave() ); // syst_Vii_OS;


      

      RET.VK_TM = Br_sigma_Vii_tm_extr.ave() + ( Br_sigma_Vii_tm_extr - Br_sigma_Vii_tm_extr.ave())*sqrt( 1.0 + pow( syst_tm_Vk/Br_sigma_Vii_tm_extr.err(),2));
      RET.AK_TM = Br_sigma_Aii_tm_extr.ave() + ( Br_sigma_Aii_tm_extr - Br_sigma_Aii_tm_extr.ave())*sqrt( 1.0 + pow( syst_tm_Ak/Br_sigma_Aii_tm_extr.err(),2));
      RET.V0_TM = Br_sigma_V0_tm_extr.ave() + ( Br_sigma_V0_tm_extr - Br_sigma_V0_tm_extr.ave())*sqrt( 1.0 + pow( syst_tm_V0/Br_sigma_V0_tm_extr.err(),2));
      RET.A0_TM = Br_sigma_A0_tm_extr.ave() + ( Br_sigma_A0_tm_extr - Br_sigma_A0_tm_extr.ave())*sqrt( 1.0 + pow( syst_tm_A0/Br_sigma_A0_tm_extr.err(),2));

      RET.VK_OS = Br_sigma_Vii_OS_extr.ave() + ( Br_sigma_Vii_OS_extr - Br_sigma_Vii_OS_extr.ave())*sqrt( 1.0 + pow( syst_OS_Vk/Br_sigma_Vii_OS_extr.err(),2));
      RET.AK_OS = Br_sigma_Aii_OS_extr.ave() + ( Br_sigma_Aii_OS_extr - Br_sigma_Aii_OS_extr.ave())*sqrt( 1.0 + pow( syst_OS_Ak/Br_sigma_Aii_OS_extr.err(),2));
      RET.V0_OS = Br_sigma_V0_OS_extr.ave() + ( Br_sigma_V0_OS_extr - Br_sigma_V0_OS_extr.ave())*sqrt( 1.0 + pow( syst_OS_V0/Br_sigma_V0_OS_extr.err(),2));
      RET.A0_OS = Br_sigma_A0_OS_extr.ave() + ( Br_sigma_A0_OS_extr - Br_sigma_A0_OS_extr.ave())*sqrt( 1.0 + pow( syst_OS_A0/Br_sigma_A0_OS_extr.err(),2));
    

      cout<<"Ensemble: "<<ls_data_tm_VKVK.Tag[iens]<<", sigma: "<<s<<" completed!"<<endl<<flush;
     
   
     
    
    
      cout<<endl;
          
    }
    
  }
    


  return RET;

}
