#include "../include/Num_integrate_l4_decay_rate.h"

using namespace std;


const double r_mu = pow(0.10565837/0.493677,2); //squared ratio m_mu / m_K
const double r_el = pow( 0.000510998950/0.493677,2); //squared ratio m_el / mK
const double Gf = 1.1663787*1e-5; //Fermi constant in GeV^-2
const double Vus=  0.2243; //Vus CKM matrix element
const double Gamma = 5.317*1e-17; // Total decay width of the K^+ in GeV
const double alpha = 1/137.04; //alpha constant
const double eps_mumu= 1e-10;
const double eps_ee = 1e-10;
const double MkPh=  0.493677; //GeV
const double fkPh = 0.155 ;//GeV
double cut =  pow(0.140/MkPh,2); //pow(0.020/MkPh,2); //pow(0.140/MkPh,2); //pow(0.020/MkPh,2); //0.001; // 
double error_goal_MC_ee = 1000*1.0e-2*(Gamma/pow(Gf*Vus*alpha,2))*8.0*1.0e-8;
double error_goal_MC_mumu= 1000*1.0e-2*(Gamma/pow(Gf*Vus*alpha,2))*1.0e-8;


Decay_Rate_Integration_Result Num_Integrate_Decay_Rate(const vector<function<double(double, double)>> &H1, const vector<function<double(double, double)>> &H2  , const vector<function<double(double, double)>> &FA, const vector<function<double(double, double)>> &FV, distr_t &m_distr, distr_t & fp_distr, bool UseJack) {

  

  
  if( (H1.size() != H2.size()) || (H1.size() != FA.size()) || (H1.size() != FV.size()) || (H1.size() != (unsigned)m_distr.size()) || (H1.size() != (unsigned)fp_distr.size())) crash("Sizes of form factors in Decay_Rate_Integration_Result are different.");

  int Njacks= H1.size();
  Decay_Rate_Integration_Result Res(UseJack);

  double eps;

  //relative accuracy in integration
  Res.eps_rel_mumu= eps_mumu;
  Res.eps_rel_ee = eps_ee;

  //loop over jackknives
  for(int ijack=0; ijack< Njacks; ijack++) {
    //Compute decay rates for a given jackknife

    double m= m_distr.distr[ijack];
    double fp= fp_distr.distr[ijack];
    //set to physical to make test
    m=MkPh;
    fp=fkPh;
    

    //define double differential decay rate
    auto Rt_diff = [&](double xk, double xq, double l, double ll) -> double {

      double Int= ptrate(xk,xq, l, ll, m, fp) + H1[ijack](xk,xq)*kern1(xk, xq, l, ll, m,fp) + H2[ijack](xk,xq)*kern2(xk, xq, l, ll, m,fp) + FA[ijack](xk,xq)*kernA(xk,xq,l,ll,m,fp) + FV[ijack](xk,xq)*kernV(xk,xq,l,ll,m,fp) + pow(H1[ijack](xk,xq),2)*kern11(xk,xq,l,ll,m) + pow(H2[ijack](xk,xq),2)*kern22(xk,xq,l,ll,m) + pow(FA[ijack](xk,xq),2)*kernAA(xk,xq,l,ll,m) + pow(FV[ijack](xk,xq),2)*kernVV(xk,xq,l,ll,m) + H1[ijack](xk,xq)*H2[ijack](xk,xq)*kern12(xk,xq,l,ll,m) + H1[ijack](xk,xq)*FA[ijack](xk,xq)*kernA1(xk,xq,l,ll,m);
      double jacobian= 4.0*xk*xq;
      return Int*jacobian;
    };

    cout<<"#####Quad integration#####"<<endl;
    //e+e-

    double rl=r_mu;
    double rll=r_el;

    /*
    for(double xk=2*sqrt(rll);xk<1-sqrt(rl); xk+=0.001)
      for(double xq=sqrt(rl);xq<1-xk;xq+=0.001) {
	cout<<"kern1("<<xk<<","<<xq<<"): "<<kern1(xk,xq,rl,rll,m,fp)<<endl;
	cout<<"kern2("<<xk<<","<<xq<<"): "<<kern2(xk,xq,rl,rll,m,fp)<<endl;
	cout<<"kernA("<<xk<<","<<xq<<"): "<<kernA(xk,xq,rl,rll,m,fp)<<endl;
	cout<<"kernV("<<xk<<","<<xq<<"): "<<kernV(xk,xq,rl,rll,m,fp)<<endl;
	cout<<"kern11("<<xk<<","<<xq<<"): "<<kern11(xk,xq,rl,rll,m)<<endl;
	cout<<"kern22("<<xk<<","<<xq<<"): "<<kern22(xk,xq,rl,rll,m)<<endl;
	cout<<"kernAA("<<xk<<","<<xq<<"): "<<kernAA(xk,xq,rl,rll,m)<<endl;
	cout<<"kernVV("<<xk<<","<<xq<<"): "<<kernVV(xk,xq,rl,rll,m)<<endl;
	cout<<"kern12("<<xk<<","<<xq<<"): "<<kern12(xk,xq,rl,rll,m)<<endl;
	cout<<"kernA1("<<xk<<","<<xq<<"): "<<kernA1(xk,xq,rl,rll,m)<<endl;
	cout<<"log+("<<xk<<","<<xq<<"): "<<logplus(xk,xq,rl)<<endl;
	cout<<"log-("<<xk<<","<<xq<<"): "<<logminus(xk,xq,rl)<<endl;

      }
    */
    eps=eps_ee;

    auto Fxk= [&](double xk) -> double {
      double inf = sqrt(rl);
      double sup = 1-xk;
      auto g = [&](double xq) {
	return Rt_diff(xk, xq, rl, rll);
        };
      return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, inf, sup, 5, eps);
    };

    double error_ee;
    double res_ee= (pow(alpha*Gf*Vus,2)/Gamma)*boost::math::quadrature::gauss_kronrod<double,15>::integrate(Fxk, sqrt(cut), 1-sqrt(rl), 5, eps, &error_ee);
    Res.Jack_Distr_Int_Quad_ee.distr.push_back(res_ee );
    cout<<"Relative error achieved in Quad e+e- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_ee<<endl;
    cout<<"Branching ratio Quad e+e-: "<<res_ee<<endl;
    cout<<"################"<<endl;
    
    //mu+mu-

    rl=r_el;
    rll=r_mu;
    eps=eps_mumu;
    

    double error_mm;
    double res_mm = (pow(alpha*Gf*Vus,2)/Gamma)*boost::math::quadrature::gauss_kronrod<double,15>::integrate(Fxk, 2*sqrt(rll), 1-sqrt(rl), 5, eps, &error_mm);
    Res.Jack_Distr_Int_Quad_mumu.distr.push_back(res_mm);
    cout<<"Relative error achieved in Quad mu+mu- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_mm<<endl;
    cout<<"Branching ratio Quad  mu+mu-: "<<res_mm<<endl;
    cout<<"################"<<endl;
    //perform MonteCarlo Integration
    cout<<"#####Monte Carlo integration#####"<<endl;
    
    //define lambda function for Monte Carlo Integration


    auto Rt_diff_MC = [&](Vfloat const &par) -> double {
      if(par.size() != 2) crash("Rt_diff_MC in Num_Integrate_Decay_Rate only accepts vector<double> of size 2. "+to_string( (signed)par.size())+" provided.");
      double xk= par[0];
      double xq= par[1];
      if(xq +xk >= 1) return 0;
      
      return Rt_diff(xk,xq,rl,rll);
    };



    //e+e-

    rl=r_mu;
    rll=r_el;
    VPfloat bounds_ee{ {sqrt(cut),1.0-sqrt(rl)}, { sqrt(rl),1.0-sqrt(cut)}};

    
    boost::math::quadrature::naive_monte_carlo<double, decltype(Rt_diff_MC)> mc_ee(Rt_diff_MC, bounds_ee, error_goal_MC_ee,  true, thread::hardware_concurrency() - 1);

    
    future<double> task_ee = mc_ee.integrate();
    while (task_ee.wait_for(chrono::seconds(1)) != future_status::ready)
      {
	// do nothing
      }

    double res_ee_MC = (pow(alpha*Gf*Vus,2)/Gamma)*task_ee.get();
    Res.Jack_Distr_Int_MonteCarlo_ee.distr.push_back(res_ee_MC);
    Res.Stat_err_MonteCarlo_ee.push_back( (pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_ee/res_ee_MC);
    cout<<"Relative error achieved in Monte Carlo e+e- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_ee/res_ee_MC<<endl;
    cout<<"Branching ratio Monte Carlo e+e-: "<<res_ee_MC<<endl;
    cout<<"################"<<endl;


    //mu+mu-

    rl=r_el;
    rll=r_mu;

    VPfloat bounds_mumu{ {2.0*sqrt(rll), 1-sqrt(rl)},  {sqrt(rl) , 1.0- 2.0*sqrt(rll) }};
    
    boost::math::quadrature::naive_monte_carlo<double, decltype(Rt_diff_MC)> mc_mumu(Rt_diff_MC, bounds_mumu, error_goal_MC_mumu,  true, thread::hardware_concurrency() - 1);

    
    future<double> task_mumu = mc_mumu.integrate();
    while (task_mumu.wait_for(chrono::seconds(1)) != future_status::ready)
      {
	// do nothing
      }

    double res_mumu_MC = (pow(alpha*Gf*Vus,2)/Gamma)*task_mumu.get();
    Res.Jack_Distr_Int_MonteCarlo_mumu.distr.push_back(res_mumu_MC);
    Res.Stat_err_MonteCarlo_mumu.push_back(error_goal_MC_mumu/res_mumu_MC);
    cout<<"Relative error achieved in Monte Carlo mu+mu- for jack: "<<ijack<<" is: "<<(pow(alpha*Gf*Vus,2)/Gamma)*error_goal_MC_mumu/res_mumu_MC<<endl;
    cout<<"Branching ratio Monte Carlo e+e-: "<<res_mumu_MC<<endl;
    cout<<"################"<<endl;


  }


  //compute expectation value and std err
  Res.Int_Quad_val_ee= Res.Jack_Distr_Int_Quad_ee.ave();
  Res.Int_Quad_val_mumu= Res.Jack_Distr_Int_Quad_mumu.ave();
  Res.Int_MonteCarlo_val_ee = Res.Jack_Distr_Int_MonteCarlo_ee.ave();
  Res.Int_MonteCarlo_val_mumu= Res.Jack_Distr_Int_MonteCarlo_mumu.ave();
  Res.Int_Quad_err_ee= Res.Jack_Distr_Int_Quad_ee.err();
  Res.Int_Quad_err_mumu= Res.Jack_Distr_Int_Quad_mumu.err();
  Res.Int_MonteCarlo_err_ee = Res.Jack_Distr_Int_MonteCarlo_ee.err();
  Res.Int_MonteCarlo_err_mumu= Res.Jack_Distr_Int_MonteCarlo_mumu.err();


  return Res;
}
