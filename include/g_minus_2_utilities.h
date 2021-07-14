#ifndef  __g_minus_2_utilities__
#define  __g_minus_2_utilities__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "LatInfo.h"
#include "T_min.h"
#include "input.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/finite_difference.hpp>


using namespace std;

bool Is_perfect_square(int x);
double kernel_K(double t,double MV, string channel);
double Vdual(double t, double m_rho, double Edual, double Rdual);
double Gamma_rpp(double omega, double g_rho_pipi, double Mpi);
double Gamma_rpp_der(double omega, double g_rho_pipi, double Mpi);
double h(double omega, double g_rho_pipi, double Mpi);
double h_prime(double omega, double g_rho_pipi, double Mpi);
double A_pipi_0(double m_rho, double g_rho_pipi, double Mpi);
double Re_A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi);
double Im_A_pipi(double omega, double g_rho_pipi, double Mpi);
Pfloat A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi);
double F_pi_GS_mod(double omega, double m_rho, double g_rho_pipi, double Mpi);
double cot_d_11(double k, double m_rho, double g_rho_pipi, double Mpi);
double d_11(double k, double m_rho, double g_rho_pipi, double Mpi);
double cot_d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi);
double d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi);
int degeneracy(int m);
double Zeta_function_laplacian(double z);
double Zeta_function_laplacian_Luscher(double z);
double tan_phi(double z);
double phi(double z);
double phi_der(double z);
double Find_pipi_energy_lev(int L, int n, double m_rho, double g_rho_pipi, double Mpi);
double V_pipi(double t, int L, double m_rho, double g_rho_pipi, double Mpi);
double V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL);








#endif
