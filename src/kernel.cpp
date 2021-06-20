#include "../include/kernel.h"
using namespace std;

//variables:
//x is normalized square root of photon invariant mass
//y is normalized square root of W invariant mass
//rl and rl2 are the two normalized lepton masses
//mk is the meson mass

//this two logarithms appears in the various kernel functions
double log_plus_minus(double x, double y, double rl){
  double result=log((sqrt(pow(x,4)-2*pow(x,2)*(pow(y,2)+1)+pow((pow(y,2)-1),2))+rl/pow(y,2)*(-sqrt(-2*(pow(x,2)+1)*pow(y,2)+pow((pow(x,2)-1),2)+pow(y,4))+pow(x,2)-1)-pow(x,2)+pow(y,2)+rl-1)/(-sqrt(pow(x,4)-2*pow(x,2)*(pow(y,2)+1)+pow((pow(y,2)-1),2))+rl/pow(y,2)*(sqrt(-2*(pow(x,2)+1)*pow(y,2)+pow((pow(x,2)-1),2)+pow(y,4))+pow(x,2)-1)-pow(x,2)+pow(y,2)+rl-1));
    return result;  
}

double logminus(double x, double y, double rl){
  // double result=log(-sqrt(pow(x,4)-2*pow(x,2)*(pow(y,2)+1)+pow((pow(y,2)-1),2))+rl/pow(y,2)*(sqrt(-2*(pow(x,2)+1)*pow(y,2)+pow((pow(x,2)-1),2)+pow(y,4))+pow(x,2)-1)-pow(x,2)+pow(y,2)+rl-1);
  //return result;
  return 0;
}



//point-like contribution to the rate
double ptrate(double x, double y, double rl, double rl2, double mk, double fk){
    double result=pow(fk,2)*pow(mk,3)*rl/96/pow(M_PI,3)/pow(x,2)*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*(2/(pow(y,2)-1)*(-pow(x,2)*pow(y,2)+pow(x,2)+pow(y,4)-2*pow(y,2)*rl-2*(rl-1)*rl+1)*(log_plus_minus(x,y,rl))+sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(y,2)-rl)*((pow(x,2)*pow(y,2)-pow(x,2)*rl-2*pow(y,4)+4*pow(y,2)*rl-2)/pow(y,2)/pow((pow(y,2)-1),2)+2*(rl-1)*(pow(x,2)+2*rl)/(pow((pow(y,2)-1),2)*rl-pow(x,2)*(rl-1)*(pow(y,2)-rl))));
    return result;
}




//interaction kernels
double kernV(double x, double y, double rl, double rl2, double mk, double fk){
    double result=fk*pow(mk,4)*rl/(48*pow(M_PI,3)*pow(x,2))*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*((pow(x,2)*(pow(y,2)-2*rl+1)-pow((pow(y,2)-1),2))*(log_plus_minus(x,y,rl))+sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(x,2)+pow(y,2)-1)*(pow(y,2)-rl)/pow(y,2));
    return result;
}

double kernA(double x, double y, double rl, double rl2, double mk, double fk){
    double result=fk*pow(mk,4)*rl/(48*pow(M_PI,3)*pow(x,2)*(pow(y,2)-1))*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*(pow((pow(y,2)-1),2)*(-pow(x,2)-pow(y,2)-2*rl+1)*(log_plus_minus(x,y,rl))+(pow(y,2)-1)*sqrt(-2*(pow(x,2)+1)*pow(y,2)+pow((pow(x,2)-1),2)+pow(y,4))*(pow(y,2)-rl)*(pow(x,2)+2*pow(y,2)+rl-1)/pow(y,2));
    return result;
}

double kern1(double x, double y, double rl, double rl2, double mk, double fk){
    double result=fk*pow(mk,4)*rl/(96*pow(M_PI,3))*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*(4*(pow(y,2)+rl-2)*(log_plus_minus(x,y,rl))-sqrt(-2*(pow(x,2)+1)*pow(y,2)+pow((pow(x,2)-1),2)+pow(y,4))*(pow(y,2)-rl)*(rl*(-pow(x,2)+3*pow(y,2)+1)+pow(y,2)*(pow(x,2)+5*pow(y,2)-9))/(pow(y,4)*(pow(y,2)-1)));
    return result;
}

double kern2(double x, double y, double rl, double rl2, double mk, double fk){
    double result=fk*pow(mk,4)*rl/(96*pow(M_PI,3)*pow((pow(y,2)-1),2))*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*(2*(pow(y,2)-1)*(pow(y,2)-pow(rl,2))*(log_plus_minus(x,y,rl))-sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(y,2)-rl)*(rl*(pow(x,2)-3*pow(y,2)-1)+pow(y,2)*(-pow(x,2)+pow(y,2)+3))/pow(y,2));
    return result;    
}



//SD kernels

double kernVV(double x, double y, double rl, double rl2, double mk){
    double result=pow(mk,5)/(96*pow(M_PI,3)*pow(x,4)*pow(y,2))*pow(sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2)),3)*(pow(x,2)+2*rl2)*sqrt(1-4*rl2/pow(x,2))*pow((pow(y,2)-rl),2);
    return result;
}

double kernAA(double x, double y, double rl, double rl2, double mk){
    double result=pow(mk,5)/(576*pow(M_PI,3)*pow(x,4)*pow(y,4))*sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(x,4)+pow(x,2)*(4*pow(y,2)-2)+pow((pow(y,2)-1),2))*(pow(x,2)+2*rl2)*sqrt(1-4*rl2/pow(x,2))*pow((pow(y,2)-rl),2)*(2*pow(y,2)+rl);
    return result;
}

double kern11(double x, double y, double rl, double rl2, double mk){
    double result=pow(mk,5)/(576*pow(M_PI,3)*pow(x,2)*pow(y,6))*sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(x,2)+2*rl2)*sqrt(1-4*rl2/pow(x,2))*pow((pow(y,2)-rl),2)*(2*pow(y,4)*(5*pow(x,2)+rl-1)+pow(y,2)*(2*(pow(x,2)-2)*rl+pow((pow(x,2)-1),2))+2*pow((pow(x,2)-1),2)*rl+pow(y,6));
    return result;
}

double kern22(double x, double y, double rl, double rl2, double mk){
    double result=pow(mk,5)*rl/(384*pow(M_PI,3)*pow((pow(y,2)-1),2)*pow(y,2))*pow(sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2)),3)*sqrt(1-4*rl2/pow(x,2))*(2*rl2/pow(x,2)+1)*pow((pow(y,2)-rl),2);
    return result;
}

double kernA1(double x, double y, double rl, double rl2, double mk){
    double result=-pow(mk,5)/(96*pow(M_PI,3)*pow(x,2)*pow(y,4))*sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2))*(pow(x,2)+pow(y,2)-1)*(pow(x,2)+2*rl2)*sqrt(1-4*rl2/pow(x,2))*pow((pow(y,2)-rl),2)*(2*pow(y,2)+rl);
    return result;
}

double kern12(double x, double y, double rl, double rl2, double mk){
    double result=-pow(mk,5)*rl/(192*pow(M_PI,3)*pow(x,2)*(pow(y,2)-1)*pow(y,4))*pow(sqrt(pow((pow(x,2)-pow(y,2)+1),2)-4*pow(x,2)),3)*(pow(x,2)+2*rl2)*sqrt(1-4*rl2/pow(x,2))*pow((pow(y,2)-rl),2);
    return result;    
}
