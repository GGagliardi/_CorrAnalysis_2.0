#ifndef __numerics__
#define __numerics__

#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <map>
#include <string>
#include <functional>
#include <stdarg.h>
#include <numeric>
#include <sstream>
#include <cassert>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <bits/stdc++.h> 
#include <boost/algorithm/string.hpp> 
#include <ctime>
#include <iomanip>
#include <math.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnMinimize.h>

using namespace std;

typedef boost::minstd_rand base_generator_type;
typedef vector<int> Vint;
typedef vector<vector<int>> VVint;
typedef pair<int,int> Pint;
typedef pair<double,double> Pfloat;
typedef vector<pair<double,double>> VPfloat;
typedef vector<vector<pair<double,double>>> VVPfloat;
typedef vector<pair<int,int>> VPint;
typedef vector<double> Vfloat;
typedef vector<vector<double>> VVfloat;
typedef vector<long double> Vdouble;
typedef vector<vector<pair<int,int>>> VVPint;
typedef vector<vector<vector<vector<double>>>> VVVVfloat;
typedef vector<vector<vector<double>>> VVVfloat;
typedef vector<vector<vector<pair<double,double>>>> VVVPfloat;
typedef vector<vector<vector<vector<pair<double,double>>>>> VVVVPfloat;



void D(int k);
long long int ipow(int a,int n);
double fpow(double a, int n);
long long int fact(int n);
long long int BinomialCoeff(int n,int k);
double RatioPol(Vfloat &x_num, Vfloat &x_den, Vfloat &Num, Vfloat &Den, Vint &NumPow, Vint &DenPow, string MODE_POL);
void crash(string Message);
double eps(int k);
double F(double x, double a, int alpha);
void derivative(Vfloat &RES, Vfloat &INPUT, string MODE);
double FTAN(double x1, int t, int NT);
double FTAN_SYMM(double x1, int t, int NT);
double Root_Brent(double R, int nt, int NT);
double DoConstantFit(Vfloat &data, Vfloat &err);
void Print_To_File(const vector<string>& row_id, const VVfloat &data, string Path, string MODE, string Header);







template <typename T>
string to_string_with_precision(T a_value, const int n)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

template <typename T>
void printV(const vector<T> &A, string B, bool mode) {
  int size = A.size();
  cout.precision(10);
  cout << B << endl;
  if(mode ==0 )for(int i = 0;i < size; i++) cout<< A[i] << " ";
  else for(int i=0;i<size;i++) cout<<i<<"  "<<A[i]<<endl;
  cout << endl;
  return;

}

template <typename T>

T Compute_scalar_product( const vector<T> & A,const vector<T> &B) {

T res =0;

if(A.size() != B.size()) crash("In compute scalar product, size of the two vectors do not coincide");

for (unsigned int i=0; i < A.size(); i++) res += A[i]*B[i];

return res;

}

template <typename T>
vector<T> Multiply_vectors( const vector<T> & A, const vector<T> & B) {

  vector<T> res;
  if(A.size() != B.size()) crash("In Multiply_vectors, the size of the two vectors do not coincide");

  for(unsigned int i=0; i < A.size(); i++) res.push_back(A[i]*B[i]);

  return res;

}


template <typename T>
vector<T> Sum_vectors (const vector<T> & A,const vector<T> &B) {

 vector<T> res;
  if(A.size() != B.size()) crash("In Sum_vectors, the size of the two vectors do not coincide");

  for(unsigned int i=0; i < A.size(); i++) res.push_back(A[i]+B[i]);

  return res;


}

template <typename T>

vector<vector<T>> Sum_Vvectors( const vector<vector<T>> & A,const  vector<vector<T>> &B) {

  vector<vector<T>> res;

  if(A.size() != B.size()) crash("In Multiply_Vvectors, the size of the two vectors do not coincide");


  for(unsigned int i=0; i<A.size();i++) res.push_back( Sum_vectors(A[i], B[i]));

  return res;



}

template <typename T>
vector<T> Multiply_vector_by_scalar(const vector<T> & A, T B) {

  vector<T> res;

  for(auto & el: A) res.push_back(el*B);

  return res;


}

template <typename T>

vector<vector<T>> Multiply_Vvector_by_scalar( const vector<vector<T>> &A, T B) {

  vector<vector<T>> res;

  for( auto & el: A) res.push_back ( Multiply_vector_by_scalar( el, B));

  return res;
}

template <typename T>
void Transpose_VV(vector<vector<T>> &A) {

  if(A.size() == 0) crash("Transpose_VV called with empty matrix.");
  

  vector<vector<T>> B;
  B.resize(A[0].size());
  for(auto &B_i:B) B_i.resize(A.size());

  for(unsigned int i=0; i<A.size(); i++) {
    if(A[i].size() != B.size() ) crash("Matrix A in Transpose_VV has an invalid size");
    for(unsigned int j=0; j<A[i].size(); j++) B[j][i] = A[i][j];
  }



  A=B;

  return;



}


template <typename T>
vector<T> slicing(const vector<T>& arr, int X, int Y) {
  
    // Starting and Ending iterators 
    auto start = arr.begin() + X; 
    auto end = arr.begin() + Y + 1; 
  
    // To store the sliced vector 
    vector<T> result(Y - X + 1); 
  
    // Copy vector using copy function() 
    copy(start, end, result.begin()); 
  
    // Return the final sliced vector 
    return result; 
}


template <typename T>
vector<T> external_prod( const vector<T> &arr1, const vector<T> &arr2) {

  if (arr1.size() != 3 || arr2.size() != 3) crash("Invalid call to external_prod, vectors sizes != 3");

  vector<T> res(3);

  res[0] = arr1[1]*arr2[2] - arr1[2]*arr2[1];
  res[1] = arr1[2]*arr2[0] - arr1[0]*arr2[2];
  res[2] = arr1[0]*arr2[1] - arr1[1]*arr2[0];

  return res;


};



void cascade_resize( vector<vector<double>>& arr, const Vint& A ); 
void cascade_resize( vector<vector<vector<double>>>& arr,const Vint &A);
void cascade_resize( vector<vector<vector<vector<double>>>>& arr,const Vint &A);
void cascade_resize( vector<vector<vector<vector<vector<double>>>>>& arr,const Vint &A);


#endif


