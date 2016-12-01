#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  int dim1 = 11, dim2 = 10;
  ///A is matrix to be decomposed
  Matrix<double> A(dim1, dim2);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  ///u, s, vt are answers
  Matrix<double> q, r;
  qr( A, q, r, INPLACE);
  vector<Matrix<double>> qrs = qr(A);
  ///q_r multiplication of answers
  Matrix<double> q_r_inplace;
  dots( q_r_inplace, q, r);
  Matrix<double> q_r_noninplace;
  dots( q_r_noninplace, qrs.at(0), qrs.at(1));
  ///output error
  cout<<q<<endl;
  cout<<r<<endl;
  cout<<"error = "<<norm( A - q_r_noninplace)<<endl;
  cout<<"error = "<<norm( A - q_r_inplace)<<endl;

  return 0;
}
