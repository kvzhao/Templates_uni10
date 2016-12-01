#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  int dim1 = 100, dim2 = 120;
  ///A is matrix to be decomposed
  Matrix<double> A(dim1, dim2);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  ///u, s, vt are answers
  Matrix<double> l, q;
  lq( A, l, q, INPLACE);
  vector<Matrix<double>> lqs = lq(A);
  ///q_r multiplication of answers
  Matrix<double> l_q_inplace;
  dots( l_q_inplace, l, q);
  Matrix<double> l_q_noninplace;
  dots( l_q_noninplace, lqs.at(0), lqs.at(1));
  ///output error
  cout<<l<<endl;
  cout<<q<<endl;
  cout<<"error = "<<norm( A - l_q_noninplace)<<endl;
  cout<<"error = "<<norm( A - l_q_inplace)<<endl;

  return 0;
}
