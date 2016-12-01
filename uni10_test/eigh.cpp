#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  int dim1 = 7;
  ///A is matrix to be decomposed
  Matrix<double> A(dim1, dim1);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  ///d, u are answers
  Matrix<uni10_complex128> d, u;
  eig( A, d, u, INPLACE);
  cout<<d<<endl;
  cout<<u<<endl;

  Matrix<double> ut = u;
  ut.ctranspose();
  cout<<A;
  cout<<dot( d, u);

  return 0;
}
