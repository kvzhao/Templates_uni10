#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){
  ///create symmetric A
  Matrix<double> A( 8, 8);
  uni10_rand(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<double> At = transpose( A);
  A = A + At;
  ///get exp( A ) = U exp(d) Ut
  Matrix<double> exp_A = exph( 1.0, A);
  cout<<"exp(A)"<<exp_A;
  ///gat d, ut and check
  Matrix<double> d, u, ut;
  eigh( A, d, u, INPLACE);
  ut = dagger( u);
  ///get exp(d)
  Matrix<double> ans;
  dots( ans, u, exp_A, ut);
  cout<<"d"<<d;
  cout<<"exp(d)"<<ans;

  return 0;
}
