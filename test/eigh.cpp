#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

int main(){

  double elem_B[16] = {
    1, 2, 3, 4,
    2, 4, 9, 1,
    3, 9, 7, 3,
    4, 1, 3, 0
  };

  Matrix<double> A(4, 4);
  uni10_rng(A, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<double> B(4, 4);
  B.setElem(elem_B);
  //uni10_rng(B, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);


  Matrix<uni10_complex128> Eig, EigVec;
  
  eig(A, Eig, EigVec, INPLACE);

  //vector< Matrix<uni10_complex128> > EV = eig(B);
  vector< Matrix<double> > EV = eigh(B);

  cout << EV[0];
  cout << EV[1];

  //cout << Eig;
  //cout << EigVec;

  return 0;
}
