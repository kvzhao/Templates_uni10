#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

UniTensor<double> func1(const UniTensor<double>& A){
  UniTensor<double> B(A.bond());
  B.putBlock(A.getBlock());
  return B;
}

UniTensor<double> func2(const UniTensor<double>& A){
  return func1(A);
}


int main(){

  Matrix<double> buf(2, 3);

  uni10_rand(buf, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

  vector<Bond> bond1(2, Bond(BD_OUT, 3));
  bond1[0] = Bond(BD_IN, 2);

  UniTensor<double> A(bond1);
  A.putBlock(buf);
  cout << A;
  UniTensor<double> C(A);
  cout << C;
  set_zeros(C);
  cout << A;
  cout << C;

  C = func2(A);
  cout << C;

  int newlabel[] = {1, 0};

  vector<int> newlabel1;
  newlabel1.push_back(1);
  newlabel1.push_back(0);

  C = permute(A, newlabel, 0);
  //C = permute(A, 2);

  cout << C;
  

  return 0;
}
