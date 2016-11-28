#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

UniTensor<double> func1(const UniTensor<double>& A){
  UniTensor<double> B(A);
  return B;
}

UniTensor<double> func2(const UniTensor<double>& A){
  return func1(A);
}


int main(){

  Matrix<double> buf(3, 9);
  Matrix<double> buf2(3, 9);


  uni10_rand(buf, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);
  uni10_rand(buf2, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

  //cout << buf;
  //cout << buf2;


  vector<Bond> bonds(3, Bond(BD_OUT, 3));
  bonds[0] = Bond(BD_IN, 3);
  //UniTensor<double> B(2);
  int label[] = {1, 4, 5};
  UniTensor<double> B(bonds, label);
  B.putBlock(buf);
  cout << B;

  //cout << B.elemNum()<< std::endl;
  //cout << B.getElem()<< std::endl;
  UniTensor<double> C(B);
  int labelC[] = {2, 3, 4};

  //cout << C.elemNum()<< std::endl;
  //cout << B.getElem()<< std::endl;
  //cout << C.getElem()<< std::endl;
  //exit(0);

  C.putBlock(buf2);
  C.setLabel(labelC);

  cout << B;
  cout << C;

  //cout << B;
  //cout << C;

  UniTensor<double> E;
  //E = contract(B, C, true);
  E = otimes(B, C);
  E.printDiagram();
  return 0;
}
